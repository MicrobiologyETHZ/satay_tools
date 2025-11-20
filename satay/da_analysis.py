from pydeseq2.ds import DeseqStats
from pydeseq2.default_inference import DefaultInference
from pydeseq2.dds import DeseqDataSet
import pandas as pd
import pyranges as pr
from pathlib import Path
from datetime import datetime
import logging
import sys


def run_deseq(count_df, sample_data, condition_col, baseline, a=0.01):
    """
    Run DESeq2 analysis on count data.
    
    Args:
        count_df: DataFrame with count data (samples as rows, genes as columns)
        sample_data: DataFrame with sample metadata
        condition_col: Column name for condition/treatment in sample_data
        baseline: Baseline condition for comparisons
        a: Alpha value for significance testing
    
    Returns:
        tuple: (vst_counts DataFrame, results DataFrame)
    """
    logging.info(f"Running DESeq2 analysis with baseline: {baseline}")
    
    # Validate inputs
    if baseline not in sample_data[condition_col].values:
        raise ValueError(f"Baseline '{baseline}' not found in column '{condition_col}'")
    
    if not count_df.index.equals(sample_data.index):
        logging.warning("Count data and sample data indices don't match exactly")
    
    inference = DefaultInference(n_cpus=8)
    dds = DeseqDataSet(
        counts=count_df,
        metadata=sample_data,
        design_factors=condition_col,
        refit_cooks=True,
        inference=inference,
    )
    dds.deseq2()
    dds.vst()
    vst_counts = pd.DataFrame(
        dds.layers["vst_counts"], index=count_df.index, columns=count_df.columns
    ).T
    res_l = []
    comps = list(sample_data[condition_col].unique())
    comps.remove(baseline)    
    logging.info(f"Performing {len(comps)} comparisons against baseline '{baseline}'")
    for comp in comps:
        res = DeseqStats(
            dds, contrast=[condition_col, comp,
                           baseline], inference=inference, alpha=a
        )
        res.summary()
        res_l.append(res.results_df.assign(contrast=f"{comp}_vs_{baseline}"))
    res_df = pd.concat(res_l)
    return vst_counts, res_df


def load_filter_deseq(counts_file, sample_data_file,
                      filter=100,
                      comp_col='conc', baseline="0nMaF",
                      a=0.01,
                      sample_id_col='sample_id'
                      ):
    """Load count data and sample metadata, filter low-count genes, run DESeq2."""
    try:
        logging.info(f"Loading sample data from {sample_data_file}")
        sample_data = pd.read_csv(sample_data_file)
        
        logging.info(f"Loading count data from {counts_file}")
        df = pd.read_csv(counts_file, index_col=0)
        
        # Validate required columns exist
        if sample_id_col not in sample_data.columns:
            raise ValueError(f"Sample ID column '{sample_id_col}' not found in sample data")
        if comp_col not in sample_data.columns:
            raise ValueError(f"Comparison column '{comp_col}' not found in sample data")
        
        # Filter to samples present in both files
        metadata_samples = set(sample_data[sample_id_col].unique())
        count_samples = set(df.columns)
        
        # Find samples missing from each file
        missing_from_counts = metadata_samples - count_samples
        missing_from_metadata = count_samples - metadata_samples
        
        if missing_from_counts:
            logging.warning(f"Samples in metadata but not in count data: {sorted(missing_from_counts)}")
        if missing_from_metadata:
            logging.warning(f"Samples in count data but not in metadata: {sorted(missing_from_metadata)}")
        
        # Use only samples present in both files
        common_samples = metadata_samples & count_samples
        if not common_samples:
            raise ValueError("No common samples found between count data and metadata")
        
        logging.info(f"Using {len(common_samples)} samples present in both files")
        df = df[list(common_samples)]
        sample_data = sample_data[sample_data[sample_id_col].isin(common_samples)]
        df = df.T.sort_index()
        sample_data = sample_data.set_index(sample_id_col).sort_index()
        
        # Filter low-count genes
        initial_gene_count = df.shape[1]
        genes_to_keep = df.columns[df.sum(axis=0) >= filter]
        df = df[genes_to_keep]
        logging.info(f"Filtered genes: {initial_gene_count} -> {df.shape[1]} (removed {initial_gene_count - df.shape[1]} low-count genes)")
        
        if df.empty:
            raise ValueError("No genes remain after filtering")
        
        vst_counts, res_df = run_deseq(df, sample_data, comp_col, baseline, a=a)
        return vst_counts, res_df
        
    except FileNotFoundError as e:
        logging.error(f"File not found: {e}")
        raise
    except Exception as e:
        logging.error(f"Error in load_filter_deseq: {e}")
        raise


def merge_annotations(res_df, gff_file, ids_to_keep=["locus_tag", "gene", "product"]):
    """Merge DESeq2 results with gene annotations from GFF file."""
    try:
        logging.info(f"Loading annotations from {gff_file}")
        gff = pr.read_gff3(gff_file).as_df()
        
        # Check if required columns exist
        missing_cols = set(ids_to_keep) - set(gff.columns)
        if missing_cols:
            logging.warning(f"Missing annotation columns: {missing_cols}")
            ids_to_keep = [col for col in ids_to_keep if col in gff.columns]
        
        if not ids_to_keep:
            logging.warning("No valid annotation columns found")
            return res_df
            
        gff = gff[ids_to_keep].drop_duplicates()
        merged = res_df.merge(gff, left_on='ID', right_on='locus_tag', how='left')
        logging.info(f"Merged annotations: {len(res_df)} -> {len(merged)} rows")
        return merged
        
    except Exception as e:
        logging.error(f"Error merging annotations: {e}")
        return res_df


def gene_da(counts_file, sample_data_file, output_dir,
            filter, comp_col, baseline, a,
            gff_file='', ids_to_keep=["locus_tag", "gene", "product"]):
    """Run complete differential analysis pipeline."""
    # Set up logging
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
    
    try:
        # Convert string paths to Path objects
        counts_file = Path(counts_file)
        sample_data_file = Path(sample_data_file)
        output_dir = Path(output_dir)
        
        # Ensure output directory exists
        output_dir.mkdir(parents=True, exist_ok=True)
        
        logging.info("Starting differential analysis pipeline")
        today_str = datetime.today().strftime("%y-%m-%d")
        
        vst_counts, res_df = load_filter_deseq(counts_file, sample_data_file,
                                               filter, comp_col, baseline, a)
        logging.info(f"DESeq2 results shape: {res_df.shape}")
        
        if gff_file:
            res_df = merge_annotations(res_df, gff_file, ids_to_keep)
            logging.info(f"Results with annotations shape: {res_df.shape}")
        
        # Save results
        results_file = output_dir / f"{today_str}_{counts_file.stem}_{comp_col}_l0a{a}.csv"
        vst_file = output_dir / f"{today_str}_{counts_file.stem}_vstcounts.csv"
        
        res_df.to_csv(results_file, index=False)
        vst_counts.to_csv(vst_file)
        
        logging.info(f"Results saved to: {results_file}")
        logging.info(f"VST counts saved to: {vst_file}")
        logging.info("Differential analysis complete")
        
    except Exception as e:
        logging.error(f"Pipeline failed: {e}")
        sys.exit(1)
