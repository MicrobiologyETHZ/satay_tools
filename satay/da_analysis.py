from pydeseq2.ds import DeseqStats
from pydeseq2.default_inference import DefaultInference
from pydeseq2.dds import DeseqDataSet
import pandas as pd
import pyranges as pr
from pathlib import Path
from datetime import datetime


def run_deseq(count_df, sample_data, condition_col, baseline, a=0.01):
    """
    Make sure count_df index is sorted, as well as sample_data
    Make sure no duplicated entries in count_df
    """
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
    sample_data = pd.read_csv(sample_data_file)
    df = pd.read_csv(counts_file, index_col=0)
    df = df[sample_data[sample_id_col].unique()]
    df = df.T.sort_index()
    sample_data = sample_data.set_index(sample_id_col).sort_index()
    genes_to_keep = df.columns[df.sum(axis=0) >= filter]
    df = df[genes_to_keep]
    # check if indices are the same
    vst_counts, res_df = run_deseq(
        df, sample_data, comp_col, baseline,  a=a)
    # save
    return vst_counts, res_df


def merge_annotations(res_df, gff_file, ids_to_keep=["locus_tag", "gene", "product"]):
    gff = pr.read_gff3(gff_file).as_df()
    gff = gff[ids_to_keep].drop_duplicates()
    return res_df.merge(gff, left_on='ID', right_on='locus_tag', how='left')


def gene_da(counts_file, sample_data_file, output_dir,
            filter, comp_col, baseline, a,
            gff_file='', ids_to_keep=["locus_tag", "gene", "product"]):
    today_str = datetime.today().strftime("%y-%m-%d")
    vst_counts, res_df = load_filter_deseq(counts_file, sample_data_file,
                                           filter,
                                           comp_col, baseline, a)
    print(res_df.shape)
    if gff_file:
        res_df = merge_annotations(res_df, gff_file, ids_to_keep)
    print(res_df.shape)
    res_df.to_csv(
        output_dir/f"{today_str}_{counts_file.stem}_{comp_col}_l0a{a}.csv", index=False)
    vst_counts.to_csv(
        output_dir/f"{today_str}_{counts_file.stem}_vstcounts.csv")
