import os
import subprocess
from pathlib import Path
import argparse
import logging
import pandas as pd
from datetime import datetime
# import pyranges as pr
from typing import Tuple, Union, List
import shutil


def setup_logging(output_dir):
    os.makedirs(output_dir, exist_ok=True)
    logging.basicConfig(
        filename=Path(output_dir)/"process_samples.log",
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s"
    )


def validate_environment():
    # Validate bedtools installation
    if not shutil.which("bedtools"):
        raise RuntimeError("bedtools is not installed or not in PATH")
    if not shutil.which("samtools"):
        raise RuntimeError("samtools is not installed or not in PATH")


# def get_samples(directory):
#     """Get unique sample identifiers from the directory."""
#     try:
#         sample_dirs = [d for d in os.listdir(directory) if os.path.isdir(os.path.join(directory, d))]
#         samples = sorted(set("_".join(d.split("_")[:2]) for d in sample_dirs if d.startswith("20")))
#         logging.info(f"Identified samples: {samples}")
#         return samples
#     except Exception as e:
#         logging.error(f"Error getting samples: {e}")
#         raise


def run_command(cmd):
    """Run a subprocess command with error handling."""
    try:
        subprocess.run(cmd, check=True)
        logging.info(f"Command succeeded: {' '.join(cmd)}")
    except subprocess.CalledProcessError as e:
        logging.error(f"Command failed: {' '.join(cmd)}\nError: {e}")
        raise


def merge_bams(sample, output_dir, bam_dir, bam_suffix="Aligned.sortedByCoord.out.bam"):
    """Merge BAM files for each sample."""

    bam_files = [str(bf) for bf in Path(bam_dir).rglob(
        f"*{sample}*/*{bam_suffix}")]
    if bam_files:
        logging.info(f"Merging {', '.join(bam_files)} with samtools.")
        output_bam = Path(output_dir)/f"{sample}.bam"
        cmd = ["samtools", "merge", "-f", "-o", str(output_bam)] + \
            [str(bam) for bam in bam_files]
        run_command(cmd)
        return output_bam
    else:
        msg = (f"No bam files found for sample '{sample}' in {bam_dir} "
               f"matching '*{sample}*/*{bam_suffix}'")
        logging.error(msg)
        raise FileNotFoundError(msg)


def filter_bam(output_bam):
    """Filter BAM files and convert to BED format.
        Filter out q < 10
    """
    output_bam = output_bam
    bed_file = output_bam.with_suffix(".bed")
    cmd = ["samtools", "view", "-h", "-F",
           "2304", "-q", "10", str(output_bam)]
    try:
        with subprocess.Popen(cmd, stdout=subprocess.PIPE) as proc:
            with open(bed_file, "w") as bed_out:
                subprocess.run(["bedtools", "bamtobed"],
                               stdin=proc.stdout, stdout=bed_out, check=True)
        logging.info(f"Filtered BAM to BED: {output_bam} -> {bed_file}")
        if bed_file.stat().st_size == 0:
            raise ValueError(
                f"No reads passed filtering (primary alignments, MAPQ >= 10) for "
                f"{output_bam}; the resulting BED is empty. Verify the BAM contains "
                f"mapped reads, e.g. 'samtools view -c {output_bam}'. A failed or empty "
                f"alignment step is the usual cause."
            )
        return bed_file
    except Exception as e:
        logging.error(f"Error filtering BAM {output_bam}: {e}")
        raise


def process_orientation(filtered_bam_bed):
    """Process BED files based on orientation."""
    insertions_file = Path(filtered_bam_bed).with_suffix(".bed.insertions")
    # cmd = [script_path, str(filtered_bam_bed)]
    try:
        # Read the file assuming whitespace-delimited columns without a header.
        df = pd.read_table(filtered_bam_bed, header=None)
    except Exception as e:
        logging.error(f"Error reading file: {e}")
        raise
    # Ensure there are at least 6 columns in the input.
    if df.shape[1] < 6:
        logging.error("Error: The input file must contain at least 6 columns.")
        raise ValueError(
            "DataFrame must have at least 6 columns, but got {} columns.".format(df.shape[1]))

    # Use only the first six columns and assign column names.
    df = df.iloc[:, :6]
    df.columns = ["col1", "col2", "col3", "col4", "col5", "col6"]
    # Convert columns used for arithmetic to numeric (col2 and col3)
    df["col2"] = pd.to_numeric(df["col2"], errors="coerce")
    df["col3"] = pd.to_numeric(df["col3"], errors="coerce")

    # Compute the insertion site: if strand is "+" use col2, else use col3.
    df["insertion_site"] = df.apply(
        lambda row: row["col2"] if row["col6"] == "+" else row["col3"],
        axis=1
    )
    # Compute insertion_site + 1.
    df["insertion_site_plus_1"] = df["insertion_site"] + 1
    # Select the desired columns:
    # col1, insertion_site, insertion_site+1, col4, col5, col6
    result = df[["col1", "insertion_site",
                 "insertion_site_plus_1", "col4", "col5", "col6"]]
    # Output the result as tab-separated values to stdout.
    result.to_csv(insertions_file, sep="\t", index=False, header=False)
    return insertions_file


def sort_files(insertions_file):
    """Sort insertion files."""
    sorted_file = Path(insertions_file).with_suffix(".insertions.sorted")
    cmd = ["sort", "-k1,1", "-k2n,2", str(insertions_file)]
    try:
        with open(sorted_file, "w") as out_file:
            subprocess.run(cmd, stdout=out_file, check=True)
        logging.info(f"Sorted file: {insertions_file} -> {sorted_file}")
        return sorted_file
    except Exception as e:
        logging.error(f"Error sorting file {insertions_file}: {e}")
        raise


def merge_files(sorted_insertions_file):
    """Merge intervals mapping to same location using sorted files."""
    merged_file = sorted_insertions_file.with_suffix(".sorted.merged")
    cmd = ["bedtools", "merge", "-i",
           str(sorted_insertions_file), "-s", "-c", "1,6", "-o", "count,distinct"]
    try:
        with open(merged_file, "w") as out_file:
            subprocess.run(cmd, stdout=out_file, check=True)
        logging.info(f"Merged file: {sorted_insertions_file} -> {merged_file}")
        return merged_file
    except Exception as e:
        logging.error(f"Error merging file {sorted_insertions_file}: {e}")
        raise


def filter_insertions(merged_file):
    """Remove insertions supported by only 1 read."""
    filtered_file = merged_file.with_suffix(".merged.filtered")
    try:
        df = pd.read_table(merged_file, header=None)
        df[df[3] > 1].to_csv(filtered_file, index=False, header=False, sep='\t')
        logging.info(
            f"Filtered insertions: {merged_file} -> {filtered_file}")
        return filtered_file
    except Exception as e:
        logging.error(f"Error filtering insertions for {merged_file}: {e}")
        raise


def count_insertions_over_intervals(
    interval_files: List[Union[str, Path]],
    filtered_insertions_file: Union[str, Path]
) -> None:
    """
    Count insertions over specified genomic intervals using bedtools map.

    Parameters
    ----------
    interval_files : List[Union[str, Path]]
        List of paths to interval files (BED/GFF format)
    filtered_insertions_file : Union[str, Path]
        Path to the filtered insertions file

    Notes
    -----
    For each interval file, creates an output file with the format:
    {filtered_insertions_stem}_{interval_file_stem}.cnts
    containing the count of unique insertions and sum of insertions per interval.
    """
    # Convert to Path objects
    filtered_path = Path(filtered_insertions_file)
    # Validate filtered insertions file
    if not filtered_path.exists():
        raise FileNotFoundError(
            f"Filtered insertions file not found: {filtered_path}")
    if filtered_path.stat().st_size == 0:
        raise ValueError(f"Filtered insertions file is empty: {filtered_path}")

    # Process each interval file
    for interval_file in interval_files:
        interval_path = Path(interval_file)
        # Validate interval file
        if not interval_path.exists():
            raise FileNotFoundError(
                f"Interval file not found: {interval_path}")
        if interval_path.stat().st_size == 0:
            logging.warning(f"Interval file is empty: {interval_path}")
            continue

        # Prepare output file path
        file_name = interval_path.stem
        output_file = filtered_path.parent / \
            f"{filtered_path.name}_{file_name}.cnts"

        # Prepare bedtools command
        cmd = [
            "bedtools", "map",
            "-a", str(interval_path),
            "-b", str(filtered_path),
            "-c", "1,4",  # count column 1 and sum column 4
            "-o", "count,sum"
        ]
        logging.info(f"Processing interval file: {interval_path}")

        try:
            # Run bedtools command
            with open(output_file, "w") as out_file:
                process = subprocess.run(
                    cmd,
                    stdout=out_file,
                    stderr=subprocess.PIPE,
                    check=True,
                    text=True
                )

            # Check if output file was created and has content
            if not output_file.exists() or output_file.stat().st_size == 0:
                raise RuntimeError(
                    f"Output file is empty or was not created: {output_file}")
            logging.info(
                f"Successfully mapped intervals: {interval_path.name} with "
                f"{filtered_path.name} -> {output_file.name}"
            )

        except subprocess.CalledProcessError as e:
            error_msg = e.stderr.strip() if e.stderr else str(e)
            logging.error(
                f"bedtools mapping failed for {interval_path.name}: {error_msg}"
            )
            raise RuntimeError(f"bedtools mapping failed: {error_msg}") from e

        except Exception as e:
            logging.error(
                f"Error processing {interval_path.name}: {str(e)}"
            )
            raise


def map_sample(sample, bam_dir, output_dir, interval_files):
    """Process a sample"""
    try:
        setup_logging(output_dir)
        logging.info(f"Starting processing sample {sample}")
        logging.info("Merging bams")
        output_bam = merge_bams(sample, output_dir, bam_dir)
        logging.info("Filtering bams")
        bed_file = filter_bam(output_bam)
        logging.info("Processing orientation")
        insertions_file = process_orientation(bed_file)
        logging.info("Sorting")
        sorted_file = sort_files(insertions_file)
        logging.info("Merging sorted")
        merged_file = merge_files(sorted_file)
        logging.info("Filtering insertions")
        filtered_file = filter_insertions(merged_file)
        logging.info("Counting insertions")
        count_insertions_over_intervals(interval_files, filtered_file)
        logging.info(f"Sample {sample} processing completed successfully")
    except Exception as e:
        logging.error(f"Sample processing failed: {e}")
        raise

# This would be to merge output from multiple samples


def merge_cnts_files(
    interval_file: Union[str, Path],
    counts_dir: Union[str, Path],
    name: str = "",
    format: str = 'gff'
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Merge count files from multiple samples into consolidated transposon and read count matrices.

    Parameters
    ----------
    interval_file : str or Path
        Path to the interval file (GFF or BED format)
    counts_dir : str or Path
        Directory containing the count files to merge
    name : str, optional
        Prefix for output files, by default ""
    format : str, optional
        Format of the interval file ('gff' or 'bed'), by default 'gff'

    Returns
    -------
    Tuple[pd.DataFrame, pd.DataFrame]
        Two DataFrames containing merged transposon counts and read counts respectively

    """
    # Input validation
    if format not in ['gff', 'bed']:
        raise ValueError("Format must be either 'gff' or 'bed'")

    interval_path = Path(interval_file)
    counts_path = Path(counts_dir)

    if not interval_path.exists():
        raise FileNotFoundError(f"Interval file not found: {interval_path}")
    if not counts_path.exists():
        raise FileNotFoundError(f"Counts directory not found: {counts_path}")

    # Find count files
    file_name = interval_path.stem
    cnt_files = list(counts_path.rglob(f"*_{file_name}*.cnts"))

    if not cnt_files:
        raise ValueError(
            f"No count files found matching pattern '*_{file_name}*.cnts' in {counts_dir}")

    logging.info(f"Found {len(cnt_files)} count files to process")

    # Initialize lists for DataFrames
    tn_dfs: List[pd.DataFrame] = []
    read_dfs: List[pd.DataFrame] = []

    def process_id_column(df: pd.DataFrame, file_format: str) -> pd.Series:
        """Helper function to process ID column based on file format"""
        if file_format == 'bed':
            return df[[0, 1, 2]].astype(str).agg('_'.join, axis=1)
        else:  # gff format
            id_col = df[df.columns[8]]
            try:
                return (id_col.str.split("locus_tag=", expand=True)[1]
                        .str.split(";", expand=True)[0])
            except (IndexError, KeyError):
                raise ValueError(
                    "Invalid GFF format: missing or malformed locus_tag field")

    # Process each count file
    for file in cnt_files:
        try:
            sample_id = file.stem.split(".bed")[0].replace(".bam", '')
            logging.info(f"Processing sample: {sample_id}")
            # Read and validate input file
            ndf = pd.read_table(file, header=None)
            required_cols = [ndf.shape[1] - 2,
                             ndf.shape[1] - 1]  # Last two columns
            if any(col < 0 for col in required_cols):
                raise ValueError(
                    f"File {file} does not have the expected number of columns")

            # Rename columns
            ndf = ndf.rename(columns={
                ndf.columns[required_cols[0]]: "tn_per_gene",
                ndf.columns[required_cols[1]]: "read_per_gene"
            })

            # Process ID column
            ndf['ID'] = process_id_column(ndf, format)

            # Clean and transform data
            ndf = ndf[["ID", "tn_per_gene", "read_per_gene"]].set_index('ID')
            ndf["read_per_gene"] = pd.to_numeric(
                ndf["read_per_gene"].replace({".": "0"}),
                errors='coerce'
            ).fillna(0)

            # Append transformed data
            tn_dfs.append(ndf[["tn_per_gene"]].rename(
                columns={"tn_per_gene": sample_id}).T)
            read_dfs.append(
                ndf[["read_per_gene"]].rename(
                    columns={"read_per_gene": sample_id}).T
            )

        except Exception as e:
            logging.error(f"Error processing file {file}: {str(e)}")
            raise

    # Merge and convert data
    tn_df = pd.concat(tn_dfs).fillna(0).astype(int).sort_index()
    read_df = pd.concat(read_dfs).fillna(0).astype(int).sort_index()

    # Save output files
    today_str = datetime.today().strftime("%Y-%m-%d")
    output_prefix = f"{today_str}_{name}" if name else today_str
    print(output_prefix)
    output_tn = counts_path / f"{output_prefix}_transposon_counts.csv"
    output_reads = counts_path / f"{output_prefix}_read_counts.csv"

    tn_df.T.to_csv(output_tn)
    read_df.T.to_csv(output_reads)

    logging.info(f"Saved transposon counts to: {output_tn}")
    logging.info(f"Saved read counts to: {output_reads}")

    return tn_df, read_df
