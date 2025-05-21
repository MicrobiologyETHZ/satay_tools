#!/usr/bin/env python3
"""
Script to create tiered test datasets for transposon mutagenesis analysis:
1. Medium dataset: Extract reads mapping to chromosome IV from 6 BAM files
2. Small dataset: Subsample medium dataset to 5%
3. Tiny dataset: Subsample specified FASTQ files to 0.05%
"""

import os
import random
import argparse
import pysam
import gzip
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq


def extract_chrom_iv_reads(bam_file, output_dir, chrom_name="IV"):
    """
    Extract reads mapping to chromosome IV and output as FASTQ

    Args:
        bam_file (str): Path to input BAM file
        output_dir (str): Directory to save output FASTQ
        chrom_name (str): Chromosome name to extract (default: "IV")

    Returns:
        str: Path to output FASTQ file
    """
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)

    # Get base name for output file
    base_name = os.path.basename(bam_file).replace(".bam", "")
    output_file = os.path.join(
        output_dir, f"{base_name}_chr{chrom_name}.fastq.gz")

    print(f"Extracting chromosome {chrom_name} reads from {bam_file}...")

    # Open BAM file
    bamfile = pysam.AlignmentFile(bam_file, "rb")

    # Check if chromosome exists in the BAM file
    try:
        chrom_id = None
        for ref in bamfile.references:
            if chrom_name in ref:
                chrom_id = ref
                break

        if chrom_id is None:
            print(
                f"Warning: Chromosome {chrom_name} not found in {bam_file}. Checking for numeric format...")
            # Try numeric format
            for ref in bamfile.references:
                if ref.endswith(chrom_name) or ref == chrom_name:
                    chrom_id = ref
                    break

        if chrom_id is None:
            print(
                f"Error: Chromosome {chrom_name} not found in {bam_file}. Available references: {bamfile.references}")
            return None

    except Exception as e:
        print(f"Error checking chromosome: {e}")
        return None

    # Open output file
    with gzip.open(output_file, 'wt') as out_handle:
        # Iterate through reads mapped to chromosome IV
        count = 0
        for read in bamfile.fetch(chrom_id):
            if read.is_unmapped:
                continue

            # Extract read sequence and quality
            seq = read.query_sequence
            if seq is None:  # Skip reads without sequence
                continue

            qual = ''.join([chr(q+33) for q in read.query_qualities]
                           ) if read.query_qualities else '*' * len(seq)

            # Format as FASTQ and write
            header = f"@{read.query_name}"
            if read.is_read1:
                header += "/1"
            elif read.is_read2:
                header += "/2"

            out_handle.write(f"{header}\n{seq}\n+\n{qual}\n")
            count += 1

            # Print progress
            if count % 100000 == 0:
                print(f"Processed {count} reads...")

    print(
        f"Extracted {count} reads from chromosome {chrom_name} to {output_file}")
    return output_file


def subsample_fastq(fastq_file, output_dir, sample_rate, seed=42):
    """
    Subsample a FASTQ file

    Args:
        fastq_file (str): Path to input FASTQ file (can be gzipped)
        output_dir (str): Directory to save output FASTQ
        sample_rate (float): Proportion of reads to keep (0.0-1.0)
        seed (int): Random seed for reproducibility

    Returns:
        str: Path to output FASTQ file
    """
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)

    # Get base name for output file
    base_name = os.path.basename(fastq_file).replace(
        ".fastq.gz", "").replace(".fq.gz", "")
    output_file = os.path.join(
        output_dir, f"{base_name}_subsample_{int(sample_rate*100)}pct.fastq.gz")

    print(f"Subsampling {fastq_file} at {sample_rate*100}% rate...")

    # Set random seed for reproducibility
    random.seed(seed)

    # Open input file
    if fastq_file.endswith(".gz"):
        in_handle = gzip.open(fastq_file, "rt")
    else:
        in_handle = open(fastq_file, "r")

    # Open output file
    with gzip.open(output_file, 'wt') as out_handle:
        # Read through FASTQ file using SeqIO
        count = 0
        kept = 0

        for record in SeqIO.parse(in_handle, "fastq"):
            count += 1

            # Decide whether to keep this read
            if random.random() < sample_rate:
                # Write to output file
                SeqIO.write(record, out_handle, "fastq")
                kept += 1

            # Print progress
            if count % 100000 == 0:
                print(f"Processed {count} reads, kept {kept}...")

    in_handle.close()
    print(
        f"Subsampled {kept}/{count} reads ({kept/count*100:.2f}%) to {output_file}")
    return output_file


def create_tiered_datasets(bam_files, output_base_dir, tiny_fastq_files=None, chrom_name="IV"):
    """
    Create tiered test datasets

    Args:
        bam_files (list): List of paths to input BAM files
        output_base_dir (str): Base directory for all output files
        tiny_fastq_files (list, optional): List of paths to FASTQ files for tiny dataset
        chrom_name (str): Chromosome name to extract
    """
    # Create directories for each tier
    medium_dir = os.path.join(output_base_dir, "medium_dataset")
    small_dir = os.path.join(output_base_dir, "small_dataset")
    tiny_dir = os.path.join(output_base_dir, "tiny_dataset")

    os.makedirs(medium_dir, exist_ok=True)
    os.makedirs(small_dir, exist_ok=True)
    os.makedirs(tiny_dir, exist_ok=True)

    # Step 1: Create medium dataset - extract chromosome IV reads from all BAM files
    medium_fastqs = []
    for bam_file in bam_files:
        fastq_file = extract_chrom_iv_reads(bam_file, medium_dir, chrom_name)
        if fastq_file:
            medium_fastqs.append(fastq_file)

    # Step 2: Create small dataset - subsample medium dataset to 5%
    small_fastqs = []
    for fastq_file in medium_fastqs:
        small_fastq = subsample_fastq(fastq_file, small_dir, 0.05)
        small_fastqs.append(small_fastq)

    # Step 3: Create tiny dataset - subsample provided FASTQ files to 0.05%
    if tiny_fastq_files:
        print(
            f"Creating tiny dataset from {len(tiny_fastq_files)} provided FASTQ files")
        for fastq_file in tiny_fastq_files:
            if os.path.exists(fastq_file):
                # Direct 0.05% sampling
                subsample_fastq(fastq_file, tiny_dir, 0.0005)
            else:
                print(f"Warning: FASTQ file not found: {fastq_file}")
    else:
        print("No FASTQ files provided for tiny dataset. Skipping tiny dataset creation.")


def main():
    parser = argparse.ArgumentParser(
        description='Create tiered test datasets for transposon mutagenesis analysis')
    parser.add_argument('--bam_files', nargs='+', required=True,
                        help='Space-separated list of input BAM files')
    parser.add_argument('--output_dir', required=True,
                        help='Base directory for output datasets')
    parser.add_argument('--chromosome', default="NC_001136.10",
                        help='Chromosome name to extract (default: NC_001136.10 (IV))')
    parser.add_argument('--tiny_fastq_files', nargs='+',
                        help='FASTQ files to use for tiny dataset generation')

    args = parser.parse_args()

    create_tiered_datasets(args.bam_files, args.output_dir,
                           args.tiny_fastq_files, args.chromosome)

    print("Tiered datasets creation complete.")
    print(
        f"- Medium dataset: {os.path.join(args.output_dir, 'medium_dataset')}")
    print(f"- Small dataset: {os.path.join(args.output_dir, 'small_dataset')}")
    print(f"- Tiny dataset: {os.path.join(args.output_dir, 'tiny_dataset')}")


if __name__ == "__main__":
    main()
