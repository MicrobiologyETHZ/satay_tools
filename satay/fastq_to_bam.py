#!/usr/bin/env python3
"""
Script to map FASTQ files to the yeast genome using STAR aligner.
Takes a directory of FASTQ files and aligns them to a specified yeast genome reference.

"""

import os
import subprocess
import logging
from datetime import datetime
import sys
import re
import gzip
import sys
import tempfile
from pathlib import Path
import shutil


def setup_logger(log_dir):
    """Set up the logger to write to both console and file"""
    if not os.path.exists(log_dir):
        os.makedirs(log_dir)

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    log_file = os.path.join(log_dir, f"yeast_mapping_{timestamp}.log")

    # Create logger
    logger = logging.getLogger("yeast_mapper")
    logger.setLevel(logging.INFO)

    # Create handlers
    file_handler = logging.FileHandler(log_file)
    console_handler = logging.StreamHandler(sys.stdout)

    # Create formatter
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    file_handler.setFormatter(formatter)
    console_handler.setFormatter(formatter)

    # Add handlers to logger
    logger.addHandler(file_handler)
    logger.addHandler(console_handler)

    return logger


def find_fastq_files(fastq_dir, suffix=None):
    """Find all FASTQ files in the given directory"""
    fastq_files = []
    extensions = [suffix] if suffix else [
        ".fastq", ".fq", ".fastq.gz", ".fq.gz"]
    for ext in extensions:
        fastq_files.extend(Path(fastq_dir).rglob(f"*{ext}"))
    return sorted(fastq_files)


def get_sample_name(fastq_file):
    """Extract sample name from FASTQ file path"""
    base_name = os.path.basename(fastq_file)

    # Remove extensions
    sample_name = re.sub(r'\.fastq(\.gz)?$|\.fq(\.gz)?$', '', base_name)

    # Remove common sequencing-related suffixes
    sample_name = re.sub(r'_R[12]$|_R[12]_001$|_[12]$', '', sample_name)

    return sample_name


def find_paired_files(fastq_files):
    """Group FASTQ files into pairs if they appear to be paired-end reads"""
    paired_files = {}

    for fastq_file in fastq_files:
        base_name = os.path.basename(fastq_file)

        # Look for common paired-end indicators
        r1_match = re.search(r'_(R1|1)(_001)?\.f(ast)?q(\.gz)?$', base_name)
        r2_match = re.search(r'_(R2|2)(_001)?\.f(ast)?q(\.gz)?$', base_name)

        if r1_match:
            sample_name = get_sample_name(fastq_file)
            if sample_name not in paired_files:
                paired_files[sample_name] = {'R1': None, 'R2': None}
            paired_files[sample_name]['R1'] = fastq_file
        elif r2_match:
            sample_name = get_sample_name(fastq_file)
            if sample_name not in paired_files:
                paired_files[sample_name] = {'R1': None, 'R2': None}
            paired_files[sample_name]['R2'] = fastq_file
        else:
            # Treat as single-end
            sample_name = get_sample_name(fastq_file)
            if sample_name not in paired_files:
                paired_files[sample_name] = {'R1': fastq_file, 'R2': None}

    return paired_files


def run_star_alignment(fastq_file1, fastq_file2, output_dir, genome_dir, sample_name, threads, logger,
                       limit_bam_sort_ram=2_000_000_000):
    """Run STAR alignment for single or paired-end reads"""

    # Coerce path-like args to str so they can be passed to STAR and joined for logging
    fastq_file1 = str(fastq_file1)
    if fastq_file2 is not None:
        fastq_file2 = str(fastq_file2)

    # Create sample output directory
    sample_output_dir = os.path.join(output_dir, sample_name)
    if not os.path.exists(sample_output_dir):
        os.makedirs(sample_output_dir)

    # Construct STAR command
    star_cmd = [
        "STAR",
        "--runThreadN", str(threads),
        "--genomeDir", genome_dir,
        "--outSAMtype", "BAM", "SortedByCoordinate",
        "--outFileNamePrefix", os.path.join(sample_output_dir,
                                            f"{sample_name}_"),
        "--limitBAMsortRAM", str(limit_bam_sort_ram),
        "--readFilesIn", fastq_file1
    ]

    # Add second read file if paired-end
    if fastq_file2:
        star_cmd.append(fastq_file2)

    # If input is gzipped, add appropriate parameter
    if str(fastq_file1).endswith('.gz'):
        star_cmd.extend(["--readFilesCommand", "zcat"])

    # Log the command
    logger.info(
        f"Running STAR alignment for {sample_name}: {' '.join(star_cmd)}")

    # Run STAR
    try:
        subprocess.run(star_cmd, check=True)
        logger.info(f"STAR alignment completed successfully for {sample_name}")

        # Create a convenience symlink to the final BAM (overwrite if re-running)
        final_bam = os.path.join(
            sample_output_dir, f"{sample_name}_Aligned.sortedByCoord.out.bam")
        renamed_bam = os.path.join(output_dir, f"{sample_name}.bam")
        if os.path.lexists(renamed_bam):
            logger.info(f"Replacing existing symlink/file: {renamed_bam}")
            os.remove(renamed_bam)
        os.symlink(final_bam, renamed_bam)
        logger.info(f"Created symlink: {renamed_bam}")

        return renamed_bam
    except subprocess.CalledProcessError as e:
        logger.error(f"STAR alignment failed for {sample_name}: {e}")
        return None


def check_star_installed():
    """Check if STAR is installed and available"""
    try:
        subprocess.run(["STAR", "--version"],
                       stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        return True
    except FileNotFoundError:
        return False


def create_genome_index(genome_fasta, genome_dir, threads, logger):
    logger.info(f"Creating STAR genome index in {genome_dir}")

    # Coerce path-like args to str for string ops and subprocess/logging
    genome_fasta = str(genome_fasta)
    genome_dir = str(genome_dir)

    # Create genome directory if it doesn't exist
    if not os.path.exists(genome_dir):
        logger.error("Genome directory does not exist")
        sys.exit(1)

    is_gzipped = genome_fasta.endswith('.gz')
    tmp_fasta = None
    star_success = False

    try:
        if is_gzipped:
            logger.info(f"Decompressing {genome_fasta} to a temporary file")
            fd, tmp_fasta = tempfile.mkstemp(suffix=".fa")
            with gzip.open(genome_fasta, 'rb') as f_in, os.fdopen(fd, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
            fasta_for_star = tmp_fasta
        else:
            fasta_for_star = genome_fasta

        # For yeast genome, these parameters work well
        # Adjusted for smaller genome size
        star_cmd = [
            "STAR",
            "--runMode", "genomeGenerate",
            "--runThreadN", str(threads),
            "--genomeDir", genome_dir,
            "--genomeFastaFiles", fasta_for_star,  # Use the appropriate file
            # Reduced for smaller genome (default: 14)
            "--genomeSAindexNbases", "10",
            "--limitGenomeGenerateRAM", "16000000000"  # 16GB limit
        ]

        cmd_str = ' '.join(star_cmd)
        logger.info(f"Running STAR genome generation: {cmd_str}")

        # Run STAR command and capture output for logging
        process = subprocess.run(
            star_cmd,
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
        )

        # Log some of the output for debugging purposes
        logger.info("STAR completed with output:")
        # Log last 10 lines of stdout
        for line in process.stdout.splitlines()[-10:]:
            logger.info(f"STAR stdout: {line}")

        # Check specifically for successful index creation
        if not verify_star_index(genome_dir, logger):
            logger.error("STAR completed but did not create a valid index")
            star_success = False
        else:
            logger.info("STAR genome index creation completed successfully")
            star_success = True

    except subprocess.CalledProcessError as e:
        logger.error(
            f"STAR genome index creation failed with code {e.returncode}")
        if e.stdout:
            logger.error(f"STAR stdout: {e.stdout}")
        if e.stderr:
            logger.error(f"STAR stderr: {e.stderr}")
        star_success = False

    except Exception as e:
        logger.error(f"Error during genome index creation: {e}")
        star_success = False

    finally:
        if tmp_fasta and os.path.exists(tmp_fasta):
            try:
                os.remove(tmp_fasta)
                logger.info(f"Removed temporary file {tmp_fasta}")
            except Exception as cleanup_error:
                logger.error(
                    f"Failed to remove temporary file: {cleanup_error}")

        return star_success


def verify_star_index(genome_dir, logger):
    """
    Verify that STAR actually created a valid index.

    Args:
        genome_dir (str): Directory containing the genome index
        logger: Logger object for logging messages

    Returns:
        bool: True if valid index exists, False otherwise
    """
    import os

    required_files = [
        "Genome",
        "SA",
        "SAindex",
        "chrLength.txt",
        "chrNameLength.txt",
        "chrName.txt",
        "chrStart.txt",
        "genomeParameters.txt"
    ]

    missing_files = []
    for file in required_files:
        file_path = os.path.join(genome_dir, file)
        if not os.path.exists(file_path):
            missing_files.append(file)

    if missing_files:
        logger.info(
            f"Missing required index files: {', '.join(missing_files)}")
        return False

    logger.info("Verified STAR index was created successfully")
    return True


def map_fastq_to_bam(fastq_dir, output_dir, genome_fasta, threads,
                     suffix=[], single_end=True,
                     limit_bam_sort_ram=2_000_000_000):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Set up logger
    logger = setup_logger(output_dir)

    # Check if STAR is installed
    if not check_star_installed():
        logger.error(
            "STAR aligner not found. Please install STAR and add it to your PATH.")
        sys.exit(1)

    # Verify genome index
    genome_dir = Path(genome_fasta).parent
    if not verify_star_index(genome_dir, logger):
        if Path(genome_fasta).is_file():
            logger.info(
                f"STAR genome index not found or incomplete in {genome_dir}")
            logger.info(
                f"Will create index using FASTA file: {genome_fasta}")
            index_done = create_genome_index(
                genome_fasta, genome_dir, threads, logger)
            if not index_done:
                logger.error("Failed to create genome index. Exiting.")
                sys.exit(1)
        else:
            logger.error(
                f"Genome FASTA file not found: {genome_fasta}")
            sys.exit(1)
    logger.info(f"START genome index found in {genome_dir}")
    # Find FASTQ files
    logger.info(f"Searching for FASTQ files in {fastq_dir} with rglob")
    fastq_files = find_fastq_files(fastq_dir, suffix=suffix)

    if not fastq_files:
        logger.error(f"No FASTQ files found in {fastq_dir}")
        sys.exit(1)

    logger.info(f"Found {len(fastq_files)} FASTQ files")

    # Process FASTQ files based on pairing mode
    if single_end:
        # Force single-end mode
        logger.info("Processing all files as single-end reads")
        for fastq_file in fastq_files:
            sample_name = get_sample_name(fastq_file)
            run_star_alignment(str(fastq_file), None, str(output_dir),
                               str(genome_dir), sample_name, threads, logger,
                               limit_bam_sort_ram=limit_bam_sort_ram)
    else:
        logger.error('Paired end not implemented')
        # Auto-detect or use paired mode
        # paired_files = find_paired_files(fastq_files)
        # logger.info(f"Grouped files into {len(paired_files)} sample(s)")

        # for sample_name, files in paired_files.items():
        #     if args.paired and (files['R1'] is None or files['R2'] is None):
        #         logger.warning(
        #             f"Skipping {sample_name} - missing R1 or R2 file for paired-end mode")
        #         continue

        #     if files['R2'] is not None:
        #         logger.info(f"Processing {sample_name} as paired-end")
        #         run_star_alignment(files['R1'], files['R2'], args.output_dir,
        #                            args.genome_dir, sample_name, args.threads, logger)
        #     else:
        #         logger.info(f"Processing {sample_name} as single-end")
        #         run_star_alignment(files['R1'], None, args.output_dir,
        #                            args.genome_dir, sample_name, args.threads, logger)

    logger.info("STAR alignment process completed.")
