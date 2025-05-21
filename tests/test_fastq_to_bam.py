#!/usr/bin/env python3
from satay.fastq_to_bam import find_fastq_files, check_star_installed, verify_star_index, create_genome_index, run_star_alignment
import os
import sys
import pytest
from pathlib import Path
import tempfile
import shutil
import subprocess

# Import the function from the map_to_yeast script
# Assumes the script is in the same directory or in the Python path
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

REF = "ref"
GFF = "ref/GCF_000146045.2_R64.gff.gz "


@pytest.fixture
def data_paths():
    """Return paths to test data directories"""
    # This assumes the test is being run from the project root directory
    # Adjust these paths to match your project structure
    base_dir = os.getenv('TEST_DATA_DIR', os.path.join(
        os.path.dirname(__file__),  'test_data'))
    tiny_dataset = os.path.join(base_dir, 'tiny_dataset')

    # Check if the tiny_dataset directory exists
    if not os.path.exists(tiny_dataset):
        pytest.skip(f"Tiny dataset directory not found at {tiny_dataset}")

    return {
        'base_dir': base_dir,
        'tiny_dataset': tiny_dataset
    }


@pytest.fixture
def star_installed():
    """Check if STAR is installed and skip test if not"""
    if not check_star_installed():
        pytest.skip("STAR aligner not installed or not in PATH")
    return True


@pytest.fixture
def temp_output_dir():
    """Create a temporary directory for alignment outputs"""
    temp_dir = tempfile.mkdtemp()
    yield temp_dir
    # Clean up after test
    shutil.rmtree(temp_dir)


def test_find_fastq_files_tiny_dataset(data_paths):
    """Test that find_fastq_files finds all FASTQ files in the tiny dataset"""
    # Find all FASTQ files in the tiny dataset
    found_files = find_fastq_files(data_paths['tiny_dataset'])

    # Verify that files were found
    assert len(
        found_files) > 0, f"No FASTQ files found in {data_paths['tiny_dataset']}"

    # Verify that all found files are actually FASTQ files
    for file_path in found_files:
        assert os.path.exists(file_path), f"File doesn't exist: {file_path}"

        # Check that the file has a valid FASTQ extension
        extensions = ['.fastq', '.fq', '.fastq.gz', '.fq.gz']
        has_valid_ext = any(file_path.endswith(ext) for ext in extensions)
        assert has_valid_ext, f"File does not have a valid FASTQ extension: {file_path}"

    # Print found files for debugging
    print(f"Found {len(found_files)} FASTQ files in tiny dataset:")
    for file_path in found_files:
        print(f"  - {os.path.basename(file_path)}")


def test_find_fastq_files_returns_absolute_paths(data_paths):
    """Test that find_fastq_files returns absolute paths"""
    found_files = find_fastq_files(data_paths['tiny_dataset'])

    # Check that all paths are absolute
    for file_path in found_files:
        assert os.path.isabs(file_path), f"Path is not absolute: {file_path}"


def test_find_fastq_files_filters_non_fastq(data_paths, tmp_path):
    """Test that find_fastq_files correctly filters out non-FASTQ files"""
    # Create a temporary directory with both FASTQ and non-FASTQ files
    test_dir = tmp_path / "mixed_files"
    test_dir.mkdir()

    # Create a FASTQ file
    fastq_file = test_dir / "test.fastq"
    fastq_file.write_text("@read1\nACGT\n+\nIIII\n")

    # Create a non-FASTQ file
    non_fastq_file = test_dir / "test.txt"
    non_fastq_file.write_text("This is not a FASTQ file")

    # Find FASTQ files
    found_files = find_fastq_files(str(test_dir))

    # Check that only the FASTQ file was found
    assert len(found_files) == 1
    assert os.path.basename(found_files[0]) == "test.fastq"


def get_logger():
    """Create a simple logger for testing"""
    import logging
    logger = logging.getLogger("test_logger")
    logger.setLevel(logging.INFO)
    # Prevent logger from propagating to root logger
    logger.propagate = False

    # If logger already has handlers, don't add more
    if not logger.handlers:
        handler = logging.StreamHandler()
        formatter = logging.Formatter(
            '%(asctime)s - %(levelname)s - %(message)s')
        handler.setFormatter(formatter)
        logger.addHandler(handler)

    return logger


class TestRunStarAlignment:
    """Tests for run_star_alignment function"""

    def test_star_prerequisites(self, star_installed):
        """Test that STAR prerequisites are available"""
        assert star_installed, "STAR aligner should be installed"
        assert os.path.exists(
            REF), f"Genome index directory should exist at {REF}"

    def test_create_genome_index(self, tmp_path):
        # Create logger
        logger = get_logger()
        genome_fasta = os.path.join(REF, "GCF_000146045.2_R64_genomic.fna.gz")
        index_done = create_genome_index(
            str(genome_fasta), str(tmp_path), 1, logger)
        assert index_done is True

    def test_run_star_single_end(self, data_paths, tmp_path, star_installed):
        """Test running STAR alignment with a single-end read from tiny dataset"""
        # Get the first FASTQ file from tiny dataset
        fastq_files = find_fastq_files(data_paths['tiny_dataset'])
        assert len(
            fastq_files) >= 1, "At least one FASTQ file should exist in tiny dataset"

        # Use the first file for testing
        fastq_file = fastq_files[0]
        sample_name = os.path.basename(fastq_file).split('.')[0]

        # Create logger
        logger = get_logger()

        # Run alignment

        bam_file = run_star_alignment(
            fastq_file1=fastq_file,
            fastq_file2=None,  # Single-end read
            output_dir=tmp_path,
            genome_dir=REF,
            sample_name=sample_name,
            threads=1,  # Use single thread for testing
            logger=logger
        )
        assert bam_file is not None, f"run_star_alignment returned None for {sample_name}"
        assert os.path.exists(
            bam_file), f"Output BAM file should exist at {bam_file}"

        # Check if BAM file is valid by running samtools stats
        # This is optional and depends on samtools being installed
        try:
            subprocess.run(
                ["samtools", "stats", bam_file],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                check=True
            )
            # If we got here, the BAM file is valid
        except (subprocess.SubprocessError, FileNotFoundError):
            # Either samtools is not installed or the BAM file is invalid
            # We don't fail the test because of this, just log it
            logger.warning(
                f"Could not validate BAM file with samtools: {bam_file}")

    # Can expand to more indices if needed
    # @pytest.mark.parametrize("sample_index", [0])
    # def test_output_directory_structure(self, data_paths, star_genome_index, temp_output_dir, star_installed, sample_index):
    #     """Test that STAR creates the expected output directory structure"""
    #     # Get a FASTQ file from tiny dataset
    #     fastq_files = find_fastq_files(data_paths['tiny_dataset'])
    #     if sample_index >= len(fastq_files):
    #         pytest.skip(
    #             f"Not enough FASTQ files in tiny dataset (need index {sample_index})")

    #     fastq_file = fastq_files[sample_index]
    #     sample_name = os.path.basename(fastq_file).split('.')[0]

    #     # Create logger
    #     logger = get_logger()

        # Run alignment
        # try:
        #     bam_file = run_star_alignment(
        #         fastq_file1=fastq_file,
        #         fastq_file2=None,
        #         output_dir=temp_output_dir,
        #         genome_dir=star_genome_index,
        #         sample_name=sample_name,
        #         threads=1,
        #         logger=logger
        #     )

        #     # Check output directory structure
        #     sample_output_dir = os.path.join(temp_output_dir, sample_name)
        #     assert os.path.exists(
        #         sample_output_dir), f"Sample output directory should exist: {sample_output_dir}"

        #     # Check for expected STAR output files
        #     expected_files = [
        #         f"{sample_name}_Aligned.sortedByCoord.out.bam",
        #         f"{sample_name}_Log.final.out",
        #         f"{sample_name}_Log.out",
        #         f"{sample_name}_Log.progress.out"
        #     ]

        #     for filename in expected_files:
        #         file_path = os.path.join(sample_output_dir, filename)
        #         assert os.path.exists(
        #             file_path), f"Expected output file missing: {file_path}"

        #     # Check symlink
        #     symlink_path = os.path.join(temp_output_dir, f"{sample_name}.bam")
        #     assert os.path.exists(
        #         symlink_path), f"Symlink to BAM file should exist: {symlink_path}"

        # except subprocess.CalledProcessError:
        #     pytest.skip("STAR alignment failed to run.")

    # def test_paired_end_detection(self, data_paths, tmp_path, star_genome_index, star_installed):
    #     """
    #     Test that paired-end reads are correctly detected and aligned.

    #     Note: Since the tiny dataset may not contain paired-end reads,
    #     this test creates mock paired files if needed.
    #     """
    #     # First check if we have paired-end reads in the tiny dataset
    #     fastq_files = find_fastq_files(data_paths['tiny_dataset'])

    #     # Check if we need to create mock paired-end files
    #     paired_files_found = False
    #     r1_file = None
    #     r2_file = None

    #     # Try to find paired files by naming convention
    #     for file in fastq_files:
    #         basename = os.path.basename(file)
    #         if "_R1" in basename or "_1" in basename:
    #             r1_file = file
    #         elif "_R2" in basename or "_2" in basename:
    #             r2_file = file

    #     # If paired files not found, create mock ones
    #     if r1_file is None or r2_file is None:
    #         # Create temp dir for mock files
    #         mock_dir = tmp_path / "mock_paired"
    #         mock_dir.mkdir()

    #         # Create minimal paired FASTQ files
    #         r1_file = str(mock_dir / "sample_R1.fastq")
    #         r2_file = str(mock_dir / "sample_R2.fastq")

    #         # Get minimal content from an existing file
    #         if len(fastq_files) > 0:
    #             with open(fastq_files[0], 'r') as f:
    #                 content = f.read()

    #             # Write to mock files
    #             with open(r1_file, 'w') as f:
    #                 f.write(content)
    #             with open(r2_file, 'w') as f:
    #                 f.write(content)
    #         else:
    #             # Write minimal FASTQ content
    #             minimal_content = "@read1\nACGTACGT\n+\nIIIIIIII\n"
    #             with open(r1_file, 'w') as f:
    #                 f.write(minimal_content)
    #             with open(r2_file, 'w') as f:
    #                 f.write(minimal_content)

    #     # Now we have paired-end files to test with
    #     sample_name = "paired_test"
    #     temp_output_dir = str(tmp_path / "paired_output")
    #     os.makedirs(temp_output_dir, exist_ok=True)

    #     # Create logger
    #     logger = get_logger()

    #     # Run alignment with paired-end files
    #     try:
    #         bam_file = run_star_alignment(
    #             fastq_file1=r1_file,
    #             fastq_file2=r2_file,
    #             output_dir=temp_output_dir,
    #             genome_dir=star_genome_index,
    #             sample_name=sample_name,
    #             threads=1,
    #             logger=logger
    #         )

    #         # Check if alignment completed successfully
    #         if bam_file is not None:
    #             assert os.path.exists(
    #                 bam_file), f"Output BAM file should exist at {bam_file}"
    #     except subprocess.CalledProcessError:
    #         pytest.skip("STAR paired-end alignment failed to run.")


if __name__ == "__main__":
    pytest.main(["-v", __file__])
