#!/usr/bin/env python3
from satay.process_satay_bams import merge_bams, filter_bam
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


@pytest.fixture
def data_paths():
    """Return paths to test data directories"""
    # This assumes the test is being run from the project root directory
    # Adjust these paths to match your project structure
    base_dir = os.getenv('TEST_DATA_DIR', os.path.join(
        os.path.dirname(__file__),  'test_data'))
    small_dataset = os.path.join(base_dir, 'small_dataset')
    filter_control_bed = os.path.join(small_dataset, '20190221.A-1_noaF.bed')
    # Check if the tiny_dataset directory exists
    if not os.path.exists(small_dataset):
        pytest.skip(f"Small dataset directory not found at {small_dataset}")

    return {
        'base_dir': base_dir,
        'small_dataset': small_dataset,
        'filter_controld_bed': filter_control_bed
    }


@pytest.fixture
def control_bed_file(test_data_path):
    """Return the path to the control BED file"""
    control_bed = test_data_path / "control.bed"
    assert control_bed.exists(
    ), f"Control BED file {control_bed} does not exist"
    return control_bed


def compare_bed_files(file1, file2):
    """Compare two BED files, returning True if they have the same content"""
    with open(file1, 'r') as f1, open(file2, 'r') as f2:
        # Read and sort the lines from both files (to handle different ordering)
        lines1 = sorted(f1.readlines())
        lines2 = sorted(f2.readlines())

        # Compare lines
        return lines1 == lines2


def test_merge_bams_integration(data_paths, tmp_path):
    """Integration test that actually runs samtools merge (optional)"""
    # Skip if samtools is not installed
    try:
        subprocess.run(["samtools", "--version"],
                       capture_output=True, check=True)
    except (subprocess.SubprocessError, FileNotFoundError):
        pytest.skip("samtools not available")

    # Create output directory
    output_dir = tmp_path / "merged"
    output_dir.mkdir()

    # Sample name to test
    sample = "20190221.A-1_noaF"  # Adjust to match your actual sample name
    output_bam = merge_bams(sample, str(output_dir),
                            str(data_paths['small_dataset']))

    # Verify the BAM file was created
    print(output_bam)
    assert output_bam.exists()

    # Optional: Verify the BAM file is valid by running samtools view
    try:
        result = subprocess.run(
            ["samtools", "view", "-H", str(output_bam)],
            capture_output=True,
            check=True,
            text=True
        )
        # Check that we got a header
        assert "@HD" in result.stdout
    except subprocess.SubprocessError:
        pytest.fail("Generated BAM file appears to be invalid")


def test_filter_bam_integration(data_paths, tmp_path):
    """Integration test that runs the actual filter_bam function and compares with control"""
    # Skip if samtools or bedtools are not installed
    try:
        subprocess.run(["samtools", "--version"],
                       capture_output=True, check=True)
        subprocess.run(["bedtools", "--version"],
                       capture_output=True, check=True)
    except (subprocess.SubprocessError, FileNotFoundError):
        pytest.skip("samtools or bedtools not available")

    # Find a test BAM file
    test_bam = Path(data_paths['small_dataset'])/'20190221.A-1_noaF.bam'

    # Copy the test BAM file to the temp directory (to avoid modifying original)

    temp_bam = tmp_path / test_bam.name
    shutil.copy(test_bam, temp_bam)

    # Run the filter_bam function
    result_bed = filter_bam(temp_bam)
    # Verify the BED file was created
    assert result_bed.exists(), f"Output BED file {result_bed} was not created"
    control_bed_file = data_paths['filter_controld_bed']
    # Compare with control BED file
    assert compare_bed_files(result_bed, control_bed_file), \
        f"Output BED file {result_bed} does not match control {control_bed_file}"
