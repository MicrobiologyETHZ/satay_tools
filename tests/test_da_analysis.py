#!/usr/bin/env python3
import os
from pathlib import Path
import pytest

pytest.importorskip("pydeseq2")
from satay.da_analysis import gene_da


@pytest.fixture
def data_paths():
    base_dir = os.getenv('TEST_DATA_DIR', os.path.join(
        os.path.dirname(__file__), 'test_data'))
    counts = Path(base_dir) / "test_merged_counts.txt"
    metadata = Path(base_dir) / "test-metadata.csv"
    if not counts.exists() or not metadata.exists():
        pytest.skip("DESeq test fixtures not found")
    return {"counts": counts, "metadata": metadata}


def test_gene_da_smoke(data_paths, tmp_path):
    """Run the full DESeq differential-analysis pipeline end-to-end and check outputs."""
    gene_da(
        counts_file=data_paths["counts"],
        sample_data_file=data_paths["metadata"],
        output_dir=tmp_path,
        filter=100,
        comp_col="conc",
        baseline="0",
        a=0.05,
        n_cpus=1,
    )

    results = list(tmp_path.glob("*_conc_*.csv"))
    vst = list(tmp_path.glob("*_vstcounts.csv"))
    assert results, "DESeq results CSV was not created"
    assert vst, "VST counts CSV was not created"

    import pandas as pd
    res = pd.read_csv(results[0])
    assert len(res) > 0, "DESeq results table is empty"
    for col in ("log2FoldChange", "padj", "contrast"):
        assert col in res.columns, f"missing expected column '{col}' in results"
