Tutorial
========

This tutorial demonstrates how to process SATAY (Saturation Transposon Analysis in Yeast) 
data using the test dataset included with satay-tools.


Setting Up Your Workspace
-------------------------

For detailed installation instructions, see :doc:`../installation`.

Before running the tutorial, make sure your activate `satay-tools` environment, `cd` into satay_tools, and create an output directory for your results:

.. code-block:: bash

   conda activate satay-tools
   cd satay_tools
   output_dir=~/test_out # Change this to desired location
   mkdir -p $output_dir


Pipeline Overview
-----------------

The SATAY analysis pipeline consists of four main steps:

1. **align**: Map reads to reference genome using STAR aligner
2. **map**: Identify transposon insertion sites from aligned reads
3. **merge**: Combine insertion counts across samples into count matrices
4. **analyze**: Perform statistical analysis to identify fitness-altering mutations

Step 1: Align Reads
-------------------

Align SATAY sequencing reads to the reference genome:

.. code-block:: bash

   satay align \
     -f tests/test_data/medium_dataset/ \
     -o $output_dir \
     -g ref/GCF_000146045.2_R64_genomic.fna.gz

.. note::

   On **macOS**, aligning gzipped FASTQ files fails because of a STAR
   limitation (see :doc:`../installation`). Decompress the test FASTQ files
   first (``gunzip tests/test_data/medium_dataset/*.fastq.gz``) and point
   ``-f`` at the uncompressed files, or run this step on Linux.

Parameters
^^^^^^^^^^

* ``-f, --fastq-dir``: Directory containing FASTQ files (can be gzipped). Reads are processed as single-end.
* ``-o, --output-dir``: Output directory for BAM files
* ``-g, --genome-fasta``: Reference genome in FASTA format (GCF_000146045.2 is included with satay-tools)

Optional parameters:

* ``-t, --threads``: Number of threads for alignment (default: 4)
* ``-r, --limit-bam-sort-ram``: Max RAM in bytes for STAR BAM sorting (default: 2000000000). Increase if STAR reports a BAM-sorting RAM error.

.. note::

   Only single-end reads are supported. On the first run, if a STAR genome
   index is not already present next to the FASTA, satay-tools builds one
   automatically (tuned for the small yeast genome, with a 16 GB index-build
   RAM cap). The index is written alongside the reference and reused on
   subsequent runs.

Outputs
^^^^^^^

The align step produces:

* ``{sample}.bam``: Merged and sorted BAM files for each sample


Step 2: Map Insertion Sites
-----------------------------

Identify transposon insertion sites from aligned BAM files:

.. code-block:: bash

   satay map \
     -b $output_dir \
     -o $output_dir \
     -s 20190221.A-2_noaF \
     -a ref/GCF_000146045.2.genes.gff.gz

Parameters
^^^^^^^^^^

* ``-b, --bam-dir``: Directory containing BAM files from align step
* ``-o, --output-dir``: Output directory for insertion site files
* ``-s, --sample-name``: Sample identifier to process. Must be a part of the BAM file name
* ``-a, --gff``: Annotation file(s) in GFF or BED format. May be given multiple times to count over several interval sets.

Outputs
^^^^^^^

The map step produces:


* ``{sample}_*.cnts``: Per-gene insertion and read counts
* ``{sample}.bed``: Filtered alignments in BED format
* ``{sample}.bed.insertions.sorted.merged.filtered``: High-confidence insertion sites
* ``process_samples.log``: Processing log file

Step 3: Merge Counts
--------------------

Combine insertion site data from multiple samples into one count matrix:

.. code-block:: bash

   satay merge \
     -d $output_dir \
     -a ref/GCF_000146045.2.genes.gff.gz \
     -n test1

Parameters
^^^^^^^^^^

* ``-d, --counts-dir``: Directory containing count files from map step
* ``-a, --gff``: Name of the annotation file that was used for counting
* ``-n, --name``: Name prefix for output files

Optional parameters:

* ``--format``: Format of the annotation file, ``gff`` or ``bed`` (default: ``gff``)

Outputs
^^^^^^^

The merge step produces two count matrices:

* ``{date}_test1_transposon_counts.csv``: Number of unique transposon insertions per gene per sample
* ``{date}_test1_read_counts.csv``: Total read depth per gene per sample

These matrices have genes as rows and samples as columns, suitable for downstream statistical analysis.

Step 4: Statistical Analysis
-----------------------------

Perform differential abundance analysis to identify genes with altered fitness between conditions:

.. code-block:: bash

   satay analyze \
     -f tests/test_data/test_merged_counts.txt \
     -s tests/test_data/test-metadata.csv \
     -o $output_dir \
     -c conc \
     -b "0"

Parameters
^^^^^^^^^^

* ``-f, --counts-file``: Merged count matrix from merge step
* ``-s, --sample_data``: Sample metadata file (CSV format)
* ``-o, --output-dir``: Output directory for analysis results
* ``-c, --comp-col``: Column name in metadata that defines the contrasts (experimental condition)
* ``-b, --baseline``: Baseline/reference value in ``--comp-col``; all other values are compared to it

Optional parameters:

* ``-a, --gff``: GFF file used to annotate results with gene names (optional)
* ``-l, --filter``: Drop genes/intervals with fewer than this many counts across samples (default: 100)
* ``--alpha``: FDR cutoff for DESeq2 (default: 0.05)
* ``--sample-id-col``: Column in the metadata holding sample IDs; must match the count matrix columns (default: ``sample_id``)
* ``--ids``: GFF fields to keep when annotating (default: ``locus_tag gene``)
* ``-t, --threads``: Number of CPUs for DESeq2 inference (default: min(8, available cores))

Sample Metadata Format
^^^^^^^^^^^^^^^^^^^^^^^

The metadata file is a CSV with one row per sample. It must contain a sample-ID
column (named ``sample_id`` by default, configurable with ``--sample-id-col``)
whose values match the column names of the count matrix, plus the column passed
to ``--comp-col`` holding the condition/grouping variable. For example:

.. code-block:: text

   sample_id,conc
   20190221.A-1_noaF,0
   20190221.A-1_4nMaF,4
   20190221.A-2_noaF,0
   20190221.A-2_4nMaF,4

Here ``--comp-col conc`` defines the contrasts and ``--baseline 0`` sets the
reference level, so each non-zero concentration is compared against ``0``.

Outputs
^^^^^^^

The analyze step produces:

* Differential abundance results table

Key columns in the results table:

* ``log2FoldChange``: Effect size (positive = enriched, negative = depleted)
* ``padj``: Adjusted p-value (FDR-corrected)
* ``baseMean``: Average normalized count across samples
