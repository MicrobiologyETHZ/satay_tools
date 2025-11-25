Tutorial
========

This tutorial demonstrates how to process SATAY (Saturation Transposon Analysis in Yeast) 
data using the test dataset included with satay-tools.


Setting Up Your Workspace
------------

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

Parameters
^^^^^^^^^^

* ``-f, --fastq-dir``: Directory containing FASTQ files (can be gzipped)
* ``-o, --output-dir``: Output directory for BAM files
* ``-g, --genome``: Reference genome in FASTA format (GCF_000146045.2 is included with satay-tools)

Optional parameters:

* ``--threads``: Number of threads for alignment (default: 4)

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
* ``-s, --sample``: Sample identifier to process. Must be a part of `bam` file name
* ``-a, --annotation``: Gene annotation file in GFF format

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
     -d test_out/ \
     -a ref/GCF_000146045.2.genes.gff.gz \
     -n test1

Parameters
^^^^^^^^^^

* ``-d, --data-dir``: Directory containing count files from map step
* ``-a, --annotation``: Gene annotation file that was used for counting
* ``-n, --name``: Name prefix for output files

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
     -o test_out \
     -c conc \
     -b "0"

Parameters
^^^^^^^^^^

* ``-f, --counts-file``: Merged count matrix from merge step
* ``-s, --sample-info``: Sample metadata file (CSV format)
* ``-o, --output-dir``: Output directory for analysis results
* ``-c, --condition``: Column name in metadata for experimental condition
* ``-b, --baseline``: Baseline/reference condition for comparison

Sample Metadata Format
^^^^^^^^^^^^^^^^^^^^^^^

The metadata file should be a CSV with at minimum 2 columns:

Where:

* ``sample`` column: Sample identifiers matching column names in count matrix
* ``treatment`` column: Experimental condition (or any grouping variable)

Outputs
^^^^^^^

The analyze step produces:

* Differential abundance results table

Key columns in the results table:

* ``log2FoldChange``: Effect size (positive = enriched, negative = depleted)
* ``padj``: Adjusted p-value (FDR-corrected)
* ``baseMean``: Average normalized count across samples
