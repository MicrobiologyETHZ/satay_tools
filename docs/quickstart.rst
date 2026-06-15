Quick Start Guide
=================

This guide will help you get started with SATAY Tools for analyzing transposon insertion data.


1. **Prepare your data**: Ensure you have quality-controlled FASTQ files from your SATAY experiment


2. **Align reads**: Map FASTQ files to reference genome

.. code-block:: bash

   # Align FASTQ files to genome and generate BAM files
   satay align -f /path/to/fastq_dir -o /path/to/output_dir -g genome.fasta

  
3. **Map insertions**: Identify genomic location of transposon insertions and count reads supporting each insertion. Generates a file with transposon insertions and read counts per genome interval (i.e. CDS).

.. code-block:: bash

   # Call transposon insertions from BAM files
   satay map -b /path/to/bam_dir -o /path/to/output_dir -s sample_name -a annotations.gff


4. **Merge counts**: Combine transposon/read counts data from multiple samples.
   This writes ``{date}_{experiment_name}_transposon_counts.csv`` and
   ``{date}_{experiment_name}_read_counts.csv``.

.. code-block:: bash

   # Merge count files from multiple samples
   satay merge -d /path/to/counts_dir -a annotations.gff -n experiment_name


5. **Analyze**: Perform differential abundance analysis to identify significant changes in insertion frequency/ abundance between treatments. ``--counts-file`` is one of the count matrices from the merge step, and ``--sample_data`` is a CSV with a sample-ID column matching the matrix columns plus a condition column (see the :doc:`tutorials/index` for the format).

.. code-block:: bash

   # Perform differential analysis
   satay analyze -f {date}_experiment_name_transposon_counts.csv -s sample_data.csv -o /path/to/output_dir -c condition_column -b baseline_condition


