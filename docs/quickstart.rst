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


4. **Merge counts**: Combine transposon/read counts data from multiple samples

.. code-block:: bash

   # Merge count files from multiple samples
   satay merge -d /path/to/counts_dir -a annotations.gff -n experiment_name


5. **Analyze**: Perform differential abundance analysis to identify significant changes in insertion frequency/ abundance between treatments

.. code-block:: bash

   # Perform differential analysis
   satay analyze -f merged_counts.txt -s sample_data.txt -o /path/to/output_dir -c condition_column -b baseline_condition


