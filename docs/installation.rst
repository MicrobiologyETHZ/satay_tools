Installation
============


.. code-block:: bash

   git clone https://github.com/MicrobiologyETHZ/satay_tools.git
   cd satay_tools
   conda env create -f satay_environment.yaml
   conda activate satay-tools
   pip install --no-deps .

   satay --help


.. note::

   **macOS users:** the ``satay align`` step runs the STAR aligner with
   ``--readFilesCommand`` to read **gzipped** FASTQ files on the fly. STAR's
   macOS (conda) builds cannot spawn this decompression subprocess, so aligning
   gzipped FASTQ files fails on macOS with
   ``EXITING: ... Failed spawning readFilesCommand``. Work around this by either:

   * decompressing the FASTQ files first (e.g. ``gunzip *.fastq.gz``) and
     passing the uncompressed files to ``satay align``, or
   * running the ``align`` step on Linux.

   Uncompressed input works on macOS, and the other pipeline steps
   (``map``, ``merge``, ``analyze``) are unaffected.
