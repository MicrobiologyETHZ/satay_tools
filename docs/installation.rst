Installation
============

Prerequisites
-------------

`git <https://git-scm.com/>`_ and a `conda <https://docs.conda.io/>`_
distribution (e.g. Miniconda or Miniforge) must be installed first.

.. code-block:: bash

   git clone https://github.com/MicrobiologyETHZ/satay_tools.git
   cd satay_tools
   conda env create -f satay_environment.yaml
   conda activate satay-tools
   pip install --no-deps .

   satay --help


.. note::

   **The** ``align`` **step requires Linux.** It runs the STAR aligner, which
   does not work on macOS (STAR fails to read the input reads). Run ``align`` on
   a Linux machine. The other steps (``map``, ``merge`` and ``analyze``) work on
   both macOS and Linux.
