# SATAY Tools

A Python package for analyzing SATAY (Saturation Transposon Analysis in Yeast) transposon insertion sequencing data. This tool provides utilities for processing transposon mutagenesis datasets, including read mapping, insertion site identification, and downstream analysis of yeast fitness data.

## Features

- Process FASTQ files from transposon sequencing experiments
- Map transposon insertion sites to yeast genome
- Count and analyze insertion frequencies
- Perform differential analysis of insertion patterns

## Installation

```bash
   git clone https://github.com/MicrobiologyETHZ/satay_tools.git
   cd satay_tools
   conda env create -f satay_environment.yaml
   conda activate satay-tools
   pip install --no-deps .

```

## Usage

```bash

satay --help

```