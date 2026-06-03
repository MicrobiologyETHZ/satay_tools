import click
import os
from pathlib import Path
from typing import List, Optional, Tuple, Union
from satay.fastq_to_bam import map_fastq_to_bam
from satay.process_satay_bams import map_sample, merge_cnts_files
from satay.da_analysis import gene_da


@click.group()
def satay():
    """SATAY: A tool for analyzing transposon insertion sequencing data for yeast."""
    pass


@satay.command()
@click.option('--fastq-dir', '-f',
              type=click.Path(exists=True, dir_okay=True, path_type=Path),
              required=True,
              help='Directory containing FASTQ files')
@click.option('--output-dir', '-o',
              type=click.Path(exists=True, dir_okay=True, path_type=Path),
              required=True,
              help='Output directory')
@click.option('--genome-fasta', '-g',
              type=click.Path(exists=True, path_type=Path),
              required=True,
              help='Genome fasta')
@click.option('--threads', '-t',
              type=int,
              default=4,
              show_default=True,
              help='Number of threads to run with')
@click.option('--limit-bam-sort-ram', '-r',
              type=int,
              default=2_000_000_000,
              show_default=True,
              help='Max RAM in bytes for STAR BAM sorting (STAR --limitBAMsortRAM). '
                   'Increase if STAR reports a BAM-sorting RAM error.')
# TODO add option for different bam extension
def align(fastq_dir: Path,
          output_dir: Path,
          genome_fasta: Path,
          threads: int,
          limit_bam_sort_ram: int,
          ):
    """Perform read alignment against the reference genome using STAR"""
    click.echo(genome_fasta)
    map_fastq_to_bam(fastq_dir=fastq_dir,
                     output_dir=output_dir, genome_fasta=genome_fasta,
                     threads=threads, limit_bam_sort_ram=limit_bam_sort_ram)
    

@satay.command()
@click.option('--bam-dir', '-b',
              type=click.Path(exists=True, dir_okay=True, path_type=Path),
              required=True,
              help='Directory containing bam files')
@click.option('--output-dir', '-o',
              type=click.Path(exists=True, dir_okay=True, path_type=Path),
              required=True,
              help='Output directory')
@click.option('--sample-name', '-s',
              type=str,
              required=True,
              help='Name of the sample being processed')
@click.option('--gff', '-a',
              type=click.Path(exists=True, path_type=Path),
              multiple=True,
              required=True,
              help='Annotation file(s), ex. gff or bed intervals of interest')
# TODO add option for different bam extension
def map(bam_dir: Path,
        sample_name: str,
        output_dir: Path,
        gff: Tuple[Path, ...],
        ):
    """Identify genomic locations of transposon insertions and compute read counts at each site"""
    click.echo(f"Processing files for sample: {sample_name}")
    click.echo(f"Output directory: {output_dir}")
    map_sample(bam_dir=bam_dir, sample=sample_name,
               output_dir=output_dir, interval_files=gff)


@satay.command()
@click.option('--counts-dir', '-d',
              type=click.Path(exists=True, dir_okay=True, path_type=Path),
              required=True,
              help='Directory containing count files to merge')
@click.option('--gff', '-a',
              type=str,
              required=True,
              help='Name of annotation file used to count # of tns/reads')
@click.option('--name', '-n',
              type=str,
              default='',
              help='Name of the experiment')
@click.option('--format',
              type=click.Choice(['gff', 'bed'], case_sensitive=False),
              default='gff',
              help='Format of the annotation file (gff or bed)')
def merge(counts_dir: Path,
          gff: str,
          format: str,
          name: str,
          ):
    """Merge multiple count files into a single matrix."""
    click.echo(f"Merging files in directory: {counts_dir}")
    click.echo(f"Annotation file name: {gff}")
    click.echo(f"Using format: {format}")
    merge_cnts_files(interval_file=gff, counts_dir=counts_dir,
                     format=format, name=name)


@satay.command()
@click.option('--counts-file', '-f',
              type=click.Path(exists=True, path_type=Path),
              required=True,
              help='Merged counts file')
@click.option('--sample_data', '-s',
              type=click.Path(exists=True, path_type=Path),
              required=True,
              help='Sample data file')
@click.option('--gff', '-a',
              default='',
              help='gff file (optional)')
@click.option('--output-dir', '-o',
              type=click.Path(exists=True, dir_okay=True, path_type=Path),
              required=True,
              help='Output directory')
@click.option('--filter', '-l',
              type=int,
              default=100,
              help='Filter out genes/intervals with less than [100] counts across samples')
@click.option('--alpha',
              type=float,
              default=0.05,
              help='FDR cutoff for DESeq')
@click.option('--comp-col', '-c',
              type=str,
              required=True,
              help='Column in sample data used to define contrasts')
@click.option('--baseline', '-b',
              type=str,
              required=True,
              help='Baseline value, all other options will be compared to this one')
@click.option('--ids',
              multiple=True,
              default=["locus_tag", "gene"],
              help='gff fields to keep')
@click.option('--threads', '-t',
              type=int,
              default=min(8, os.cpu_count() or 1),
              show_default=True,
              help='Number of CPUs for DESeq2 inference')
@click.option('--sample-id-col',
              type=str,
              default='sample_id',
              show_default=True,
              help='Column in sample data holding sample IDs (must match count matrix columns)')
def analyze(counts_file: Path, sample_data: Path,
            output_dir: Path,
            comp_col: str,
            baseline: str,
            filter: int,
            gff: Union[str, Path],
            alpha: float,
            ids: List,
            threads: int,
            sample_id_col: str
            ):
    """Run differential abundance analysis (DESeq2) on the merged count matrix."""
    gene_da(counts_file, sample_data, output_dir,
            filter, comp_col, baseline, alpha,
            gff, list(ids), sample_id_col=sample_id_col, n_cpus=threads)


def main():
    """Main entry point for the CLI."""
    try:
        satay()
    except Exception as e:
        click.echo(f"Error: {str(e)}", err=True)
        raise click.Abort()


if __name__ == '__main__':
    main()
