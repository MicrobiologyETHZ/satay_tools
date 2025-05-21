import pyranges as pr
import pandas as pd
import typing
from typing import Union, List, Dict
from pathlib import Path
"""
Part 1: Count transposons over different bed intervals, generate 9 count files

1. Generate bed file of size 500 across genome
2. Calculate intersection with tn posistions
3. ouptut count file
4. Do this 9 times

"""
# Chromosome names for <> genome, release version <>

chr_names = """Chromosome	Genbank ID	RefSeq ID	Length (bp)
Chromosome I	BK006935.2	NC_001133.9	230218
Chromosome II	BK006936.2	NC_001134.8	813184
Chromosome III	BK006937.2	NC_001135.5	316620
Chromosome IV	BK006938.2	NC_001136.10	1531933
Chromosome V	BK006939.2	NC_001137.3	576874
Chromosome VI	BK006940.2	NC_001138.5	270161
Chromosome VII	BK006941.2	NC_001139.9	1090940
Chromosome VIII	BK006934.2	NC_001140.6	562643
Chromosome IX	BK006942.2	NC_001141.2	439888
Chromosome X	BK006943.2	NC_001142.9	745751
Chromosome XI	BK006944.2	NC_001143.9	666816
Chromosome XII	BK006945.2	NC_001144.5	1078177
Chromosome XIII	BK006946.2	NC_001145.3	924431
Chromosome XIV	BK006947.3	NC_001146.8	784333
Chromosome XV	BK006948.2	NC_001147.6	1091291
Chromosome XVI	BK006949.2	NC_001148.4	948066
Chromosome Mito	AJ011856.1	NC_001224.1	85779""".split('\n')

chr_names = [ch.split("\t") for ch in chr_names if 'Mito' not in ch]
chr_map = {f"chr{c[0].split()[1]}": int(c[3]) for c in chr_names[1:]}


def chr_range(chr_name: str, chr_len: int, start: int = 1, size: int = 500) -> dict:
    """
    Given a chromosome name and length, generate ranges of a given size, starting at a given position
    :param chr_name: chromosome name
    :param chr_len: chromosome length
    :param start: starting position for the intervals
    :param size: size of the intervals
    :return: dictionary of intervals
    """
    starts = list(range(start,  chr_len, size))
    ends = list(range(start+size-1,  chr_len, size)) + [chr_len]
    chrs = [chr_name]*len(starts)
    return {'Chromosome': chrs, "Start": starts, "End": ends}


def create_ranges(chr_map: dict, start: int = 1, size: int = 500):
    """
    Given all a dict of chromosome names and length, create ranges of given size for each, concat into a df
    :param chr_map: dict of {name: length}
    :param start: starting position for the intervals
    :param size: size of the intervals
    :return: ranges dataframe
    """
    genomic_ranges = []
    for chr, chr_len in chr_map.items():
        chr_ranges_dict = chr_range(chr, chr_len, start, size)
        chr_ranges = pr.from_dict(chr_ranges_dict)
        genomic_ranges.append(chr_ranges)
    ranges_df = pr.concat(genomic_ranges)
    return ranges_df


def get_tn_positions(bed_file: Union[str, Path]):
    """
    reads 'bed' file provided by the collaborators -> not actually bed, sep = ' '
    and some files start with 'track name ....' -> not even commented out ...
    :param bed_file:
    :return: pyranges dataframe
    """
    df = pd.read_table(bed_file, names="Chromosome Start End Strand Score".split(), sep=' ',
                       comment="t")
    gr = pr.PyRanges(df)
    return gr


def get_tn_sites_for_all_samples(sample_files: List[Path]) -> Dict:
    """
    Take list of sample files and convert to dictionary of genomic ranges
    :param sample_files:
    :return: dictionary of tn insertion sites for each sample
    """
    tn_sites = {}
    for file in sample_files:
        name = file.parent.name.replace("_", "-")
        print(name)
        tn_positions = get_tn_positions(file)
        tn_sites[name] = tn_positions
    return tn_sites


def count_tns(tn_sites: Dict, chr_gr):
    """
    Take dictionary of tn insertion sites and genomic ranges and calculate overlaps
    :param tn_sites:
    :param chr_gr:
    :return:
    """
    return pr.count_overlaps(tn_sites, chr_gr)


def count_tns_over_custom_range(chr_map, start, size, sample_files, out_file):
    # Create custom range
    custom_range = create_ranges(chr_map, start, size)

    # Get all the transposon sites
    tn_sites = get_tn_sites_for_all_samples(sample_files)

    # Count overlaps
    overlaps = count_tns(tn_sites, custom_range)
    overlaps.to_csv(out_file)
    return overlaps

"""

Part 2: Count reads over different intervals, generate 9 count files

1. Generate gtf/saf file of size 500 across genome

GeneID	Chr	Start	End	Strand
497097	chr1	3204563	3207049	-
497097	chr1	3411783	3411982	-
497097	chr1	3660633	3661579	-


2. featureCounts the reads
3. ouptut count file
4. Do this 9 times

"""




def main():
    pass

if __name__ == "__main__":
    data_dir = Path("/Users/ansintsova/git_repos/parfenova_satay/data/")
    out_dir = data_dir/"02_22_tn_counts"
    samples = [data_dir/"15776-1_noaF/15776-1_noaF.bam.bed", data_dir/"15776-1_4nMaF/15776-1_4nMaF.bam.bed"]
    starts = [1] + list(range(50, 500, 50))
    size = 500
    for start in starts:
        overlaps = count_tns_over_custom_range(chr_map, start, size, samples, out_dir/f"15776-1-4aF-0aF-{size}-{start}.csv")

