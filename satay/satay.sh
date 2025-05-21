#!/bin/bash
cdir=/nfs/cds-peta/exports/biol_micro_cds_gr_sunagawa/scratch/ansintsova/satay/notebooks
cd /nfs/cds-peta/exports/biol_micro_cds_gr_sunagawa/scratch/Projects_NCCR/other/lheistinger/bam
#samples=$(ls -d  20*| awk -F'_' '{print $1 "_" $2}' | sort | uniq)
#echo $samples

# Merging the bams
#for i in $samples; do echo $i; bam_files=$(ls $i*/*Aligned.sortedByCoord.out.bam); echo samtools merge -o meged_bams/$i.bam $bam_files; samtools merge merged_bams/$i.bam $bam_files; done


#Filtering the bams

#for i in merged_bams/*.bam; do echo $i; samtools view -h -F 256,272 -q 10 $i |bedtools bamtobed > $i.bed; done 

#Process based on orientation

#for i in merged_bams/*.bed; do echo $i; $cdir/process_bed.sh $i > $i.insertions; done

# Sort

#for i in merged_bams/*.insertions; do echo $i; sort -k1,1 -k2n,2 $i > $i.sorted; done

# Merge

#for i in merged_bams/*.sorted; do echo $i; bedtools merge -i $i  -s -c 1,6 -o count,distinct > $i.merged; done

# Remove insertions supported by only 1 read

for i in merged_bams/*.sorted.merged; do echo $i;  awk '$4 > 1' $i > $i.filtered; done
#

#Map
interval_files=$(ls ../Int*Off*.bed)

for i in $interval_files
do
  
  file_name_with_ext=$(basename "$i")
  file_name="${file_name_with_ext%.*}"
  echo $file_name
  for j in merged_bams/*merged.filtered
    do
      echo "$j"_"$file_name".cnts
      bedtools map -a $i -b $j -c 1,4 -o count,sum > "$j"_"$file_name".cnts
    done
done


gff_file=/nfs/cds-peta/exports/biol_micro_cds_gr_sunagawa/scratch/Projects_NCCR/other/lheistinger/GCF_000146045.2.genes.gff
gff_name=GCF_000146045
for j in merged_bams/*merged.filtered
do
    echo "$j"_"$gff_name".cnts
    bedtools map -a $gff_file -b $j -c 1,4 -o count,sum > "$j"_"$gff_name".cnts
done