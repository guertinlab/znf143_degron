#! /usr/bin/bash

#SBATCH --job-name=pro_dreg.sh     # name for job
#SBATCH -N 1                  
#SBATCH -n 1                 
#SBATCH -c 32                  
#SBATCH -p general           
#SBATCH --qos=general       
#SBATCH --mem=32G               
#SBATCH --mail-type=ALL 
#SBATCH --mail-user=guertin@uchc.edu
#SBATCH -o pro_dreg.sh_%j.out
#SBATCH -e pro_dreg.sh_%j.err

module load samtools/1.16.1
module load genometools/1.5.10
module load ucsc_genome/2012.05.22
module load rust
module load bowtie2
module load bedtools
sizes=hg38.chrom.sizes
genome=hg38.fa
ncore=16
read_size=62
seqOutBias=/home/FCAM/mguertin/software/seqOutBias
table=hg38_62.4.2.2.tbl
tallymer=hg38.tal_62.gtTxt.gz

module load samtools/1.16.1
module load genometools/1.5.10
module load ucsc_genome/2012.05.22
module load rust
module load bowtie2
module load bedtools
sizes=hg38.chrom.sizes
genome=hg38.fa
ncore=32

for i in *sorted.bam
do
    name=$(echo $i | awk -F"/" '{print $NF}' | awk -F".sorted.bam" '{print $1}')
    samtools sort -@ $ncore -o ${name}.name.bam $i
done


x=`ls *name*bam | grep -v us`
name=ZNF143_combined
seqOutBias scale $table $x --no-scale --stranded --bed-stranded-positive --bed=$name.bed --out-split-pairends --only-paired --tail-edge --read-size=$read_size --tallymer=$tallymer

#I messed up with the $name variable
mv HEK_CloneZD29_30min_control_rep1.name_minus_PE1.bigWig HEK_CloneZD29_merged_minus_PE1.bigWig
mv HEK_CloneZD29_30min_control_rep1.name_minus_PE2.bigWig HEK_CloneZD29_merged_minus_PE2.bigWig
mv HEK_CloneZD29_30min_control_rep1.name_plus_PE1.bigWig HEK_CloneZD29_merged_plus_PE1.bigWig
mv HEK_CloneZD29_30min_control_rep1.name_plus_PE2.bigWig HEK_CloneZD29_merged_plus_PE2.bigWig
