#! /usr/bin/bash

#SBATCH --job-name=pro_ucsc.sh     # name for job
#SBATCH -N 1                  
#SBATCH -n 1                 
#SBATCH -c 32                  
#SBATCH -p general           
#SBATCH --qos=general       
#SBATCH --mem=32G               
#SBATCH --mail-type=ALL 
#SBATCH --mail-user=guertin@uchc.edu
#SBATCH -o pro_ucsc.sh_%j.out
#SBATCH -e pro_ucsc.sh_%j.err

module load samtools/1.16.1
module load genometools/1.5.10
module load ucsc_genome/2012.05.22
module load rust
module load bowtie2
module load bedtools
sizes=hg38.chrom.sizes
genome=hg38.fa
ncore=32

for i in *.bam
do
    name=$(echo $i | awk -F"/" '{print $NF}' | awk -F".bam" '{print $1}')
    samtools sort -@ $ncore -n -o ${name}.sorted.bam $i
    bedtools bamtobed -mate1 -i ${name}.sorted.bam -bedpe > ${name}_bed12.bed
    awk '$1==$4 {print $0}' ${name}_bed12.bed | awk '{OFS="\t";} {print $1, $2, $6, $9}' | awk '$1!="." && $4=="+" && $3>$2 && (($3 - $2)<1000) {print $0}' | sort -k1,1 -k2,2n > ${name}_read_span.plus.bed
    awk '$1==$4 {print $0}' ${name}_bed12.bed | awk '{OFS="\t";} {print $1, $2, $6, $9}' | awk '$1!="." && $4=="-" && $3>$2 && (($3 - $2)<1000) {print $0}' | sort -k1,1 -k2,2n > ${name}_read_span.minus.bed
    genomeCoverageBed -bg -i ${name}_read_span.plus.bed -g $sizes > ${name}.plus.bedGraph
    genomeCoverageBed -bg -i ${name}_read_span.minus.bed -g $sizes > ${name}.minus.bedGraph
    #from previous 
    depth=`samtools view -c -F 260 ${name}.sorted.bam`
    scaled=$(bc <<< "scale=3 ; 10000000 / $depth / 2") 
    echo $scaled
    awk -v scaled="$scaled" '{OFS="\t";} {print $1, $2, $3, $4*scaled}' ${name}.plus.bedGraph > ${name}_normalized.plus.bedGraph
    awk -v scaled="$scaled" '{OFS="\t";} {print $1, $2, $3, $4*scaled}' ${name}.minus.bedGraph > ${name}_normalized.minus.bedGraph
    wigToBigWig -clip ${name}_normalized.plus.bedGraph $sizes ${name}_normalized.plus.bigWig
    wigToBigWig -clip ${name}_normalized.minus.bedGraph $sizes ${name}_normalized.minus.bigWig
done
