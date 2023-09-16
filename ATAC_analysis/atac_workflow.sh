#! /usr/bin/bash

#SBATCH --job-name=cutadapt_atac.sh     # name for job
#SBATCH -N 1                  
#SBATCH -n 1                 
#SBATCH -c 32                  
#SBATCH -p general           
#SBATCH --qos=general       
#SBATCH --mem=32G               
#SBATCH --mail-type=ALL 
#SBATCH --mail-user=guertin@uchc.edu
#SBATCH -o cutadapt_atac.sh_%j.out
#SBATCH -e cutadapt_atac.sh_%j.err

#module load fastq-pair
module load samtools/1.16.1
module load genometools/1.5.10
module load ucsc_genome/2012.05.22
module load rust
module load cutadapt
module load bowtie2
module load bedtools
sizes=/home/FCAM/mguertin/ZNF143_PRO/hg38.chrom.sizes
genome=hg38.fa
ncore=16
read_size=62
seqOutBias=/home/FCAM/mguertin/software/seqOutBias
table=/home/FCAM/mguertin/ZNF143_PRO/hg38_62.4.2.2.tbl
tallymer=/home/FCAM/mguertin/ZNF143_PRO/hg38.tal_62.gtTxt.gz
genome=/home/FCAM/mguertin/ZNF143_PRO/hg38.fa
ncore=32

for i in *PE1.fastq.gz
do
    name=$(echo $i | awk -F"/" '{print $NF}' | awk -F"_PE1" '{print $1}')
    echo $name
    echo unzipping $i
    gunzip $i
    echo unzipping ${name}_PE2.fastq.gz
    gunzip ${name}_PE2.fastq.gz
    cutadapt -a CTGTCTCTTATACACATCT -A CTGTCTCTTATACACATCT -j $ncore -m 10 -O 1 -o ${name}_PE1_no_adapt.fastq -p ${name}_PE2_no_adapt.fastq ${name}_PE1.fastq ${name}_PE2.fastq 2>&1 | tee ${name}_cutadapt.log
done

gunzip chrM.fa.gz
bowtie2-build chrM.fa chrM

#Align to chrM first and remove aligned reads
for i in *_PE1_no_adapt.fastq
do
    name=$(echo $i | awk -F"/" '{print $NF}' | awk -F"_PE1" '{print $1}')
    echo $name
    bowtie2 -p 8 -x chrM -U ${name}_PE1_no_adapt.fastq | samtools sort -n - | samtools fastq -f 0x4 - > ${name}_PE1.chrM.fastq 2>&1 | tee ${name}_chrM_alignment.log
done

for i in *_PE1.chrM.fastq
do
    name=$(echo $i | awk -F"/" '{print $NF}' | awk -F"_PE1" '{print $1}')
    echo $name
    /home/FCAM/mguertin/software/fastq_pair ${name}_PE1.chrM.fastq ${name}_PE2_no_adapt.fastq 2>&1 | tee ${name}_fastq_pair.log
done

genome_index=/home/FCAM/mguertin/ZNF143_PRO/hg38

for i in *_PE1_no_adapt.fastq
do
    name=$(echo $i | awk -F"/" '{print $NF}' | awk -F"_PE1" '{print $1}')
    echo $name
    bowtie2 -p $ncore --maxins 800 -x $genome_index -1 ${name}_PE1.chrM.fastq.paired.fq -2 ${name}_PE2_no_adapt.fastq.paired.fq | samtools view -bS -f 0x2 - | samtools sort -@ $ncore -n - | samtools fixmate -m - - | samtools sort - | samtools markdup -s -r - ${name}.hg38.bam 2>&1 | tee ${name}_bowtie2_hg38.log
done
readLength=62

for i in *.hg38.bam
do
    name=$(echo $i | awk -F"/" '{print $NF}' | awk -F".hg38.bam" '{print $1}')
    echo $name
    seqOutBias scale $table ${name}.hg38.bam --tallymer=$tallymer --no-scale --custom-shift=4,-4 --read-size=${readLength} 2>&1 | tee ${name}_seqOutBias.log
    grep -v "random" ${name}.hg38_not_scaled.bed | grep -v "chrUn" | grep -v "chrEBV" | grep -v "alt" | sort -k1,1 -k2,2n > ${name}_tmp.txt && mv ${name}_tmp.txt ${name}_not_scaled.bed 
    rm ${name}.hg38_not_scaled.bed
done

gzip *q

module load macs2

mkdir temp_macs
macs2 callpeak --call-summits -t *.hg38.bam -n ZNF143_degron_ATAC -g hs -q 0.01 --keep-dup all -f BAM --nomodel --shift -100 --extsize 200 --tempdir temp_macs

wget https://github.com/Boyle-Lab/Blacklist/raw/master/lists/hg38-blacklist.v2.bed.gz
gunzip hg38-blacklist.v2.bed.gz
blacklist=hg38-blacklist.v2.bed
name=ZNF143_degron_ATAC

grep -v "random" ${name}_summits.bed | grep -v "chrUn" | grep -v "chrEBV" | grep -v "chrM" | grep -v "alt" | intersectBed -v -a stdin -b $blacklist > ${name}_tmp.txt

mv ${name}_tmp.txt ${name}_summits.bed

slopBed -b 200 -i ${name}_summits.bed -g $sizes > ${name}_summit_window.bed

sort -k1,1 -k2,2n ${name}_summit_window.bed > ${name}_summit_window_sorted.bed

peaks=${name}_summit_window_sorted.bed

for i in *_not_scaled.bed
do
    nm=$(echo $i | awk -F"/" '{print $NF}' | awk -F"_not_scaled.bed" '{print $1}')
    sort -k1,1 -k2,2n $i > ${nm}_sorted.bed
    mapBed -null '0' -a $peaks -b ${nm}_sorted.bed > ${nm}_peak_counts.txt 
done

for i in *_peak_counts.txt
do
  name=$(echo $i | awk -F"/" '{print $NF}' | awk -F"_peak_counts.txt" '{print $1}')
  awk '{print $NF}' ${name}_peak_counts.txt > ${name}_peak_counts_only.txt
  echo $name | cat - ${name}_peak_counts_only.txt > ${name}_peak_counts.txt
  rm ${name}_peak_counts_only.txt
done

echo -e "chr\tstart\tend\tname\tqvalue" | cat - $peaks | paste -d'\t' - *peak_counts.txt > Combined_ATAC_peak_counts.txt
