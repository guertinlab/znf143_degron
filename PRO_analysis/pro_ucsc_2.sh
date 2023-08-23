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


#conda install -c bioconda ucsc-bigwigmerge

for i in *rep1.sorted.bam
do
    nm=$(echo $i | awk -F"/" '{print $NF}' | awk -F".sorted.bam" '{print $1}')
    name=$(echo $nm | awk -F"_rep1" '{print $1}')
    echo $name
    reps=$(ls ${name}_rep*normalized.plus.bigWig | wc -w | bc)
    echo $reps
    filesPlus=$(ls ${name}_rep*normalized.plus.bigWig)
    echo $filesPlus
    bigWigMerge ${name}_rep*normalized.plus.bigWig tmpPlus.bg
    bigWigMerge ${name}_rep*normalized.minus.bigWig tmpMinus.bg
    scaleall=$(bc <<< "scale=4 ; 1.0 / $reps")
    echo scale:
    echo $scaleall
    awk -v scaleall="$scaleall" '{OFS="\t";} {print $1, $2, $3, $4*scaleall}' tmpPlus.bg > ${name}_normalized.plus.bedGraph
    awk -v scaleall="$scaleall" '{OFS="\t";} {print $1, $2, $3, $4*scaleall}' tmpMinus.bg > ${name}_normalized.minus.bedGraph
    rm tmpPlus.bg
    rm tmpMinus.bg
    wigToBigWig ${name}_normalized.plus.bedGraph $sizes ${name}.plus.bigWig
    wigToBigWig ${name}_normalized.minus.bedGraph $sizes ${name}.minus.bigWig
    awk -v var="$name" 'BEGIN {  print "browser position chr11:5,289,521-5,291,937"; print "track type=bedGraph name=plus\"" var "\" description=\"" var "plus_bedGraph\" visibility=full autoScale=on alwaysZero=on color=255,0,0"}  { print $0}' ${name}_normalized.plus.bedGraph > ${name}_header_normalized.plus.bedGraph
    gzip ${name}_header_normalized.plus.bedGraph
    awk -v var="$name" 'BEGIN {  print "browser position chr11:5,289,521-5,291,937"; print "track type=bedGraph name=minus\"" var "\" description=\"" var "minus_bedGraph\" visibility=full autoScale=on alwaysZero=on color=0,0,255"}  { print $0}' ${name}_normalized.minus.bedGraph > ${name}_header_normalized.minus.bedGraph
    gzip ${name}_header_normalized.minus.bedGraph
done

