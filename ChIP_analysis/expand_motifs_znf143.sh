genome=/Users/guertinlab/meds5420_2023/genomes/hg38.fa
sizes=/Users/guertinlab/meds5420_2023/genomes/hg38.chrom.sizes

awk '{OFS="\t";} {if($6 == "+") print $1,$2-14,$3+2,$4,$5,$6; else print $1,$2-3,$3+13,$4,$5,$6}' mast_ZNF143_motif_in_peaks_round1.bed > expanded_mast_ZNF143_motif_in_peaks_round1.bed
fastaFromBed -s -fi $genome -bed expanded_mast_ZNF143_motif_in_peaks_round1.bed > expanded_mast_ZNF143_motif_in_peaks_round1.fasta

awk '{OFS="\t";} {if($6 == "+") print $1,$2-4,$3+10,$4,$5,$6; else print $1,$2-11,$3+3,$4,$5,$6}' mast_ZNF143_motif_in_peaks_round2_RC.bed > expanded_mast_ZNF143_motif_in_peaks_round2_RC.bed
fastaFromBed -s -fi $genome -bed expanded_mast_ZNF143_motif_in_peaks_round2_RC.bed > expanded_mast_ZNF143_motif_in_peaks_round2_RC.fasta

cat expanded_mast_ZNF143_motif_in_peaks_round1.bed expanded_mast_ZNF143_motif_in_peaks_round2_RC.bed > motifs_in_rounds12.bed
slopBed -b 50 -i over40peakIntens_ZNF143peaksChIP.bed -g $sizes > over40peakIntens_ZNF143peaksChIP_101window.bed
intersectBed -v -a over40peakIntens_ZNF143peaksChIP_101window.bed -b motifs_in_rounds12.bed > without_motifs_12.bed

fastaFromBed -fi $genome -bed without_motifs_12.bed -fo without_motifs_12.fasta

meme -oc ZNF143_no_12.meme_output -nmotifs 1 -objfun classic -csites 20000 -searchsize 0 -minw 10 -maxw 20 -revcomp -dna -markov_order 3 -maxsize 100000000 without_motifs_12.fasta
tomtom -no-ssc -oc motif3.tomtom_output -verbosity 1 -min-overlap 5 -dist ed -evalue -eps -thresh 5.0 ZNF143_no_12.meme_output/meme.txt ZNF143_final_PSWM.txt

meme-get-motif -id TCCCACGGCACGGGAACTCC -rc ZNF143_no_12.meme_output/meme.txt > motif3_pswm.txt
#difficult to eye-ball alignment
tomtom -no-ssc -oc motif3_rc.tomtom_output -verbosity 1 -min-overlap 5 -dist ed -evalue -eps -thresh 5.0 motif3_pswm.txt ZNF143_final_PSWM.txt
#this stringency keeps the motif confined to the U RNA regulatory regions
mast -mt 0.000001 -hit_list -best motif3_pswm.txt without_motifs_12.fasta > mast_ZNF143_PSWM_in_peaks_round3.txt
Rscript /Users/guertinlab/rscripts/parse_mast_to_coordinates.R mast_ZNF143_PSWM_in_peaks_round3.txt
fastaFromBed -s -fi $genome -bed mast_ZNF143_PSWM_in_peaks_round3.bed -fo mast_ZNF143_PSWM_in_peaks_round3.fasta
awk '{OFS="\t";} {if($6 == "+") print $1,$2-17,$3-7,$4,$5,$6; else print $1,$2+6,$3+16,$4,$5,$6}' mast_ZNF143_PSWM_in_peaks_round3.bed  > expanded_mast_ZNF143_PSWM_in_peaks_round3.bed 
fastaFromBed -s -fi $genome -bed expanded_mast_ZNF143_PSWM_in_peaks_round3.bed > expanded_mast_ZNF143_PSWM_in_peaks_round3.fasta 


#round 4
intersectBed -v -a without_motifs_12.bed -b expanded_mast_ZNF143_PSWM_in_peaks_round3.bed > without_motifs_123.bed
#make wider window
slopBed -b 70 -i without_motifs_123.bed -g $sizes > without_motifs_123_wide.bed
fastaFromBed -fi $genome -bed without_motifs_123_wide.bed -fo without_motifs_123_wide.fasta
meme -oc ZNF143_no_123.meme_output -nmotifs 1 -objfun classic -csites 20000 -searchsize 0 -minw 5 -maxw 15 -revcomp -dna -markov_order 3 -maxsize 100000000 without_motifs_123_wide.fasta
#all the motifs are on the flanks, as expected.
meme-get-motif -id RRACTACAWYTCCCA -rc ZNF143_no_123.meme_output/meme.txt > motif4_pswm.txt
mast -mt 0.0005 -hit_list -best motif4_pswm.txt without_motifs_123_wide.fasta > mast_ZNF143_PSWM_in_peaks_round4.txt
Rscript /Users/guertinlab/rscripts/parse_mast_to_coordinates.R mast_ZNF143_PSWM_in_peaks_round4.txt
fastaFromBed -s -fi $genome -bed mast_ZNF143_PSWM_in_peaks_round4.bed -fo mast_ZNF143_PSWM_in_peaks_round4.fasta
awk '{OFS="\t";} {if($6 == "+") print $1,$2-15,$3,$4,$5,$6; else print $1,$2-1,$3+14,$4,$5,$6}' mast_ZNF143_PSWM_in_peaks_round4.bed  > expanded_mast_ZNF143_PSWM_in_peaks_round4.bed
fastaFromBed -s -fi $genome -bed expanded_mast_ZNF143_PSWM_in_peaks_round4.bed > expanded_mast_ZNF143_PSWM_in_peaks_round4.fasta

#round 5
intersectBed -wa -a without_motifs_123_wide.bed -b expanded_mast_ZNF143_PSWM_in_peaks_round4.bed > peaks_with_motif4.bed
intersectBed -v -a without_motifs_123.bed -b peaks_with_motif4.bed > without_motifs_1234.bed

slopBed -b 100 -i without_motifs_1234.bed -g $sizes > without_motifs_1234_wide.bed
fastaFromBed -fi $genome -bed without_motifs_1234_wide.bed -fo without_motifs_1234_wide.fasta

#this replaces the central 70 bases with 50 Ns. the number of N is not important, as long as it is greater than the max motif width
awk '!/^>/ { mid = int(length($0) / 2); $0 = substr($0, 1, mid - 35) "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN" substr($0, mid + 36); } 1' without_motifs_1234_wide.fasta > without_motifs_1234_wide_NNN.fasta
meme -oc ZNF143_no_1234.meme_output -nmotifs 1 -objfun classic -csites 20000 -searchsize 0 -minw 5 -maxw 20 -revcomp -dna -markov_order 3 -maxsize 100000000 without_motifs_1234_wide_NNN.fasta
meme-get-motif -id CRRDGCATTVTGGGWA ZNF143_no_1234.meme_output/meme.txt > motif5_pswm.txt
ceqlogo -i motif5_pswm.txt -m CRRDGCATTVTGGGWA -o round5_motif.eps
mast -mt 0.0005 -hit_list -best motif5_pswm.txt without_motifs_1234_wide.fasta > mast_ZNF143_PSWM_in_peaks_round5.txt
Rscript /Users/guertinlab/rscripts/parse_mast_to_coordinates.R mast_ZNF143_PSWM_in_peaks_round5.txt
fastaFromBed -s -fi $genome -bed mast_ZNF143_PSWM_in_peaks_round5.bed -fo mast_ZNF143_PSWM_in_peaks_round5.fasta
awk '{OFS="\t";} {if($6 == "+") print $1,$2-5,$3+9,$4,$5,$6; else print $1,$2-10,$3+4,$4,$5,$6}' mast_ZNF143_PSWM_in_peaks_round5.bed  > expanded_mast_ZNF143_PSWM_in_peaks_round5.bed
fastaFromBed -s -fi $genome -bed expanded_mast_ZNF143_PSWM_in_peaks_round5.bed > expanded_mast_ZNF143_PSWM_in_peaks_round5.fasta

#round 6
intersectBed -wa -a without_motifs_1234_wide.bed -b expanded_mast_ZNF143_PSWM_in_peaks_round5.bed > peaks_with_motif5.bed
intersectBed -v -a without_motifs_1234.bed -b peaks_with_motif5.bed > without_motifs_12345.bed

slopBed -b 50 -i without_motifs_12345.bed -g $sizes > without_motifs_12345_wide.bed
fastaFromBed -fi $genome -bed without_motifs_12345_wide.bed -fo without_motifs_12345_wide.fasta
awk '!/^>/ { mid = int(length($0) / 2); $0 = substr($0, 1, mid - 25) "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN" substr($0, mid + 26); } 1' without_motifs_12345_wide.fasta > without_motifs_12345_wide_NNN.fasta
meme -oc ZNF143_no_12345.meme_output -nmotifs 4 -objfun classic -csites 20000 -searchsize 0 -minw 5 -maxw 20 -revcomp -dna -markov_order 3 -maxsize 100000000 without_motifs_12345_wide_NNN.fasta
streme -oc ZNF143_no_12345.streme_output --maxw 29 --nmotifs 10 --p without_motifs_12345_wide_NNN.fasta


#make note of who these are and we want to see if the queried composite motif scores are also really low or not central:
without_motifs_12345.bed


