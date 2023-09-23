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
intersectBed -v -a without_motifs_123.bed -b expanded_mast_ZNF143_PSWM_in_peaks_round4.bed > without_motifs_1234.bed
slopBed -b 70 -i without_motifs_1234.bed -g $sizes > without_motifs_1234_wide.bed
fastaFromBed -fi $genome -bed without_motifs_1234_wide.bed -fo without_motifs_1234_wide.fasta
#meme -oc ZNF143_no_1234.meme_output -nmotifs 1 -objfun classic -csites 20000 -searchsize 0 -minw 5 -maxw 20 -revcomp -dna -markov_order 3 -maxsize 100000000 without_motifs_1234_wide.fasta

streme -oc ZNF143_no_1234.streme_output --nmotifs 1 --p without_motifs_1234_wide.fasta
mast -mt 0.0005 -hit_list -best ZNF143_no_1234.streme_output/streme.txt without_motifs_1234_wide.fasta > mast_ZNF143_PSWM_in_peaks_round5.txt
Rscript /Users/guertinlab/rscripts/parse_mast_to_coordinates.R mast_ZNF143_PSWM_in_peaks_round5.txt
tr '+-' '-+' < mast_ZNF143_PSWM_in_peaks_round5.bed > mast_ZNF143_PSWM_in_peaks_round5_rc.bed
fastaFromBed -s -fi $genome -bed mast_ZNF143_PSWM_in_peaks_round5_rc.bed -fo mast_ZNF143_PSWM_in_peaks_round5.fasta
awk '{OFS="\t";} {if($6 == "+") print $1,$2-9,$3+10,$4,$5,$6; else print $1,$2-11,$3+8,$4,$5,$6}' mast_ZNF143_PSWM_in_peaks_round5_rc.bed > expanded_mast_ZNF143_PSWM_in_peaks_round5_rc.bed
fastaFromBed -s -fi $genome -bed expanded_mast_ZNF143_PSWM_in_peaks_round5_rc.bed > expanded_mast_ZNF143_PSWM_in_peaks_round5_rc.fasta

#round 6
intersectBed -v -a without_motifs_1234.bed -b expanded_mast_ZNF143_PSWM_in_peaks_round5_rc.bed > without_motifs_12345.bed
#slopBed -b -30 -i without_motifs_12345.bed -g $sizes > without_motifs_12345_narrow.bed
fastaFromBed -fi $genome -bed without_motifs_12345.bed -fo without_motifs_12345.fasta
streme -oc ZNF143_no_12345.streme_output --nmotifs 2 --p without_motifs_12345.fasta
meme-get-motif -id 2-CACTTCCGGGG -rc ZNF143_no_12345.streme_output/streme.txt > motif6_pswm.txt
ceqlogo -i motif6_pswm.txt -m  2-CACTTCCGGGG -o round6_motif.eps
tomtom -no-ssc -oc motif6.tomtom_output -verbosity 1 -min-overlap 5 -dist ed -evalue -eps -thresh 1.0 motif6_pswm.txt ZNF143_final_PSWM.txt
mast -mt 0.0005 -hit_list -best motif6_pswm.txt without_motifs_12345.fasta > mast_ZNF143_PSWM_in_peaks_round6.txt
Rscript /Users/guertinlab/rscripts/parse_mast_to_coordinates.R mast_ZNF143_PSWM_in_peaks_round6.txt
fastaFromBed -s -fi $genome -bed mast_ZNF143_PSWM_in_peaks_round6.bed -fo mast_ZNF143_PSWM_in_peaks_round6.fasta
awk '{OFS="\t";} {if($6 == "+") print $1,$2-13,$3+6,$4,$5,$6; else print $1,$2-7,$3+12,$4,$5,$6}' mast_ZNF143_PSWM_in_peaks_round6.bed  > expanded_mast_ZNF143_PSWM_in_peaks_round6.bed
fastaFromBed -s -fi $genome -bed expanded_mast_ZNF143_PSWM_in_peaks_round6.bed > expanded_mast_ZNF143_PSWM_in_peaks_round6.fasta

#round 7
intersectBed -v -a without_motifs_12345.bed -b expanded_mast_ZNF143_PSWM_in_peaks_round6.bed > without_motifs_123456.bed
slopBed -b 100 -i without_motifs_123456.bed -g $sizes > without_motifs_123456_wide.bed
fastaFromBed -fi $genome -bed without_motifs_123456_wide.bed -fo without_motifs_123456_wide.fasta
streme -oc ZNF143_no_123456.streme_output --maxw 29 --nmotifs 5 --p without_motifs_123456_wide.fasta
tomtom -no-ssc -oc motif7.tomtom_output -verbosity 1 -min-overlap 5 -dist ed -evalue -eps -thresh 1.0 ZNF143_no_123456.streme_output/streme.txt ZNF143_final_PSWM.txt

meme-get-motif -id 4-RACTACAHTTCCCAGMAKSCHH -rc ZNF143_no_123456.streme_output/streme.txt > motif7_pswm.txt
ceqlogo -i motif7_pswm.txt -m 4-RACTACAHTTCCCAGMAKSCHH -o round7_motif.eps
mast -mt 0.0005 -hit_list -best motif7_pswm.txt without_motifs_123456_wide.fasta > mast_ZNF143_PSWM_in_peaks_round7.txt
Rscript /Users/guertinlab/rscripts/parse_mast_to_coordinates.R mast_ZNF143_PSWM_in_peaks_round7.txt
