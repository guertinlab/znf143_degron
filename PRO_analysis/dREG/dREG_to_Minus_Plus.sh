awk '{OFS="\t";} {print $1,$2,$6}' out.dREG.peak.full.bed > out.dREG.peak.minus.bed
awk '{OFS="\t";} {print $1,$6,$3}' out.dREG.peak.full.bed > out.dREG.peak.plus.bed
