#!/usr/bin/env Rscript

# Edited from 'parse_mast_to_coordinates.R' to allow for 
# parsing of MAST results containing multiple different motifs

args = commandArgs(trailingOnly=TRUE)

file = args[[1]]
filename = paste(strsplit(file, '.txt')[[1]][1], '.bed', sep='')
print(filename)

if(length(args) == 1){
  numMotifs=1
} else {
  numMotifs = as.numeric(args[2])
}

print(numMotifs)

parse.mast <- function(file, motif.num = 1) {
  mast.data = read.table(file, colClasses = c('character','character','character','character', 'integer','integer','numeric','numeric'))
  chrom = vector(mode="character", length = nrow(mast.data))
  start = vector(mode="integer", length = nrow(mast.data))
  end = vector(mode="integer", length = nrow(mast.data))
  strand.motif = vector(mode="integer", length = nrow(mast.data))
  score = vector(mode="numeric", length = nrow(mast.data))
  pval = vector(mode="numeric", length = nrow(mast.data))

  for (j in 1:nrow(mast.data)) {
    chrom[j] = strsplit(mast.data[j,1], ":")[[1]][1]
    start[j] = as.numeric(strsplit(strsplit(mast.data[j,1], ":")[[1]][2], "-")[[1]][1]) + mast.data[j,5]
    end[j] = as.numeric(strsplit(strsplit(mast.data[j,1], ":")[[1]][2], "-")[[1]][1]) + mast.data[j,6]
    strand.motif[j] = as.character(mast.data[j,2])
    score[j] = mast.data[j,7]
    pval[j] = mast.data[j,8]
  }
  
  arg = data.frame(cbind(chrom, start, end, score, pval, strand.motif))
  arg = arg[arg$strand.motif == paste('+', motif.num, sep='') | arg$strand.motif ==  paste('-', motif.num, sep=''),]

  strand = vector(mode="character", length = nrow(arg))
  motif = vector(mode="character", length = nrow(arg))
  for (j in 1:nrow(arg)) {
    strand[j] = strsplit(as.character(arg[j,6]), "")[[1]][1]
    motif[j] = strsplit(as.character(arg[j,6]), "")[[1]][2]
  }
  res = cbind(arg[,c(T,T,T,T,T,F)], strand, motif)
  colnames(res) = c('chr', 'start', 'end', 'score', 'pval', 'strand', 'motif')
  res$start = as.numeric(as.character(res$start))
  res$end = as.numeric(as.character(res$end))

  return(res)
}

finalBed = data.frame(matrix(nrow = 0, ncol = 7))
colnames(finalBed) = c('chr', 'start', 'end', 'score', 'pval', 'strand', 'motif')
for (m in 1:numMotifs){
  print(paste0("Processing motif ",m))
  finalBed = rbind(finalBed,parse.mast(file, m))
}


write.table(finalBed, file=filename, sep='\t', quote=F, row.names=F, col.names=F)

