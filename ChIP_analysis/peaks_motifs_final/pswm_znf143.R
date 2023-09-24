source('https://raw.githubusercontent.com/guertinlab/seqOutBias/master/docs/R/seqOutBias_functions.R') 
library(ggseqlogo)

pswm.func.2 <- function(x.ligation, out = 'outfilename', ncols = 29) {
  col.matrix <- matrix(unlist(strsplit(as.character(x.ligation), '')), ncol = ncols, byrow = TRUE)
  
  a.nuc <- colSums(col.matrix == "A")
  t.nuc <- colSums(col.matrix == "T")
  c.nuc <- colSums(col.matrix == "C")
  g.nuc <- colSums(col.matrix == "G")
  
  pswm <- cbind(a.nuc, c.nuc, g.nuc, t.nuc)
  
  pswm <- pswm / rowSums(pswm)
  
  outfile <- file(paste0(out, '.txt'))
  on.exit(close(outfile))
  writeLines(c(
    "MEME version 4", 
    "ALPHABET= ACGT", 
    "strands: + -", 
    " ",
    "Background letter frequencies (from uniform background):",
    "A 0.30000 C 0.20000 G 0.20000 T 0.30000", 
    paste("MOTIF", out), 
    " ", 
    paste0("letter-probability matrix: alength= 4 w= ", ncols)
  ), outfile)
  
  write.table(pswm, file = paste0(out, '.txt'), append = TRUE, quote = FALSE, row.names = FALSE, col.names = FALSE)
  
  return(pswm)
}


x = read.table('inferred_mast_ZNF143_PSWM_in_peaks_29mer.fasta', comment.char = ">", header =FALSE)
znf143pswm = pswm.func.2(x[,1], "ZNF143")
write.table(znf143pswm, file = "znf143_denovo_pswm_ACGT.txt", sep="\t", quote = FALSE, col.names =FALSE, row.names=FALSE)
