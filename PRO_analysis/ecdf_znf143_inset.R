library(latticeExtra)
library(DESeq2)
library(lattice)
library(dplyr)
library(ggplot2)
library(limma)
library(bigWig)
library(gplots)
library(RColorBrewer)
library(viridis)
source('https://raw.githubusercontent.com/guertinlab/seqOutBias/master/docs/R/seqOutBias_functions.R')
source('https://raw.githubusercontent.com/mjg54/znf143_pro_seq_analysis/master/docs/ZNF143_functions.R')
source('https://raw.githubusercontent.com/guertinlab/genex/main/ChIP_analysis/cdf_functions.R')

# use these versions of cdf.deseq.df and bedTools.Closest
cdf.deseq.df <- function(genes = gene.file, chip.peaks = chip.peaks, cat = "Repressed", opt.str="") {
  bed.tss.activated = get.tss(genes[genes[,5] == cat,])
  bed.tss.unchanged = get.tss(genes[genes[,5] == paste0("Matched to ", cat),])
  
  # note 'functionstring' in bedTools.closest
  act.distance = bedTools.closest(bed1 = bed.tss.activated, bed2 = chip.peaks[,c(1:3)], opt.string = '-D a')
  unreg.distance = bedTools.closest(bed1 = bed.tss.unchanged, bed2 = chip.peaks[,c(1:3)], opt.string = '-D a')

  df.up.can = cbind(act.distance[,c(4, 10)], paste0(opt.str,cat))
  df.un.can = cbind(unreg.distance[,c(4, 10)], paste0(opt.str,"Matched to ", cat))

  colnames(df.up.can) = c(colnames(df.up.can)[1:2], 'status')
  colnames(df.un.can) = c(colnames(df.up.can)[1:2], 'status')

  df.all = rbind(df.up.can, df.un.can)
  df.all$status = factor(df.all$status, levels = c(paste0(opt.str,cat), paste0(opt.str,"Matched to ", cat)))
  return(df.all)
}

bedTools.closest <- function(functionstring="/usr/local/bin/closestBed",bed1,bed2,opt.string="") {
  
  options(scipen =99) # not to use scientific notation when writing out
  
  #write bed formatted dataframes to tempfile
  write.table(bed1,file= 'a.file.bed', quote=F,sep="\t",col.names=F,row.names=F)
  write.table(bed2,file= 'b.file.bed', quote=F,sep="\t",col.names=F,row.names=F)
  
  # create the command string and call the command using system()
  command1=paste('sort -k1,1 -k2,2n', 'a.file.bed', '> a.file.sorted.bed')
  cat(command1,"\n")
  try(system(command1))
  command2=paste('sort -k1,1 -k2,2n', 'b.file.bed', '> b.file.sorted.bed')
  cat(command2,"\n")
  try(system(command2))
  
  command=paste(functionstring,opt.string,"-a",'a.file.sorted.bed',"-b",'b.file.sorted.bed',">",'out.file.bed',sep=" ")
  cat(command,"\n")
  try(system(command))
  
  res=read.table('out.file.bed',sep ="\t", header=F, comment.char='')
  
  command3=paste('rm', 'a.file.bed', 'b.file.bed', 'a.file.sorted.bed', 'b.file.sorted.bed', 'out.file.bed')
  cat(command3,"\n")
  try(system(command3))
  
  colnames(res) = c(colnames(bed1), colnames(bed2), 'dis' )
  return(res)
}

plot_cdf <- function(df.all, tf="quantile", cat = "Repressed", col.lines = c("#ce228e", "grey60", "#2290cf","grey90"), line.type = c(1), cex = 1, abline=0) {
pdf(paste0(tf, "_CDF_", cat, ".pdf"), width=6.2, height=3.83) 
         print(ecdfplot(~log(abs(dis), base = 10), groups = status, data = df.all,
         auto.key = list(lines=TRUE, points=FALSE, cex = cex),
         col = col.lines,
         aspect = 1,
         scales=list(relation="free",alternating=c(1,1,1,1)),
         ylab = 'Cumulative Distribution Function',
         xlab = expression('log'[10]~'ZNF143 Distance from TSS'),
                                        #index.cond = list(c(2,1)),
         between=list(y=1.0),
         type = 'a',
         xlim = c(0,7),
         lwd=2,
         lty=line.type,
         par.settings = list(superpose.line = list(col = col.lines, lwd=3), strip.background=list(col="grey85")),
         panel = function(...) {
             panel.abline(v= log(440, base=10), lty =2) # variable line location
             panel.ecdfplot(...)
         }))
    dev.off()
}



#motifs = read.table(file="funcZNF143_closestPSWM.bedGraph", sep="\t", header=FALSE)
#gene.file=read.table(file = "HEK_ZNF143_gene_classes.bed", sep="\t", header=FALSE)

#df.rep.high = cdf.deseq.df(genes = gene.file, chip.peaks=motifs, cat = "Repressed")
#df.act.high = cdf.deseq.df(genes = gene.file, chip.peaks=motifs, cat = "Activated")

#plot_cdf(df.rep.high, tf = "ZNF143_motif_Repressed", col.lines = c("#ce228e", "grey60", "#2290cf","grey90"), cat = cat, line.type = c(1), cex = 0.75)
#plot_cdf(df.act.high, tf = "ZNF143_motif_Activated", col.lines = c("#ce228e", "grey60", "#2290cf","grey90"), cat = cat, line.type = c(1), cex = 0.75)

#x = rbind(df.rep.high, df.act.high)

#plot_cdf(x, tf = "ZNF143_motif_Both", col.lines = c( "#2290cf", "grey60","#ce228e", "grey90"), cat = cat, line.type = c(1), cex = 0.75)


#motifs = read.table(file="expanded_MAST_motifs_closestBindingSite_shift1.bed", sep="\t", header=FALSE)
#gene.file=read.table(file = "HEK_ZNF143_gene_annotations.bed", sep="\t", header=FALSE)

#tss.interval.bed = merge(potential.tss2, largest.interval.expr.bed, by.x="gene", by.y="gene")



xx = read.table("ZNF143dTag_PRO_DESeq_genes_categorized.bed", sep="\t", header=FALSE)
yy = read.table("HEK_ZNF143_gene_annotations_new_infCoords.bed", sep="\t", header=FALSE)
test = merge(xx, yy, by.x="V4", by.y="V4")[,c(2, 3, 4, 1, 11, 16)]

gene.file = test[test[,5] == "Activated" | test[,5] == "Repressed" | test[,5] == "Matched to Activated" | test[,5] == "Matched to Repressed",] 

df.rep.high = cdf.deseq.df(genes = gene.file, chip.peaks=motifs, cat = "Repressed")
df.act.high = cdf.deseq.df(genes = gene.file, chip.peaks=motifs, cat = "Activated")

plot_cdf(df.rep.high, tf = "ZNF143_motif_Repressed", col.lines = c("#ce228e", "grey60", "#2290cf","grey90"), cat = cat, line.type = c(1), cex = 0.75)
plot_cdf(df.act.high, tf = "ZNF143_motif_Activated", col.lines = c("#ce228e", "grey60", "#2290cf","grey90"), cat = cat, line.type = c(1), cex = 0.75)

x = rbind(df.rep.high, df.act.high)

plot_cdf(x, tf = "ZNF143_motif_Both", col.lines = c( "#2290cf", "grey60","#ce228e", "grey90"), cat = cat, line.type = c(1), cex = 0.75)

#inset
match = ecdf(abs(df.rep.high$dis)[df.rep.high$status == 'Matched to Repressed'])
rep = ecdf(abs(df.rep.high$dis)[df.rep.high$status == 'Repressed'])
match.y = seq(0, 20000, by=20)
rep.y = seq(0, 20000, by=20)

spl = smooth.spline(rep.y, rep(rep.y) - match(match.y))
pred = predict(spl)
pred1 = predict(spl, deriv=1)

print('the distance that CDFs are parallel or converge')
print(rep.y[min(which(pred1$y<=0)) - 1])

pdf("empirical_distance_determination_ZNF143.pdf", width=3.83, height=3.83) 
plot(rep.y, rep(rep.y) - match(match.y),
     xlim = c(0,2000),
     cex=0.7,
     xlab = 'ZNF143 distance from TSS (bp)',
     ylab = 'Repressed Genes CDF - Matched genes CDF')
     abline(v = rep.y[min(which(pred1$y<=0)) - 1], col =2, lty =2)
      lines(pred[[1]], pred[[2]], col = 'blue')
dev.off()


pdf("empirical_distance_determination_ZNF143_wide.pdf", width=3.83, height=3.83) 
plot(rep.y, rep(rep.y) - match(match.y),
     xlim = c(0,10000),
     cex=0.7,
     xlab = 'ZNF143 distance from TSS (bp)',
     ylab = 'Repressed Genes CDF - Matched genes CDF')
     abline(v = rep.y[min(which(pred1$y<=0)) - 1], col =2, lty =2)
      lines(pred[[1]], pred[[2]], col = 'blue')
dev.off()
