library(bigWig)                                        
library(DESeq2)
library(lattice)
library(dplyr)
library(ggplot2)
library(limma)

get.raw.counts.dREG <- function(df.prefix, path.to.bigWig, file.prefix = 'HEK') {
    vec.names = c()
    df.plus = read.table(paste0(df.prefix, '.plus.bed'))
    #dREG gives summits outside range for some reason
    df.plus[,2] = df.plus[,2] - 11
    df.plus[,3] = df.plus[,3] + 100
    df.plus[,2][df.plus[,2] < 0] = 1
    df.minus = read.table(paste0(df.prefix, '.minus.bed'))
    df.minus[,3] = df.minus[,3] + 11
    df.minus[,2] = df.minus[,2] -100
    df.minus[,2][df.minus[,2] < 0] = 1
    print('plus')
    print(head(df.plus))
    print('minus')
    print(head(df.minus))
    inten.df=data.frame(matrix(ncol = 0, nrow = nrow(df.plus)))
    for (mod.bigWig in Sys.glob(file.path(path.to.bigWig, 
                       paste(file.prefix, "*_plus_PE1.bigWig", sep ='')))) {
        factor.name = strsplit(strsplit(mod.bigWig, 
                        "/")[[1]][length(strsplit(mod.bigWig, "/")[[1]])], '_plus')[[1]][1]
        print(factor.name)
        vec.names = c(vec.names, factor.name)
        loaded.bw.plus = load.bigWig(mod.bigWig)
        print(mod.bigWig)
        print(paste(path.to.bigWig,'/',factor.name, '_minus_PE1.bigWig', sep=''))
        loaded.bw.minus = load.bigWig(paste(path.to.bigWig,'/',factor.name, 
                        '_minus_PE1.bigWig', sep=''))
        mod.inten.plus = bed.region.bpQuery.bigWig(loaded.bw.plus, df.plus, abs.value = TRUE)
        mod.inten.minus = bed.region.bpQuery.bigWig(loaded.bw.minus, df.minus, abs.value = TRUE)
        inten.df = cbind(inten.df, mod.inten.plus + mod.inten.minus)
    }
    colnames(inten.df) = vec.names
    r.names = paste(df.plus[,1], ':', df.plus[,2] + 11, '-', df.plus[,2] + 12, sep='')
    row.names(inten.df) = r.names
    return(inten.df)
}




source('https://raw.githubusercontent.com/guertinlab/seqOutBias/master/docs/R/seqOutBias_functions.R')
source('https://raw.githubusercontent.com/mjg54/znf143_pro_seq_analysis/master/docs/ZNF143_functions.R')

get.raw.counts.interval <- function(df, path.to.bigWig, file.prefix = 'M') {
    vec.names = c()
    inten.df=data.frame(matrix(ncol = 0, nrow = nrow(df)))
    
    for (mod.bigWig in Sys.glob(file.path(path.to.bigWig, paste(file.prefix, "*plus_PE1.bigWig", sep ='')))) {
        factor.name = strsplit(strsplit(mod.bigWig, "/")[[1]][length(strsplit(mod.bigWig, "/")[[1]])], '_plus')[[1]][1]
        print(factor.name)
        vec.names = c(vec.names, factor.name)
        loaded.bw.plus = load.bigWig(mod.bigWig)
        print(mod.bigWig)
        print(paste(path.to.bigWig,'/',factor.name, '_minus.bigWig', sep=''))
        loaded.bw.minus = load.bigWig(paste(path.to.bigWig,'/',factor.name, '_minus_PE1.bigWig', sep=''))
        mod.inten = bed6.region.bpQuery.bigWig(loaded.bw.plus, loaded.bw.minus, df)
        inten.df = cbind(inten.df, mod.inten)
    }
    colnames(inten.df) = vec.names
    r.names = paste(df[,1], ':', df[,2], '-', df[,3],'_', df[,4], sep='')
    row.names(inten.df) = r.names
    return(inten.df)
}


get.raw.counts.genic <- function(df.prefix, path.to.bigWig, file.prefix = 'HEK') {
    vec.names = c()
    df.plus = read.table(df.prefix)
    #dREG gives summits outside range for some reason
    df.plus[,2] = df.plus[,2] - 11
    df.plus[,3] = df.plus[,3] + 200
    df.plus[,2][df.plus[,2] < 0] = 1
    df.minus = read.table(df.prefix)
    df.minus[,3] = df.minus[,3] + 11
    df.minus[,2] = df.minus[,2] -200
    df.minus[,2][df.minus[,2] < 0] = 1
    print('plus')
    print(head(df.plus))
    print('minus')
    print(head(df.minus))
    inten.df=data.frame(matrix(ncol = 0, nrow = nrow(df.plus)))
    for (mod.bigWig in Sys.glob(file.path(path.to.bigWig, 
                       paste(file.prefix, "*_plus_PE1.bigWig", sep ='')))) {
        factor.name = strsplit(strsplit(mod.bigWig, 
                        "/")[[1]][length(strsplit(mod.bigWig, "/")[[1]])], '_plus')[[1]][1]
        print(factor.name)
        vec.names = c(vec.names, factor.name)
        loaded.bw.plus = load.bigWig(mod.bigWig)
        print(mod.bigWig)
        print(paste(path.to.bigWig,'/',factor.name, '_minus_PE1.bigWig', sep=''))
        loaded.bw.minus = load.bigWig(paste(path.to.bigWig,'/',factor.name, 
                        '_minus_PE1.bigWig', sep=''))
        mod.inten.plus = bed.region.bpQuery.bigWig(loaded.bw.plus, df.plus, abs.value = TRUE)
        mod.inten.minus = bed.region.bpQuery.bigWig(loaded.bw.minus, df.minus, abs.value = TRUE)
        inten.df = cbind(inten.df, mod.inten.plus + mod.inten.minus)
    }
    colnames(inten.df) = vec.names
    r.names = paste(df.plus[,1], ':', df.plus[,2] + 11, '-', df.plus[,2] + 12, sep='')
    row.names(inten.df) = r.names
    return(inten.df)
}


inferred.coords=read.table('/Users/guertinlab/Desktop/ZNF143_ChIP/bigWigs_for_dREG/HEK_ZNF143_gene_annotations_new_infCoords.bed', header =FALSE, sep ="\t")

counts.df.genes = abs(get.raw.counts.interval(inferred.coords, "/Users/guertinlab/Desktop/ZNF143_ChIP/bigWigs_for_dREG/counts/", file.prefix = "H"))
#load from Jinhong

sf.dREG = estimateSizeFactorsForMatrix(counts.df.genes)

load('~/Desktop/ZNF143_ChIP/actDownMatchedRepDown.Rdata')

write.table(actDownMatchedRepDown[,c(7:9)][actDownMatchedRepDown$status == "Activated",], file = "ZNF143_genic_Activated.bed", quote = FALSE, col.names =FALSE, row.names=FALSE, sep = "\t")  
write.table(actDownMatchedRepDown[,c(7:9)][actDownMatchedRepDown$status == "Matched to Activated",], file = "ZNF143_genic_MatchedActivated.bed", quote = FALSE, col.names =FALSE, row.names=FALSE, sep = "\t")  
write.table(actDownMatchedRepDown[,c(7:9)][actDownMatchedRepDown$status == "Repressed",], file = "ZNF143_genic_Repressed.bed", quote = FALSE, col.names =FALSE, row.names=FALSE, sep = "\t")  


direc = '/Users/guertinlab/Desktop/ZNF143_ChIP/bigWigs_for_dREG/counts'

counts.genic.act = get.raw.counts.genic('ZNF143_genic_Activated.bed', direc, file.prefix = 'HEK')

sample.conditions = factor(sapply(strsplit(colnames(counts.genic.act), '_'), '[', 4), levels=c("control","dTAGV1"))
rep = factor(sapply(strsplit(colnames(counts.genic.act), 'rep'), '[', 2))


deseq.counts.table = DESeqDataSetFromMatrix(countData = counts.genic.act,
                colData = cbind.data.frame(sample.conditions, rep), 
                design = ~ sample.conditions)

sizeFactors(deseq.counts.table) = sf.dREG

deseq.counts.table 

dds.act = DESeq(deseq.counts.table)

act.df = results(dds.act)

nrow(act.df[act.df$padj < 0.1 & !is.na(act.df$padj) & act.df$log2FoldChange < 0,])
                                        #8/28

system(paste("awk '!seen[$0]++' ", "ZNF143_genic_MatchedActivated.bed >ZNF143_genic_MatchedActivated_unique.bed"))
counts.genic.matched = get.raw.counts.genic('ZNF143_genic_MatchedActivated_unique.bed', direc, file.prefix = 'HEK')

sample.conditions = factor(sapply(strsplit(colnames(counts.genic.matched), '_'), '[', 4), levels=c("control","dTAGV1"))
rep = factor(sapply(strsplit(colnames(counts.genic.matched), 'rep'), '[', 2))


deseq.counts.table = DESeqDataSetFromMatrix(countData = counts.genic.matched,
                colData = cbind.data.frame(sample.conditions, rep), 
                design = ~ sample.conditions)

deseq.counts.table 

dds.matched = DESeq(deseq.counts.table)

matched.df = results(dds.matched)

nrow(matched.df[matched.df$padj < 0.1 & !is.na(matched.df$padj) & matched.df$log2FoldChange < 0,])
matched.df
                                        #0/30


counts.genic.rep = get.raw.counts.genic('ZNF143_genic_Repressed.bed', direc, file.prefix = 'HEK')

sample.conditions = factor(sapply(strsplit(colnames(counts.genic.rep), '_'), '[', 4), levels=c("control","dTAGV1"))
rep = factor(sapply(strsplit(colnames(counts.genic.rep), 'rep'), '[', 2))


deseq.counts.table = DESeqDataSetFromMatrix(countData = counts.genic.rep,
                colData = cbind.data.frame(sample.conditions, rep), 
                design = ~ sample.conditions)

sizeFactors(deseq.counts.table) = sf.dREG

deseq.counts.table 

dds.rep = DESeq(deseq.counts.table)

rep.df = results(dds.rep)

nrow(rep.df[rep.df$padj < 0.1 & !is.na(rep.df$padj) & rep.df$log2FoldChange < 0,])
rep.df


rep.df[rep.df$padj < 0.1 & !is.na(rep.df$padj) & rep.df$log2FoldChange < 0,]





setwd('/Users/guertinlab/Desktop/ZNF143_ChIP/ZNF143_dREG_peak_calling_on_Aug_23_2023_1250_PM_ARCHIVE')

direc = '/Users/guertinlab/Desktop/ZNF143_ChIP/bigWigs_for_dREG/counts'

counts.dreg = get.raw.counts.dREG('ZNF143_total.dREG.peak', 
                                  direc, file.prefix = 'HEK')


colnames(counts.dreg) = sapply(strsplit(colnames(counts.dreg), 'HEK_CloneZD29_30min_'), '[', 2)
#colnames(counts.dreg) = sapply(strsplit(colnames(counts.dreg), '_pro'), '[', 1)

save(counts.dreg, file = 'znf143.counts.dreg.Rdata')



sample.conditions = factor(sapply(strsplit(colnames(counts.dreg), '_'), '[', 1), levels=c("control","dTAGV1"))
rep = factor(sapply(strsplit(colnames(counts.dreg), 'rep'), '[', 2))


deseq.counts.table = DESeqDataSetFromMatrix(countData = counts.dreg,
                colData = cbind.data.frame(sample.conditions, rep), 
                design = ~ sample.conditions)

sizeFactors(deseq.counts.table) = sf.dREG

deseq.counts.table 

dds = DESeq(deseq.counts.table)
#dds
#save(dds, file = "genex_ATAC_dds.R")

normalized.counts = counts(dds, normalized=TRUE)
    
rld = rlog(dds, blind=TRUE)
#save(rld, file = "genex_ATAC_rld.R")

# plot principle components 

pca.plot = plotPCA(rld, intgroup="sample.conditions", returnData=TRUE)

#pca.plot 

plotPCAlattice(pca.plot, file = 'PCA_plot_30min_ZNF143_degradation-dREG.pdf')
        
DE.results = results(dds)

head(DE.results)

DE.results.lattice = 
    categorize.deseq.df(DE.results, 
                        fdr = 0.1, log2fold = 0.0, treat = 'ZNF143_degron_30min')


ma.plot.lattice(DE.results.lattice, filename = '30min_ZNF143_degradation dREG', 
        title.main = "Differential Bidirectional TXN")   

repressed.all = DE.results.lattice[DE.results.lattice$response == 
                                             'ZNF143_degron_30min Repressed',]

activated.all = DE.results.lattice[DE.results.lattice$response == 
                                             'ZNF143_degron_30min Activated',]

chr = sapply(strsplit(rownames(activated.all), ":"), "[", 1)
rnge = sapply(strsplit(rownames(activated.all), ":"), "[", 2)
start = as.numeric(sapply(strsplit(rnge, "-"), "[", 1)) - 50
end = as.numeric(sapply(strsplit(rnge, "-"), "[", 2)) + 50

write.table(cbind(chr, start, end), file = "dREG_ZNF143_degron_activated.bed", quote = FALSE,
col.names =FALSE, row.names=FALSE, sep = "\t")


chr = sapply(strsplit(rownames(repressed.all), ":"), "[", 1)
rnge = sapply(strsplit(rownames(repressed.all), ":"), "[", 2)
start = as.numeric(sapply(strsplit(rnge, "-"), "[", 1)) - 50
end = as.numeric(sapply(strsplit(rnge, "-"), "[", 2)) + 50

write.table(cbind(chr, start, end), file = "dREG_ZNF143_degron_repressed.bed", quote = FALSE,
col.names =FALSE, row.names=FALSE, sep = "\t")

chr = sapply(strsplit(rownames(DE.results.lattice), ":"), "[", 1)
rnge = sapply(strsplit(rownames(DE.results.lattice), ":"), "[", 2)
start = as.numeric(sapply(strsplit(rnge, "-"), "[", 1)) - 50
end = as.numeric(sapply(strsplit(rnge, "-"), "[", 2)) + 50

write.table(cbind(chr, start, end, DE.results.lattice), file = "dREG_ZNF143_degron_changes.bed", quote = FALSE,
col.names =FALSE, row.names=FALSE, sep = "\t")


chr = sapply(strsplit(rownames(DE.results.lattice), ":"), "[", 1)
rnge = sapply(strsplit(rownames(DE.results.lattice), ":"), "[", 2)
start = as.numeric(sapply(strsplit(rnge, "-"), "[", 1)) 
end = as.numeric(sapply(strsplit(rnge, "-"), "[", 2)) 
write.table(cbind(chr, start, end, DE.results.lattice), file = "dREG_ZNF143_degron_changes_summit.bed", quote = FALSE, col.names =FALSE, row.names=FALSE, sep = "\t")
