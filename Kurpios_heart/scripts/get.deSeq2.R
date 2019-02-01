##
## deSeq2 ...

require(bigWig)

load("data-counts.RData")

library("DESeq2")

FDR.lv<- 0.01

## Add specific genes to the plot.
addlab <- function(GENEID, deRes, genes, ...) {
 idx<-sapply(GENEID, function(GENEID) {
  #ig <- which(genes[,4] == GENEID)
  ig <- which(genes[,7] == GENEID)
  io <- ig[which.min(deRes$padj[ig])]
  if(NROW(io)>0) {
        text(deRes$baseMean[io], deRes$log2FoldChange[io], labels= GENEID, cex= 1, pos= 3, ...)
        points(deRes$baseMean[io], deRes$log2FoldChange[io], col="blue", cex=1.5)
        io
  }
 })
 print(idx)
 idx <- unlist(idx); idx <- idx[!is.null(idx)]
 return(data.frame(Gene= genes[idx,7], AveExpr= deRes$baseMean[idx], logFC= deRes$log2FoldChange[idx], adj.P.Val= deRes$padj[idx]))
}

condition <- c(rep("WT", 4), rep("MUT", 4))
replicate <- factor(rep(c(1:4), 2))
Design <- data.frame(colnames(counts), condition, replicate)

## Workaround for deSeq2 issues ...
#setClassUnion("characterORNULL", c("character", "NULL"))
#setClassUnion("DataTableORNULL", c("DataTable", "NULL"))


## Create DESeq2 object.
dds <- DESeqDataSetFromMatrix(countData= counts, colData= Design, design= ~condition+replicate)

dds <- DESeq(dds)
res_wt_mut <- results(dds, contrast=c("condition","WT","MUT"))


 
print(paste("Number of changes: ", sum(res_wt_mut$padj < FDR.lv, na.rm=TRUE))) ## Number of transcripts.
#print(paste("Number of changes: ", sum(res_p_d$padj < FDR.lv, na.rm=TRUE))) ## Number of transcripts.

 pdf(paste("WT-MUT.pdf", sep=""))
  plotMA(res_wt_mut, alpha= FDR.lv, cex=0.75, ylab= "Log Fold-Change", ylim=c(-5, 5))
  addlab(c("Pitx2"), res_wt_mut, bodies) ## ENSMUSG00000028023.16 is PITX2
  #addlab(c("ENSMUSG00000028023.16"), res_wt_mut, bodies) ## ENSMUSG00000028023.16 is PITX2
 dev.off()



## Write out genes.
PX <- cbind(bodies, res_wt_mut)
PX <- PX[order(PX$padj),]
write.table(PX, "WT-MUT-changed.genes.tsv", sep="\t", row.names=FALSE, col.names=T, quote=FALSE)
#write.table(PX, "WT-MUT-changed.genes.tsv", sep="\t", row.names=T, col.names=T, quote=FALSE)


