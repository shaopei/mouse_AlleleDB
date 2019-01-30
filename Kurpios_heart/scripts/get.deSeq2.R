##
## deSeq2 ...

require(bigWig)

load("data-counts.RData")

library("DESeq2")

FDR.lv<- 0.01

## Add specific genes to the plot.
addlab <- function(GENEID, deRes, genes, ...) {
 idx<-sapply(GENEID, function(GENEID) {
  ig <- which(genes[,4] == GENEID)
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

condition <- c(rep("LZ", 4), rep("P", 4), rep("D", 4))
replicate <- factor(rep(c(1:4), 3))
Design <- data.frame(colnames(counts), condition, replicate)

## Workaround for deSeq2 issues ...
setClassUnion("characterORNULL", c("character", "NULL"))
setClassUnion("DataTableORNULL", c("DataTable", "NULL"))

## Create DESeq2 object.
dds <- DESeqDataSetFromMatrix(countData= counts, colData= Design, design= ~condition+replicate)

dds <- DESeq(dds)
res_lz_p <- results(dds, contrast=c("condition","LZ","P"))
res_p_d  <- results(dds, contrast=c("condition","P","D"))

 
print(paste("Number of changes: ", sum(res_lz_p$padj < FDR.lv, na.rm=TRUE))) ## Number of transcripts.
print(paste("Number of changes: ", sum(res_p_d$padj < FDR.lv, na.rm=TRUE))) ## Number of transcripts.

 pdf(paste("results/LZ-P.pdf", sep=""))
  plotMA(res_lz_p, alpha= FDR.lv, cex=0.75, ylab= "Log Fold-Change", ylim=c(-5, 5))
  addlab(c("Adam3", "Zfy1", "Adam18", "Stag3"), res_lz_p, bodies) ##
 dev.off()

 pdf(paste("results/P-D.pdf", sep=""))
  plotMA(res_p_d, alpha= FDR.lv, cex=0.75, ylab= "Log Fold-Change", ylim=c(-5, 5))
  addlab(c("Adam3", "Zfy1", "Adam18", "Stag3"), res_p_d, bodies) ##
 dev.off()

## Write out genes.
PX <- cbind(bodies, res_lz_p)
PX <- PX[order(PX$padj),]
write.table(PX, "results/LZ-P-changed.genes.tsv", sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)

PX <- cbind(bodies, res_p_d)
PX <- PX[order(PX$padj),]
write.table(PX, "results/P-D-changed.genes.tsv", sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)




