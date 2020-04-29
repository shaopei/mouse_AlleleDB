##
##
##
require(bigWig)
require(cluster)

load("data-rpkms_wPHDHET.RData")
load("data-counts_wPHDHET.RData")

indx_counts <- rowSums(counts>10) >= dim(counts)[2]  #every sample has at least 10 reads
indx_trxSize<- (tus[,3]-tus[,2])>10000  # to get a robost signal
indx <- indx_counts & indx_trxSize
rpkm <- rpkm[indx,]
write.table(rpkm, file = "rpkms.txt", quote = F, sep = "\t")

yb.sig.pal <- function(n, scale=10) {
 ints<- c(0:(n-1))/(n-1)   ## Linear scale from 0:1 x N values.
 ints<- 1/(1+exp(scale*(0.5-ints)))## Transfer to sigmoidal scale.
 b<- min(ints)
 m<- 2*b/(n-1)
 ints<- ints+(m*c(0:(n-1)) -b)## Transform by linear function to fill colors out to maxes.

 ## Transfer to colorspace.
 # Yellow: 255, 255, 0
 # White:  255, 255, 255
 # Blue:   0, 0, 255
  hinge.point<-0.1
 
 YW <- ints[ints < hinge.point] *2
 WB <- (ints[ints >= hinge.point]-hinge.point) *2
 YW[YW<0] <- 0; WB[WB>1] <- 1
 c(rgb(1, 1 , YW), rgb(1-WB, 1-WB, 1))
}

drawCor <- function(var1, var2, Labels, hcMethod="average") {
        rpkm_df <- as.matrix(rpkm)

        cond <- as.double(as.factor(var1))
        spec <- as.double(as.factor(var2))
        labs <- Labels

        cc <- cor(rpkm_df, method="spearman")

        pal2 <- c("#567C34", "#D14C29", "#567C34", "#69D270", "#7073C8", "#557571", "#CD9537", "#C0CB88")
        pal1 <- c("#B65BCB", "#C0513A", "#84CA54", "#92C2AF", "#4D4639", "#7B7EB5", "#BDA04D", "#B1517B")
 	pal3 <- c("#E03CE9", "#17B92B", "#E6350D", "#6FD2F0", "#F9F77F", "#5B6C0C", "#68003D", "#310F08")

        ## Print dendrogram and heatmap with latticeExtra.
         library(latticeExtra)
         hc1 <- hclust(dist(cc, method = "euclidean"),  method=hcMethod)
         hc1 <- as.dendrogram(hc1)
         ord.hc1 <- order.dendrogram(hc1)
         hc2 <- reorder(hc1, cond[ord.hc1])
         ord.hc2 <- order.dendrogram(hc2)

        cols <- yb.sig.pal(300, scale=3)
        brks <- seq(min(cc), max(cc), length.out=NROW(cols)+1)

         pl <- levelplot((cc)[ord.hc2, ord.hc2], at= brks, col.regions= cols, xlab="", ylab="", #rev(cm.colors(100)),  # #c("white", "yellow", "blue") # c("#E9F231", "#B1EC2C", "#5DBDEF")
                 colorkey = list(space="left", labels=list(cex=1.5)),
                 scales = list(x= list(rot=90, cex=1.00, labels=labs[ord.hc2]), y=list(draw=FALSE)), #scales = list(x = list(rot = 90)),
                 legend = list(
                        right = list(fun = dendrogramGrob,
                                 args = list(x = hc2, ord = ord.hc2, side = "right", #lwd=2,
                                 size = 7, size.add = 0.5,
                                 add = list(rect = list(col = "transparent", fill = pal3[c(7, 1, 8)][cond])),
                                 type = "rectangle")),
                        top = list(fun = dendrogramGrob,
                                 args = list(x = hc2, ord = ord.hc2, side = "top", #lwd=2,
                                 size = 1, size.add = 0.5,
                                 add = list(rect = list(col = "transparent", fill = pal3[2:6][spec])),
                                 type = "rectangle"))
                                 ))
         print(pl)
}

pdf("Kurpios_heart-correlation-matrix-complete.pdf")
#drawCor(rep(replicate, 3), c(rep(c(stage[1]), 3), rep(c(stage[2]), 3), rep(c(stage[3]), 3)), colnames(rpkm))
drawCor(c(rep(c(stage[1]), 4), rep(c(stage[2]), 4), rep(c(stage[3]), 2)), c(rep(c(stage[1]), 4), rep(c(stage[2]), 4), rep(c(stage[3]), 2)), colnames(rpkm), "complete")
dev.off()

pdf("Kurpios_heart-correlation-matrix-single.pdf")
#drawCor(rep(replicate, 3), c(rep(c(stage[1]), 3), rep(c(stage[2]), 3), rep(c(stage[3]), 3)), colnames(rpkm))
drawCor(c(rep(c(stage[1]), 4), rep(c(stage[2]), 4), rep(c(stage[3]), 2)), c(rep(c(stage[1]), 4), rep(c(stage[2]), 4), rep(c(stage[3]), 2)), colnames(rpkm), "single")
dev.off()

pdf("Kurpios_heart-correlation-matrix-average.pdf")
#drawCor(rep(replicate, 3), c(rep(c(stage[1]), 3), rep(c(stage[2]), 3), rep(c(stage[3]), 3)), colnames(rpkm))
drawCor(c(rep(c(stage[1]), 4), rep(c(stage[2]), 4), rep(c(stage[3]), 2)), c(rep(c(stage[1]), 4), rep(c(stage[2]), 4), rep(c(stage[3]), 2)), colnames(rpkm), "average")
dev.off()


# PCA
rpkm_wPHDHET <- rpkm
pca <- prcomp(rpkm_wPHDHET)
cols <- c(rep("red",4), rep("blue", 4), rep("black", 2))
pch <- c(rep(19,4), rep(17,4), rep(15,2))

data.frame(colnames(rpkm_wPHDHET), cols, pch)
summary(pca)
pdf("PC1.PC2_wPHDHET.pdf", width = 9, height = 8)
par(mar=c(6.1, 7.1, 2.1, 2.1)) #d l u r 5.1, 4.1, 4.1, 2.1
par(mgp=c(3,1,0))
par(cex.lab=2.2, cex.axis=2.2)
plot(y= pca$rotation[,1], x= pca$rotation[,2], col=cols, pch=pch, xlab="PC2", ylab="PC1", cex=1.5, las=1)
legend("topright",col=c("red","blue", "black"),pch=c(19,17,15),
       legend=c("WT","MUT","PHDHET"), bty="n", cex=1.5)
#pairs(pca$rotation[,1:5], col=cols, pch=pch)
dev.off()

pdf("PC1.PC2_wPHDHET_wLabel.pdf", width = 9, height = 8)
par(mar=c(6.1, 7.1, 2.1, 2.1)) #d l u r 5.1, 4.1, 4.1, 2.1
par(mgp=c(3,1,0))
par(cex.lab=2.2, cex.axis=2.2)
plot(y= pca$rotation[,1], x= pca$rotation[,2], col=cols, pch=pch, xlab="PC2", ylab="PC1", cex=1.5, las=1, label=TRUE)
legend("topright",col=c("red","blue", "black"),pch=c(19,17,15),
       legend=c("WT","MUT","PHDHET"), bty="n", cex=1.5)
text(y= pca$rotation[,1], x= pca$rotation[,2] , labels=colnames(rpkm_wPHDHET), pos=3)
#pairs(pca$rotation[,1:5], col=cols, pch=pch)
dev.off()


#pca <- prcomp(rpkm_df[,1:8], center=FALSE, scale=FALSE) ## UNT
#pca <- prcomp(rpkm_df[,9:15], center=FALSE, scale=FALSE) ## PI
#pca <- prcomp(rpkm_df, scale=FALSE, center=FALSE) ## ALL
pca <- prcomp(rpkm_df, scale=TRUE, center=TRUE)
cols <- c(rep("red",3), rep("green",3), rep("blue", 3), rep("dark red", 3), rep("dark green", 3), rep("dark blue", 3), "black", "black")#, "pink", "pink")
pch <- c(rep(19,9), rep(6,9), 9, 24)#, 19, 6)
data.frame(colnames(rpkm_df), cols, pch)

summary(pca) # Prints variance summary for all principal components.

pdf("PC1.PC2.pdf")
plot(y= pca$rotation[,1], x= pca$rotation[,2], col=cols, pch=pch, xlab="PC2", ylab="PC1")
pairs(pca$rotation[,1:5], col=cols, pch=pch)
dev.off()
