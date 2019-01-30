##
##
##
require(bigWig)
require(cluster)

load("data-rpkms.RData")
indx <- (tus[,3]-tus[,2])>10000
pkm <- rpkm[indx,]

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
 YW <- ints[ints < 0.5] *2
 WB <- (ints[ints >= 0.5]-0.5) *2
 YW[YW<0] <- 0; WB[WB>1] <- 1
 c(rgb(1, 1, YW), rgb(1-WB, 1-WB, 1))
}

drawCor <- function(var1, var2, Labels) {
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
         hc1 <- hclust(dist(cc, method = "euclidean"),  method="average")
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


pdf("prophaseI-correlation-matrix.pdf")
drawCor(rep(replicate, 3), c(rep(c(stage[1]), 3), rep(c(stage[2]), 3), rep(c(stage[3]), 3)), colnames(rpkm))
dev.off()
}


