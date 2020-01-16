library(scales)
# MDS plot #
setwd("../data")
#load results from biohpc #
load("MDSoutput_update_v1.Rdata")

#### 1. case group ####
par(mar = c(4, 4, 1, 1))
ord.real = real.nmds[["points"]]; ord.wgan = wgan.nmds[["points"]]; ord.norta = NorTA.nmds[["points"]]; 
xrange = c(min(ord.real[, 1],ord.wgan[, 1], ord.norta[, 1] ) - 0.4, 
           max(ord.real[, 1],ord.wgan[, 1], ord.norta[, 1] ) + 0.4)
yrange = c(min(ord.real[, 2],ord.wgan[, 2], ord.norta[, 2] ) - 0.4, 
           max(ord.real[, 2],ord.wgan[, 2], ord.norta[, 2] ) + 0.4)
plot(1, 1, type="n", main=" ",
     xlim = xrange, ylim = yrange, 
     xlab = "nMDS1", ylab = "nMDS2")
points(NorTA.nmds$points[, 1], NorTA.nmds$points[, 2], pch=20, col = alpha(rgb(0,1,0,0.5), 0.1))
points(real.nmds$points[, 1], real.nmds$points[, 2], pch=20, col = "blue")
points(wgan.nmds$points[, 1], wgan.nmds$points[, 2], pch=20, col = alpha(rgb(1,0,0,0.5), 0.3))
# add legend
legend("bottomright", 
       legend=c("Real data","MB-GAN","NorTA"), 
       col=c("blue",rgb(1,0,0,0.5),rgb(0,1,0,0.5)), pch=20,
       cex = 0.8,bty = "n")


#### 2. ctrl group ####
ord.real = real.nmds.ctrl[["points"]];
ord.wgan = wgan.nmds.ctrl[["points"]];
ord.norta = NorTA.nmds.ctrl[["points"]]; 
xrange = c(min(ord.real[, 1],ord.wgan[, 1], ord.norta[, 1] ) - 0.4, 
           max(ord.real[, 1],ord.wgan[, 1], ord.norta[, 1] ) + 0.4)
yrange = c(min(ord.real[, 2],ord.wgan[, 2], ord.norta[, 2] ) - 0.4, 
           max(ord.real[, 2],ord.wgan[, 2], ord.norta[, 2] ) + 0.4)
plot(1, 1, type="n", main=" ",
     xlim = xrange, ylim = yrange, 
     xlab = "nMDS1", ylab = "nMDS2")
points(NorTA.nmds.ctrl$points[, 1], NorTA.nmds.ctrl$points[, 2], pch=20, col = alpha(rgb(0,1,0,0.5), 0.1))
points(real.nmds.ctrl$points[, 1], real.nmds.ctrl$points[, 2], pch=20, col = "blue")
points(wgan.nmds.ctrl$points[, 1], wgan.nmds.ctrl$points[, 2] , pch=20, col = alpha(rgb(1,0,0,0.5), 0.3))
# add legend
legend("bottomright", 
       legend=c("Real data","MB-GAN","NorTA"), 
       col=c("blue",rgb(1,0,0,0.5),rgb(0,1,0,0.5)), pch=20,
       cex = 0.8,bty = "n")

