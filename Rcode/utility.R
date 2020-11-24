clr <- function(x.f, base=exp(1), tol=.Machine$double.eps, ...) {
  nzero <- (x.f >= tol)
  LOG <- log(ifelse(nzero, x.f, 1), base)
  ifelse(nzero, LOG - mean(LOG)/mean(nzero), 0.0)
}

##### second-order comparison #####


##### spearman correlation ####

vis.spearman = function(Xmat1, Xmat2, Xmat3, Xmat4, taxa.name, plot.hist = FALSE, nbreak = 20){
  Xmat1 = as.matrix(Xmat1); 
  Xmat2 = as.matrix(Xmat2); 
  Xmat3 = as.matrix(Xmat3); 
  Xmat4 = as.matrix(Xmat4)
  
  # 1. extract the taxa for plot #
  Xmat1.plt = Xmat1[ , colnames(Xmat1) %in% taxa.name]
  Xmat2.plt = Xmat2[ , colnames(Xmat2) %in% taxa.name]
  Xmat3.plt = Xmat3[ , colnames(Xmat3) %in% taxa.name]
  Xmat4.plt = Xmat4[ , colnames(Xmat4) %in% taxa.name]
  
  # 2. compare the correlation structure #
  phi1 = cor(Xmat1.plt)
  phi2 = cor(Xmat2.plt)
  phi3 = cor(Xmat3.plt)
  phi4 = cor(Xmat4.plt)
  col.range = range(c(phi1, phi2, phi3, phi4 ))
  
  # correlation plot #
  par(mfcol = c(1, 4))
  corrplot(phi1, type = "upper", diag = F, method = "ellipse",is.corr = T,
           main = "Real Data", col=colorRampPalette(c("blue","white","red"))(200),
           mar = c(0, 0, 1, 0),tl.pos='n', cl.lim = c(min(col.range), max(col.range)))
  corrplot(phi2, type = "upper", diag = F, method = "ellipse",is.corr = T,
           main = "MB-GAN", col=colorRampPalette(c("blue","white","red"))(200),
           mar = c(0, 0, 1, 0),tl.pos='n', cl.lim = c(min(col.range), max(col.range)))
  corrplot(phi3, type = "upper", diag = F, method = "ellipse",is.corr = T,
           main = "Normal-To-Anything", col=colorRampPalette(c("blue","white","red"))(200),
           mar = c(0, 0, 1, 0),tl.pos='n', cl.lim = c(min(col.range), max(col.range)))
  corrplot(phi4, type = "upper", diag = F, method = "ellipse",is.corr = T,
           main = "metaSPARSim", col=colorRampPalette(c("blue","white","red"))(200),
           mar = c(0, 0, 1, 0),tl.pos='n', cl.lim = c(min(col.range), max(col.range)))
  par(mfrow = c(1, 1)) 
  
  
  corrplot(phi1, type = "upper", diag = F, method = "ellipse",is.corr = T,
           main = " ", tl.cex = 0.5,tl.pos = "td",col=colorRampPalette(c("blue","white","red"))(200),
           mar = c(1, 0, 0, 6),cl.lim = c(min(col.range), max(col.range)))
  # 
  #xrange = c(0, max(phi1.v, phi2.v, phi3.v) + 0.05)
  if(plot.hist){
    # histogram of th entries #
    phi1.v = phi1[upper.tri(phi1)]
    phi2.v = phi2[upper.tri(phi2)]
    phi3.v = phi3[upper.tri(phi3)]
    phi4.v = phi4[upper.tri(phi4)]
    
    xrange = range(phi1.v, phi2.v, phi3.v, phi4.v)
    # hist 2 is for the real data 
    hist2 = hist(phi1.v, breaks = nbreak, plot=FALSE )
    hist1 = hist(phi2.v, breaks = hist2$breaks, plot=FALSE)
    hist3 = hist(phi3.v, breaks = hist2$breaks, plot=FALSE)
    hist4 = hist(phi4.v, breaks = hist2$breaks, plot=FALSE)
    
    yhist1 = max(hist1$density)
    yhist2 = max(hist2$density)
    yhist3 = max(hist3$density)
    yhist4 = max(hist4$density)
    
    par(mar = c(4, 5, 2, 1))
    par(mfrow = c(1, 3))
    plot(hist2, xlim = xrange, 
         col="blue", freq = F, main = "Real Data vs MB-GAN",
         ylim = c(0, max(yhist1, yhist2, yhist3)+1),
         xlab = "Spearman Correlation", 
         ylab = "Density", cex.lab=2, cex.main=2)
    plot(hist1, col=rgb(1,0,0,0.5),xlim = xrange, add = T, freq = F)
    
    plot(hist2, xlim = xrange, col="blue", freq = F, main = "Real Data vs NorTA",
         ylim = c(0, max(yhist1, yhist2, yhist3)+1),
         xlab = "Spearman Correlation", 
         ylab = "Density", cex.lab=2, cex.main=2)
    plot(hist3, col=rgb(0,0.5,0,0.5),xlim = xrange, add = T, freq = F)
    
    plot(hist2, xlim = xrange, col="blue", freq = F, main = "Real Data vs metaSPARSim",
         ylim = c(0, max(yhist1, yhist2, yhist3)+1),
         xlab = "Spearman Correlation", 
         ylab = "Density", cex.lab=2, cex.main=2)
    plot(hist4, col=rgb(0.12, 0.56, 1,0.5),xlim = xrange, add = T, freq = F)
    par(mfrow = c(1, 1))
    # ks.test #
    ks1 = ks.test(phi1.v, phi2.v)
    ks2 = ks.test(phi1.v, phi3.v)
    ks3 = ks.test(phi1.v, phi4.v)
    
    pval.1 = ks1$p.value
    pval.2 = ks2$p.value
    pval.3 = ks3$p.value
    cat("mbgan pval = ", pval.1, "\n norta p-val = ", pval.2, "\n metap p-val = ", pval.3)
  }
  
}

compare.spearman.hist = function(Xmat1, Xmat2, Xmat3, 
                                 threshold = 0.1, nbreak = 10){
  if(any(colnames(data.matrix(Xmat1)) != colnames(data.matrix(Xmat2))) ){
    stop("The column names of 2 matrices should be identical.")
  }else if(any(data.matrix(Xmat1) > 1) | any(data.matrix(Xmat2) < 0) ){
    stop("True matrix should have values ranging from 0 to 1.")
  }
  # remove the taxa with all zeros across all the samples #
  nreal = nrow(Xmat1); nsimu = nrow(Xmat2); 
  zero.real = apply(Xmat1, 2, function(x){sum(x == 0)}) == nreal
  zero.simu = apply(Xmat2, 2, function(x){sum(x == 0)}) == nsimu
  zero.simu2 = apply(Xmat3, 2, function(x){sum(x == 0)}) == nsimu
  
  
  # extract the index of taxa to be considered next #
  zero.rm.idx = zero.real | zero.simu |zero.simu2
  zero.count.by.taxa = apply(Xmat1, 2, function(x){sum(x == 0)})
  
  extr.idx = (zero.count.by.taxa < nreal * threshold) & (!zero.rm.idx)
  
  Y.real.extr = Xmat1[, extr.idx]
  Y.simu.extr = Xmat2[, extr.idx]
  Y.simu.extr2 = Xmat3[, extr.idx]
  
  cor.simu = cor(Y.simu.extr, method = "spearman")
  cor.simu2 = cor(Y.simu.extr2, method = "spearman")
  cor.real = cor(Y.real.extr, method = "spearman")
  
  cor.simu.upper = cor.simu[upper.tri(cor.simu)]
  cor.simu.upper2 = cor.simu2[upper.tri(cor.simu2)]
  cor.real.upper = cor.real[upper.tri(cor.real)]
  
  cor.simu.v = as.vector(cor.simu.upper)
  cor.simu.v2 = as.vector(cor.simu.upper2)
  cor.real.v = as.vector(cor.real.upper)
  
  xrange = c(min(cor.simu.v,cor.simu.v2, cor.real.v)-0.05, max(cor.simu.v,cor.simu.v2, cor.real.v)+0.05)
  #par(mfrow = c(2, 1))
  hist2 = hist(cor.real.v, breaks = nbreak, plot=FALSE)
  
  hist1 = hist(cor.simu.v, breaks = nbreak, plot=FALSE)
  #lines1 = lines(density(cor.simu.v))
  yhist1 = max(hist1$density)
  
  #lines2 = lines(density(cor.real.v))
  yhist2 = max(hist2$density)
  hist3 = hist(cor.simu.v2, breaks = nbreak, plot=FALSE)
  yhist3 = max(hist3$density)
  
  plot(hist2, xlim = xrange, 
       col="blue", freq = F, main = " ",
       ylim = c(0, max(yhist1, yhist2, yhist3)),
       xlab = "Spearman Correlation", 
       ylab = "Density")
  plot(hist1, col=rgb(1,0,0,0.5),xlim = xrange, add = T, freq = F)
  
  plot(hist2, xlim = xrange, col="blue", freq = F, main = " ",
       ylim = c(0, max(yhist1, yhist2, yhist3)),
       xlab = "Spearman Correlation", 
       ylab = "Density")
  plot(hist3, col=rgb(0,1,0,0.5),xlim = xrange, add = T, freq = F)
  
  #lines(density(cor.simu.v), col="blue")
  #lines(lines(density(cor.real.v)),col="red")
  #par(mfrow = c(1, 1))
  
  
}


scatter.spearman = function(Xmat1, Xmat2, Xmat3, Xmat4, 
                            taxa.name,plot.hist = F,
                            nbreak = 10){
  
  Xmat1.sub = Xmat1[, colnames(Xmat1) %in% taxa.name]
  Xmat2.sub = Xmat2[, colnames(Xmat2) %in% taxa.name]
  Xmat3.sub = Xmat3[, colnames(Xmat3) %in% taxa.name]
  Xmat4.sub = Xmat4[, colnames(Xmat4) %in% taxa.name]
  
  cor.simu = cor(Xmat2.sub, method = "spearman")
  cor.simu2 = cor(Xmat3.sub, method = "spearman")
  cor.simu3 = cor(Xmat4.sub, method = "spearman")
  cor.real = cor(Xmat1.sub, method = "spearman")
  
  cor.simu.upper = cor.simu[upper.tri(cor.simu)]
  cor.simu.upper2 = cor.simu2[upper.tri(cor.simu2)]
  cor.simu.upper3 = cor.simu3[upper.tri(cor.simu3)]
  cor.real.upper = cor.real[upper.tri(cor.real)]
  
  cor.simu.v = as.vector(cor.simu.upper)
  cor.simu.v2 = as.vector(cor.simu.upper2)
  cor.simu.v3 = as.vector(cor.simu.upper3)
  cor.real.v = as.vector(cor.real.upper)
  
  # scatter plot #
  # compare the correlation structure #
  phi1 = cor.real.v
  phi2 = cor.simu.v
  phi3 = cor.simu.v2
  phi4 = cor.simu.v3
  col.range = range(c(phi1, phi2, phi3, phi4 ))
  
  lm1 = lm(as.vector(phi1)~as.vector(phi2))
  lm1.info = summary(lm1) 
  
  par(mar = c(4, 6, 3, 2))
  par(pty="s")
  par(mfrow = c(1, 3))
  plot(phi1, phi2,xlim = col.range, ylim = col.range, pch = 16,
       xlab = "Real data",
       ylab = "MB-GAN simulated data",
       main = sprintf("R\UB2 = %3.4f \nMSE = %3.4f", 
                      lm1.info$adj.r.squared, 
                      mean(lm1.info$residuals^2)),
       cex.lab=1.5, cex.main=1.5)
  
  lm2 = lm(as.vector(phi1)~as.vector(phi3))
  lm2.info = summary(lm2) 
  lm2.info$adj.r.squared
  
  # scatter plot #
  plot(phi1, phi3,xlim = col.range, ylim = col.range, pch = 16,
       xlab = "Real data",
       ylab = "NorTA simulated data",
       main = sprintf("R\UB2 = %3.4f \nMSE = %3.4f", 
                      lm2.info$adj.r.squared, 
                      mean(lm2.info$residuals^2)),
       cex.lab=1.5, cex.main=1.5)
  
  lm3 = lm(as.vector(phi1)~as.vector(phi4))
  lm3.info = summary(lm3) 
  lm3.info$adj.r.squared
  
  # scatter plot #
  plot(phi1, phi4,xlim = col.range, ylim = col.range, pch = 16,
       xlab = "Real data",
       ylab = "metaSPARSim simulated data",
       cex.lab=2, cex.main=2,
       main = sprintf("R\UB2 = %3.4f \nMSE = %3.4f", 
                      lm3.info$adj.r.squared, 
                      mean(lm3.info$residuals^2)),
       cex.lab=1.5, cex.main=1.5)
  
  par(mfrow = c(1, 1))
  if(plot.hist){
    # histogram #
    xrange = c(min(cor.simu.v,cor.simu.v2,cor.simu.v3, cor.real.v)-0.05, max(cor.simu.v,cor.simu.v2,cor.simu.v3,cor.real.v)+0.05)
    
    hist2 = hist(cor.real.v, breaks = nbreak, plot=FALSE)
    hist1 = hist(cor.simu.v, breaks = hist2$breaks, plot=FALSE)
    hist3 = hist(cor.simu.v2, breaks = hist2$breaks, plot=FALSE)
    hist4 = hist(cor.simu.v3, breaks = hist2$breaks, plot=FALSE)
    yhist1 = max(hist1$density)
    yhist2 = max(hist2$density)
    yhist3 = max(hist3$density)
    yhist4 = max(hist4$density)
    
    plot(hist2, xlim = xrange, 
         col="blue", freq = F, main = " ",
         ylim = c(0, max(yhist1, yhist2, yhist3)),
         xlab = "Spearman Correlation", 
         ylab = "Density")
    plot(hist1, col=rgb(1,0,0,0.5),xlim = xrange, add = T, freq = F)
    
    plot(hist2, xlim = xrange, col="blue", freq = F, main = " ",
         ylim = c(0, max(yhist1, yhist2, yhist3)),
         xlab = "Spearman Correlation", 
         ylab = "Density")
    plot(hist3, col=rgb(0,1,0,0.5),xlim = xrange, add = T, freq = F)
    
    plot(hist2, xlim = xrange, col="blue", freq = F, main = " ",
         ylim = c(0, max(yhist1, yhist2, yhist3)),
         xlab = "Spearman Correlation", 
         ylab = "Density")
    plot(hist4, col=rgb(1,1,0,0.5),xlim = xrange, add = T, freq = F)
  }
}



scatter.spearman.2mat = function(Xmat1, Xmat2,
                                 taxa.name,plot.hist = F,
                                 nbreak = 10){
  
  Xmat1.sub = Xmat1[, colnames(Xmat1) %in% taxa.name]
  Xmat2.sub = Xmat2[, colnames(Xmat2) %in% taxa.name]
  
  cor.simu = cor(Xmat2.sub, method = "spearman")
  cor.real = cor(Xmat1.sub, method = "spearman")
  
  cor.simu.upper = cor.simu[upper.tri(cor.simu)]
  cor.real.upper = cor.real[upper.tri(cor.real)]
  
  cor.simu.v = as.vector(cor.simu.upper)
  cor.real.v = as.vector(cor.real.upper)
  
  # scatter plot #
  # compare the correlation structure #
  phi1 = cor.real.v
  phi2 = cor.simu.v
  col.range = range(c(phi1, phi2))
  
  lm1 = lm(as.vector(phi1)~as.vector(phi2))
  lm1.info = summary(lm1) 
  
  par(pty="s")
  plot(phi1, phi2,xlim = col.range, ylim = col.range, pch = 16,
       xlab = "Spearman correlation of the real data",
       ylab = "Spearman correlation of the \n perturbed data")
  legend("topleft", bty="n", 
         legend=sprintf("R\UB2 = %3.4f \nMSE = %3.4f", 
                        lm1.info$adj.r.squared, 
                        mean(lm1.info$residuals^2))
  )
  
  if(plot.hist){
    # histogram #
    xrange = c(min(cor.simu.v,cor.real.v)-0.05, max(cor.simu.v,cor.real.v)+0.05)
    
    hist2 = hist(cor.real.v, breaks = nbreak, plot=FALSE)
    hist1 = hist(cor.simu.v, breaks = hist2$breaks, plot=FALSE)
    yhist1 = max(hist1$density)
    yhist2 = max(hist2$density)
    
    plot(hist2, xlim = xrange, 
         col="blue", freq = F, main = " ",
         ylim = c(0, max(yhist1, yhist2)),
         xlab = "Spearman Correlation", 
         ylab = "Density")
    plot(hist1, col=rgb(1,0,0,0.5),xlim = xrange, add = T, freq = F)
  }
}


#### proportionality ####

vis.proportionality = function(Xmat1, Xmat2, Xmat3, Xmat4, taxa.name, plot.hist = FALSE, nbreak = 20){
  require(compositions, quietly = T)
  
  Xmat1 = as.matrix(Xmat1); Xmat2 = as.matrix(Xmat2); Xmat3 = as.matrix(Xmat3); Xmat4 = as.matrix(Xmat4)
  # get the clr transformation #
  # 1. impute the 0s in the compositional by the second largest number #
  impute.val1 = unique(sort(as.vector(Xmat1)))[2]
  Xmat1[Xmat1 == 0] = impute.val1
  
  impute.val2 = unique(sort(as.vector(Xmat2)))[2]
  Xmat2[Xmat2 == 0] = impute.val2
  
  impute.val3 = unique(sort(as.vector(Xmat3)))[2]
  Xmat3[Xmat3 == 0] = impute.val3
  
  impute.val4 = unique(sort(as.vector(Xmat4)))[2]
  Xmat4[Xmat4 == 0] = impute.val4
  # 2. compute the variations #
  class(Xmat1) = 'acomp'
  Xmat1.vlr = variation(Xmat1)
  
  class(Xmat2) = 'acomp'
  Xmat2.vlr = variation(Xmat2)
  
  class(Xmat3) = 'acomp'
  Xmat3.vlr = variation(Xmat3)
  
  class(Xmat4) = 'acomp'
  Xmat4.vlr = variation(Xmat4)
  
  # 3. apply the clr transformation #
  Xmat1.clr = t(apply(Xmat1, 1, clr))
  Xmat2.clr = t(apply(Xmat2, 1, clr))
  Xmat3.clr = t(apply(Xmat3, 1, clr))
  Xmat4.clr = t(apply(Xmat4, 1, clr))
  
  # 4. compute the variance of clr transformated matrices #
  Xmat1.var = apply(Xmat1.clr, 2, var)
  Xmat2.var = apply(Xmat2.clr, 2, var)
  Xmat3.var = apply(Xmat3.clr, 2, var)
  Xmat4.var = apply(Xmat4.clr, 2, var)
  
  # 5. compute phi #
  Xmat1.phi = sweep(Xmat1.vlr, 2, Xmat1.var, FUN="/")
  Xmat2.phi = sweep(Xmat2.vlr, 2, Xmat2.var, FUN="/")
  Xmat3.phi = sweep(Xmat3.vlr, 2, Xmat3.var, FUN="/")
  Xmat4.phi = sweep(Xmat4.vlr, 2, Xmat4.var, FUN="/")
  
  # 6. extract the taxa for plot #
  Xmat1.phi.sub = Xmat1.phi[rownames(Xmat1.phi) %in% taxa.name, colnames(Xmat1.phi) %in% taxa.name]
  Xmat2.phi.sub = Xmat2.phi[rownames(Xmat2.phi) %in% taxa.name, colnames(Xmat2.phi) %in% taxa.name]
  Xmat3.phi.sub = Xmat3.phi[rownames(Xmat3.phi) %in% taxa.name, colnames(Xmat3.phi) %in% taxa.name]
  Xmat4.phi.sub = Xmat4.phi[rownames(Xmat4.phi) %in% taxa.name, colnames(Xmat4.phi) %in% taxa.name]
  
  # 7. compare the correlation structure #
  phi1 = Xmat1.phi.sub
  phi2 = Xmat2.phi.sub
  phi3 = Xmat3.phi.sub
  phi4 = Xmat4.phi.sub
  col.range = range(c(phi1, phi2, phi3, phi4 ))
  
  # correlation plot #
  par(mfcol = c(2, 2))
  corrplot(phi1, type = "full", diag = F, method = "color",is.corr = F,
           main = "Real Data", 
           mar = c(1, 0, 1, 0),tl.pos='n', cl.lim = c(min(col.range), max(col.range)),
           col= colorRampPalette(c("red", "white"))(30))
  corrplot(phi2, type = "full", diag = F, method = "color",is.corr = F,
           main = "MB-GAN", 
           mar = c(1, 0, 1, 0),tl.pos='n', cl.lim = c(min(col.range), max(col.range)),
           col= colorRampPalette(c("red", "white"))(30))
  corrplot(phi3, type = "full", diag = F, method = "color",is.corr = F,
           main = "Normal-To-Anything", 
           mar = c(1, 0, 1, 0),tl.pos='n', cl.lim = c(min(col.range), max(col.range)),
           col= colorRampPalette(c("red", "white"))(30))
  corrplot(phi4, type = "full", diag = F, method = "color",is.corr = F,
           main = "metaSPARSim", 
           mar = c(1, 0, 1, 0),tl.pos='n', cl.lim = c(min(col.range), max(col.range)),
           col= colorRampPalette(c("red", "white"))(30))
  par(mfrow = c(1, 1))
  # histogram of th entries #
  phi1.v = as.vector(phi1[phi1!=0])
  phi2.v = as.vector(phi2[phi2!=0])
  phi3.v = as.vector(phi3[phi3!=0])
  phi4.v = as.vector(phi4[phi4!=0])
  
  #xrange = c(0, max(phi1.v, phi2.v, phi3.v) + 0.05)
  if(plot.hist){
    xrange = range(phi1.v, phi2.v, phi3.v, phi4.v)
    # hist 2 is for the real data 
    hist2 = hist(phi1.v, breaks = nbreak, plot=FALSE )
    hist1 = hist(phi2.v, breaks = hist2$breaks, plot=FALSE)
    hist3 = hist(phi3.v, breaks = hist2$breaks, plot=FALSE)
    hist4 = hist(phi4.v, breaks = hist2$breaks, plot=FALSE)
    
    yhist1 = max(hist1$density)
    yhist2 = max(hist2$density)
    yhist3 = max(hist3$density)
    yhist4 = max(hist4$density)
    
    plot(hist2, xlim = xrange, 
         col="blue", freq = F, main = "Real Data vs MB-GAN",
         ylim = c(0, max(yhist1, yhist2, yhist3)+0.1),
         xlab = "Proportionality", 
         ylab = "Density")
    plot(hist1, col=rgb(1,0,0,0.5),xlim = xrange, add = T, freq = F)
    
    plot(hist2, xlim = xrange, col="blue", freq = F, main = "Real Data vs NorTA",
         ylim = c(0, max(yhist1, yhist2, yhist3)+0.1),
         xlab = "Proportionality", 
         ylab = "Density")
    plot(hist3, col=rgb(0,0.5,0,0.5),xlim = xrange, add = T, freq = F)
    
    plot(hist2, xlim = xrange, col="blue", freq = F, main = "Real Data vs metaSPARSim",
         ylim = c(0, max(yhist1, yhist2, yhist3)+0.1),
         xlab = "Proportionality", 
         ylab = "Density")
    plot(hist4, col=rgb(0.12, 0.56, 1,0.5),xlim = xrange, add = T, freq = F)
    
    
  }
  
}


scatter.proportionality = function(Xmat1, Xmat2, Xmat3,Xmat4, taxa.name){
  require(compositions, quietly = T)
  
  Xmat1 = as.matrix(Xmat1); Xmat2 = as.matrix(Xmat2); Xmat3 = as.matrix(Xmat3);Xmat4 = as.matrix(Xmat4)
  # get the clr transformation #
  # 1. impute the 0s in the compositional by the second largest number #
  impute.val1 = unique(sort(as.vector(Xmat1)))[2]
  Xmat1[Xmat1 == 0] = impute.val1
  
  impute.val2 = unique(sort(as.vector(Xmat2)))[2]
  Xmat2[Xmat2 == 0] = impute.val2
  
  impute.val3 = unique(sort(as.vector(Xmat3)))[2]
  Xmat3[Xmat3 == 0] = impute.val3
  
  impute.val4 = unique(sort(as.vector(Xmat4)))[2]
  Xmat4[Xmat4 == 0] = impute.val4
  
  # 2. compute the variations #
  class(Xmat1) = 'acomp'
  Xmat1.vlr = variation(Xmat1)
  
  class(Xmat2) = 'acomp'
  Xmat2.vlr = variation(Xmat2)
  
  class(Xmat3) = 'acomp'
  Xmat3.vlr = variation(Xmat3)
  
  class(Xmat4) = 'acomp'
  Xmat4.vlr = variation(Xmat4)
  
  # 3. apply the clr transformation #
  Xmat1.clr = t(apply(Xmat1, 1, clr))
  Xmat2.clr = t(apply(Xmat2, 1, clr))
  Xmat3.clr = t(apply(Xmat3, 1, clr))
  Xmat4.clr = t(apply(Xmat4, 1, clr))
  
  # 4. compute the variance of clr transformated matrices #
  Xmat1.var = apply(Xmat1.clr, 2, var)
  Xmat2.var = apply(Xmat2.clr, 2, var)
  Xmat3.var = apply(Xmat3.clr, 2, var)
  Xmat4.var = apply(Xmat4.clr, 2, var)
  
  # 5. compute phi #
  Xmat1.phi = sweep(Xmat1.vlr, 2, Xmat1.var, FUN="/")
  Xmat2.phi = sweep(Xmat2.vlr, 2, Xmat2.var, FUN="/")
  Xmat3.phi = sweep(Xmat3.vlr, 2, Xmat3.var, FUN="/")
  Xmat4.phi = sweep(Xmat4.vlr, 2, Xmat4.var, FUN="/")
  
  # 6. extract the taxa for plot #
  Xmat1.phi.sub = Xmat1.phi[rownames(Xmat1.phi) %in% taxa.name, colnames(Xmat1.phi) %in% taxa.name]
  Xmat2.phi.sub = Xmat2.phi[rownames(Xmat2.phi) %in% taxa.name, colnames(Xmat2.phi) %in% taxa.name]
  Xmat3.phi.sub = Xmat3.phi[rownames(Xmat3.phi) %in% taxa.name, colnames(Xmat3.phi) %in% taxa.name]
  Xmat4.phi.sub = Xmat4.phi[rownames(Xmat4.phi) %in% taxa.name, colnames(Xmat4.phi) %in% taxa.name]
  
  # 7. compare the correlation structure #
  phi1 = Xmat1.phi.sub
  phi2 = Xmat2.phi.sub
  phi3 = Xmat3.phi.sub
  phi4 = Xmat4.phi.sub
  col.range = range(c(phi1, phi2, phi3,phi4 ))
  
  lm1 = lm(as.vector(phi1)~as.vector(phi2))
  lm1.info = summary(lm1) 
  
  # scatter plot #
  par(mar = c(4, 6, 3, 2))
  par(pty="s")
  par(mfrow = c(1, 3))
  plot(phi1, phi2,xlim = col.range, ylim = col.range, pch = 16,
       xlab = "Real data",
       ylab = "MB-GAN simulated data", 
       main = sprintf("R\UB2 = %3.4f \nMSE = %3.4f", 
                      lm1.info$adj.r.squared, 
                      mean(lm1.info$residuals^2)),
       cex.lab=1.5, cex.main=1.5)
  
  lm2 = lm(as.vector(phi1)~as.vector(phi3))
  lm2.info = summary(lm2) 
  
  # scatter plot #
  plot(phi1, phi3,xlim = col.range, ylim = col.range, pch = 16,
       xlab = "Real data",
       ylab = "NorTA simulated data", 
       main = sprintf("R\UB2 = %3.4f \nMSE = %3.4f", 
                      lm2.info$adj.r.squared, 
                      mean(lm2.info$residuals^2)),
       cex.lab=1.5, cex.main=1.5)
  
  
  lm3 = lm(as.vector(phi1)~as.vector(phi4))
  lm3.info = summary(lm3) 
  
  # scatter plot #
  plot(phi1, phi4,xlim = col.range, ylim = col.range, pch = 16,
       xlab = "Real data",
       ylab = "metaSPARSim simulated data", 
       main = sprintf("R\UB2 = %3.4f \nMSE = %3.4f", 
                      lm3.info$adj.r.squared, 
                      mean(lm3.info$residuals^2)),
       cex.lab=1.5, cex.main=1.5)
  par(mfrow = c(1, 1))
}



vis.proportionality.2mat = function(Xmat1, Xmat2, taxa.name, nbreak = 20){
  require(compositions, quietly = T)
  
  Xmat1 = as.matrix(Xmat1); Xmat2 = as.matrix(Xmat2)
  # get the clr transformation #
  # 1. impute the 0s in the compositional by the second largest number #
  impute.val1 = unique(sort(as.vector(Xmat1)))[2]
  Xmat1[Xmat1 == 0] = impute.val1
  
  impute.val2 = unique(sort(as.vector(Xmat2)))[2]
  Xmat2[Xmat2 == 0] = impute.val2
  
  # 2. compute the variations #
  class(Xmat1) = 'acomp'
  Xmat1.vlr = variation(Xmat1)
  
  class(Xmat2) = 'acomp'
  Xmat2.vlr = variation(Xmat2)
  
  # 3. apply the clr transformation #
  Xmat1.clr = t(apply(Xmat1, 1, clr))
  Xmat2.clr = t(apply(Xmat2, 1, clr))
  
  # 4. compute the variance of clr transformated matrices #
  Xmat1.var = apply(Xmat1.clr, 2, var)
  Xmat2.var = apply(Xmat2.clr, 2, var)
  
  # 5. compute phi #
  Xmat1.phi = sweep(Xmat1.vlr, 2, Xmat1.var, FUN="/")
  Xmat2.phi = sweep(Xmat2.vlr, 2, Xmat2.var, FUN="/")
  
  # 6. extract the taxa for plot #
  Xmat1.phi.sub = Xmat1.phi[rownames(Xmat1.phi) %in% taxa.name, colnames(Xmat1.phi) %in% taxa.name]
  Xmat2.phi.sub = Xmat2.phi[rownames(Xmat2.phi) %in% taxa.name, colnames(Xmat2.phi) %in% taxa.name]
  
  # 7. compare the correlation structure #
  phi1 = Xmat1.phi.sub
  phi2 = Xmat2.phi.sub
  col.range = range(c(phi1, phi2))
  
  # scatter plot #
  lm1 = lm(as.vector(phi1)~as.vector(phi2))
  lm1.info = summary(lm1) 
  
  # scatter plot #
  par(pty="s")
  plot(phi1, phi2,xlim = col.range, ylim = col.range, pch = 16,
       xlab = "Proportionality of the real data",
       ylab = "Proportionality of the \n perturbed data")
  legend("topleft", bty="n", 
         legend=sprintf("R\UB2 = %3.4f \nMSE = %3.4f", 
                        lm1.info$adj.r.squared, 
                        mean(lm1.info$residuals^2))
  )
  
  
  # histogram of th entries #
  phi1.v = as.vector(phi1[phi1!=0])
  phi2.v = as.vector(phi2[phi2!=0])
  
  #xrange = c(0, max(phi1.v, phi2.v, phi3.v) + 0.05)
  xrange = c(0,  max(phi2.v) -4)
  # hist 2 is for the real data 
  hist2 = hist(phi1.v, breaks = nbreak, plot=FALSE )
  hist1 = hist(phi2.v, breaks = hist2$breaks, plot=FALSE)
  
  yhist1 = max(hist1$density)
  yhist2 = max(hist2$density)
  
  plot(hist2, xlim = xrange, 
       col="blue", freq = F, main = " ",
       ylim = c(0, max(yhist1, yhist2)+0.1),
       xlab = "Proportionality", 
       ylab = "Density")
  plot(hist1, col=rgb(1,0,0,0.5),xlim = xrange, add = T, freq = F)
  
}


##### first-order comparison #####
get.sparsity = function(x){sum(x == 0)/length(x)}

get.shannon.entropy = function(target) {
  #freq = table(target)/length(target)
  #freq = target
  vec = target
  #drop 0 to avoid NaN resulting from log
  vec = vec[vec>0]
  #compute entropy
  -sum(vec * log(vec))
}

compare.abundance.boxplot = function(sparsity.real, 
                                     df1, df2, df3, df4,
                                     threshold1, threshold2, y.hi = 1){
  # subset the data #
  check.idx = which(sparsity.real <= threshold2 & sparsity.real>=threshold1)
  tmp1 = as.matrix(df1[, check.idx])
  tmp2 = as.matrix(df2[, check.idx])
  tmp3 = as.matrix(df3[, check.idx])
  tmp4 = as.matrix(df4[, check.idx])
  
  # t-test for p-value 
  my_comparisons = list( c("Real", "MB-GAN"), c("Real", "NorTA"),  c("Real", "metaSPARSim"))
  ctrl.abd.df = data.frame(Type = rep(c("Real", "MB-GAN", "NorTA","metaSPARSim"),
                                      c(prod(dim(tmp1)), prod(dim(tmp2)), prod(dim(tmp3)),
                                        prod(dim(tmp4)))),
                           Abundance = c(as.vector(tmp1), as.vector(tmp2),as.vector(tmp3),
                                         as.vector(tmp4)), stringsAsFactors = F)
  ctrl.abd.df$Type = factor(ctrl.abd.df$Type, levels = c("Real", "MB-GAN", "NorTA","metaSPARSim")) 
  # generate plot #
  if(is.null(y.hi)){
    ctrl.abd = ggplot(ctrl.abd.df, aes(x=Type, y=Abundance, color = Type)) + 
      geom_boxplot(outlier.size = 0.5) + 
      xlab(" ") +
      scale_color_manual(values=c("black", rgb(1,0,0,0.5), rgb(0,0.5,0,0.5),  rgb(0.12, 0.56, 1,0.5)))+
      theme(legend.position="none",
            axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))+
      stat_compare_means(comparisons = my_comparisons,  method = "wilcox.test")
    
  }else{
    ctrl.abd = ggplot(ctrl.abd.df, aes(x=Type, y=Abundance, color = Type)) + 
      geom_boxplot(outlier.size = 0.5) + 
      xlab(" ") + ylim(0, y.hi)+
      scale_color_manual(values=c("black", rgb(1,0,0,0.5), rgb(0,0.5,0,0.5),  rgb(0.12, 0.56, 1,0.5)))+
      theme(legend.position="none",
            axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))+
      stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "wilcox.test")
  }
  
  plot(ctrl.abd)
}

#### data preturbation ####
perturb.taxon = function(x, zero.rm = TRUE, sd.ratio = 0.1){
  zero.idx = x == 0
  x.var = ifelse(zero.rm, var(x[!zero.idx]), var(x))
  x.noise = x + abs(rnorm(length(x), mean = 0, sd = sd.ratio* x.var))
  if(zero.rm){
    x.noise[zero.idx] = 0
  }
  return(x.noise)
}

constrain.sample = function(x){
  x.norm = x
  zero.idx = x == 0
  nonzero.x = x[!zero.idx]
  nonzero.x = nonzero.x/sum(nonzero.x)
  x.norm[!zero.idx] = nonzero.x
  return(x.norm)
}


##### first-order comparison #####
get.sparsity = function(x){sum(x == 0)/length(x)}

get.shannon.entropy = function(target) {
  #freq = table(target)/length(target)
  #freq = target
  vec = target
  #drop 0 to avoid NaN resulting from log
  vec = vec[vec>0]
  #compute entropy
  -sum(vec * log(vec))
}

compare.abundance.boxplot = function(sparsity.real, 
                                     df1, df2, df3, df4,
                                     threshold1, threshold2, y.hi = 1){
  # subset the data #
  check.idx = which(sparsity.real <= threshold2 & sparsity.real>=threshold1)
  tmp1 = as.matrix(df1[, check.idx])
  tmp2 = as.matrix(df2[, check.idx])
  tmp3 = as.matrix(df3[, check.idx])
  tmp4 = as.matrix(df4[, check.idx])
  
  # t-test for p-value 
  my_comparisons = list( c("Real", "MB-GAN"), c("Real", "NorTA"),  c("Real", "metaSPARSim"))
  ctrl.abd.df = data.frame(Type = rep(c("Real", "MB-GAN", "NorTA","metaSPARSim"),
                                      c(prod(dim(tmp1)), prod(dim(tmp2)), prod(dim(tmp3)),
                                        prod(dim(tmp4)))),
                           Abundance = c(as.vector(tmp1), as.vector(tmp2),as.vector(tmp3),
                                         as.vector(tmp4)), stringsAsFactors = F)
  ctrl.abd.df$Type = factor(ctrl.abd.df$Type, levels = c("Real", "MB-GAN", "NorTA","metaSPARSim")) 
  # generate plot #
  if(is.null(y.hi)){
    ctrl.abd = ggplot(ctrl.abd.df, aes(x=Type, y=Abundance, color = Type)) + 
      geom_boxplot(outlier.size = 0.5) + 
      xlab(" ") +
      scale_color_manual(values=c("black", rgb(1,0,0,0.5), rgb(0,0.5,0,0.5),  rgb(0.12, 0.56, 1,0.5)))+
      theme(legend.position="none",
            axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))+
      stat_compare_means(comparisons = my_comparisons,  method = "wilcox.test")
    
  }else{
    ctrl.abd = ggplot(ctrl.abd.df, aes(x=Type, y=Abundance, color = Type)) + 
      geom_boxplot(outlier.size = 0.5) + 
      xlab(" ") + ylim(0, y.hi)+
      scale_color_manual(values=c("black", rgb(1,0,0,0.5), rgb(0,0.5,0,0.5),  rgb(0.12, 0.56, 1,0.5)))+
      theme(legend.position="none",
            axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))+
      stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "wilcox.test")
  }
  
  plot(ctrl.abd)
}





DA.test = function(Xmat1, Xmat2, z1, z2, method = c("KW"), sig.level = 0.05, top = 10, subsample = FALSE, seed = 123){
  Xmat1 = as.matrix(Xmat1); Xmat2 = as.matrix(Xmat2)
  if(any(Xmat1 > 1) | any(Xmat2 > 1)){
    stop("Input should be compositional data.")
  }
  if(any(colnames(Xmat1) !=  colnames(Xmat2)) ){
    stop("Column names of 2 input compositional matrices should be the same.")
  }
  # taxa names #
  bug.names = colnames(Xmat1)
  
  # perform KW test #
  pval.from.real = as.numeric(apply(Xmat1, 2, KW.test, z1))
  if(!subsample){
    pval.from.simu = as.numeric(apply(Xmat2, 2, KW.test, z2))
    
    # select the taxa with smallest p-values on the real data #
    adj.pval.from.real = p.adjust(pval.from.real, "BH")
    pval.order = order(adj.pval.from.real)
    pval.top = pval.order[1:top]
    taxa.top.real = bug.names[pval.top]
    real.top.rank = data.frame(taxa = taxa.top.real, pvalue = sort(adj.pval.from.real)[1:top], stringsAsFactors = F)
    
    # compare on the simulated data #
    adj.pval.from.simu = p.adjust(pval.from.simu, "BH")
    pval.simu.rank = rank(adj.pval.from.simu)
    names(pval.simu.rank) = bug.names
    names(adj.pval.from.simu)  = bug.names
    simu.top.rank =  data.frame(taxa = taxa.top.real, ResultRank = as.numeric(pval.simu.rank[taxa.top.real]),
                                pvalue = adj.pval.from.simu[taxa.top.real], stringsAsFactors = F)
    
  }else{
    set.seed(seed)
    ng1 = sum(z1 == 1); ng0 = sum(z1 == 0)
    g1.idx = sample(1:1000, ng1); g0.idx =  sample(1:1000, ng0)
    Xmat2.sub = rbind(Xmat2[g1.idx, ], Xmat2[1000+ g0.idx, ])
    pval.from.simu = as.numeric(apply(Xmat2.sub, 2, KW.test, z1))
  }
  
  # adjuest P-values #
  adj.pval.from.real = p.adjust(pval.from.real, "BH")
  adj.pval.from.simu = p.adjust(pval.from.simu, "BH")
  
  # select DA bugs #
  res.by.real = which(adj.pval.from.real < sig.level)
  sel.real.nam = bug.names[res.by.real]
  res.by.simu = which(adj.pval.from.simu < sig.level)
  sel.simu.nam = bug.names[res.by.simu]
  
  #
  if(!subsample){
    res.list = list(RealRank = real.top.rank,
                    SimuRank = simu.top.rank, 
                    RealDataResult = sel.real.nam, 
                    SimuDataResult = sel.simu.nam, 
                    OverlapResult = intersect(sel.simu.nam, sel.real.nam), 
                    OverlapNumber = length(intersect(sel.simu.nam, sel.real.nam)))
  }else{
    res.list = list(RealDataResult = sel.real.nam, 
                    SimuDataResult = sel.simu.nam, 
                    OverlapResult = intersect(sel.simu.nam, sel.real.nam), 
                    OverlapNumber = length(intersect(sel.simu.nam, sel.real.nam)))
  }
  
  
  return(res.list)
}


KW.test = function(x, z){
  kw.res = kruskal.test(x = x, g = z)
  kw.pval = kw.res$p.value
  return(kw.pval)
}


DA.resample = function(Xmat1, Xmat2, z1, z2, sig.level = 0.05,seed = 123, nresample = 100){
  Xmat1 = as.matrix(Xmat1); Xmat2 = as.matrix(Xmat2)
  if(any(Xmat1 > 1) | any(Xmat2 > 1)){
    stop("Input should be compositional data.")
  }
  if(any(colnames(Xmat1) !=  colnames(Xmat2)) ){
    stop("Column names of 2 input compositional matrices should be the same.")
  }
  # taxa names #
  bug.names = colnames(Xmat1)
  
  # perform KW test #
  ng1 = sum(z1 == 1); ng0 = sum(z1 == 0)
  set.seed(seed)
  resample.list = vector("list", nresample)
  for(ii in 1:nresample){
    if(ii %% 10==0) {
      # Print on the screen some message
      cat(paste0("iteration: ", ii, "\n"))
    }
    g1.idx = sample(1:1000, ng1); g0.idx =  sample(1:1000, ng0)
    Xmat2.sub = rbind(Xmat2[g1.idx, ], Xmat2[1000+ g0.idx, ])
    pval.from.simu = as.numeric(apply(Xmat2.sub, 2, KW.test, z1))
    adj.pval.from.simu = p.adjust(pval.from.simu, "BH")
    res.by.simu = which(adj.pval.from.simu < sig.level)
    sel.simu.nam = bug.names[res.by.simu]
    resample.list[[ii]] = sel.simu.nam
  }
  overlap.num = length(Reduce(intersect, resample.list))
  overlap.nam = Reduce(intersect, resample.list)
  uniq = lapply(resample.list, function(x){length(x) - overlap.num})
  
  # KW test on the real data #
  pval.from.real = as.numeric(apply(Xmat1, 2, KW.test, z1))
  adj.pval.from.real = p.adjust(pval.from.real, "BH")
  res.by.real = which(adj.pval.from.real < sig.level)
  sel.real.nam = bug.names[res.by.real]
  
  inters.with.real = intersect(sel.real.nam, overlap.nam)
  
  res.list = list(RealDataResult = sel.real.nam, 
                  SimuDataIntersect = inters.with.real, 
                  UniqueNumber = uniq)
  return(res.list)
}

KW.boxplot = function(Xmat1, Xmat2, z1, z2, top = 10, seed = 123, nresample = 100){
  require(cowplot)
  theme_set(theme_cowplot())
  require(reshape2)
  Xmat1 = as.matrix(Xmat1); Xmat2 = as.matrix(Xmat2)
  if(any(Xmat1 > 1) | any(Xmat2 > 1)){
    stop("Input should be compositional data.")
  }
  if(any(colnames(Xmat1) !=  colnames(Xmat2)) ){
    stop("Column names of 2 input compositional matrices should be the same.")
  }
  # taxa names #
  bug.names = colnames(Xmat1)
  
  # KW test on the real data #
  pval.from.real = as.numeric(apply(Xmat1, 2, KW.test, z1))
  adj.pval.from.real = p.adjust(pval.from.real, "BH")
  
  # select the taxa with smallest p-values #
  pval.order = order(adj.pval.from.real)
  pval.top = pval.order[1:top]
  taxa.top.real = bug.names[pval.top]
  real.top.rank = data.frame(taxa = taxa.top.real, 
                             pvalue = sort(adj.pval.from.real)[1:top], 
                             stringsAsFactors = F)
  
  
  # perform KW test #
  ng1 = sum(z1 == 1); ng0 = sum(z1 == 0)
  set.seed(seed)
  rank.mat = matrix(0, nrow = nresample, ncol = top)
  for(ii in 1:nresample){
    if(ii %% 10==0) {
      # Print on the screen some message
      cat(paste0("iteration: ", ii, "\n"))
    }
    rand.idx = sample(1:1000, 200); 
    
    Xg1 = rbind(Xmat1[z1 == 0, ], Xmat2[rand.idx, ])
    Xg2 = rbind(Xmat1[z1 == 1, ], Xmat2[rand.idx + 1000, ])
    
    Xmat.new = rbind(Xg1, Xg2)
    z.new = rep(c(0, 1), c(nrow(Xg1), nrow(Xg2)) )
    pval.from.simu = as.numeric(apply(Xmat.new, 2, KW.test, z.new))
    adj.pval.from.simu = p.adjust(pval.from.simu, "BH")
    pval.simu.rank = rank(adj.pval.from.simu)
    names(pval.simu.rank) = bug.names
    
    current.rank = pval.simu.rank[taxa.top.real]
    rank.mat[ii, ] = current.rank
    
  }
  
  colnames(rank.mat) = taxa.top.real
  
  box.df = melt(rank.mat)
  box.df = box.df[, -1]
  colnames(box.df) = c("Taxa", "Rank")
  gplt = ggplot(data = box.df, aes(x = Taxa, y = Rank))+
    geom_boxplot() + theme(legend.position = 'none',axis.text.x = element_text(angle = 45,
                                                                               hjust = 1,
                                                                               size = 9))+
    geom_hline(yintercept=10, linetype = "dashed")+
    ggtitle("")+ylab("Rank of P-values")
  plot(gplt)
  #boxplot(Rank~Taxa, data = box.df,srt = 45, pos = 1)
  return(list(gplt, rank.mat))
}


