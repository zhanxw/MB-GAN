##### second-order comparison #####
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

vis.spearman.mat = function(Xmat1, 
                            Xmat2,
                            shape = "ellipse",
                            threshold = 0.1, lim.num = 50){
  if(any(colnames(data.matrix(Xmat1)) != colnames(data.matrix(Xmat2))) ){
    stop("The column names of 2 matrices should be identical.")
  }else if(any(data.matrix(Xmat1) > 1) | any(data.matrix(Xmat2) < 0) ){
    stop("True matrix should have values ranging from 0 to 1.")
  }
  # remove the taxa with all zeros across all the samples #
  nreal = nrow(Xmat1); nsimu = nrow(Xmat2)
  zero.real = apply(Xmat1, 2, function(x){sum(x == 0)}) == nreal
  zero.simu = apply(Xmat2, 2, function(x){sum(x == 0)}) == nsimu
  
  # extract the index of taxa to be considered next #
  zero.rm.idx = zero.real | zero.simu
  zero.count.by.taxa = apply(Xmat1, 2, function(x){sum(x == 0)})
  
  extr.idx = (zero.count.by.taxa < nreal * threshold) & (!zero.rm.idx)
  
  Y.real.extr = Xmat1[, extr.idx]
  Y.simu.extr = Xmat2[, extr.idx]
  cor.simu = cor(Y.simu.extr, method = "spearman")
  cor.real = cor(Y.real.extr, method = "spearman")
  
  # compare the correlation structure #
  cor.simu.plot = cor.simu
  cor.true.plot = cor.real
  
  if(sum(extr.idx) > lim.num){
    stop("Too many species selected. Try a smaller threshould.")
  }else{
    
    corrplot(cor.simu.plot, type = "upper", diag = F, method = shape,
             mar = c(1, 0, 1, 0),tl.pos='n')
    corrplot(cor.true.plot, type = "upper", diag = F, method = shape,
             mar = c(1, 0, 1, 0),tl.pos='n')
    
  }
  
  # a spearman correlation test #
  cor.simu.v = as.vector(cor.simu[upper.tri(cor.simu)])
  cor.real.v = as.vector(cor.real[upper.tri(cor.real)])
  
  cor.res = cor(cor.simu.v, cor.real.v, method = c("spearman"))
  cat("Spearman-Correlation between the two correlation matrices is ", round(cor.res, 2))
  dist.test = ks.test(cor.simu.v, cor.real.v)
  cat("\n Test of distribution: p-value = ", dist.test$p.value)
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
                                     df1, df2, df3, 
                                     threshold1, threshold2, y.hi = 1){
  # subset the data #
  check.idx = which(sparsity.real <= threshold2 & sparsity.real>=threshold1)
  tmp1 = as.matrix(df1[, check.idx])
  tmp2 = as.matrix(df2[, check.idx])
  tmp3 = as.matrix(df3[, check.idx])
  # t-test for p-value 
  my_comparisons = list( c("Real", "MB-GAN"), c("Real", "NorTA") )
  ctrl.abd.df = data.frame(Type = rep(c("Real", "MB-GAN", "NorTA"),
                                      c(prod(dim(tmp1)), prod(dim(tmp2)),
                                        prod(dim(tmp3)))),
                           Abundance = c(as.vector(tmp1), as.vector(tmp2),
                                         as.vector(tmp3)), stringsAsFactors = F)
  ctrl.abd.df$Type = factor(ctrl.abd.df$Type, levels = c("Real", "MB-GAN", "NorTA")) 
  # generate plot #
  if(is.null(y.hi)){
    ctrl.abd = ggplot(ctrl.abd.df, aes(x=Type, y=Abundance, color = Type)) + 
      geom_boxplot(outlier.size = 0.5) + 
      xlab(" ") +
      scale_color_manual(values=c("blue", rgb(1,0,0,0.5), rgb(0,1,0,0.5)))+
      theme(legend.position="none",
            axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))+
      stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "wilcox.test")
  }else{
    ctrl.abd = ggplot(ctrl.abd.df, aes(x=Type, y=Abundance, color = Type)) + 
      geom_boxplot(outlier.size = 0.5) + 
      xlab(" ") + ylim(0, y.hi)+
      scale_color_manual(values=c("blue", rgb(1,0,0,0.5), rgb(0,1,0,0.5)))+
      theme(legend.position="none",
            axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))+
      stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "wilcox.test")
  }
  
  plot(ctrl.abd)
}





