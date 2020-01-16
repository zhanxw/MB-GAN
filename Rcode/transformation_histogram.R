source("utility.R")
library(ggplot2)
library(scales)
library(cowplot)
theme_set(theme_cowplot())

setwd("../data/")
# load the real data #
load("Data_NielsenHB_real.Rdata", verbose = T)

# case group #
abd.raw = case.real.spec
abd.trans = log((1 + 1000*abd.raw) / (1 + abd.raw))
abd.trans[is.na(abd.trans)] = 0

abd.raw.v = apply(abd.raw, 2, mean)
abd.trans.v = apply(abd.trans, 2, mean)

case.abd.df = data.frame(Abundance = c(abd.raw.v, abd.trans.v), 
                         Type = rep(c("Original", "Transformed"), 
                                    each = length(abd.raw.v)) )


# 
ggplot(subset(case.abd.df, Type == "Original"), aes(x=Abundance)) + 
  geom_histogram( bins = 30)+ylim(0, 620)+
  xlab("Original abundance")

ggplot(subset(case.abd.df, Type == "Transformed"), aes(x=Abundance)) + 
  geom_histogram( bins = 30)+ylim(0, 620)+
  xlab("Transformed abundance")

# get the transformation function #
x = seq(0, 1, by = 0.0001)
x.trans = log(1 + 1000*x ) / log(1 + x )
x.trans = log((1 + 1000*x ) / (1+x))
par(mar = c(5, 7, 2, 2))
plot(x, x.trans, type = "l",
     ylab = expression(paste("log ", frac("1+1000x", "1+x"))), 
     main = "Abundance transformation")




