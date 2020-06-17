# 50-fold MI: run on computing cluster

library(mice)
#setwd("...")
load("data_incl_censor.RData")

set.seed(2333)
dat2.imp=dat2[!is.na(dat2$bin3yr),-c(1:2)]
imp=mice(dat2.imp, method='cart', printFlag=FALSE, m=50, seed=2333)
imp.dat=complete(imp, action="long")
save.image("imputed_50fold_assess6_incl_censor.RData")