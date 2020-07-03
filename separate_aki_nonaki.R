# separating AKI/non-AKI subgroups

# rf and lasso
library(ROCR)
library(mice)
library(randomForest)
library(glmnet)
library(Metrics)
library(gplots)

load("imputed_50fold_assess6_incl_censor.RData")

# create the RF formula using all/subset of predictors
rf.formula="bin3yr~"
n.var=0
# set exclude.idx to obtain results with selected biomarkers (free + top 5)
# set exclude.idx to NULL will include all predictors
exclude.idx=c(64:68,70:71,34,36,38,41:49,52:53,55:61)
#exclude.idx=NULL
for (j in 4:length(names(dat2))){
  if (! (j %in% exclude.idx)){ 
    rf.formula=paste(rf.formula,names(dat2)[j],sep="")
    if (j < length(names(dat2)))
      rf.formula=paste(rf.formula,"+",sep="")
    n.var=n.var+1
  } 
}

set.seed(2333)
MDacc=MDgini=matrix(NA,nrow=n.var,ncol=imp$m) # RF: mean decrease accuracy/Gini for each predictor at each lambda
auc=auc.lasso=rep(NA,imp$m) # RF&lasso: model AUC at each lambda
lam=exp(seq(-8,-2,length.out=60)) # range of lambda to consider
nonzero=matrix(0,nrow=n.var+2,ncol=length(lam)) # Lasso: number of time each predictor remains active
nonzero.optim=rep(0,n.var+2) # Lasso: number of time each predictor remains active at the *optimal* lambda
for (j in 1:imp$m){
  # repeat the procedure for each imputed dataset
  imp.j=imp.dat[imp.dat$.imp==j & imp.dat$V0_AKIN_stage!="non-AKI",-1]
  train.idx=sample(1:nrow(imp.j),round(nrow(imp.j)*0.8),replace=F)
  training=imp.j[train.idx,-1]
  test=imp.j[-train.idx,-1]
  
  #rf
  rf=randomForest(as.formula(rf.formula),data=training,importance=TRUE)
  MDacc[,j]=rf$importance[,3]/rf$importanceSD[,3]
  MDgini[,j]=rf$importance[,4]
  pred.test=predict(rf,test,type="prob")
  pred=prediction(pred.test[,2],test[,1])
  #roc.perf=performance(pred, measure = "tpr", x.measure = "fpr")
  #plot(roc.perf, main="ROC curve")
  #abline(a=0, b= 1)
  auc.perf = performance(pred, measure = "auc")
  auc[j]=auc.perf@y.values
  
  #lasso
  predictors.imp=training[,-c(1,exclude.idx-2)]
  outcome=training$bin3yr
  lasso.fit=glmnet(model.matrix( ~ ., predictors.imp),outcome,family = "binomial",
                   alpha = 1, lambda=lam)
  cv.lasso=cv.glmnet(x=model.matrix( ~ ., predictors.imp), y=outcome, family = "binomial",
                     type.measure="auc", nfolds=5, alpha=1, lambda=lam)# CV to select the optimal lambda
  #plot(cv.lasso, main="AUC vs lambda (training)")
  lambda.star=cv.lasso$lambda.min
  nonzero=nonzero+as.matrix(coef(lasso.fit)[-(1:2),]!=0)
  nonzero.optim=nonzero.optim+
    as.numeric(coef(lasso.fit)[-(1:2),which(lam==lambda.star)]!=0)
  lasso.pred=predict(lasso.fit, newx=model.matrix( ~ ., test[,-c(1,exclude.idx-2)]), newy=test[,1],
                     family="binomial")
  lasso.auc=apply(lasso.pred, 2, auc, actual=test[,1])
  auc.lasso[j]=lasso.auc[which(lasso.fit$lambda==lambda.star)]
}

rf.auc.aki=mean(unlist(auc)) # average AUC from 50 MI datasets
par(mfrow=c(1,2)) # RF importance plot
dotchart(sort(rowMeans(MDacc),decreasing = FALSE),
         labels=(names(test)[-c(1,exclude.idx-2)])[order(rowMeans(MDacc),decreasing = FALSE)],
         pch=1, main="Mean Decrease Accuracy, AKI")
dotchart(sort(rowMeans(MDgini),decreasing = FALSE),
         labels=(names(test)[-c(1,exclude.idx-2)])[order(rowMeans(MDgini),decreasing = FALSE)],
         pch=1, main="Mean Decrease Gini, AKI")
MDacc.aki=MDacc
MDgini.aki=MDgini

lasso.auc.aki=mean(auc.lasso)
temp.rownames=rownames(coef(lasso.fit))[-(1:2)]
ord.rownames=rep(NA,length(temp.rownames))
for (j in 1:length(temp.rownames)) # create the row labels in the lasso var importance plots
  ord.rownames[j]=paste("(",nonzero.optim[j],") ",temp.rownames[j],sep="")
heatmap.2(matrix(as.numeric(nonzero[order(nonzero.optim,decreasing = TRUE),length(lam):1]),ncol=ncol(nonzero)), 
          scale="none", Colv = FALSE, Rowv = FALSE, labCol = "", xlab="log(lambda)",
          dendrogram = "none", trace = "none", key="FALSE",
          margins = c(5,0), lwid = c(2,8),
          col= colorRampPalette(colors=c("white","gray"))(11),#col=c("white","gray"),
          labRow = ord.rownames[order(nonzero.optim,decreasing = TRUE)],
          offsetRow = -60) # variable importance ranked at the *optimal* lambda
sum.nonzero=rowSums(nonzero)/length(lam)/imp$m*100
ord.sum.rownames=rep(NA,length(temp.rownames))
for (j in 1:length(temp.rownames))
  ord.sum.rownames[j]=paste("(",round(sum.nonzero[j],digits=2),") ",temp.rownames[j],sep="")
heatmap.2(matrix(as.numeric(nonzero[order(sum.nonzero,decreasing = TRUE),length(lam):1]),ncol=ncol(nonzero)), 
          scale="none", Colv = FALSE, Rowv = FALSE, labCol = "", xlab="log(lambda)",
          dendrogram = "none", trace = "none", key="FALSE",
          margins = c(5,0), lwid = c(2,8),
          col= colorRampPalette(colors=c("white","gray"))(11),#col=c("white","gray"),
          labRow = ord.sum.rownames[order(sum.nonzero,decreasing = TRUE)],
          offsetRow = -62) # variable importance ranked by the sum across *all* lambdas

# non-AKI, essentially just repeating what's done for the AKI group
set.seed(2333)
MDacc=MDgini=matrix(NA,nrow=n.var,ncol=imp$m)
auc=auc.lasso=rep(NA,imp$m)
nonzero=matrix(0,nrow=n.var+2,ncol=length(lam))
nonzero.optim=rep(0,n.var+2)
for (j in 1:imp$m){
  imp.j=imp.dat[imp.dat$.imp==j & imp.dat$V0_AKIN_stage=="non-AKI",-1]
  train.idx=sample(1:nrow(imp.j),round(nrow(imp.j)*0.8),replace=F)
  training=imp.j[train.idx,-1]
  test=imp.j[-train.idx,-1]
  
  #rf
  rf=randomForest(as.formula(rf.formula),data=training,importance=TRUE)
  MDacc[,j]=rf$importance[,3]/rf$importanceSD[,3]
  MDgini[,j]=rf$importance[,4]
  pred.test=predict(rf,test,type="prob")
  pred=prediction(pred.test[,2],test[,1])
  auc.perf = performance(pred, measure = "auc")
  auc[j]=auc.perf@y.values
  
  #lasso
  predictors.imp=training[,-c(1,exclude.idx-2)]
  outcome=training$bin3yr
  lasso.fit=glmnet(model.matrix( ~ ., predictors.imp),outcome,family = "binomial",
                   alpha = 1, lambda=lam)
  cv.lasso=cv.glmnet(x=model.matrix( ~ ., predictors.imp), y=outcome, family = "binomial",
                     type.measure="auc", nfolds=5, alpha=1, lambda=lam)
  lambda.star=cv.lasso$lambda.min
  nonzero=nonzero+as.matrix(coef(lasso.fit)[-(1:2),]!=0)
  nonzero.optim=nonzero.optim+
    as.numeric(coef(lasso.fit)[-(1:2),which(lam==lambda.star)]!=0)
  lasso.pred=predict(lasso.fit, newx=model.matrix( ~ ., test[,-c(1,exclude.idx-2)]), newy=test[,1],
                     family="binomial")
  lasso.auc=apply(lasso.pred, 2, auc, actual=test[,1])
  auc.lasso[j]=lasso.auc[which(lasso.fit$lambda==lambda.star)]
}
rf.auc.nonaki=mean(unlist(auc))
dotchart(sort(rowMeans(MDacc),decreasing = FALSE),
         labels=(names(test)[-c(1,exclude.idx-2)])[order(rowMeans(MDacc),decreasing = FALSE)],
         pch=1, main="Mean Decrease Accuracy, non-AKI")
dotchart(sort(rowMeans(MDgini),decreasing = FALSE),
         labels=(names(test)[-c(1,exclude.idx-2)])[order(rowMeans(MDgini),decreasing = FALSE)],
         pch=1, main="Mean Decrease Gini, non-AKI")
lasso.auc.nonaki=mean(auc.lasso)
temp.rownames=rownames(coef(lasso.fit))[-(1:2)]
ord.rownames=rep(NA,length(temp.rownames))
for (j in 1:length(temp.rownames))
  ord.rownames[j]=paste("(",nonzero.optim[j],") ",temp.rownames[j],sep="")
heatmap.2(matrix(as.numeric(nonzero[order(nonzero.optim,decreasing = TRUE),length(lam):1]),ncol=ncol(nonzero)), 
          scale="none", Colv = FALSE, Rowv = FALSE, labCol = "", xlab="log(lambda)",
          dendrogram = "none", trace = "none", key="FALSE",
          margins = c(5,0), lwid = c(2,8),
          col= colorRampPalette(colors=c("white","gray"))(11),#col=c("white","gray"),
          labRow = ord.rownames[order(nonzero.optim,decreasing = TRUE)],
          offsetRow = -60)
sum.nonzero=rowSums(nonzero)/length(lam)/imp$m*100
ord.sum.rownames=rep(NA,length(temp.rownames))
for (j in 1:length(temp.rownames))
  ord.sum.rownames[j]=paste("(",round(sum.nonzero[j],digits=2),") ",temp.rownames[j],sep="")
heatmap.2(matrix(as.numeric(nonzero[order(sum.nonzero,decreasing = TRUE),length(lam):1]),ncol=ncol(nonzero)), 
          scale="none", Colv = FALSE, Rowv = FALSE, labCol = "", xlab="log(lambda)",
          dendrogram = "none", trace = "none", key="FALSE",
          margins = c(5,0), lwid = c(2,8),
          col= colorRampPalette(colors=c("white","gray"))(11),#col=c("white","gray"),
          labRow = ord.sum.rownames[order(sum.nonzero,decreasing = TRUE)],
          offsetRow = -60)
