#This script is for 10-fold cross-validation. It is modified from huangclin_randomforest_cv10_gappedpeak_select5_allhs.R,
#but use epitensor selfrun results of roadmap data, instead of peak count

library(randomForest)
library(ROCR)

huangh3k4me3 = read.delim("huangclin_h3k4me3_gplength.txt",row.names=1)
colnames(huangh3k4me3) = sub('E(\\d+)','E\\1h3k4me3',colnames(huangh3k4me3))
huangh3k9ac = read.delim("huangclin_h3k9ac_gplength.txt",row.names=1)
colnames(huangh3k9ac) = sub('E(\\d+)','E\\1h3k9ac',colnames(huangh3k9ac))
huangh3k27me3 = read.delim("huangclin_h3k27me3_gplength.txt",row.names=1)
colnames(huangh3k27me3) = sub('E(\\d+)','E\\1h3k27me3',colnames(huangh3k27me3))
huangh2az = read.delim("huangclin_h2a.z_gplength.txt",row.names=1)
colnames(huangh2az) = sub('E(\\d+)','E\\1h2az',colnames(huangh2az))
huang=cbind(huangh3k4me3,huangh3k9ac,huangh3k27me3,huangh2az)

hsh3k4me3 = read.delim("hs_h3k4me3_gplength_select.txt",row.names=1)
colnames(hsh3k4me3) = sub('E(\\d+)','E\\1h3k4me3',colnames(hsh3k4me3))
hsh3k9ac = read.delim("hs_h3k9ac_gplength_select.txt",row.names=1)
colnames(hsh3k9ac) = sub('E(\\d+)','E\\1h3k9ac',colnames(hsh3k9ac))
hsh3k27me3 = read.delim("hs_h3k27me3_gplength_select.txt",row.names=1)
colnames(hsh3k27me3) = sub('E(\\d+)','E\\1h3k27me3',colnames(hsh3k27me3))
hsh2az = read.delim("hs_h2a.z_gplength_select.txt",row.names=1)
colnames(hsh2az) = sub('E(\\d+)','E\\1h2az',colnames(hsh2az))
hs=cbind(hsh3k4me3,hsh3k9ac,hsh3k27me3,hsh2az)

overlap = as.vector(intersect(rownames(huang), rownames(hs)))
length(overlap)
huang = huang[!(rownames(huang) %in% overlap), ] #349 genes left
hs = hs[!(rownames(hs) %in% overlap), ] #800 genes left

exac = read.table("exac_list_pli0.1_elof10.txt",header =T)
huang = huang[! (rownames(huang) %in% exac$Gene),] #287 genes
exac2 = read.table("exac_list3.txt",header =T)
hs = hs[! (rownames(hs) %in% exac2$Gene),] #717 genes


epitensor = read.delim("epitensor_count.txt",row.names=1)
colnames(epitensor) = sub('E(\\d+)','E\\1tensor',colnames(epitensor))
hitensor = epitensor[match(rownames(huang),rownames(epitensor)),]
which(as.vector(is.na(hitensor)))
huang = cbind(huang,hitensor)
hstensor = epitensor[match(rownames(hs),rownames(epitensor)),]
which(as.vector(is.na(hstensor)))
hs = cbind(hs,hstensor)

hs$type = 0
huang$type = 1
hs$type = as.factor(hs$type)
huang$type = as.factor(huang$type)


nrandom = 100 #the randomization times to split 90/10
epi.rf.pr = list()
label = list()
times = 0
for (i in 1:nrandom){
  times = times + 1
  
  ind = sample(nrow(hs), round(nrow(hs)*0.9,digits = 0), replace=F)
  ind2 = sample(nrow(huang), round(nrow(huang)*0.9,digits = 0), replace=F)
  
  hstrain = hs[ind,]
  hstest = hs[-ind,]
  huangtrain = huang[ind2,]
  huangtest = huang[-ind2,]
  train = rbind(hstrain,huangtrain)
  test = rbind(hstest,huangtest)
  
  epi.rf = randomForest(type~.,data=train, ntree = 2000)
  epi.rf.pr[[times]] = predict(epi.rf,type="prob",newdata=test)[,2]
  label[[times]] = test$type
}

epi.rf.pred = prediction(epi.rf.pr, label)
epi.rf.perf = performance(epi.rf.pred,"tpr","fpr")
plot(epi.rf.perf,main="ROC Curve of HIS and HS",col=2,lwd=2, avg='threshold', spread.estimate='stddev')
abline(a=0,b=1,lwd=2,lty=2,col="gray")
summary(as.numeric(cbind(performance(epi.rf.pred,"auc")@y.values)))
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.7821  0.8596  0.8879  0.8806  0.9098  0.9583
text(0.7,0.08,"mean AUC = 0.88")
text(0.7,0.2,"median AUC = 0.89")

###for final plot###
epi.rf.pred = prediction(epi.rf.pr, label)
epi.rf.perf = performance(epi.rf.pred,"tpr","fpr")
par(mar=c(5.1,4.5,1,2.1))
plot(epi.rf.perf,col=2,lwd=2, avg='threshold', spread.estimate='stddev',cex.axis =1.2,cex.lab = 1.5)
abline(a=0,b=1,lwd=2,lty=2,col="gray")
summary(as.numeric(cbind(performance(epi.rf.pred,"auc")@y.values)))
text(0.8,0.08,"mean AUC = 0.88",cex = 1.3)
####################

epi.rf.f = performance(epi.rf.pred,"f")
plot(epi.rf.f)
