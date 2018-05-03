#It is modified from huangclin_randomforest_gappedpeak_select5_predict_allhs.R,
#but use epitensor selfrun results of roadmap data, instead of peak count

library(randomForest)

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

otherh3k4me3 = read.delim("other_h3k4me3_gplength.txt",row.names=1)
colnames(otherh3k4me3) = sub('E(\\d+)','E\\1h3k4me3',colnames(otherh3k4me3))
otherh3k9ac = read.delim("other_h3k9ac_gplength.txt",row.names=1)
colnames(otherh3k9ac) = sub('E(\\d+)','E\\1h3k9ac',colnames(otherh3k9ac))
otherh3k27me3 = read.delim("other_h3k27me3_gplength.txt",row.names=1)
colnames(otherh3k27me3) = sub('E(\\d+)','E\\1h3k27me3',colnames(otherh3k27me3))
otherh2az = read.delim("other_h2a.z_gplength.txt",row.names=1)
colnames(otherh2az) = sub('E(\\d+)','E\\1h2az',colnames(otherh2az))
other=cbind(otherh3k4me3,otherh3k9ac,otherh3k27me3,otherh2az)

overlap = as.vector(intersect(rownames(huang), rownames(hs)))
length(overlap)
huang = huang[!(rownames(huang) %in% overlap), ] #349 genes left
hs = hs[!(rownames(hs) %in% overlap), ] #800 genes left

sum(rownames(other) %in% rownames(huang)) #test for overlap, should be zero
sum(rownames(other) %in% rownames(hs))

exac = read.table("exac_list_pli0.1_elof10.txt",header =T)
moved = huang[ (rownames(huang) %in% exac$Gene),]
huang = huang[! (rownames(huang) %in% exac$Gene),] #287 genes left
exac2 = read.table("exac_list3.txt",header =T)
moved = rbind(moved,hs[ (rownames(hs) %in% exac2$Gene),]) #145 genes
hs = hs[! (rownames(hs) %in% exac2$Gene),] #717 genes left

other = rbind(other,moved)

epitensor = read.delim("epitensor_count.txt",row.names=1)
colnames(epitensor) = sub('E(\\d+)','E\\1tensor',colnames(epitensor))
hitensor = epitensor[match(rownames(huang),rownames(epitensor)),]
which(as.vector(is.na(hitensor)))
huang = cbind(huang,hitensor)
######
hstensor = epitensor[match(rownames(hs),rownames(epitensor)),]
which(as.vector(is.na(hstensor)))
hs = cbind(hs,hstensor)
######
epitensor2 = read.delim("epitensor_count_predict.txt",row.names=1)
colnames(epitensor2) = sub('E(\\d+)','E\\1tensor',colnames(epitensor2))
epitensor2 = rbind(epitensor2,epitensor)
othertensor = epitensor2[match(rownames(other),rownames(epitensor2)),]
which(as.vector(is.na(othertensor)))
other = cbind(other,othertensor)


hs$type = 0
huang$type = 1
hs$type = as.factor(hs$type)
huang$type = as.factor(huang$type)

nrandom = 30 #the randomization times to obtain the same number of hs as the huangclin set
epi.rf.pr = list()
hipredicted = list()
hspredicted = list()
times = 0
for (i in 1:nrandom){
  train = rbind(hs,huang)
  
  times = times + 1
  epi.rf = randomForest(type~.,data=train, ntree = 2000)
  epi.rf.pr[[times]] = predict(epi.rf,type="prob",newdata=other)[,2]
  hipredicted[[times]] = predict(epi.rf,type="prob",newdata=huang)[,2]
  hspredicted[[times]] = predict(epi.rf,type="prob",newdata=hs)[,2]
}

predicted = matrix(unlist(epi.rf.pr),ncol = nrandom)
rownames(predicted) = rownames(other)
write.table(apply(predicted,1,mean), "predicted_prob_allhs_epitensor.txt",quote = F,sep = "\t",col.names = F)

hipre = matrix(unlist(hipredicted),ncol = nrandom)
hspre = matrix(unlist(hspredicted),ncol = nrandom)
summary(apply(hipre,1,mean))
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.2810  0.7626  0.8368  0.8300  0.9101  0.9759 
summary(apply(hspre,1,mean))
#     Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
#3.333e-05 2.423e-02 4.945e-02 6.818e-02 9.785e-02 3.108e-01 
quantile(apply(hspre,1,mean),probs = seq(0,1,0.1))
#          0%          10%          20%          30%          40%          50%          60%          70%          80% 
#3.333333e-05 9.006667e-03 1.847333e-02 2.846000e-02 3.817000e-02 4.945000e-02 6.272000e-02 8.506667e-02 1.156633e-01 
#90%         100% 
#  1.629033e-01 3.108167e-01 
tail(sort(apply(hspre,1,mean)))
#[1] 0.2645833 0.2696667 0.2735667 0.2973500 0.3035500 0.3108167
quantile(apply(hipre,1,mean),probs = seq(0,1,0.1))
#       0%       10%       20%       30%       40%       50%       60%       70%       80%       90%      100% 
#0.2810000 0.7046100 0.7449067 0.7785233 0.8102633 0.8368333 0.8691500 0.9007100 0.9189767 0.9389900 0.9759333 
head(sort(apply(hipre,1,mean)))
#[1] 0.2810000 0.6124333 0.6478333 0.6485167 0.6514833 0.6617167

par(mfrow=c(1,2))
hist(apply(hipre,1,mean),breaks = 40,xlim = c(0,1),ylim = c(0,13),freq = F,col=rgb(0,0,1,1/4),main='Training Genes',xlab = 'Predicted Probability of Being HIS',lwd=2)
hist(apply(hspre,1,mean),breaks = 20,xlim = c(0,1),ylim = c(0,13),freq = F,col=rgb(1,0,0,1/4),add = T,lwd=2)
hist(apply(predicted,1,mean),breaks = 40,freq= F,ylim =c(0,2.5),xlim = c(0,1),xlab = "Predicted Probability of Being HIS",main = "All Other Genes",col = 'red',lwd=2)

rownames(hipre) = rownames(huang)
write.table(apply(hipre,1,mean), "predicted_prob_hi_allhs_epitensor.txt",quote = F,sep = "\t",col.names = F)
rownames(hspre) = rownames(hs)
write.table(apply(hspre,1,mean), "predicted_prob_hs_allhs_epitensor.txt",quote = F,sep = "\t",col.names = F)

plot(hipre[,19],hipre[,9])
plot(hipre[,15],hipre[,4])

plot(hspre[,7],hspre[,13])
plot(hspre[,21],hspre[,14])

plot(predicted[,5],predicted[,16])
plot(predicted[,11],predicted[,25])
