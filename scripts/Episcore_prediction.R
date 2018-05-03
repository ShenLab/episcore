
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

hs.pre = hs

pli = read.csv("~/Desktop/Episcore_working/enrichment/data/exac_pli_combined_addtarget.csv", header = TRUE)
hs.pli = pli[match(rownames(hs), pli$EnsemblID), ]
length(which(hs.pli$pLI>0.2))

hs.pliH0.2 = hs.pli[which(hs.pli$pLI>0.2),]

hs = hs.pre[!(rownames(hs.pre) %in% hs.pliH0.2$EnsemblID), ]
hs.out = hs.pre[rownames(hs.pre) %in% hs.pliH0.2$EnsemblID, ]
other = rbind(other,hs.out)

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
write.table(apply(predicted,1,mean), "predicted_prob_otherwithoutpli>0.2_epitensor.txt",quote = F,sep = "\t",col.names = F)

hipre = matrix(unlist(hipredicted),ncol = nrandom)
hspre = matrix(unlist(hspredicted),ncol = nrandom)

rownames(hipre) = rownames(huang)
write.table(apply(hipre,1,mean), "predicted_prob_hiwithoutpli>0.2_epitensor.txt",quote = F,sep = "\t",col.names = F)
rownames(hspre) = rownames(hs)
write.table(apply(hspre,1,mean), "predicted_prob_hswithoutpli>0.2_epitensor.txt",quote = F,sep = "\t",col.names = F)

other.mean = as.data.frame(apply(predicted,1,mean))
colnames(other.mean) = "Episcore"
other.mean$batch = "testing"
hs.mean = as.data.frame(apply(hspre,1,mean))
colnames(hs.mean) = "Episcore"
hs.mean$batch = "HStraining"
hi.mean = as.data.frame(apply(hipre,1,mean))
colnames(hi.mean) = "Episcore"
hi.mean$batch = "HItraining"

all = rbind(hi.mean, hs.mean, other.mean)
write.table(all, "predicted_prob_allgenes_epitensor.txt",quote = F,sep = "\t",col.names = F)
