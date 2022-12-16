library(randomForest)
library(ggplot2)
#----------------------------------
## load data ##

#load HepG2
load("RBP_eCLIP_RMBase_HEPG2_GSM908331_matrix_mean_model_abs.rda") 

#or, load K562
#load("RBP_eCLIP_RMBase_K562_Parental_PMID30297871_matrix_mean_model_abs.rda")

#or, load histone
#load("RBP_eCLIP_RMBase_HEPG2_GSM908331_matrix_mean_model_Histone_abs.rda")
#load("RBP_eCLIP_RMBase_K562_Parental_PMID30297871_matrix_mean_model_Histone_abs.rda")

#----------------------------------

tr.siz = 5000 # random sample 5000 sites

se = sample(1:nrow(pos))[1:tr.siz] # sampling positive sitess
tr.pos = pos[se,]
te.pos = pos[-se,]

se = sample(1:nrow(neg))[1:tr.siz] # sampling negative sitess
tr.neg = neg[se,]
te.neg = neg[-se,]

myy = c(rep("Y", nrow(tr.pos)), rep("N", nrow(tr.neg))) # generate training data
tr.dat = rbind(tr.pos, tr.neg)
tr.dat = data.frame(myy=factor(myy), tr.dat)

myy = c(rep("Y", nrow(te.pos)), rep("N", nrow(te.neg))) # generate predicting data
te.dat = rbind(te.pos, te.neg)
te.dat = data.frame(myy=factor(myy), te.dat)

#------------------------------
## fit the random forest model ##
myfit <- randomForest(myy~.,importance=T,data=tr.dat) # training

tmp = predict(myfit, te.dat[,-1], type="prob") # predicting

# calculate AUC
res = data.frame(mod=te.dat[,1], prob=tmp)
thr = (1:99)*0.01
yy =  xx =  rep(0, length(thr))
for(j in 1:length(thr))
{
  aa = sum(res[,"prob.Y"]>=thr[j] & res[,1]=="Y")
  bb = sum(res[,"prob.Y"]<thr[j] & res[,1]=="Y")
  cc = sum(res[,"prob.Y"]>=thr[j] & res[,1]=="N")
  dd = sum(res[,"prob.Y"]<thr[j] & res[,1]=="N")
  yy[j] = aa/(aa+bb)
  xx[j] = cc/(cc+dd)
}
xx = c(1, xx, 0)
yy = c(1, yy, 0)
tmp1 = tmp2 = rep(0,100)
for(j in 1:100)
{
  tmp1[j] = xx[j]-xx[j+1]
  tmp2[j] = (yy[j+1]+yy[j])/2	
}
myauc = sum(tmp1*tmp2)

df<-data.frame(x=xx,y=yy) # data for AUC plot

imp=as.data.frame(myfit$importance) # relative importance

#---------------------------------
# AUC plot
ggplot(df,aes(x=x,y=y))+
  geom_line()+
  coord_equal()+
  labs(x="False positive",y="False negative",title = paste0("AUC=",signif(myauc,2)))+
  theme_classic()

#---------------------------------
# relative importance plot

imp<-data.frame(RBP=row.names(imp),imp)
imp<-imp[order(imp$MeanDecreaseGini,decreasing = T),]
imp$RBP<-factor(as.character(imp$RBP),levels = as.character(imp$RBP))

ggplot(imp,aes(x=RBP,y=MeanDecreaseGini))+
  geom_bar(stat = "identity")+
  labs(x="RBP",y="Relative importance")+
  theme_classic()+
  theme(axis.text.x = element_text(angle=90,hjust=1,vjust=0.5))
