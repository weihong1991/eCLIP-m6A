library(randomForest)
library(foreach)
library(doParallel)
library(e1071)
library(rpart)
library(party)

# parallel setting: 12 core
registerDoParallel(12)

#----------------------------------
## load data ##

#load HepG2
load("RBP_eCLIP_RMBase_HEPG2_GSM908331_matrix_mean_model_abs.rda") 

#or, load K562
#load("RBP_eCLIP_RMBase_K562_Parental_PMID30297871_matrix_mean_model_abs.rda")

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

#----------------------------------
# turning for ntree

ntree<-c(100,200,300,400,500,600,800,1000,1500,2000,3000,4000,5000)

turning_ntree<-foreach(xxxx=ntree, .packages="randomForest")%dopar%{
  #------------------------------
  ## using all featues
  
  myfit <- randomForest(myy~., importance=T, data=tr.dat,ntree=xxxx)
  
  out<-lapply(list(train=tr.dat,predict=te.dat), function(aa.dat){
    tmp = predict(myfit, aa.dat[,-1], type="prob")
    res = data.frame(mod=aa.dat[,1], prob=tmp)
    
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
    myauc
    
    #++++++++++++ AUC plot 
    
    df<-data.frame(x=xx,y=yy)
    #ggplot(df,aes(x=x,y=y))+
    #  geom_line()
    
    imp = myfit$importance
    
    oob<-myfit$err.rate[,1]
    return(list(myauc=myauc,df=df,imp=imp,oob=oob))
  })
}
saveRDS(turning,paste0("/mount/ictr1/chenglab/whong/m6a/RBP/random100/ntree_100/ntree_",i),compress = "gzip")

  
#----------------------------------
# turning for mtry
data<-readRDS(paste0("/mount/ictr1/chenglab/whong/m6a/RBP/random100/all_prediction/",i))
tr.dat=data$tr.dat
te.dat=data$te.dat
mtry<-c(6,8,10,11,12,14,16,18,20,22,24,26,30,35,40,45)

turning<-foreach(xxxx=mtry, .packages="randomForest")%dopar%{
  #------------------------------
  ## using all featues
  
  myfit <- randomForest(myy~., importance=T, data=tr.dat,mtry=xxxx)
  
  tmp = predict(myfit, te.dat[,-1], type="prob")
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
  myauc
  
  #++++++++++++ AUC plot 
  
  df<-data.frame(x=xx,y=yy)
  #ggplot(df,aes(x=x,y=y))+
  #  geom_line()
  
  imp = myfit$importance
  
  oob<-myfit$err.rate[,1]
  return(list(myauc=myauc,df=df,imp=imp,oob=oob))
}
saveRDS(turning,paste0("/mount/ictr1/chenglab/whong/m6a/RBP/random100/mtry_100/mtry_",i),compress = "gzip")


