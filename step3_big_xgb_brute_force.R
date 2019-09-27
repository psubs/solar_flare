args=commandArgs(trailingOnly=TRUE) 

library(parallel)
library(xgboost)

parcores<-as.numeric(args[1])

set.seed(11)

tg<-CJ(max_depth=c(2,4,6),
	subsample=0.5,
	min_child_weight=1,
	eta=0.05,
	gamma=c(10,20),
        colsample_bytree=c(0.5),
        colsample_bylevel=c(0.5),
	colsample_bynode=c(0.5))

tunepars<-names(tg)
grpvars <- c("id", "fid", "fold", "grp", "class_label")


### read in data ###
datafiles<-paste0("data/ssdat_window_", 1:369)
for(i in 1:369){
  di<-fread(datafiles[i])
  if(i==1){
  d<-di
  } else {
  d<-rbindlist(list(d,di))
  }
}

################

d[,class_label:=factor(class_label, levels=c("0", "1"), 
	labels=c("no_flare", "big_flare"))]

cnvrt_num<-melt(d[,lapply(.SD, class)][,id:="id"],
     id.vars="id",variable.factor=F)[value=="character"][!variable %in% grpvars,variable]

d[,paste(cnvrt_num):=lapply(.SD, as.numeric),.SDcols=cnvrt_num]
trainset<-d[fold!=4]
testset<-d[fold==4]

     
f1<-function(pred, label,cut=0.5){
        tp<-sum(label* (pred > cut))
        fp<-sum((1-label)* (pred > cut))
        fn<-sum(label* (pred <= cut))
        prec<-tp / (tp + fp)
        rec<-tp / (tp + fn)
        return(2*prec*rec / (prec + rec))
      }

f1_5 <- function(preds, dtrain) {
	labels <- getinfo(dtrain, "label")
#	preds<-1/(1 + exp(-preds)) 
	tp<-sum(labels * 
		    (preds > 0.5)) 
	fp<-sum((1-labels) * 
		    (preds > 0.5)) 
	fn<-sum(labels * 
		    (preds <= 0.5)) 
        prec<-tp / (tp + fp)
        rec<-tp / (tp + fn)
        f1<-2*prec*rec / (prec + rec)	
        return(list(metric = c("f1-score"), 
		    value = c(f1)))
}


tg<-CJ(max_depth=c(2,4,6),
	subsample=0.5,
	min_child_weight=1,
	eta=0.05,
	gamma=c(10,20),
        colsample_bytree=c(0.5),
        colsample_bylevel=c(0.5),
	colsample_bynode=c(0.5))



for(j in 1:nrow(tg)){
  resj<-list()
  impj<-list()
  for(i in 1:3){
    tstv<-unlist(tg[j,])
    tstl<-lapply(1:length(tstv), function(i) unname(tstv[i]))
    names(tstl)<-names(tstv)
    tstl$objective="binary:logistic"
    tgj<-tg[j][,ktmp:=1]
    dtest<-xgb.DMatrix(as.matrix(trainset[fold==i,
    	-c(grpvars),with=FALSE]),
    	label=trainset[fold==i,as.numeric(class_label=="big_flare")])
    dtrain<-xgb.DMatrix(as.matrix(trainset[fold!=i,
    	-c(grpvars),with=FALSE]),
    	label=trainset[fold!=i,as.numeric(class_label=="big_flare")])
    model1 <- xgb.train(param=tstl,
                       	data=dtrain, 
  			nthread=parcores,
                        maximize=TRUE,
                        feval=f1_5, 
                        nrounds=1000,
			verbose=1,
                	early_stopping_rounds=200, 
                     	watchlist=list(train=dtrain,eval=dtest)) 
  
    resj[[i]]<-data.table(
          flare=predict(model1, newdata=dtest,
          ntreelimit=model1$best_ntree_limit),
          obs=getinfo(dtest,'label'),
          best_iter=model1$best_iteration,
          fold=i,
          index=trainset[fold==i,which=TRUE],
          ktmp=1)	 
  
  
    if(length(model1$best_ntreelimit) > 0)	resj[[i]][,ntreelimit:=model1$best_ntreelimit]
    impj[[i]]<- xgb.importance(feature_names=dimnames(dtrain)[[2]],
    		               model=model1)[,`:=`(ktmp=1)]
  
  }
     
res<-rbindlist(resj)
imp<-rbindlist(lapply(1:3, function(i) { impj[[i]][,fold:=i] }))
fwrite(res, file=paste0("window_res.", j, ".tsv"),sep="\t")
fwrite(imp, file=paste0("window_imp.", j, ".tsv"),sep="\t")
outj<-merge(res,tgj,by="ktmp",all=TRUE,allow.cartesian=TRUE)[,ktmp:=NULL]
outimpj<-merge(tgj,imp,by="ktmp",all=TRUE,allow.cartesian=TRUE)[,ktmp:=NULL]

  if(j==1){
   out<-outj
   outimp<-outimpj 
  }else{
   out<-rbindlist(list(out, outj),use.names=T, fill=T)
   outimp<-rbindlist(list(outimp, outimpj),use.names=T,fill=T )
  }

}


print(out[,.(f1=f1(flare, obs,cut=0.35)),by=tunepars])

fwrite(out, 
	file=paste0("results/fullmodel_window.rprt.tsv"),
	sep="\t")
fwrite(outimp, 
	file=paste0("results/fullmodel_window.imp.tsv"),
	sep="\t")


