library(parallel)
library(xgboost)

source("make_sumdt_window_func.R")
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
outimp<-fread("results/fullmodel_window.imp.tsv")
outimp[,imp:=Gain/sum(Gain),by=c("fold", tunepars)]
csimp<-outimp[,.(mean_imp=mean(imp)),by=.(Feature)][order(-mean_imp)]
usevars_window<-csimp[mean_imp >=0.001,Feature]
grpvars <- c("id", "fid", "fold", "grp", "class_label")


cts<-seq(0.05,0.5,0.05)
outimp[,.(cts, f1=sapply(cts, function(x) f1(flare, obs,cut=x))),by=tunepars]


outimp[,imp:=Gain/sum(Gain),by=c("fold", tunepars)]
csimp<-outimp[,.(mean_imp=mean(imp)),by=.(Feature)][order(-mean_imp)]
usevars_window<-csimp[mean_imp >=0.001,Feature]


datafiles<-paste0("data/ssdat_window_", )
for(i in 1:369){
  di<-fread(datafiles[i], select=c(grpvars, usevars))
  if(i==1){
  d<-di
  } else {
  d<-rbindlist(list(d,di))
  }
}



#### 

tg<-CJ(max_depth=c(2,4,8),
	subsample=0.5,
	min_child_weight=1,
	eta=0.01,
	gamma=c(0,10,30),
        colsample_bynode=c(0.1, 0.3, 0.5, 0.9))
tunepars<-names(tg)

trainset_select<-d[fold!=4,c(grpvars, usevars_window),with=F]
testset_select<-d[fold==4,c(grpvars, usevars_window),with=F]


for(j in 2:nrow(tg)){
  resj<-list()
  impj<-list()
  for(i in 1:3){
    tstv<-unlist(tg[j,])
    tstl<-lapply(1:length(tstv), function(i) unname(tstv[i]))
    names(tstl)<-names(tstv)
    tstl$objective="binary:logistic"
    tgj<-tg[j][,ktmp:=1]
    dtest<-xgb.DMatrix(as.matrix(trainset_select[fold==i,
    	-c(grpvars),with=FALSE]),
    	label=trainset_select[fold==i,as.numeric(class_label=="big_flare")])
    dtrain<-xgb.DMatrix(as.matrix(trainset_select[fold!=i,
    	-c(grpvars),with=FALSE]),
    	label=trainset_select[fold!=i,as.numeric(class_label=="big_flare")])
    model1 <- xgb.train(param=tstl,
                       	data=dtrain, 
  			nthread=parcores,
                        maximize=TRUE,
                        feval=f1_5, 
                        nrounds=10000,
                        verbose=1,
                 	early_stopping_rounds=200, 
                      	watchlist=list(train=dtrain,eval=dtest)) 
    resj[[i]]<-data.table(
          flare=predict(model1, newdata=dtest,
          ntreelimit=model1$best_ntree_limit),
          obs=getinfo(dtest,'label'),
          best_iter=model1$best_iteration,
          fold=i,
          index=trainset_select[fold==i,which=TRUE],
          ktmp=1)	 
  
  
       if(length(model1$best_ntreelimit) > 0)	resj[[i]][,ntreelimit:=model1$best_ntreelimit]
       impj[[i]]<- xgb.importance(feature_names=dimnames(dtrain)[[2]],
       		                  model=model1)[,`:=`(ktmp=1)]
  
  }
     
res<-rbindlist(resj)
imp<-rbindlist(lapply(1:3, function(i) { impj[[i]][,fold:=i] }))
fwrite(res, file=paste0("window_res_select.", j, ".tsv"),sep="\t")
fwrite(imp, file=paste0("window_imp_select.", j, ".tsv"),sep="\t")
outj<-merge(res,tgj,by="ktmp",all=TRUE,allow.cartesian=TRUE)[,ktmp:=NULL]
outimpj<-merge(tgj,imp,by="ktmp",all=TRUE,allow.cartesian=TRUE)[,ktmp:=NULL]

  if(j==2){
   out<-outj
   outimp<-outimpj 
  }else{
   out<-rbindlist(list(out, outj),use.names=T, fill=T)
   outimp<-rbindlist(list(outimp, outimpj),use.names=T,fill=T )
  }

}


#####
###
## Select best tune and fit final model.
####
###

oob_setting<-out[,.(f1=f1(flare,obs,0.35)),
		 	 by=c(tunepars)][order(-f1)][1]
tstv<-unlist(oob_setting[1,-c("f1"),with=F])
tstl<-lapply(1:length(tstv), function(i) unname(tstv[i]))
names(tstl)<-names(tstv)
tstl$objective="binary:logistic"
        
     dtestf<-xgb.DMatrix(as.matrix(testset_select[fold==4,
     	-c(grpvars),with=FALSE]))
     dtrainf<-xgb.DMatrix(as.matrix(trainset_select[fold!=4,
     	-c(grpvars),with=FALSE]),
     	label=trainset_select[fold!=4,as.numeric(class_label=="big_flare")])
     model_f <- xgb.train(param=tstl,
                       	  data=dtrainf, 
		          nthread=2,
                          maximize=TRUE,
                          feval=f1_5, 
                          nrounds=merge(oob_setting, out,by=tunepars)[,floor(mean(ntreelimit))],
                          verbose=1,
     	                  watchlist=list(train=dtrainf)) 

submit_dt<- data.table(Id = testset_select[,id],
                  ClassLabel = as.numeric(predict(model_f, dtestf) > 0.35),
		  flare=predict(model_f, dtestf))
fwrite(submit_dt, file=paste0("oobs/modfull_window.csv"))
fwrite(submit_dt[,.(Id,ClassLabel=as.numeric(flare > 0.35))], file=paste0("oobs/modfull_window_submit.csv"))


