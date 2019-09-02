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
outimp<-fread("results/fullmodel_window.imp.tsv")
outimp[,imp:=Gain/sum(Gain),by=c("fold", tunepars)]
csimp<-outimp[,.(mean_imp=mean(imp)),by=.(Feature)][order(-mean_imp)]
usevars_window<-csimp[mean_imp >=0.001,Feature]
grpvars <- c("id", "fid", "fold", "grp", "class_label")



allcovs<-tolower(fread(paste0("../data/fold1_all.tsv"),nrow=1,header=F))[-1]
#csimp[,base:=sub(paste0(allcovs,collapse="|"), "\\1", Feature)]
library(stringr)
pat<-paste0(allcovs,collapse="|")
csimp[,base:=str_match(Feature, pat)[,1]]
csimp[,window:=tstrsplit(Feature, base)[[2]]]
csimp[,stat:=tstrsplit(window, "t10")[[1]]]
csimp[,window:=gsub("^.*t10_", "", window)]
csimp[,stat:=gsub("^_|_$", "", stat)]
chk<-csimp[Feature %in% usevars_window][order(-mean_imp)]
chk[,.N,by=base]
chk[,rnk:=1:.N]


csimp[,base:=sub(pat, "\\1", Feature[grepl(pat, Feature)])]
csimp[,base:=grep(paste0(allcovs,collapse="|"), Feature)]



#d1<-fread("big_window_195.tsv",select=c(grpvars, usevars_window))
#d2<-fread("big_window_d196.369.tsv",select=c(grpvars, usevars_window))
#d<-rbindlist(list(d1, d2), use.names=T, fill=T)
#fwrite(d, file="big_window.tsv",sep="\t")
d<-fread("big_window.tsv")

rm(d1,d2)

d[,class_label:=factor(class_label, levels=c("0", "1"), 
	labels=c("no_flare", "big_flare"))]

cnvrt_num<-melt(d[,lapply(.SD, class)][,id:="id"],
     id.vars="id",variable.factor=F)[value=="character"][!variable %in% grpvars,variable]

d[,paste(cnvrt_num):=lapply(.SD, as.numeric),.SDcols=cnvrt_num]
trainset<-d[fold!=4]
testset<-d[fold==4]

#if(tgf!="default"){
#  tg<-fread(tgf)
#  names(tg)<-unname(unlist(fread(tgf.header,header=FALSE,nrow=1)))
#  } else {
#  tg<-setDT(thresh_code$grid(
#		x=trainset[,-c("fraud"),with=FALSE], 
#		y=trainset[,fraud], len=3))
#  }

#load("fold10rep10.RData")
     
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


#tg<-CJ(max_depth=c(2,4,6),
#	subsample=0.5,
#	min_child_weight=1,
#	eta=0.05,
#	gamma=c(10,20),
#        colsample_bytree=c(0.5),
#        colsample_bylevel=c(0.5),
#	colsample_bynode=c(0.5))
#
##tg2<-CJ(max_depth=c(2, 6),
##	subsample=c(0.5, 0.8),
##	min_child_weight=1,
##	eta=0.05,
##	gamma=c(1, 10, 20),
##        alpha=c(10),
##        lambda=c(0,10),
##        colsample_bynode=c(0.1, 0.2, 0.5))
#
#
#for(j in 1:nrow(tg)){
#  resj<-list()
#  impj<-list()
#  for(i in 1:3){
#    tstv<-unlist(tg[j,])
#    tstl<-lapply(1:length(tstv), function(i) unname(tstv[i]))
#    names(tstl)<-names(tstv)
#    tstl$objective="binary:logistic"
#    tgj<-tg[j][,ktmp:=1]
#       dtest<-xgb.DMatrix(as.matrix(trainset[fold==i,
#       	-c(grpvars),with=FALSE]),
#       	label=trainset[fold==i,as.numeric(class_label=="big_flare")])
#       dtrain<-xgb.DMatrix(as.matrix(trainset[fold!=i,
#       	-c(grpvars),with=FALSE]),
#       	label=trainset[fold!=i,as.numeric(class_label=="big_flare")])
#       model1 <- xgb.train(param=tstl,
#                       	data=dtrain, 
#  			nthread=parcores,
#                          maximize=TRUE,
#                          feval=f1_5, 
#                          nrounds=1000,verbose=1,
#  	early_stopping_rounds=200, 
#       	watchlist=list(train=dtrain,eval=dtest)) 
#     
#       
#  
#       resj[[i]]<-data.table(
#             flare=predict(model1, newdata=dtest,
#  	   ntreelimit=model1$best_ntree_limit),
#             obs=getinfo(dtest,'label'),
#  	   best_iter=model1$best_iteration,
#             fold=i,
#             index=trainset[fold==i,which=TRUE],
#             ktmp=1)	 
#  
#  
#       if(length(model1$best_ntreelimit) > 0)	resj[[i]][,ntreelimit:=model1$best_ntreelimit]
#       impj[[i]]<- xgb.importance(feature_names=dimnames(dtrain)[[2]],
#       		      model=model1)[,`:=`(
#       		# index=names(foldlist)[i],
#  		 ktmp=1)]
#  
#  }
#     
#res<-rbindlist(resj)
#imp<-rbindlist(lapply(1:3, function(i) { impj[[i]][,fold:=i] }))
#fwrite(res, file=paste0("window_res.", j, ".tsv"),sep="\t")
#fwrite(imp, file=paste0("window_imp.", j, ".tsv"),sep="\t")
#outj<-merge(res,tgj,by="ktmp",all=TRUE,allow.cartesian=TRUE)[,ktmp:=NULL]
#outimpj<-merge(tgj,imp,by="ktmp",all=TRUE,allow.cartesian=TRUE)[,ktmp:=NULL]
#
#  if(j==1){
#   out<-outj
#   outimp<-outimpj 
#  }else{
#   out<-rbindlist(list(out, outj),use.names=T, fill=T)
#   outimp<-rbindlist(list(outimp, outimpj),use.names=T,fill=T )
#  }
#
#}
#
#
#print(out[,.(f1=f1(flare, obs,cut=0.35)),by=tunepars])
#
#fwrite(out, 
#	file=paste0("results/fullmodel_window.rprt.tsv"),
#	sep="\t")
#fwrite(outimp, 
#	file=paste0("results/fullmodel_window.imp.tsv"),
#	sep="\t")
#
#cts<-seq(0.05,0.5,0.05)
#out[,.(cts, f1=sapply(cts, function(x) f1(flare, obs,cut=x))),by=tunepars]
#
#
#outimp[,imp:=Gain/sum(Gain),by=c("fold", tunepars)]
#csimp<-outimp[,.(mean_imp=mean(imp)),by=.(Feature)][order(-mean_imp)]
#usevars_window<-csimp[mean_imp >=0.001,Feature]
#
#
##### 
#
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
                          nrounds=10000,verbose=1,
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
       		      model=model1)[,`:=`(
       		# index=names(foldlist)[i],
  		 ktmp=1)]
  
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

















oob_setting<-out[,.(f1=f1(flare,obs,0.35)),
		 	 by=c(tunepars)][order(-f1)][1]
tstv<-unlist(oob_setting[1,-c("f1"),with=F])
tstl<-lapply(1:length(tstv), function(i) unname(tstv[i]))
names(tstl)<-names(tstv)
tstl$objective="binary:logistic"
#tgj<-tg[j][,ktmp:=1]
        
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

merge(oob_setting, out,by=tunepars)[,f1(flare, obs, cut=0.4)]

fwrite(submit_dt[,.(Id,ClassLabel=as.numeric(flare > 0.4))], file=paste0("oobs/modfull_window_submit2.csv"))


#submit_df.to_csv("/gpfs/group/asb17/default/dan/dan_psubs_solar2.csv",index=False)
##
###rbindlist(resj)[,f1(flare, obs,cut=0.3),by=fold] 
###rbindlist(resj)[,f1(flare, obs,cut=0.3)] 
##
#
#
#outimp[,imp:=Gain/sum(Gain),by=c("fold", tunepars)]
#impsum<-outimp[,.(imp=mean(imp)),by=Feature][order(-imp)]
#impsum[,csimp:=cumsum(imp)/sum(imp)]
#
#
#allcovs<-tolower(fread(paste0("../data/fold1_all.tsv"),nrow=1,header=F))[-1]
#chk<-impsum[imp >=0.005][order(-imp)]
#chk[,base:=sub(paste0(allcovs,collapse="|"), "\\1", Feature)]
#pat<-paste0(allcovs,collapse="|")
#library(stringr)
#chk[,base:=str_match(Feature, pat)[,1]]
#chk[,t:=tstrsplit(Feature, base)[[1]]]
#chk[,m:=tstrsplit(Feature, base)[[2]]]
#chk[,.N,by=base]
#
#chk[,base:=sub(pat, "\\1", Feature[grepl(pat, Feature)])]
#chk[,base:=grep(paste0(allcovs,collapse="|"), Feature)]
#
#
#
#usevars<-c(grpvars, impsum[imp >=0.005,Feature])
#toptrain<-trainset[,usevars,with=F]
#toptest<-testset[,usevars,with=F]
#
#tg3<-CJ(max_depth=c(2,6,8),
#	subsample=0.5,
#	min_child_weight=1,
#	eta=0.07,
#	gamma=c(10,20,50),
#        alpha=c(10,30,50),
#        lambda=c(10,30,50),
#        colsample_bytree=c(0.1,0.5),
#	colsample_bynode=c(0.2))
#tg3[-c(1:17),alpha:=0]
#tg3[-c(1:17),lambda:=0]
#tg3<-tg3[,unique(.SD)]
#for(j in 1:nrow(tg3)){
#  resj<-list()
#  impj<-list()
#  for(i in 1:3){
#    tstv<-unlist(tg3[j,])
#    tstl<-lapply(1:length(tstv), function(i) unname(tstv[i]))
#    names(tstl)<-names(tstv)
#    tstl$objective="binary:logistic"
#    tgj<-tg3[j][,ktmp:=1]
#       dtest<-xgb.DMatrix(as.matrix(toptrain[fold==i,
#       	-c(grpvars),with=FALSE]),
#       	label=toptrain[fold==i,as.numeric(class_label=="big_flare")])
#       dtrain<-xgb.DMatrix(as.matrix(toptrain[fold!=i,
#       	-c(grpvars),with=FALSE]),
#       	label=toptrain[fold!=i,as.numeric(class_label=="big_flare")])
#       model1 <- xgb.train(param=tstl,
#                       	data=dtrain, 
#  			nthread=parcores,
#                          maximize=TRUE,
#                          feval=f1_5, 
#                          nrounds=2000,verbose=1,
#  	early_stopping_rounds=200, 
#       	watchlist=list(train=dtrain,eval=dtest)) 
#     
#       
#  
#       resj[[i]]<-data.table(
#             flare=predict(model1, newdata=dtest,
#  	   ntreelimit=model1$best_ntree_limit),
#             obs=getinfo(dtest,'label'),
#  	   best_iter=model1$best_iteration,
#             fold=i,
#             index=trainset[fold==i,which=TRUE],
#             ktmp=1)	 
#  
#  
#       if(length(model1$best_ntreelimit) > 0)	resj[[i]][,ntreelimit:=model1$best_ntreelimit]
#       impj[[i]]<- xgb.importance(feature_names=dimnames(dtrain)[[2]],
#       		      model=model1)[,`:=`(
#       		# index=names(foldlist)[i],
#  		 ktmp=1)]
#  
#  }
#     
#res<-rbindlist(resj)
#imp<-rbindlist(lapply(1:3, function(i) { impj[[i]][,fold:=i] }))
#
#outj<-merge(res,tgj,by="ktmp",all=TRUE,allow.cartesian=TRUE)[,ktmp:=NULL]
#outimpj<-merge(tgj,imp,by="ktmp",all=TRUE,allow.cartesian=TRUE)[,ktmp:=NULL]
# fwrite(outj[,.(tg3row=j,
#		cut=seq(0.05,0.95,0.05),
#		f1=sapply(seq(0.05,0.95,0.05), function(ct) f1(flare,obs,ct))),by=tunepars] ,
#	file=paste0("slimtune_",j,".csv"))
#  if(j==1){
#   out<-outj
#   outimp<-outimpj 
#  }else{
#   out<-rbindlist(list(out, outj),use.names=T, fill=T)
#   outimp<-rbindlist(list(outimp, outimpj),use.names=T,fill=T )
#  }
#
#}
#
#fwrite(out, 
#	file=paste0("results/slimmodel.rprt.tsv"),
#	sep="\t")
#fwrite(outimp, 
#	file=paste0("results/slimmodel.imp.tsv"),
#	sep="\t")
#
#
#
#
#oob_setting<-out[,.(f1=f1(flare,obs,0.35)),
#		 	 by=c(tunepars)][order(-f1)][1]
#tstv<-unlist(oob_setting[1,-c("f1"),with=F])
#tstl<-lapply(1:length(tstv), function(i) unname(tstv[i]))
#names(tstl)<-names(tstv)
#tstl$objective="binary:logistic"
##tgj<-tg[j][,ktmp:=1]
#        
#     dtestf<-xgb.DMatrix(as.matrix(testset[fold==4,
#     	-c(grpvars),with=FALSE]))
#     dtrainf<-xgb.DMatrix(as.matrix(trainset[fold!=4,
#     	-c(grpvars),with=FALSE]),
#     	label=trainset[fold!=4,as.numeric(class_label=="big_flare")])
#     model_f <- xgb.train(param=tstl,
#                     	data=dtrainf, 
#			nthread=2,
#                        maximize=TRUE,
#                        feval=f1_5, 
#                        nrounds=merge(oob_setting, out,by=tunepars)[,floor(mean(ntreelimit))],
#                       verbose=1,
#     	watchlist=list(train=dtrainf)) 
#
#submit_dt<- data.table(Id = testset[,id],
#                  ClassLabel = as.numeric(predict(model_f, dtestf) > 0.35),
#		  flare=predict(model_f, dtestf))
#fwrite(submit_dt, file=paste0("oobs/slimmod.csv"))
#
#
#
#
##allpars<-as.character(getModelInfo("xgbTree")[[1]]$parameters$parameter)
##   cutdts<-lapply(seq(.5,.99,.01), function(cutval){
##    p<-preddt[obs==1,
##    	.(tp=sum(fraud > cutval),fn=sum(fraud < cutval)),
##    	by=c(allpars,"rep")]
##    n<-preddt[obs!=1,
##    	.(tn=sum(fraud < cutval),fp=sum(fraud > cutval)),
##    	by=c(allpars,"rep")]
##  
##    s<-merge(n, p,by=c(allpars,"rep"))
##    s[,profit:=fp*-25 + tp*5 + fn*-5]
##    s[,cut:=cutval]
##  })
##
##} else {
##
## fwrite(tg, #	file=paste0(modelname, ".",tgf,".tsv"), 
##	sep="\t")
##
##}
#
