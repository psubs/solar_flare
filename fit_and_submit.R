
args=commandArgs(trailingOnly=T)
library(parallel)
library(xgboost)

setwd("/gpfs/group/asb17/default/BD2019_solar/dan/")
source("solar_flare/make_sumdt_window_func.R")
#args<-c("2", "2")
parcores<-as.numeric(args[1])
tgind<-as.numeric(args[2])

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
usevars_window<-csimp[mean_imp >=.001,Feature]
grpvars <- c("id", "fid", "fold", "grp", "class_label")


out<-fread("results/fullmodel_window.rprt.tsv")
cts<-seq(0.05,0.5,0.05)
out[,.(cts, f1=sapply(cts, function(x) f1(flare, obs,cut=x))),by=tunepars]


outimp[,imp:=Gain/sum(Gain),by=c("fold", tunepars)]
csimp<-outimp[,.(mean_imp=mean(imp)),by=.(Feature)][order(-mean_imp)]
usevars_window<-csimp[mean_imp >=0.001,Feature]

d<-fread("window_select_data.tsv")
d[,class_label:=factor(class_label, levels=c("0", "1"), 
	labels=c("no_flare", "big_flare"))]

cnvrt_num<-melt(d[,lapply(.SD, class)][,id:="id"],
     id.vars="id",variable.factor=F)[value=="character"][!variable %in% grpvars,variable]

d[,paste(cnvrt_num):=lapply(.SD, as.numeric),.SDcols=cnvrt_num]


tg<-CJ(max_depth=c(2,4,8),
	subsample=0.5,
	min_child_weight=1,
	eta=0.01,
	gamma=c(0,10,30),
        colsample_bynode=c(0.1, 0.3, 0.5, 0.9))
tunepars<-names(tg)

trainset_select<-d[fold!=4,c(grpvars, usevars_window),with=F]
testset_select<-d[fold==4,c(grpvars, usevars_window),with=F]


out<-fread(paste0("window_res_select.", tgind, ".tsv"))
out<-merge(out, tg[tgind][,ktmp:=1],by="ktmp")


tstv<-unlist(out[1,tunepars,with=F])
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
                  nrounds=out[,floor(mean(ntreelimit))],
                  verbose=1,
	          watchlist=list(train=dtrainf)) 
xgb.save(model_f, fname=paste0("modfull_window_pt.",tgind,".model"))
submit_dt<- data.table(Id = testset_select[,id],
                  ClassLabel = as.numeric(predict(model_f, dtestf) > 0.35),
		  flare=predict(model_f, dtestf))

fwrite(submit_dt, file=paste0("oobs/modfull_window_pt.",tgind,".csv"))
fwrite(submit_dt[,.(Id,ClassLabel=as.numeric(flare > 0.35))], 
       file=paste0("oobs/modfull_window_pt.",tgind,".submit.csv"))







