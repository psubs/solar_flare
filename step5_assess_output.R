
source("solar_flare/make_sumdt_window_func.R")
library(ggplot2)
library(data.table)
library(stringr)
library(xgboost)
library(parallelDist)
library(stringr)
## data read-in ##

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

d<-fread("window_select_data.tsv")
d[,class_label:=factor(class_label, levels=c("0", "1"), 
	labels=c("no_flare", "big_flare"))]
grpvars <- c("id", "fid", "fold", "grp", "class_label")
cnvrt_num<-melt(d[,lapply(.SD, class)][,id:="id"],
     id.vars="id",variable.factor=F)[value=="character"][!variable %in% grpvars,variable]

d[,paste(cnvrt_num):=lapply(.SD, as.numeric),.SDcols=cnvrt_num]


#### 
trainset_select<-d[fold!=4,c(grpvars, usevars_window),with=F]
testset_select<-d[fold==4,c(grpvars, usevars_window),with=F]


tg<-CJ(max_depth=c(2,4,8),
	subsample=0.5,
	min_child_weight=1,
	eta=0.01,
	gamma=c(0,10,30),
        colsample_bynode=c(0.1, 0.3, 0.5, 0.9))
tunepars<-names(tg)

for(i in 1:nrow(tg)){
  outi<-fread(paste0("window_res_select.", i, ".tsv"))
  outi<-merge(outi, tg[i][,ktmp:=1],by="ktmp")[,tgind:=i]
  if(i==1){
  out<-outi
  } else {
  out<-rbindlist(list(out,outi))
  }
  print(i)
}

## Select best tune and fit final model.
f1sum<-out[,.(tgind=unique(tgind),
	      f1=f1(flare,obs,0.35),
	      tss=tss(flare, obs, 0.35),
	      gs=gs(flare, obs, 0.35),
	      prec=prec(flare, obs, 0.35),
	      rec=rec(flare, obs, 0.35)),
		 	 by=c(tunepars)]#[order(-f1)]


tgbst<-f1sum[max_depth < 8][1,tgind]

for(i in 1:3){
  mod<-xgb.load(paste0("cv.tgind",tgbst, ".fold", i,".model"))
  dtest<-xgb.DMatrix(as.matrix(trainset_select[fold==i,
  	-c(grpvars),with=FALSE]),
  	label=trainset_select[fold==i,as.numeric(class_label=="big_flare")])
  dtrain<-xgb.DMatrix(as.matrix(trainset_select[fold!=i,
  	-c(grpvars),with=FALSE]),
  	label=trainset_select[fold!=i,as.numeric(class_label=="big_flare")])
  cvshapi<-data.table(predict(mod, newdata=dtest,approxcontrib=TRUE,predcontrib=TRUE))
  print("shap values calculated..")
  outf=out[tgind==tgbst & fold==i, -tunepars, with=F][,
           id:=1:.N,by=fold][,
	   `:=`(fid=paste0("f",fold,"_",id),tgind=NULL,ktmp=NULL,index=NULL)]
  cvshapi<-cbind(outf, cvshapi)
  levs<-c(0, 0.01, 0.2, 0.35, 0.6, 0.9, 1)
    for(s in 1:length(levs)){
      shw<-cvshapi[order(abs(flare - levs[s]))][1:2,fid]
      cvshapi[fid %in% shw, lev:=s]
      print(shw)
    }
    if(i==1){
      cvshap<-cvshapi 
    #  longf<-longfi
    }else{
      cvshap<-rbindlist(list(cvshap, cvshapi)) 
    #  longf<-rbindlist(list(longf, longfi))
    }
  print(i)
}

longf<-fread(paste0("imp_full.tsv"))
longf<-longf[fid %in% cvshap[!is.na(lev),fid]][order(fold, id,time)][,V1:=NULL]
setnames(longf, old=names(longf), new=tolower(names(longf)))
cvshap[,lapply(.SD, function(x) mean(abs(x))),.SDcols=shapcols]


outf=out[tgind==tgbst , -tunepars, with=F][,
         id:=1:.N,by=fold][,
	   `:=`(fid=paste0("f",fold,"_",id),tgind=NULL,ktmp=NULL,index=NULL)]
fwrite(outf, file="ieee/cvout.tsv", sep="\t")
fwrite(cvshap, file="ieee/cvshap.tsv",sep="\t")


mod<-xgb.load(paste0("modfull_window_pt.",tgbst, ".model"))
dtest<-xgb.DMatrix(as.matrix(testset_select[,
	-c(grpvars),with=FALSE]))
tsshap<-data.table(predict(mod, newdata=dtest,approxcontrib=TRUE,predcontrib=TRUE))
tsout<-fread("oobs/modfull_window_pt.14.csv")
fwrite(cbind(tsout, tsshap),file="ieee/testpred_shap.tsv", sep="\t")

cvshap<-fread("ieee/cvshap.tsv")
outf<-fread("ieee/cvout.tsv")
shapcols<-setdiff(names(cvshap),c(names(outf), "lev", "BIAS"))
cvcols<-setdiff(names(cvshap), c(shapcols))


#library(parallelDist)
library(cluster)
library(SHAPforxgboost)
tstshap<-fread("ieee/testpred_shap.tsv")
plt_tst<-tstshap[order(flare)][seq(1,.N,18)]#[,.N] # 15
plt_cv<-cvshap[order(flare)][seq(1,.N,20)]#[,.N]
fwrite(plt_cv, file="ieee/plt_subset_shapcv.tsv", sep="\t")
fwrite(plt_tst, file="ieee/plt_subset_shaptst.tsv", sep="\t")

plt_cv<-fread("plt_subset_shapcv.tsv")
plt_tst<-fread("plt_subset_shaptst.tsv")
cvcols<-c(names(plt_cv)[1:7], "BIAS", "lev")
shapcols<-setdiff(names(plt_cv), cvcols)
global_shap<-function(shapdt, shapcols){
  global<-melt(shapdt[,lapply(.SD, function(x) mean(abs(x))),.SDcols=shapcols][,id:="var"], 
		value.name="global_shap",
		id.vars="id")[,id:=NULL][order(-global_shap)]
  return(global)
}
allcovs<-fread("allcovsbase.tsv")[,unlist(covs)]
pat<-paste0(allcovs,collapse="|")

decomp<-function(dd){
  d<-copy(dd)
  d[,base:=str_match(variable, pat)[,1]]
  d[,window:=tstrsplit(variable, base)[[2]]]
  d[,stat:=tstrsplit(window, "t10")[[1]]]
  d[,window:=gsub("^.*t10_", "", window)]
  d[,stat:=gsub("^_|_$", "", stat)]
  return(d)
}

std1 <- function(x){
   return ((x - min(x, na.rm = T))/(max(x, na.rm = T) - min(x, na.rm = T)))
}

gcv<-global_shap(plt_cv, shapcols)
gtst<-global_shap(plt_tst, shapcols)

gcv<-decomp(gcv)[,rnk:=1:.N]
gtst<-decomp(gtst)[,rnk:=1:.N]

gcv[,mean(global_shap),by=base][order(-V1)]

pltcols<-intersect(shapcols,gcv[1:20,variable])

gcvl<-melt(plt_cv[,c(pltcols, "flare", "fid","obs"),with=F],
           measure.vars=pltcols,value.name="shap")

gcvl[,stdshap:=std1(shap),by="variable"]
















cvshap[order(flare)]
calls<-cvshap[flare > 0.35]
calls<-calls[order(flare)][seq(1,.N,5)]
#callsc<-calls[, lapply(.SD, function(x) c(scale(x))),
#	     .SDcols=shapcols]
#calls<-cbind(calls[,cvcols,with=F], callsc)
#bigd<-parDist(as.matrix(calls[,shapcols,with=F]))
bigd<-as.matrix(dist(as.matrix(calls[,shapcols,with=F])))
bigd<-cbind(calls[,cvcols,with=F],as.data.table(bigd))
colnames(bigd)[(length(cvcols) + 1):ncol(bigd)]<-bigd[,fid]

#fwrite(bigd, file="ieee/distance_matrix_selected_calls.tsv", sep="\t") 
#bigd<-fread("ieee/distance_matrix_selected_calls.tsv")

#hc1<-hclust(as.dist(bigd[,-c(cvcols),with=F]), 
#	    method="complete")
#
#save(hc1, "hc1.RData")
#load("hc1.RData")
#plot(hc1)
rect.hclust(hc1,k=10)
kcl<-cutree(hc1, k=3)
ss<-merge(cvshap, data.table(fid=names(kcl), clust=kcl), by="fid",sort=F)
ssm<-ss[,lapply(.SD, mean),by=.(clust),.SDcols=shapcols]
msh<-melt(ssm, id.vars="clust",value.name="shap")
msh[,diff:=max(shap) - min(shap),by=variable]
msh<-msh[order(-diff,-shap)]
dcast(, variable~clust, value.var="shap")
tmp<-dcast(msh, variable~clust,value.var="shap")
setnames(tmp, old=paste0(1:3), new=paste0("c_", 1:3))
tmp[,diff12:=c_1 - c_2]
tmp[,diff3_12:=c_3 - (c_2 + c_1)/2]
tmp[order(-abs(diff3_12)),rnk3_12:=1:.N]
tmp[order(-abs(diff12)),rnk12:=1:.N]


bigd[fid %in% bigd[clust==1,fid], 
ss[clust==3,ss[clust==3,fid],with=F]
subm<-as.matrix(ss[clust==3,ss[clust==3,fid],with=F])
rownames(subm)<-colnames(subm)
isSymmetric(subm)
sum(subm[upper.tri(subm)])
for(i in 1:10){
	for(j in 1:10){
 if(subm[i,j] !=subm[j,i]) print(i)
	}
}


chk<-calls[,shapcols,with=F]
wss <- function(d) {
  sum(scale(d, scale = FALSE)^2)
}

wrap <- function(i, hc, x) {
  cl <- cutree(hc, i)
  spl <- split(x, cl)
  wss <- sum(sapply(spl, wss))
  wss
}

1 - (sapply(1:10, wrap, hc=hc1, x=chk)) / wrap(1, hc1, chk)








res1<-agnes(as.dist(bigd[,-cvcols,with=F]),method="complete")
m <- c( "average", "single", "complete", "ward")
names(m) <- c( "average", "single", "complete", "ward")
res<-sapply(m, function(x) agnes(as.dist(bigd[,-cvcols,with=F]),method=x))
for(i in 1:length(m)){
  

	
}

longd<-melt(cvshap, id.vars=cvcols,value.name="shap")
imp<-longd[,.(imp=mean(abs(shap))),by=variable][order(-imp)][,rnk:=1:.N]


ggplot(imp[1:20], aes(x=variable, y=imp)) + geom_col() + coord_flip()

std1 <- function(x){
            return ((x - min(x, na.rm = T))/(max(x, na.rm = T) - min(x, na.rm = T)))
}

plotd<-merge(imp[rnk < 10], longd, by="variable")[order(rnk)]
plotd[,stdfval:=std1(shap),by=variable]
data_long<-copy(plotd)
devtools::install_github("liuyanguu/SHAPforxgboost")
library(SHAPforxgboost)
dilute<-x_bound<-NULL
  # check number of observations
  N_features <- setDT(data_long)[,uniqueN(variable)]
  if (is.null(dilute)) dilute = FALSE

  nrow_X <- nrow(data_long)/N_features # n per feature
  if (dilute!=0){
    # if nrow_X <= 10, no dilute happens
    dilute <- ceiling(min(nrow_X/10, abs(as.numeric(dilute)))) # not allowed to dilute to fewer than 10 obs/feature
    set.seed(1234)
    data_long <- data_long[sample(nrow(data_long),
                                  min(nrow(data_long)/dilute, nrow(data_long)/2))] # dilute
  }

  x_bound <- if (is.null(x_bound)) max(abs(data_long$shap))*1.1 else as.numeric(abs(x_bound))
  plot1 <- ggplot(data = data_long)+
    coord_flip(ylim = c(-x_bound, x_bound)) +
    # sina plot:
    ggforce::geom_sina(aes(x = variable, y = shap, color = stdfval),
              method = "counts", maxwidth = 0.7, alpha = 0.7) +
    # print the mean absolute value:
    geom_text(data = unique(data_long[, c("variable", "mean_value")]),
              aes(x = variable, y=-Inf, label = sprintf(label_format, mean_value)),
              size = 3, alpha = 0.7,
              hjust = -0.2,
              fontface = "bold") + # bold
    # # add a "SHAP" bar notation
    # annotate("text", x = -Inf, y = -Inf, vjust = -0.2, hjust = 0, size = 3,
    #          label = expression(group("|", bar(SHAP), "|"))) +
    scale_color_gradient(low="#FFCC33", high="#6600CC",
                         breaks=c(0,1), labels=c(" Low","High "),
                         guide = guide_colorbar(barwidth = 12, barheight = 0.3)) +
    theme_bw() +
    theme(axis.line.y = element_blank(),
          axis.ticks.y = element_blank(), # remove axis line
          legend.position="bottom",
          legend.title=element_text(size=10),
          legend.text=element_text(size=8),
          axis.title.x= element_text(size = 10)) +
    geom_hline(yintercept = 0) + # the vertical line
    # reverse the order of features, from high to low
    # also relabel the feature using `label.feature`
    scale_x_discrete(limits = rev(levels(data_long$variable)),
                     labels = label.feature(rev(levels(data_long$variable))))+
    labs(y = "SHAP value (impact on model output)", x = "", color = "Feature value  ")





#out[tgind==tgbst, -tunepars, with=F][,
#           id:=1:.N,by=fold][
#	   `:=`(fid=paste0("f",fold,"_",id),tgind=NULL,ktmp=NULL,index=NULL)],
#














#shapint<-predict(mod, newdata=dtest,approxcontrib=TRUE, predinteraction=TRUE)





for(i in 1:3){
    tstv<-unlist(tg[j,])
    tstl<-lapply(1:length(tstv), function(i) unname(tstv[i]))
    names(tstl)<-names(tstv)
    tstl$objective="binary:logistic"
    tgj<-tg[j][,ktmp:=1]
    model1 <- xgb.train(param=tstl,
                       	data=dtrain, 
  			nthread=parcores,
                        maximize=TRUE,
                        feval=f1_5, 
                        nrounds=10000,
                        verbose=1,
                 	early_stopping_rounds=200, 
                      	watchlist=list(train=dtrain,eval=dtest)) 
    xgb.save(model1, fname=paste0("hldt_f",i,".tg",tgind,".model"))
    shapval<-predict(model1, newdata=dtest, predcontrib=TRUE, approxcontrib=TRUE)
    mshap<-melt(data.table(shapval)[1==1][,id:=1:.N],id.vars="id")[,fold:=i][,fid:=paste0("f",fold,"_",id)]
    mshap<-merge(mshap, 
		 mshap[variable!="BIAS",
		       .(imp=mean(value)),by=variable][order(-abs(imp))][,rnk:=1:.N][],
		 by="variable",all.x=T)
    p<-predict(model1, dtest)[1] 
    pr<-mshap[id==1][,sum(value)]
    exp(pr) / (1 + exp(pr))


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


res<-fread(paste0("window_res_select.", tgind, ".tsv"))
imp<-fread(paste0("window_imp_select.", tgind, ".tsv"))
imp<-merge(imp[1==1][,ktmp:=1], tg[tgind][,ktmp:=1],all=TRUE, by="ktmp")[,ktmp:=NULL]
imp[,imp:=Gain/sum(Gain),by=c("fold", tunepars)]
csimp<-imp[,.(mean_imp=mean(imp)),by=.(Feature)][order(-mean_imp)]


res[,.(f1=f1(flare, obs, 0.35),
       prec=prec(flare, obs, 0.35), 
       rec=rec(flare, obs, 0.35))]

###
## original long format data, ld => long_d
#### 
for(i in 1:4){
  ldi<-fread(paste0("/gpfs/group/asb17/default/BD2019_solar/data/fold",i,"_all.tsv"))
  ldi[,fold:=i]
  ldi[,fid:=paste0("f",fold,"_",id)]
  if(i==1){
  ld<-ldi
  } else {
  ld<-rbindlist(list(ld,ldi),
	       use.names=TRUE,fill=TRUE)
  }
  print(i)
}
ld[,V1:=NULL]
onames <- copy(names(ld))
setnames(ld, old=onames, new=tolower(onames))
keys<-c("fold", "fid", "id")
resp<-c("class_label")
covs<-setdiff(names(ld), c(keys, resp))


###
## Featurized ts data wd => wide_d
###
wd<-fread("window_select_data.tsv")

allcovs<-tolower(fread(paste0("../data/fold1_all.tsv"),nrow=1,header=F))[-1]
pat<-paste0(allcovs,collapse="|")
csimp[,base:=str_match(Feature, pat)[,1]]
csimp[,window:=tstrsplit(Feature, base)[[2]]]
csimp[,stat:=tstrsplit(window, "t10")[[1]]]
csimp[,window:=gsub("^.*t10_", "", window)]
csimp[,stat:=gsub("^_|_$", "", stat)]

res[,`:=`(fid=paste0("f", fold, "_", 1:.N), id=1:.N )],by=fold]

shw<-merge(res[,.(flare, obs, fold, fid, id)],
      ld[fold!=4]
      , by=c("fid", "id", "fold"))[order(fold, id, time)]

shw[,uniqueN(fid)==res[,.N]]
shw[,.N,by=.(class_label, obs)]

set.seed(123)
id_sample<-c(shw[flare < quantile(flare, 0.1),sample(fid, 2)],
	     shw[flare > 0.3 & flare < 0.5, sample(fid, 2)],)


for(ct in c(0, 0.01, 0.2, 0.35, 0.6, 0.9, 1)){
 shwi<-res[order(abs(flare - ct))][1:2, .(fid, flare, obs)]
 if(ct==0){
   
 } else{
 
 }
}
shw[flare > 0.1 & flare < 0.35]
shw[flare > 0.35 & flare < 0.8]
shw[flare > 0.35 & flare < 0.8]














