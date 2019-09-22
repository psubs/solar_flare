allcovs<-fread("ieee/allcovsbase.tsv")[,unlist(covs)]
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

global_shap<-function(shapdt, shapcols){
  global<-melt(shapdt[,lapply(.SD, function(x) mean(abs(x))),.SDcols=shapcols][,id:="var"], 
		value.name="global_shap",
		id.vars="id")[,id:=NULL][order(-global_shap)]
  return(global)
}

source("solar_flare/make_sumdt_window_func.R")
library(ggplot2)
library(data.table)
library(stringr)
library(xgboost)
library(parallelDist)
library(stringr)
## data read-in ##



d<-fread("window_select_data.tsv")
d[,class_label:=factor(class_label, levels=c("0", "1"), 
		       	labels=c("no_flare", "big_flare"))]
grpvars <- c("id", "fid", "fold", "grp", "class_label")
cnvrt_num<-melt(d[,lapply(.SD, class)][,id:="id"],
		   id.vars="id",
		   variable.factor=F)[value=="character"][!variable %in% grpvars,variable]

d[,paste(cnvrt_num):=lapply(.SD, as.numeric),.SDcols=cnvrt_num]
d[,(paste0("q_", shapcols)):=lapply(.SD, function(x) rank(x) / .N),.SDcols=shapcols]


cvshap<-fread("ieee/cvshap.tsv")
tstshap<-fread("ieee/testpred_shap.tsv")
setnames(tstshap, old="Id", new="id")
tstshap<-tstshap[,`:=`(fold=4)][,fid:=paste0("f",fold,"_",id)]
allshap<-rbindlist(list(cvshap,tstshap),use.names=T,fill=T)

cvcols<-c(names(allshap)[1:7], "BIAS", "lev", "ClassLabel", "Id")
shapcols<-setdiff(names(allshap), cvcols)

gcv<-global_shap(allshap, shapcols)
#gtst<-global_shap(tstshap, shapcols)

gcv<-decomp(gcv)[,rnk:=1:.N]
#gtst<-decomp(#gtst)[,rnk:=1:.N]

gcv[,base_shap_sum:=sum(abs(global_shap)),by=base]
gcv[,window_shap_sum:=sum(abs(global_shap)),by=window]
gcv[,stat_shap_sum:=sum(abs(global_shap)),by=stat]

#gtst[,base_shap_sum:=sum(abs(global_shap)),by=base]
#gtst[,stat_shap_sum:=sum(abs(global_shap)),by=stat]
#gtst[,window_shap_sum:=sum(abs(global_shap)),by=window]
shaplong<-melt(allshap[,c(shapcols, "flare", "fid","obs", "fold"),with=F],
           measure.vars=shapcols,value.name="shap")

qcols<-paste0("q_", shapcols)
dlong<-melt(d,id.vars=setdiff(names(d),shapcols))[,-qcols,with=F]
dlong[,stdfvalue:=std1(value),by=variable]
dqlong<-melt(d,id.vars=setdiff(names(d),qcols),value.name="qvalue")[,-shapcols,with=F]
dqlong[,variable:=gsub("q_", "", variable)]
dlong<-merge(dlong, dqlong, by=intersect(names(dlong), names(dqlong)))

bigl<-merge(shaplong, dlong, by=c("fid", "fold", "variable"))
#shaplong<-melt(allshap,id.vars=setdiff(names(allshap),shapcols))
bigl<-merge(bigl, gcv, by="variable")
setnames(bigl, old="value", new="rfvalue")
bigl[,base_shap_id:=sum(shap),by=.(fid, base)]
bigl[,window_shap_id:=sum(shap),by=.(fid, window)]
bigl[,stat_shap_id:=sum(shap),by=.(fid, stat)]




#fwrite(bigl, file="big_plot_shap_long.tsv", sep="\t")

top20<-gcv[rnk < 21,as.character(variable)]
plotd<-bigl[variable %in% top20]

pfid<-c()
for(v in top20){
  setnames(d, old=v, new="v_tmp")
  pfid<-c(pfid, d[order(v_tmp)][seq(1,.N,length.out=500),fid])
  print(v) 
  setnames(d, old="v_tmp", new=v)
}
pfid<-unique(pfid)
pfid2<-allshap[order(flare)][seq(1,.N,length.out=15000-length(pfid)),fid]
pfid2<-setdiff(pfid2, pfid)
allshap[fid %in% pfid,quantile(flare,seq(0,1,.1))]
allshap[fid %in% c(pfid,pfid2),quantile(flare,seq(0,1,.1))]
allshap[,quantile(flare,seq(0,1,.1))]


#num_skip<-ceiling(allshap[,.N]/10^4)
#showid<-allshap[order(flare)][seq(1,.N,num_skip)][,fid]

showid<-c(pfid, pfid2)
fwrite(bigl[fid %in% showid],file="sub_plot_shap_long.tsv",sep="\t")

top20<-gcv[rnk < 21,as.character(variable)]
plotd<-bigl[fid %in% showid][variable %in% top20]

#fwrite(plotd,file="sub_plot_shap_long_top20.tsv",sep="\t")
library(SHAPforxgboost)
library(RColorBrewer)
plotd<-fread(file="sub_plot_shap_long_top20.tsv")
setnames(plotd,old=c("shap","global_shap"), new=c("value","mean_value"))
levs<-plotd[order(rnk),head(.SD,1),by=variable][,as.character(variable)]
plotd[,variable:=factor(variable, levels=levs)]
plotd[,rfvalue:=as.numeric(rfvalue)]
shap.plot.summary(data_long=plotd)
source("solar_flare/shap_plot_funcs.R")
cls<-brewer.pal(8, "Set1")[c(2,2,5,4,4)]
cls<-brewer.pal(8, "Dark2")[c(1,3,3, 3, 2)]
cls<-c(brewer.pal(8, "Paired")[3], brewer.pal(8, "Dark2")[c(3,3,3, 2)])
       
me.shap.plot.summary(data_long=plotd,cols=cls)

combos<-CJ(xvar=plotd[,unique(variable)], 
   yvar=plotd[,unique(variable)])[xvar!=yvar]
combos[,ind:=1:.N]
for(i in 1:combos[,.N]){
#for(i in (i+1):combos[,.N]){
  print(
  #shap.plot.dependence(data_long=plotd, 
  #		     x=as.character(combos[i,xvar]),
  #		     #y="shrgt45_min_t10_5",
  #		     color_feature=as.character(combos[i,yvar]))#,size0=2) #+ 
  me.dep.plot(plotd,combos[i,xvar],combos[i,yvar],quantx=T)
  )
  readline(prompt="Press [enter] to continue")

}
## 58



me.dep2.plot(plotd, 
	     x="totpot_max_t10_5", 
	     color_feature="epsx_min_t10_5", 
	     shape_feature="r_value_min_t10_5")

x="totpot_max_t10_5"; 
color_feature="epsx_min_t10_5"; 
shape_feature="r_value_max_t10_5"
shape_feature="absnjzh_sd_t10_5"
me.dep2.plot(plotd, x, color_feature, shape_feature)




















shap.plot.dependence(data_long=plotd, 
		     x="r_value_max_l1_t10_5",
		     #y="shrgt45_min_t10_5",
		     color_feature="shrgt45_min_t10_5")#,size0=2) #+ 
#shap.dependenc.plot(
plt_tst<-tstshap[order(flare)][seq(1,.N,18)]#[,.N] # 15
plt_cv<-allshap[order(flare)][seq(1,.N,20)]#[,.N]









