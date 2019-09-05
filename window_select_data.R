
library(data.table)
source("solar_flare/make_sumdt_window_func.R")
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


datafiles<-paste0("data/ssdat_window_", 1:369)
for(i in 1:369){
  di<-fread(datafiles[i], select=c(grpvars, usevars_window))
  if(i==1){
  d<-di
  } else {
  d<-rbindlist(list(d,di))
  }
 print(i)
}

fwrite(d, file="window_select_data.tsv", sep="\t")
