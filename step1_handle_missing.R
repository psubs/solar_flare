library(zoo)
library(xtable)

setwd("/gpfs/group/asb17/default/BD2019_solar/data")
for(i in 1:4){
  di<-fread(paste0("fold",i,"_all.tsv"))
  di[,fold:=i]
  di[,fid:=paste0("f",fold,"_",id)]
  if(i==1){
  d<-di
  } else {
  d<-rbindlist(list(d,di),
	       use.names=TRUE,fill=TRUE)
  }
  print(i)
}

#summaryf <-mean
d[,V1:=NULL]
onames <- copy(names(d))
setnames(d, old=onames, new=tolower(onames))
keys<-c("fold", "fid", "id")
resp<-c("class_label")
covs<-setdiff(names(d), c(keys, resp))

miss <- d[,lapply(.SD, function(x) mean(is.na(x)))]
dat<-copy(d)
rm(d)

impvars <- melt(miss, 
		id.vars="fold",
		variable.factor=FALSE)[variable!="class_label"][value > 0,unique(variable)]


imputeit<-function(dat, impvars){

d<-copy(dat)
d<-d[order(fold, id, fid, time)]
missdf<-list()
for(nam in impvars){
  setnames(d, old=nam, new="tmpvar")
  d[, value_imp := na.locf(tmpvar, na.rm = F), by = keys]
  d[, value_imp := na.locf(value_imp, fromLast = T), by = keys] #Fix when T=0 is missing
#  d[fid %in% d[is.na(value_imp),fid]][,c(keys, "time", "tmpvar", "value_imp"),with=F]
  if(d[is.na(value_imp),.N] > 0){
    impc <- d[,median(tmpvar,na.rm=T)]
    missdf[[nam]]<-d[is.na(tmpvar),c(keys, "time"),with=F][,`:=`(val=impc, var=paste0(nam))]
    d[is.na(value_imp),value_imp:=impc] 
  }
  d[,tmpvar:=value_imp]
  setnames(d, old="tmpvar", new=nam)
  d[,value_imp:=NULL]
  print(nam) 
 }
return(list(missdf=missdf, d=d))
}
dd<-imputeit(dat, impvars)
d<-dd$d
missdf<-dd$missdf
rm(dd)


missdt<-rbindlist(missdf)
missdt<-dcast(missdt, fold + fid + id + time ~ var, value.var="val")

missdt<-merge(missdt, d[,c(keys, "time", "class_label"),with=F],by=c(keys, "time"))
missdt[,.(sum(class_label), mean(class_label),uniqueN(fid)),by=fold]
print(xtable(rbindlist(list(
       d[,.(sum(class_label), as.character(round(mean(class_label,na.rm=T),3)), uniqueN(fid)),by=fold],
       d[,.(fold="overall", sum(class_label,na.rm=T), as.character(round(mean(class_label,na.rm=T),3)), uniqueN(fid))]
))),include.rownames=F)



d[,lapply(.SD, function(x) mean(is.na(x)))]
fwrite(d, file="../dan/imp_full.tsv", sep="\t")





