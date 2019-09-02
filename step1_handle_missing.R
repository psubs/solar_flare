library(zoo)

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

#nam<-impvars[1]

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
    #impc <- d[,median(tmpvar,na.rm=T)]
    impc <- -9 
    #missdf[[nam]]<-d[is.na(tmpvar),c(keys, "time"),with=F][,`:=`(val=impc, var=paste0(nam))]
    d[is.na(value_imp),value_imp:=impc] 
  }
  d[,tmpvar:=value_imp]
#  impc <- merge(impc, d[is.na(tmpvar),.(fid,time)], by="fid")
#  d <- merge(d, impc, by=c("fid","time"),all.x=TRUE)
#  d[is.na(tmpvar),tmpvar:=m]
#  d[,`:=`(m=NULL, s=NULL)]
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

d[,lapply(.SD, function(x) mean(is.na(x)))]
# fwrite(d, file="../dan/imp_full.tsv", sep="\t")
dat<-fread("../dan/imp_full.tsv", sep="\t")

# renans code
#    data[, value_imp := na.locf(value, na.rm = F), by = byvars]
#    data[, value_imp := na.locf(value_imp, fromLast = T), by = byvars] #Fix when T=0 is missing
#    # 3. Compute summary statistics
#    data[, vmax := max(value_imp), by = byvars]
#    data[, vmin := min(value_imp), by = byvars]
#    data[, vmed := median(value_imp), by = byvars]
#    data[, vmean := mean(value_imp), by = byvars]
#    data[, vvar := var(value_imp), by = byvars]
#    # 4. Fit linear model for each id/variable
#    data[, intercept := glm(formula = value_imp ~ time)$coefficients[1], by = byvars]
#    data[, slope := glm(formula = value_imp ~ time)$coefficients[2], by = byvars]

#### standardize ###
#d[,(allcovs):=lapply(.SD, function(x) c(scale(x,center=TRUE,scale=TRUE))),.SDcols=allcovs]


### covariate-wise sumstats ###

make_sumdt<-function(dat, grpvars, allcovs){
  d<-copy(dat)
  
  fillin_func<-function(x) ifelse(is.infinite(x) | is.na(x), -99, x)

  for(nam in allcovs){
    setnames(d, old=nam, new="tmpvar")
    #sumdt <- d[,(paste0(nam, c("_sd", "_mean")))`:=`c(sd(tmpvar), mean(tmpvar)),by=grpvars]
    #print("cov set 1")
    if(d[,class(tmpvar)]=="character"){
     d[,tmpvar:=as.numeric(tmpvar)] 
    }
    sumdt <- d[,.(tmpvar_sd=sd(tmpvar),
  		tmpvar_mean=mean(tmpvar),
  		tmpvar_q1 = quantile(tmpvar, .1),
  		tmpvar_q9 = quantile(tmpvar, .9),
  		tmpvar_q5 =quantile(tmpvar, .5), 
  		tmpvar_min =min(tmpvar), 
  		tmpvar_max =max(tmpvar),
  		tmpvar_range = max(tmpvar) - min(tmpvar)
  		),by=grpvars]
    
    sumdt[,`:=`(tmpvar_q9overq1 = tmpvar_q9/tmpvar_q1,
  #	      tmpvar_q5overq1 = tmpvar_q5/tmpvar_q1,
  #	      tmpvar_q9overq5 = tmpvar_q9/tmpvar_q5,
                tmpvar_mean_over_sd = tmpvar_mean / tmpvar_sd,
  #              tmpvar_qrange_over_range = (tmpvar_q9 - tmpvar_q1) / tmpvar_range,
                tmpvar_qrange_over_sd = (tmpvar_q9 - tmpvar_q1) / tmpvar_sd
               )] 
    sumdt<-sumdt[,lapply(.SD, fillin_func)] 
    
    sumdt_l5 <- d[,tail(.SD, 5), by=grpvars]
    sumdt_l5 <- sumdt_l5[,.(tmpvar_sd=sd(tmpvar),
  		tmpvar_mean=mean(tmpvar),
 # 		tmpvar_q1 = quantile(tmpvar, .1),
 # 		tmpvar_q9 = quantile(tmpvar, .9),
 # 		tmpvar_q5 =quantile(tmpvar, .5), 
  		tmpvar_min =min(tmpvar), 
  		tmpvar_max =max(tmpvar),
  		tmpvar_range = max(tmpvar) - min(tmpvar)
  		),by=grpvars]
    sumdt_l5[,`:=`(
  #	       tmpvar_q9overq1 = tmpvar_q9/tmpvar_q1,
  #            tmpvar_q5overq1 = tmpvar_q5/tmpvar_q1,
  #            tmpvar_q9overq5 = tmpvar_q9/tmpvar_q5,
                tmpvar_mean_over_sd = tmpvar_mean / tmpvar_sd#,
  #              tmpvar_qrange_over_range = (tmpvar_q9 - tmpvar_q1) / tmpvar_range,
  #              tmpvar_qrange_over_sd = (tmpvar_q9 - tmpvar_q1) / tmpvar_sd
               )] 
    sumdt_l5<-sumdt_l5[,lapply(.SD, fillin_func)] 

    sumdt_f5 <- d[,head(.SD, 5), by=grpvars]
    sumdt_f5 <- sumdt_f5[,.(tmpvar_sd=sd(tmpvar),
  		tmpvar_mean=mean(tmpvar),
 # 		tmpvar_q1 = quantile(tmpvar, .1),
 # 		tmpvar_q9 = quantile(tmpvar, .9),
 # 		tmpvar_q5 =quantile(tmpvar, .5), 
  		tmpvar_min =min(tmpvar), 
  		tmpvar_max =max(tmpvar),
  		tmpvar_range = max(tmpvar) - min(tmpvar)
  		),by=grpvars]
    sumdt_f5[,`:=`(
  #	       tmpvar_q9overq1 = tmpvar_q9/tmpvar_q1,
  #            tmpvar_q5overq1 = tmpvar_q5/tmpvar_q1,
  #            tmpvar_q9overq5 = tmpvar_q9/tmpvar_q5,
                tmpvar_mean_over_sd = tmpvar_mean / tmpvar_sd#,
  #              tmpvar_qrange_over_range = (tmpvar_q9 - tmpvar_q1) / tmpvar_range,
  #              tmpvar_qrange_over_sd = (tmpvar_q9 - tmpvar_q1) / tmpvar_sd
               )] 

    sumdt_f5<-sumdt_f5[,lapply(.SD, fillin_func)] 
    
    m<-melt(sumdt, id.vars=grpvars)    
    m5<-melt(sumdt_l5, id.vars=grpvars)    
    mf5<-melt(sumdt_f5, id.vars=grpvars)    
    m<-merge(m, m5, by=c(grpvars, "variable"), suffixes=c("", "_l5"),all=T)
    m<-merge(m, mf5, by=c(grpvars, "variable"), suffixes=c("", "_f5"),all=T)

    m[!is.na(value_l5),all_over_l5:=value/value_l5] 
    m[!is.na(value_f5),all_over_f5:=value/value_f5] 
    m[!is.na(value_f5) & !is.na(value_l5),l5_over_f5:=value_l5/value_f5] 
    
   # ss<-m[,lapply(.SD , median),
   #       by=.(variable,class_label),
   #       .SDcols=grep("value", names(m),value=T)]
   # dcast(ss, variable ~class_label, value.vars=grep("value", names(ss),value=T))
    #m[,lapply(.SD ,function(x) mean(is.na(x))),
    #  by="variable",
    #  .SDcols=grep("value|_over_", names(m),value=T)]
    mw<-dcast(m, class_label + fold + grp + id + fid ~ variable, 
	      value.var=grep("value|_over_", names(m),value=T))
    rmvars<-melt(mw[,lapply(.SD ,function(x) mean(is.na(x))),
	 .SDcols=grep("tmpvar", names(mw),value=T)][,
	 var:="var"],
    id.vars="var",variable.factor=F)[value==1,variable]
    
    mw<-mw[,-c(rmvars),with=F]
    mw<- mw[,lapply(.SD, fillin_func)]
   
    ### time-based ### 

    maxt<-d[,.SD[which.max(tmpvar)],by=grpvars][,.(fid, tmpvar_maxtime=time)]   
    mint<-d[,.SD[which.min(tmpvar)],by=grpvars][,.(fid, tmpvar_mintime=time)]   

    d[,`:=`(l1_tmp=tmpvar - shift(tmpvar, 1, type="lag")#,
           # l10_tmp=tmpvar - shift(tmpvar, 10, type="lag")
	    ),
        by=grpvars]
   

    tb1l<-d[,.SD[which.max(l1_tmp)],by=grpvars][,
	.(fid, tmpvar_maxlag1=l1_tmp,tmpvar_maxlag1t=time)]
    
    tbdiff<-dcast(d[time==0 | time==59], fid ~ time, value.var="tmpvar")[,
	.(fid, tmpvar_last_minus_first=`59`-`0`)]

    tb<-merge(maxt, mint, by="fid")
    tb<-merge(tb, tb1l, by="fid")
    tb<-merge(tb, tbdiff, by="fid") 
   
    tb<-tb[,lapply(.SD, fillin_func)] 
    ### regression slope / int  and z-score ## 
    lmcoef<-d[,.(tmpvar_intercept = lm(formula = tmpvar ~ time)$coefficients[1],
         tmpvar_slope = lm(formula = tmpvar ~ time)$coefficients[2],       
         tmpvar_slope_z = coef(summary(lm(formula = tmpvar ~ time)))[,"t value"][2]), 
       by = grpvars]
     
    lmcoef<-lmcoef[,lapply(.SD, fillin_func)] 
    
    mw<-merge(mw, lmcoef, by=grpvars)
    mw<-merge(mw, tb, by="fid")
    
    setnames(mw, 
	     old = names(mw), 
	     new=gsub("tmpvar", nam, names(mw))) 
    setnames(mw, 
	     old = names(mw), 
	     new=gsub("^value_", "", names(mw))) 
    
    setnames(d, old="tmpvar", new=nam) 
    d[,l1_tmp:=NULL]
#   if(!exists("sumdt_all")){}
   if(nam==allcovs[1]){
    sumdt_all <- mw 
   } else {
    sumdt_all <- merge(sumdt_all, mw, by=grpvars) 
   }
   print(nam)
#   fwrite(sumdt_all, file="../dan/sumdt_all.tsv", sep="\t")
  }
  return(sumdt_all)
}






n_fid<-dat[,length(unique(fid))]
k<-1000
n_grp<-floor(n_fid  / k)
grpid<-rep(1:n_grp, times=c(rep(k, n_grp-1), n_fid - k*(n_grp-1)))
chk<-dat[,.(fid=unique(fid))]
chk[,grp:=grpid]
dat<-merge(dat, chk, by=c("fid"))

grpvars <- c("id", "fid", "fold", "grp", "class_label")
allcovs <- setdiff(names(dat), c("time", grpvars))

for(i in 1:length(unique(dat[,grp]))){
  di<-make_sumdt(dat[grp==i], grpvars, allcovs)
  fwrite(di, file=paste0("../dan/data/ssdat_",i),sep="\t")
  print(i)
  rm(di)
  print(i / length(unique(dat[,grp])))
}

i<-1
di<-dat[grp==1]
dii<-make_sumdt(di, grpvars, allcovs)





sapply(1:length(allcovs), function(i) length(grep(allcovs[i], names(sumdt_all))))



### sample covariates to mix datasets ###
mod_names<-lapply(1:10, function(i) c())
mod_names<-lapply(1:10, function(i) c(grpvars))
for(i in 1:length(allcovs)){
        av<-grep(allcovs[i], names(sumdt_all),value=T)
        assgn<-sample(1:10, length(av), replace=T)
	naml<-lapply(1:10, function(i) av[which(assgn==i)])
	naml2<-lapply(1:10, function(i){
          if(length(naml[[i]]) < 6){
	    c(naml[[i]], sample(setdiff(av, naml[[i]]), 6 - length(naml[[i]])))
	  } else {
	   naml[[i]] 
	  }
        })
      mod_names<-lapply(1:10, function(i) c(mod_names[[i]], naml2[[i]]))
}








sumstat_covs <-setdiff(names(sumdt_all), grpvars)
nuniq <- melt(sumdt_all[,lapply(.SD, uniqueN),.SDcols=sumstat_covs][, id:="id"], id.vars="id")
nzv <- melt(sumdt_all[,lapply(.SD, sd),.SDcols=sumstat_covs][,id:="id"], 
	    id.vars="id")


### standardize the covariate-wise sumstats
sumdt_all[,(sumstat_covs):=lapply(.SD, function(x) 
	   c(scale(x, center=TRUE, scale=TRUE))), 
	  .SDcols=sumstat_covs]

sumdt_all1 <-sumdt_all[,-c(grep("mean_over_sd", names(sumdt_all), value=TRUE)),with=FALSE]
sumdt_all1[,fold:=tstrsplit(fid, "_")[[1]]]
 melt(sumdt_all1[1:3], id.vars=grpvars)[is.na(value)]#[40]#[,tstrsplit(variable, "_")[[3]]]
for(i in 1:4){
  fwrite(sumdt_all1[fold==paste0("f",i)], file=paste0("fold",i,"_sumstatcovs.tsv"),sep="\t")
  print(i)
}


for(i in 1:4){
  di<-fread(paste0("fold", i, "_sumstatcovs.tsv"))
  if(i==1){
  d<-di
  } else {
  d<-rbindlist(list(d,di))
  }
  print(i)
}
dtw <- d; rm(d)




nzv[is.na(value) | value < 1e-2][,variable]
d[,.(sum(class_label)),by=fold]

### look at difference between last few observations (closest to solar event ) ####
### and the mean, median, q9, q1, min, max ##
## which time has max(x), which has min(x) ## 
### time based summary stats ####

d<-d[order(fid, time)]
for(nam in allcovs){
  setnames(d, old=nam, new="tmpvar")
  d[time >=50, .(tmpvar_l10_sd=sd(tmpvar),
                 tmpvar_l10_mean=mean(tmpvar),
		 tmpvar_l10_max=quantile(tmpvar, 0.1),
		 tmpvar_l10_min=quantile(tmpvar, 0.1),
  d[time >=55, .(tmpvar_l5_sd=sd(tmpvar),
                 tmpvar_l5_mean=mean(tmpvar),
		 tmpvar_l5_max=quantile(tmpvar, 0.1),
		 tmpvar_l5_min=quantile(tmpvar, 0.1),
  d[time >59, .(tmpvar_l1_sd=sd(tmpvar),
                 tmpvar_l1_mean=mean(tmpvar),
		 tmpvar_l1_max=quantile(tmpvar, 0.1),
		 tmpvar_l1_min=quantile(tmpvar, 0.1),
  d[time ==0, .(tmpvar_l1_sd=sd(tmpvar),
                 tmpvar_l1_mean=mean(tmpvar),
		 tmpvar_l1_max=quantile(tmpvar, 0.1),
		 tmpvar_l1_min=quantile(tmpvar, 0.1),
  setnames(sumdt, old = names(sumdt), new=gsub("tmpvar", nam, names(sumdt))) 
  setnames(d, old = "tmpvar", new=nam) 
  
}
    


### covariate interaction summary stats #### 




### ####

p <-fread("tgpreds.csv")
f1 <- function(pred, label){
  tp <- sum((pred > 0.5) * label)
  fp <- sum((pred > 0.5) * (1-label))
  fn <- sum((pred < 0.5) * label)
  prec <- tp/(tp + fp)
  rec <- tp/(tp + fn)
  2*prec*rec / (prec + rec)

}
p[,fold:=tstrsplit(fid, "_")[[1]]]
f1r<-p[,.(f1(pred, label)), by=tgi][order(-V1)]
f1r<-p[tgi==18,.(f1(pred, label)), by=fold]



              tmpvar_q9over91 = tmpvar_q9/tmpvar_q1,
              tmpvar_q9over91 = tmpvar_q9/tmpvar_q1,

  sumdt <- d[,.(tmpvar_sd=sd(tmpvar)),by=grpvars]
  sumdt <- d[,.(tmpvar_mean=mean(tmpvar)),by=grpvars]
  sumdt <- d[,.(tmpvar_sd=sd(tmpvar)),by=grpvars]
  sumdt <- d[,.(tmpvar_sd=sd(tmpvar)),by=grpvars]
  sumdt <- d[,.(tmpvar_sd=sd(tmpvar)),by=grpvars]
  sumdt <- d[,.(tmpvar_sd=sd(tmpvar)),by=grpvars]
   
    
}


for(lag in c(1:10, 20, 30, 40, 50, 60)){
  for(nam in allcovars){
  setnames(d, old=nam, new="tmpvar")
  d[,tmpvar]
   
    
  }
}



























#ds = d[fold==1]
#ds[,(paste0(covs, "_ct")):=lapply(.SD, function(x){ 
#				     cut(x, breaks=10, include.lowest=T)}),
#	.SDcols=covs]
#grpvars <-"time_ct" 
#
#ds[,lapply(.SD, function(x){quantile(x, seq(0,1,.1))})]
#quantmap<-ds[,
#	     lapply(.SD, function(x){quantile(x, seq(0,1,.1),na.rm=T)}),
#	  .SDcols=covs]
#ds2<-d[fold==2]
#for(i in names(quantmap)){
#  print(i)
#  setnames(ds2, old=i, new="tmpvar")
#  ds2[,(paste0(i, "_ct")):=cut(tmpvar, breaks=quantmap[[i]], include.lowest=T)]
#  setnames(ds2, new=i, old=c("tmpvar"))
#}
#catmap_fit<-function(d, f, grpvars,sdcols){
#  ds[,lapply(.SD, f, na.rm=T), .SDcols=sdcols, by=grpvars] 
#
#}
#
#
#library(matrixStats)
#
#d[,V1:=NULL]
#covs <- setdiff(names(d), c("fold", "time", "fid", "class_label", "id"))
#chk<-d[,covs,with=FALSE][,lapply(.SD, function(x) as.numeric(is.na(x) | x==0))] 
#chk[,`:=`(chksum=rowSums(.SD, na.rm=TRUE)),.SDcols=covs]
#chk3<-chk[chksum==3,which=TRUE]
#chk2<-chk[chksum==2,which=TRUE]
#chk23<-chk[chksum==23,which=TRUE]
#chk24<-chk[chksum==24,which=TRUE]
#chk25<-chk[chksum==25,which=TRUE]
#
#
#chk9<-d[,covs,with=FALSE][,lapply(.SD, function(x) x==-9999)] 
#chk9[,`:=`(chksum=rowSums(.SD, na.rm=TRUE)),.SDcols=covs]
#d[chk3]
#d[chk]
#d[chk3]
#
#
#
#
