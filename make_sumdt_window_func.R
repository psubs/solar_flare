
make_sumdt<-function(dat, grpvars, allcovs){
  d<-copy(dat)
  
  fillin_func<-function(x) ifelse(is.infinite(x) | is.na(x), -99, x)
  d[,window:=paste0("t10_",floor(time/10))]

  for(nam in allcovs){
    setnames(d, old=nam, new="tmpvar")
    #sumdt <- d[,(paste0(nam, c("_sd", "_mean")))`:=`c(sd(tmpvar), mean(tmpvar)),by=grpvars]
    #print("cov set 1")
    if(d[,class(tmpvar)]=="character"){
     d[,tmpvar:=as.numeric(tmpvar)] 
    }

    d[,l1_shift:=tmpvar - shift(tmpvar, 1, type="lag"),by=grpvars]
    sumdt_window <- d[,.(tmpvar_sd=sd(tmpvar),
  		tmpvar_mean=mean(tmpvar),
  		tmpvar_q5 =quantile(tmpvar, .5), 
  		tmpvar_min =min(tmpvar), 
  		tmpvar_max =max(tmpvar),
  		tmpvar_range = max(tmpvar) - min(tmpvar)
  		),by=c(grpvars, "window")]
    
    sumdt_window<-sumdt_window[,lapply(.SD, fillin_func)] 
    
   
    lmcoef_window<-d[,.(tmpvar_intercept = lm(formula = tmpvar ~ time)$coefficients[1],
         tmpvar_slope = lm(formula = tmpvar ~ time)$coefficients[2],       
         tmpvar_slope_z = coef(summary(lm(formula = tmpvar ~ time)))[,"t value"][2]), 
       by = c(grpvars, "window")]
     
    lmcoef_window<-lmcoef_window[,lapply(.SD, fillin_func)] 
    

    sumdt_shift<-d[,.(tmpvar_sd_l1=sd(l1_shift,na.rm=T),
  		tmpvar_mean_l1=mean(l1_shift,na.rm=T),
  		tmpvar_q1_l1= quantile(l1_shift, .1,na.rm=T),
  		tmpvar_q9_l1= quantile(l1_shift, .9,na.rm=T),
  		tmpvar_q5_l1=quantile(l1_shift, .5,na.rm=T), 
  		tmpvar_min_l1=min(l1_shift,na.rm=T), 
  		tmpvar_max_l1=max(l1_shift,na.rm=T),
  		tmpvar_range_l1= max(l1_shift,na.rm=T) - min(l1_shift,na.rm=T)
  		),by=c(grpvars)]
    sumdt_shift<-sumdt_shift[,lapply(.SD, fillin_func)] 

    sumdt_shift_window<-d[,.(tmpvar_sd_l1=sd(l1_shift,na.rm=T),
  		tmpvar_mean_l1=mean(l1_shift,na.rm=T),
  		tmpvar_q5_l1=quantile(l1_shift, .5,na.rm=T), 
  		tmpvar_min_l1=min(l1_shift,na.rm=T), 
  		tmpvar_max_l1=max(l1_shift,na.rm=T),
  		tmpvar_range_l1= max(l1_shift,na.rm=T) - min(l1_shift,na.rm=T)
  		),by=c(grpvars, "window")]
    sumdt_shift_window<-sumdt_shift_window[,lapply(.SD, fillin_func)] 


    sumdt_window<-merge(sumdt_window, lmcoef_window, by=c(grpvars, "window"))
    sumdt_window<-merge(sumdt_window, sumdt_shift_window, by=c(grpvars, "window"))

    sumdt_w<-dcast(sumdt_window, id + fid + fold + grp + class_label ~ window, 
		   value.var=setdiff(names(sumdt_window), c("window", grpvars)))
     
    
    setnames(sumdt_w, 
	     old = names(sumdt_w), 
	     new=gsub("tmpvar", nam, names(sumdt_w))) 
    
    setnames(d, old="tmpvar", new=nam) 
    d[,l1_shift:=NULL]
#   if(!exists("sumdt_all")){}
   if(nam==allcovs[1]){
    sumdt_all <- sumdt_w 
   } else {
    sumdt_all <- merge(sumdt_all, sumdt_w, by=grpvars) 
   }
   print(nam)
#   fwrite(sumdt_all, file="../dan/sumdt_all.tsv", sep="\t")
  }
  return(sumdt_all)
}
