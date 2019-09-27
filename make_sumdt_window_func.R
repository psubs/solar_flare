

f1<-function(pred, label,cut=0.5){
        tp<-sum(label* (pred > cut))
        fp<-sum((1-label)* (pred > cut))
        fn<-sum(label* (pred <= cut))
        prec<-tp / (tp + fp)
        rec<-tp / (tp + fn)
        return(2*prec*rec / (prec + rec))
      }

prec<-function(pred, label,cut=0.5){
        tp<-sum(label* (pred > cut))
        fp<-sum((1-label)* (pred > cut))
        fn<-sum(label* (pred <= cut))
        prec<-tp / (tp + fp)
        return(prec)
      }

fpr<-function(pred, label,cut=0.5){
        fp<-sum((1-label)* (pred > cut))
        tn<-sum((1-label)* (pred <= cut))
        fpr<-fp / (fp + tn)
        return(fpr)
      }

tss<-function(pred, label,cut=0.5){
        tss<-rec(pred, label, cut) - fpr(pred, label, cut)
        return(tss)
      }

hss1<-function(pred, label, cut=0.5){
  hss1<-rec(pred, label, cut) *(2 - 1/prec(pred, label, cut))
  return(hss1) 
 }

gs<-function(pred, label, cut=0.5){ 
 p<-sum(label)
 n<-sum(1-label)
 fp<-sum((1-label)* (pred > cut))
 tn<-sum((1-label)* (pred <= cut))
 tp<-sum((label)* (pred > cut))
 fn<-sum(label*(pred <=cut))
 ch<-((tp + fp) * (tp + fn)) / (p + n) 
 gs<-(tp - ch) / (tp + fp + fn - ch)
 return(gs)
}
logit<-function(x){
 return(log(x / (1-x)))
}

invlogit<-function(x){
 return(exp(x) / (1 + exp(x)))
}

rec<-function(pred, label,cut=0.5){
        tp<-sum(label* (pred > cut))
        fp<-sum((1-label)* (pred > cut))
        fn<-sum(label* (pred <= cut))
        rec<-tp / (tp + fn)
        return(rec)
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
