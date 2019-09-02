

setwd("/gpfs/group/asb17/default/BD2019_solar/dan")
source('make_sumdt_window_func.R')
args=commandArgs(trailingOnly=T)

i<-as.numeric(args[1])
print("reading data...")
dat<-fread("imp_full.tsv", sep="\t")

n_fid<-dat[,length(unique(fid))]
k<-1000
n_grp<-floor(n_fid  / k)
grpid<-rep(1:n_grp, times=c(rep(k, n_grp-1), n_fid - k*(n_grp-1)))
chk<-dat[,.(fid=unique(fid))]
chk[,grp:=grpid]
dat<-merge(dat, chk, by=c("fid"))

grpvars <- c("id", "fid", "fold", "grp", "class_label")
allcovs <- setdiff(names(dat), c("time", grpvars))

print("making featurized ts data...")
di<-make_sumdt(dat[grp==i], grpvars, allcovs)

fwrite(di, file=paste0("data/ssdat_window_",i),sep="\t")

#for(i in 1:length(unique(dat[,grp]))){
#  di<-make_sumdt(dat[grp==i], grpvars, allcovs)
#  fwrite(di, file=paste0("../dan/data/ssdat_",i),sep="\t")
#  print(i)
#  rm(di)
#  print(i / length(unique(dat[,grp])))
#}

