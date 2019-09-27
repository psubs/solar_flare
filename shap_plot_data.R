allcovs<-fread("ieee/allcovsbase.tsv")[,unlist(covs)]
pat<-paste0(allcovs,collapse="|")



global_shap<-function(shapdt, shapcols,abs=TRUE,byvars=NULL){
 
  if(abs){
    aggf<-function(x){ mean(abs(x)) }
  } else{
    aggf<-function(x){ mean(x) }
  }
  
 if(!missing(byvars)){
          global<-melt(shapdt[,lapply(.SD, function(x) aggf(x)),
                       .SDcols=shapcols,by=byvars][,id:="var"], 
                	value.name="global_shap",
                	id.vars=c("id", byvars))[,id:=NULL][order(-global_shap)]
	}else{
          global<-melt(shapdt[,lapply(.SD, function(x) aggf(x)),
		       .SDcols=shapcols][,id:="var"], 
        		value.name="global_shap",
        		id.vars="id")[,id:=NULL][order(-global_shap)]
        	
	}
  return(global)
}



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

#global_shap<-function(shapdt, shapcols){
#  global<-melt(shapdt[,lapply(.SD, function(x) mean(abs(x))),.SDcols=shapcols][,id:="var"], 
#		value.name="global_shap",
#		id.vars="id")[,id:=NULL][order(-global_shap)]
#  return(global)
#}

source("solar_flare/make_sumdt_window_func.R")
library(ggplot2)
library(data.table)
library(stringr)
library(xgboost)
library(parallelDist)
library(stringr)
## data read-in ##

cvshap<-fread("ieee/cvshap.tsv")
tstshap<-fread("ieee/testpred_shap.tsv")
setnames(tstshap, old="Id", new="id")
tstshap<-tstshap[,`:=`(fold=4)][,fid:=paste0("f",fold,"_",id)]
allshap<-rbindlist(list(cvshap,tstshap),use.names=T,fill=T)

cvcols<-c(names(allshap)[1:7], "BIAS", "lev", "ClassLabel", "Id")
shapcols<-setdiff(names(allshap), cvcols)


d<-fread("window_select_data.tsv")
d[,class_label:=factor(class_label, levels=c("0", "1"), 
		       	labels=c("no_flare", "big_flare"))]
grpvars <- c("id", "fid", "fold", "grp", "class_label")
cnvrt_num<-melt(d[,lapply(.SD, class)][,id:="id"],
		   id.vars="id",
		   variable.factor=F)[value=="character"][!variable %in% grpvars,variable]

d[,paste(cnvrt_num):=lapply(.SD, as.numeric),.SDcols=cnvrt_num]
d[,(paste0("q_", shapcols)):=lapply(.SD, function(x) rank(x) / .N),.SDcols=shapcols]



gcv<-global_shap(allshap, shapcols)
#gtst<-global_shap(tstshap, shapcols)

gcv<-decomp(gcv)[,rnk:=1:.N]
#gtst<-decomp(#gtst)[,rnk:=1:.N]

gcv[,.(list(unique(variable))),by=.(base,window)]
gcv[grep("l1", variable),vtype:="var"]
gcv[grep("sd|range|slope", variable),vtype:="var"]
gcv[is.na(vtype),vtype:="loc"]


#gcv[,base_shap_sum:=sum(abs(global_shap)),by=base]
#gcv[,window_shap_sum:=sum(abs(global_shap)),by=window]
#gcv[,stat_shap_sum:=sum(abs(global_shap)),by=stat]
#gcv[,.(Nvar=uniqueN(variable)),by=.(base, base_shap_sum)][order(-base_shap_sum)]
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
bigl[,base_window_shap_id:=sum(shap),by=.(fid, base, window)]
bigl[,vtype_shap_id:=sum(shap),by=.(fid, vtype)]
bigl[,base_vtype_shap_id:=sum(shap),by=.(fid, base, vtype)]




bigl[,global_base_shap:=mean(abs(base_shap_id)),by=base]
bigl[,global_window_shap:=mean(abs(window_shap_id)),by=window]
bigl[,global_stat_shap:=mean(abs(stat_shap_id)),by=stat]
bigl[,global_base_window_shap:=mean(abs(base_window_shap_id)),by=.(base,window)]
bigl[,global_vtype_shap:=mean(abs(vtype_shap_id)),by=.(vtype)]
bigl[,global_base_vtype_shap:=mean(abs(base_vtype_shap_id)),by=.(base,vtype)]


chk<-bigl[,.(fid, base_shap_sum, base_shap_id, base)][,unique(.SD)]
chk[,.(base_shap_sum, mean(abs(base_shap_id))),by=base]

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
source("solar_flare/shap_plot_funcs.R")

plotd<-fread(file="sub_plot_shap_long.tsv")
setnames(plotd,old=c("shap","global_shap"), new=c("value","mean_value"))
levs<-plotd[order(rnk),head(.SD,1),by=variable][,as.character(variable)]
plotd[,variable:=factor(variable, levels=levs)]
plotd[,rfvalue:=as.numeric(rfvalue)]
#shap.plot.summary(data_long=plotd)
#cls<-brewer.pal(8, "Set1")[c(2,2,5,4,4)]
#cls<-brewer.pal(8, "Dark2")[c(1,3,3, 3, 2)]
cls<-c(brewer.pal(8, "Paired")[3], brewer.pal(8, "Dark2")[c(3,3,3, 2)])



#plotdw<-plotd[,.(fid, base, base_shap_id, global_base_shap)][,unique(.SD)]
#plotds<-plotd[,.(fid, base, base_shap_id, global_base_shap)][,unique(.SD)]

plotdb<-plotd[,.(fid, base, base_shap_id, global_base_shap)][,unique(.SD)]
blevs<-plotdb[order(-global_base_shap),head(.SD,1),by=base][,base]
plotdb[,base:=factor(base, levels=blevs)]

base_g<-me.shap.plot.summary(dat=plotdb,
		     cols=cls,
		     variable_name="base",
		     shap_name="base_shap_id", 
		     global_shap_name="global_base_shap",
		     show_color=F)
base_g + 
	annotate("text",
                x = 11 + 0.25, y=-3.9,label="Global Importance",
                size = 3.1, alpha = 0.7,
                hjust = -0.2,
                fontface = "bold")  + 
        labs(title="Aggregated SHAP by SHARP Parameter") + 
        theme(plot.title=element_text(hjust=0.5,size=12))
ggsave("global_shap_base.png",width=3,height=4.2,dpi=300)


plotdb<-plotd[,.(fid, window, window_shap_id, global_window_shap)][,unique(.SD)]
blevs<-plotdb[order(-global_window_shap),head(.SD,1),by=window][,window]
wlabs<-c("w5: 50-59", "w4: 40-49", "w3: 30-39", "w2: 20-29", "w0: time 0-9", "w1: time 10-19")
plotdb[,window:=factor(window, levels=blevs,labels=wlabs)]

window_g<-me.shap.plot.summary(dat=plotdb,
		     cols=cls,
		     variable_name="window",
		     shap_name="window_shap_id", 
		     global_shap_name="global_window_shap",
		     show_color=F)
window_g + 
	annotate("text",
                x = 6 + 0.25, y=-7,label="Global Importance",
                size = 3.1, alpha = 0.7,
                hjust = -0.2,
                fontface = "bold")  + 
        labs(title="Aggregated SHAP by Time-Series Window") + 
        theme(plot.title=element_text(hjust=0.5,size=12))
ggsave("global_shap_window.png",width=4,height=4,dpi=300)


plotdb<-plotd[,.(fid, bw=paste0(base, "_w", window), base_window_shap_id, global_base_window_shap)][,unique(.SD)]
blevs<-plotdb[order(-global_base_window_shap),head(.SD,1),by=bw][,bw]
#wlabs<-c("time 50-59", "time 40-49", "time 30-39", "time 20-29", "time 0-9", "time 10-19")
plotdb[,bw:=factor(bw, levels=blevs)]

bw_g<-me.shap.plot.summary(dat=plotdb,
		     cols=cls,
		     variable_name="bw",
		     shap_name="base_window_shap_id", 
		     global_shap_name="global_base_window_shap",
		     show_color=F)
bw_g + 
	annotate("text",
                x = 6 + 0.25, y=-7,label="Global Importance",
                size = 3.1, alpha = 0.7,
                hjust = -0.2,
                fontface = "bold")  + 
        labs(title="Aggregated SHAP by \nSHARP Parameter and Time-Series Window") + 
        theme(plot.title=element_text(hjust=0.5,size=12))
ggsave("global_base_shap_window.png",width=4,height=4,dpi=300)


plotdb<-plotd[,.(fid, vtype, base, bv=paste0(base, "_", vtype), base_vtype_shap_id, global_base_vtype_shap)][,unique(.SD)]
blevs<-plotdb[order(-global_base_vtype_shap),head(.SD,1),by=bv][,bv]
#wlabs<-c("time 50-59", "time 40-49", "time 30-39", "time 20-29", "time 0-9", "time 10-19")
plotdb[,bv:=factor(bv, levels=blevs)]

bv_g<-me.shap.plot.summary(dat=plotdb,
		     cols=cls,
		     variable_name="bv",
		     shap_name="base_vtype_shap_id", 
		     global_shap_name="global_base_vtype_shap",
		     show_color=F) 
bv_g + 
	annotate("text",
                x = 6 + 0.25, y=-7,label="Global Importance",
                size = 3.1, alpha = 0.7,
                hjust = -0.2,
                fontface = "bold")  + 
        labs(title="Aggregated SHAP by \nSHARP Parameter and Variability/Location Type of Derived Statistic") + 
        theme(plot.title=element_text(hjust=0.5,size=12))
ggsave("global_base_shap_vtype.png",width=4,height=4,dpi=300)


combos<-CJ(xvar=plotd[,unique(variable)], 
   yvar=plotd[,unique(variable)])[xvar!=yvar]
combos[,ind:=1:.N]
combos<-combos[grep("r_value_mean_", xvar)][grep("r_value_mean_",yvar)]
#combos<-combos[!grep("r_value", xvar)][!grep("r_value",yvar)]
#combos<-combos[xvar=="totpot_max_t10_5"]

dw<-dcast(plotd,fid + fold + flare + obs + id + grp ~ variable, value.var="value")
kk<-list()
for(k in 2:15){
 kk[[k]]<-kmeans(dw[,-c("fid", "fold", "flare", "obs", "id", "grp"),with=F],centers=k) 
 print(kk[[k]]$betweenss / kk[[k]]$totss)
}

dw[,cluster:=kk[[13]]$cluster]
dw[,.(mean(flare),.N),by=.(cluster,fold)][order(-V1)]

clusters <- hclust(dist(scale(as.matrix(dw[,-c("fid", "fold", "flare", "obs", "id", "grp"),with=F]))), method = "ward.D")
kcl<-cutree(clusters, k=10)
dw[,hclust10:=kcl]
dw<-dw[,clusterid:=clusters$order][rank(clusterid),]
dw[,clusterid:=NULL]
dw[,id:=.I]
dw[,mflare:=mean(flare),by=hclust10]
dws<-dw[,.(gsub("NA", "-", paste0(round(mean(obs),2)),
#				  "(", .N, ")")
)),
	by=.(mflare,hclust10,fold)]
dws<-dcast(dws, mflare + hclust10 ~ fold, value.var=c("V1"))

ch<-melt(dw, id.vars=c("fid", "fold", "flare", "obs", "hclust10"),measure.vars=plotd[,unique(as.character(variable))])
library(class)

gcv_clust<-global_shap(dw, shapcols=plotd[,as.character(unique(variable))], 
		       byvars=c("mflare", "hclust10"))[order(hclust10,-global_shap)]

gcv_clust[,rnk:=1:.N,by=.(hclust10)]

allcovs<-fread("allcovsbase.tsv")[,unlist(covs)]
pat<-paste0(allcovs,collapse="|")
gcv_clust<-decomp(gcv_clust)

plotd<-merge(dw[,.(fid,mflare,hclust10)],plotd,by="fid")
bases<-dcast(plotd[,.(fid,base,hclust10,mflare,base_shap_id)][,unique(.SD)], fid + hclust10+mflare ~ base, value.var="base_shap_id")
gcvb<-global_shap(bases, shapcols=intersect(names(bases),allcovs),
		  byvars=c("hclust10","mflare"))

gcvb<-global_shap(bases, shapcols=intersect(names(bases),allcovs),abs=F,
		  byvars=c("hclust10","mflare"))

blevs<-gcvb[order(-mflare,-global_shap),unique(as.character(variable))]
gcvb[,variable:=factor(variable, levels=blevs)]
clevs<-gcvb[order(-mflare)][,unique(hclust10)]
clabs<-gcvb[order(-mflare)][,as.character(round(unique(mflare),2))]
clabs<-paste(paste0("Cluster", 1:10, ":"), clabs)

gcvb[,hclust10:=factor(hclust10,levels=clevs,labels=clabs)]
clustcols<-brewer.pal(11,"Paired")
ggplot(gcvb, 
       aes(x=variable, y=global_shap,fill=hclust10)) + 
  geom_bar(stat="identity",position="dodge") + 
#  facet_wrap(~hclust10,scales="free_x") + 
  scale_fill_manual(values=clustcols[1:10],
		    #breaks=clevs,
		    labels=clabs,
		    name="Cluster ID:\nAvg Predicted Probability \nof MSF\nwithin Cluster") + 
  theme_bw() + 
  labs(x="SHARP Parameter", 
       title="Avg Magnitude of SHARP on predicted log-odds of MSF",
       y="Avg Magnitude of SHARP" )
ggsave("subcluster_breakdownBdir.png")
ggplot(gcvb, 
       aes(x=hclust10, y=global_shap,fill=variable)) + 
  geom_bar(stat="identity",position="dodge") + 
#  facet_wrap(~hclust10,scales="free_x") + 
  scale_fill_manual(values=clustcols,name="SHARP parameter") + theme_bw() + 
  labs(x="Cluster ID: (Avg Predicted Probability of Flare within Cluster)",
       title="Avg influence of SHARP on predicted log-odds of MSF by Cluster",
       y="Avg influence of SHARP" )
ggsave("subcluster_breakdown_dir.png")
ggsave("subcluster_breakdown_dir.pdf")


plot(x=0:59,y=rep(0,60),pch="",ylim=c(-.3,.3),
     xaxt='n',yaxt='n', ann=F,axes=F)
rect(xleft=seq(0,59,1), 
     ybottom=rep(-.1, 60),
     xright=seq(0,59,1)+1,
     ytop=rep(0.1,60),
     col=rep(clustcols[1:6],each=10),ylab="")
text(x=seq(5, 55, 10),y=rep(-0.2,10),labels=paste0("window",0:5),cex=0.8)
text(x=c(seq(0.5, 50.5,10),59.5),y=rep(-0.2,10),
     labels=paste0("t=",c(seq(0,50,10),59)),
		   cex=0.8)
points(x=c(seq(0.5, 50.5,10),59.5), y=
points(x=0:59,y=rep(0,60),pch="|")
abline(h=0)
text(x=
#text(x = 0.5, y = 0.5, '{', srt = 90, cex = 8, family = 'Helvetica Neue UltraLight')

plotd[,.(list(unique(stat)),uniqueN(variable)),by=base][order(-V2)]
plotd[,.(list(unique(paste0(stat, "_w",window))),uniqueN(variable)),by=base][order(-V2)]



shapobs_long<-merge(dw[mflare > 0.1,.(cid=1:.N,hclust10,mflare,fid)], 
                    plotd[,.(fid,base,base_shap_id)][,unique(.SD)],
                    by=c("fid"))
shapobs_long<-shapobs_long[order(cid)]
id<-"cid"
value<-"base_shap_id"
variable<-"base"
zoom_in_location = NULL;
y_parent_limit = NULL;
y_zoomin_limit = NULL;   # c(a;b) to limit the y-axis in zoom-in
zoom_in = TRUE;  # default zoom in to zoom_in_location
zoom_in_group = NULL
blevs<-gcvb[order(-mflare,-global_shap),unique(as.character(variable))]
shapobs_long[,base:=factor(base, levels=blevs)]

p <- ggplot(shapobs_long, aes_string(x = id, y = value , fill = variable)) +
    geom_col(width =1, alpha = 0.9) +
    # geom_area() +
    labs(fill = 'Feature', x = 'Observation',
         y = 'SHAP values by feature:\n (Contribution to the base value)') +
    geom_hline(yintercept = 0, col = "gray40") +
    theme_bw() + 
#    theme(legend.position="none") + 
    scale_fill_manual(values=clustcols,name="SHARP parameter")
    





#gcv_sub<-gcv_clust[rnk < 10]
#incl<-gcv_sub[order(-mflare, rnk),unique(as.character(variable))]
#gcv_clust[,variable:=factor(variable, levels=incl)]
#gcv[,hclust10
#ggplot(gcv_clust[variable %in% incl][mflare > 0.1], 
#       aes(x=variable, y=global_shap,fill=factor(hclust10))) + 
#  geom_bar(stat="identity",position="dodge") #+ 
#  facet_wrap(~hclust10,scales="free_x") + 
#  scale_fill_manual(breaks=)
#
#
#
#
#chk<-gcv_clust[rnk < 20]
#allv<-chk[,as.character(unique(variable))]
#for(k





for(i in 1:combos[,.N]){
#for(i in (i+1):combos[,.N]){
  print(
  me.dep.plot(plotd,combos[i,xvar],combos[i,yvar],quantx=T) + 
	  facet_wrap(~fold,labeller="label_both") 
  )
  #)
  #readline(prompt="Press [enter] to continue")
  #print(
#  b<-me.dep.plot(plotd,combos[i,yvar],combos[i,xvar],quantx=T) + 
#	  facet_wrap(~fold,labeller="label_both") 
  readline(prompt="Press [enter] to continue")
  
#  print(g<-grid.arrange(a,b) )
  
}
## 58
ggsave("interaction_epsx_totusjh_ind45.png")
ggsave("interaction_totpot_epsx_ind82.png")
ggsave("interaction_xr_max_totusjh_ind100.png")

















plotd<-fread(file="sub_plot_shap_long_top20.tsv")
setnames(plotd,old=c("shap","global_shap"), new=c("value","mean_value"))
levs<-plotd[order(rnk),head(.SD,1),by=variable][,as.character(variable)]
plotd[,variable:=factor(variable, levels=levs)]
plotd[,rfvalue:=as.numeric(rfvalue)]
#shap.plot.summary(data_long=plotd)
#cls<-brewer.pal(8, "Set1")[c(2,2,5,4,4)]
#cls<-brewer.pal(8, "Dark2")[c(1,3,3, 3, 2)]
cls<-c(brewer.pal(8, "Paired")[3], brewer.pal(8, "Dark2")[c(3,3,3, 2)])
top20<-plotd[,.(variable, mean_value)][order(-mean_value)][,unique(.SD)][1:20][,as.character(variable)]
me.shap.plot.summary(dat=plotd[variable %in% top20],cols=cls) + 
    annotate("text",
                    x = 20 + 0.25, y=-1.3,label="Global Importance",
                    size = 3.1, alpha = 0.7,
                    hjust = -0.2,
                    fontface = "bold")  + 
            labs(title="SHAP: Top 20 Final Variables") + 
            theme(plot.title=element_text(hjust=0.5,size=12)) + 
	    coord_cartesian(clip="off") + coord_flip()
ggsave("global_shap.png",width=4,height=4,dpi=300)





xvars<-plotd[,as.character(unique(variable))]
xvars<-intersect(xvars, "epsx_min_t10_5")
for(x in xvars){
  for(i in plotd[variable!=x,as.character(unique(variable))]){
    for(j in plotd[variable!=x][variable!=i,as.character(unique(variable))]){
          print(
          me.dep2.plot(plotd, 
          	     x=x, 
          	     color_feature=j, 
          	     shape_feature=i,wrapfold=T)
          )
            readline(prompt="Press [enter] to continue")
    }
  }
}

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









