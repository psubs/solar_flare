#data_long<-copy(plotd)
#x_bound = NULL;
#dilute = FALSE; scientific = FALSE;
#my_format = NULL
#cls<-brewer.pal(3,"Dark2")[2:3]
#cls<-cls2<-c("#FFCC33", "#6600CC")
#display.brewer.all(colorblindFriendly=T)
#cls<-brewer.pal(12,"Paired")[10:11]


#dat=plotd;
#shap_name="value";
#variable_name="variable";
#global_shap_name="mean_value";
#cols=NULL;
#x_bound = NULL;
#qval=TRUE;
#dilute = FALSE; 
#scientific = FALSE;
#my_format = NULL
#
#
#dat=plotdb;
#cols=cls;
#variable_name="base";
#shap_name="base_shap_id"; 
#global_shap_name="global_base_shap";
#show_color=F


me.shap.plot.summary <- function(dat,
				 variable_name="variable",
                                 shap_name="value",
				 global_shap_name="mean_value",
				 show_color=TRUE,
                                 cols=NULL,
                                 x_bound = NULL,
			         qval=TRUE,
                                 dilute = FALSE, 
				 scientific = FALSE,
                                 my_format = NULL){
  data_long<-copy(dat)
  setnames(data_long, old=shap_name, new="shap_name_tmp") 
  setnames(data_long, old=global_shap_name, new="global_name_tmp") 
  if(variable_name!="variable") setnames(data_long, old=variable_name, new="variable")

  if(missing(cols)) cols<-c("#FFCC33", "#6600CC")
  if (scientific){label_format = "%.1e"} else {label_format = "%.3f"}
  if (!is.null(my_format)) label_format <- my_format
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

  x_bound <- if (is.null(x_bound)) max(abs(data_long$shap_name_tmp))*1.1 else as.numeric(abs(x_bound))

  if(show_color){
  
    plot1 <- ggplot(data = data_long)+
      coord_flip(ylim = c(-x_bound, x_bound)) +
      # sina plot:
      ggforce::geom_sina(aes(x = variable, 
			     y = shap_name_tmp, 
			     color = if(qval) qvalue else stdfvalue),
                method = "counts", 
		maxwidth = 0.7, 
		size=0.01,
		alpha = 0.7) +
      # print the mean absolute value:
      geom_text(data = unique(data_long[, c("variable", "global_name_tmp")]),
                aes(x = variable, y=-Inf, label = sprintf(label_format, global_name_tmp)),
                size = 2.1, alpha = 0.7,
                hjust = -0.2,
                fontface = "bold") + # bold
      # # add a "SHAP" bar notation
      # annotate("text", x = -Inf, y = -Inf, vjust = -0.2, hjust = 0, size = 3,
      #          label = expression(group("|", bar(SHAP), "|"))) +
      #scale_color_gradient(low=cls[1], high=cls[2],
      #                     breaks=c(0,1), labels=c(" Low","High "),
      #                     guide = guide_colorbar(barwidth = 12, barheight = 0.3)) +
      #scale_color_distiller(palette="Spectral",
      #                     guide = guide_colorbar(barwidth = 12, barheight = 0.3)) +
      scale_color_gradientn(colours=cols,
                           breaks=c(0,.2,.8,1), #labels=c("Low","High "),
                           labels=c("0", "0.2", "0.8", "1"),
                           guide = guide_colorbar(barwidth = 4, 
						  barheight = 0.3,
						  label.vjust=1)) + # barwidth=12
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
                       labels =SHAPforxgboost:::label.feature(rev(levels(data_long$variable))))+
      labs(y = "SHAP value (impact on log-odds of MSF)", 
	   x = "", 
	   color = if(qval) "Quantile Value\n of Feature" else "Feature value  ")
    
  
  } else {
  
    plot1 <- ggplot(data = data_long)+
      coord_flip(ylim = c(-x_bound, x_bound)) +
      # sina plot:
      ggforce::geom_sina(aes(x = variable, 
			     y = shap_name_tmp), 
			     method = "counts", 
			     maxwidth = 0.7, 
			     alpha = 0.2,
			     size=0.01,
			     colour="#D95F02") +
      # print the mean absolute value:
      geom_text(data = unique(data_long[, c("variable", "global_name_tmp")]),
                aes(x = variable, y=-Inf, label = sprintf(label_format, global_name_tmp)),
                size = 2.1, alpha = 0.7,
                hjust = -0.2,
                fontface = "bold") + # bold
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
                       labels =SHAPforxgboost:::label.feature(rev(levels(data_long$variable))))+
      labs(y = "SHAP value (impact on log-odds of MSF)", x = "")
    
  
  
  
  
  
  }
  #setnames(data_long, new=shap_name, old="shap_name_tmp") 
  return(plot1)
}

#data_long=copy(plotd); 
#
#
#x="epsx_min_t10_5"; 
#y = NULL;
#y<-x
#color_feature = "shrgt45_min_t10_5"; 
#data_int = NULL;
#
#x<-"r_value_max_t10_5"
#color_feature<-"epsx_min_t10_5"
#me.dep.plot(plotd,x,color_feature,quantx=F)
#dat<-copy(plotd)
#x<-combos[i,as.character(xvar)]
#color_feature<-combos[i,as.character(yvar)]
me.dep.plot<-function(dat, x, color_feature,quantx=T){
  data_long<-copy(dat)
  x<-as.character(x)
  color_feature<-as.character(color_feature)
  y<-x
  cast_vars=intersect(names(data_long),c("qvalue", "value", "rfvalue"))
  dilute = FALSE; smooth = TRUE; size0 = NULL; add_hist = FALSE
  scientific=T
  label_format<-function(x) if((x > 1e-3 && x < 10^3) || x==0){
	  sprintf("%.2f", x)
  } else {
	  sprintf("%.1e", x)
  }
  
  
  dw<-data_long[variable %in% c(y, x, color_feature)]
  if("hclust10" %in% names(data_long)){
    dw<-dcast(dw, fid + fold + hclust10 ~ variable, value.var=cast_vars)
  } else{
    dw<-dcast(dw, fid + fold ~ variable, value.var=cast_vars)
  }
  pvars<-unique(unlist(sapply(cast_vars, 
  			function(i) paste0(i, "_", c(x, y, color_feature)))))
  setnames(dw, old=paste0("rfvalue_", color_feature), new=paste0("rfvalue_cf"))
  setnames(dw, old=paste0("rfvalue_", x), new="xtmp")
  
  labsh<-dw[,quantile(rfvalue_cf,c(0, 0.2, 0.5, 0.8, 1))]
  labsh<-sapply(labsh, label_format)#label_format(labs)
  labsh<-paste0(labsh, ": ","(", names(labsh),")" )
  
  if(quantx){
    labsv<-dw[,quantile(xtmp,c(0, 0.2, 0.5, 0.8, 1))]
    labsv<-sapply(labsv, label_format)#label_format(labs)
    labsv<-paste0(labsv, "\n", "(", names(labsv), ")")
  } 
 # else{
 #   labsv<-dw[,quantile(xtmp,c(0, 0.2, 0.5, 0.8, 1))]
 #   labsv<-sapply(labsv, label_format)#label_format(labs)
 #   labsv<-paste0(labsv, "\n", "(", names(labs), ")")
 # } 
  
  
  setnames(dw, new=paste0("rfvalue_", color_feature), old=paste0("rfvalue_cf"))
  setnames(dw, new=paste0("rfvalue_", x), old="xtmp")
  p<-ggplot(dw, aes_string(x=paste0("qvalue_",x), 
  		           y=paste0("value_", y),
  		           color=paste0("qvalue_", color_feature))) + 
             geom_point(size=0.6,alpha=0.5) + 
  	     scale_color_gradientn(colours=cols,
                                   breaks=c(0,.2,0.5, .8,1), #labels=c("Low","High "),
                                   labels=labsh, limits=c(0,1),
  				   name=paste0(color_feature, "\nRaw value (quantile)")) 
  if(quantx){
   p<-p + scale_x_continuous(breaks=c(0, .2, .5, .8, 1),
  			      labels=labsv,
			      limits=c(0,1)) + 
           theme_bw() + 
  	   labs(y = paste0("SHAP value: ", y), 
  		x = paste0(x,"\nRaw value (quantile)"),
		title=paste0("Dependence plot: \n", x, "vs. \n", color_feature)) +
           theme(plot.title=element_text(hjust=0.5))
  } else {
   p<-p + theme_bw() + 
  	  labs(y = paste0("SHAP value: ", y), 
  	       x = paste0(x,"\nRaw value"),
	       title=paste0("Dependence plot: \n", x, "vs. \n", color_feature)) +
          theme(plot.title=element_text(hjust=0.5))
  }
  return(p)
}
#me.dep2.plot(plotd, 
#	     x="totpot_max_t10_5", 
#	     color_feature="epsx_min_t10_5", 
#	     shape_feature="r_value_min_t10_5")
#
#x="totpot_max_t10_5"; 
#color_feature="epsx_min_t10_5"; 
#shape_feature="r_value_max_t10_5"
#shape_feature="absnjzh_sd_t10_5"
#me.dep2.plot(plotd, x, color_feature, shape_feature)

me.dep2.plot<-function(dat, x, color_feature,shape_feature,quantx=T,wrapfold=F){

  data_long<-copy(dat)
  x<-as.character(x)
  color_feature<-as.character(color_feature)
  shape_feature<-as.character(shape_feature)
  y<-x
  cast_vars=c("qvalue", "value", "rfvalue")
  dilute = FALSE; smooth = TRUE; size0 = NULL; add_hist = FALSE
  scientific=T
  label_format<-function(x) if((x > 1e-3 && x < 10^3) || x==0){
	  sprintf("%.2f", x)
  } else {
	  sprintf("%.1e", x)
  }
  
  
  dw<-data_long[variable %in% c(y, x, color_feature,shape_feature)]
  dw<-dcast(dw, fid + fold ~ variable, value.var=cast_vars)
  pvars<-unique(unlist(sapply(cast_vars, 
  			function(i) paste0(i, "_", c(x, y, color_feature)))))
  setnames(dw, old=paste0("rfvalue_", color_feature), new=paste0("rfvalue_cf"))
  setnames(dw, old=paste0("rfvalue_", x), new="xtmp")
  setnames(dw, old=paste0("qvalue_", shape_feature), new="qvalue_sf")

  dw[,sf:=cut(qvalue_sf,c(0,.2,.8,1),include.lowest=T,labels=c("low", "med", "high"))]
  dw[,sf:=factor(sf, levels=c("low", "med", "high"))]
  labdt<-dw[,lapply(.SD, function(x) paste0(round(min(x),2),"< ", 
				     shape_feature, " <", 
				     round(max(x),2))),
     by=sf,
     .SDcols=c(paste0("rfvalue_",shape_feature))]#[,paste0("rfvalue_",shape_feature),with=F]
  labssh<-labdt[order(sf)][[2]]
  breakssh<-labdt[order(sf)][[1]]

  labsh<-dw[,quantile(rfvalue_cf,c(0, 0.2, 0.5, 0.8, 1))]
  labsh<-sapply(labsh, label_format)#label_format(labs)
  labsh<-paste0(labsh, ": ","(", names(labsh),")" )
  
  if(quantx){
    labsv<-dw[,quantile(xtmp,c(0, 0.2, 0.5, 0.8, 1))]
    labsv<-sapply(labsv, label_format)#label_format(labs)
    labsv<-paste0(labsv, "\n", "(", names(labs), ")")
  } 
 # else{
 #   labsv<-dw[,quantile(xtmp,c(0, 0.2, 0.5, 0.8, 1))]
 #   labsv<-sapply(labsv, label_format)#label_format(labs)
 #   labsv<-paste0(labsv, "\n", "(", names(labs), ")")
 # } 
  
  
  setnames(dw, new=paste0("rfvalue_", color_feature), old=paste0("rfvalue_cf"))
  setnames(dw, new=paste0("rfvalue_", x), old="xtmp")
  setnames(dw, new=paste0("qvalue_", shape_feature), old=paste0("qvalue_sf"))
  
  p<-ggplot(dw, aes_string(x=paste0("qvalue_",x), 
  		           y=paste0("value_", y),
  		           color=paste0("qvalue_", color_feature),
			   shape="sf",size="sf")) + 
             geom_point(alpha=0.9) + 
  	     scale_color_gradientn(colours=cols,
                                   breaks=c(0,.2,0.5, .8,1), #labels=c("Low","High "),
                                   labels=labsh, 
				   limits=c(0,1),
  				   name=paste0(color_feature, "\nRaw value (quantile)")) +

  	     scale_shape_manual(values=c(2, 4, 1),
                                breaks=breakssh, 
                                labels=labssh,
  				name=paste0(shape_feature)) +

  	     scale_size_manual(values=c(1, 5, 10),
                                breaks=breakssh, 
                                labels=labssh,
  				name=paste0(shape_feature)) 
  if(quantx){
   p<-p + scale_x_continuous(breaks=c(0, .2, .5, .8, 1),
  			      labels=labsv,
			      limits=c(0,1)) + 
           theme_bw() + 
  	   labs(y = paste0("SHAP value: ", y), 
  		x = paste0(x,"\nRaw value (quantile)"),
		title=paste0("Dependence plot: \n", x, "vs. \n", color_feature)) +
           theme(plot.title=element_text(hjust=0.5))
  } else {
   p<-p + theme_bw() + 
  	  labs(y = paste0("SHAP value: ", y), 
  	       x = paste0(x,"\nRaw value"),
	       title=paste0("Dependence plot: \n", x, "vs. \n", color_feature)) +
          theme(plot.title=element_text(hjust=0.5))
  }
  if(wrapfold){
    p<-p + facet_wrap(~fold) 
  }
  return(p)
}

                                    #guide = guide_colorbar(barwidth = 12, barheight = 0.3)) 
#shap.plot.dependence(plotd, x)
me.shap.plot.dependence<-
function (data_long, x, y = NULL, color_feature = NULL, data_int = NULL, 
    dilute = FALSE, smooth = TRUE, size0 = NULL, add_hist = FALSE) 
{
    if (is.null(y)) 
        y <- x
    data0 <- data_long[variable == y, .(variable, shap_name_tmp)]
    data0$x_feature <- data_long[variable == x, rfvalue]
    if (!is.null(color_feature)) 
        data0$color_value <- data_long[variable == color_feature, 
            rfvalue]
    if (!is.null(data_int)) 
        data0$int_value <- data_int[, x, y]
    nrow_X <- nrow(data0)
    if (is.null(dilute)) 
        dilute = FALSE
    if (dilute != 0) {
        dilute <- ceiling(min(nrow(data0)/10, abs(as.numeric(dilute))))
        set.seed(1234)
        data0 <- data0[sample(nrow(data0), min(nrow(data0)/dilute, 
            nrow(data0)/2))]
    }
    if (x == "dayint") {
        data0[, `:=`(x_feature, as.Date(data0[, x_feature], format = "%Y-%m-%d", 
            origin = "1970-01-01"))]
    }
    if (is.null(size0)) 
        size0 <- if (nrow(data0) < 1000L) 
            1
        else 0.4
    plot1 <- ggplot(data = data0, aes(x = x_feature, y = if (is.null(data_int)) 
        shap_name_tmp # chk
    else int_value, color = if (!is.null(color_feature)) 
        color_value
    else NULL)) + geom_point(size = size0, alpha = if (nrow(data0) < 
        1000L) 
        1
    else 0.6) + labs(y = if (is.null(data_int)) 
        paste0("SHAP value for ", label.feature(y))
    else paste0("SHAP interaction values for\n", label.feature(x), 
        " and ", label.feature(y)), x = label.feature(x), color = if (!is.null(color_feature)) 
        paste0(label.feature(color_feature), "\n", "(Feature value)")
    else NULL) + scale_color_gradient(low = "#FFCC33", high = "#6600CC", 
        guide = guide_colorbar(barwidth = 10, barheight = 0.3)) + 
        theme_bw() + theme(legend.position = "bottom", legend.title = element_text(size = 10), 
        legend.text = element_text(size = 8))
    if (smooth) {
        plot1 <- plot1 + geom_smooth(method = "loess", color = "red", 
            size = 0.4, se = F)
    }
    plot1 <- plot.label(plot1, show_feature = x)
    if (add_hist) {
        plot1 <- ggExtra::ggMarginal(plot1, type = "histogram", 
            bins = 50, size = 10, color = "white")
    }
    plot1
}























#shap.plot.dependence(data_long=plotd, 
#		     x="r_value_max_l1_t10_5",
#		    # y="shrgt45_min_t10_5",
#		     color_feature="shrgt45_min_t10_5",add_hist=T)#,size0=2) #+ 
#coord_cartesian(ylim=c(-1,1)) 
##dcast(plotd[variable %in% c("r_value_max_l1_t10_5", "shrgt45_min_t10_5")]
## ggplot(plotd, aes(x=,y=)) + geom_point()
#x="r_value_max_l1_t10_5";
#y ="shrgt45_min_t10_5";
#color_feature = "shrgt45_min_t10_5";
#data_int = NULL;  # if supply; will plot SHAP interaction values
#dilute = FALSE;
#smooth = TRUE;
#size0 = NULL;
#add_hist = FALSE
#data_long<-copy(plotd)
#plotd[,rfvalue:=as.numeric(rfvalue)]
