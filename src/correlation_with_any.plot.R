## #!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

message(sprintf("cmd: Rscript correlation_with_any.plot.R longmethyl.forplot.txt [plot_title query_title target_title]"))

library(ggplot2)
## geom_bin2d() may be a better way ##
# https://stackoverflow.com/questions/50331320/how-to-move-tick-marks-and-labels-at-right-left-end-of-tiles-in-geom-tile-ggplot
# https://stackoverflow.com/questions/10710463/how-to-force-the-x-axis-tick-marks-to-appear-at-the-end-of-bar-in-heatmap-graph
library(gridExtra)
library(grid)
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(scales)
library(ggthemes)
library(extrafont)
library(plyr)
library(stringr)

corrdata = args[1]

str2bool = function(input_str)
{
    if(input_str == "yes"){
        input_str = TRUE
    }else{
        input_str = FALSE
    }
    return(input_str)
}

if(length(args)>=2){
    plotlabel = args[2]
}else{
    plotlabel = "CG-HG002"
}

if(length(args)>=3){
    query_label = args[3]  # query label
}else{
    query_label = "Query methylation frequency"
}

if(length(args)>=4){
    target_label = args[4]  # target label
}else{
    target_label = "Target methylation frequency"
}

# font_import()
# fonts()


classify_one_type <- function(rmet, sranges){
  for(i in 2:(length(sranges))){
    if(rmet<=sranges[i]){return(i-1)}
  }
}

classify_rtype <- function(Rmet, sranges){
  return(sapply(Rmet, classify_one_type, sranges))
}

generate_rtype_matrix <- function(Rmet1, Rmet2, rescale_num=10){
  scale_range <- 1 / rescale_num
  sranges <- seq(0, 1, scale_range)
  
  res_mat <- matrix(data=0, nrow=length(sranges)-1, ncol=length(sranges)-1, 
                    dimnames = list(sranges[2:length(sranges)], 
                                    sranges[2:length(sranges)]))
  # res_mat <- matrix(data=0, nrow=length(sranges)+1, ncol=length(sranges)+1, 
  #                   dimnames = list(c(sranges, '1+'), c(sranges, '1+')))
  rtype1 = classify_rtype(Rmet1, sranges)
  rtype2 = classify_rtype(Rmet2, sranges)
  for(i in 1:length(rtype1)){
    tmp = res_mat[rtype1[i], rtype2[i]]
    res_mat[rtype1[i], rtype2[i]] = tmp + 1
  }
  res_mat
}

plot_rmet_heatmap <- function(rmet_df_file, hunit=40, 
                              data_label="", 
                              name_query="",
                              name_target="", 
                              pearcor=0.0, 
                              setted_maxvalue = 0, 
                              titlehjust=-0.3,
                              istitley=TRUE,
                              istitlex=TRUE,
                              istitle=TRUE){
  rmet_df <- read.table(rmet_df_file, 
                        header = T, sep = '\t', stringsAsFactors = F)
  rmet_target <- rmet_df$rmet_target
  rmet_query <- rmet_df$rmet_query
  if(setted_maxvalue==0){
    # hunit=40
    rmet_mat_10 <- as.data.frame(generate_rtype_matrix(rmet_query, rmet_target, hunit))
    rmet_mat_10[rmet_mat_10==0] = 1
    sim.df <- cbind(ID=as.numeric(rownames(rmet_mat_10)), rmet_mat_10)
    sim.df <- melt(sim.df, id.vars = 'ID')
    sim.df$variable <- as.numeric(as.character(sim.df$variable))
    sim.df$logvalue <- log10(sim.df$value)
    sim.df$label <- as.character(round(sim.df$value, 2))
    
    maxvalue = max(sim.df$value)
  }else{
    maxvalue = setted_maxvalue
  }
  message(sprintf("maxvalue: %d", maxvalue))
  # if(maxvalue < 1000000){
  #   maxvalue = 1000000
  # }
  
  # RdYlGn
  # RdYlBu
  rf <- colorRampPalette(rev(brewer.pal(11,'RdYlBu')))
  r <- rf(32)
  
  if(pearcor<=0){
    pearcor = cor(rmet_target, rmet_query, method='pearson')
  }
  pearcor_fmt <- format(round(pearcor, 4), nsmall = 4)
  
  p <- ggplot(rmet_df, aes(rmet_target, rmet_query)) + 
    geom_bin2d(bins=hunit) +
    theme_bw() + 
    scale_fill_gradientn(colors = r, 
                         limit = c(1, maxvalue), 
                         space = "Lab", name="count", trans="log10", 
                         labels=comma) + 
    scale_x_continuous(breaks = seq(0, 1, 0.2)) +
    scale_y_continuous(breaks = seq(0, 1, 0.2)) +
    theme(text=element_text(size=9, family = "Arial"), 
          axis.title.x = element_text(size = 9),
          axis.title.y = element_text(size = 9), 
          axis.text=element_text(size=8), 
          legend.text = element_text(size=6), 
          legend.title = element_text(size=8), 
          legend.key.height = unit(0.5, "cm"), 
          legend.key.width = unit(0.3, "cm"),
          plot.title = element_text(hjust = titlehjust)) +
    ylab(name_query) + xlab(name_target)
  
  if(istitle){
    p <- p + ggtitle(bquote(paste(.(data_label), ",  ", italic("r"), " = ", .(pearcor_fmt))))
  }
  if(!istitlex){
    p <- p + theme(axis.title.x = element_blank())
  }
  if(!istitley){
    p <- p + theme(axis.title.y = element_blank())
  }
  reslist <- list("plot"=p, "pearson"=pearcor_fmt)
  return(reslist)
}


# plot =======
setted_hunit = 40

r_qt <- plot_rmet_heatmap(corrdata, 
                         hunit = setted_hunit,
                         data_label = plotlabel, 
                         name_target = str_replace_all(target_label, "_", " "),
                         name_query = str_replace_all(query_label, "_", " "),
                         titlehjust = -0.0, 
                         pearcor = 0.0, 
                         istitlex = TRUE, 
                         istitley = TRUE)
p_qt <- r_qt$plot


ppi= 300
png(paste(corrdata, "_plot.png", sep=""),
    width = 12,
    height = 9, units = "cm", res=ppi)
p_qt
dev.off()







