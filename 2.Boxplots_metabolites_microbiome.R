# Read the SE object containing all the metabolites from Liver/Serum/Fecal matrices and the microbiome data.
# Generate pdf file containing the Log ion counts for all metabolites from different matrices and another pdf file containing Log abundance of microbiomes.

library(ggpubr)
library(rstatix)
library(gridExtra)
library(tidyverse)
library(autonomics)
library(maplet)
library(magrittr)
library(grid)
library(gridExtra)
library(data.table)
library(broom)


Ds <- D_rmerged_ALL


rowData(D_rmerged_ALL)


rowData(Ds)$tissue= sub(".*\\.", "", rownames(Ds))
Ds %<>% mt_modify_filter_features(filter=tissue!="Microbiome") 

## get stat table for CTR vs PRN
stat_id_ctr_PRN <- metadata(Ds)$results %>%
  purrr::map_lgl(~"stats" %in% .x$fun && .x$output$name == "Gctr_prot") %>%
  which()
stat_tab_ctr_PRN <- metadata(Ds)$results[[stat_id_ctr_PRN]]$output$table  %>% select(var,term,statistic,p.value)
stat_tab_ctr_PRN$group1 <- "CTR"
stat_tab_ctr_PRN$group2 <- "PRN"


## get stat table for CTR vs NEAA
stat_id_ctr_NEAA <- metadata(Ds)$results %>%
  purrr::map_lgl(~"stats" %in% .x$fun && .x$output$name == "Gctr_NEAA") %>%
  which()
stat_tab_ctr_NEAA <- metadata(Ds)$results[[stat_id_ctr_NEAA]]$output$table %>% select(var,term,statistic,p.value)
stat_tab_ctr_NEAA$group1 <- "CTR"
stat_tab_ctr_NEAA$group2 <- "NEAA"

## get stat table for PRN vs NEAA
stat_id_PRN_NEAA <- metadata(Ds)$results %>%
  purrr::map_lgl(~"stats" %in% .x$fun && .x$output$name == "prot_NEAA") %>%
  which()
stat_tab_PRN_NEAA <- metadata(Ds)$results[[stat_id_PRN_NEAA]]$output$table %>% select(var,term,statistic,p.value)
stat_tab_PRN_NEAA$group1 <- "PRN"
stat_tab_PRN_NEAA$group2 <- "NEAA"

stat_tab <- rbind(stat_tab_ctr_PRN, stat_tab_ctr_NEAA) %>% rbind(stat_tab_PRN_NEAA)

stat_tab %<>% mutate(p.sig.if = case_when(
                                p.value <=0.001 ~ "***",
                                p.value > 0.001 & p.value <=0.01 ~ "**",
                                p.value>0.01 & p.value <=0.05 ~ "*"
                               
                                )
                    )



stat_tab[stat_tab$var == "Xanthine.Serum",]





#Ds <- Ds[-which(Ds$tissue == "Microbiome"),]


df <- data.frame(x=c(0,1),y=c(1,1))
empty_plot <- ggplot(df) + geom_point(aes(x,y),color='black',  shape=21) + xlim(0, 1) + ylim(0, 1) +  xlab(" ")+ylab(" ")+ggtitle(" ") +
  theme(panel.background = element_blank(), panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=.5),aspect.ratio=.8)

try(dev.off(),silent = TRUE)
pdf("plot_metabolites_color.pdf",height=12,width=12)
par(mfrow=c(4,3))

rowData(Ds)$name <- make.names(rowData(Ds)$name)
ix=sort(rowData(Ds)$name,index.return=TRUE)$ix

colour=rep(NA,length(rowData(Ds)$tissue))
colour[which(rowData(Ds)$tissue=="Serum")]="#F4B183"
colour[which(rowData(Ds)$tissue=="Fecal")]="#525252"
colour[which(rowData(Ds)$tissue=="Liver")]="#6868AE"
# for (i in seq(dim(D_rmerged2)[1])){
plot_list <- list()
my_comparisons = list(c("CTR", "NEAA"), c("NEAA", "PRN"), c("CTR", "PRN"))
k=-1
j=1
for (i in ix){
  k=k+1
  l=k%%3
  lmet=rowData(Ds)$name[i]
  
 # if(lmet = "Citric.acid.Liver"){
  cat("working on", i,lmet,"\n")
  pdata=assay(Ds)[i,]
  
  
  
  if(grepl("Fecal",lmet)&(l!=0)){
   
    p <- empty_plot
    plot_list[[j]] <- p
    
    j <- j+1
    k=k+1
    l=k%%3
  } 
  if(grepl("Fecal",lmet)&(l!=0)){
    
    p <- empty_plot
    plot_list[[j]] <- p
    
    j <- j+1
    k=k+1
    l=k%%3
  } 
  if(grepl("Liver",lmet)&(l!=1)){
    
    p <- empty_plot
    plot_list[[j]] <- p
    
    j <- j+1
    k=k+1
    l=k%%3
  } 
  if(grepl("Liver",lmet)&(l!=1)){
    
    p <- empty_plot
    plot_list[[j]] <- p
    
    j <- j+1
    k=k+1
    l=k%%3
  } 
  if(grepl("Serum",lmet)&(l!=2)){
    
    p <- empty_plot
    plot_list[[j]] <- p
    
    j <- j+1
    k=k+1
    l=k%%3
  } 
  if(grepl("Serum",lmet)&(l!=2)){
    
    p <- empty_plot
    plot_list[[j]] <- p
    
    j <- j+1
    k=k+1
    l=k%%3
  } 
 # xdata= Ds$Group
 # ggdata <- cbind(pdata,xdata)
 #  boxplot(pdata~xdata,notch=FALSE,xlab="Group",ylab="Log.Concentration",main=lmet, col=colour[i])
  
 # xdata=data.frame(Group =  Ds$Group)
  xdata=data.frame(Group =  case_when(
    Ds$Group == "Control" ~"CTR",   
    Ds$Group == "Reduced NEAA" ~ "NEAA",  
    Ds$Group == "Reduced Protein" ~ "PRN"))
  
 
  ggdata <- cbind(pdata,xdata) %>% dplyr::mutate(met = .[,1])
  ymax <- max(ggdata$met, na.rm = TRUE)
  
  ggdata$Group <- factor(ggdata$Group , levels=c("CTR", "PRN", "NEAA"))
  
  # stat.test <- 
  #   lm(met ~ Group, ggdata ) %>% summary(.) %>% tidy()
  
  # stat.test1<- ggdata %>% 
  #   t_test(met ~ Group,p.adjust.method = "none" ) 
  
  bon_sig <- .05/437
  
  stat.test2 <- stat_tab[stat_tab$var == lmet,] 
  
   stat.test <- ggdata %>% 
     t_test(met ~ Group,p.adjust.method = "bonferroni" ) %>%
    add_xy_position()
   
   
   
   
  
   
  # 
  
  p <- ggboxplot(ggdata, x = "Group", y = "met",fill = colour[i], palette = "jco", outlier.size=1.1, outlier.shape = 21 ) +
    xlab("")+ ylab ("Log ion count") + ggtitle(lmet) + theme(aspect.ratio=.8) +
    stat_pvalue_manual(stat.test2, hide.ns = TRUE, label = "p.sig.if",y.position = max(ggdata$met)+0.12, tip.length = 0.01 , size = 8.5,
                       bracket.shorten = 0.05, step.increase = .13 ) +stat_boxplot(geom = "errorbar", width = .3) + 
    border() + scale_y_continuous(expand = c(0.1,0))+theme(text = element_text(size = 11), plot.title = element_text(face="bold"))+
    font("xy.text", size = 10)
 
  
  plot_list[[j]] <- p
  
  j <- j+1
  
  #}
   
}


ggsave(
  filename = "Boxplots_all_metabolites.pdf",
  plot = marrangeGrob(plot_list, nrow=4, ncol=3,layout_matrix = matrix(1:9, 3, 3, TRUE)),
  width = 10, height = 12
)

dev.off()



######################################################################################################

#Microbiom

######################################################################################

Db <- D_rmerged_ALL

rowData(Db)$tissue= sub(".*\\.", "", rownames(Db))
Db %<>% mt_modify_filter_features(filter=tissue=="Microbiome") 

## get stat table for CTR vs PRN
stat_id_ctr_PRN <- metadata(Db)$results %>%
  purrr::map_lgl(~"stats" %in% .x$fun && .x$output$name == "Gctr_prot") %>%
  which()
stat_tab_ctr_PRN <- metadata(Db)$results[[stat_id_ctr_PRN]]$output$table  %>% select(var,term,statistic,p.value)
stat_tab_ctr_PRN$group1 <- "CTR"
stat_tab_ctr_PRN$group2 <- "PRN"


## get stat table for CTR vs NEAA
stat_id_ctr_NEAA <- metadata(Db)$results %>%
  purrr::map_lgl(~"stats" %in% .x$fun && .x$output$name == "Gctr_NEAA") %>%
  which()
stat_tab_ctr_NEAA <- metadata(Db)$results[[stat_id_ctr_NEAA]]$output$table %>% select(var,term,statistic,p.value)
stat_tab_ctr_NEAA$group1 <- "CTR"
stat_tab_ctr_NEAA$group2 <- "NEAA"

## get stat table for PRN vs NEAA
stat_id_PRN_NEAA <- metadata(Db)$results %>%
  purrr::map_lgl(~"stats" %in% .x$fun && .x$output$name == "prot_NEAA") %>%
  which()
stat_tab_PRN_NEAA <- metadata(Db)$results[[stat_id_PRN_NEAA]]$output$table %>% select(var,term,statistic,p.value)
stat_tab_PRN_NEAA$group1 <- "PRN"
stat_tab_PRN_NEAA$group2 <- "NEAA"

stat_tab <- rbind(stat_tab_ctr_PRN, stat_tab_ctr_NEAA) %>% rbind(stat_tab_PRN_NEAA)

stat_tab %<>% mutate(p.sig.if = case_when(
  p.value <=0.001 ~ "***",
  p.value > 0.001 & p.value <=0.01 ~ "**",
  p.value>0.01 & p.value <=0.05 ~ "*"
  
)
)
stat_tab$var = sub(".Microbiome", "", stat_tab$var)


df <- data.frame(x=c(0,1),y=c(1,1))
empty_plot <- ggplot(df) + geom_point(aes(x,y),color='black',  shape=21) + xlim(0, 1) + ylim(0, 1) +  xlab(" ")+ylab(" ")+ggtitle(" ") +
  theme(panel.background = element_blank(), panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=.5),aspect.ratio=.8)

try(dev.off(),silent = TRUE)

rowData(Db)$name <- make.names(rowData(Db)$name)
ix=sort(rowData(Db)$name,index.return=TRUE)$ix

colour="yellow"

plot_list <- list()
my_comparisons = list(c("CTR", "NEAA"), c("NEAA", "PRN"), c("CTR", "PRN"))

j=1
for (i in ix){
 
  lmet=rowData(Db)$name[i]
  
  # if(lmet = "Citric.acid.Liver"){
  cat("working on", i,lmet,"\n")
  pdata=assay(Db)[i,]
  
  
  xdata=data.frame(Group =  case_when(
    Ds$Group == "Control" ~"CTR",   
    Ds$Group == "Reduced NEAA" ~ "NEAA",  
    Ds$Group == "Reduced Protein" ~ "PRN"))
  
  
  ggdata <- cbind(pdata,xdata) %>% dplyr::mutate(met = .[,1])
  ymax <- max(ggdata$met, na.rm = TRUE)
  
  ggdata$Group <- factor(ggdata$Group , levels=c("CTR", "PRN", "NEAA"))
  
   bon_sig <- .05/437
  
  stat.test <- stat_tab[stat_tab$var == lmet,] 
  
  p <- ggboxplot(ggdata, x = "Group", y = "met",fill = colour, palette = "jco", outlier.size=1.1, outlier.shape = 21 ) +
    xlab("")+ ylab ("Log abundance") + ggtitle(lmet) + theme(aspect.ratio=.8) +
    stat_pvalue_manual(stat.test, hide.ns = TRUE, label = "p.sig.if",y.position = max(ggdata$met)+0.12, tip.length = 0.01 , size = 8.5,
                       bracket.shorten = 0.05, step.increase = .13 ) +stat_boxplot(geom = "errorbar", width = .3) + 
    border() + scale_y_continuous(expand = c(0.1,0))+theme(text = element_text(size = 11), plot.title = element_text(face="bold"))+
    font("xy.text", size = 10)
  
  
  plot_list[[j]] <- p
  
  j <- j+1
  
  #}
  
}


ggsave(
  filename = "Boxplots_all_microbiome.pdf",
  plot = marrangeGrob(plot_list, nrow=4, ncol=3,layout_matrix = matrix(1:9, 3, 3, TRUE)),
  width = 10, height = 12
)

dev.off()

