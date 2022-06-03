#Figure3
getwd()
EBSAI_dd=read.csv("EBSAI_fullMH.csv",header=TRUE)
View(EBSAI_dd)
library(dplyr)
detach(package:plyr)
library(tidyverse)
#devtools::install_github("norment/normentR")
library(normentR)
library(janitor)
#devtools::install_github("eliocamp/ggnewscale")
library(ggnewscale)

#EBSAI_dd=read.csv("/Users/ingridspies/Documents/Carolyn_Cod/EBSAI_fullMH.csv",header=TRUE)

Ofst=data.frame(cbind(EBSAI_dd$X[which(EBSAI_dd$EBSAI_fst_dxy.upper_outlier=="HighFstOutlier")],
                      EBSAI_dd$fst[which(EBSAI_dd$EBSAI_fst_dxy.upper_outlier=="HighFstOutlier")],
                      EBSAI_dd$chr[which(EBSAI_dd$EBSAI_fst_dxy.upper_outlier=="HighFstOutlier")]))
colnames(Ofst)=c("X","fst","chr")

Odxy=data.frame(cbind(EBSAI_dd$X[which(EBSAI_dd$EBSAI_dxy_upper5=="TRUE"&EBSAI_dd$fst>0.023)],
                      EBSAI_dd$weightedAverageDxy[which(EBSAI_dd$EBSAI_dxy_upper5=="TRUE"&EBSAI_dd$fst>0.023)]),rep(.075,length(EBSAI_dd$X[which(EBSAI_dd$EBSAI_dxy_upper5=="TRUE"&EBSAI_dd$fst>0.023)])))
colnames(Odxy)=c("X","dxy","val")

axis_set <- EBSAI_dd %>% 
 group_by(chr) %>% 
 summarize(center = mean(X))

axis_set3 <- Ofst %>% 
 group_by(chr) %>% 
 summarize(center = mean(X))

manhplot1 <- ggplot(EBSAI_dd, aes(x = X, y = EBSAI_fst_dxy.weightedAverageFST, 
                                  color = as_factor(chr)),size=.05)  + 
 ggtitle("Eastern Bering Sea vs. Aleutian Islands")+
 geom_point(alpha = 0.5,size=.5) +
 scale_x_continuous(label = axis_set$chr, breaks = axis_set$center) +
 scale_y_continuous(expand = c(0,0), limits = c(0.0, 0.1)) +
 scale_color_manual(values = rep(c("gray72", "grey30"), unique(length(axis_set$chr)))) +
 new_scale_color()+
 geom_point(data=Ofst,aes(x=X,y=fst,color = as_factor(chr)),alpha=0.5,size=0.8)+
 scale_color_manual(values = rep(c("red", "darkred"), unique(length(axis_set3$chr)))) +
 scale_size_continuous(range = c(0,3)) +
 ylab(expression(paste("Weighted Average ", italic(F)[ST])))+
 labs(x = "Linkage Group",size=12) + 
 theme_minimal() +
 theme( 
  legend.position = "none",
  panel.border = element_blank(),
  panel.grid.major.x = element_blank(),
  panel.grid.minor.x = element_blank(),
  axis.text.x = element_text(angle = 0, size = 10, vjust = 9)
 )
print(manhplot1)