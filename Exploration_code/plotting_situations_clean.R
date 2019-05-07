### Situation plotting clean ####

rm(list = ls())
setwd("/Users/arnsteinvestre/Desktop/Studier/19_vår/STK-MAT2011/3_Files/FINAL_NOR/analysis")

PATH1 = "/Users/arnsteinvestre/Desktop/Studier/19_vår/STK-MAT2011/3_Files/FINAL_NOR/computation/RecordsAnalysis/AIS_with_stats.csv"
PATH2 = "/Users/arnsteinvestre/Desktop/Studier/19_vår/STK-MAT2011/3_Files/FINAL_NOR/computation/RecordsAnalysis/AIS_summary_stats.csv"


library(data.table)
library(fasttime)
library(tidyr)
library(dplyr)
library(ggplot2)
library(gridExtra)


AIS_sum <- fread(PATH2
                 ,check.names=TRUE
                 ,showProgress=TRUE
                 ,nThread=3
                 ,stringsAsFactors=TRUE
)
AIS_sum

AIS_NOR <- fread(PATH1
                 ,check.names=TRUE
                 ,showProgress=TRUE
                 ,nThread=3
                 ,stringsAsFactors=TRUE
)

# RELOAD COPY #
AIS <- AIS_NOR

table1 <- AIS[,c("ID_1","MT_datetime","LAT_1","LON_1","COG_1","SOG_1","dSOG_1","dCOG_1","Distance_to_CPA_1","Time_to_CPA","CPA_distance","Category_1","Length_1","Heading_1","Time_1","Time_datetime_1", "Merge_times", "Situations_1"),]
table1 <- table1[,Ship:=1,]

table2 <- AIS[,c("ID_2","MT_datetime","LAT_2","LON_2","COG_2","SOG_2","dSOG_2","dCOG_2","Distance_to_CPA_2","Time_to_CPA","CPA_distance","Category_2","Length_2","Heading_2","Time_2","Time_datetime_2", "Merge_times", "Situations_2"),]
table2 <- table2[,Ship:=2,]

colnames(table1) <- c("ID","MT_datetime","LAT","LON","COG","SOG","dSOG","dCOG","Distance_to_CPA","Time_to_CPA","CPA_distance","Category","Length","Heading","Time","Time_datetime", "Merge_times", "Situations","Ship")
colnames(table2) <- c("ID","MT_datetime","LAT","LON","COG","SOG","dSOG","dCOG","Distance_to_CPA","Time_to_CPA","CPA_distance","Category","Length","Heading","Time","Time_datetime", "Merge_times", "Situations","Ship")

AIS_final <- rbind(table1,table2)
AIS_final

AIS <- AIS_final

AIS <- as.data.table(AIS %>% group_by(Situations) %>% mutate(Intervals = (Time - min(Time)) %/% 600))
AIS[,ColorCode:= Intervals %% 2,]
Time_tresh = 3600
DMC = 1.1132e+5

rename <- function(x){
    if (x < 10) {
        return(name <- paste('000',x,'plot.png',sep=''))
    }
    if (x < 100 && x >= 10) {
        return(name <- paste('00',x,'plot.png', sep=''))
    }
    if (x < 1000 && x >= 100) {
        return(name <- paste('0', x,'plot.png', sep=''))
    }
    if (x >= 1000) {
        return(name <- paste(x,'plot.png', sep=''))
    }
}
### Plotting for separate situations

i = 2494
plotframe = AIS[Situations == i,,]
IDs = unique(plotframe[,ID,])
setkey(plotframe,Time)
plotframe <- plotframe[(ID == IDs[1])&(abs(Time - min(Time[ID == IDs[1]])) < 10),ID:= 2,]
plotframe <- plotframe[(ID == IDs[2])&(abs(Time - min(Time[ID == IDs[2]])) < 10),ID:= 2,]
plotframe <- plotframe[,,]
plotA <- ggplot(data=plotframe,aes(x=LON,y=LAT,colour=as.factor(ID), alpha = as.factor(ColorCode), shape = as.factor(ID), size = as.factor(ID))) + geom_point() +
    ggtitle(paste("Situation number ",i, "\nExample: Overtaking")) + 
    theme(panel.background = element_rect(fill = "grey95", colour = "grey60",size = 2, linetype = "solid"),
          legend.text = element_text(size=12),legend.title = element_text(size=15, face = "bold", hjust = 0.3),
          plot.title = element_text(face="bold",size=12)) + 
    scale_color_manual(name = "ID", values = c("green2", "red3", "darkblue")) + 
    scale_shape_manual(name = "Shape", values = c(18,16,16)) +
    scale_size_manual(name = "Size", values = c(6,1.5,1.5)) +
    scale_alpha_manual(name = "Time", values = c(0.8,0.05)) + 
    guides(size = FALSE,
           alpha = FALSE,
           shape = FALSE,
           color = FALSE) + 
    annotate("point",x=plotframe[AIS_sum[Situation==i,Yield_start_record,]*2,LON,],y=plotframe[AIS_sum[Situation==i,Yield_start_record,]*2,LAT,],shape=15,size=4,colour="brown") + 
    annotate("point",x=plotframe[AIS_sum[Situation==i,Yield_start_record,]*2+1,LON,],y=plotframe[AIS_sum[Situation==i,Yield_start_record,]*2+1,LAT,],shape=15,size=4,colour="brown") +
    annotate("point",x=plotframe[AIS_sum[Situation==i,Yield_finish_record,]*2,LON,],y=plotframe[AIS_sum[Situation==i,Yield_finish_record,]*2,LAT,],shape=15,size=4,colour="blue") + 
    annotate("point",x=plotframe[AIS_sum[Situation==i,Yield_finish_record,]*2-1,LON,],y=plotframe[AIS_sum[Situation==i,Yield_finish_record,]*2+1,LAT,],shape=15,size=4,colour="blue") +
    annotate("text",x=plotframe[AIS_sum[Situation==i,Yield_start_record,]*2,LON,],y=plotframe[AIS_sum[Situation==i,Yield_start_record,]*2,LAT,],label="t[1]",parse=TRUE,vjust=2) + 
    annotate("text",x=plotframe[AIS_sum[Situation==i,Yield_start_record,]*2+1,LON,],y=plotframe[AIS_sum[Situation==i,Yield_start_record,]*2+1,LAT,],label="t[1]",parse=TRUE,vjust=-1) +
    annotate("text",x=plotframe[AIS_sum[Situation==i,Yield_finish_record,]*2,LON,],y=plotframe[AIS_sum[Situation==i,Yield_finish_record,]*2,LAT,],label="t[F]",parse=TRUE,vjust=-1) + 
    annotate("text",x=plotframe[AIS_sum[Situation==i,Yield_finish_record,]*2-1,LON,],y=plotframe[AIS_sum[Situation==i,Yield_finish_record,]*2+1,LAT,],label="t[F]",parse=TRUE,vjust=2) +
    annotate("point",x=plotframe[1,LON,],y=plotframe[1,LAT,],shape=18,size=6,colour="green") + 
    annotate("point",x=plotframe[2,LON,],y=plotframe[2,LAT,],shape=18,size=6,colour="green")
plotA

i = 2495
plotframe = AIS[Situations == i,,]
IDs = unique(plotframe[,ID,])
setkey(plotframe,Time)
plotframe <- plotframe[(ID == IDs[1])&(abs(Time - min(Time[ID == IDs[1]])) < 10),ID:= 2,]
plotframe <- plotframe[(ID == IDs[2])&(abs(Time - min(Time[ID == IDs[2]])) < 10),ID:= 2,]
plotB <- ggplot(data=plotframe,aes(x=LON,y=LAT,colour=as.factor(ID), alpha = as.factor(ColorCode), shape = as.factor(ID), size = as.factor(ID))) + geom_point() +
    ggtitle(paste("Situation number ",i, "\nExample: Crossing")) + 
    theme(panel.background = element_rect(fill = "grey95", colour = "grey60",size = 2, linetype = "solid"),
          legend.text = element_text(size=12),legend.title = element_text(size=15, face = "bold", hjust = 0.3),
          plot.title = element_text(face="bold",size=12)) + 
    scale_color_manual(name = "ID", values = c("green2", "red3", "darkblue")) + 
    scale_shape_manual(name = "Shape", values = c(18,16,16)) +
    scale_size_manual(name = "Size", values = c(6,1.5,1.5)) +
    scale_alpha_manual(name = "Time", values = c(1,0.02)) + 
    guides(size = FALSE,
           alpha = FALSE,
           shape = FALSE,
           color = FALSE) + 
    annotate("point",x=plotframe[AIS_sum[Situation==i,Yield_start_record,]*2,LON,],y=plotframe[AIS_sum[Situation==i,Yield_start_record,]*2,LAT,],shape=15,size=4,colour="brown") + 
    annotate("point",x=plotframe[AIS_sum[Situation==i,Yield_start_record,]*2+1,LON,],y=plotframe[AIS_sum[Situation==i,Yield_start_record,]*2+1,LAT,],shape=15,size=4,colour="brown") +
    annotate("point",x=plotframe[AIS_sum[Situation==i,Yield_finish_record,]*2,LON,],y=plotframe[AIS_sum[Situation==i,Yield_finish_record,]*2,LAT,],shape=15,size=4,colour="blue") + 
    annotate("point",x=plotframe[AIS_sum[Situation==i,Yield_finish_record,]*2-1,LON,],y=plotframe[AIS_sum[Situation==i,Yield_finish_record,]*2+1,LAT,],shape=15,size=4,colour="blue") +
    annotate("text",x=plotframe[AIS_sum[Situation==i,Yield_start_record,]*2,LON,],y=plotframe[AIS_sum[Situation==i,Yield_start_record,]*2,LAT,],label="t[1]",parse=TRUE,vjust=-1,hjust=-1) + 
    annotate("text",x=plotframe[AIS_sum[Situation==i,Yield_start_record,]*2+1,LON,],y=plotframe[AIS_sum[Situation==i,Yield_start_record,]*2+1,LAT,],label="t[1]",parse=TRUE,vjust=-1) +
    annotate("text",x=plotframe[AIS_sum[Situation==i,Yield_finish_record,]*2,LON,],y=plotframe[AIS_sum[Situation==i,Yield_finish_record,]*2,LAT,],label="t[F]",parse=TRUE,vjust=-1) + 
    annotate("text",x=plotframe[AIS_sum[Situation==i,Yield_finish_record,]*2-1,LON,],y=plotframe[AIS_sum[Situation==i,Yield_finish_record,]*2+1,LAT,],label="t[F]",parse=TRUE,vjust=-1) +
    annotate("point",x=plotframe[1,LON,],y=plotframe[1,LAT,],shape=18,size=6,colour="green") + 
    annotate("point",x=plotframe[2,LON,],y=plotframe[2,LAT,],shape=18,size=6,colour="green")
plotB


i = 15
plotframe = AIS[Situations == i,,]
IDs = unique(plotframe[,ID,])
setkey(plotframe,Time)
plotframe <- plotframe[(ID == IDs[1])&(abs(Time - min(Time[ID == IDs[1]])) < 10),ID:= 2,]
plotframe <- plotframe[(ID == IDs[2])&(abs(Time - min(Time[ID == IDs[2]])) < 10),ID:= 2,]
plotC <- ggplot(data=plotframe,aes(x=LON,y=LAT,colour=as.factor(ID), alpha = as.factor(ColorCode), shape = as.factor(ID), size = as.factor(ID))) + geom_point() +
    ggtitle(paste("Situation number ",i, "\nExample: Head-on")) + 
    theme(panel.background = element_rect(fill = "grey95", colour = "grey60",size = 2, linetype = "solid"),
          legend.text = element_text(size=12),legend.title = element_text(size=15, face = "bold", hjust = 0.3),
          plot.title = element_text(face="bold",size=12)) + 
    scale_color_manual(name = "ID", values = c("green2", "red3", "darkblue")) + 
    scale_shape_manual(name = "Shape", values = c(18,16,16)) +
    scale_size_manual(name = "Size", values = c(6,1.5,1.5)) +
    scale_alpha_manual(name = "Time", values = c(1,0.02)) + 
    guides(size = FALSE,
           alpha = FALSE,
           shape = FALSE,
           color = FALSE) + 
    annotate("point",x=plotframe[ID == IDs[1]][AIS_sum[Situation==i,Yield_start_record,],LON,],y=plotframe[ID == IDs[1]][AIS_sum[Situation==i,Yield_start_record,],LAT,],shape=15,size=4,colour="brown") + 
    annotate("point",x=plotframe[ID == IDs[2]][AIS_sum[Situation==i,Yield_start_record,],LON,],y=plotframe[ID == IDs[2]][AIS_sum[Situation==i,Yield_start_record,],LAT,],shape=15,size=4,colour="brown") +
    annotate("point",x=plotframe[ID == IDs[1]][AIS_sum[Situation==i,Yield_finish_record,],LON,],y=plotframe[ID == IDs[1]][AIS_sum[Situation==i,Yield_finish_record,],LAT,],shape=15,size=4,colour="blue") + 
    annotate("point",x=plotframe[ID == IDs[2]][AIS_sum[Situation==i,Yield_finish_record,],LON,],y=plotframe[ID == IDs[2]][AIS_sum[Situation==i,Yield_finish_record,],LAT,],shape=15,size=4,colour="blue") +
    annotate("text",x=plotframe[ID == IDs[1]][AIS_sum[Situation==i,Yield_start_record,],LON,],y=plotframe[ID == IDs[1]][AIS_sum[Situation==i,Yield_start_record,],LAT,],label="t[1]",parse=TRUE,vjust=2) + 
    annotate("text",x=plotframe[ID == IDs[2]][AIS_sum[Situation==i,Yield_start_record,],LON,],y=plotframe[ID == IDs[2]][AIS_sum[Situation==i,Yield_start_record,],LAT,],label="t[1]",parse=TRUE,vjust=-1,hjust=2) +
    annotate("text",x=plotframe[ID == IDs[1]][AIS_sum[Situation==i,Yield_finish_record,],LON,],y=plotframe[ID == IDs[1]][AIS_sum[Situation==i,Yield_finish_record,],LAT,],label="t[F]",parse=TRUE,vjust=2,hjust=-1) + 
    annotate("text",x=plotframe[ID == IDs[2]][AIS_sum[Situation==i,Yield_finish_record,],LON,],y=plotframe[ID == IDs[2]][AIS_sum[Situation==i,Yield_finish_record,],LAT,],label="t[F]",parse=TRUE,vjust=-1) +
    annotate("point",x=plotframe[ID == IDs[1]][1,LON,],y=plotframe[ID == IDs[1]][1,LAT,],shape=18,size=6,colour="green") + 
    annotate("point",x=plotframe[ID == IDs[2]][1,LON,],y=plotframe[ID == IDs[2]][1,LAT,],shape=18,size=6,colour="green")
plotC


i = 18
plotframe = AIS[Situations == i,,]
IDs = unique(plotframe[,ID,])
setkey(plotframe,Time)
plotframe <- plotframe[(ID == IDs[1])&(abs(Time - min(Time[ID == IDs[1]])) < 10),ID:= 2,]
plotframe <- plotframe[(ID == IDs[2])&(abs(Time - min(Time[ID == IDs[2]])) < 10),ID:= 2,]
plotD <- ggplot(data=plotframe,aes(x=LON,y=LAT,colour=as.factor(ID), alpha = as.factor(ColorCode), shape = as.factor(ID), size = as.factor(ID))) + geom_point() +
    ggtitle(paste("Situation number ",i, "\nExample: Not distinct")) + 
    theme(panel.background = element_rect(fill = "grey95", colour = "grey60",size = 2, linetype = "solid"),
          legend.text = element_text(size=12),legend.title = element_text(size=15, face = "bold", hjust = 0.3),
          plot.title = element_text(face="bold",size=12)) + 
    scale_color_manual(name = "ID", values = c("green2", "red3", "darkblue")) + 
    scale_shape_manual(name = "Shape", values = c(18,16,16)) +
    scale_size_manual(name = "Size", values = c(6,1.5,1.5)) +
    scale_alpha_manual(name = "Time", values = c(1,0.02)) + 
    guides(size = FALSE,
           alpha = FALSE,
           shape = FALSE,
           color = FALSE) + 
    annotate("point",x=plotframe[AIS_sum[Situation==i,Yield_start_record,]*2,LON,],y=plotframe[AIS_sum[Situation==i,Yield_start_record,]*2,LAT,],shape=15,size=4,colour="brown") + 
    annotate("point",x=plotframe[AIS_sum[Situation==i,Yield_start_record,]*2+1,LON,],y=plotframe[AIS_sum[Situation==i,Yield_start_record,]*2+1,LAT,],shape=15,size=4,colour="brown") +
    annotate("point",x=plotframe[AIS_sum[Situation==i,Yield_finish_record,]*2,LON,],y=plotframe[AIS_sum[Situation==i,Yield_finish_record,]*2,LAT,],shape=15,size=4,colour="blue") + 
    annotate("point",x=plotframe[AIS_sum[Situation==i,Yield_finish_record,]*2-1,LON,],y=plotframe[AIS_sum[Situation==i,Yield_finish_record,]*2+1,LAT,],shape=15,size=4,colour="blue") +
    annotate("text",x=plotframe[AIS_sum[Situation==i,Yield_start_record,]*2,LON,],y=plotframe[AIS_sum[Situation==i,Yield_start_record,]*2,LAT,],label="t[1]",parse=TRUE,vjust=-1) + 
    annotate("text",x=plotframe[AIS_sum[Situation==i,Yield_start_record,]*2+1,LON,],y=plotframe[AIS_sum[Situation==i,Yield_start_record,]*2+1,LAT,],label="t[1]",parse=TRUE,vjust=0,hjust=-1) +
    annotate("text",x=plotframe[AIS_sum[Situation==i,Yield_finish_record,]*2,LON,],y=plotframe[AIS_sum[Situation==i,Yield_finish_record,]*2,LAT,],label="t[F]",parse=TRUE,vjust=-1) + 
    annotate("text",x=plotframe[AIS_sum[Situation==i,Yield_finish_record,]*2-1,LON,],y=plotframe[AIS_sum[Situation==i,Yield_finish_record,]*2+1,LAT,],label="t[F]",parse=TRUE,hjust=-1) +
    annotate("point",x=plotframe[1,LON,],y=plotframe[1,LAT,],shape=18,size=6,colour="green") + 
    annotate("point",x=plotframe[2,LON,],y=plotframe[2,LAT,],shape=18,size=6,colour="green")
plotD


plotFIN <- grid.arrange(plotA,plotB,plotC,plotD,ncol=2)
plotFIN





