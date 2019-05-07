### Plotting Distributions ###

rm(list = ls())
setwd("/Users/arnsteinvestre/Desktop/Studier/19_vår/STK-MAT2011/3_Files/FINAL_NOR/analysis")

PATHA= "~/code/DATA_DNV/ais_small.csv"
PATH0 = "/Users/arnsteinvestre/Desktop/Studier/19_vår/STK-MAT2011/3_Files/FINAL_NOR/computation/RecordsAnalysis/AIS_Databank_tight.csv"
PATH1 = "/Users/arnsteinvestre/Desktop/Studier/19_vår/STK-MAT2011/3_Files/FINAL_NOR/computation/RecordsAnalysis/AIS_summary_stats.csv"

library(data.table)
library(fasttime)
library(tidyr)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(GGally)
library(scales)
library(stargazer)

AIS_sum <- fread(PATH1
                 ,check.names=TRUE
                 ,showProgress=TRUE
                 ,nThread=3
                 ,stringsAsFactors=TRUE
)
AIS_sum
AIS <- AIS_sum

### Plotting for all ###
AIS_save = AIS
AIS = AIS_save
# Seconds spent
# 
AIS_save = AIS

plotframe1 = AIS
meanval1A = mean(plotframe1$Seconds_spent)

plot1A <- ggplot(data=AIS,mapping=aes(x=Seconds_spent)) + 
    geom_histogram(aes(y = ..density..),colour="black",fill="grey", position="dodge") +
    geom_density(alpha=.6, fill="royalblue") +
    labs(x="Seconds (s)", y="Density") +
    geom_vline(aes(xintercept=meanval1A),linetype="dashed",size=1) +
    ggtitle("Seconds spent on maneuver\n(Whole data set)") +
    theme_minimal() + 
    theme(plot.title = element_text(face="bold",size=15),axis.title = element_text(size=10),
          legend.position = "none",legend.text = element_blank(),legend.title = element_blank(),
          axis.text = element_text(size=10)) +
    annotate("text",x=meanval1A+200, y=0.00125, label=paste("Mean: ",format(meanval1A, digits=0),"sec"),size=5,hjust="left",vjust="bottom") 

plot1A

# Distance to CPA at yield

filter2 = ((AIS[,Distance_at_yield_1,] < 20000)|(AIS[,Distance_at_yield_2,] < 20000))&((AIS[,Distance_at_yield_1,] > 0)|(AIS[,Distance_at_yield_2,] > 0))
plotframe2 = AIS[filter2,,]
meanval2A = mean((AIS[filter2,Distance_at_yield_1] + AIS[filter2,Distance_at_yield_2])/2)

plot2A <- ggplot(data=plotframe2,mapping=aes(x=(Distance_at_yield_1+Distance_at_yield_2)/2)) + 
    geom_histogram(aes(y = ..density..),colour="black",fill="grey", position="dodge") +
    geom_density(alpha=.6, fill="royalblue") +
    labs(x="Distance (m)", y="Density") +
    geom_vline(aes(xintercept=meanval2A),linetype="dashed",size=1) +
    ggtitle("Distance to CPA at yield \n(Whole data set)") +
    theme_minimal() + 
    theme(plot.title = element_text(face="bold",size=15),axis.title = element_text(size=10),
          legend.position = "none",legend.text = element_blank(),legend.title = element_blank(),
          axis.text = element_text(size=10)) + xlim(0,40000) +
    annotate("text",x=meanval2A+3000, y=1.25e-4, label=paste("Mean: ",format(meanval2A, digits=0),"m"),size=5,hjust="left",vjust="bottom") 

plot2A

# Passing Distance

filter3 = (AIS[,Passing_distance,] < 10000)&(AIS[,Passing_distance,] > 0)
plotframe3 = AIS[filter3,,]
meanval3A = mean(plotframe3[,Passing_distance,])

plot3A <- ggplot(data=plotframe3,mapping=aes(x=Passing_distance)) + 
    geom_histogram(aes(y = ..density..),colour="black",fill="grey", position="dodge") +
    geom_density(alpha=.6, fill="royalblue") +
    labs(x="Distance (m)", y="Density") +
    geom_vline(aes(xintercept=meanval3A),linetype="dashed",size=1) +
    ggtitle("Passing distance \n(Whole data set)") +
    theme_minimal() + 
    theme(plot.title = element_text(face="bold",size=15),axis.title = element_text(size=10),
          legend.position = "none",legend.text = element_blank(),legend.title = element_blank(),
          axis.text = element_text(size=10)) + xlim(0,7500) +
    annotate("text",x=meanval3A+1000, y=7e-4, label=paste("Mean: ",format(meanval3A, digits=0),"m"),size=5,hjust="left",vjust="bottom") 

plot3A

# Approach speed
# 
# 
filter4 = (AIS[,Approach_speed,] < 1000000000)&(AIS[,Approach_speed,] > 0)
plotframe4 = AIS[filter4,,]
meanval4A = mean(plotframe4[,Approach_speed,])

plot4A <- ggplot(data=plotframe4,mapping=aes(x=Approach_speed)) + 
    geom_histogram(aes(y = ..density..),colour="black",fill="grey", position="dodge") +
    geom_density(alpha=.6, fill="royalblue") +
    labs(x="Speed (m/s)", y="Density") +
    geom_vline(aes(xintercept=meanval4A),linetype="dashed",size=1) +
    ggtitle("Approach speed\n(Whole data set)") +
    theme_minimal() + 
    theme(plot.title = element_text(face="bold",size=15),axis.title = element_text(size=10),
          legend.position = "none",legend.text = element_blank(),legend.title = element_blank(),
          axis.text = element_text(size=10)) +
    annotate("text",x=meanval4A+15, y=0.05, label=paste("Mean: \n",format(meanval4A, digits=3,nsmall=1),"m/s"),size=5,hjust="left",vjust="top") 

plot4A

# Max COG change
colnames(AIS)


filter5 = ((AIS[,Max_course_change_1,] < 100)|(AIS[,Max_course_change_2,] < 100))&((AIS[,Max_course_change_1,] > 0)|(AIS[,Max_course_change_2,] > 0))
plotframe5 = AIS[filter5,,]
meanval5A = mean((plotframe5[,Max_course_change_1] + plotframe5[,Max_course_change_2])/2)

plot5A <- ggplot(data=plotframe5,mapping=aes(x=(Max_course_change_1+Max_course_change_2)/2)) + 
    geom_histogram(aes(y = ..density..),colour="black",fill="grey", position="dodge") +
    geom_density(alpha=.6, fill="royalblue") +
    labs(x="Course change (°)", y="Density") +
    geom_vline(aes(xintercept=meanval5A),linetype="dashed",size=1) +
    ggtitle("Max course change \n(Whole data set)") +
    theme_minimal() + 
    theme(plot.title = element_text(face="bold",size=15),axis.title = element_text(size=10),
          legend.position = "none",legend.text = element_blank(),legend.title = element_blank(),
          axis.text = element_text(size=10)) +
    annotate("text",x=meanval5A+10, y=0.035, label=paste("Mean: ",format(meanval5A, digits=3),"°"),size=5,hjust="left",vjust="bottom") 

plot5A

# Max SOG change


filter6 = ((AIS[,Max_speed_change_1,] < 20)|(AIS[,Max_speed_change_2,] < 20))&((AIS[,Max_speed_change_1,] > 0)|(AIS[,Max_speed_change_2,] > 0))
plotframe6 = AIS[filter6,,]
meanval6A = mean((plotframe6[,Max_speed_change_1] + plotframe6[,Max_speed_change_2])/2)

plot6A <- ggplot(data=plotframe6,mapping=aes(x=(Max_speed_change_1+Max_speed_change_2)/2)) + 
    geom_histogram(aes(y = ..density..),colour="black",fill="grey", position="dodge") +
    geom_density(alpha=.6, fill="royalblue") +
    labs(x="Speed change (m/s)", y="Density") +
    geom_vline(aes(xintercept=meanval6A),linetype="dashed",size=1) +
    ggtitle("Max speed change \n(Whole data set)") +
    theme_minimal() + 
    theme(plot.title = element_text(face="bold",size=15),axis.title = element_text(size=10),
          legend.position = "none",legend.text = element_blank(),legend.title = element_blank(),
          axis.text = element_text(size=10)) +
    annotate("text",x=meanval6A+2, y=0.35, label=paste("Mean: ",format(meanval6A, digits=3),"m/s"),size=5,hjust="left",vjust="bottom") 

plot6A

grid.arrange(plot1,plot2,plot3,plot4,plot5,plot6,ncol=2)

## Plotting for overtaking

AIS = AIS_save

AIS = AIS[COLREG == 1,,]
# Seconds spent

plotframe1 = AIS
meanval1B = mean(plotframe1$Seconds_spent)

plot1B <- ggplot(data=AIS,mapping=aes(x=Seconds_spent)) + 
    geom_histogram(aes(y = ..density..),colour="black",fill="grey", position="dodge",binwidth=200) +
    geom_density(alpha=.5, fill="limegreen") +
    labs(x="Seconds (s)", y="Density") +
    geom_vline(aes(xintercept=meanval1B),linetype="dashed",size=1) +
    ggtitle("Seconds spent on maneuver\n(Overtaking)") +
    theme_minimal() + 
    theme(plot.title = element_text(face="bold",size=15),axis.title = element_text(size=10),
          legend.position = "none",legend.text = element_blank(),legend.title = element_blank(),
          axis.text = element_text(size=10)) +
    annotate("text",x=meanval1B+200, y=0.0007, label=paste("Mean: ",format(meanval1B, digits=0),"sec"),size=5,hjust="left",vjust="bottom") 

plot1B

# Distance to CPA at yield

filter2 = ((AIS[,Distance_at_yield_1,] < 20000)|(AIS[,Distance_at_yield_2,] < 20000))&((AIS[,Distance_at_yield_1,] > 0)|(AIS[,Distance_at_yield_2,] > 0))
plotframe2 = AIS[filter2,,]
meanval2B = mean((AIS[filter2,Distance_at_yield_1] + AIS[filter2,Distance_at_yield_2])/2)

plot2B <- ggplot(data=plotframe2,mapping=aes(x=(Distance_at_yield_1+Distance_at_yield_2)/2)) + 
    geom_histogram(aes(y = ..density..),colour="black",fill="grey", position="dodge",binwidth=2000) +
    geom_density(alpha=.5, fill="limegreen") +
    labs(x="Distance (m)", y="Density") +
    geom_vline(aes(xintercept=meanval2B),linetype="dashed",size=1) +
    ggtitle("Distance to CPA at yield \n(Overtaking)") +
    theme_minimal() + 
    theme(plot.title = element_text(face="bold",size=15),axis.title = element_text(size=10),
          legend.position = "none",legend.text = element_blank(),legend.title = element_blank(),
          axis.text = element_text(size=10)) +
    annotate("text",x=meanval2B+2000, y=8e-5, label=paste("Mean: ",format(meanval2B, digits=0),"m"),size=5,hjust="left",vjust="bottom") 

plot2B

# Passing Distance

filter3 = (AIS[,Passing_distance,] < 10000)&(AIS[,Passing_distance,] > 0)
plotframe3 = AIS[filter3,,]
meanval3B = mean(plotframe3[,Passing_distance,])

plot3B <- ggplot(data=plotframe3,mapping=aes(x=Passing_distance)) + 
    geom_histogram(aes(y = ..density..),colour="black",fill="grey", position="dodge",binwidth=500) +
    geom_density(alpha=.5, fill="limegreen") +
    labs(x="Distance (m)", y="Density") +
    geom_vline(aes(xintercept=meanval3B),linetype="dashed",size=1) +
    ggtitle("Passing distance \n(Overtaking)") +
    theme_minimal() + 
    theme(plot.title = element_text(face="bold",size=15),axis.title = element_text(size=10),
          legend.position = "none",legend.text = element_blank(),legend.title = element_blank(),
          axis.text = element_text(size=10)) +
    annotate("text",x=meanval3B+1000, y=5e-4, label=paste("Mean: ",format(meanval3B, digits=0),"m"),size=5,hjust="left",vjust="bottom") 

plot3B

# Approach speed

filter4 = (AIS[,Approach_speed,] < 1000000000)&(AIS[,Approach_speed,] > 0)
plotframe4 = AIS[filter4,,]
meanval4B = mean(plotframe4[,Approach_speed,])

plot4B <- ggplot(data=plotframe4,mapping=aes(x=Approach_speed)) + 
    geom_histogram(aes(y = ..density..),colour="black",fill="grey", position="dodge",binwidth=3) +
    geom_density(alpha=.5, fill="limegreen") +
    labs(x="Speed (m/s)", y="Density") +
    geom_vline(aes(xintercept=meanval4B),linetype="dashed",size=1) +
    ggtitle("Approach speed\n(Overtaking)") +
    theme_minimal() + 
    theme(plot.title = element_text(face="bold",size=15),axis.title = element_text(size=10),
          legend.position = "none",legend.text = element_blank(),legend.title = element_blank(),
          axis.text = element_text(size=10)) +
    annotate("text",x=meanval4B+10, y=0.04, label=paste("Mean: \n",format(meanval4B, digits=3,nsmall=1),"m/s"),size=5,hjust="left",vjust="bottom") 

plot4B

# Max COG change

filter5 = ((AIS[,Max_course_change_1,] < 100)|(AIS[,Max_course_change_2,] < 100))&((AIS[,Max_course_change_1,] > 0)|(AIS[,Max_course_change_2,] > 0))
plotframe5 = AIS[filter5,,]
meanval5B = mean((plotframe5[,Max_course_change_1] + plotframe5[,Max_course_change_2])/2)

plot5B <- ggplot(data=plotframe5,mapping=aes(x=(Max_course_change_1+Max_course_change_2)/2)) + 
    geom_histogram(aes(y = ..density..),colour="black",fill="grey", position="dodge",binwidth=7) +
    geom_density(alpha=.5, fill="limegreen") +
    labs(x="Course change (°)", y="Density") +
    geom_vline(aes(xintercept=meanval5B),linetype="dashed",size=1) +
    ggtitle("Max course change \n(Overtaking)") +
    theme_minimal() + 
    theme(plot.title = element_text(face="bold",size=15),axis.title = element_text(size=10),
          legend.position = "none",legend.text = element_blank(),legend.title = element_blank(),
          axis.text = element_text(size=10)) +
    annotate("text",x=meanval5B+10, y=0.027, label=paste("Mean: ",format(meanval5B, digits=3),"°"),size=5,hjust="left",vjust="bottom") 

plot5B

# Max SOG change


filter6 = ((AIS[,Max_speed_change_1,] < 20)|(AIS[,Max_speed_change_2,] < 20))&((AIS[,Max_speed_change_1,] > 0)|(AIS[,Max_speed_change_2,] > 0))
plotframe6 = AIS[filter6,,]
meanval6B = mean((plotframe6[,Max_speed_change_1] + plotframe6[,Max_speed_change_2])/2)

plot6B <- ggplot(data=plotframe6,mapping=aes(x=(Max_speed_change_1+Max_speed_change_2)/2)) + 
    geom_histogram(aes(y = ..density..),colour="black",fill="grey", position="dodge",binwidth = 1.5) +
    geom_density(alpha=.5, fill="limegreen") +
    labs(x="Speed change (m/s)", y="Density") +
    geom_vline(aes(xintercept=meanval6B),linetype="dashed",size=1) +
    ggtitle("Max speed change \n(Overtaking)") +
    theme_minimal() + 
    theme(plot.title = element_text(face="bold",size=15),axis.title = element_text(size=10),
          legend.position = "none",legend.text = element_blank(),legend.title = element_blank(),
          axis.text = element_text(size=10)) +
    annotate("text",x=meanval6B+2, y=0.15, label=paste("Mean: ",format(meanval6B, digits=3),"m/s"),size=5,hjust="left",vjust="bottom") 

plot6B

grid.arrange(plot1,plot2,plot3,plot4,plot5,plot6,ncol=2)


## Plotting for crossing
AIS = AIS_save

AIS = AIS[COLREG == 2,,]
# Seconds spent

plotframe1 = AIS
meanval1C = mean(plotframe1$Seconds_spent)

plot1C <- ggplot(data=AIS,mapping=aes(x=Seconds_spent)) + 
    geom_histogram(aes(y = ..density..),colour="black",fill="grey", position="dodge",binwidth=200) +
    geom_density(alpha=.5, fill="yellow3") +
    labs(x="Seconds (s)", y="Density") +
    geom_vline(aes(xintercept=meanval1C),linetype="dashed",size=1) +
    ggtitle("Seconds spent on maneuver\n(Crossing)") +
    theme_minimal() + 
    theme(plot.title = element_text(face="bold",size=15),axis.title = element_text(size=10),
          legend.position = "none",legend.text = element_blank(),legend.title = element_blank(),
          axis.text = element_text(size=10)) +
    annotate("text",x=meanval1C+200, y=0.001, label=paste("Mean: ",format(meanval1C, digits=0),"sec"),size=5,hjust="left",vjust="bottom") 

plot1C

# Distance to CPA at yield

filter2 = ((AIS[,Distance_at_yield_1,] < 20000)|(AIS[,Distance_at_yield_2,] < 20000))&((AIS[,Distance_at_yield_1,] > 0)|(AIS[,Distance_at_yield_2,] > 0))
plotframe2 = AIS[filter2,,]
meanval2C = mean((AIS[filter2,Distance_at_yield_1] + AIS[filter2,Distance_at_yield_2])/2)

plot2C <- ggplot(data=plotframe2,mapping=aes(x=(Distance_at_yield_1+Distance_at_yield_2)/2)) + 
    geom_histogram(aes(y = ..density..),colour="black",fill="grey", position="dodge",binwidth=2000) +
    geom_density(alpha=.5, fill="yellow3") +
    labs(x="Distance (m)", y="Density") +
    geom_vline(aes(xintercept=meanval2C),linetype="dashed",size=1) +
    ggtitle("Distance to CPA at yield \n(Crossing)") +
    theme_minimal() + 
    theme(plot.title = element_text(face="bold",size=15),axis.title = element_text(size=10),
          legend.position = "none",legend.text = element_blank(),legend.title = element_blank(),
          axis.text = element_text(size=10)) +
    annotate("text",x=meanval2C+5000, y=1.1e-4, label=paste("Mean: ",format(meanval2C, digits=0),"m"),size=5,hjust="left",vjust="bottom") 

plot2C

# Passing Distance

filter3 = (AIS[,Passing_distance,] < 10000)&(AIS[,Passing_distance,] > 0)
plotframe3 = AIS[filter3,,]
meanval3C = mean(plotframe3[,Passing_distance,])

plot3C <- ggplot(data=plotframe3,mapping=aes(x=Passing_distance)) + 
    geom_histogram(aes(y = ..density..),colour="black",fill="grey", position="dodge",binwidth=500) +
    geom_density(alpha=.5, fill="yellow3") +
    labs(x="Distance (m)", y="Density") +
    geom_vline(aes(xintercept=meanval3C),linetype="dashed",size=1) +
    ggtitle("Passing distance \n(Crossing)") +
    theme_minimal() + 
    theme(plot.title = element_text(face="bold",size=15),axis.title = element_text(size=10),
          legend.position = "none",legend.text = element_blank(),legend.title = element_blank(),
          axis.text = element_text(size=10)) +
    annotate("text",x=meanval3C+1000, y=5.5e-4, label=paste("Mean: ",format(meanval3C, digits=0),"m"),size=5,hjust="left",vjust="bottom") 

plot3C

# Approach speed

filter4 = (AIS[,Approach_speed,] < 1000000000)&(AIS[,Approach_speed,] > 0)
plotframe4 = AIS[filter4,,]
meanval4C = mean(plotframe4[,Approach_speed,])

plot4C <- ggplot(data=plotframe4,mapping=aes(x=Approach_speed)) + 
    geom_histogram(aes(y = ..density..),colour="black",fill="grey", position="dodge",binwidth=3) +
    geom_density(alpha=.5, fill="yellow3") +
    labs(x="Speed (m/s)", y="Density") +
    geom_vline(aes(xintercept=meanval4C),linetype="dashed",size=1) +
    ggtitle("Approach speed\n(Crossing)") +
    theme_minimal() + 
    theme(plot.title = element_text(face="bold",size=15),axis.title = element_text(size=10),
          legend.position = "none",legend.text = element_blank(),legend.title = element_blank(),
          axis.text = element_text(size=10)) +
    annotate("text",x=meanval4C+10, y=0.04, label=paste("Mean: \n",format(meanval4C, digits=3,nsmall=1),"m/s"),size=5,hjust="left",vjust="bottom") 

plot4C

# Max COG change

filter5 = ((AIS[,Max_course_change_1,] < 100)|(AIS[,Max_course_change_2,] < 100))&((AIS[,Max_course_change_1,] > 0)|(AIS[,Max_course_change_2,] > 0))
plotframe5 = AIS[filter5,,]
meanval5C = mean((plotframe5[,Max_course_change_1] + plotframe5[,Max_course_change_2])/2)

plot5C <- ggplot(data=plotframe5,mapping=aes(x=(Max_course_change_1+Max_course_change_2)/2)) + 
    geom_histogram(aes(y = ..density..),colour="black",fill="grey", position="dodge",binwidth=7) +
    geom_density(alpha=.5, fill="yellow3") +
    labs(x="Course change (°)", y="Density") +
    geom_vline(aes(xintercept=meanval5C),linetype="dashed",size=1) +
    ggtitle("Max course change \n(Crossing)") +
    theme_minimal() + 
    theme(plot.title = element_text(face="bold",size=15),axis.title = element_text(size=10),
          legend.position = "none",legend.text = element_blank(),legend.title = element_blank(),
          axis.text = element_text(size=10)) +
    annotate("text",x=meanval5C+10, y=0.025, label=paste("Mean: ",format(meanval5C, digits=3,nsmall=1),"°"),size=5,hjust="left",vjust="bottom") 

plot5C

# Max SOG change


filter6 = ((AIS[,Max_speed_change_1,] < 20)|(AIS[,Max_speed_change_2,] < 20))&((AIS[,Max_speed_change_1,] > 0)|(AIS[,Max_speed_change_2,] > 0))
plotframe6 = AIS[filter6,,]
meanval6C = mean((plotframe6[,Max_speed_change_1] + plotframe6[,Max_speed_change_2])/2)

plot6C <- ggplot(data=plotframe6,mapping=aes(x=(Max_speed_change_1+Max_speed_change_2)/2)) + 
    geom_histogram(aes(y = ..density..),colour="black",fill="grey", position="dodge",binwidth = 1.5) +
    geom_density(alpha=.5, fill="yellow3") +
    labs(x="Speed change (m/s)", y="Density") +
    geom_vline(aes(xintercept=meanval6C),linetype="dashed",size=1) +
    ggtitle("Max speed change \n(Crossing)") +
    theme_minimal() + 
    theme(plot.title = element_text(face="bold",size=15),axis.title = element_text(size=10),
          legend.position = "none",legend.text = element_blank(),legend.title = element_blank(),
          axis.text = element_text(size=10)) +
    annotate("text",x=meanval6C+2, y=0.2, label=paste("Mean: ",format(meanval6C, digits=3),"m/s"),size=5,hjust="left",vjust="bottom") 

plot6C

grid.arrange(plot1,plot2,plot3,plot4,plot5,plot6,ncol=2)

## Plotting for head-to-head
AIS = AIS_save

AIS = AIS[COLREG == 3,,]
# Seconds spent

filter1 = AIS[,Seconds_spent,] < 350000
plotframe1 = AIS[filter1,]
meanval1D = mean(plotframe1$Seconds_spent)

plot1D <- ggplot(data=plotframe1,mapping=aes(x=Seconds_spent)) + 
    geom_histogram(aes(y = ..density..),colour="black",fill="grey", position="dodge",binwidth=200) +
    geom_density(alpha=.5, fill="red") +
    labs(x="Seconds (s)", y="Density") +
    geom_vline(aes(xintercept=meanval1D),linetype="dashed",size=1) +
    ggtitle("Seconds spent on maneuver\n(Head-on)") +
    theme_minimal() + 
    theme(plot.title = element_text(face="bold",size=15),axis.title = element_text(size=10),
          legend.position = "none",legend.text = element_blank(),legend.title = element_blank(),
          axis.text = element_text(size=10)) + xlim(0,3500) +
    annotate("text",x=meanval1D+300, y=0.0025, label=paste("Mean: ",format(meanval1D, digits=0),"sec"),size=5,hjust="left",vjust="bottom") 

plot1D

# Distance to CPA at yield

filter2 = ((AIS[,Distance_at_yield_1,] < 20000)|(AIS[,Distance_at_yield_2,] < 20000))&((AIS[,Distance_at_yield_1,] > 0)|(AIS[,Distance_at_yield_2,] > 0))
plotframe2 = AIS[filter2,,]
meanval2D = mean((AIS[filter2,Distance_at_yield_1] + AIS[filter2,Distance_at_yield_2])/2)

plot2D <- ggplot(data=plotframe2,mapping=aes(x=(Distance_at_yield_1+Distance_at_yield_2)/2)) + 
    geom_histogram(aes(y = ..density..),colour="black",fill="grey", position="dodge",binwidth=2000) +
    geom_density(alpha=.5, fill="red") +
    labs(x="Distance (m)", y="Density") +
    geom_vline(aes(xintercept=meanval2D),linetype="dashed",size=1) +
    ggtitle("Distance to CPA at yield \n(Head-on)") +
    theme_minimal() + 
    theme(plot.title = element_text(face="bold",size=15),axis.title = element_text(size=10),
          legend.position = "none",legend.text = element_blank(),legend.title = element_blank(),
          axis.text = element_text(size=10) ) + xlim(0,30000) +
    annotate("text",x=meanval2D+5000, y=3.3e-4, label=paste("Mean: ",format(meanval2D, digits=0),"m"),size=5,hjust="left",vjust="bottom") 

plot2D

# Passing Distance

filter3 = (AIS[,Passing_distance,] < 10000)&(AIS[,Passing_distance,] > 0)
plotframe3 = AIS[filter3,,]
meanval3D = mean(plotframe3[,Passing_distance,])

plot3D <- ggplot(data=plotframe3,mapping=aes(x=Passing_distance)) + 
    geom_histogram(aes(y = ..density..),colour="black",fill="grey", position="dodge",binwidth=500) +
    geom_density(alpha=.5, fill="red") +
    labs(x="Distance (m)", y="Density") +
    geom_vline(aes(xintercept=meanval3D),linetype="dashed",size=1) +
    ggtitle("Passing distance \n(Head-on)") +
    theme_minimal() + 
    theme(plot.title = element_text(face="bold",size=15),axis.title = element_text(size=10),
          legend.position = "none",legend.text = element_blank(),legend.title = element_blank(),
          axis.text = element_text(size=10)) +
    annotate("text",x=meanval3D+1000, y=9e-4, label=paste("Mean: ",format(meanval3D, digits=0),"m"),size=5,hjust="left",vjust="bottom") 

plot3D

# Approach speed

filter4 = (AIS[,Approach_speed,] < 1000000000)&(AIS[,Approach_speed,] > 0)
plotframe4 = AIS[filter4,,]
meanval4D = mean(plotframe4[,Approach_speed,])

plot4D <- ggplot(data=plotframe4,mapping=aes(x=Approach_speed)) + 
    geom_histogram(aes(y = ..density..),colour="black",fill="grey", position="dodge",binwidth=3) +
    geom_density(alpha=.5, fill="red") +
    labs(x="Speed (m/s)", y="Density") +
    geom_vline(aes(xintercept=meanval4D),linetype="dashed",size=1) +
    ggtitle("Approach speed\n(Head-on)") +
    theme_minimal() + 
    theme(plot.title = element_text(face="bold",size=15),axis.title = element_text(size=10),
          legend.position = "none",legend.text = element_blank(),legend.title = element_blank(),
          axis.text = element_text(size=10)) +
    annotate("text",x=meanval4D+10, y=0.05, label=paste("Mean: \n",format(meanval4D, digits=3,nsmall=1),"m/s"),size=5,hjust="left",vjust="bottom") 

plot4D

# Max COG change

filter5 = ((AIS[,Max_course_change_1,] < 100)|(AIS[,Max_course_change_2,] < 100))&((AIS[,Max_course_change_1,] > 0)|(AIS[,Max_course_change_2,] > 0))
plotframe5 = AIS[filter5,,]
meanval5D = mean((plotframe5[,Max_course_change_1] + plotframe5[,Max_course_change_2])/2)

plot5D <- ggplot(data=plotframe5,mapping=aes(x=(Max_course_change_1+Max_course_change_2)/2)) + 
    geom_histogram(aes(y = ..density..),colour="black",fill="grey", position="dodge",binwidth=7) +
    geom_density(alpha=.5, fill="red") +
    labs(x="Course change (°)", y="Density") +
    geom_vline(aes(xintercept=meanval5D),linetype="dashed",size=1) +
    ggtitle("Max course change \n(Head-on)") +
    theme_minimal() + 
    theme(plot.title = element_text(face="bold",size=15),axis.title = element_text(size=10),
          legend.position = "none",legend.text = element_blank(),legend.title = element_blank(),
          axis.text = element_text(size=10)) + xlim(0,100) + 
    annotate("text",x=meanval5D+10, y=0.1, label=paste("Mean: ",format(meanval5D, digits=3,nsmall=1),"°"),size=5,hjust="left",vjust="bottom") 

plot5D

# Max SOG change


filter6 = ((AIS[,Max_speed_change_1,] < 20)|(AIS[,Max_speed_change_2,] < 20))&((AIS[,Max_speed_change_1,] > 0)|(AIS[,Max_speed_change_2,] > 0))
plotframe6 = AIS[filter6,,]
meanval6D = mean((plotframe6[,Max_speed_change_1] + plotframe6[,Max_speed_change_2])/2)

plot6D <- ggplot(data=plotframe6,mapping=aes(x=(Max_speed_change_1+Max_speed_change_2)/2)) + 
    geom_histogram(aes(y = ..density..),colour="black",fill="grey", position="dodge") +
    geom_density(alpha=.5, fill="red") +
    labs(x="Speed change (m/s)", y="Density") +
    geom_vline(aes(xintercept=meanval6D),linetype="dashed",size=1) +
    ggtitle("Max speed change \n(Head-on)") +
    theme_minimal() + 
    theme(plot.title = element_text(face="bold",size=15),axis.title = element_text(size=10),
          legend.position = "none",legend.text = element_blank(),legend.title = element_blank(),
          axis.text = element_text(size=10)) +
    annotate("text",x=meanval6D+3, y=2, label=paste("Mean: ",format(meanval6D, digits=3),"m/s"),size=5,hjust="left",vjust="bottom") 

plot6D

grid.arrange(plot1,plot2,plot3,plot4,plot5,plot6,ncol=2)

### Updated plot combinations

grid.arrange(plot1A,plot1B,plot1C,plot1D,ncol=2)
grid.arrange(plot2A,plot2B,plot2C,plot2D,ncol=2)
grid.arrange(plot3A,plot3B,plot3C,plot3D,ncol=2)
grid.arrange(plot4A,plot4B,plot4C,plot4D,ncol=2)
grid.arrange(plot5A,plot5B,plot5C,plot5D,ncol=2)
grid.arrange(plot6A,plot6B,plot6C,plot6D,ncol=2)
