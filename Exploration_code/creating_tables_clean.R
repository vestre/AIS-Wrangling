### Creating tables ###

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

### Groupings ###
#2 Chemical tankers
#3 Gas tankers
#4 Bulk carriers
#14 Crude oil tankers
#15 Oil product tankers

#5 General cargo ships
#6 Container ships 
#7 Ro-Ro cargo ships
#8 Refrigerated cargo ships 

#10 Offshore supply ships
#11 Other service offshore vessels

#12 Other activities

#13 Fishing vessels

#16 Passenger ships

#17 Cruise ships

### Summary statistics ###

AIS[,Cat1 := Category_1,]
AIS[,Cat2 := Category_2,]

AIS[,Category_1 := NULL,]
AIS[,Category_2 := NULL,]

AIS[,Category_1 := "NULL",]
AIS[,Category_2 := "NULL",]

AIS[Cat1 %in% c(2,3,4,14,15),Category_1:= "A",]
AIS[Cat2 %in% c(2,3,4,14,15),Category_2:= "A",]

AIS[Cat1 %in% c(10,11),Category_1:= "B",]
AIS[Cat2 %in% c(10,11),Category_2:= "B",]

AIS[Cat1 %in% c(5,6,7,8),Category_1:= "C",]
AIS[Cat2 %in% c(5,6,7,8),Category_2:= "C",]

AIS[Cat1 == 13,Category_1:= "D",]
AIS[Cat2 == 13,Category_2:= "D",]

AIS[Cat1 == 16,Category_1:= "E",]
AIS[Cat2 == 16,Category_2:= "E",]

AIS[Cat1 == 17,Category_1:= "F",]
AIS[Cat2 == 17,Category_2:= "F",]

setkey(AIS,Category_1,Category_2)

unique(AIS[,Category_1,])
unique(AIS[,Category_2,])

# Mulighet 3b) Table triangle

table1 <- with(AIS, table(Category_1, Category_2) + table(Category_2, Category_1) )
data.table(table1)

# stargazer

{r, echo=TRUE, results='asis'}
stargazer(format(table1, quote=FALSE, justify="right"), type="latex")

### Mean Seconds in yield ###

AIS_mean1 <- data.table(AIS %>% group_by(Category_1,Category_2) %>%
                            summarize(Sec_spent1 = mean(Seconds_spent), Count1 = n()) )

AIS_mean2 <- data.table(AIS %>% group_by(Category_2,Category_1) %>%
                            summarize(Sec_spent2 = mean(Seconds_spent), Count2 = n()) )

AIS_mean <- AIS_mean1[,c("Category_1","Category_2"),]
AIS_mean <- AIS_mean2[,c("Category_1","Category_2"),]
AIS_mean
AIS_mean2
AIS_mean1
colnames(AIS_mean2) <- c("Category_1","Category_2","Sec_spent2","Count2")

AIS_mean <- merge(AIS_mean1,AIS_mean2,all.x = TRUE,all.y = TRUE)
setkey(AIS_mean,Category_2)
AIS_mean

AIS_mean <- AIS_mean %>% replace_na(list(Sec_spent1 = 0, Sec_spent2 = 0,Count1 = 0, Count2 = 0))

AIS_mean <- data.table(AIS_mean %>% group_by(Category_1,Category_2) %>% 
                           mutate(Sec_spent_mean = (Sec_spent1*Count1 + Sec_spent2*Count2)/(Count1+Count2)))

AIS_table <- dcast(AIS_mean,Category_1 ~ Category_2, value.var = "Sec_spent_mean")
t(as.matrix(AIS_table))


{r, echo=TRUE, results='asis'}
stargazer(format(AIS_table, quote=FALSE, justify="right"), type="latex")

### Mean Distance to CPA ###

AIS_mean1 <- data.table(AIS %>% group_by(Category_1,Category_2) %>%
                            summarize(mean_yield_dist_1 = mean((Distance_at_yield_1 + Distance_at_yield_2)/2), Count1 = n()) )

AIS_mean2 <- data.table(AIS %>% group_by(Category_2,Category_1) %>%
                            summarize(mean_yield_dist_2 = mean((Distance_at_yield_1 + Distance_at_yield_2)/2), Count2 = n()) )

colnames(AIS_mean2) <- c("Category_1","Category_2","mean_yield_dist_2","Count2")

AIS_mean <- merge(AIS_mean1,AIS_mean2,all.x = TRUE,all.y = TRUE)

AIS_mean <- AIS_mean %>% replace_na(list(mean_yield_dist_1 = 0, mean_yield_dist_2 = 0,Count1 = 0, Count2 = 0))

AIS_mean <- data.table(AIS_mean %>% group_by(Category_1,Category_2) %>% 
                           mutate(yield_dist_mean = (mean_yield_dist_1*Count1 + mean_yield_dist_2*Count2)/(Count1+Count2)))

AIS_mean

AIS_table <- dcast(AIS_mean,Category_1 ~ Category_2, value.var = "yield_dist_mean")

AIS_table


{r, echo=TRUE, results='asis'}
stargazer(format(AIS_table, quote=FALSE, justify="right"), type="latex")


### Mean Passing Distance ###

AIS_mean1 <- data.table(AIS %>% group_by(Category_1,Category_2) %>%
                            summarize(PassDist1 = mean(Passing_distance), Count1 = n()) )

AIS_mean2 <- data.table(AIS %>% group_by(Category_2,Category_1) %>%
                            summarize(PassDist2 = mean(Passing_distance), Count2 = n()) )

AIS_mean <- AIS_mean1[,c("Category_1","Category_2"),]
AIS_mean <- AIS_mean2[,c("Category_1","Category_2"),]

colnames(AIS_mean2) <- c("Category_1","Category_2","PassDist2","Count2")

AIS_mean <- merge(AIS_mean1,AIS_mean2,all.x = TRUE,all.y = TRUE)

AIS_mean

AIS_mean <- AIS_mean %>% replace_na(list(PassDist1 = 0, PassDist2 = 0,Count1 = 0, Count2 = 0))

AIS_mean <- data.table(AIS_mean %>% group_by(Category_1,Category_2) %>% 
                           mutate(PassDist = (PassDist1*Count1 + PassDist2*Count2)/(Count1+Count2)))

mean(c(AIS_mean$PassDist1,AIS_mean$PassDist2))

AIS_table <- dcast(AIS_mean,Category_1 ~ Category_2, value.var = "PassDist")

{r, echo=TRUE, results='asis'}
stargazer(format(AIS_table, quote=FALSE, justify="right"), type="latex")


### Mean Approach speed ###

AIS_mean1 <- data.table(AIS %>% group_by(Category_1,Category_2) %>%
                            summarize(App_speed1 = mean(Approach_speed), Count1 = n()) )

AIS_mean2 <- data.table(AIS %>% group_by(Category_2,Category_1) %>%
                            summarize(App_speed2 = mean(Approach_speed), Count2 = n()) )


colnames(AIS_mean2) <- c("Category_1","Category_2","App_speed2","Count2")

AIS_mean <- merge(AIS_mean1,AIS_mean2,all.x = TRUE,all.y = TRUE)

AIS_mean <- AIS_mean %>% replace_na(list(App_speed1 = 0, App_speed2 = 0,Count1 = 0, Count2 = 0))

AIS_mean <- data.table(AIS_mean %>% group_by(Category_1,Category_2) %>% 
                           mutate(App_speed_mean = (App_speed1*Count1 + App_speed2*Count2)/(Count1+Count2)))

AIS_mean

AIS_table <- dcast(AIS_mean,Category_1 ~ Category_2, value.var = "App_speed_mean")

AIS_table

{r, echo=TRUE, results='asis'}
stargazer(format(AIS_table, quote=FALSE, justify="right"), type="latex")


