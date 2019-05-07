### Plotting to map final ###

rm(list = ls())
setwd("/Users/arnsteinvestre/Desktop/Studier/19_vår/STK-MAT2011/3_Files/FINAL_NOR/analysis")

PATH = "~/code/DATA_DNV/ais_small.csv"

library(data.table)
library(fasttime)
library(tidyr)
library(zoo)
library(lubridate)
library(dplyr)
library(chron)

library(tidyverse)
library(mapdata)
library(maps)
library(stringr)
library(viridis)
library(ggplot2)
library(ggmap)

library(OpenStreetMap)
library(rJava)

library(cowplot)

AIS_NOR <- fread(PATH,
                 select=c("mmsi",
                          "date_time_utc",
                          "lon",
                          "lat",
                          "sog",
                          "cog",
                          "true_heading",
                          "nav_status",
                          "length",
                          "width",
                          "RISK_Norwegian_Main_Vessel_Category_ID",
                          "RISK_Norwegian_Main_Vessel_Category_Name")
                 ,check.names=TRUE
                 ,showProgress=TRUE
                 ,nThread=3
                 ,stringsAsFactors=TRUE
)



# RELOAD COPY #
AIS <- AIS_NOR


colnames(AIS)[colnames(AIS) == "date_time_utc"] <- "Time"
colnames(AIS)[colnames(AIS) == "mmsi"] <- "ID"
colnames(AIS)[colnames(AIS) == "nav_status"] <- "Status"
colnames(AIS)[colnames(AIS) == "lon"] <- "LON"
colnames(AIS)[colnames(AIS) == "lat"] <- "LAT"
colnames(AIS)[colnames(AIS) == "true_heading"] <- "Heading"
colnames(AIS)[colnames(AIS) == "sog"] <- "SOG"
colnames(AIS)[colnames(AIS) == "cog"] <- "COG"
colnames(AIS)[colnames(AIS) == "length"] <- "Length"
colnames(AIS)[colnames(AIS) == "width"] <- "Width"
colnames(AIS)[colnames(AIS) == "RISK_Norwegian_Main_Vessel_Category_ID"] <- "CatID"
colnames(AIS)[colnames(AIS) == "RISK_Norwegian_Main_Vessel_Category_Name"] <- "CatName"



# FORMATTING TIME #
AIS[,Time := fastPOSIXct(AIS[,Time,]),]
AIS[,TimeAsInt := as.integer(AIS[,Time,]),]
typeof(AIS[,TimeAsInt,][1])

#### Plotting
AIS <- AIS[,Interval := (TimeAsInt - min(TimeAsInt)) %/% 20 ,]

LAT1C =  55 ; LAT2C = 75
LON1C = 0 ; LON2C = 32

DMC = 1.1132e+5

map <- openmap(c(LAT2C,LON1C), c(LAT1C,LON2C), zoom = NULL,
               type = c("osm", "stamen-toner", "stamen-terrain","stamen-watercolor", "esri","esri-topo")[3],
               mergeTiles = TRUE) 

### If plot in latlon projection ###
## OSM CRS :: "+proj=merc +a=6378137 +b=6378137 +lat_ts=0.0 +lon_0=0.0 +x_0=0.0 +y_0=0 +k=1.0 +units=m +nadgrids=@null +no_defs"

map.latlon <- openproj(map, projection = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
#map.latlon <- openproj(map, projection = "+proj=merc +ellps=WGS84 +datum=WGS84 +no_defs +b=6378137 +lon_0=0.0 +x_0=0.0")

mytheme <- theme(
                #plot.title = element_text(face = "bold",size = 15),
                 plot.title = element_blank(),
                 panel.background = element_rect(colour = NA),
                 plot.background = element_rect(colour = NA),
                 axis.title = element_text(size = 10),
                 axis.title.y = element_text(angle=90,vjust =2),
                 axis.title.x = element_text(vjust = -0.2),
                 panel.border = element_rect(colour = "black", fill=NA, size=5),
                 axis.line = element_line(colour="grey",size=0.5),
                 axis.line.x.top = element_line(colour="grey",size=0.5),
                 axis.line.y.right =element_line(colour="grey",size=0.5))

# Only  map latlon
mapplot1 <- autoplot(map.latlon) + mytheme + 
    labs(title = "Ships in Norway, Interval 1",x = "LAT", y="LON")

mapplot1

# Only map mercator
mapplot2 <- autoplot(map) + mytheme + 
    labs(title = "Ships in Norway, Interval 1",x = "LAT", y="LON")

mapplot2

# Filters for data
filter1 = (AIS[,Interval] == 1500)&(AIS[,Length] != 0)&(AIS[,Width] != 0)&(!is.na(AIS[,Width]))&(!is.na(AIS[,Length]))
filter2 <- (AIS[,LAT] > LAT1C & AIS[,LAT] < LAT2C)&(AIS[,LON] > LON1C & AIS[,LON] < LON2C)

# Creating data table
AISforplot <- AIS[filter1&filter2,,]

maxWid <- max(AISforplot[,Width,])
maxLength <- max(AISforplot[,Length,])

## Plotting vessels in latlon ##
matplot3 <- autoplot(map.latlon) +
    labs(title = "Vessels in Norway",x = "LON", y="LAT") + mytheme +
    geom_point(data=AISforplot,
               aes(x=LON,y=LAT)) +
    geom_segment(data=AISforplot,
                 aes(xend=(LON),
                     yend=(LAT),
                     x=(LON) + 3000*SOG*1852*sin(COG*3.14/180)/(3600*DMC),
                     y=(LAT) + 3000*SOG*1852*cos(COG*3.14/180)/(3600*DMC))
                 ,col="red",alpha=0.5) + 
    guides(size=FALSE)
matplot3

## Plotting vessels in mercator ##

filter1 = (AIS[,Interval] == 1500)&(AIS[,Length] != 0)&(!is.na(AIS[,Length])&(AIS[,SOG,]<50))
filter2 <- (AIS[,LAT] > LAT1C & AIS[,LAT] < LAT2C)&(AIS[,LON] > LON1C & AIS[,LON] < LON2C)

AISforplot <- AIS[filter1&filter2,,]

AISforplot[,x:=NULL,]
AISforplot[,y:=NULL,]
AISforplot

# Converting latlons to mercator coordinatesd WGS
points2 <- as.data.frame( 
    OpenStreetMap::projectMercator( lat = AISforplot[,LAT,], 
                                    long = AISforplot[,LON,] ) 
)
points2
AISforplot <- cbind(AISforplot,points2)

# Final plot of vessels in mercator
matplot4 <- autoplot(map) +
    labs(title = "Vessels in Norway",x = "LON", y="LAT") + mytheme +
    geom_point(data=AISforplot,
               aes(x=x,y=y)) +
    geom_segment(data=AISforplot,
                 aes(xend=(x),
                     yend=(y),
                     x=(x) + 300000000*SOG*1852*sin(COG*3.14/180)/(3600*DMC),
                     y=(y) + 300000000*SOG*1852*cos(COG*3.14/180)/(3600*DMC))
                 ,col="red",alpha=0.5) + 
    guides(size=FALSE)
matplot4

# Same plot w/o scale

matplot5 <- autoplot(map) +
    labs(title = "Vessels in Norway",x = NULL, y=NULL) + mytheme +
    geom_point(data=AISforplot,
               aes(x=x,y=y)) +
    geom_segment(data=AISforplot,
                 aes(xend=(x),
                     yend=(y),
                     x=(x) + 300000000*SOG*1852*sin(COG*3.14/180)/(3600*DMC),
                     y=(y) + 300000000*SOG*1852*cos(COG*3.14/180)/(3600*DMC))
                 ,col="red",alpha=0.5) + 
    scale_x_continuous(breaks = NULL) + 
    scale_y_continuous(breaks = NULL) +
    guides(size=FALSE)
matplot5

# Same plot w/ border and scales (final)

matplot6 <- autoplot(map) +
    labs(title = NULL,x = "LON", y="LAT") +
    geom_point(data=AISforplot,
               aes(x=x,y=y)) +
    geom_segment(data=AISforplot,
                 aes(xend=(x),
                     yend=(y),
                     x=(x) + 300000000*SOG*1852*sin(COG*3.14/180)/(3600*DMC),
                     y=(y) + 300000000*SOG*1852*cos(COG*3.14/180)/(3600*DMC))
                 ,col="red",alpha=0.5) + 
    #scale_x_continuous(breaks = waiver()) + 
    #scale_y_continuous(breaks = waiver()) +
    scale_x_continuous(breaks = round(seq(min(AISforplot$x),max(AISforplot$x),length.out=5),1), labels = round(seq(LON1C,LON2C,length.out=5),2) ,name="LON") + 
    scale_y_continuous(breaks = round(seq(min(AISforplot$y),max(AISforplot$y),length.out=6),1), labels = round(seq(LAT1C,LAT2C,length.out=6),2) ,name="LAT") +
    guides(size=FALSE) + 
    theme_classic() +
    mytheme
matplot6
