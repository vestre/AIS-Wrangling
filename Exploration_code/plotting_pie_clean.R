### Plot pie chart ###

library(ggplot2)
library(scales)

antall <- c(262,651,481)

antallframe <- data.frame(Number=antall)
rownames(antallframe) <- c("Overtaking","Crossing","Head-on")
antallframe

blank_theme <- theme_minimal()+
    theme(
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.border = element_blank(),
        panel.grid=element_blank(),
        axis.ticks = element_blank(),
        plot.title=element_text(size=14, face="bold")
    )

plot2 <- ggplot(antallframe,aes(x="",y=Number,fill=rownames(antallframe))) + 
    geom_bar(width = 1, stat = "identity", alpha=0.4,colour="black") + 
    coord_polar("y", start=0) + 
    blank_theme +
    #scale_fill_brewer(palette = "Set1",direction = 1,name="COLREGs") +
    scale_fill_manual(values = c("yellow3","#E41A1C","#4DAF4A"), name="Type") + 
    geom_text(aes(y = Number/2 + c(0, cumsum(Number)[-length(Number)]), 
                  label = percent(Number/sum(Number))), size=10,fontface=1) + 
    theme(axis.text.x=element_blank(), plot.title = element_text(face="bold",size=25),
          legend.title = element_text(face="bold",size=20), legend.text = element_text(size=15)) + 
    ggtitle("Instances of situations by COLREGs type")
plot2


library("RColorBrewer")

brewer.pal(12, "Set1")
