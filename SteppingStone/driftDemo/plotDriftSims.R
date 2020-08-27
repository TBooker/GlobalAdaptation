rm(list = ls())

x<-read.csv("driftedPop.csv", header = 0)
names(x) <- c("generation","p","population")
library(ggplot2)
png("quickDriftPlot.png", res = 300, units = "in", width = 6, height = 5)
ggplot(data = x, aes(x = generation, y = p, col = as.factor(population)))+
  geom_line()+
  theme_bw()+
  scale_x_continuous("Generation")+
  scale_y_continuous("Allele Frequency", limits = c(0,1))+
  scale_colour_brewer("Population", palette = "Set1")
dev.off()
