rm(list=ls())

library(ggplot2)
library(reshape2)

pal<-c( "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

### Neutral Sims
newt <-read.csv('/media/booker/mammal/work/FishersWaveFst/simulations/FstInNeutralSims2.csv',na.string = 'Na')

newtMelt <- melt(newt, id = c('meanP'))


newtMelt[ which(newtMelt$value == max(na.omit(newtMelt)$value)) ,]

ggplot(data = newtMelt, aes(x= value, fill = variable))+
  geom_density(alpha = 0.5)+
  scale_fill_manual(values = pal)+
#  facet_grid(~variable, scales  ='free_x')
#  scale_x_continuous(limits = c(0,1))+
  theme_bw()

### Selected sims

byGen_c0<-read.csv('~/work/FishersWaveFst/DFEsampling/N5K/summary.c0.wc.csv')
byGen_c0$c <- 0
byGen_c0.0001<-read.csv('~/work/FishersWaveFst/DFEsampling/N5K/summary.c1e-4.wc.csv')
byGen_c0.0001$c <- 0.0001

byGen <- rbind(byGen_c0, byGen_c0.0001)
byGen<- byGen[byGen$s < 0.02,]

byGen$s <- format(byGen$s, scientific = FALSE)

byGen$Ua <- factor(byGen$Ua, levels = levels(as.factor(byGen$Ua)), labels = c(expression(italic(U[a])*' = 0.1'), expression(italic(U[a])*' = 1'), expression(italic(U[a])*' = 10')))
byGen$distance <- factor(byGen$distance, levels = levels(as.factor(byGen$distance)), labels = c(expression(italic(d)*' = 100'),expression(italic(d)*' = 200'), expression(italic(d)*' = 300'), expression(italic(d)*' = 400')))

outputPlot<- ggplot(data =byGen , aes(x = s, y = fstOutliers, col = as.factor(c), fill = as.factor(c)))+
  geom_violin(alpha = 0.2, adjust = 2, width = 1)+
  stat_summary(stat = mean, aes(fill = as.factor(c), col = as.factor(c)), geom = 'point',shape = 3, size = 4, position=position_dodge(1))+
  facet_grid(Ua ~ distance, scales = 'free_y', labeller = label_parsed)+
  scale_color_manual('Recombination\nFraction',values = pal)+
  scale_fill_manual('Recombination\nFraction',values = pal)+
  scale_y_continuous(expression('Number of Outliers ('*italic(F[st])*' > 0.6)'))+
  scale_x_discrete(expression(italic(hat(s[a]))))+
  theme_bw()+
  theme(
    axis.text.y =  element_text(size = 12, color = 'black'),
    strip.text.x =  element_text(size = 14.5, color = 'black'),
    strip.text.y =  element_text(size = 14.5, color = 'black'),
    axis.text.x =  element_text(size = 12, color = 'black'),
    legend.title = element_text(size = 15, color = 'black'),
    legend.text = element_text(size = 12, color = 'black'),
    axis.title.x = element_text(size = 15, color = 'black'),
    axis.title.y = element_text(size = 15,  vjust = 1.5, color = 'black'),
    panel.grid = element_line(),         # All grid lines
  )


png('~/work/FishersWaveFst/Plots/number_of_outliers_steppingStone.png', width = 950, height = 450)
print(outputPlot)
dev.off()
