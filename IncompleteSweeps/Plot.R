rm(list = ls())
library(ggplot2) 
library(ggpubr)
library(wesanderson)
cols <- c( wes_palette("Rushmore1", 5)[3], wes_palette("Rushmore1", 5)[5])


options(scipen = 99)
TommyTheme <-theme_bw()+
  theme(
    axis.title.x = element_text(size=14,angle=0, color = 'black'),
    axis.title.y = element_text(size=20,angle=90,vjust=0.5, color = 'black'),
    axis.text.x = element_text(size=12,angle=0, color = 'black'),
    axis.text.y = element_text(size=12,angle=0, color = 'black'),
    legend.title = element_text(size =15),
    legend.text = element_text(size =13, face = 'italic'),
    strip.text.y = element_text(size = 15)
  )

base_breaks <- function(n = 10){
  function(x) {
    axisTicks(log10(range(x, na.rm = TRUE)), log = TRUE, n = n)
  }
}

pal <- c( "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

### A small script to plot the output of model for the number of intermediate frequnecy sweeps present in any popuation at a given point in time

num0.01 <- read.csv('~/UBC/FishersWaveFst/MKultra/Number_of_sweeps_pA0.01.csv', header = F)
colnames( num0.01 ) <- c('Nes', 'Num') 
num0.01$pA <- 0.01

num0.001 <- read.csv('~/UBC/FishersWaveFst/MKultra/Number_of_sweeps_pA0.001.csv')
colnames( num0.001 ) <- c('Nes', 'Num')
num0.001$pA <- 0.001

num0.0001 <- read.csv('~/UBC/FishersWaveFst/MKultra/Number_of_sweeps_pA0.0001.csv')
num0.0001
colnames( num0.0001 ) <- c('Nes', 'Num')
num0.0001$pA <- 0.0001

num <- rbind(num0.001, num0.0001)
num$pA<- factor(num$pA, levels = c(0.001, 0.0001), labels = c(expression(italic(p[a])*' = 0.001'), expression(italic(p[a])*' = 0.0001')))

temp <- list.files(path='~/UBC/FishersWaveFst/MKultra/',pattern="segSites.csv", full.names= T)
myfiles = lapply(temp, read.csv)
segSites = as.data.frame(do.call(rbind, myfiles))
segSites <- segSites[segSites$M == 1,]
segSites$pA <- as.factor(segSites$pA)
segSites$pA<- factor(segSites$pA, levels = c(0.001, 0.0001), labels = c(expression(italic(p[a])*' = 0.001'), expression(italic(p[a])*' = 0.0001')))



Nump<-ggplot(data  = num , aes(x = Nes*2, y = Num, col = pA))+
  geom_line(lwd = 2, alpha  = 0.6)+
  stat_summary(data  = segSites, geom = 'point', fun.y = mean, aes(x = Nes*2, y = sweeps, group = pA), shape = 4, size = 3)+
#  stat_summary(fun.y=mean,geom="line",lwd=2,aes(group=1))
  
  scale_y_log10('Number of\nIncomplete\nSweeps', limits = c(20,1000), breaks = c(20,50,100,200,500,1000))+
  scale_x_log10(expression(italic('2'*N[e]*bar(s[a]))), breaks = c(10,20,50,100,200,500), limits = c(4,500))+
  scale_color_manual('', values = pal, labels = scales::parse_format())+
  guides(col = guide_legend(reverse = TRUE), shape = F)+
  theme_bw()+
  theme(axis.title.y = element_text(angle = 90, vjust = 0.5))+
  TommyTheme

print(Nump)

MK0.01 <- read.csv('~/UBC/FishersWaveFst/MKultra/MKalphaInt_pA0.01.csv', header = F)
colnames( MK0.01 ) <- c('Nes', 'alpha')
MK0.01$pA <- 0.01

MK0.001 <- read.csv('~/UBC/FishersWaveFst/MKultra/MKalphaInt_pA0.001.csv')

colnames( MK0.001 ) <- c('Nes', 'alpha')


MK0.001$pA <- 0.001

MK0.0001 <- read.csv('~/UBC/FishersWaveFst/MKultra/MKalphaInt_pA0.0001.csv')
colnames( MK0.0001 ) <- c('Nes', 'alpha')
MK0.0001$pA <- 0.0001

MK <- rbind(MK0.001, MK0.0001)
MK$pA<- as.factor(MK$pA)
MK$pA<- factor(MK$pA, levels = c(0.001, 0.0001), labels = c(expression(italic(p[a])*' = 0.001'), expression(italic(p[a])*' = 0.0001')))


MKp<- ggplot(data  = MK , aes(x = Nes*2, y = alpha, col = pA))+
  geom_line(lwd = 1.5)+
  scale_y_continuous(expression(alpha), breaks = c(0, 0.25, 0.5, 0.75, 1.0), limits = c(0,1))+
  scale_x_log10(expression(italic('2'*N[e]*bar(s[a]))), breaks = c(10,20,50,100,200,500), limits = c(4,500))+
  scale_color_manual(expression(italic(p[a])), values = pal)+
  guides(col = guide_legend(reverse = TRUE), shape = F)+
  theme_bw()+
  theme(axis.title.y = element_text(angle = 0, vjust = 0.5))+
  TommyTheme



zz<- ggarrange(Nump, MKp, nrow = 2, ncol = 1, align = 'v', common.legend = T, legend = 'right', labels = 'AUTO')
zz
png('~/UBC/FishersWaveFst/MKultra/NumSweeps_v_alpha.png', height = 400, width = 550)
print(zz)
dev.off()



