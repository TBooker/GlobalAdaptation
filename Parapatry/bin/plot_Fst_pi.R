rm(list=ls())
library(ggplot2)

pal<-c(  "#000000", "#E69F00", "#56B4E9")

glob <- read.csv('~/work/FishersWaveFst/slim/N10000/analysisFiles/shat_0.02_M_1_pA_0.0001.fstPi.csv')
glob$lab <- 'Global Adaptation'
glob$M <- 1
glob$shat <- 0.02
glob$generation <- 1
globTop <- glob[order(glob$WEIGHTED_FST, decreasing = T),]
globTops <- globTop[1:100,]
glob_sig <- glob[glob$WEIGHTED_FST > 0.46,]
glob_sig$sig <- 'sig'
glob_nonsig <- glob[glob$WEIGHTED_FST <= 0.46,]
glob_nonsig$sig <- 'nonsig'

neu <- read.csv('~/work/FishersWaveFst/slim/N10000/analysisFiles/pi_shat_0.0_M_1_pA_0.0_.csv')
neu$lab <- 'Neutral'
neu$M <- 1
neu$shat <- 0.02
neu$generation <- 1
neu_sig <- neu[neu$WEIGHTED_FST > 0.46,]
neu_sig$sig <- 'sig'
neu_nonsig <- neu[neu$WEIGHTED_FST <= 0.46,]
neu_nonsig$sig <- 'nonsig'

antag <- read.csv('~/work/FishersWaveFst/slim/N10000/analysisFiles/M1_s0.01.fstPi.csv')
antag$lab <- 'Local Adaptation'
local <- antag[antag$generation > 11000,]
l_sig <- local[local$WEIGHTED_FST > 0.46,]
l_sig$sig <- 'sig'
l_nonsig <- local[local$WEIGHTED_FST <= 0.46,]
l_nonsig$sig <- 'nonsig'
l_sig <- local[local$WEIGHTED_FST > 0.46,]
l_sig$sig <- 'sig'



pal <- setNames(pal, c('Neutral', 'Local Adaptation', 'Global Adaptation'))

ggplot(data = neu_nonsig, aes(x = all_pi/0.01, y = WEIGHTED_FST))+
  geom_point(color = '#000000', alpha = 0.2)+
  geom_point(data = neu_sig, aes(x = all_pi/0.01, y = WEIGHTED_FST),color = '#000000')+
  geom_point(data = l_nonsig, aes(x = all_pi/0.01, y = WEIGHTED_FST), color = '#E69F00', alpha = 0.5)+
  geom_point(data = l_sig, aes(x = all_pi/0.01, y = WEIGHTED_FST),color = '#E69F00')+
  geom_point(data = glob_nonsig, aes(x = all_pi/0.01, y = WEIGHTED_FST), color = '#56B4E9', alpha = 0.1)+
  geom_point(data = glob_sig, aes(x = all_pi/0.01, y = WEIGHTED_FST),color = '#56B4E9')+
  geom_hline(yintercept = 0.46, lty = 2) +
  scale_x_continuous(expression(pi/pi[0]), limits = c(0,2))+
  scale_y_continuous(expression(F[ST]), limits = c(0,1))+
  theme_bw()+
  theme(
    axis.text.y =  element_text(size = 15, color = 'black'),
    strip.text.x =  element_text(size = 17.5, color = 'black'),
    axis.text.x =  element_text(size = 15, color = 'black'),
    legend.title = element_blank(),
    legend.text = element_text(size = 15, color = 'black'),
    axis.title.x = element_text(size = 20, color = 'black'),
    axis.title.y = element_text(size = 20,  vjust = 1.5, color = 'black'),
    panel.grid = element_line(),         # All grid lines
    panel.grid.major = element_line(),   # Major grid lines
    panel.grid.minor = element_line(),   # Minor grid lines
    panel.grid.major.y = element_blank(), # Horizontal major grid lines
    panel.grid.minor.y = element_blank(),  # Vertical major grid lines
    panel.grid.major.x = element_blank(), # Horizontal major grid lines
    panel.grid.minor.x = element_blank()  # Vertical major grid lines
  )



glob_sig <- globTops[globTops$WEIGHTED_FST > 0.46,]
glob_sig$sig <- 'sig'
glob_nonsig <- globTops[globTops$WEIGHTED_FST <= 0.46,]
glob_nonsig$sig <- 'nonsig'

sig<- rbind(neu_sig, l_sig, glob_sig)
nonsig<- rbind(neu_nonsig, l_nonsig, glob_nonsig)

outty <- ggplot(data = nonsig, aes(x = all_pi/0.01, y = WEIGHTED_FST, col = lab))+
  geom_point(alpha = 0.5, size = 3)+
  geom_point(data =sig, aes(x = all_pi/0.01, y = WEIGHTED_FST, col = lab), size = 3)+
  geom_hline(yintercept = 0.46, lty = 2) +
  scale_x_continuous(expression(pi[T]/pi[0]), limits = c(0,2))+
  scale_y_continuous(expression(italic(F[ST])), limits = c(0,1))+
  scale_color_manual(values = pal)+
  theme_bw()+
  theme(
    axis.text.y =  element_text(size = 15, color = 'black'),
    strip.text.x =  element_text(size = 17.5, color = 'black'),
    axis.text.x =  element_text(size = 15, color = 'black'),
    legend.title = element_blank(),
    legend.text = element_text(size = 15, color = 'black'),
    axis.title.x = element_text(size = 20, color = 'black'),
    axis.title.y = element_text(size = 20,  vjust = 1.5, color = 'black'),
    panel.grid = element_line(),         # All grid lines
    panel.grid.major = element_line(),   # Major grid lines
    panel.grid.minor = element_line(),   # Minor grid lines
    panel.grid.major.y = element_blank(), # Horizontal major grid lines
    panel.grid.minor.y = element_blank(),  # Vertical major grid lines
    panel.grid.major.x = element_blank(), # Horizontal major grid lines
    panel.grid.minor.x = element_blank()  # Vertical major grid lines
  )


pdf('~/work/FishersWaveFst/Plots/Fst_V_diversity.pdf', width = 6.5, height = 5)
outty
dev.off()
?png
png('~/work/FishersWaveFst/Plots/Fst_V_diversity.png', width = 650, height = 500)
outty
dev.off()
tiff('~/work/FishersWaveFst/Plots/Fst_V_diversity.tiff', width = 650, height = 500)
outty
dev.off()
