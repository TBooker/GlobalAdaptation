rm(list=ls())

library(ggplot2)
options(scipen=999)
cbbPalette <- c( "#0072B2", "#D55E00", "#CC79A7")

s02M1p001<-read.csv('~/work/FishersWaveFst/slim/N10000/test/SNPs_shat_0.02_M_1_pA_0.001.summary.csv')
s02M1p0001<-read.csv('~/work/FishersWaveFst/slim/N10000/test/SNPs_shat_0.02_M_1_pA_0.0001.summary.csv')
s02M10p001<-read.csv('~/work/FishersWaveFst/slim/N10000/test/SNPs_shat_0.02_M_10_pA_0.001.summary.csv')
s02M10p0001<-read.csv('~/work/FishersWaveFst/slim/N10000/test/SNPs_shat_0.02_M_10_pA_0.0001.summary.csv')

s01M1p001<-read.csv('~/work/FishersWaveFst/slim/N10000/test/SNPs_shat_0.01_M_1_pA_0.001.summary.csv')
s01M1p0001<-read.csv('~/work/FishersWaveFst/slim/N10000/test/SNPs_shat_0.01_M_1_pA_0.0001.summary.csv')
s01M10p001<-read.csv('~/work/FishersWaveFst/slim/N10000/test/SNPs_shat_0.01_M_10_pA_0.001.summary.csv')
s01M10p0001<-read.csv('~/work/FishersWaveFst/slim/N10000/test/SNPs_shat_0.01_M_10_pA_0.0001.summary.csv')

s005M1p001<-read.csv('~/work/FishersWaveFst/slim/N10000/test/SNPs_shat_0.005_M_1_pA_0.001.summary.csv')
s005M1p0001<-read.csv('~/work/FishersWaveFst/slim/N10000/test/SNPs_shat_0.005_M_1_pA_0.0001.summary.csv')
s005M10p001<-read.csv('~/work/FishersWaveFst/slim/N10000/test/SNPs_shat_0.005_M_10_pA_0.001.summary.csv')
s005M10p0001<-read.csv('~/work/FishersWaveFst/slim/N10000/test/SNPs_shat_0.005_M_10_pA_0.0001.summary.csv')

s001M1p001<-read.csv('~/work/FishersWaveFst/slim/N10000/test/SNPs_shat_0.001_M_1_pA_0.001.summary.csv')
s001M1p0001<-read.csv('~/work/FishersWaveFst/slim/N10000/test/SNPs_shat_0.001_M_1_pA_0.0001.summary.csv')
s001M10p001<-read.csv('~/work/FishersWaveFst/slim/N10000/test/SNPs_shat_0.001_M_10_pA_0.001.summary.csv')
s001M10p0001<-read.csv('~/work/FishersWaveFst/slim/N10000/test/SNPs_shat_0.001_M_10_pA_0.0001.summary.csv')


M1 <-rbind(s02M1p001,
           s02M1p0001,
           s01M1p001,
           s01M1p0001,
           s005M1p001,
           s005M1p0001,
           s001M1p001,
           s001M1p0001
)


M10 <-rbind(s02M10p001,
            s02M10p0001,
            s01M10p001,
            s01M10p0001,
            s005M10p001,
            s005M10p0001,
            s001M10p001,
            s001M10p0001)

M1$M <- factor(M1$M,levels = c(1), labels = c(expression(italic(N[e]*'m = 1'))))
M10$M <- factor(M10$M,levels = c(10), labels = c(expression(italic(N[e]*'m = 10'))))
M <- rbind(M1,M10)

M







outty <- ggplot(dat = M, aes(x = as.factor(shat*20000), y = prop, fill = as.factor(pA), col = as.factor(pA)))+
  geom_violin(alpha = 0.2, adjust = 2)+
  stat_summary(fun.y = mean, aes(fill = as.factor(pA), col = as.factor(pA)), geom = 'point',shape = 3, size = 4, position=position_dodge(0.9))+
  #  scale_color_manual(alpha = 0.3)+
  facet_grid(~as.factor(M), labeller = label_parsed)+
  guides(fill = FALSE)+
  scale_color_manual(expression(italic(p[a])), values = cbbPalette)+
  scale_fill_manual(values = cbbPalette)+
  scale_y_sqrt("Proportion of Regions\nContaining Top-Candidates",  expand = c(0.001, 0.001), limits = c(0,0.001))+
  scale_x_discrete(expression(italic('2'*N[e] * bar(s[a]))))+
  theme_bw()+
  theme(
    axis.text.y =  element_text(size = 15, color = 'black'),
    strip.text.x =  element_text(size = 17.5, color = 'black'),
    axis.text.x =  element_text(size = 15, color = 'black'),
    legend.title = element_text(size = 20, color = 'black'),
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

outty


png('~/work/FishersWaveFst/Plots/Proportion_of_top_candidates.png', width = 650, height = 300)
print(outty)
dev.off()
