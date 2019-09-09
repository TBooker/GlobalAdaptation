rm(list=ls())

library(ggplot2)
library(reshape2)

options(scipen=999)
cbbPalette <- c( "#0072B2", "#D55E00", "#CC79A7")

s02M1p001<-read.csv('~/work/FishersWaveFst/slim/N10000/analysisFiles/shat_0.02_M_1_pA_0.001.w10000.SEG.csv')
s02M1p0001<-read.csv('~/work/FishersWaveFst/slim/N10000/analysisFiles/shat_0.02_M_1_pA_0.0001.w10000.SEG.csv')
s02M10p001<-read.csv('~/work/FishersWaveFst/slim/N10000/analysisFiles/shat_0.02_M_10_pA_0.001.w10000.SEG.csv')
s02M10p0001<-read.csv('~/work/FishersWaveFst/slim/N10000/analysisFiles/shat_0.02_M_10_pA_0.0001.w10000.SEG.csv')

s01M1p001<-read.csv('~/work/FishersWaveFst/slim/N10000/analysisFiles/shat_0.01_M_1_pA_0.001.w10000.SEG.csv')
s01M1p0001<-read.csv('~/work/FishersWaveFst/slim/N10000/analysisFiles/shat_0.01_M_1_pA_0.0001.w10000.SEG.csv')
s01M10p001<-read.csv('~/work/FishersWaveFst/slim/N10000/analysisFiles/shat_0.01_M_10_pA_0.001.w10000.SEG.csv')
s01M10p0001<-read.csv('~/work/FishersWaveFst/slim/N10000/analysisFiles/shat_0.01_M_10_pA_0.0001.w10000.SEG.csv')

s005M1p001<-read.csv('~/work/FishersWaveFst/slim/N10000/analysisFiles/shat_0.005_M_1_pA_0.001.w10000.SEG.csv')
s005M1p0001<-read.csv('~/work/FishersWaveFst/slim/N10000/analysisFiles/shat_0.005_M_1_pA_0.0001.w10000.SEG.csv')
s005M10p001<-read.csv('~/work/FishersWaveFst/slim/N10000/analysisFiles/shat_0.005_M_10_pA_0.001.w10000.SEG.csv')
s005M10p0001<-read.csv('~/work/FishersWaveFst/slim/N10000/analysisFiles/shat_0.005_M_10_pA_0.0001.w10000.SEG.csv')

s001M1p001<-read.csv('~/work/FishersWaveFst/slim/N10000/analysisFiles/shat_0.001_M_1_pA_0.001.w10000.SEG.csv')
s001M1p0001<-read.csv('~/work/FishersWaveFst/slim/N10000/analysisFiles/shat_0.001_M_1_pA_0.0001.w10000.SEG.csv')
s001M10p001<-read.csv('~/work/FishersWaveFst/slim/N10000/analysisFiles/shat_0.001_M_10_pA_0.001.w10000.SEG.csv')
s001M10p0001<-read.csv('~/work/FishersWaveFst/slim/N10000/analysisFiles/shat_0.001_M_10_pA_0.0001.w10000.SEG.csv')


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


outty <- ggplot(dat = M, aes(x = as.factor(shat*20000), y = o99999, fill = as.factor(pA), col = as.factor(pA)))+
  geom_violin(alpha = 0.2, adjust = 2)+
  #stat_summary(fun.y = mean, aes(group = as.factor(pA), col = as.factor(pA)), geom = 'line',  position=position_dodge(0.9), alpha = 0.3)+
  stat_summary(fun.y = mean, aes(fill = as.factor(pA), col = as.factor(pA)), geom = 'point',shape = 3, size = 4, position=position_dodge(0.9))+
  #  scale_color_manual(alpha = 0.3)+
  facet_grid(~as.factor(M), labeller = label_parsed)+
  guides(fill = FALSE)+
  scale_color_manual(expression(italic(p[a])), values = cbbPalette)+
  scale_fill_manual(values = cbbPalette)+
  scale_y_sqrt("Proportion of Regions Containing Outliers",  expand = c(0.001, 0.001), limits = c(0,0.0051),breaks = c(0,0.0001,0.0005,0.001,0.002,0.003,0.004,0.005))+
  scale_x_discrete(expression(italic('2'*N[e] * hat(s[a]))))+
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

png('~/work/FishersWaveFst/Plots/Proportion_of_outliers.png', width = 650, height = 500)
print(outty)
dev.off()

outtyIRS <- ggplot(dat = M, aes(x = as.factor(shat*20000), y = o99999_irs, fill = as.factor(pA), col = as.factor(pA)))+
  geom_violin(alpha = 0.2, adjust = 2)+
  stat_summary(stat = mean, aes(fill = as.factor(pA), col = as.factor(pA)), geom = 'point',shape = 3, size = 4, position=position_dodge(0.9))+
  #  scale_color_manual(alpha = 0.3)+
  facet_grid(~as.factor(M), labeller = label_parsed)+
  guides(fill = FALSE)+
#  ggtitle('Windows linked to incomplete range-wide sweeps')+
  scale_color_manual(expression(italic(p[a])), values = cbbPalette)+
  scale_fill_manual(values = cbbPalette)+
  scale_y_sqrt("Proportion of Regions Containing Outliers",  expand = c(0.001, 0.001), limits = c(0,0.0051),breaks = c(0,0.0001,0.0005,0.001,0.002,0.003,0.004,0.005))+
  scale_x_discrete(expression(italic('2'*N[e] * hat(s[a]))))+
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

png('~/work/FishersWaveFst/Plots/Proportion_of_outliers_IRS.png', width = 650, height = 500)
print(outtyIRS)
dev.off()


outtyIRSratio <- ggplot(dat = M, aes(x = as.factor(shat*20000), y = o99999_irs/o99999, fill = as.factor(pA), col = as.factor(pA)))+
  geom_boxplot(alpha = 0.2, adjust = 2)+
  stat_summary(stat = mean, aes(fill = as.factor(pA), col = as.factor(pA)), geom = 'point',shape = 3, size = 4, position=position_dodge(0.9))+
  #  scale_color_manual(alpha = 0.3)+
  facet_grid(~as.factor(M), labeller = label_parsed)+
  guides(fill = FALSE)+
  scale_color_manual(expression(italic(p[a])), values = cbbPalette)+
  scale_fill_manual(values = cbbPalette)+
#  scale_y_sqrt("Proportion of Regions Containing Outliers",  expand = c(0.001, 0.001), limits = c(0,0.0051))+
  scale_x_discrete(expression(italic('2'*N[e] * hat(s[a]))))+
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
png('~/work/FishersWaveFst/Plots/Proportion_of_outliers_IRS_ratio.png', width = 650, height = 500)
print(outtyIRSratio)
dev.off()

outty2 <- ggplot(dat = M1, aes(x = as.factor(shat*20000), y = o99999, fill = as.factor(pA), col = as.factor(pA)))+
  geom_violin(alpha = 0.2, adjust = 2, width = .5, , position=position_dodge(0.5))+
  stat_summary(stat = mean, aes(fill = as.factor(pA), col = as.factor(pA)), geom = 'point',shape = 3, size = 4, position=position_dodge(0.5))+
  #  scale_color_manual(alpha = 0.3)+
  guides(fill = FALSE)+
  scale_color_manual(expression(italic(p[a])), values = cbbPalette)+
  scale_fill_manual(values = cbbPalette)+
  scale_y_sqrt("Proportion of Regions Containing Outliers",  expand = c(0.001, 0.001))+
  scale_x_discrete(expression(italic('2'*N[e] * hat(s[a]))))+
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

outty2

png('~/work/FishersWaveFst/Plots/Proportion_of_outliers_forTalk.png', width = 650, height = 500)
print(outty2)
dev.off()

col_names <- c("shat","w_size","M", "pA","rep","o999","o9999","o999_irs","o9999_irs")

M_melted <- melt(M,id = col_names)

M_melted$variable <- factor(M_melted$variable, levels = levels(M_melted$variable), labels = c( expression('All ' * 'Windows'), expression('Linked '*'to '*'Incomplete '*'Sweep') ))



incompleteComparison <- ggplot(dat = M_melted, aes(x = as.factor(shat*20000), y = value, fill = as.factor(pA), col = as.factor(pA)))+
  geom_violin(alpha = 0.2, adjust = 2)+
  stat_summary(stat = mean, aes(fill = as.factor(pA), col = as.factor(pA)), geom = 'point',shape = 3, size = 4, position=position_dodge(0.9))+
  #  scale_color_manual(alpha = 0.3)+
  facet_grid(as.factor(M)~variable, labeller = label_parsed)+
#  facet_grid(~as.factor(M), labeller = label_parsed)+
  guides(fill = FALSE)+
  scale_color_manual(expression(italic(p[a])), values = cbbPalette)+
  scale_fill_manual(values = cbbPalette)+
  scale_y_sqrt("Proportion of Regions Containing Outliers",  expand = c(0.001, 0.001), limits = c(0,0.0051),breaks = c(0,0.0001,0.0005,0.001,0.002,0.003,0.004,0.005))+
  scale_x_discrete(expression(italic('2'*N[e] * hat(s[a]))))+
  theme_bw()+
  theme(
    axis.text.y =  element_text(size = 15, color = 'black'),
    strip.text.y =  element_text(size = 15, color = 'black'),
    strip.text.x =  element_text(size = 12, color = 'black'),
    axis.text.x =  element_text(size = 15, color = 'black'),
    legend.title = element_text(size = 17, color = 'black'),
    legend.text = element_text(size = 15, color = 'black'),
    axis.title.x = element_text(size = 17, color = 'black'),
    axis.title.y = element_text(size = 17,  vjust = 1.5, color = 'black'),
    panel.grid = element_line(),         # All grid lines
    panel.grid.major = element_line(),   # Major grid lines
    panel.grid.minor = element_line(),   # Minor grid lines
    
    panel.grid.major.y = element_blank(), # Horizontal major grid lines
    panel.grid.minor.y = element_blank(),  # Vertical major grid lines
    panel.grid.major.x = element_blank(), # Horizontal major grid lines
    panel.grid.minor.x = element_blank()  # Vertical major grid lines
  )

png('~/work/FishersWaveFst/Plots/Proportion_of_outliers_incompleteComparison.png', width = 650, height = 500)
print(incompleteComparison)
dev.off()
