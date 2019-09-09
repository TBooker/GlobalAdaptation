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

#M1$M <- factor(M1$M,levels = c(1), labels = c(expression(italic(N[e]*'m = 1'))))
#M10$M <- factor(M10$M,levels = c(10), labels = c(expression(italic(N[e]*'m = 10'))))
M_snps <- rbind(M1,M10)
M_snps$source <- 'SNPs'




ws02M1p001<-read.csv('~/work/FishersWaveFst/slim/N10000/analysisFiles/shat_0.02_M_1_pA_0.001.w10000.SEG.csv')
ws02M1p0001<-read.csv('~/work/FishersWaveFst/slim/N10000/analysisFiles/shat_0.02_M_1_pA_0.0001.w10000.SEG.csv')
ws02M10p001<-read.csv('~/work/FishersWaveFst/slim/N10000/analysisFiles/shat_0.02_M_10_pA_0.001.w10000.SEG.csv')
ws02M10p0001<-read.csv('~/work/FishersWaveFst/slim/N10000/analysisFiles/shat_0.02_M_10_pA_0.0001.w10000.SEG.csv')

ws01M1p001<-read.csv('~/work/FishersWaveFst/slim/N10000/analysisFiles/shat_0.01_M_1_pA_0.001.w10000.SEG.csv')
ws01M1p0001<-read.csv('~/work/FishersWaveFst/slim/N10000/analysisFiles/shat_0.01_M_1_pA_0.0001.w10000.SEG.csv')
ws01M10p001<-read.csv('~/work/FishersWaveFst/slim/N10000/analysisFiles/shat_0.01_M_10_pA_0.001.w10000.SEG.csv')
ws01M10p0001<-read.csv('~/work/FishersWaveFst/slim/N10000/analysisFiles/shat_0.01_M_10_pA_0.0001.w10000.SEG.csv')

ws005M1p001<-read.csv('~/work/FishersWaveFst/slim/N10000/analysisFiles/shat_0.005_M_1_pA_0.001.w10000.SEG.csv')
ws005M1p0001<-read.csv('~/work/FishersWaveFst/slim/N10000/analysisFiles/shat_0.005_M_1_pA_0.0001.w10000.SEG.csv')
ws005M10p001<-read.csv('~/work/FishersWaveFst/slim/N10000/analysisFiles/shat_0.005_M_10_pA_0.001.w10000.SEG.csv')
ws005M10p0001<-read.csv('~/work/FishersWaveFst/slim/N10000/analysisFiles/shat_0.005_M_10_pA_0.0001.w10000.SEG.csv')

ws001M1p001<-read.csv('~/work/FishersWaveFst/slim/N10000/analysisFiles/shat_0.001_M_1_pA_0.001.w10000.SEG.csv')
ws001M1p0001<-read.csv('~/work/FishersWaveFst/slim/N10000/analysisFiles/shat_0.001_M_1_pA_0.0001.w10000.SEG.csv')
ws001M10p001<-read.csv('~/work/FishersWaveFst/slim/N10000/analysisFiles/shat_0.001_M_10_pA_0.001.w10000.SEG.csv')
ws001M10p0001<-read.csv('~/work/FishersWaveFst/slim/N10000/analysisFiles/shat_0.001_M_10_pA_0.0001.w10000.SEG.csv')


wM1 <-rbind(ws02M1p001,
           ws02M1p0001,
           ws01M1p001,
           ws01M1p0001,
           ws005M1p001,
           ws005M1p0001,
           ws001M1p001,
           ws001M1p0001
)


wM10 <-rbind(ws02M10p001,
            ws02M10p0001,
            ws01M10p001,
            ws01M10p0001,
            ws005M10p001,
            ws005M10p0001,
            ws001M10p001,
            ws001M10p0001)

#M1$M <- factor(M1$M,levels = c(1), labels = c(expression(italic(N[e]*'m = 1'))))
#M10$M <- factor(M10$M,levels = c(10), labels = c(expression(italic(N[e]*'m = 10'))))
M_windows <- rbind(wM1,wM10)
M_windows$source <- 'windows'
M_windows$prop <- M_windows$o99999

keeps <- c( "shat", "pA","M","prop","source")  

M_windows <- subset(M_windows, select = keeps)
M_snps <- subset(M_snps, select = keeps)

M <-rbind( M_snps, M_windows)
M$source<- levels(as.factor(M$source))
M$M <- as.factor(M$M)


M2<- rbind(M_windows, M_snps)   


M2$M <- factor(M2$M,levels = c(1, 10), labels = c(expression(italic(N[e]*'m = 1')), expression(italic(N[e]*'m = 10'))))
M2$source <- factor(M2$source,levels = c('SNPs', 'windows'), labels = c(expression('SNP-based'), expression('Window-based')))

outty <- ggplot(dat = M2, aes(x = as.factor(shat*20000),y = prop, fill = as.factor(pA)))+
  geom_violin(alpha = 0.2, adjust = 2)+
  stat_summary(fun.y = mean, aes(fill = as.factor(pA), col = as.factor(pA)), geom = 'point',shape = 3, size = 4, position=position_dodge(0.9))+
  facet_grid(as.factor(M)~source, labeller = label_parsed)+
  scale_y_sqrt("Proportion of regions\ncontaining outliers/top-candidates", breaks = c(0,0.0001,0.0005,0.001,0.002,0.003,0.004,0.005))+
  scale_x_discrete(expression(italic('2'*N[e] * bar(s[a]))))+
  scale_color_manual(expression(italic(p[a])), values = cbbPalette)+
  scale_fill_manual(expression(italic(p[a])),values = cbbPalette)+
  theme_bw()+
  theme(
    axis.text.y =  element_text(size = 15, color = 'black'),
    strip.text.x =  element_text(size = 17.5, color = 'black'),
    strip.text.y =  element_text(size = 17.5, color = 'black'),
    axis.text.x =  element_text(size = 15, color = 'black'),
    legend.title = element_text(size = 20, color = 'black'),
    legend.text = element_text(size = 15, color = 'black'),
    axis.title.x = element_text(size = 20, color = 'black'),
    axis.title.y = element_text(size = 20,  vjust = 1.5, color = 'black'),
    panel.grid = element_line(),         # All grid lines
    panel.grid.major = element_line(),   # Major grid lines
    panel.grid.minor = element_line(),   # Minor grid lines
    
  #  panel.grid.major.y = element_blank(), # Horizontal major grid lines
   # panel.grid.minor.y = element_blank(),  # Vertical major grid lines
    panel.grid.major.x = element_blank(), # Horizontal major grid lines
    panel.grid.minor.x = element_blank()  # Vertical major grid lines
  )



png('~/work/FishersWaveFst/Plots/Proportion_of_both.png', width = 650, height = 650)
print(outty)
dev.off()



