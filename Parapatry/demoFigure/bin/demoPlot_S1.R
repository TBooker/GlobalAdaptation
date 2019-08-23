rm(list = ls() )

library(ggplot2)
library(ggpubr)


pal <- c( "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

trace <- read.csv('~/UBC/FishersWaveFst/DemoFigure/freqSlice.M1_p0.0001_4.csv')

#sub <- xx[(xx$WEIGHTED_FST>0.65)&(xx$gen<5000),]
#sub$id <- factor(sub$id)
yy<- trace#[xx$id %in% sub$id,]

yy[(yy$gen < yy$firstGen)&(is.na(yy$q) ) ,]$q <- 0
yy[(yy$gen < yy$firstGen)&(is.na(yy$q1) ) ,]$q1 <- 0
yy[(yy$gen < yy$firstGen)&(is.na(yy$q2) ) ,]$q2 <- 0

yy[(yy$gen > yy$firstGen)&(is.na(yy$q) ) ,]$q <- 1
yy[(yy$gen > yy$firstGen)&(is.na(yy$q1) ) ,]$q1 <- 1
yy[(yy$gen > yy$firstGen)&(is.na(yy$q2) ) ,]$q2 <- 1

slice <- yy[yy$BIN_START <2130001,]
slice<- slice[order(slice$REP_x),]
slice$id <- paste(slice$REP_x,slice$BIN_START, sep = '_')
slice$POS <- slice$BIN_START + (slice$BIN_END-slice$BIN_START)/2
slice$CUM_POS <- slice$POS + (slice$REP_x * 2000000) - 2000000 


manhattan1 <- read.csv('~/UBC/FishersWaveFst/DemoFigure/M1_p0.0001_4_files/gen5000.M1_p0.0001_4.csv.gz')
manhattan1 <- manhattan1[manhattan1$BIN_START <2130001,]
manhattan1<- manhattan1[order(manhattan1$REP),]
manhattan1$id <- paste(manhattan1$REP,manhattan1$BIN_START)
manhattan1$POS <- manhattan1$BIN_START + (manhattan1$BIN_END-manhattan1$BIN_START)/2
manhattan1$CUM_POS <- manhattan1$POS + (manhattan1$REP * 2000000) - 2000000 
slice1 <- slice[slice$gen == 5000,]

manhattan2 <- read.csv('~/UBC/FishersWaveFst/DemoFigure/M1_p0.0001_4_files/gen5500.M1_p0.0001_4.csv.gz')
manhattan2 <- manhattan2[manhattan2$BIN_START <2130001,]
manhattan2<- manhattan2[order(manhattan2$REP),]
manhattan2$id <- paste(manhattan2$REP,manhattan2$BIN_START)
manhattan2$POS <- manhattan2$BIN_START + (manhattan2$BIN_END-manhattan2$BIN_START)/2
manhattan2$CUM_POS <- manhattan2$POS + (manhattan2$REP * 2000000) - 2000000 
slice2 <- slice[slice$gen == 5500,]


manhattan3 <- read.csv('~/UBC/FishersWaveFst/DemoFigure/M1_p0.0001_4_files/gen6000.M1_p0.0001_4.csv.gz')
manhattan3 <- manhattan3[manhattan3$BIN_START <2130001,]
manhattan3<- manhattan3[order(manhattan3$REP),]
manhattan3$id <- paste(manhattan3$REP,manhattan3$BIN_START)
manhattan3$POS <- manhattan3$BIN_START + (manhattan3$BIN_END-manhattan3$BIN_START)/2
manhattan3$CUM_POS <- manhattan3$POS + (manhattan3$REP * 2000000) - 2000000 
slice3 <- slice[slice$gen == 6000,]


manhattan4 <- read.csv('~/UBC/FishersWaveFst/DemoFigure/M1_p0.0001_4_files/gen6500.M1_p0.0001_4.csv.gz')
manhattan4 <- manhattan4[manhattan4$BIN_START <2130001,]
manhattan4<- manhattan4[order(manhattan4$REP),]
manhattan4$id <- paste(manhattan4$REP,manhattan4$BIN_START)
manhattan4$POS <- manhattan4$BIN_START + (manhattan4$BIN_END-manhattan4$BIN_START)/2
manhattan4$CUM_POS <- manhattan4$POS + (manhattan4$REP * 2000000) - 2000000 
slice4 <- slice[slice$gen == 6500,]

manhattan5 <- read.csv('~/UBC/FishersWaveFst/DemoFigure/M1_p0.0001_4_files/gen7000.M1_p0.0001_4.csv.gz')
manhattan5 <- manhattan5[manhattan5$BIN_START <2130001,]
manhattan5<- manhattan5[order(manhattan5$REP),]
manhattan5$id <- paste(manhattan5$REP,manhattan5$BIN_START)
manhattan5$POS <- manhattan5$BIN_START + (manhattan5$BIN_END-manhattan5$BIN_START)/2
manhattan5$CUM_POS <- manhattan5$POS + (manhattan5$REP * 2000000) - 2000000 
slice5 <- slice[slice$gen == 7000,]


sel1plot<- ggplot(data = manhattan1, aes(x = CUM_POS/1e6, y = WEIGHTED_FST))+
  geom_point(aes(colour = cut(WEIGHTED_FST, c(-Inf, 0.46, Inf)), size = cut(WEIGHTED_FST, c(-Inf, 0.46, Inf)))) +
  scale_x_continuous('Position in Genome', limits = c(0, max(manhattan1$CUM_POS)/1e6))+
  scale_y_continuous(expression(italic(F[ST])), limits = c(0,0.75))+
  geom_hline(aes(yintercept = 0.46), lty = 2, alpha = 0.5)+
  scale_color_manual(values = c("black","red"))+
  scale_size_manual(values = c(0.7,2.7))+
  ggtitle('Generation 5000')+
  guides(colour = F, size = F)+
  theme(
    axis.ticks.y  = element_blank(),
    axis.ticks.x  = element_blank(),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size=13,colour = "black"),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size=15)
  )

sel2plot<- ggplot(data = manhattan2, aes(x = CUM_POS/1e6, y = WEIGHTED_FST))+
  geom_point(aes(colour = cut(WEIGHTED_FST, c(-Inf, 0.46, Inf)), size = cut(WEIGHTED_FST, c(-Inf, 0.46, Inf)))) +
  scale_x_continuous('Position in Genome', limits = c(0, max(manhattan2$CUM_POS)/1e6))+
  scale_y_continuous(expression(italic(F[ST])), limits = c(0,0.75))+
  geom_hline(aes(yintercept = 0.46), lty = 2, alpha = 0.5)+
  scale_color_manual(values = c("black","red"))+
  scale_size_manual(values = c(0.7,2.7))+
  ggtitle('Generation 5500')+
  guides(colour = F, size = F)+
  theme(
    axis.ticks.y  = element_blank(),
    axis.ticks.x  = element_blank(),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size=13,colour = "black"),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size=15)
  )

sel3plot<- ggplot(data = manhattan3, aes(x = CUM_POS/1e6, y = WEIGHTED_FST))+
  geom_point(aes(colour = cut(WEIGHTED_FST, c(-Inf, 0.46, Inf)), size = cut(WEIGHTED_FST, c(-Inf, 0.46, Inf)))) +
  scale_x_continuous('Position in Genome', limits = c(0, max(manhattan3$CUM_POS)/1e6))+
  scale_y_continuous(expression(italic(F[ST])), limits = c(0,0.75))+
  geom_hline(aes(yintercept = 0.46), lty = 2, alpha = 0.5)+
  scale_color_manual(values = c("black","red"))+
  scale_size_manual(values = c(0.7,2.7))+
  ggtitle('Generation 6000')+
  guides(colour = F, size = F)+
  theme(
    axis.ticks.y  = element_blank(),
    axis.ticks.x  = element_blank(),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size=13,colour = "black"),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size=15)
  )

sel4plot<- ggplot(data = manhattan4, aes(x = CUM_POS/1e6, y = WEIGHTED_FST))+
  geom_point(aes(colour = cut(WEIGHTED_FST, c(-Inf, 0.46, Inf)), size = cut(WEIGHTED_FST, c(-Inf, 0.46, Inf)))) +
  scale_x_continuous('Position in Genome', limits = c(0, max(manhattan4$CUM_POS)/1e6))+
  scale_y_continuous(expression(italic(F[ST])), limits = c(0,0.75))+
  geom_hline(aes(yintercept = 0.46), lty = 2, alpha = 0.5)+
  scale_color_manual(values = c("black","red"))+
  scale_size_manual(values = c(0.7,2.7))+
  ggtitle('Generation 6500')+
  guides(colour = F, size = F)+
  theme(
    axis.ticks.y  = element_blank(),
    axis.ticks.x  = element_blank(),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size=13,colour = "black"),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size=15)
  )

sel5plot<- ggplot(data = manhattan5, aes(x = CUM_POS/1e6, y = WEIGHTED_FST))+
  geom_point(aes(colour = cut(WEIGHTED_FST, c(-Inf, 0.46, Inf)), size = cut(WEIGHTED_FST, c(-Inf, 0.46, Inf)))) +
  scale_x_continuous('Position in Genome', limits = c(0, max(manhattan5$CUM_POS)/1e6))+
  scale_y_continuous(expression(italic(F[ST])), limits = c(0,0.75))+
  geom_hline(aes(yintercept = 0.46), lty = 2, alpha = 0.5)+
  scale_color_manual(values = c("black","red"))+
  scale_size_manual(values = c(0.7,2.7))+
  ggtitle('Generation 7000')+
  guides(colour = F, size = F)+
  theme(
    axis.ticks.y  = element_blank(),
    axis.ticks.x  = element_blank(),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size=13,colour = "black"),
    axis.title.x = element_text(size=15),
    axis.title.y = element_text(size=15)
  )



png(file="~/UBC/FishersWaveFst/DemoFigure/plots/RedManhattan.png",width=2000,height=450,res=200)
print(sel3plot)
dev.off()


a<-ggplot(data = yy, aes( x = gen, y = WEIGHTED_FST, group = id, col = id))+
  geom_vline(xintercept= 6000, lty =3 , alpha = 0.5, lty = 0.8)+
  geom_line(lwd = 1.4)+
  # geom_point(data = fst, aes( x = gen, y = WEIGHTED_FST), lwd = 5, col = 'red')+
  scale_y_continuous(expression(italic(F[ST])), limits = c(0,0.75))+
  scale_x_continuous(expression(italic('Time →')), limits = c(5000,7000))+
  guides(color = F)+
  scale_color_manual(values = pal)+
  geom_hline(yintercept= 0.46, lty = 2, alpha = 0.5, col = 'red')+
  guides(colour = F)+
  theme(
    axis.ticks.y  = element_blank(),
    axis.ticks.x  = element_blank(),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    plot.title = element_text(hjust = 0.5),
  #  axis.text.x = element_blank(),
    axis.text.y = element_text(size=13,colour = "black"),
    axis.title.x = element_text(size=15),
    axis.title.y = element_text(size=15)
  )

a

nam <- c("CHROM" ,       "BIN_START",    "BIN_END"   ,   "N_VARIANTS",   "WEIGHTED_FST", "MEAN_FST"   ,
         "REP_x"  ,      "gen"       ,   "id" ,          "MutId"     ,   "MutGen",       "firstGen"    ,
         "POS"     ,     "q"          ,   "s"   ,"REP_y"       )

yyy <- melt(yy, id = nam)
str(yyy)
b<-ggplot(data = yyy, aes( x = gen, y = value,  col = id, lty = variable))+
  geom_vline(xintercept= 6000, lty =3 , alpha = 0.5, lty = 0.8)+
  geom_line(lwd = 1.4)+
  # geom_point(data = fst, aes( x = gen, y = WEIGHTED_FST), lwd = 5, col = 'red')+
  scale_y_continuous(expression(italic("Allele Frequency")), limits = c(0,1.0))+
  scale_x_continuous(expression(italic('Time →')), limits = c(5000,7000))+
  guides(color = F)+
  scale_color_manual(values = pal)+
  guides(colour = F, lty = F)+
  theme(
    axis.ticks.y  = element_blank(),
    axis.ticks.x  = element_blank(),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size=13,colour = "black"),
    axis.title.x = element_text(size=15),
    axis.title.y = element_text(size=15)
  )

b


skyline<-ggarrange(sel1plot,
                   sel2plot,
                   sel3plot,
                   sel4plot,
                   sel5plot,
                   nrow = 5,
                   ncol = 1,
                   labels = "" )

zx<-ggarrange(manhattan.plot,                                                 
              ggarrange( b, a,  ncol = 2, labels = c("", "")), 
              nrow = 2, 
              labels = "" )
png(file="~/UBC/FishersWaveFst/DemoFigure/plots/FigureS1.png",width=2000,height=450*5,res=200)
print(skyline)
dev.off()

png(file="~/UBC/FishersWaveFst/DemoFigure/plots/combined_3panel.png",width=2000,height=900,res=200)
print(zx)
dev.off()


