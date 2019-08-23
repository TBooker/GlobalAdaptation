rm(list = ls() )

library(ggplot2)
library(reshape2)
library(ggpubr)

pal <- c( "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

trace <- read.csv('~/UBC/FishersWaveFst/DemoFigure/gen6000.M1_p0.0001.gz')

yy<- trace#[xx$id %in% sub$id,]

yy[(yy$gen < yy$firstGen)&(is.na(yy$q) ) ,]$q <- 0
yy[(yy$gen < yy$firstGen)&(is.na(yy$q1) ) ,]$q1 <- 0
yy[(yy$gen < yy$firstGen)&(is.na(yy$q2) ) ,]$q2 <- 0

yy[(yy$gen > yy$firstGen)&(is.na(yy$q) ) ,]$q <- 1
yy[(yy$gen > yy$firstGen)&(is.na(yy$q1) ) ,]$q1 <- 1
yy[(yy$gen > yy$firstGen)&(is.na(yy$q2) ) ,]$q2 <- 1


manhattan <- read.csv('~/UBC/FishersWaveFst/DemoFigure/M1_p0.0001_4_files/gen6000.M1_p0.0001.csv.gz')
manhattan <- manhattan[manhattan$BIN_START <2130001,]
manhattan<- manhattan[order(manhattan$REP),]
manhattan$id <- paste(manhattan$REP,manhattan$BIN_START)
manhattan$POS <- manhattan$BIN_START + (manhattan$BIN_END-manhattan$BIN_START)/2
manhattan$CUM_POS <- manhattan$POS + (manhattan$REP * 2000000) - 2000000 

manhattan.new = manhattan[seq(1, nrow(manhattan), 10), ]

slice <- yy[yy$BIN_START <2130001,]
slice<- slice[order(slice$REP_x),]
slice$id <- paste(slice$REP_x,slice$BIN_START, sep = '_')
slice$POS <- slice$BIN_START + (slice$BIN_END-slice$BIN_START)/2
slice$CUM_POS <- slice$POS + (slice$REP_x * 2000000) - 2000000 
slice <- slice[slice$gen == 6000,]

maxFst<- manhattan.new[manhattan.new$WEIGHTED_FST == max(manhattan.new$WEIGHTED_FST), ]

manhattan.plot<- ggplot( manhattan.new, aes(x = CUM_POS/1e6, y = WEIGHTED_FST))+
  geom_point()+
  geom_point(data = slice, aes(x = CUM_POS/1e6, y = WEIGHTED_FST, col = id, group = id), size  = 4)+
  scale_color_manual(values = pal)+
  scale_y_continuous(expression(italic(F[ST])), limits = c(0,0.75))+
  scale_x_continuous('Position in Genome (Mbp)')+
  geom_hline(yintercept= 0.46, lty = 2, alpha = 0.5, col = 'red')+
  guides(colour = F)+
  theme(
    axis.ticks.y  = element_blank(),
   # axis.ticks.x  = element_blank(),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    plot.title = element_text(hjust = 0.5),
  #  axis.text.x = element_blank(),
    axis.text.y = element_text(size=13,colour = "black"),
    axis.title.x = element_text(size=15),
    axis.title.y = element_text(size=15)
  )
  





sel1plot<- ggplot(data = manhattan.new, aes(x = CUM_POS/1e6, y = WEIGHTED_FST))+
  geom_point(aes(colour = cut(WEIGHTED_FST, c(-Inf, 0.46, Inf)), size = cut(WEIGHTED_FST, c(-Inf, 0.46, Inf)))) +
  scale_x_continuous('Position in Genome', limits = c(0, max(manhattan$CUM_POS)/1e6))+
  scale_y_continuous(expression(italic(F[ST])), limits = c(0,0.75))+
  annotate('point', x = maxFst$CUM_POS/1e6, y = maxFst$WEIGHTED_FST+0.05, shape = 8)+
  geom_hline(aes(yintercept = 0.46), lty = 2, alpha = 0.5)+
  scale_color_manual(values = c("black","red"))+
  scale_size_manual(values = c(0.7,2.7))+
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
print(sel1plot)
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
  #  axis.ticks.x  = element_blank(),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    plot.title = element_text(hjust = 0.5),
    #axis.text.x = element_blank(),
  axis.text.x = element_text(size=13,colour = "black"),
  axis.text.y = element_text(size=13,colour = "black"),
  axis.title.x = element_text(size=15),
    axis.title.y = element_text(size=15)
  )



nam <- c("CHROM" ,       "BIN_START",    "BIN_END"   ,   "N_VARIANTS",   "WEIGHTED_FST", "MEAN_FST"   ,
         "REP_x"  ,      "gen"       ,   "id" ,          "MutId"     ,   "MutGen",       "firstGen"    ,
         "POS"     ,     "q"          ,   "s"   ,"REP_y"       )

yyy <- melt(yy, id = nam)
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
    #axis.ticks.x  = element_blank(),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_blank(),
#    axis.text.x = element_text(size=13,colour = "black"),
    axis.text.y = element_text(size=13,colour = "black"),
    axis.title.x = element_text(size=15),
    axis.title.y = element_text(size=15)
  )




zx<-ggarrange(manhattan.plot,                                                 
              ggarrange( b, a,  ncol = 2, labels = c("", "")), 
              nrow = 2, 
              labels = "" )

png(file="combined_3panel.png",width=2000,height=900,res=200)
print(zx)
dev.off()

png(file="combined_3panel.jpg",width=2000,height=900,res=200)
print(zx)
dev.off()

zy<-ggarrange( b, a,  ncol = 1, nrow = 2, labels = c("A", "B"), align = 'v')
    
png(file="2panel.png",width=2000,height=900,res=200)
print(zy)
dev.off()

