library(data.table)
library(ggplot2)

data <- fread(file="data-parser/alldata-static.csv", header=T)
setnames(data, 1, "algorithm")
setnames(data, "numberGO", "GO")

accuracy_levels <-  c(0.1, 0.01, 0.001, 0.0001, 0.00001)

DT <- unique(data[,.(numGOacc0,numGOacc1, numGOacc2, numGOacc3, numGOacc4,  GO), by=.(algorithm,problem, run)])
# calculate PR per accuracy level
DT[, PR0:=numGOacc0/GO, by=list(algorithm,problem)]
DT[, PR1:=numGOacc1/GO, by=list(algorithm,problem)]
DT[, PR2:=numGOacc2/GO, by=list(algorithm,problem)]
DT[, PR3:=numGOacc3/GO, by=list(algorithm,problem)]
DT[, PR4:=numGOacc4/GO, by=list(algorithm,problem)]

#All PR together 
DTAll <- DT[,.(PR0,PR1,PR2,PR3,PR4),by=.(algorithm)]


suppressWarnings(suppressMessages(library(reshape2)))
DTAll_melted <- melt(DTAll, id.vars = c("algorithm"))
finalMeanRanking <- DTAll_melted[, .(meanPR = mean(value)), by=.(algorithm)]

print(finalMeanRanking)
results.classic <- copy(finalMeanRanking)
fwrite(finalMeanRanking, "results/Analysis1-classic-mean-PR-all-acc.csv", row.names = FALSE)

## Create appropriate graphs

bp <- ggplot(DTAll_melted,aes(x=factor(algorithm), y=c(value) ) ) + 
  geom_boxplot(notch=FALSE,aes(fill = factor(algorithm) ) ) + 
  stat_summary(fun.y="mean", geom="point", shape=23, size=3, fill="black") + 
  theme(axis.text.x  = element_text(angle=45, vjust=1,hjust=1) ) 	+ 
  ggtitle("All Accuracy levels") 	+ 
  theme(axis.title.x = element_blank() ) 	+ 
  theme(plot.background = element_blank() )  +
  ylab("Peak Ratio in all benchmark functions") +
  guides( fill=guide_legend(title="Algorithms") )

print(bp)
cairo_pdf("figs/bp-pr-acc-sum.pdf",bg = "transparent")
print(bp)
dev.off()

#Accuracy level 1
bp <- ggplot(DT,aes(x=factor(algorithm), y=PR0) ) +
  geom_boxplot(notch=FALSE,aes(fill = factor(algorithm) ) ) +
  stat_summary(fun.y="mean", geom="point", shape=23, size=3, fill="black") +
  theme(axis.text.x  = element_text(angle=45, vjust=1,hjust=1)) +
  ggtitle("Accuracy level 1.0e-1") +
  theme(axis.title.x = element_blank() ) +
  theme(plot.background = element_blank() )  +
  ylab("Peak Ratio in all benchmark functions") +
  guides( fill=guide_legend(title="Algorithms") )

print(bp)
cairo_pdf("figs/bp-pr-acc1.pdf",bg = "transparent")
print(bp)
dev.off()

#Accuracy level 2
bp <- ggplot(DT,aes(x=factor(algorithm), y=PR1) ) +
  geom_boxplot(notch=FALSE,aes(fill = factor(algorithm) ) ) +
  stat_summary(fun.y="mean", geom="point", shape=23, size=3, fill="black") +
  theme(axis.text.x  = element_text(angle=45, vjust=1,hjust=1) )  +
  ggtitle("Accuracy level 1.0e-2")  +
  theme(axis.title.x = element_blank() )  +
  theme(plot.background = element_blank() )  +
  ylab("Peak Ratio in all benchmark functions") +
  guides( fill=guide_legend(title="Algorithms") )

print(bp)
cairo_pdf("figs/bp-pr-acc2.pdf",bg = "transparent")
print(bp)
dev.off()

#Accuracy level 3
bp <- ggplot(DT,aes(x=factor(algorithm), y=PR2) )  +
  geom_boxplot(notch=FALSE,aes(fill = factor(algorithm) ) )  +
  stat_summary(fun.y="mean", geom="point", shape=23, size=3, fill="black") +
  theme(axis.text.x  = element_text(angle=45, vjust=1,hjust=1) )  +
  ggtitle("Accuracy level 1.0e-3")  +
  theme(axis.title.x = element_blank() )  +
  theme(plot.background = element_blank() )  +
  ylab("Peak Ratio in all benchmark functions") +
  guides( fill=guide_legend(title="Algorithms") )

print(bp)
cairo_pdf("figs/bp-pr-acc3.pdf",bg = "transparent")
print(bp)
dev.off()

#Accuracy level 4
bp <- ggplot(DT,aes(x=factor(algorithm), y=PR3) )  +
  geom_boxplot(notch=FALSE,aes(fill = factor(algorithm) ) )  +
  stat_summary(fun.y="mean", geom="point", shape=23, size=3, fill="black") +
  theme(axis.text.x  = element_text(angle=45, vjust=1,hjust=1) )  +
  ggtitle("Accuracy level 1.0e-4")  +
  theme(axis.title.x = element_blank() )  +
  theme(plot.background = element_blank() )  +
  ylab("Peak Ratio in all benchmark functions") +
  guides( fill=guide_legend(title="Algorithms") )

print(bp)
cairo_pdf("figs/bp-pr-acc4.pdf",bg = "transparent")
print(bp)
dev.off()

#Accuracy level 5
bp <- ggplot(DT,aes(x=factor(algorithm), y=PR4) )  +
  geom_boxplot(notch=FALSE,aes(fill = factor(algorithm) ) )  +
  stat_summary(fun.y="mean", geom="point", shape=23, size=3, fill="black") +
  theme(axis.text.x  = element_text(angle=45, vjust=1,hjust=1) )  +
  ggtitle("Accuracy level 1.0e-5")  +
  theme(axis.title.x = element_blank() )  +
  theme(plot.background = element_blank() )  +
  ylab("Peak Ratio in all benchmark functions") +
  guides( fill=guide_legend(title="Algorithms") )

print(bp)
cairo_pdf("figs/bp-pr-acc5.pdf",bg = "transparent")
print(bp)
dev.off()

