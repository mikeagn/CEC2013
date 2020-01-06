library(data.table)
library(ggplot2)

accuracy_levels <-  c(0.1, 0.01, 0.001, 0.0001, 0.00001)
nruns <- 50
nproblems <- 20

#DT <- fread(file="data-parser/alldata-static.csv", sep=",", header=T)
DT <- fread(file="final-data-gecco-static.csv", sep=",", header=T)
setnames(DT,1,"algorithm")
setnames(DT,"numberGO","GO")

DT[, SO:=max(counter), by=list(algorithm,problem,run)]

F1 <- function(precision, recall){ 
	if( (precision==0) && (recall==0) ) return(0)
	return(2*(precision * recall)/(precision + recall)) 
}
### F1 produces a lot of NaNs if recall==0 and precision==0, set to 0 manually 

# Calculate precision, recall 
# Given a set of solutions SO and a set of all global optima GO

# TP = |relevant documents ^ retrieved documents|
# FP = |retrieved documents| - TP
# FN = |relevant documents| - TP
# TN
# precision = |relevant documents ^ retrieved documents|/|retrieved documents|
# recall = |relevant documents ^ retrieved documents|/|relevant documents|
# relevant documents = GO # 
# retrieved documents = SO
# TGOSOX (X <- acc) = |relevant documents ^ retrieved documents| == numGOaccX
# F1 = 2*(precision * recall)/(precision + recall)

# acc0
DT[, precision0:=numGOacc0/SO, by=list(algorithm,problem,run)]
DT[, recall0:=numGOacc0/GO, by=list(algorithm,problem,run)]
DT[, F10:=F1(precision0, recall0), by=list(algorithm,problem,run)]

# acc1
DT[, precision1:=numGOacc1/SO, by=list(algorithm,problem,run)]
DT[, recall1:=numGOacc1/GO, by=list(algorithm,problem,run)]
DT[, F11:=F1(precision1, recall1), by=list(algorithm,problem,run)]

# acc2
DT[, precision2:=numGOacc2/SO, by=list(algorithm,problem,run)]
DT[, recall2:=numGOacc2/GO, by=list(algorithm,problem,run)]
DT[, F12:=F1(precision2, recall2), by=list(algorithm,problem,run)]

# acc3
DT[, precision3:=numGOacc3/SO, by=list(algorithm,problem,run)]
DT[, recall3:=numGOacc3/GO, by=list(algorithm,problem,run)]
DT[, F13:=F1(precision3, recall3), by=list(algorithm,problem,run)]

# acc4
DT[, precision4:=numGOacc4/SO, by=list(algorithm,problem,run)]
DT[, recall4:=numGOacc4/GO, by=list(algorithm,problem,run)]
DT[, F14:=F1(precision4, recall4), by=list(algorithm,problem,run)]

DTF1All <- unique(DT[,.(F10,F11,F12,F13,F14), by=.(algorithm,problem, run)])
DTAll <- unique(DT[,.(F10,F11,F12,F13,F14,
                      precision0,recall0,
                      precision1,recall1,
                      precision2,recall2,
                      precision3,recall3,
                      precision4,recall4), by=.(algorithm,problem, run)])

DTF1All_melted <- melt(DTF1All, id.vars = c("algorithm","problem", "run"))
finalMeanRanking <- DTF1All_melted[, .(meanF1 = mean(value)), by=.(algorithm)]

print(finalMeanRanking)
results.staticF1 <- copy(finalMeanRanking)
fwrite(finalMeanRanking, "results/Analysis2-staticF1-mean-F1-all-acc.csv", row.names = FALSE)

fbeta <- function(beta, precision, recall){
  result = (1.0 + beta*beta)*(precision*recall)/(beta*beta*precision+recall)
  return(result)
}

# boxplot of F1 acc 0
bp <- ggplot(DTAll)+ 
  geom_boxplot(notch=F, aes(x=factor(algorithm), y=F10,fill = factor(algorithm) ) ) + 
  theme(text = element_text(size=20)) +
  stat_summary(aes(x=factor(algorithm), y=F10), fun.y = "mean", geom = "point", shape= 23, size= 3, fill= "red") +
  theme(plot.background = element_blank() )  + 
  theme(legend.position="none") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  ylab("F1 value Accuracy level 1.0e-1") +
  xlab(NULL) +
  ylim(0,1)
print(bp);
pdf("figs/boxed-boxplot-algorithm-F10.pdf",family="URWHelvetica",bg = "transparent", width=7, height=5)
print(bp) 
dev.off()

# boxplot of precision acc 0
bp <- ggplot(DTAll)+ 
  geom_boxplot(notch=F, aes(x=factor(algorithm), y=precision0,fill = factor(algorithm) ) ) + 
  theme(text = element_text(size=20)) +
  stat_summary(aes(x=factor(algorithm), y=precision0), fun.y = "mean", geom = "point", shape= 23, size= 3, fill= "red") +
  theme(plot.background = element_blank() )  + 
  theme(legend.position="none") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  ylab("Precision Accuracy level 1.0e-1") +
  xlab(NULL)  +
  ylim(0,1)
print(bp);
pdf("figs/boxed-boxplot-algorithm-precision0.pdf",family="URWHelvetica",bg = "transparent", width=7, height=5)
print(bp) 
dev.off()

# boxplot of recall acc 0
bp <- ggplot(DTAll)+ 
  geom_boxplot(notch=F, aes(x=factor(algorithm), y=recall0,fill = factor(algorithm) ) ) + 
  theme(text = element_text(size=20)) +
  stat_summary(aes(x=factor(algorithm), y=recall0), fun.y = "mean", geom = "point", shape= 23, size= 3, fill= "red") +
  theme(plot.background = element_blank() )  + 
  theme(legend.position="none") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  ylab("Recall Accuracy level 1.0e-1") +
  xlab(NULL)  +
  ylim(0,1)
print(bp);
pdf("figs/boxed-boxplot-algorithm-recall0.pdf",family="URWHelvetica",bg = "transparent", width=7, height=5)
print(bp) 
dev.off()

# boxplot of F1 acc 1
bp <- ggplot(DTAll)+ 
  geom_boxplot(notch=F, aes(x=factor(algorithm), y=F11,fill = factor(algorithm) ) ) + 
  theme(text = element_text(size=20)) +
  stat_summary(aes(x=factor(algorithm), y=F11), fun.y = "mean", geom = "point", shape= 23, size= 3, fill= "red") +
  theme(plot.background = element_blank() )  + 
  theme(legend.position="none") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  ylab("F1 value Accuracy level 1.0e-2") +
  xlab(NULL) +
  ylim(0,1)
print(bp);
pdf("figs/boxed-boxplot-algorithm-F11.pdf",family="URWHelvetica",bg = "transparent", width=7, height=5)
print(bp) 
dev.off()

# boxplot of precision acc 1
bp <- ggplot(DTAll)+ 
  geom_boxplot(notch=F, aes(x=factor(algorithm), y=precision1,fill = factor(algorithm) ) ) + 
  theme(text = element_text(size=20)) +
  stat_summary(aes(x=factor(algorithm), y=precision1), fun.y = "mean", geom = "point", shape= 23, size= 3, fill= "red") +
  theme(plot.background = element_blank() )  + 
  theme(legend.position="none") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  ylab("Precision Accuracy level 1.0e-2") +
  xlab(NULL)  +
ylim(0,1)
print(bp);
pdf("figs/boxed-boxplot-algorithm-precision1.pdf",family="URWHelvetica",bg = "transparent", width=7, height=5)
print(bp) 
dev.off()

# boxplot of recall acc 1
bp <- ggplot(DTAll)+ 
  geom_boxplot(notch=F, aes(x=factor(algorithm), y=recall1,fill = factor(algorithm) ) ) + 
  theme(text = element_text(size=20)) +
  stat_summary(aes(x=factor(algorithm), y=recall1), fun.y = "mean", geom = "point", shape= 23, size= 3, fill= "red") +
  theme(plot.background = element_blank() )  + 
  theme(legend.position="none") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  ylab("Recall Accuracy level 1.0e-2") +
  xlab(NULL)  +
ylim(0,1)
print(bp);
pdf("figs/boxed-boxplot-algorithm-recall1.pdf",family="URWHelvetica",bg = "transparent", width=7, height=5)
print(bp) 
dev.off()

# boxplot of F1 acc 2
bp <- ggplot(DTAll)+ 
  geom_boxplot(notch=F, aes(x=factor(algorithm), y=F12,fill = factor(algorithm) ) ) + 
  theme(text = element_text(size=20)) +
  stat_summary(aes(x=factor(algorithm), y=F12), fun.y = "mean", geom = "point", shape= 23, size= 3, fill= "red") +
  theme(plot.background = element_blank() )  + 
  theme(legend.position="none") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  ylab("F1 value Accuracy level 1.0e-3") +
  xlab(NULL) +
  ylim(0,1)
print(bp);
pdf("figs/boxed-boxplot-algorithm-F12.pdf",family="URWHelvetica",bg = "transparent", width=7, height=5)
print(bp) 
dev.off()

# boxplot of precision acc 2
bp <- ggplot(DTAll)+ 
  geom_boxplot(notch=F, aes(x=factor(algorithm), y=precision2,fill = factor(algorithm) ) ) + 
  theme(text = element_text(size=20)) +
  stat_summary(aes(x=factor(algorithm), y=precision2), fun.y = "mean", geom = "point", shape= 23, size= 3, fill= "red") +
  theme(plot.background = element_blank() )  + 
  theme(legend.position="none") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  ylab("Precision Accuracy level 1.0e-3") +
  xlab(NULL)  +
ylim(0,1)
print(bp);
pdf("figs/boxed-boxplot-algorithm-precision2.pdf",family="URWHelvetica",bg = "transparent", width=7, height=5)
print(bp) 
dev.off()

# boxplot of recall acc 2
bp <- ggplot(DTAll)+ 
  geom_boxplot(notch=F, aes(x=factor(algorithm), y=recall2,fill = factor(algorithm) ) ) + 
  theme(text = element_text(size=20)) +
  stat_summary(aes(x=factor(algorithm), y=recall2), fun.y = "mean", geom = "point", shape= 23, size= 3, fill= "red") +
  theme(plot.background = element_blank() )  + 
  theme(legend.position="none") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  ylab("Recall Accuracy level 1.0e-3") +
  xlab(NULL)  +
ylim(0,1)
print(bp);
pdf("figs/boxed-boxplot-algorithm-recall2.pdf",family="URWHelvetica",bg = "transparent", width=7, height=5)
print(bp) 
dev.off()

# boxplot of F1 acc 3
bp <- ggplot(DTAll)+ 
  geom_boxplot(notch=F, aes(x=factor(algorithm), y=F13,fill = factor(algorithm) ) ) + 
  theme(text = element_text(size=20)) +
  stat_summary(aes(x=factor(algorithm), y=F13), fun.y = "mean", geom = "point", shape= 23, size= 3, fill= "red") +
  theme(plot.background = element_blank() )  + 
  theme(legend.position="none") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  ylab("F1 value Accuracy level 1.0e-4") +
  xlab(NULL) +
  ylim(0,1)
print(bp);
pdf("figs/boxed-boxplot-algorithm-F13.pdf",family="URWHelvetica",bg = "transparent", width=7, height=5)
print(bp) 
dev.off()

# boxplot of precision acc 2
bp <- ggplot(DTAll)+ 
  geom_boxplot(notch=F, aes(x=factor(algorithm), y=precision3,fill = factor(algorithm) ) ) + 
  theme(text = element_text(size=20)) +
  stat_summary(aes(x=factor(algorithm), y=precision3), fun.y = "mean", geom = "point", shape= 23, size= 3, fill= "red") +
  theme(plot.background = element_blank() )  + 
  theme(legend.position="none") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  ylab("Precision Accuracy level 1.0e-4") +
  xlab(NULL)  +
ylim(0,1)
print(bp);
pdf("figs/boxed-boxplot-algorithm-precision3.pdf",family="URWHelvetica",bg = "transparent", width=7, height=5)
print(bp) 
dev.off()

# boxplot of recall acc 3
bp <- ggplot(DTAll)+ 
  geom_boxplot(notch=F, aes(x=factor(algorithm), y=recall3,fill = factor(algorithm) ) ) + 
  theme(text = element_text(size=20)) +
  stat_summary(aes(x=factor(algorithm), y=recall3), fun.y = "mean", geom = "point", shape= 23, size= 3, fill= "red") +
  theme(plot.background = element_blank() )  + 
  theme(legend.position="none") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  ylab("Recall Accuracy level 1.0e-4") +
  xlab(NULL)  +
ylim(0,1)
print(bp);
pdf("figs/boxed-boxplot-algorithm-recall3.pdf",family="URWHelvetica",bg = "transparent", width=7, height=5)
print(bp) 
dev.off()


# boxplot of F1 acc 4
bp <- ggplot(DTAll)+ 
  geom_boxplot(notch=F, aes(x=factor(algorithm), y=F14,fill = factor(algorithm) ) ) + 
  theme(text = element_text(size=20)) +
  stat_summary(aes(x=factor(algorithm), y=F14), fun.y = "mean", geom = "point", shape= 23, size= 3, fill= "red") +
  theme(plot.background = element_blank() )  + 
  theme(legend.position="none") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  ylab("F1 value Accuracy level 1.0e-5") +
  xlab(NULL) +
  ylim(0,1)
print(bp);
pdf("figs/boxed-boxplot-algorithm-F14.pdf",family="URWHelvetica",bg = "transparent", width=7, height=5)
print(bp) 
dev.off()

# boxplot of precision acc 4
bp <- ggplot(DTAll)+ 
  geom_boxplot(notch=F, aes(x=factor(algorithm), y=precision4,fill = factor(algorithm) ) ) + 
  theme(text = element_text(size=20)) +
  stat_summary(aes(x=factor(algorithm), y=precision4), fun.y = "mean", geom = "point", shape= 23, size= 3, fill= "red") +
  theme(plot.background = element_blank() )  + 
  theme(legend.position="none") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  ylab("Precision Accuracy level 1.0e-5") +
  xlab(NULL)  +
ylim(0,1)
print(bp);
pdf("figs/boxed-boxplot-algorithm-precision4.pdf",family="URWHelvetica",bg = "transparent", width=7, height=5)
print(bp) 
dev.off()

# boxplot of recall acc 4
bp <- ggplot(DTAll)+ 
  geom_boxplot(notch=F, aes(x=factor(algorithm), y=recall4,fill = factor(algorithm) ) ) + 
  theme(text = element_text(size=20)) +
  stat_summary(aes(x=factor(algorithm), y=recall4), fun.y = "mean", geom = "point", shape= 23, size= 3, fill= "red") +
  theme(plot.background = element_blank() )  + 
  theme(legend.position="none") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  ylab("Recall Accuracy level 1.0e-5") +
  xlab(NULL)  +
ylim(0,1)
print(bp);
pdf("figs/boxed-boxplot-algorithm-recall4.pdf",family="URWHelvetica",bg = "transparent", width=7, height=5)
print(bp) 
dev.off()

# Overall accuracy levels
bp <- ggplot(DTF1All_melted,aes(x=factor(algorithm), y=c(value) ) ) + 
  geom_boxplot(notch=FALSE,aes(fill = factor(algorithm) ) ) + 
  stat_summary(fun.y="mean", geom="point", shape=23, size=3, fill="black") + 
  theme(axis.text.x  = element_text(angle=45, vjust=1,hjust=1) ) 	+ 
  ggtitle("All Accuracy levels") 	+ 
  theme(axis.title.x = element_blank() ) 	+ 
  theme(plot.background = element_blank() )  +
  ylab("F1 value in all benchmark functions") +
  guides( fill=guide_legend(title="Algorithms") )

print(bp)
pdf("figs/F1-Overall-Acc.pdf",family="URWHelvetica",bg = "transparent")
print(bp)
dev.off()

