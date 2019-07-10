library(data.table)
library(ggplot2)
library(progress)

nruns <- 50
accuracy_levels <-  c(0.1, 0.01, 0.001, 0.0001, 0.00001)
maxRL = c(5E4,5E4,5E4,5E4,5E4, 
          2E5,2E5,
          4E5,4E5, 
          2E5,2E5,2E5,2E5,
          4E5,4E5,4E5,4E5,4E5,4E5,4E5)
budget <- data.table(FEs=maxRL)
budget[,problem:=paste0("prob",.I)]

DT <- fread(file="data-parser/alldata-dynamic.csv", sep=",", header=T)
setnames(DT, 1, "algorithm")
setnames(DT, "numberGO", "GO")
# Ensure that values are sorted correctly
setkeyv(DT, c("algorithm","problem", "run", "fes"))
DT[, SO:= as.integer(popsize)]

algorithms <- unique(DT[,algorithm])
probNames <- unique(DT[,problem])
nproblems <- length(probNames)

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

# compute dynamic F1/AUC, takes preselected data for alg, prob, run (right Riemann sum)
dynamicF1 <- function(subDT, acc, probName, algName) {
  lastTime=0
  f1sum=0
  f1part=0
  
  subDTrows = nrow(subDT)	
  if( subDTrows > 0 ) {
    
    GO = subDT[1,GO]
    for( lineNo in 1:subDTrows ) {  
      
      # find appropriate accuracy value
      if( acc==0.1 ) {
        numGOacc = subDT[lineNo,numGOacc0]
      } else if( acc==0.01 ) {
        numGOacc = subDT[lineNo,numGOacc1]
      } else if( acc==0.001 ) {
        numGOacc = subDT[lineNo,numGOacc2]
      } else if( acc==0.0001 ) {
        numGOacc = subDT[lineNo,numGOacc3]
      } else if( acc==0.00001 ) {
        numGOacc = subDT[lineNo,numGOacc4]
      }
      
      # compute temporary F1
      precision = numGOacc / subDT[lineNo,SO]
      recall = numGOacc / GO
      f1part = F1( precision, recall )
      f1sum = f1sum + ((subDT[lineNo, fes] - lastTime)) * f1part
      lastTime = subDT[lineNo, fes]  
    }
    # we have to add the last interval from lastTime to the max amount of evals
    maxFes = maxRL[which(probNames==probName)]
    f1sum = f1sum + (maxFes - lastTime) * f1part
    return(abs(f1sum) / maxFes)
  } else {
    #print(paste("WARNING! empty data for alg ", algName, ", prob ", probName, ", run ", runNo))
    return(0)
  }
}

# compute dynamic F1/AUC, takes preselected data for alg, prob, run (left Riemann sum)
dynamicF1v2 <- function(subDT, acc, probName, algName, probBudget) {
  lastTime=0
  f1sum=0
  f1part=0
  
  subDTrows = nrow(subDT)	
  if( subDTrows > 0 ) {
    
    GO = subDT[1,GO]
    for( lineNo in 1:subDTrows ) {  
      
      # find appropriate accuracy value
      if( acc==0.1 ) {
        numGOacc = subDT[lineNo,numGOacc0]
      } else if( acc==0.01 ) {
        numGOacc = subDT[lineNo,numGOacc1]
      } else if( acc==0.001 ) {
        numGOacc = subDT[lineNo,numGOacc2]
      } else if( acc==0.0001 ) {
        numGOacc = subDT[lineNo,numGOacc3]
      } else if( acc==0.00001 ) {
        numGOacc = subDT[lineNo,numGOacc4]
      }
      
      # compute temporary F1
      precision = numGOacc / subDT[lineNo,SO]
      recall = numGOacc / GO
      f1part = F1( precision, recall )
      if (lineNo == subDTrows) {
        f1sum <- f1sum + f1part * (probBudget - subDT[lineNo, fes]) / probBudget
      } else {
        f1sum <- f1sum + f1part * (subDT[lineNo+1, fes] - subDT[lineNo, fes]) / probBudget
      }
    }
    return(f1sum)
  } else {
    #print(paste("WARNING! empty data for alg ", algName, ", prob ", probName, ", run ", runNo))
    return(0)
  }
}

# set up data frame for detailed results 
results <- data.frame(algorithm=character(), problem=character(), run=integer(), acc1=double(), 
                      acc2=double(), acc3=double(), acc4=double(), acc5=double(),
                      stringsAsFactors=FALSE) 

pb <- progress_bar$new(
  format = "(:spin) [:bar] :percent",
  total = length(algorithms) * length(probNames), clear = FALSE, width = 60)

aucresults <- data.table()
# algorithm loop: sum everything up
for( alg in algorithms ) {
  counter = 0
  aucSum = 0

  for( prob in probNames ) {
    probBudget <- budget[prob == problem]$FEs
    for( runNo in 0:(nruns-1) ) {
      accSums = rep(0, length(accuracy_levels))
      # same data for all accuracy levels
      subDT = subset(DT, algorithm==alg & problem==prob & run==runNo )   
      
      accNum = 1
      for( accu in accuracy_levels ) {
        counter = counter + 1
        tmp = dynamicF1v2(subDT, acc=accu, probName=prob, algName=alg, probBudget = probBudget)
        aucSum = aucSum + tmp
        accSums[accNum] = accSums[accNum] + tmp
        accNum = accNum + 1
      }	
      results[nrow(results)+1,] = c(alg, prob, runNo, accSums[1], accSums[2], accSums[3], accSums[4], accSums[5])
    }
    pb$tick()
  }
  aucAvg = aucSum / counter
  aucresults <- rbind(aucresults, data.table(algorithm=alg, AUCAvg=aucAvg))
}
print(aucresults)
results.dynF1 <- copy(aucresults)
fwrite(aucresults, "results/Analysis3-dynamicF1-mean-AUC-all-acc.csv", quote=FALSE, row.names = FALSE)

DTResults <- data.table(results)
DTResults[, acc1:= as.numeric(acc1)]
DTResults[, acc2:= as.numeric(acc1)]
DTResults[, acc3:= as.numeric(acc1)]
DTResults[, acc4:= as.numeric(acc1)]
DTResults[, acc5:= as.numeric(acc1)]
fwrite(DTResults, "results/Analysis3-dynamicF1-data.csv", quote=FALSE, row.names = FALSE)

accuracy_levels <-  c(0.1, 0.01, 0.001, 0.0001, 0.00001)

DTResults[, dynF1PR1:=mean(acc1), by=list(algorithm,problem)]
DTResults[, dynF1PR2:=mean(acc2), by=list(algorithm,problem)]
DTResults[, dynF1PR3:=mean(acc3), by=list(algorithm,problem)]
DTResults[, dynF1PR4:=mean(acc4), by=list(algorithm,problem)]
DTResults[, dynF1PR5:=mean(acc5), by=list(algorithm,problem)]

#Accuracy level 1
bp <- ggplot(DTResults,aes(x=factor(algorithm), y=problem) ) +
  geom_tile(aes(fill = dynF1PR1 ) ) +
  theme(axis.text.x  = element_text(angle=45, vjust=1,hjust=1)) +
  ggtitle("Accuracy level 1.0e-1") +
  theme(axis.title.x = element_blank() ) +
  theme(plot.background = element_blank() )  +
  ylab("Benchmark functions") +
  xlab(NULL) +
  theme_bw() 
print(bp)
cairo_pdf("figs/dynF1-acc1.pdf",bg = "transparent")
print(bp)
dev.off()

#Accuracy level 2
bp <- ggplot(DTResults,aes(x=factor(algorithm), y=problem) ) +
  geom_tile(aes(fill = dynF1PR2 ) ) +
  theme(axis.text.x  = element_text(angle=45, vjust=1,hjust=1)) +
  ggtitle("Accuracy level 1.0e-2") +
  theme(axis.title.x = element_blank() ) +
  theme(plot.background = element_blank() )  +
  ylab("Benchmark function") +
  xlab(NULL) +
  theme_bw() 
print(bp)
cairo_pdf("figs/dynF1-acc2.pdf",bg = "transparent")
print(bp)
dev.off()

#Accuracy level 3
bp <- ggplot(DTResults,aes(x=factor(algorithm), y=problem) ) +
  geom_tile(aes(fill = dynF1PR3 ) ) +
  theme(axis.text.x  = element_text(angle=45, vjust=1,hjust=1)) +
  ggtitle("Accuracy level 1.0e-3") +
  theme(axis.title.x = element_blank() ) +
  theme(plot.background = element_blank() )  +
  ylab("Benchmark functions") +
  xlab(NULL) +
  theme_bw() 
print(bp)
cairo_pdf("figs/dynF1-acc3.pdf",bg = "transparent")
print(bp)
dev.off()

#Accuracy level 4
bp <- ggplot(DTResults,aes(x=factor(algorithm), y=problem) ) +
  geom_tile(aes(fill = dynF1PR4 ) ) +
  theme(axis.text.x  = element_text(angle=45, vjust=1,hjust=1)) +
  ggtitle("Accuracy level 1.0e-4") +
  theme(axis.title.x = element_blank() ) +
  theme(plot.background = element_blank() )  +
  ylab("Benchmark functions") +
  xlab(NULL) +
  theme_bw() 
print(bp)
cairo_pdf("figs/dynF1-acc4.pdf",bg = "transparent")
print(bp)
dev.off()

#Accuracy level 5
bp <- ggplot(DTResults,aes(x=factor(algorithm), y=problem) ) +
  geom_tile(aes(fill = dynF1PR5 ) ) +
  theme(axis.text.x  = element_text(angle=45, vjust=1,hjust=1)) +
  ggtitle("Accuracy level 1.0e-5") +
  theme(axis.title.x = element_blank() ) +
  theme(plot.background = element_blank() )  +
  ylab("Benchmark functions") +
  xlab(NULL) +
  theme_bw() 
print(bp)
cairo_pdf("figs/dynF1-acc5.pdf",bg = "transparent")
print(bp)
dev.off()

