## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
AsSlides <- TRUE

## ----setup2, include=FALSE,eval=FALSE,echo=FALSE-------------------------
## library(ShortRead)
## 
## fqSample <- FastqSampler("~/Downloads/ENCFF001NQP.fastq.gz",n=10^6)
## temp <- yield(fqSample)
## 
## 
## writeFastq(fastqSample,file = "~/Projects/Software/Github/RUBioconductor_Introduction/r_course/Data/sampled_ENCFF000CXH.fastq.gz",mode = "w")

## ----shortreada,include=FALSE--------------------------------------------
library(ShortRead)

## ----shortread-----------------------------------------------------------
library(ShortRead)

## ----mycRep1Reads,cache=TRUE---------------------------------------------
fqSample <- FastqSampler("~/Downloads/ENCFF001NQP.fastq.gz",n=10^6)
fastq <- yield(fqSample)

## ----mycRep1ReadsShortReadQ,cache=TRUE,dependson="mycRep1Reads"----------
fastq

## ----mycRep1ReadsAccessor,cache=TRUE,dependson="mycRep1Reads"------------
readSequences <- sread(fastq)
readQuality <- quality(fastq)
readIDs <- id(fastq)
readSequences

## ----mycRep1ReadsQScores,cache=TRUE,dependson="mycRep1Reads"-------------
readQuality <- quality(fastq)
readQualities <- alphabetScore(readQuality)
readQualities[1:10]

## ----mycRep1ReadsQScoresPlot,cache=TRUE,dependson="mycRep1ReadsQScores",fig.height=3,fig.width=8----
library(ggplot2)
toPlot <- data.frame(ReadQ=readQualities)
ggplot(toPlot,aes(x=ReadQ))+geom_histogram()+theme_minimal()

## ----mycRep1ReadsAlpFreq,cache=TRUE,dependson="mycRep1Reads"-------------
readSequences <- sread(fastq)
readSequences_AlpFreq <- alphabetFrequency(readSequences)
readSequences_AlpFreq[1:3,]

## ----mycRep1ReadsAlpFreqSum,cache=TRUE,dependson="mycRep1ReadsAlpFreq"----
summed__AlpFreq  <- colSums(readSequences_AlpFreq)
summed__AlpFreq[c("A","C","G","T","N")]

## ----mycRep1ReadsAlpByCycle,cache=TRUE,dependson="mycRep1ReadsAlpFreq"----
readSequences_AlpbyCycle <- alphabetByCycle(readSequences)
readSequences_AlpbyCycle[1:4,1:10]

## ----mycRep1ReadsAlpByCyclePlot,cache=TRUE,dependson="mycRep1ReadsAlpFreq"----
AFreq <- readSequences_AlpbyCycle["A",]
CFreq <- readSequences_AlpbyCycle["C",]
GFreq <- readSequences_AlpbyCycle["G",]
TFreq <- readSequences_AlpbyCycle["T",]
toPlot <- data.frame(Count=c(AFreq,CFreq,GFreq,TFreq),
                     Cycle=rep(1:36,4),
                     Base=rep(c("A","C","G","T"),each=4))


## ----mycRep1ReadsAlpByCyclePlot2,cache=TRUE,eval=FALSE,dependson="mycRep1ReadsAlpByCyclePlot",fig.height=4,fig.width=8----
## 
## ggplot(toPlot,aes(y=Count,x=Cycle,colour=Base))+geom_line()+
##   theme_bw()

## ----mycRep1ReadsAlpByCyclePlot3,cache=TRUE,echo=FALSE,dependson="mycRep1ReadsAlpByCyclePlot",fig.height=4,fig.width=8----

ggplot(toPlot,aes(y=Count,x=Cycle,colour=Base))+geom_line()+ylim(150000,400000)+
  theme_bw()

## ----mycRep1ReadsQByCycle,cache=TRUE,dependson="mycRep1ReadsAlpFreq"-----
qualAsMatrix <- as(readQuality,"matrix")
qualAsMatrix[1:2,]

## ----mycRep1ReadsQByCyclePlot,cache=TRUE,dependson="mycRep1ReadsQByCycle",fig.width=8,fig.height=4----
boxplot(qualAsMatrix[1:1000,])

## ----out,eval=FALSE------------------------------------------------------
## fqStreamer <- FastqStreamer("~/Downloads/ENCFF001NQP.fastq.gz",
##                             n=100000)

## ----out1,eval=FALSE-----------------------------------------------------
## TotalReads <- 0
## TotalReadsFilt <- 0
## while (length(fq <- yield(fqStreamer))>0) {
##     TotalReads <- TotalReads+length(fq)
##     filt1 <- fq[alphabetScore(fq) > 300 ]
##     filt2 <- filt1[alphabetFrequency(sread(filt1))[,"N"] < 10]
##     TotalReadsFilt <- TotalReadsFilt+length(filt2)
##     writeFastq(filt2,"filtered_ENCFF001NQP.fastq.gz",mode="a")
## }
## TotalReads
## TotalReadsFilt

## ----fa1q, include=FALSE-------------------------------------------------
library(BSgenome.Mmusculus.UCSC.mm10)

## ----fa1, echo=TRUE------------------------------------------------------
library(BSgenome.Mmusculus.UCSC.mm10)
BSgenome.Mmusculus.UCSC.mm10

## ----fa2,cache=FALSE,echo=TRUE-------------------------------------------
mainChromosomes <- paste0("chr",c(1:19,"X","Y","M"))
mainChrSeq <- lapply(mainChromosomes,
                     function(x)BSgenome.Mmusculus.UCSC.mm10[[x]])
names(mainChrSeq) <- mainChromosomes
mainChrSeqSet <- DNAStringSet(mainChrSeq)
mainChrSeqSet

## ----fa3, echo=TRUE,eval=FALSE-------------------------------------------
## writeXStringSet(mainChrSeqSet,
##                 "BSgenome.Mmusculus.UCSC.mm10.mainChrs.fa")

## ----sampleTable1, echo=TRUE,eval=FALSE----------------------------------
## FileName <- "filtered_ENCFF001NQP.fastq.gz"
## SampleName <- "Myc_Mel_1"
## sampleTable <- data.frame(FileName,SampleName)
## write.table(sampleTable,file="sampleTable.txt",sep="\t",
##             quote=FALSE,row.names = FALSE)

## ----sampleTable1s, echo=FALSE,eval=TRUE---------------------------------
FileName <- "filtered_ENCFF001NQP.fastq.gz"
SampleName <- "Myc_Mel_1" 
data.frame(FileName,SampleName)

## ----bsgenome, echo=TRUE,eval=FALSE--------------------------------------
## library(QuasR)
## qAlign("sampleTable.txt","BSgenome.Mmusculus.UCSC.mm10.mainChrs.fa")

## ----index, echo=TRUE,eval=FALSE-----------------------------------------
## library(Rsubread)
## buildindex("BSgenome.Mmusculus.UCSC.mm10.mainChrs",
##            "BSgenome.Mmusculus.UCSC.mm10.mainChrs.fa")

## ----rsalign, echo=TRUE,eval=FALSE---------------------------------------
## align("BSgenome.Mmusculus.UCSC.mm10.mainChrs",
##       "filtered_ENCFF001NQP.fastq.gz",
##       output_file="filtered_ENCFF001NQP.bam",
##       type="dna",phredOffset = 64)

## ----sortindex, echo=TRUE,eval=FALSE-------------------------------------
## library(Rsamtools)
## sortBam("filtered_ENCFF001NQP.bam","SR_Myc_Mel_rep1")
## indexBam("SR_Myc_Mel_rep1.bam")

## ----mapped1, echo=TRUE,eval=FALSE---------------------------------------
## mappedReads <- idxstatsBam("../../../SR_Myc_Mel_rep1.bam")
## TotalMapped <- sum(mappedReads[,"mapped"])
## ggplot(mappedReads,aes(x=seqnames,y=mapped))+
##   geom_bar(stat="identity")

## ----mapped, echo=FALSE,eval=TRUE,fig.width=4,fig.height=4---------------
mappedReads <- idxstatsBam("../../../SR_Myc_Mel_rep1.bam")
TotalMapped <- sum(mappedReads[,"mapped"])
suppressPackageStartupMessages(library(ggplot2))
ggplot(mappedReads,aes(x=seqnames,y=mapped))+geom_bar(stat="identity")+coord_flip()

## ----coverage, echo=TRUE,eval=FALSE--------------------------------------
## forBigWig <- coverage("SR_Myc_Mel_rep1.bam")
## forBigWig

## ----bw, echo=TRUE,eval=FALSE--------------------------------------------
## library(rtracklayer)
## export.bw(forBigWig,con="SR_Myc_Mel_rep1.bw")

## ----weightedCover, echo=TRUE,eval=FALSE---------------------------------
## forBigWig <- coverage("SR_Myc_Mel_rep1.bam",
##                       weight = (10^6)/TotalMapped)
## forBigWig
## export.bw(forBigWig,con="SR_Myc_Mel_rep1_weighted.bw")

