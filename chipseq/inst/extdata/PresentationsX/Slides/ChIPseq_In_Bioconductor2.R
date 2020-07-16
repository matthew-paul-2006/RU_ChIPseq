## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
AsSlides <- TRUE

## ----mycQCdwdwshowL,include=FALSE----------------------------------------
library(ChIPQC)

## ----eval=F--------------------------------------------------------------
## QCresult <- ChIPQCsample(reads="/pathTo/myChIPreads.bam",
##                          genome="mm10",
##                          blacklist = "/pathTo/mm10_Blacklist.bed")

## ----mycQC,cache=TRUE,eval=FALSE-----------------------------------------
## library(ChIPQC)
## chipqc_MycMel_rep1 <- ChIPQCsample("SR_Myc_Mel_rep1.bam",
##                          annotation = "mm10",
##                          blacklist = toBlkList,
##                          chromosomes = paste0(c("chr",1:10)))
## class(chipqc_MycMel_rep1)

## ----mycQCshowLa,echo=FALSE,eval=TRUE------------------------------------
library(ChIPQC)
load(file="rep1.RData")
class(chipqc_MycMel_rep1)

## ----mycQCshowL,include=FALSE--------------------------------------------
library(ChIPQC)
load(file="rep1.RData")

## ----mycQCshow,eval=TRUE-------------------------------------------------
chipqc_MycMel_rep1

## ----mycQCshowd2,cache=TRUE,eval=FALSE,include=FALSE---------------------
## setwd("~/Projects/Results/chipseq/testRun/BAMs/")
## bamsToQC <- c("Sorted_Myc_Ch12_1.bam","Sorted_Myc_Ch12_2.bam",
##              "Sorted_Myc_MEL_1.bam","Sorted_Myc_MEL_2.bam",
##              "Sorted_Input_MEL.bam","Sorted_Input_Ch12.bam")
## myQC <- bplapply(bamsToQC,ChIPQCsample,
##         annotation = "mm10",
##         blacklist = toBlkList,
##         chromosomes = paste0("chr",1:10))
## names(myQC) <- bamsToQC

## ----mycQCshow2,cache=TRUE,eval=FALSE------------------------------------
## bamsToQC <- c("Sorted_Myc_Ch12_1.bam","Sorted_Myc_Ch12_2.bam",
##              "Sorted_Myc_MEL_1.bam","Sorted_Myc_MEL_2.bam",
##              "Sorted_Input_MEL.bam","Sorted_Input_Ch12.bam")
## myQC <- bplapply(bamsToQC,ChIPQCsample,
##         annotation = "mm10",
##         blacklist = toBlkList,
##         chromosomes = paste0("chr",1:10))
## names(myQC) <- bamsToQC

## ----qcmetricsA,include=FALSE--------------------------------------------
load(file="myQCnoPeaks.RData")

## ----qcmetrics,cache=FALSE,eval=TRUE-------------------------------------
QCmetrics(myQC)

## ----qcmetridedecs,cache=FALSE,eval=TRUE,fig.width=6,fig.height=4--------
plotCC(myQC,facetBy = "Sample")

## ----qcmetridecs,cache=FALSE,eval=TRUE-----------------------------------
myMeta <- data.frame(Sample= names(myQC),
                     Tissue=c("Ch12","Ch12","MEL","MEL","MEL","Ch12"),
                     Antibody=c(rep("Myc",4),rep("Input",2)))
myMeta

## ----qcmetricsede,cache=FALSE,eval=TRUE,fig.width=6,fig.height=3---------
plotCC(myQC,facetBy = "Tissue",addMetaData = myMeta,
       colourBy="Antibody")

## ----qcmetricsrf,cache=FALSE,eval=TRUE,fig.width=6,fig.height=3----------
plotCC(myQC,facetBy = "Tissue",addMetaData = myMeta,
       colourBy="Antibody")+theme_bw()+
  ggtitle("ChIPQC results")

## ----fig.width=6,fig.height=2,warning=FALSE,message=FALSE----------------
plotSSD(myQC)+xlim(0,5)

## ----fig.width=6,fig.height=3,warning=FALSE,message=FALSE----------------
plotSSD(myQC)+xlim(0.2,0.4)

## ----fig.width=6,fig.height=3,warning=FALSE,message=FALSE----------------
myFlags <- flagtagcounts(myQC)
myFlags["DuplicateByChIPQC",]/myFlags["Mapped",]

## ----warning=FALSE,message=FALSE,fig.width=8,fig.height=4----------------
p <- plotRegi(myQC)

## ----warning=FALSE,fig.width=12,fig.height=6-----------------------------
p

## macs2 callpeak -t Sorted_Myc_MEL_1.bam

## ----fig.height=5, fig.width=15,eval=FALSE-------------------------------
## myChIP <- "Sorted_Myc_MEL_1.bam"
## myControl <- "Sorted_Input_MEL.bam"
## 
## macsCommand <- paste0("macs2 callpeak -t ", myChIP,
##                       " -n ", "Mel_Rep1",
##                       " –-outdir ","PeakDirectory",
##                       " -c ", myControl)
## system(macsCommand)

## ----eval=T,echo=F,  warning=FALSE,collapse=T----------------------------
macsPeaks <- "../Slides/data/Mel1_peaks.xls"
macsPeaks_DF <- read.delim(macsPeaks,sep="",comment="#")
macsPeaks_DF[1:3,]

## ----eval=T,echo=F,  warning=FALSE,collapse=T----------------------------
macsPeaks_DF <- read.delim(macsPeaks,comment="",stringsAsFactors = F)
macsPeaks_DF[1:8,]

## ----eval=T,echo=T, eval=T, echo=T, warning=FALSE------------------------
macsPeaks_DF <- read.delim(macsPeaks,comment.char="#")
macsPeaks_DF[1:2,]

## ----eval=T,echo=T, eval=T, echo=T, warning=FALSE------------------------
library(GenomicRanges)
macsPeaks_GR <- GRanges(
 seqnames=macsPeaks_DF[,"chr"],
 IRanges(macsPeaks_DF[,"start"],
         macsPeaks_DF[,"end"]
 )
)
macsPeaks_GR

## ----eval=T,echo=T, eval=T, echo=T, warning=FALSE------------------------
seqnames(macsPeaks_GR)
ranges(macsPeaks_GR)

## ----eval=T,echo=T, eval=T, echo=T, warning=FALSE------------------------
mcols(macsPeaks_GR) <- macsPeaks_DF[,c("abs_summit", "fold_enrichment")]
macsPeaks_GR

## ----eval=T,echo=T, eval=T, echo=T, warning=FALSE------------------------
library(rtracklayer)
blkList <- import.bed("mm10-blacklist.bed")
macsPeaks_GR <- macsPeaks_GR[!macsPeaks_GR %over% blkList] 

## ----include=FALSE-------------------------------------------------------
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(GenomeInfoDb)
library(ChIPseeker)


## ----eval=F,echo=T, eval=T, echo=T, warning=FALSE,tidy=T,message=FALSE----
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(GenomeInfoDb)
library(ChIPseeker)


## ----eval=T,echo=T, message=FALSE,messages=FALSE, eval=T, echo=T, warning=FALSE----
peakAnno <- annotatePeak(macsPeaks_GR, tssRegion=c(-500, 500), 
                         TxDb=TxDb.Mmusculus.UCSC.mm10.knownGene, 
                         annoDb="org.Mm.eg.db")
class(peakAnno)

## ----eval=T,echo=T, message=F,messages=F, eval=T, echo=T, warning=FALSE,tidy=T----
peakAnno

## ----eval=T,echo=T, eval=T, echo=T, warning=FALSE,tidy=T-----------------
peakAnno_GR <- as.GRanges(peakAnno)
peakAnno_DF <- as.data.frame(peakAnno)

## ----eval=T,echo=T, eval=T, echo=T, warning=FALSE,tidy=T-----------------
peakAnno_GR[1,]

## ---- eval=T, echo=T, fig.height=5, fig.width=15, warning=FALSE, tidy=T----
plotAnnoBar(peakAnno)

## ----eval=T,echo=T, eval=F, echo=T, warning=FALSE,fig.height=5, fig.width=15,tidy=T----
## plotDistToTSS(peakAnno)

## ---- eval=T, echo=T, fig.height=5, fig.width=15, warning=FALSE, tidy=T----
upsetplot(peakAnno, vennpie=F)

