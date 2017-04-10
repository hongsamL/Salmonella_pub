#!/usr/bin/env Rscript

library(seqinr)
library(Biostrings)
library(stringr)
library(ape)
library(reshape2)

#--------------Pairwise SNPS CODE (without counting N)---------------------

num_snps <- function(seq1,seq2){
  #turn seq1 and seq2 into strings
  s1<-unlist(getSequence(seq1,as.string = TRUE))
  s2<-unlist(getSequence(seq2,as.string = TRUE))
  #locations of bases with "N" for seq1 and seq2
  n1<-str_locate_all(s1,'n')[[1]][,1]
  n2<-str_locate_all(s2,'n')[[1]][,1]
  # calculate non-overlapping "N" and remove from total pairwise distance
  outer <-union(n1,n2)
  inner <-intersect(n1,n2)
  remain <- setdiff(outer,inner)
  return(neditAt(s1,s2)-length(remain))
}

# input list of ids in fasta file "seqnames", fasta DNA alignment, and output file name "fh"
# writes a csv file with pairwise snp count
pairwise_snps <- function(seqnames,alignment,fh){
  write.table(matrix(c('seq1','seq2','numSNP'),ncol=3,byrow=TRUE),file=fh,sep=",",row.names = FALSE,col.names = FALSE)
  while(length(seqnames)>1){
    index<-seqnames[1]
    seqnames<-seqnames[-1]
    for(seq in seqnames) {
      data<-c(index,seq,num_snps(alignment[index],alignment[seq]))
      data <- matrix(data,ncol=3,byrow=TRUE)
      #data <- as.table(data)
      write.table(data,file=fh,append = TRUE,sep=",",row.names = FALSE,col.names = FALSE)
    }
  }
  return(c("Pairwise SNP Counts in file ", fh))
  }

#-------------------------------------------------------------------
args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)==1) {
  # default output file
  args[2] = paste("snps_",args[1],".csv",sep="")
}

aln <- read.fasta(args[1],as.string = TRUE)
# creates a vector with all ids in fasta file
seqs_names <-getName(aln)
#Calculate pairwise snps
pairwise_snps(seqs_names,aln,args[2])


