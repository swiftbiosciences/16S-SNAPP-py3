#!/usr/bin/env Rscript
## Swift Biosciences 16S snapp workflow
## Author Benli Chai & Sukhinder Sandhu 20200502

## run dada2 for 16S amplicon sequences
#Supply the arguments from the command line

args <- commandArgs(TRUE)
if (length(args) != 3){
    stop("Usage: run_dada2.R inputDir workDir trimLength");
}

library(dada2)
path <- normalizePath(args[1]) #inputDir
wd <- normalizePath(args[2]) #workdir
READLEN <- as.numeric(args[3]) #the expected read length after the primer is trimmed
log <- args[4] #the log file to append process stats

fnFs <- sort(list.files(path, pattern="_R1_001.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq.gz", full.names = TRUE))

sample.names <- sapply(strsplit(basename(fnFs), "_L001"), `[`, 1)

filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(READLEN, READLEN), maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE, compress=TRUE, multithread=TRUE)

errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

#Sample inference with standard all sample pooling
dadaFs <- dada(filtFs, err=errF, multithread=TRUE, pool=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE, pool=TRUE)

mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE, justConcatenate=TRUE)

#get asv sequence and count table and write to a tab-delimited file
seqtab <- makeSequenceTable(mergers)
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
write.csv(t(seqtab.nochim), file=paste(wd, 'asv_seqNcount.csv', sep='/'))

#Get process stats and write to a tab-delimited file
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("Primer_trimmed", "Filtered", "Denoised_FWD", "Denoised_REV", "Merged", "Non_chimera")
rownames(track) <- sample.names
write.table(track, file=paste(wd, 'DADA2_summary.csv', sep="/"), sep=" ") #save processing stats to a file
