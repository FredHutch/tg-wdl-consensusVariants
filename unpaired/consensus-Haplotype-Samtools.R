#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# test if there is one argument: if not, return an error
if (length(args) != 3) {
  stop("3 arguments must be provided, GATK file, SAMtools file and base file name.n", call.=FALSE)
}

GATKfile <- args[1]
SAMfile <- args[2]
baseName <- args[3]

## Load Libraries
library(dplyr)

print("Pull down and process GATK variants")
G <- read.delim(file = GATKfile, row.names = NULL, 
                  header=FALSE, stringsAsFactors = FALSE, skip = 1)
annHead <- read.delim(file = GATKfile, row.names = NULL, 
                        nrows = 1, header=FALSE, stringsAsFactors = FALSE)
extraHead <- c("QUAL1", "DP.GATK", "CHR", "POS", "END", "REF", 
               "ALT", "QUAL.GATK", "FILTER", "INFO.GATK", "FORMAT.GATK", "MAGIC.GATK")
colnames(G) <- c(annHead[1,], extraHead)
G[, "QUAL1"]<- NULL
G <- G %>% filter(Otherinfo %in% c("0", "0.5", "1"))
G$DP.GATK <- as.integer(G$DP.GATK)
G$VariantID <- paste(G$CHR, G$Gene.refGene, G$POS, G$REF, G$ALT, sep = "-")
G$AD.GATK <- as.integer(gsub("^[^:]+:[^,]+,|:.*", "", G$MAGIC.GATK));
G$VAF.GATK <- 100*as.numeric(G$AD.GATK)/as.numeric(G$DP.GATK);
G <- Filter(function(y)!all(y == "."), G);
G[G == ""] <- NA
G <- G %>% dplyr::mutate_if(is.factor, as.character);
G <- G %>% mutate_if(is.integer, as.character);


print("Pull down and process SAMTOOLS variants")
S <- read.delim(file = SAMfile, row.names=NULL, 
                  header = FALSE, skip = 1)
annHead <- read.delim(file = SAMfile, row.names = NULL, 
                        nrows = 1, header=FALSE, stringsAsFactors = FALSE)
extraHead <- c("QUAL1", "DP.SAM", "CHR", "POS", "END", "REF", 
               "ALT", "QUAL.SAM", "FILTER", "INFO.SAM", "FORMAT.SAM", "MAGIC.SAM")
colnames(S) <- c(annHead[1,], extraHead)
#DP4 = Number of 1) forward ref alleles; 2) reverse ref; 3) forward non-ref; 4) reverse non-ref alleles,
# used in variant calling. Sum can be smaller than DP because low-quality bases are not counted.
S[,"QUAL1"]<- NULL
S <- S %>% filter(Otherinfo %in% c("0", "0.5", "1"))
S$VariantID <- paste(S$CHR, S$Gene.refGene, S$POS, S$REF, S$ALT, sep = "-")
a <- gsub(".*:", "",S$MAGIC.SAM);
REF <- as.numeric(gsub(",.*$","",a));
S$AD.SAM <- as.numeric(gsub("^.*,","",a));
S$DP.SAM <- REF+S$AD.SAM;
S$VAF.SAM <- 100*S$AD.SAM/S$DP.SAM;
S <- Filter(function(y)!all(y == "."), S);
S[S == ""] <- NA
S <- S %>% mutate_if(is.factor, as.character);
S <- S %>% mutate_if(is.integer, as.character);

## Deal with flubbed annotation issues where the same variant is now annotated two different ways!  ARGH!
commonCols <- c("VariantID","CHR", "POS", "REF", "ALT", "Chr", "Start", "End", "Ref", "Alt")
GATKCols <- colnames(G)[!colnames(G) %in% colnames(S)]
SAMCols <- colnames(S)[!colnames(S) %in% colnames(G)]

GMerge <- G %>% select(c(commonCols, GATKCols))
SMerge <- S %>% select(c(commonCols, SAMCols))

variants <- full_join(GMerge, SMerge)
variants$Type <- ifelse(nchar(variants$REF) == nchar(variants$ALT), "SNV", "INDEL");

reannotate <- left_join(variants, G)
reannotate <- left_join(reannotate, S)

reannotate$VAF.DIFF <- reannotate$VAF.GATK-reannotate$VAF.SAM;
reannotate$Confidence  <- ifelse(!is.na(reannotate$AD.GATK) & !is.na(reannotate$AD.SAM), "conftier1",
                  ifelse(is.na(reannotate$AD.GATK) & !is.na(reannotate$AD.SAM) & reannotate$Type == "SNV", "conftier2",
                         ifelse(!is.na(reannotate$AD.GATK) & is.na(reannotate$AD.SAM) & reannotate$Type == "INDEL", "conftier2", 
                                "conftier3")))
write.table(reannotate, file = paste0(baseName, ".consensus.tsv"), sep = "\t",
            row.names = F)

