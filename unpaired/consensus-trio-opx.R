#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# test if there is one argument: if not, return an error
if (length(args) != 4) {
  stop("4 arguments must be provided, GATK file, SAMtools file, Mutect file and the base file name.n", call.=FALSE)
}

GATKfile <- args[1]
SAMfile <- args[2]
MutectFile <- args[3]
baseName <- args[4]

## Load Libraries
library(tidyverse)

print("Process GATK variants")
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
#G <- Filter(function(y)!all(y == "."), G);
G[G == ""] <- NA
G$MAGIC.GATK <- NULL
G <- G %>% dplyr::mutate_if(is.factor, as.character);
G <- G %>% mutate_if(is.integer, as.character);


print("Process SAMTOOLS variants")
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
#S <- Filter(function(y)!all(y == "."), S);
S[S == ""] <- NA
S$MAGIC.SAM <- NULL
S <- S %>% mutate_if(is.factor, as.character);
S <- S %>% mutate_if(is.integer, as.character);


print("Process Mutect variants")

Mu <- read.delim(file = MutectFile, row.names=NULL, 
                header = FALSE, skip = 1)
annHead <- read.delim(file = MutectFile, row.names = NULL, 
                      nrows = 1, header=FALSE, stringsAsFactors = FALSE)
extraHead <- c("col1", "DP.Mu", "CHR", "POS", "END", "REF", 
               "ALT", "col2", "FILTER.Mu", "INFO.Mu", "FORMAT.Mu", "MAGIC.Mu")
colnames(Mu) <- c(annHead[1,], extraHead)
Mu[,"col1"]<- NULL
Mu <- Mu %>% filter(Otherinfo %in% c("0", "0.5", "1"))
Mu$VariantID <- paste(Mu$CHR, Mu$Gene.refGene, Mu$POS, Mu$REF, Mu$ALT, sep = "-")
Mu$DP.Mu <- as.integer(Mu$DP.Mu)
Mu$VariantID <- paste(Mu$CHR, Mu$Gene.refGene, Mu$POS, Mu$REF, Mu$ALT, sep = "-")
Mu$AD.Mu <- as.integer(gsub("^[^:]+:[^,]+,|:.*", "", Mu$MAGIC.Mu));
Mu$VAF.Mu <- 100*as.numeric(Mu$AD.Mu)/as.numeric(Mu$DP.Mu);
#Mu <- Filter(function(y)!all(y == "."), Mu);
Mu[Mu == ""] <- NA
Mu$MAGIC.Mu <- NULL
Mu <- Mu %>% dplyr::mutate_if(is.factor, as.character);
Mu <- Mu %>% mutate_if(is.integer, as.character);

## Deal with flubbed annotation issues where the same variant is now annotated two different ways!  ARGH!
#commonCols <- intersect(intersect(colnames(G),colnames(S)),colnames(Mu))
commonCols <- c("CHR", "POS", "REF", "ALT")
GMerge <- G %>% select(c(commonCols), AD.GATK, DP.GATK)
SMerge <- S %>% select(c(commonCols), AD.SAM, DP.SAM)
MuMerge <- Mu %>% select(c(commonCols), FILTER.Mu,AD.Mu, DP.Mu)

variants <- full_join(full_join(GMerge, SMerge), MuMerge)
variants$Type <- ifelse(nchar(variants$REF) == nchar(variants$ALT), "SNV", "INDEL");


filtered <- variants %>% 
  filter(!ExonicFunc.refGene %in% c("synonymous SNV", ".", "unknown")) %>%
  filter(FILTER.Mu == "PASS" & is.na(AD.SAM)==F & is.na(AD.GATK)==F |  # all three callers
           FILTER.Mu == "PASS" & is.na(AD.GATK)==F |  # pass for Mutect, and called by GATK too - NPM1 falls here
           is.na(AD.SAM)==F  & is.na(AD.GATK)==F | # regardless of mutect, if called by both others
           grepl("germline", FILTER.Mu) ==T  & is.na(AD.GATK)==F | # GATK, and germline, any samtools, this includes DNMT3A R882, FLT3 insertions and many NPM1 insertions
           grepl("slippage", FILTER.Mu) ==T  & is.na(AD.GATK)==F | # GATK, and contains slippage filter for Mutect, any samtools - FLT3 ITD gets lost here too
           FILTER.Mu == "clustered_events" & is.na(AD.SAM)==F & is.na(AD.GATK)==F # Mutect clustered events, and both other callers
  )  %>% filter(grepl("base_qual", FILTER.Mu) == F & grepl("map_qual", FILTER.Mu) == F) # regardless of other reasons, remove things that mutect filtered due to base or mapping quality concerns

filtered$ExAC_ALL[grepl("^\\.$",filtered$ExAC_ALL)] <- 0
filtered$ExAC_ALL <- as.numeric(filtered$ExAC_ALL)
filtered <- filtered %>% filter(ExAC_ALL < 0.01 | cosmic70 != ".")  %>% select(-ExAC_ALL)
# UWLM filter:
#.03 minimum variant frequency and a minimum of 5 variant reads for SNVs and 0.01 minimum variant frequency and a minimum of 4 variant reads for indels.
filtered <- filtered %>% filter(AD.SAM >= 4 | AD.Mu >= 4 | AD.GATK >= 4) %>% unique()


# take all the common annotation columns from the GATK and Strelka data and make a giant full join that has all the annotations from any variant called in either dataset
sharedCols <- intersect(intersect(colnames(G),colnames(S)),colnames(Mu))
annotations <- full_join(G %>% select(all_of(sharedCols)), S %>% select(all_of(sharedCols))) %>% unique()

# now reapply those annotations to the filtered variants list, with the exception of AF_popmax since that's been reformatted in the variants list. 
reannotate <- left_join(filtered, annotations)
reannotate$molecular_id <-baseName

write.table(reannotate, file = paste0(baseName, ".consensus.filtered.tsv"), sep = "\t",
            row.names = F)

