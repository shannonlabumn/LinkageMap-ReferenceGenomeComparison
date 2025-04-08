library(vcfR)
library(updog)
library(dplyr)
library(tidyr)

DATA.DIR <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(DATA.DIR)
# Run it one at a time for each vcf dataset
vcf <- read.vcfR("NGSEPXJDY_merged204060AD.vcf")
vcf <- read.vcfR("NGSEPZm1-merged204060AD.vcf")
vcf <- read.vcfR("NGSEPPG-merged204060AD.vcf")

# Shown for ZhongmuNo.1 genome based dataset
extract.info(vcf, "CHROM")

ADmat <- extract.gt(vcf, "AD")
refmat <- strsplit(ADmat, ",")
refmat <- sapply(refmat, "[[", 1)
refmat <- matrix(refmat, nrow=nrow(ADmat))
refmat <- apply(refmat, 2, as.numeric)

DPmat <- extract.gt(vcf, "DP")
DPmat <- matrix(DPmat, nrow=nrow(ADmat))
DPmat <- apply(DPmat, 2, as.numeric)

altmat <- DPmat - refmat
altmat <- strsplit(ADmat, ",")
altmat <- sapply(altmat, "[[", 2)
altmat <- matrix(altmat, nrow=nrow(ADmat))
altmat <- apply(altmat, 2, as.numeric)

colnames(refmat) <- colnames(ADmat)
colnames(altmat) <- colnames(ADmat)
rownames(refmat) <- rownames(ADmat)
rownames(altmat) <- rownames(ADmat)
sizemat <- refmat + altmat

mout <- multidog(refmat=refmat,
                 sizemat = sizemat,
                 model="f1",
                 ploidy=4,
                 p1_id = "RSY_S170",
                 p2_id = "3988_S171",
                 nc=32)
dim(mout$snpdf)
mout_cleaned <- filter_snp(mout, prop_mis < 0.075 & bias > exp(-1.5) & bias < exp(1.5) )


dim(mout_cleaned$snpdf)
dim(mout_cleaned$inddf)

genomat1 <- format_multidog(mout_cleaned, varname = "geno")
head(genomat1)
dim(genomat1)
colnames(genomat1)
class(genomat1)

# Prepare files to import in Mappoly 
genomat1 <- as.data.frame(genomat1)
genomat1$P1 <- mout$snpdf$p1geno
genomat1$P2 <- mout$snpdf$p2geno
marks.1 <- tibble::rownames_to_column(genomat1, "snp_name")
marks.2 <- marks.1 %>% relocate("P1", .before = "102_S1") %>% relocate("P2", .before = "102_S1")
marks.3 <- marks.2 %>% separate(col = 1, into = c("sequence", "sequence_position"), remove = FALSE, sep = "_") %>% relocate("P1", .after = "snp_name") %>% relocate("P2", .after = "P1")
marks.3$sequence <- as.factor(marks.3$sequence)

# For XJDY
marks.3$sequence <- recode_factor(marks.3$sequence,
                             chr1.1 = "1", chr2.1 = "2", chr3.1 = "3",
                             chr4.1 = "4", chr5.1 = "5", chr6.1 = "6",
                             chr7.1 = "7", chr8.1 = "8")

# For Zm1
marks.3$sequence <- recode_factor(marks.3$sequence,
                                  Chr1 = "1", Chr2 = "2", Chr3 = "3",
                                  Chr4 = "4", Chr5 = "5", Chr6 = "6",
                                  Chr7 = "7", Chr8 = "8")
# For pangenome
marks.3$sequence <- recode_factor(marks.3$sequence,
                                  `ZhongmuNo.1:Chr1` = "1", `ZhongmuNo.1:Chr2` = "2", 
                                  `ZhongmuNo.1:Chr3` = "3", `ZhongmuNo.1:Chr4` = "4", 
                                  `ZhongmuNo.1:Chr5` = "5", `ZhongmuNo.1:Chr6` = "6",
                                  `ZhongmuNo.1:Chr7` = "7", `ZhongmuNo.1:Chr8` = "8")

# remove contigs
rem <- grep("contig", marks.3$snp_name)
marks.4 <- marks.3[-c(rem[1]:rem[length(rem)]),]

# export files ready to import in mappoly
write.csv(marks.4, file = "NGSEPXJDY_merged204060ADupdog.csv", quote = FALSE, row.names = FALSE)
write.csv(marks.4, file = "NGSEPZm1_merged204060ADupdog.csv", quote = FALSE, row.names = FALSE)
write.csv(marks.4, file = "NGSEPPG_merged204060ADupdog.csv", quote = FALSE, row.names = FALSE)

# Manually check order and position of markers in Pangenome file



