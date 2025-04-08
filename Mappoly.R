library(mappoly)
library(dplyr)
library(tidyr)
library(vcfR)

# Set the working directory
DATA.DIR <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(DATA.DIR)

# Read the data, this is the output of updog (run this script for all three output files from updog)
# This code is shown for ZhongmuNo.1 based dataset only
data <- read_geno_csv(file.in  = "NGSEPZm1_merged204060ADupdog.csv", ploidy = 4)

# Filter data
dat <- filter_missing(input.data = data, type = "marker", 
                      filter.thres = 0.05, inter = FALSE)
dat <- filter_missing(input.data = dat, type = "individual", 
                      filter.thres = 0.05, inter = FALSE)
pval.bonf <- 0.05/dat$n.mrk
mrks.chi.filt <- filter_segregation(dat, chisq.pval.thres = pval.bonf, inter = FALSE)

seq.init <- make_seq_mappoly(mrks.chi.filt)
seq.red  <- elim_redundant(seq.init)
seq.init <- make_seq_mappoly(seq.red)

png("seq.init.png")
plot(seq.init)
dev.off()

length(seq.init$seq.mrk.names)
length(seq.red$unique.seq$seq.mrk.names)
class(seq.red)
go <- get_genomic_order(input.seq = seq.init)
png("GenomicOrder.png")
plot(go)
dev.off()

# Two-point analysis
## local computation
n.cores = parallel::detectCores() - 1
all.rf.pairwise <- est_pairwise_rf(input.seq = seq.init,
                                   ncpus = n.cores, verbose=TRUE, memory.warning = TRUE)

mat <- rf_list_to_matrix(input.twopt = all.rf.pairwise)

go <- get_genomic_order(seq.init)
s.o <- make_seq_mappoly(go)
png("Rmatrix.png")
plot(mat,fact = 4)
dev.off()

# Grouping
g <- group_mappoly(input.mat = mat,
                   expected.groups = 8,
                   comp.mat = TRUE, 
                   inter=FALSE,
                   LODweight=TRUE,
                   verbose=TRUE)

g # Groups doesn't match the genomic information, using reference genome sequennce and position information to group and order markers

## Making groups only with genomic information (tends to bring more markers, but also discards scaffolds)
LGS.genomic = vector("list", 8)
for (i in 1:8){
  tempseq1 = make_seq_mappoly(dat, arg = paste0("seq",i), genomic.info = 1)
  mrks = intersect(tempseq1$seq.mrk.names, seq.init$seq.mrk.names)
  tempseq = make_seq_mappoly(dat, arg = mrks)
  temptpt = make_pairs_mappoly(all.rf.pairwise, input.seq = tempseq)
  rffilt = rf_snp_filter(input.twopt = temptpt, diagnostic.plot = FALSE)
  temptpt2 = make_pairs_mappoly(temptpt, input.seq = rffilt)
  LGS.genomic[[as.numeric(unique(rffilt$chrom))]] = list(seq = rffilt, tpt = temptpt2)
}


# Phasing----------------
LGS <- LGS.genomic
class(LGS[[1]]$tpt) 

MAPs <- vector("list", 8)
for (i in c(1:8)) {
  MAPs[[i]] <- est_rf_hmm_sequential(input.seq = LGS[[i]]$seq,
                                     start.set = 3,
                                     thres.twopt = 10,
                                     thres.hmm = 50,
                                     extend.tail = 30,
                                     twopt = LGS[[i]]$tpt,
                                     verbose = TRUE,
                                     phase.number.limit = 20,
                                     sub.map.size.diff.limit = 5)
}

MAPs_summary <- summary_maps(MAPs)
write.csv(MAPs_summary, file="SummaryMAPsNGSEPZm1_merged204060ADupdog.csv",
          row.names=FALSE, quote=FALSE)
          
png("MAPs.png")
plot_map_list(MAPs, horiz = TRUE, col = "ggstyle")
dev.off()

# Update Maps
MAPs.up <- vector("list", 8)
for(i in 1:8) {
  MAPs.up[[i]] <- est_full_hmm_with_global_error(MAPs[[i]], error = 0.05, verbose = TRUE)
}


MAPsUp_summary <- summary_maps(MAPs.up)
write.csv(MAPsUp_summary, file="SummaryMAPsUp_NGSEPZm1_merged204060ADupdog.csv",
          row.names=FALSE, quote=FALSE)

# Export the final map
export_map_list(MAPs.up, file="FinalMAPsUp_NGSEPZm1_merged204060ADupdog.csv")

# Manually check the genetic distance between marker pairs and 
# remove one of the SNP markers where distance is more than 10 cM
# Update the original csv file 

data <- read.csv("NGSEPZm1_merged204060ADupdog.csv", header=TRUE)
 
# Read the file with single column containing list of markers to be removed                 
Rmlist <- read.table("RmMarkers.txt",sep="\t", header=TRUE) 

data2 <- data[-which(data$snp_name %in% Rmlist$Marker.Name),]

write.csv(data2, "Updtd_NGSEPZm1_merged204060ADupdog.csv",
          quote=FALSE, row.names=FALSE)

# and run this script again from start usinng updated file.
# AAdding 
# Ref Alt allele information in the final map
vcf <- read.vcfR(file="NGSEPZm1_merged204060AD.vcf")
head(vcf)
chrom <- getFIX(vcf)[,1]

# For ZhongmuNo.1 genome
chrom <- paste("Chr", chrom, sep="")
position <- getFIX(vcf)[,2]
chrom_pos <- paste(chrom, position, sep="_")
Ref <- getFIX(vcf)[,4]
Alt <- getFIX(vcf)[,5]

MapFile <- read.csv("FinalMAPsUp_NGSEPZm1_merged204060ADupdog.csv",
                    header=TRUE)
MapFile

new.info <- data.frame(Marker.Name=chrom_pos, Ref.Allele=Ref, Alt.Allele=Alt)

# Uppercase Chr (may be needed for pangenome based file)
# MapFile$Marker.Name <- toupper(MapFile$Marker.Name) 


New.mapfile <- merge.data.frame(MapFile, new.info, by=c("Marker.Name"), all=FALSE)
New.mapfile <- New.mapfile[match(MapFile$Marker.Name, New.mapfile$Marker.Name),]

New.mapfile$Ref.Allele.x <- New.mapfile$Ref.Allele.y

New.mapfile$Alt.Allele.x <- New.mapfile$Alt.Allele.y
New.mapfile <- New.mapfile[,c(1:17)]

# Count number of dosages in Dose columns
New.mapfile %>% count(Dose.in.P1, sort = TRUE, name= "Dose.in.P1.freq")
New.mapfile %>% count(Dose.in.P2, sort = TRUE, name= "Dose.in.P2.freq")

write.csv(New.mapfile, "FinalMAPsUp_NGSEPZm1_merged204060ADupdog_RefAlt.csv")









