#!/bin/bash -l  

# check the quality of raw reads
module load fastqc/0.11.9
for i in $PWD/*R1_001.fastq.gz
do
fastqc -o FastQC ${i}
done

# clean the reads and drop low quallity reads
module load trimmomatic/0.39
for infile in $PWD/*_R1_001.fastq.gz
do
outfile=$infile\_trim.fastq.gz
java -jar Trimmomatic-0.39/trimmomatic-0.39.jar SE -phred33 $infile $outfile ILLUMINACLIP:TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
done

# Check the quality again after trimming the reads

# Align the reads to reference genome, we used 3 different reference genome
# Subsequent steps were performed simultanoeuly for 3 genomes

# generate the index for the genomes
module load bwa/0.7.17

bwa index final_homo1_genome.fasta
bwa index ZhongmuNo.1_genome.fasta
bwa index Zm1Zm4XinJiangDaYeRSYCADL.fasta

for infile in $PWD/*R1_001.fastq.gz_trim.fastq.gz
do
bwa mem Zm1Zm4XinJiangDaYeRSYCADL.fasta $infile -M  > $infile.PG.bam
bwa mem final_homo1_genome.fasta $infile -M  > $infile.XJDY.bam
bwa mem ZhongmuNo.1_genome.fasta $infile -M  > $infile.Zm1.bam
done

# Sort the reads
module load picard/2.25.6
for i in $PWD/*.bam
do
java -jar /common/software/install/migrated/picard/2.25.6/picard.jar SortSam MAX_RECORDS_IN_RAM=1000000 SO=coordinate CREATE_INDEX=true I=${i} O=${i}_sorted.bam 
done

# Add read groups for all three datsets
module load picard/2.25.6
for i in $(ls | grep ".bam$" | cut -d"_" -f1,2 | uniq)
do
java -jar $PICARD_DIR/picard.jar AddOrReplaceReadGroups \
I=${i}_R1_001.fastq.gz_trim.fastq.gz.PG.bam_sorted.bam \
O=${i}_R1_001_PG_sortedRG.bam \
RGID=${i} \
RGPL=illumina \
RGLB=lib \
RGPU=barcode \
RGSM=${i}
done

for i in $(ls | grep ".bam$" | cut -d"_" -f1,2 | uniq)
do
java -jar $PICARD_DIR/picard.jar AddOrReplaceReadGroups \
I=${i}_R1_001.fastq.gz_trim.fastq.gz.XJDY.bam_sorted.bam \
O=${i}_R1_001_XJDY_sortedRG.bam \
RGID=${i} \
RGPL=illumina \
RGLB=lib \
RGPU=barcode \
RGSM=${i}
done

for i in $(ls | grep ".bam$" | cut -d"_" -f1,2 | uniq)
do
java -jar $PICARD_DIR/picard.jar AddOrReplaceReadGroups \
I=${i}_R1_001.fastq.gz_trim.fastq.gz.Zm1.bam_sorted.bam \
O=${i}_R1_001_Zm1_sortedRG.bam \
RGID=${i} \
RGPL=illumina \
RGLB=lib \
RGPU=barcode \
RGSM=${i}
done

# Call the variants (example is shown for Zm1 dataset only, repeat it for PG and XJDY reference genomes based bam files as well 
# to get NGSEP_bwaZm1.vcf, NGSEP_bwaPG.vcf and NGSEP_bwaXJDY.vcf files
java -jar NGSEPcore_4.2.1/NGSEPcore.jar MultisampleVariantsDetector -maxAlnsPerStartPos 100 -maxBaseQS 30 -ploidy 4 -psp  -r ZhongmuNo.1_genome.fasta -o NGSEP_bwaZm1.vcf 102_S1_R1_001.fastq.gz_trim.fastq.gz.Zm1.Zm1.bam_sortedRG.bam 104_S2_R1_001.fastq.gz_trim.fastq.gz.Zm1.bam_sortedRG.bam 111_S3_R1_001.fastq.gz_trim.fastq.gz.Zm1.bam_sortedRG.bam 112_S4_R1_001.fastq.gz_trim.fastq.gz.Zm1.bam_sortedRG.bam 113_S5_R1_001.fastq.gz_trim.fastq.gz.Zm1.bam_sortedRG.bam 114_S6_R1_001.fastq.gz_trim.fastq.gz.Zm1.bam_sortedRG.bam 116_S7_R1_001.fastq.gz_trim.fastq.gz.Zm1.bam_sortedRG.bam 117_S8_R1_001.fastq.gz_trim.fastq.gz.Zm1.bam_sortedRG.bam 118_S9_R1_001.fastq.gz_trim.fastq.gz.Zm1.bam_sortedRG.bam 119_S10_R1_001.fastq.gz_trim.fastq.gz.Zm1.bam_sortedRG.bam 11_S11_R1_001.fastq.gz_trim.fastq.gz.Zm1.bam_sortedRG.bam 120_S12_R1_001.fastq.gz_trim.fastq.gz.Zm1.bam_sortedRG.bam 121_S13_R1_001.fastq.gz_trim.fastq.gz.Zm1.bam_sortedRG.bam 122_S14_R1_001.fastq.gz_trim.fastq.gz.Zm1.bam_sortedRG.bam 124_S15_R1_001.fastq.gz_trim.fastq.gz.Zm1.bam_sortedRG.bam 125_S16_R1_001.fastq.gz_trim.fastq.gz.Zm1.bam_sortedRG.bam 126_S17_R1_001.fastq.gz_trim.fastq.gz.Zm1.bam_sortedRG.bam 130_S18_R1_001.fastq.gz_trim.fastq.gz.Zm1.bam_sortedRG.bam 131_S19_R1_001.fastq.gz_trim.fastq.gz.Zm1.bam_sortedRG.bam 132_S20_R1_001.fastq.gz_trim.fastq.gz.Zm1.bam_sortedRG.bam 135_S21_R1_001.fastq.gz_trim.fastq.gz.Zm1.bam_sortedRG.bam 136_S22_R1_001.fastq.gz_trim.fastq.gz.Zm1.bam_sortedRG.bam 138_S23_R1_001.fastq.gz_trim.fastq.gz.Zm1.bam_sortedRG.bam 139_S24_R1_001.fastq.gz_trim.fastq.gz.Zm1.bam_sortedRG.bam 140_S25_R1_001.fastq.gz_trim.fastq.gz.Zm1.bam_sortedRG.bam 143_S26_R1_001.fastq.gz_trim.fastq.gz.Zm1.bam_sortedRG.bam 145_S27_R1_001.fastq.gz_trim.fastq.gz.Zm1.bam_sortedRG.bam 146_S28_R1_001.fastq.gz_trim.fastq.gz.Zm1.bam_sortedRG.bam 147_S29_R1_001.fastq.gz_trim.fastq.gz.Zm1.bam_sortedRG.bam 148_S30_R1_001.fastq.gz_trim.fastq.gz.Zm1.bam_sortedRG.bam 149_S31_R1_001.fastq.gz_trim.fastq.gz.Zm1.bam_sortedRG.bam 14_S32_R1_001.fastq.gz_trim.fastq.gz.Zm1.bam_sortedRG.bam 150_S33_R1_001.fastq.gz_trim.fastq.gz.Zm1.bam_sortedRG.bam 152_S34_R1_001.fastq.gz_trim.fastq.gz.Zm1.bam_sortedRG.bam 153_S35_R1_001.fastq.gz_trim.fastq.gz.Zm1.bam_sortedRG.bam 154_S36_R1_001.fastq.gz_trim.fastq.gz.Zm1.bam_sortedRG.bam 155_S37_R1_001.fastq.gz_trim.fastq.gz.Zm1.bam_sortedRG.bam 156_S38_R1_001.fastq.gz_trim.fastq.gz.Zm1.bam_sortedRG.bam 157_S39_R1_001.fastq.gz_trim.fastq.gz.Zm1.bam_sortedRG.bam 158_S40_R1_001.fastq.gz_trim.fastq.gz.Zm1.bam_sortedRG.bam 159_S41_R1_001.fastq.gz_trim.fastq.gz.Zm1.bam_sortedRG.bam 15_S42_R1_001.fastq.gz_trim.fastq.gz.Zm1.bam_sortedRG.bam 160_S43_R1_001.fastq.gz_trim.fastq.gz.Zm1.bam_sortedRG.bam 161_S44_R1_001.fastq.gz_trim.fastq.gz.Zm1.bam_sortedRG.bam 162_S45_R1_001.fastq.gz_trim.fastq.gz.Zm1.bam_sortedRG.bam 163_S46_R1_001.fastq.gz_trim.fastq.gz.Zm1.bam_sortedRG.bam 164_S47_R1_001.fastq.gz_trim.fastq.gz.Zm1.bam_sortedRG.bam 165_S48_R1_001.fastq.gz_trim.fastq.gz.Zm1.bam_sortedRG.bam 166_S49_R1_001.fastq.gz_trim.fastq.gz.Zm1.bam_sortedRG.bam 168_S50_R1_001.fastq.gz_trim.fastq.gz.Zm1.bam_sortedRG.bam 169_S51_R1_001.fastq.gz_trim.fastq.gz.Zm1.bam_sortedRG.bam 170_S52_R1_001.fastq.gz_trim.fastq.gz.Zm1.bam_sortedRG.bam 171_S53_R1_001.fastq.gz_trim.fastq.gz.Zm1.bam_sortedRG.bam 172_S54_R1_001.fastq.gz_trim.fastq.gz.Zm1.bam_sortedRG.bam 173_S55_R1_001.fastq.gz_trim.fastq.gz.Zm1.bam_sortedRG.bam 176_S57_R1_001.fastq.gz_trim.fastq.gz.Zm1.bam_sortedRG.bam 178_S58_R1_001.fastq.gz_trim.fastq.gz.Zm1.bam_sortedRG.bam 179_S59_R1_001.fastq.gz_trim.fastq.gz.Zm1.bam_sortedRG.bam 181_S60_R1_001.fastq.gz_trim.fastq.gz.Zm1.bam_sortedRG.bam 182_S61_R1_001.fastq.gz_trim.fastq.gz.Zm1.bam_sortedRG.bam 183_S62_R1_001.fastq.gz_trim.fastq.gz.Zm1.bam_sortedRG.bam 184_S63_R1_001.fastq.gz_trim.fastq.gz.Zm1.bam_sortedRG.bam 185_S64_R1_001.fastq.gz_trim.fastq.gz.Zm1.bam_sortedRG.bam 186_S65_R1_001.fastq.gz_trim.fastq.gz.Zm1.bam_sortedRG.bam 187_S66_R1_001.fastq.gz_trim.fastq.gz.Zm1.bam_sortedRG.bam 188_S67_R1_001.fastq.gz_trim.fastq.gz.Zm1.bam_sortedRG.bam 190_S68_R1_001.fastq.gz_trim.fastq.gz.Zm1.bam_sortedRG.bam 191_S69_R1_001.fastq.gz_trim.fastq.gz.Zm1.bam_sortedRG.bam 192_S70_R1_001.fastq.gz_trim.fastq.gz.Zm1.bam_sortedRG.bam 193_S71_R1_001.fastq.gz_trim.fastq.gz.Zm1.bam_sortedRG.bam 194_S72_R1_001.fastq.gz_trim.fastq.gz.Zm1.bam_sortedRG.bam 195_S73_R1_001.fastq.gz_trim.fastq.gz.Zm1.bam_sortedRG.bam 197_S74_R1_001.fastq.gz_trim.fastq.gz.Zm1.bam_sortedRG.bam 198_S75_R1_001.fastq.gz_trim.fastq.gz.Zm1.bam_sortedRG.bam 200_S76_R1_001.fastq.gz_trim.fastq.gz.Zm1.bam_sortedRG.bam 201_S77_R1_001.fastq.gz_trim.fastq.gz.Zm1.bam_sortedRG.bam 202_S78_R1_001.fastq.gz_trim.fastq.gz.Zm1.bam_sortedRG.bam 203_S79_R1_001.fastq.gz_trim.fastq.gz.Zm1.bam_sortedRG.bam 204_S80_R1_001.fastq.gz_trim.fastq.gz.Zm1.bam_sortedRG.bam 205_S81_R1_001.fastq.gz_trim.fastq.gz.Zm1.bam_sortedRG.bam 206_S82_R1_001.fastq.gz_trim.fastq.gz.Zm1.bam_sortedRG.bam 208_S83_R1_001.fastq.gz_trim.fastq.gz.Zm1.bam_sortedRG.bam 210_S84_R1_001.fastq.gz_trim.fastq.gz.Zm1.bam_sortedRG.bam 211_S85_R1_001.fastq.gz_trim.fastq.gz.Zm1.bam_sortedRG.bam 212_S86_R1_001.fastq.gz_trim.fastq.gz.Zm1.bam_sortedRG.bam 213_S87_R1_001.fastq.gz_trim.fastq.gz.Zm1.bam_sortedRG.bam 214_S88_R1_001.fastq.gz_trim.fastq.gz.Zm1.bam_sortedRG.bam 216_S89_R1_001.fastq.gz_trim.fastq.gz.Zm1.bam_sortedRG.bam 217_S90_R1_001.fastq.gz_trim.fastq.gz.Zm1.bam_sortedRG.bam 218_S91_R1_001.fastq.gz_trim.fastq.gz.Zm1.bam_sortedRG.bam 219_S92_R1_001.fastq.gz_trim.fastq.gz.Zm1.bam_sortedRG.bam 220_S93_R1_001.fastq.gz_trim.fastq.gz.Zm1.bam_sortedRG.bam 223_S94_R1_001.fastq.gz_trim.fastq.gz.Zm1.bam_sortedRG.bam 226_S95_R1_001.fastq.gz_trim.fastq.gz.Zm1.bam_sortedRG.bam 228_S96_R1_001.fastq.gz_trim.fastq.gz.Zm1.bam_sortedRG.bam 229_S97_R1_001.fastq.gz_trim.fastq.gz.Zm1.bam_sortedRG.bam 231_S98_R1_001.fastq.gz_trim.fastq.gz.Zm1.bam_sortedRG.bam 232_S99_R1_001.fastq.gz_trim.fastq.gz.Zm1.bam_sortedRG.bam 236_S100_R1_001.fastq.gz_trim.fastq.gz.Zm1.bam_sortedRG.bam 237_S101_R1_001.fastq.gz_trim.fastq.gz.Zm1.bam_sortedRG.bam 238_S102_R1_001.fastq.gz_trim.fastq.gz.Zm1.bam_sortedRG.bam 239_S103_R1_001.fastq.gz_trim.fastq.gz.Zm1.bam_sortedRG.bam 241_S104_R1_001.fastq.gz_trim.fastq.gz.Zm1.bam_sortedRG.bam 242_S105_R1_001.fastq.gz_trim.fastq.gz.Zm1.bam_sortedRG.bam 243_S106_R1_001.fastq.gz_trim.fastq.gz.Zm1.bam_sortedRG.bam 244_S107_R1_001.fastq.gz_trim.fastq.gz.Zm1.bam_sortedRG.bam 248_S108_R1_001.fastq.gz_trim.fastq.gz.Zm1.bam_sortedRG.bam 249_S109_R1_001.fastq.gz_trim.fastq.gz.Zm1.bam_sortedRG.bam 250_S110_R1_001.fastq.gz_trim.fastq.gz.Zm1.bam_sortedRG.bam 251_S111_R1_001.fastq.gz_trim.fastq.gz.Zm1.bam_sortedRG.bam 252_S112_R1_001.fastq.gz_trim.fastq.gz.Zm1.bam_sortedRG.bam 253_S113_R1_001.fastq.gz_trim.fastq.gz.Zm1.bam_sortedRG.bam 254_S114_R1_001.fastq.gz_trim.fastq.gz.Zm1.bam_sortedRG.bam 255_S115_R1_001.fastq.gz_trim.fastq.gz.Zm1.bam_sortedRG.bam 256_S116_R1_001.fastq.gz_trim.fastq.gz.Zm1.bam_sortedRG.bam 257_S117_R1_001.fastq.gz_trim.fastq.gz.Zm1.bam_sortedRG.bam 258_S118_R1_001.fastq.gz_trim.fastq.gz.Zm1.bam_sortedRG.bam 259_S119_R1_001.fastq.gz_trim.fastq.gz.Zm1.bam_sortedRG.bam 261_S120_R1_001.fastq.gz_trim.fastq.gz.Zm1.bam_sortedRG.bam 263_S121_R1_001.fastq.gz_trim.fastq.gz.Zm1.bam_sortedRG.bam 264_S122_R1_001.fastq.gz_trim.fastq.gz.Zm1.bam_sortedRG.bam 267_S123_R1_001.fastq.gz_trim.fastq.gz.Zm1.bam_sortedRG.bam 268_S124_R1_001.fastq.gz_trim.fastq.gz.Zm1.bam_sortedRG.bam 269_S125_R1_001.fastq.gz_trim.fastq.gz.Zm1.bam_sortedRG.bam 270_S126_R1_001.fastq.gz_trim.fastq.gz.Zm1.bam_sortedRG.bam 271_S127_R1_001.fastq.gz_trim.fastq.gz.Zm1.bam_sortedRG.bam 272_S128_R1_001.fastq.gz_trim.fastq.gz.Zm1.bam_sortedRG.bam 275_S129_R1_001.fastq.gz_trim.fastq.gz.Zm1.bam_sortedRG.bam 276_S130_R1_001.fastq.gz_trim.fastq.gz.Zm1.bam_sortedRG.bam 278_S131_R1_001.fastq.gz_trim.fastq.gz.Zm1.bam_sortedRG.bam 281_S132_R1_001.fastq.gz_trim.fastq.gz.Zm1.bam_sortedRG.bam 282_S133_R1_001.fastq.gz_trim.fastq.gz.Zm1.bam_sortedRG.bam 283_S134_R1_001.fastq.gz_trim.fastq.gz.Zm1.bam_sortedRG.bam 284_S135_R1_001.fastq.gz_trim.fastq.gz.Zm1.bam_sortedRG.bam 285_S136_R1_001.fastq.gz_trim.fastq.gz.Zm1.bam_sortedRG.bam 286_S137_R1_001.fastq.gz_trim.fastq.gz.Zm1.bam_sortedRG.bam 28_S138_R1_001.fastq.gz_trim.fastq.gz.Zm1.bam_sortedRG.bam 295_S139_R1_001.fastq.gz_trim.fastq.gz.Zm1.bam_sortedRG.bam 296_S140_R1_001.fastq.gz_trim.fastq.gz.Zm1.bam_sortedRG.bam 2_S141_R1_001.fastq.gz_trim.fastq.gz.Zm1.bam_sortedRG.bam 30_S142_R1_001.fastq.gz_trim.fastq.gz.Zm1.bam_sortedRG.bam 34_S143_R1_001.fastq.gz_trim.fastq.gz.Zm1.bam_sortedRG.bam 3988_S171_R1_001.fastq.gz_trim.fastq.gz.Zm1.bam_sortedRG.bam 39_S144_R1_001.fastq.gz_trim.fastq.gz.Zm1.bam_sortedRG.bam 44_S145_R1_001.fastq.gz_trim.fastq.gz.Zm1.bam_sortedRG.bam 5_S146_R1_001.fastq.gz_trim.fastq.gz.Zm1.bam_sortedRG.bam 68_S147_R1_001.fastq.gz_trim.fastq.gz.Zm1.bam_sortedRG.bam 6_S148_R1_001.fastq.gz_trim.fastq.gz.Zm1.bam_sortedRG.bam 70_S149_R1_001.fastq.gz_trim.fastq.gz.Zm1.bam_sortedRG.bam 71_S150_R1_001.fastq.gz_trim.fastq.gz.Zm1.bam_sortedRG.bam 73_S151_R1_001.fastq.gz_trim.fastq.gz.Zm1.bam_sortedRG.bam 76_S152_R1_001.fastq.gz_trim.fastq.gz.Zm1.bam_sortedRG.bam 77_S153_R1_001.fastq.gz_trim.fastq.gz.Zm1.bam_sortedRG.bam 84_S154_R1_001.fastq.gz_trim.fastq.gz.Zm1.bam_sortedRG.bam 85_S155_R1_001.fastq.gz_trim.fastq.gz.Zm1.bam_sortedRG.bam 86_S156_R1_001.fastq.gz_trim.fastq.gz.Zm1.bam_sortedRG.bam 87_S157_R1_001.fastq.gz_trim.fastq.gz.Zm1.bam_sortedRG.bam 89_S158_R1_001.fastq.gz_trim.fastq.gz.Zm1.bam_sortedRG.bam 8_S159_R1_001.fastq.gz_trim.fastq.gz.Zm1.bam_sortedRG.bam 91_S160_R1_001.fastq.gz_trim.fastq.gz.Zm1.bam_sortedRG.bam 92_S161_R1_001.fastq.gz_trim.fastq.gz.Zm1.bam_sortedRG.bam 93_S162_R1_001.fastq.gz_trim.fastq.gz.Zm1.bam_sortedRG.bam 94_S163_R1_001.fastq.gz_trim.fastq.gz.Zm1.bam_sortedRG.bam 95_S164_R1_001.fastq.gz_trim.fastq.gz.Zm1.bam_sortedRG.bam 96_S165_R1_001.fastq.gz_trim.fastq.gz.Zm1.bam_sortedRG.bam RSY_S170_R1_001.fastq.gz_trim.fastq.gz.Zm1.bam_sortedRG.bam >& NGSEP_bwaZm1.log

# Filter the varaints -1
java -jar NGSEPcore_4.2.1/NGSEPcore.jar VCFFilter -i NGSEP_bwaZm1.vcf -o NGSEP_bwaZm1-1.vcf -q 40 -minMAF 0.05  -s 
java -jar NGSEPcore_4.2.1/NGSEPcore.jar VCFFilter -i NGSEP_bwaPG.vcf -o NGSEP_bwaPG-1.vcf -q 40 -minMAF 0.05  -s 
java -jar NGSEPcore_4.2.1/NGSEPcore.jar VCFFilter -i NGSEP_bwaXJDY.vcf -o NGSEP_bwaXJDY-1.vcf -q 40 -minMAF 0.05  -s 

# Filter the varaints -2, Rm_entries.txt containns the list of entries with >90% missing SNPs
java -jar NGSEPcore_4.2.1/NGSEPcore.jar VCFFilter -i NGSEP_bwaZm1-1.vcf -o NGSEP_bwaZm1-2.vcf -saf Rm_entries.txt  -fs
java -jar NGSEPcore_4.2.1/NGSEPcore.jar VCFFilter -i NGSEP_bwaPG-1.vcf -o NGSEP_bwaPG-2.vcf -saf Rm_entries.txt  -fs
java -jar NGSEPcore_4.2.1/NGSEPcore.jar VCFFilter -i NGSEP_bwaXJDY-1.vcf -o NGSEP_bwaXJDY-2.vcf -saf Rm_entries.txt  -fs


# Next read depth and missing. values filter - 3(1), for Zm1 dataset
module load bcftools/1.16-gcc-8.2.0-5d4xg4y 
java -jar NGSEPcore_4.2.1/NGSEPcore.jar VCFFilter -i  NGSEP_bwaZm1-2.vcf -o  NGSEPZm1-bwa-RD20m160.vcf -minRD 20 -m 160
java -jar NGSEPcore_4.2.1/NGSEPcore.jar VCFFilter -i  NGSEP_bwaZm1-2.vcf -o  NGSEPZm1-bwa-RD40m140.vcf -minRD 40 -m 140
java -jar NGSEPcore_4.2.1/NGSEPcore.jar VCFFilter -i  NGSEP_bwaZm1-2.vcf -o  NGSEPZm1-bwa-RD60m100.vcf -minRD 60 -m 100

bgzip NGSEPZm1-bwa-RD20m160.vcf
bgzip NGSEPZm1-bwa-RD40m140.vcf
bgzip NGSEPZm1-bwa-RD60m100.vcf
tabix NGSEPZm1-bwa-RD20m160.vcf.gz
tabix NGSEPZm1-bwa-RD40m140.vcf.gz
tabix NGSEPZm1-bwa-RD60m100.vcf.gz

# Combine files to get one  for Zm1
bcftools concat -a -d all -o NGSEPZm1_merged204060.vcf -O v NGSEPZm1-bwa-RD20m160.vcf.gz NGSEPZm1-bwa-RD40m140.vcf.gz NGSEPZm1-bwa-RD20m160.vcf.gz

# Adding AD field
java -cp NGSEPcore/NGSEPcore_4.2.2.jar ngsep.vcf.VCFGenerateADField NGSEPZm1_merged204060.vcf > NGSEPZm1_merged204060AD.vcf

# Next read depth and missing. values filter - 3(2), for PG dataset
module load bcftools/1.16-gcc-8.2.0-5d4xg4y 
java -jar NGSEPcore_4.2.1/NGSEPcore.jar VCFFilter -i  NGSEP_bwaPG-2.vcf -o  NGSEPPG-bwa-RD20m160.vcf -minRD 20 -m 160
java -jar NGSEPcore_4.2.1/NGSEPcore.jar VCFFilter -i  NGSEP_bwaPG-2.vcf -o  NGSEPPG-bwa-RD40m140.vcf -minRD 40 -m 140
java -jar NGSEPcore_4.2.1/NGSEPcore.jar VCFFilter -i  NGSEP_bwaPG-2.vcf -o  NGSEPPG-bwa-RD60m100.vcf -minRD 60 -m 100

bgzip NGSEPPG-bwa-RD20m160.vcf
bgzip NGSEPPG-bwa-RD40m140.vcf
bgzip NGSEPPG-bwa-RD60m100.vcf
tabix NGSEPPG-bwa-RD20m160.vcf.gz
tabix NGSEPPG-bwa-RD40m140.vcf.gz
tabix NGSEPPG-bwa-RD60m100.vcf.gz

# Combine files to get one  for PG
bcftools concat -a -d all -o NGSEPPG_merged204060.vcf -O v NGSEPPG-bwa-RD20m160.vcf.gz NGSEPPG-bwa-RD40m140.vcf.gz NGSEPPG-bwa-RD20m160.vcf.gz

# Adding AD field
java -cp NGSEPcore/NGSEPcore_4.2.2.jar ngsep.vcf.VCFGenerateADField NGSEPPG-maxDP200_merged204060.vcf > NGSEPPG-maxDP200_merged204060AD.vcf

# Next read depth and missing. values filter - 3(3), for XJDY dataset only
module load bcftools/1.16-gcc-8.2.0-5d4xg4y 
java -jar NGSEPcore_4.2.1/NGSEPcore.jar VCFFilter -i  NGSEP_bwaXJDY-2.vcf -o  NGSEPXJDY-bwa-RD20m160.vcf -minRD 20 -m 160
java -jar NGSEPcore_4.2.1/NGSEPcore.jar VCFFilter -i  NGSEP_bwaXJDY-2.vcf -o  NGSEPXJDY-bwa-RD40m140.vcf -minRD 40 -m 140
java -jar NGSEPcore_4.2.1/NGSEPcore.jar VCFFilter -i  NGSEP_bwaXJDY-2.vcf -o  NGSEPXJDY-bwa-RD60m100.vcf -minRD 60 -m 100

bgzip NGSEPXJDY-bwa-RD20m160.vcf
bgzip NGSEPXJDY-bwa-RD40m140.vcf
bgzip NGSEPXJDY-bwa-RD60m100.vcf
tabix NGSEPXJDY-bwa-RD20m160.vcf.gz
tabix NGSEPXJDY-bwa-RD40m140.vcf.gz
tabix NGSEPXJDY-bwa-RD60m100.vcf.gz

# Combine files to get one  for XJDY
bcftools concat -a -d all -o NGSEPXJDY_merged204060.vcf -O v NGSEPXJDY-bwa-RD20m160.vcf.gz NGSEPXJDY-bwa-RD40m140.vcf.gz NGSEPXJDY-bwa-RD20m160.vcf.gz

# Adding AD field
java -cp NGSEPcore/NGSEPcore_4.2.2.jar ngsep.vcf.VCFGenerateADField NGSEPXJDY_merged204060.vcf > NGSEPXJDY_merged204060AD.vcf

# Now we have three VCF files- NGSEPZm1_merged204060AD.vcf, NGSEPPG_merged204060AD.vcf and NGSEPXJDY_merged204060AD.vcf














