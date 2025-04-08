# Development of pangenome using 5 different genomes
# 1. Consensus genome from a autotetraploid genotype of Zhongmuno.1 cultivar (Shen et. al. 2020)
# 2. All four haplotypes of tetraploid resolved genome XinJiang DaYe cultivar (Chen et.al. 2020)
# 3. All four haplotypes of tetraploid resolved genome ZhongmuNo.4 (Long et al. 2022)
# 4. Cultivated Alfalfa at diploid level (CADL), scaffolds level genome (Legume Information System)
# 5. RegenSY27x, scaffolds level genome (Legume Information System)

# Preparing chromosome level base reference genome for minigraph, removing contigs from ZhongmuNo.1 genome
zgrep ">" ZhongmuNo.1_genome.fasta.gz | grep "^>Chr" > header
# seqkit v0.4.0
seqkit grep -f header ZhongmuNo.1_genome.fasta  -o zm-1.fasta

# Prepare files for minigraph by separating all four haplotypes in allle-aware haplotype resolved assemblies
zgrep ">" final_genome.fasta.gz | grep "^>chr" | grep ".1$" | sed  "s/>//g" > headerhap1
zgrep ">" final_genome.fasta.gz | grep "^>chr" | grep ".2$" | sed  "s/>//g" > headerhap2
zgrep ">" final_genome.fasta.gz | grep "^>chr" | grep ".3$" | sed  "s/>//g" > headerhap3
zgrep ">" final_genome.fasta.gz | grep "^>chr" | grep ".4$" | sed  "s/>//g" > headerhap4

# seqkit v0.4.0, separate fasta
seqkit grep -f headerhap1 final_genome.fasta  -o final_hap1.fasta
seqkit grep -f headerhap2 final_genome.fasta  -o final_hap2.fasta
seqkit grep -f headerhap3 final_genome.fasta  -o final_hap3.fasta
seqkit grep -f headerhap4 final_genome.fasta  -o final_hap4.fasta

# For ZhongmuNo.4 genome, extract headers
zgrep ">" zm-4.genome.fasta.gz | grep "^>chr" | grep "_1$" | sed  "s/>//g" > headerhap1
zgrep ">" zm-4.genome.fasta.gz | grep "^>chr" | grep "_2$" | sed  "s/>//g" > headerhap2
zgrep ">" zm-4.genome.fasta.gz | grep "^>chr" | grep "_3$" | sed  "s/>//g" > headerhap3
zgrep ">" zm-4.genome.fasta.gz | grep "^>chr" | grep "_4$" | sed  "s/>//g" > headerhap4

# seqkit v0.4.0, separate fasta
seqkit grep -f headerhap1 zm-4.genome.fasta  -o zm-4_hap1.fasta
seqkit grep -f headerhap2 zm-4.genome.fasta  -o zm-4_hap2.fasta
seqkit grep -f headerhap3 zm-4.genome.fasta  -o zm-4_hap3.fasta
seqkit grep -f headerhap4 zm-4.genome.fasta  -o zm-4_hap4.fasta

# CADL and RegenSY were used as such

# minigraph v0.20-r559
./minigraph/minigraph -cxggs -t16 zm-1.fasta zm-4_hap1.fasta zm-4_hap2.fasta zm-4_hap3.fasta zm-4_hap4.fasta final_hap1.fasta final_hap2.fasta final_hap3.fasta final_hap4_vg.fasta CADL.fna RSY.fna > Zm1Zm4XinJiangDaYeCADLRSY.gfa 

# Convert it to stable fasta format usinng gfatools v0.5-r253-dirty
gfatools gfa2fa -s Zm1Zm4XinJiangDaYeCADLRSY.gfa > Zm1Zm4XinJiangDaYeCADLRSY.fasta