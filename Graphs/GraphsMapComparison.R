# Load necessary packages 
library(ggplot2)
library(gridExtra)

# load the file with these columns (Markers, Chromosome, Position_cM, Position_bp, and MapType i.e. reference genome name) 
data <- read.csv("LinkageMapComparison.csv",
                 header=TRUE)
data$Chromosome <- as.factor(data$Chromosome)
data$Map_Type <- as.factor(data$Map_Type)
str(data)
data$Chromosome2 <- paste("Chr", data$Chromosome, sep="") # Add 'Chr' in front of chromosome number
data$Position_bp <- data$Position_bp/1000000 # Convert bp to MB

# Create the plot for genetic map comparison
plot1 <- ggplot(data, aes(x = Position_cM, y = Chromosome2, col = Map_Type, 
                 fill=Map_Type)) + 
  geom_col(width = 0.005, position = position_dodge(width = 0.7)) +
  scale_fill_manual(values = c("#E69F00", "#56B4E9", "#009E73")) +
  geom_point(shape=3, size=4, position = position_dodge(width = 0.7)) +
  facet_wrap(~ Chromosome2, scales = "free_x", ncol=4) + 
  labs(x = "Genetic Distance (cM)", y = NULL, color = "Map Type") + 
  theme(axis.text=element_text(size=15, face="bold"),
        axis.title=element_text(size=15, face="bold"),
        legend.position = "none",
        legend.title = element_blank(),
        strip.text = element_text(size = 15, face="bold"),
        axis.text.x=element_blank(),
        plot.title = element_text(size = 15, face="bold")) +
  scale_color_manual(values = c("#E69F00", "#56B4E9", "#009E73")) +
  scale_x_continuous(breaks=seq(0, 600, 50)) +
  coord_flip() +
  guides(color = "none") +
  ggtitle("a)")


# Physical Map Comparison
plot2 <- ggplot(data, aes(x = Position_bp, y = Chromosome2, col = Map_Type, 
                 fill=Map_Type)) + 
  geom_col(width = 0.003, position = position_dodge(width = 0.7)) +
  scale_fill_manual(values = c("#E69F00", "#56B4E9", "#009E73")) +
  geom_point(shape=3, size=4, position = position_dodge(width = 0.7)) +
  facet_wrap(~ Chromosome2, scales = "free_x", ncol=4) + 
  labs(x = "Physical Distance (MB)", y = NULL, color = "Map Type") + 
  theme(axis.text=element_text(size=15, face="bold"),
        axis.title=element_text(size=15, face="bold"),
        legend.position = "bottom",
        legend.title = element_blank(),
        strip.text = element_text(size = 15, face="bold"),
        axis.text.x=element_blank(),
        plot.title = element_text(size = 15, face="bold")) +
  scale_color_manual(values = c("#E69F00", "#56B4E9", "#009E73")) +
  scale_x_continuous(breaks=seq(0, 100, 20)) +
  coord_flip() +
  guides(color = "none") +
  ggtitle("b)")

grid.arrange(plot1, plot2, nrow=2, ncol=1)

# Add variant points in pangenome from another genomes 
data0 <- read.csv("NGSEPPG_merged204060ADupdog.csv")
data1 <- data0[,1:4]
for (i in 1:8) {
  data1[data1$Ref.Chrom==i,]$Ref.Position <- paste("Chr",i, data1[data1$Ref.Chrom==i,]$Ref.Position, sep="_")
}

data2 <- read.csv("NGSEPZm1_merged204060ADupdog.csv")
for (i in 1:8) {
  data2[data2$Ref.Chrom==i,]$Ref.Position <- paste("Chr",i, data2[data2$Ref.Chrom==i,]$Ref.Position, sep="_")
}

OtherVariants <- data1[-which(data2$Ref.Position %in% data1$Ref.Position),]
OtherVariants <- data0[which(data0$Marker.Name %in% OtherVariants$Marker.Name),]

data3 <- data[data$Map_Type=="Pangenome",]
dataPImp <- row.names(data3[which(data3$Position_bp %in% (OtherVariants$Ref.Position/1000000)),])
data$Group[as.numeric(dataPImp)] <- "Imp"
data$Map_Type <- as.character(data$Map_Type)
data$Map_Type[as.numeric(dataPImp)] <- "Pangenome Markers not in ZhongmuNo.1"
data$Map_Type <- as.factor(data$Map_Type)
which(data$Map_Type == 'Pangenome Markers not in ZhongmuNo.1')

# Map physical position against genetic distance
plot <- ggplot(data, aes(x = Position_bp, y = Position_cM, col=Map_Type)) +
  geom_point(size=1) +
  scale_color_manual(values = c("#E69F00","red","#56B4E9", "#009E73")) +
  facet_wrap(~ Chromosome2, scales = "free_x", ncol=4) +
  labs(x = "Physical Distance (MB)", y = "Genetic Distance (cM)", color = "Map Type") + 
  theme(axis.text=element_text(size=15, face="bold"),
        axis.title=element_text(size=15, face="bold"),
        legend.position = "top",
        legend.title = element_blank(),
        legend.text=element_text(size=15, face="bold"),
        strip.text = element_text(size = 15, face="bold")) +
        scale_y_continuous(breaks=seq(0, 600, 50)) +
        scale_x_continuous(breaks=seq(0, 115, 15))
plot +  guides(color = guide_legend(override.aes = list(size = 5)))
#  ggtitle("Physical vs Genetic Map Comparison")
