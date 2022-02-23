library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(ggpubr)
library(ggvenn)
library(scales)
library(ggExtra)
library(reshape2)
library(cluster)
library(parallel)
library(pbapply)
library(purrr)
library(pheatmap)
library(dendextend)
####################
# Custom Functions #
####################

VCFParsing <- function(vcf){
	# Preparing for the final part
	FinalName <- gsub(".*\\/|\\.vcf", "",vcf)
	# Getting the file Ready
	vcfFile <- read.delim(file = vcf,comment.char = "#", header = F,
			      col.names = c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "unknown")) %>%
			filter(FILTER == "PASS") %>% as_tibble()
		#filter(REF != "C" & ALT != "T") %>% filter(REF != "G" & ALT != "A") %>% as_tibble() # Filtered out the G -> A and C -> T Transitions

	if(nrow(vcfFile) == 0){
		cat("No SNPs found in",FinalName, "\n")
		return(NA)
	}
	
	tmp <- unlist(strsplit(vcfFile$FORMAT[1], ":"))
	tmp2 <- unlist(strsplit(gsub("=.*?(?=;)|=.*$","",vcfFile$INFO[1], perl = T), ";"))
	vcfFile <- vcfFile %>% select(-FORMAT) %>% separate(unknown, into = tmp, sep = ":") %>%
		mutate(INFO = gsub("(?<=;).*?=|^.*?=", "", INFO, perl = T)) %>% separate(INFO, into= tmp2, sep = ";")
	vcfFile[,c(6,8:11,14)]<- sapply(vcfFile[,c(6,8:11,14)], as.numeric)
	return(vcfFile)
}

#########################
# Preparing the colours #
#########################
colour <- c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6',
'#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3',
'#808000', '#ffd8b1', '#000075', '#808080', '#f0f0f0', '#000000')

ann_colors = list(Pathovar = c("Ancient" = "#644117", "Nonpathogenic" = colour[22], "AIEC" = "#22662A",# "EAEC" = "#39a275", <- Not directly identified in Pan-Genome
			       "EHEC" = colour[11], "EPEC" = colour[7], "ETEC" = colour[16], "ExPEC" = "#223fa9",
			       "Pathogenic" = colour[17], "STEC" = colour[9], "Shigella" = colour[20],
			       "Unknown" = "#5c5c5c"),
		  Phylogroup = c("Ancient" = "#644117", "A" = colour[4],"B1" = colour[2], "B2" = colour[1], "C" = colour[11], "Clade I" = "#a9a9a9",
			       "D" = colour[3], "E" = colour[8],"F" = colour[18], "G" = colour[5],"Shigella" = colour[20], "Unknown" = colour[22]),
		  Ancient = c("Yes" = colour[22], "No" = colour[21]))

alphaList = list(Pathovar = c("Ancient" = 1, "Nonpathogenic" = 1, "AIEC" = 1,# "EAEC" = 1, <- Not directly identified in Pan-Genome 
			       "EHEC" = 1, "EPEC" = 1, "ETEC" = 1, "ExPEC" = 1,
			       "Pathogenic" = 1, "STEC" = 1, "Shigella" = 1,
			       "Unknown" = 0.25))
# Getting the Phylogroup Data involved
phylogroupFrame <- read.delim("AllFinal_phylogroups.txt", header = F, stringsAsFactors = F)[,c(1,5)]
colnames(phylogroupFrame) <- c("Genome", "Phylogroup")
phylogroupFrame$Genome <- gsub(".fna|.fasta|.fasta.ref","",phylogroupFrame$Genome) # Removing the last bit
#which(phylogroupFrame$Genome == "991919.15")
tmp <- which(grepl("562.12824|562.29124|562.29125|562.7956", phylogroupFrame$Genome))
phylogroupFrame$Phylogroup[tmp] <- "Shigella"
phylogroupFrame$Phylogroup[c(68,234,427)] <- c("C","C","E")
phylogroupFrame[456,] <- list("O104","B1")
phylogroupFrame[457,] <- list("562.5162","Clade I")
phylogroupFrame <- phylogroupFrame %>% mutate(Phylogroup = replace(Phylogroup, grepl("Shigella", Genome), "Shigella"),
			   Phylogroup = replace(Phylogroup, grepl("Ancient", Genome), "Ancient"),
				Phylogroup = replace(Phylogroup, is.na(Phylogroup), "Unknown"),
				Phylogroup = replace(Phylogroup, grepl("cladeI", Phylogroup), "Clade I"))
phylogroupFrame$Genome[452] <- "AncientEcoli"

# Let's load in the pathovar information that I have available
patricTable <- read.csv("PatricTableEdited.csv", header = T) %>% as_tibble() %>% select(Genome.ID, Pathovar)
colnames(patricTable)[1]  <- "Genome"
virFrame <- phylogroupFrame %>% left_join(patricTable) %>% mutate(Pathovar = replace(Pathovar, Pathovar == "", "Unknown"),
						      Pathovar = replace(Pathovar, Phylogroup == "Shigella", "Shigella"),
						      Pathovar = replace(Pathovar, Phylogroup == "Ancient", "Ancient"),
						      Pathovar = replace(Pathovar, grepl("UPEC|NMEC|BSI|IPEC", Pathovar), "ExPEC"))
virFrame %>% pull(Pathovar) %>% unique()
#########################
# Reading in the Depths #
#########################
depthDf <- lapply("WholePanGeneome.tab.gz",function(f){
		tmp <- as_tibble(read.delim(f, header = F, col.names = c("Gene", "Pos", "Coverage")))
		tmp$Genome <- gsub(".*/", "", gsub(".tab.gz","",f)) 	
		tmp <- tmp %>% group_by(Genome, Gene) %>%
			summarize(MeanCoverage = mean(Coverage), sdCoverage = sd(Coverage), PercentCoverage = sum(Coverage > 0)/length(Coverage), .groups = "drop") %>%
			mutate(CV = sdCoverage/MeanCoverage)
		return(tmp)
	}) %>% bind_rows() # Can be used for multiple files in this format

depthDfNotAvg <- lapply("WholePanGeneome.tab.gz", function(f){
		tmp <- as_tibble(read.delim(f, header = F, col.names = c("Gene", "Pos", "Coverage")))
		tmp$Genome <- gsub(".*/", "", gsub(".tab.gz","",f)) 	
		return(tmp)
	}) %>% bind_rows()
#depthDf <- depthDf %>% filter(!(Genome %in% c("IAI1", "ESC_NB8751AA_AS", "ATCC11231", "ESC_LB2165AA_AS")))
#depthDf$Genome <- unname(unlist(sapply(depthDf$Genome, function(x){ifelse(grepl("^ERR",x), etecList[x],x)})))
depthDf <- depthDf %>% mutate(Genome = replace(Genome, Genome == "PanGenomeMappingFeb", "KaeroEcoli"))
#unique(depthDf$Genome)
#depthDf %>% filter(!(Genome %in% c("PanGenomeDepths","PanGenomeDepthsDup")))

############# Core Gene Presence #####
roaryOutput <- as_tibble(read.delim("gene_presence_absence.Rtab"))
colnames(roaryOutput)[2:ncol(roaryOutput)] <- gsub("X|\\.fna|\\.scaffold|\\.genome|\\.result|_genomic", "", colnames(roaryOutput)[2:ncol(roaryOutput)])
#roaryOutput <- roaryOutput %>% select(-c(IAI1, ATCC11231,ESC_LB2165AA_AS))
# Finding the Core Genome
removedgenomes <- which(colnames(roaryOutput) %in% c("562.7692","562.7736","562.7622","562.7382") | grepl("Shigella", colnames(roaryOutput)))
ecoliOnly <- roaryOutput[,-removedgenomes]
coreGenes <- rowSums(ecoliOnly[,-1]) >= floor(0.95 * ncol(ecoliOnly))
#coreGenes <- rowSums(ecoliOnly[,-1]) >= floor(0.99 * ncol(ecoliOnly))
coreGenes <- roaryOutput$Gene[coreGenes]

ancientData <- depthDf %>% filter(Genome == "KaeroEcoli", MeanCoverage >= 10, CV <=1) %>%
	pivot_wider(names_from = c(Genome), values_from=MeanCoverage) %>% select(-c(sdCoverage, PercentCoverage, CV)) #%>%
	#mutate(Genome = "10x - CV Filter")

#####################################################################
# Gene Presence table
AllPA <- roaryOutput %>% left_join(ancientData) %>% mutate(KaeroEcoli = ifelse(is.na(KaeroEcoli)|KaeroEcoli == 0, 0,1))

# Now to make a boxplot to compare our gene to everyone else
tmp <- AllPA %>% select(-Gene) %>% summarize_all(sum) %>% t() %>% as.data.frame()
tmp$Genome <- rownames(tmp)
geneCounts <- tmp %>% as_tibble() #%>% filter(!(Genome %in% c("562.7692","562.7736","562.7622","562.7382")))
colnames(geneCounts)[1] <- "Genes"
geneCounts %>% mutate(Outlier = ifelse(Genome %in% c("562.7692","562.7736","562.7622","562.7382"), "Outlier", "")) %>% filter(Genome != "KaeroEcoli") %>% 
       	ggplot(aes(y = Genes)) + geom_boxplot() + theme_bw() +
	geom_point(data = geneCounts %>% filter(Genome == "KaeroEcoli"), aes(colour = Genome, x = 0)) +
	geom_point(data = geneCounts %>% filter(Genome %in% c("562.7692","562.7736","562.7622","562.7382")), aes(x = 0), shape = 3) +
	scale_colour_manual(values = colour) + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank(), legend.position = "bottom")
ggsave("~/GeneCountsEcoli.pdf", width = 6, height = 9)

perCovMean <- depthDf %>%
       	filter(Genome == "KaeroEcoli", MeanCoverage > 0) %>% ggplot(aes(x= PercentCoverage,y = MeanCoverage)) + geom_point() + theme_bw() + scale_x_continuous(breaks = pretty_breaks(n = 10)) +
	scale_y_continuous(breaks = pretty_breaks(n = 10)) + ggtitle("Percent Coverage vs Mean Gene Coverage")
perCovMean <- ggMarginal(perCovMean, type = "densigram")
CVMean <- depthDf %>% filter(Genome == "KaeroEcoli", MeanCoverage > 0) %>% ggplot(aes(x= CV,y = MeanCoverage)) + geom_point() + theme_bw()+ scale_x_continuous(breaks = pretty_breaks(n = 10)) +
       	scale_y_continuous(breaks = pretty_breaks(n = 10))+ ggtitle("Coefficient of Variation vs Mean Gene Coverage")
CVMean <- ggMarginal(CVMean, type = "densigram")

pdf("~/CoverageMetrics.pdf", width = 9, height = 9)
perCovMean
plot.new()
CVMean
dev.off()
#####################################################################
# Let's put the gene coverage data (ie. histograms) here
ancientGenesOnly <- depthDf %>% filter(Genome == "KaeroEcoli", CV <= 1)
ancientGenesOnly$Status <- ifelse(ancientGenesOnly$Gene %in% coreGenes, "Core", "Accessory")

# Quick summary Statistics
meanSD <- ancientGenesOnly %>% filter(MeanCoverage >= 10)%>% summarize(Mean = mean(MeanCoverage), SD = sd(MeanCoverage), confInt = qnorm(0.975)*SD/sqrt(length(MeanCoverage)), high = Mean + confInt, low = Mean - confInt)
ancientGenesOnly %>% filter(MeanCoverage >= 10)%>% group_by(Status)%>% summarize(Mean = mean(MeanCoverage), SD = sd(MeanCoverage), confInt = qnorm(0.975)*SD/sqrt(length(MeanCoverage)), high = Mean + confInt, low = Mean - confInt) %>% as.data.frame()

p1 <- ancientGenesOnly %>% ggplot(aes(x = MeanCoverage)) + geom_histogram(fill = "#2e294e", colour = "black", bins = 60) +
	scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
	ylab("Genes") + theme_bw() + xlab("") + theme(axis.text.x = element_blank()) + scale_y_continuous(breaks = scales::pretty_breaks(n = 10), limits = c(0, 550))

breaks <- pretty(range(ancientGenesOnly$MeanCoverage), n = nclass.FD(ancientGenesOnly$MeanCoverage), min.n = 1)
bwidth <- breaks[2]-breaks[1]

ancientGenesOnly <- ancientGenesOnly %>% mutate(Status = factor(Status, levels = c("Core", "Accessory")))
p2 <- ancientGenesOnly  %>% ggplot(aes(x = MeanCoverage, fill = Status)) + geom_histogram(position = "identity", alpha = 0.75, colour = "black", binwidth = 1) +
	geom_rect(data = meanSD, inherit.aes = F, aes(ymin = -Inf, ymax = Inf, xmax = Inf, xmin = Mean + (2*SD)),colour = "black", fill = "black", alpha = 0.25, lty = 2) +
	scale_fill_manual(values = c(Core = "#f8333c", Accessory = "#007dba")) + theme_bw() +
	geom_vline(xintercept = 10, lty = 2) +
	xlab("Mean Read Depth") + ylab("Genes") + theme(legend.position = "bottom") + 
	scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) + scale_y_continuous(breaks = scales::pretty_breaks(n = 10), limits = c(0, 500))
p2
ggsave("~/Documents/University/EcoliPaperV2/Figures/GeneDepthHistogramWhole.pdf", width = 6, height = 4)

# A simple T test to compare the core and accessory genomes
coreMeans <- ancientGenesOnly %>% filter(MeanCoverage >= 10, Status == "Core") %>% pull(MeanCoverage)
accessMeans <- ancientGenesOnly %>% filter(MeanCoverage >= 10, Status != "Core") %>% pull(MeanCoverage)
t.test(coreMeans, accessMeans)

# Getting the estimated genome length and PCov
foundGenes <- ancientGenesOnly %>% filter(MeanCoverage >= 10) %>% pull(Gene)
depthDfNotAvg %>% filter(Gene %in% foundGenes) %>% summarize(sum(Coverage > 0)/length(Coverage)) %>% pull()

# Actually calculating the Gene Copy numbers
ancientGenesOnly <- ancientGenesOnly %>% filter(MeanCoverage >= 10) %>% mutate(CopyNumber = MeanCoverage/mean(MeanCoverage))

###########################
### The core Genome P/A ###
###########################
removedgenomes <- which(colnames(roaryOutput) %in% c("562.7692","562.7736","562.7622","562.7382"))

# Going to do this in a more sane way
coreDf <- AllPA %>% filter(Gene %in% coreGenes)
coreDf <- coreDf[,-removedgenomes] # Removing those outliers

# Preparing the data
coreData <- as.data.frame(coreDf)
rownames(coreData) <- coreData[,1]
coreData <- coreData[,-1]
genomes <- colnames(coreData)
coreData <- apply(coreData, MARGIN = 1,FUN = as.numeric)
rownames(coreData) <- genomes

coreDataSub <- coreData
maxCount <- coreDataSub %>% colSums() %>% max()
ind <- coreDataSub %>% colSums() == maxCount
coreDataSub <- coreDataSub[,!ind]

# Getting the heatmap sorted
rownames(coreDataSub)[which(rownames(coreDataSub) == "KaeroEcoli")] <- "AncientEcoli"
#rownames(coreDataSub)[which(rownames(coreDataSub) %in% c("KaeroEcoli","Lib4_3_bin","Lib4_7_bin","Lib4_8_bin"))] <- c("AncientEcoli", "Zape2.1", "Zape2.2", "Zape3")
rownamesHeatmap <- ifelse(rownames(coreDataSub) == "AncientEcoli", "AncientEcoli","")
virFrame <- as.data.frame(virFrame)

rownames(virFrame) <- virFrame$Genome
virFrame <- virFrame[rownames(coreDataSub),] %>% filter(!is.na(Genome))

rownamesHeatmap <- ifelse(rownames(coreDataSub) == "AncientEcoli", "AncientEcoli","")
#rownames(olivierTableHeatmap) <- phylogroup$Genome
#rownames(olivierTableHeatmap) <- olivierTable$Genome

# Dendrogram and order
specDist <- dist(coreDataSub, method = "binary")
clusteredSpec <- hclust(specDist, method = "ward.D2") # 
rownameOrder <- clusteredSpec$order

treeColours <- as.list(ann_colors[[1]])[virFrame[rownameOrder,]$Phylogroup] %>% unlist()
treeColours <- gsub("#ffffff", "#000000", treeColours)
dend <- as.dendrogram(clusteredSpec)
labels_colors(dend) <- treeColours

pdf("~/CVFigures/WholeCoreTree.pdf", width = 24, height = 9)
#pdf("~/Documents/University/EcoliPaperV2/AdditionalFiles/PartsofFigures/SpeciesCoreTree.pdf", width = 6, height = 8)
plot(dend, horiz = F)
legend("topright", legend = names(ann_colors[[1]]), fill = unlist(ann_colors[[1]]), title = "Phylogroup")
dev.off()

heatmapAnnotation <- virFrame[,-1]
#heatmapAnnotation <- data.frame("Phylogroup" = virFrame[,-1])
rownames(heatmapAnnotation) <- virFrame$Genome

# Making the PCoA
geneDist <- dist(coreDataSub, method = "binary")

fit <- cmdscale(geneDist, eig = T, k = 4, add = T)
coord <- fit$points %>% as_tibble()
coord$Genome <- rownames(fit$points)
contrib <- fit$eig/sum(fit$eig) * 100
coord <- coord %>% mutate(AncientLabel = ifelse(grepl("Ancient", Genome), "Ancient Ecoli", ""))
#coord <- coord %>% mutate(EAEC = ifelse(Genome %in% potentialEAEC, "Potentially EAEC", "Not EAEC"))

#olivierTableHeatmap$Genome <- rownames(olivierTableHeatmap)

corePhylo <- coord %>% left_join(virFrame) %>%
       	ggplot(aes(x = V1, y = V2, label = AncientLabel, color = Phylogroup)) +
	geom_hline(yintercept=0, lty = 2, colour = "grey90") + 
	geom_vline(xintercept=0, lty = 2, colour = "grey90") +
	geom_point() +
	geom_text_repel(show.legend = F) +
	xlab(bquote("PCoA 1 ("~.(round(contrib[1],2))~"%)")) +
	ylab(bquote("PCoA 2 ("~.(round(contrib[2],2))~"%)")) +
	theme_classic() + scale_colour_manual(values = ann_colors$Phylogroup) +
	guides(color = guide_legend(nrow = 2), shape = guide_legend(nrow = 2)) +
	theme(legend.position = "bottom")
ggsave(corePhylo, file = "~/Documents/University/EcoliPaperV2/Figures/CoreMDS.pdf", width = 6, height = 4)

# Core data counts for phylogeny
counts <- coreData %>% rowSums()
counts <- tibble(Genome = names(counts), Count = counts)
write.table(counts, sep = "\t", "../CoreGeneCounts.tab")

#################################
### Accessory Genome Analysis ###
#################################
accessDf <- AllPA %>% filter(!(Gene %in% coreGenes))
accessDf <- accessDf[,-removedgenomes] # Removing those outliers

# Preparing the data
accessData <- as.data.frame(accessDf)
rownames(accessData) <- accessData[,1]
accessData <- accessData[,-1]
genomes <- colnames(accessData)
accessData <- apply(accessData, MARGIN = 1,FUN = as.numeric)
rownames(accessData) <- genomes

# We'll be removing the unique genes now
accessDataSub <- accessData
maxCount <- accessDataSub %>% colSums() %>% max()
ind <- accessDataSub %>% colSums() < 5
accessDataSub <- accessDataSub[,!ind]

# Getting the Metadata Setup
rownames(accessDataSub)[which(rownames(accessDataSub) == "KaeroEcoli")] <- "AncientEcoli"
rownamesHeatmap <- ifelse(rownames(accessDataSub) == "AncientEcoli", "AncientEcoli","")
#phylogroupFrame <- as.data.frame(phylogroupFrame)

#rownames(phylogroupFrame) <- phylogroupFrame$Genome
virFrame <- virFrame[rownames(accessDataSub),] %>% filter(!is.na(Genome))

# Dendrogram and order
specDist <- dist(accessDataSub, method = "binary")
clusteredSpec <- hclust(specDist, method = "ward.D2") # 
rownameOrder <- clusteredSpec$order

treeColours <- as.list(ann_colors[[1]])[phylogroupFrame[rownameOrder,]$Phylogroup] %>% unlist()
treeColours <- gsub("#ffffff", "#000000", treeColours)
dend <- as.dendrogram(clusteredSpec)
labels_colors(dend) <- treeColours

#pdf("~/CVFigures/AccessTree.pdf", width = 24, height = 18)
##pdf("~/Documents/University/EcoliPaperV2/AdditionalFiles/PartsofFigures/SpeciesAccessTree.pdf", width = 6, height = 8)
#plot(dend, horiz = F, main = "Accessory Genome Clustering")
#legend("topright", legend = names(ann_colors[[1]]), fill = unlist(ann_colors[[1]]), title = "Phylogroup")
#dev.off()

# Now to cluster the trees
#geneDist <- dist(t(accessDataSub))
######################################
# If I want to use a MDS plot
geneDist <- dist(accessDataSub, method = "binary")

fit <- cmdscale(geneDist, eig = T, k = 4, add = T)
coord <- fit$points %>% as_tibble()
coord$Genome <- rownames(fit$points)
contrib <- fit$eig/sum(fit$eig) * 100
coord <- coord %>% mutate(AncientLabel = ifelse(grepl("Ancient", Genome), "Ancient Ecoli", ""))
#coord <- coord %>% mutate(EAEC = ifelse(Genome %in% potentialEAEC, "Potentially EAEC", "Not EAEC"))

#olivierTableHeatmap$Genome <- rownames(olivierTableHeatmap)

accessVir <- coord %>% left_join(virFrame) %>%
       	ggplot(aes(x = V1, y = V2, label = AncientLabel, color = Pathovar, alpha = Pathovar)) +
	geom_hline(yintercept=0, lty = 2, colour = "grey90") + 
	geom_vline(xintercept=0, lty = 2, colour = "grey90") +
	geom_point() +
	geom_text_repel(show.legend = F) +
	xlab(bquote("PCoA 1 ("~.(round(contrib[1],2))~"%)")) +
	ylab(bquote("PCoA 2 ("~.(round(contrib[2],2))~"%)")) +
	theme_classic() + scale_colour_manual(values = ann_colors$Pathovar) +
	scale_alpha_manual(values = alphaList$Pathovar) +
	guides(color = guide_legend(nrow = 2), shape = guide_legend(nrow = 2))

accessPhylo <- coord %>% left_join(virFrame) %>%
       	ggplot(aes(x = V1, y = V2, label = AncientLabel, color = Phylogroup)) +
	geom_hline(yintercept=0, lty = 2, colour = "grey90") + 
	geom_vline(xintercept=0, lty = 2, colour = "grey90") +
	geom_point() +
	geom_text_repel(show.legend = F) +
	xlab(bquote("PCoA 1 ("~.(round(contrib[1],2))~"%)")) +
	ylab(bquote("PCoA 2 ("~.(round(contrib[2],2))~"%)")) +
	theme_classic() + scale_colour_manual(values = ann_colors$Phylogroup) +
	guides(color = guide_legend(nrow = 2), shape = guide_legend(nrow = 2))

heatmapAnnotation <- virFrame[,-1]
#heatmapAnnotation <- data.frame("Phylogroup" = virFrame[,-1])
rownames(heatmapAnnotation) <- virFrame$Genome

# Making the Heatmap
pheatmap(accessDataSub[rownameOrder,], annotation_row = heatmapAnnotation,
	 color = c("#007dba", "#f8333c"), clustering_method = "ward.D2", legend = F, show_colnames = F,
	 #gaps_row = c(which(rownames(accessDataSub)[rownameOrder] == "AncientEcoli"),which(rownames(accessDataSub)[rownameOrder] == "AncientEcoli") - 1),cluster_rows = F,
	 labels_row = rownamesHeatmap[rownameOrder], fontsize_row = 5, 
	 width = 12, height = 8, border_color = NA, annotation_colors = ann_colors,
	 filename = "~/Documents/University/EcoliPaperV2/AdditionalFiles/PartsofFigures/WholeAccessoryPA.png")#, main = "Core Presence/Absence")

#######################################
### Let's do the Heterozygosity Now ###
#######################################
vcf <- VCFParsing("WholePangenomeHet.vcf")
geneLengths <- read.delim("panGenomeGeneLengths.tab", header =F, col.names = c("Gene", "Length")) %>% as_tibble()
# Getting the total gene set
FoundGenes <- AllPA %>% filter(KaeroEcoli == 1) %>% pull(Gene) 
geneLengths <- geneLengths %>% filter(Gene %in% FoundGenes)

vcfHet <- vcf %>% mutate(Hetero = ifelse(grepl("0/1|1/0", GT), T, F))
vcfPlot <- vcfHet %>% filter(!is.na(FILTER),QUAL >= 100) %>%
       	group_by(CHROM) %>% summarize(Hetero = sum(Hetero)) %>% right_join(geneLengths, by = c("CHROM" = "Gene")) %>%
	mutate(Hetero = replace(Hetero, is.na(Hetero),0)) %>% group_by(CHROM) %>% summarize(HeteroFrac = Hetero/Length) %>%
	mutate(HeteroFrac = replace(HeteroFrac, is.infinite(HeteroFrac), 0), Status = ifelse(CHROM %in% coreGenes, "Core", "Accessory")) %>%
	mutate(Status = factor(Status, levels = c("Core", "Accessory")))

vcfPlot %>% group_by(Status) %>%
       	summarize(Mean = mean(HeteroFrac), SD = sd(HeteroFrac), confInt = qnorm(0.975)*SD/sqrt(length(HeteroFrac)), 
			 high = Mean + confInt, low = Mean - confInt) %>% as.data.frame()
coreHet <-vcfPlot %>% filter(Status == "Core") %>% pull(HeteroFrac)
accessHet <-vcfPlot %>% filter(Status == "Accessory") %>% pull(HeteroFrac)
t.test(coreHet, accessHet)

coreHetNo <-vcfPlot %>% filter(HeteroFrac == 0, Status == "Core") %>% pull(CHROM) %>% length()
accessHetNo <-vcfPlot %>% filter(HeteroFrac == 0, Status == "Accessory") %>% pull(CHROM) %>% length()

hetHist <- vcfPlot %>% ggplot(aes(x = HeteroFrac, fill = Status)) +
       	geom_histogram(position = "identity", alpha = 0.75, colour = "black") +
       #	geom_vline(xintercept = vcfPlot %>% filter(HeteroFrac < 0) %>% summarize(mean(HeteroFrac)) %>% pull(), colour = "red", lty = 2) +
       	theme_bw() +
	scale_x_log10() + annotation_logticks(sides = "b") + scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) +
	ylab("Genes") + xlab("P(Heterozygous)") +
	scale_fill_manual(values = c(Core = "#f8333c", Accessory = "#007dba")) +
	theme(legend.position = "bottom") +
	annotate(geom = "text", x = 0.01, y = 60,color = "#f8333c", label = paste0("Core Genes without Heterozygous Variants: ",coreHetNo)) +
	annotate(geom = "text", x = 0.01, y = 55,color = "#007dba", label = paste0("Accessory Genes without Heterozygous Variants: ",accessHetNo)) 

# Adding the heterozygosity results to the ancientGenesOnly Data
ancientGenesOnly <- ancientGenesOnly %>% full_join(vcfPlot %>% select(-Status), by = c("Gene" = "CHROM")) 

model <- lm(data = ancientGenesOnly %>% filter(HeteroFrac > 0), log10(HeteroFrac) ~ CopyNumber*Status)
r2 <- round(summary(model)$adj.r.squared,2)
tmp <-summary(model)$fstatistic
pval <- round(pf(tmp[1],tmp[2],tmp[3], lower.tail = F),3)

copyHet <- ancientGenesOnly %>% ggplot(aes(x = CopyNumber, y = HeteroFrac, colour = Status)) +
	geom_point() +
	theme_bw() +
	geom_smooth(method = "lm") +
	scale_y_log10() + annotation_logticks(sides = "l") + scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
	xlab("Copy Number") + ylab("P(Heterozygous)") +
	scale_colour_manual(values = c(Core = "#f8333c", Accessory = "#007dba")) +
	annotate(geom = "text", x = 2, y = 3e-3, label = bquote(R[adj]^2 == .(r2))) +
	annotate(geom = "text", x = 2, y = 2e-3, label = bquote(P %~~% .(pval))) +
	theme(legend.position = "bottom") 

ggarrange(hetHist, copyHet, ncol = 1, common.legend = T, legend = "bottom", align = "hv", labels = "AUTO")

ggsave("~/Documents/University/EcoliPaperV2/AdditionalFiles/PartsofFigures/HeterozygosityHist.pdf", width = 8, height = 12)
################################################
### Now to perform the presence absence Test ###
################################################
# Getting the Gene Translations
translatedNR <- as_tibble(read.delim("BlastXAnalysis/ClassifiedGenes.tab", header = T))
colnames(translatedNR)[c(1,2,3)] <- c("Gene", "ID", "Name")
depthDfTrans <- depthDf %>% left_join(translatedNR[, c("Gene", "ID", "Name")])
depthDfTrans$Name <- sapply(1:nrow(depthDfTrans), function(x){ifelse(is.na(depthDfTrans$Name[x]),depthDfTrans$Gene[x] ,depthDfTrans$Name[x])})

# Quickly pulling out the names of those genes with high copy numbers
ancientGenesOnly <- ancientGenesOnly %>% left_join(depthDfTrans %>% select(Gene, ID, Name), by = "Gene")
highCopy <- ancientGenesOnly %>% filter(MeanCoverage > (meanSD$Mean + 2 * meanSD$SD)) %>% arrange(-CopyNumber)
highCopy %>% count(Status)
highCopy %>% summarize(Mean = mean(HeteroFrac), SD = sd(HeteroFrac)) %>% as.data.frame()
highCopy %>% write.table(file = "~/Documents/University/EcoliPaperV2/AdditionalFiles/CopyNumberTable.tab", sep = "\t", quote = F, row.names = F)

# Now to read the pan-genome blast results against K12
k12BlastHits <- read.delim("K12Blast.tab", header =F)[,1:2] %>% as_tibble() # Only care about if we got a hit, not the quality (taken care of before hand)
colnames(k12BlastHits) <- c("Gene", "K12")
k12Genes <- ancientGenesOnly %>% filter(Gene %in% k12BlastHits$Gene) # 3844 Genes....
ancientGenesOnly <- ancientGenesOnly %>% mutate(K12 = ifelse(Gene %in% k12BlastHits$Gene, T, F))
write.table(ancientGenesOnly,"~/Documents/University/EcoliPaperV2/AdditionalFiles/FoundGenes.tab", sep = "\t", row.names = F)

# Identifying the missing core genes
ancientGenes <- depthDf %>% filter(Genome == "KaeroEcoli", MeanCoverage >= 10, CV <= 1)
corePA <- coreGenes %in% ancientGenes$Gene
missingGenes <- coreGenes[!corePA]

status <- ifelse(corePA, "Core Found", "Core Missing")
status <- tibble(status) %>% count(status)
colnames(status)[1] <- "Status"
status
status %>% ggplot(aes(x = "", y = `n`, fill = Status)) + geom_col() + theme_classic() + scale_fill_manual(values = c("#f8333c", "#007dba")) +
	xlab("") + ylab("Gene Count")

tmp <- depthDfTrans %>% filter(Gene %in% missingGenes, Genome == "KaeroEcoli")
tmp <- tmp %>% left_join(translatedNR[,c("Gene", "ID", "Name")])
tmp$Name <- sapply(1:nrow(tmp), function(x){ifelse(is.na(tmp$Name[x]),tmp$Gene[x] ,tmp$Name[x])})

tmp %>% select(Gene, MeanCoverage, PercentCoverage, CV) %>% as.data.frame() %>% xtable::xtable() %>% print(file = "~/Documents/University/EcoliPaperV2/MissingGenesWhole.tex")
write.table(tmp,file = "~/Documents/University/EcoliPaperV2/AdditionalFiles/MissingCoreGenesFixedWhole.tab", row.names = F, sep = "\t", quote = F)

depthDf %>% filter(Genome == "KaeroEcoli", Gene %in% coreGenes) %>% mutate(Presence = ifelse(MeanCoverage >= 10 & CV <= 1, "Present", "Absent")) %>%
	mutate(Presence = factor(Presence, levels = c("Present", "Absent"))) %>%
	ggplot(aes(x = MeanCoverage, y = CV, colour = Presence)) + geom_point(alpha = 0.5) + theme_bw() +
	geom_hline(yintercept = 1, col = "red", lty = 2) + geom_vline(xintercept = 10, col = "red", lty = 2) + xlab("Mean Read Coverage") +
	scale_y_continuous(breaks = pretty_breaks(n = 10)) + scale_x_continuous(breaks = pretty_breaks(n = 10)) +
       	ylab("Coefficient of Variation") + theme(legend.position = "bottom") + scale_colour_manual(values = list("Present" = "#3cb44b", "Absent" = "#e6194b"))

ggsave(file = "~/Documents/University/EcoliPaperV2/Figures/CoreGeneScatterWhole.pdf", width = 8, height = 6)
#####################################
#### Let's look at virulence now ####
#####################################

virDf <- AllPA %>% left_join(depthDfTrans %>% select(Gene, Name))# %>% mutate(Name = gsub(" \\[.*|MULTISPECIES: ", "", Name))
virDf <- virDf[,-removedgenomes] # Removing those outliers

# Using the list from Erick & Olivier
virGenes <- read.delim("Virulences_noms_genes.table", header = F, col.names = c("DB", "Gene", "Name", "Type")) %>% select(Gene, Type) %>% mutate(Gene = gsub("_\\d*|.*omal_|-.*|EDL.*|FT073|E2348|9\\.8","", Gene)) %>% as_tibble() %>% distinct() %>% pull(Gene)

## Now to get virulence genes identified and mapped
#databaseVir <- as_tibble(read.delim("Escherichia_VFs_comparsionV2.csv", header = T, skip = 1)) # Trying to get pathovar information
#databaseVir[,4:ncol(databaseVir)] <- sapply(databaseVir[,4:ncol(databaseVir)], function(x){ifelse(x == "", F, T)})
#genesUnlisted <- databaseVir[,3] %>% pull()%>% strsplit(split = "/") %>% unlist() %>% unique()
##
#virGenes <- databaseVir %>% pull(Related.genes) %>% unique() %>% strsplit(split = "/") %>% unlist() %>%
#       	gsub(pattern = "-.*",replacement = "") %>% unique()
#virGenes <- virGenes[-which(virGenes == "")]

# This is for the A0 List
virGenes <- strsplit(virGenes,split = "/") %>% unlist() %>% unique()
virGenes <- virGenes[grep("paa|^cfa$|^map$", virGenes, invert = T)]
virGenes <- gsub("GI*$", "G", virGenes)
virGenes <- gsub("espX.*", "espX", virGenes)
virGenes <- gsub("espL.*", "espL", virGenes)
virGenes <- gsub("espM.*", "espM", virGenes)
virGenes <- gsub("espR.*", "espR", virGenes)
virGenes <- virGenes[-which(virGenes == "69")]

##Getting them specifically now
#t3ss <- which(grepl("^sep|^esc|^esp|^ces|^tir|^epr|^eiv|^epa|^etr|^mxi|^tss|^hcp|^vgr", virDf$Name, ignore.case = T))
#upec <- which(grepl("^hly|^sfa|^usp|^pap|^sat|^ipa|^yad|^foc|^hof|^fim", virDf$Name, ignore.case = T)) #Morales-Espinosa et al 2016
#ehec <- which(grepl("^iha|^eae|^stx|^nle|^bfp|^ehx|^esp|^etp|^kat|^ent", virDf$Name, ignore.case = T)) #Bugarel et al 2011
#eaec <- which(grepl("^agg|^astA|^tss|^aaf|^AAF|^pet", virDf$Name, ignore.case = T)) #Croxen & Finlay 2010
#etec <- which(grepl("^est|^elt|^lta|^cfa|^tib|^etp|^eat|^cyl|east", virDf$Name, ignore.case = T))
#expec<- which(grepl("^pap|^prf|^sfa|^gaf|^bma|^iha|^afa|^tsh|^ibea|^irea|^iuc|^ybt|^iro|^sita|^hlya|^cdt|^cnf|^hlyf|^clb|^sat|^pic|^trat|^ompt|^iss|^iss|^cva|^dsda|^malx", virDf$Name, ignore.case =T))

# Now to get the new simpler method involved here as well
ind <- grep("invasion|secretion|enterotoxin|cfaE_2|cfaB_2|cfaB|porcine|elf|csg|hcp|ecp|AggR|AggA|Aap", virDf$Name) %>% unique()
testOut <- sapply(virGenes, function(x){grep(x, virDf$Gene)}) %>% unlist() %>% unique()
#virData <- virDf[unique(c(testOut,ind, t3ss,upec,ehec,eaec,etec,expec)),] %>% filter(!(grepl("esterase|fimbr|Hypothetical|group_.*|DUF", Name)))
virData <- virDf[unique(c(testOut,ind)),] %>% filter(!(grepl("esterase|Hypothetical|group_.*|DUF", Name)))

# We quickly want to identify those genes which have at least two EAEC Genes
tmp <- virData %>% filter(grepl("AggR|AggA|Aap|astA", Name)) %>% select(-Gene, -Name)%>% colSums() 
potentialEAEC <- names(tmp)[tmp > 1]

# Genes found in our ancient Genome
FoundGenes <- virData %>% select(Gene,Name, KaeroEcoli) %>% filter(KaeroEcoli == 1) %>% mutate(K12 = ifelse(Gene %in% k12Genes$Gene, T, F))
write.table(FoundGenes, file ="~/Documents/University/EcoliPaperV2/AdditionalFiles/WholePanGenomeVirulenceList.tab", row.names = F, sep = "\t")

# 
#virData <- virData %>% full_join(ancientData)
virData <- virData %>% mutate(Name = gsub(" \\[.*|MULTISPECIES: |, partial","",Name)) 
virData <- virData %>% select(-c(Gene)) %>% group_by(Name) %>% summarize_all(sum) %>% group_by(Name) %>% mutate_all(as.logical)
virData[,-1] <- sapply(virData[,-1], function(x){ifelse(is.na(x),0,x)})
virData[,-1] <- sapply(virData[,-1], function(x){ifelse(x,1,0)})

virData <- as.data.frame(virData)
rownames(virData) <- virData[,1]
virData <- virData[,-1]

#st4995Genomes <- st4995Genomes[st4995Genomes %in% colnames(virData)]
#tmp <- virData[,c(st4995Genomes[-2], "KaeroEcoli")] %>% filter(KaeroEcoli == 1)
# Limiting the genes to only those which have a difference
maxCount <- virData %>% rowSums() %>% max()
ind <- virData %>% rowSums() == maxCount
virData <- virData[!ind,]
virData <- as.data.frame(t(virData))
rownames(virData)[which(rownames(virData) == "KaeroEcoli")] <- "AncientEcoli"

## Getting the Gene Virulence Clusters
virDist <- dist(t(virData))
clusteredVir <- hclust(virDist, method = "ward.D2")

# Getting ready for the heatmap
specDist <- dist(virData)
clusteredSpec <- hclust(specDist, method = "ward.D2")
rownameOrder <- clusteredSpec$order
colnameOrder <- clusteredVir$order

ancientGenesFound <- as.logical(virData["AncientEcoli",])
ancientGenesFound <- data.frame(row.names = colnames(virData), "Ancient" = ifelse(ancientGenesFound, "Yes", "No"))

treeColours <- as.list(ann_colors$Pathovar)[virFrame[rownameOrder,]$Pathovar] %>% unlist()
treeColours <- gsub("#ffffff", "#000000", treeColours)
dend <- as.dendrogram(clusteredSpec)
labels_colors(dend) <- treeColours

#pdf("~/Documents/University/EcoliPaperV2/AdditionalFiles/VirulenceWholeTree.pdf", width = 24, height = 18)
#plot(dend, horiz = T)
#dev.off()

heatmapAnnotation <- virFrame[,-1]
#heatmapAnnotation <- data.frame("Phylogroup" = virFrame[,-1])
rownames(heatmapAnnotation) <- virFrame$Genome

pheatmap(virData[rownameOrder,], annotation_row = heatmapAnnotation, annotation_col = ancientGenesFound,
	 color = c("#007dba", "#f8333c"), clustering_method = "ward.D2", legend = F, show_colnames = F,
	 #gaps_row = c(which(rownames(virData)[rownameOrder] == "AncientEcoli"),which(rownames(virData)[rownameOrder] == "AncientEcoli") - 1),cluster_rows = F,
	 labels_row = rownamesHeatmap[rownameOrder], fontsize_row = 5, 
	 width = 12, height = 8, border_color = NA, annotation_colors = ann_colors,
	 filename = "~/Documents/University/EcoliPaperV2/AdditionalFiles/VirulenceWhole.png")#, main = "Core Presence/Absence")

##########################################################
### Now to get the mean coverages of the gene families ###
##########################################################
tmp <- depthDfTrans %>% filter(Name %in% FoundGenes$Name, MeanCoverage >= 10 & CV <= 1) %>% select(Name, MeanCoverage) # Need to remove some of the proteins which have identical names but aren't actually being found.....
write.table(tmp, "VirFoundBeforeEdits.tab",row.names = F,sep = "\t") # Manual Edits to the names will be easier here
tmp <- read.delim("VirFoundAfterEdits.tab", header =T) %>% as_tibble() 
tmp %>% group_by(GeneFamily) %>% summarize(Genes = length(GeneFamily), Mean = mean(MeanCoverage), SD = sd(MeanCoverage)) %>% summarize(sum(Genes))
tmp %>% group_by(GeneFamily) %>% summarize(Genes = length(GeneFamily), Mean = mean(MeanCoverage), SD = sd(MeanCoverage)) %>% filter(Genes > 1) %>% arrange(-Genes) %>%
       xtable:::xtable() %>% print(file = "~/Documents/University/EcoliPaperV2/AdditionalFiles/VirGeneFamiliesWhole.tex")

##############################################
### Testing if a PCoA will make more Sense ###
##############################################
geneDist <- dist(virData, method = "binary")
######################################
# If I want to use a MDS plot
fit <- cmdscale(geneDist, eig = T, k = 4, add = T)
coord <- fit$points %>% as_tibble()
coord$Genome <- rownames(fit$points)
contrib <- fit$eig/sum(fit$eig) * 100
coord <- coord %>% mutate(AncientLabel = ifelse(grepl("Ancient", Genome), "Ancient Ecoli", ""))
#coord <- coord %>% mutate(EAEC = ifelse(Genome %in% potentialEAEC, "Potentially EAEC", "Not EAEC"))

#olivierTableHeatmap$Genome <- rownames(olivierTableHeatmap)

VirPhylo <- coord %>% left_join(virFrame) %>%
       	ggplot(aes(x = V1, y = V2, label = AncientLabel, color = Phylogroup)) +
	geom_hline(yintercept=0, lty = 2, colour = "grey90") + 
	geom_vline(xintercept=0, lty = 2, colour = "grey90") +
	geom_point() +
	geom_text_repel(show.legend = F) +
	xlab(bquote("PCoA 1 ("~.(round(contrib[1],2))~"%)")) +
	ylab(bquote("PCoA 2 ("~.(round(contrib[2],2))~"%)")) +
	theme_classic() + scale_colour_manual(values = ann_colors$Phylogroup) +
	guides(color = guide_legend(nrow = 3))

virVir <- coord %>% left_join(virFrame) %>%
       	ggplot(aes(x = V1, y = V2, label = AncientLabel, color = Pathovar, alpha = Pathovar)) +
	geom_hline(yintercept=0, lty = 2, colour = "grey90") + 
	geom_vline(xintercept=0, lty = 2, colour = "grey90") +
	geom_point() +
	geom_text_repel(show.legend = F) +
	xlab(bquote("PCoA 1 ("~.(round(contrib[1],2))~"%)")) +
	ylab(bquote("PCoA 2 ("~.(round(contrib[2],2))~"%)")) +
	theme_classic() + scale_colour_manual(values = ann_colors$Pathovar) +
	scale_alpha_manual(values = alphaList$Pathovar) +
	guides(color = guide_legend(nrow = 3))
phyloData <- ggarrange(accessPhylo, VirPhylo, nrow = 1, labels = "AUTO", common.legend =T, legend = "bottom", align = "hv")
virPlot <- ggarrange(accessVir, virVir, nrow = 1, labels = c("C", "D"), common.legend =T, legend = "bottom", align = "hv")

ggarrange(phyloData,virPlot, nrow =2, common.legend = F, legend = "bottom")
ggsave("~/Documents/University/EcoliPaperV2/AdditionalFiles/VirulenceWholeMDS.pdf", width = 12, height = 9)


####### AMR Analysis #######
RGI <- as_tibble(read.delim("AMRFinal.txt", header = T)) %>% filter(grepl("homolog", Model_type))
tmp <- RGI %>% dplyr:::select(Best_Hit_ARO, Drug.Class, Resistance.Mechanism) %>% count(Best_Hit_ARO, name = "Count")
RGIHomo <- RGI %>% dplyr:::select(Best_Hit_ARO, Drug.Class, Resistance.Mechanism) %>% distinct() %>% left_join(tmp)

RGIMech <- RGIHomo %>% distinct()%>% pull(Resistance.Mechanism) %>% strsplit(";") %>% unlist() %>% table() %>% data.frame()
colnames(RGIMech) <- c("Mechanism", "Genes")
RGIMech <- RGIMech %>% arrange(-Genes)
RGIMech$Mechanism <- factor(RGIMech$Mechanism, levels = RGIMech$Mechanism)

p1 <- RGIMech %>% ggplot(aes(x = Mechanism, y = Genes)) + geom_col(fill = "#2e294e") + theme_classic() +
	theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

RGIFamily <- RGIHomo %>% filter(!grepl("efflux", Resistance.Mechanism)) %>% distinct() %>% pull(Drug.Class) %>% strsplit(";") %>% unlist() %>% table() %>% data.frame()
colnames(RGIFamily) <- c("Family", "Genes")
RGIFamily <- RGIFamily %>% arrange(-Genes)
RGIFamily$Family <- factor(RGIFamily$Family, levels = RGIFamily$Family)

p2 <- RGIFamily %>% ggplot(aes(x = Family, y = Genes)) + geom_col(fill = "#2e294e") + theme_classic() +
	theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
ggarrange(p1,p2, align = "h", labels = "AUTO")

RGIHomo %>% mutate(Drug.Class = gsub("; *", ";",Drug.Class)) %>% pull(Drug.Class) %>% strsplit(split = ";") %>% unlist() %>% table() %>% sort(decreasing = T) %>% length()

RGIHomo$Gene <- RGIHomo$Best_Hit_ARO %>% gsub(pattern = "Kleb.* |pneu.* |Escherichia |coli |beta-lactamase",replacement= "")

# Note, no coverage filtering on ancientGenesOnly because they all happen to be >= 10 in this case
#ancientGenesOnly <- depthDfTrans %>% filter(MeanCoverage >= 10, CV <=1)

# Finding as many AMR genes as possible
foundAMR <- lapply(split(RGIHomo, RGIHomo$Gene), function(x){
	       ind <- grep(x[5], ancientGenesOnly$Name, ignore.case = T)
	       tmp  <- ancientGenesOnly[ind,]
	       tmp$AMRGene <- x[[5]]
	       tmp$Drug.Class <- x[[2]] 
	       tmp$Resistance.Mechanism <- x[[3]] 
	       tmp$Count <- x[[4]]
	       return(tmp)
	 }) %>% bind_rows()

# Testing if these have been found in K12

# Saving the list of Found AMR genes
foundAMR %>% select(-c(Genome)) %>% arrange(-Count) %>% mutate(K12 = ifelse(Gene %in% k12Genes$Gene, T, F )) %>%
        write.csv("~/Documents/University/EcoliPaperV2/AdditionalFiles/AMRGenesFoundWhole.csv", row.names = F, quote = F)

# Now to figure our the resistances
tmp <- foundAMR %>% select(AMRGene, Drug.Class, Resistance.Mechanism, Count) %>% distinct() %>% 
	mutate(Drug.Class = gsub("; *", ";",Drug.Class)) %>% pull(Drug.Class) %>%
	strsplit(";") %>% unlist() %>% tibble()
colnames(tmp) <- "Resistances"
withefflux <- tmp %>% count(Resistances, name = "With Efflux")

# Removing Antibiotic Effluxes
tmp <- foundAMR %>% select(AMRGene, Drug.Class, Resistance.Mechanism, Count) %>% distinct() %>% 
	mutate(Resistance.Mechanism = gsub("antibiotic efflux","",Resistance.Mechanism),
		   Drug.Class = gsub("; *", ";",Drug.Class)) %>%
	filter(Resistance.Mechanism != "") %>%
	pull(Drug.Class) %>%
	strsplit(";") %>% unlist() %>% tibble()
colnames(tmp) <- "Resistances"
noefflux <- tmp %>% count(Resistances, name = "No Efflux")
withefflux %>% full_join(noefflux) %>% arrange(-`With Efflux`, Resistances)%>% filter(!is.na(`No Efflux`)) %>%
       xtable:::xtable() %>% print(file = "~/Documents/University/EcoliPaperV2/AdditionalFiles/DrugClassesWhole.tex")


#### TEMP GC #####

# Let's take alook at GC Content and Coverage
ancientGC <- ancientGenesPan %>% inner_join(geneGC)%>% filter(Genome == "PanGenomeMappingFeb", MeanCoverage >= 0)
ancientGC$Status<- ifelse(ancientGC$Gene %in% coreGenes, "Core", "Accessory")

p1 <- ancientGC %>% ggplot(aes(y = GC, x = MeanCoverage, colour = Status)) + geom_point(alpha = 0.5) + theme_bw() +
	scale_color_manual(values = c("#007dba", "#f8333c")) +
	geom_vline(xintercept = 10, color = "red", lty = 2) + ylab("GC Content") + xlab("Mean Read Depth") +
	scale_y_continuous(breaks = scales:::pretty_breaks(n = 10)) + theme(legend.position = "bottom")
#p1 <- ggMarginal(p1, margin = "x",type = "", groupFill = T)
p1

p2 <- ancientGC %>% ggplot(aes(y = GC, x = PercentCoverage, colour = Status)) + geom_point(alpha = 0.5) +
	theme_bw() + geom_vline(xintercept = 0.9, color = "red", lty = 2) + ylab("GC Content")+ xlab("Percent Coverage") +
	scale_color_manual(values = c("#007dba", "#f8333c")) + theme(legend.position = "bottom",axis.title.y = element_blank(), axis.text.y = element_blank()) 
#p2 <- ggMarginal(p2, type = "density", groupFill = T)
p2

p12 <- ggarrange(p1,p2, ncol = 2, align = "hv", labels = "AUTO", common.legend = T, legend = "bottom")
p12

ggsave("~/Documents/University/EcoliPaperV2/Figures/GeneDistributionGC.pdf", width = 8, height = 6)
