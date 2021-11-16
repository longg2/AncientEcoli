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

#HeterozygosityParsing <- function(Sample, CoreGenes, GeneTranslations){
#	# Preparing for the final part
#	FinalName <- gsub(".*\\/|\\.vcf", "",Sample)
#	print(FinalName)
#
#	# Getting the file Ready
#	vcfFile <- read.delim(file = Sample,comment.char = "#", header = F,
#			      col.names = c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "unknown")) %>% as_tibble()
#		#filter(REF != "C" & ALT != "T") %>% filter(REF != "G" & ALT != "A") %>% as_tibble() # Filtered out the G -> A and C -> T Transitions
#
#	if(nrow(vcfFile) == 0){
#		cat("No SNPs found in",FinalName, "\n")
#		return(NA)
#	}
#	# Translating Group Names into Actual Names
#	tmp2 <- vcfFile %>% filter(grepl("group", CHROM)) # The initial translations
#	tmp3 <- GeneTranslations[tmp2$CHROM]
#	tmp2$CHROM <- tmp3
#	tmp2$CHROM <- as.character(tmp2$CHROM)
#	
#	tmp21 <- vcfFile %>% filter(grepl("group", CHROM)) # Saving the groups for future reference
#	tmp2$CHROM[tmp2$CHROM == "NULL"] <- tmp21$CHROM[tmp2$CHROM == "NULL"]
#	
#	vcfFile <- vcfFile %>% filter(!(grepl(pattern = "group.*", CHROM))) %>% bind_rows(tmp2)
#	
#	### Getting the core Genome and the hetero Count from it
#	coreVcf <- vcfFile %>% filter(CHROM %in% CoreGenes)
#	if(nrow(coreVcf) == 0){
#		cat("No Core SNPs found in",FinalName, "\n")
#		return(NA)
#	}
#	
#	tmp <- unlist(strsplit(coreVcf$FORMAT[1], ":"))
#	tmp2 <- unlist(strsplit(gsub("=.*?(?=;)|=.*$","",coreVcf$INFO[1], perl = T), ";"))
#	coreVcf <- coreVcf %>% select(-FORMAT) %>% separate(unknown, into = tmp, sep = ":") %>%
#		mutate(INFO = gsub("(?<=;).*?=|^.*?=", "", INFO, perl = T)) %>% separate(INFO, into= tmp2, sep = ";")
#	coreVcf[,c(6,8:11,14)]<- sapply(coreVcf[,c(6,8:11,14)], as.numeric)
#	return(coreVcf)	
#	#final <- coreVcf %>% summarize(Type = ifelse(grepl("1/1|0/0",GT), "Homozygote", "Heterozygote")) %>% count(Type, name = "SNPs") %>% mutate(Sample = FinalName)
#	#return(final)
#}
#
#
colour <- c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6',
'#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3',
'#808000', '#ffd8b1', '#000075', '#808080', '#f0f0f0', '#000000')

#ann_colors = list(Pathovar = c("Ancient" = colour[1], "Commensal" = colour[2],"EAEC" = colour[4], "EIEC" = colour[5], "ETEC" = colour[6], "ExPEC" = colour[7],
#			       "Hybrid ExPEC/InPEC" = colour[8], "STEC" = colour[9],"?" = colour[3]),
#		  Host = c("Environment" = colour[10], "Food" = colour[11], "Human" = colour[12], "Livestock" = colour[13], "Porc" = colour[14], "Wild animal" = colour[15],"?" = colour[3]),
#		  ST = c("1429" = colour[16], "325" = colour[17], "3630" = colour[18], "399" = colour[19], "4995" = colour[20], "Other" = colour[3]),
#		  Ancient = c("Yes" = colour[22], "No" = colour[21]))

ann_colors = list(Pathovar = c("Ancient" = colour[1], "Commensal" = "#ffffff","EAEC" = colour[2], "EIEC" = "#22662A", "ETEC" = colour[16], "ExPEC" = colour[4],
			       "Hybrid ExPEC/InPEC" = colour[3], "STEC" = colour[9],"Unknown" = colour[20], "Kostic" = colour[15]),
		  #Host = c("Environment" = colour[10], "Food" = colour[11], "Human" = colour[12], "Livestock" = colour[13], "Porc" = colour[14], "Wild animal" = colour[15],"?" = colour[3]),
		  ST = c("1429" = colour[6], "325" = colour[7], "3630" = colour[8], "399" = colour[12], "4995" = colour[18], "Other" = colour[20], "Kostic" = colour[15], "Ancient" = colour[1]),
		  Ancient = c("Yes" = colour[22], "No" = colour[21]))


##### Only If rerunning the script ####
#load("Updated04012021.RData")
#load("Updated05022021MapOrig.RData")
########### Which strains are we actually plotting? #######
#strainList <- read.table("../TreemmerPlots/RAxML_bipartitions.GubbinsTreeKraEc_trimmed_list_RTL_0.95", header = F)[,1]
#strainList <- read.table("ContigsSnippyPrunedV2.list", header =F)[,1]
#strainList <- c(strainList, "PanGenomeDedupV3", "O104")

########### Pulling out the Phylogroups #######
phylogroup <- read.csv("VirulenceHeatmapOlivier.csv") %>% as_tibble()
 
#phylogroup <- read.delim("OlivierTable.tab", header = T)
etecTranslations <- read.delim("ETECOlivierTranslations.tab", header = F, col.names= c("Lineage", "GenomeID"))
etecList <- as.list(etecTranslations$Lineage)
names(etecList) <-  etecTranslations$GenomeID
tmpERR <- phylogroup %>% filter(grepl("ERR", Genome))
tmpERR$Genome <- etecTranslations$Lineage[match(tmpERR$Genome, etecTranslations$GenomeID)]
phylogroup <- phylogroup %>%filter(!grepl("ERR", Genome)) %>% bind_rows(tmpERR)
rm(tmpERR)

# What if we want the sources instead?
#olivierTable <- as_tibble(read.delim("PhyloMetadata.tab", header = T))
#olivierTable$Assembly.barcode <- sapply(olivierTable$Assembly.barcode, function(x){ifelse(grepl("^ETEC",x), gsub("ETEC", "ETEC ",x),x)})
#colnames(olivierTable)[1] <- "Genome"
#
#olivierTable$Pathovar <- sapply(olivierTable$Source, function(x){ifelse(grepl("Blood", x), "ExPEC", ifelse(grepl("Faeces",x), "Commensal",
#													 ifelse(grepl("Mummy",x), "Ancient", "Unknown")))})
#olivierTable$Pathovar <- sapply(1:length(olivierTable$ST), function(x){ifelse(grepl("ETEC", olivierTable$ST[x]), "ETEC",olivierTable$Pathovar[x])})
#olivierTable$Source <- sapply(olivierTable$Source, function(x){ifelse(grepl("animal|pig|,|Swine|Bovine|scrofa", x), "Animal",
#									    ifelse(grepl("Homo|ND|Blood|Faeces",x), "Human",
#										   ifelse(grepl("Cur|Cor|Ground", x), "Food", 
#											  ifelse(grepl("^$",x),"Unknown",x))))})
#
#olivierTable$Genome  <- sapply(olivierTable$Genome, function(x){ifelse(grepl("^ESC_", x), paste0(x,"_AS"),x)})
#olivierTable$Source  <- sapply(olivierTable$Source, function(x){ifelse(grepl("Human ", x),"Human",x)})
#olivierTable$Source  <- sapply(olivierTable$Source, function(x){ifelse(grepl("Mummy ", x),"Ancient",x)})

########### If we're cutting the heatmap in half #######
#filteredPhylogroup <- phylogroup %>% filter(Genome %in% strainList, !grepl("Shigella|Ancient|O104", Genome))
#preparation <- filteredPhylogroup %>% count(Phylogroup) %>% full_join(filteredPhylogroup, by = "Phylogroup") %>% mutate(Weight = 102/n)
#strainListV2  <- c(sample(preparation$Genome, size = nrow(preparation) * 0.5, prob = preparation$Weight), "O104")
#
#ind <- match(strainList, phylogroup$Genome)
#phylogroupHeat <- phylogroup$Phylogroup[ind]
#phylogroupHeat[length(phylogroupHeat)] <- "Ancient"
#
########### Pulling out the Pathovars #######
#pathovars <- as_tibble(read.csv("~/Downloads/PATRIC_genome.csv"))[,c("Genome.ID", "Pathovar")]
#colnames(pathovars)[1] <- "Genome"

############ Let's get the Depths ###########
files <- list.files(c("SmallPanGenomeDepths"), full.names = T)
#files <- files[452]
#files <- c(files,list.files("MappingDepths", full.names = T))
#files <- c(files,"PanGenomeDedupV3.tab.gz")[-452]
#files <- c(files, "KrakenDepths/GenusKrakenV3.tab.gz", "KrakenDepths/FamilyKrakenV3.tab.gz")
#files <- files[grep("562.28307", files)]

op <- pboptions(type = "timer")
ncores = 10
depthDf <- do.call(bind_rows, pblapply(files, cl = ncores,function(f){
		tmp <- as_tibble(read.delim(f, header = F, col.names = c("Gene", "Pos", "Coverage")))
		tmp$Genome <- gsub(".*/", "", gsub("\\..*|_genomic","",f)) 	
		tmp <- tmp %>% group_by(Genome, Gene) %>%
			summarize(MeanCoverage = mean(Coverage), sdCoverage = sd(Coverage), PercentCoverage = sum(Coverage > 0)/length(Coverage), .groups = "drop") %>%
			mutate(CV = sdCoverage/MeanCoverage)
		return(tmp)
	}))
depthDf <- depthDf %>% filter(!(Genome %in% c("IAI1", "ESC_NB8751AA_AS", "ATCC11231", "ESC_LB2165AA_AS")))
depthDf$Genome <- unname(unlist(sapply(depthDf$Genome, function(x){ifelse(grepl("^ERR",x), etecList[x],x)})))
#unique(depthDf$Genome)
#depthDf %>% filter(!(Genome %in% c("PanGenomeDepths","PanGenomeDepthsDup")))

# Preparing for the rest of the script here
strainListV2 <- unique(depthDf$Genome)
strainListOlivier <- strainListV2[!grepl("ETEC",strainListV2)]
strainListFragmented <- c(strainListV2[grepl("ETEC",strainListV2)], "KaeroEcoli") 

############# Core Gene Presence #####
roaryOutput <- as_tibble(read.delim("SmallGenomePresenceAbsence.tab"))
colnames(roaryOutput)[2:ncol(roaryOutput)] <- gsub("X|\\.scaffold|\\.genome|\\.result|_genomic", "", colnames(roaryOutput)[2:ncol(roaryOutput)])
roaryOutput <- roaryOutput %>% select(-c(IAI1, ATCC11231,ESC_LB2165AA_AS))
# Finding the Core Genome
ecoliOnly <- roaryOutput
coreGenes <- rowSums(ecoliOnly[,-1]) >= floor(0.95 * ncol(ecoliOnly))
#coreGenes <- rowSums(ecoliOnly[,-1]) >= floor(0.99 * ncol(ecoliOnly))
coreGenes <- roaryOutput$Gene[coreGenes]


########### Figuring out how far our ancient strain is ###############
st4995 <- phylogroup[,1:2] %>% filter(ST == 4995) %>% pull(Genome)
st4995 <- st4995[st4995 %in% colnames(roaryOutput)]

ancientData <- depthDf %>% filter(Genome == "KaeroEcoli", MeanCoverage >= 10, PercentCoverage >= 0.9) %>%
	pivot_wider(names_from = c(Genome), values_from=MeanCoverage) %>% select(-c(sdCoverage, PercentCoverage, CV))
#####################################################################
# Gene Presence table

AllPA <- roaryOutput %>% left_join(ancientData) %>% mutate(KaeroEcoli = ifelse(is.na(KaeroEcoli)|KaeroEcoli == 0, 0,1))

AllPA %>% select(Gene,ATCC11229, KaeroEcoli) %>% filter(ATCC11229 == 1 | KaeroEcoli == 1) %>%
       	write.table(file = "~/Documents/University/EcoliPaperV2/AdditionalFiles/GenePAATCCAncient.tab",row.names = F,sep = "\t", quote = F)

tmp <- AllPA %>% select(Gene, ATCC11229, KaeroEcoli, ESC_VA4573AA_AS) %>%
       	mutate(ATCC11229 = ifelse(ATCC11229 == 1, Gene, NA),
		      KaeroEcoli = ifelse(KaeroEcoli == 1, Gene, NA),
		      ESC_VA4573AA_AS = ifelse(ESC_VA4573AA_AS == 1, Gene, NA))
vennList <- list("ATCC11229" = tmp$ATCC11229[!is.na(tmp$ATCC11229)],
		 "FSIS11816402" = tmp$ESC_VA4573AA_AS[!is.na(tmp$ESC_VA4573AA_AS)],
		 "Ancient" = tmp$KaeroEcoli[!is.na(tmp$KaeroEcoli)])

ggvenn(vennList, fill_color = c("#2e294e", "#f8333c", "#007dba"))
ggsave("~/Documents/University/EcoliPaperV2/Figures/VennAncientFSISATCC.pdf", width = 7, height = 7)


#####################################################################

# Going to do this in a more sane way
coreDf <- roaryOutput %>% filter(Gene %in% coreGenes)
tmp <- colnames(coreDf)
tmp <- sapply(tmp, function(x){
	       ifelse(x %in% etecTranslations$GenomeID, etecList[x], x)
		  }) %>% unlist() %>% gsub(pattern = "\\.result.*|_genomic|\\.scaffold", replacement = "") %>%
	gsub(pattern = "\\.",replacement = "-")
colnames(coreDf) <- tmp

ancientData <- depthDf %>% filter(Genome %in%  c("KaeroEcoli"), MeanCoverage >= 10, PercentCoverage >= 0.9) %>% mutate(Status = ifelse(Gene %in% coreGenes, "Core", "Accessory")) %>% filter(Status == "Core") %>%
	pivot_wider(names_from = c(Genome), values_from=MeanCoverage) %>% select(-c(sdCoverage, PercentCoverage, CV, Status))

#KosticData <- depthDf %>% filter(Genome %in% c("Lib4_3_bin","Lib4_7_bin","Lib4_8_bin"), MeanCoverage >= 1, PercentCoverage >= 0.9) %>%
#       	mutate(Status = ifelse(Gene %in% coreGenes, "Core", "Accessory")) %>%
#       	filter(Status == "Core") %>%
#	pivot_wider(names_from = c(Genome), values_from=MeanCoverage) %>%
#       	select(-c(sdCoverage, PercentCoverage, CV, Status))
#
#KosticData <- KosticData %>% group_by(Gene) %>% summarize_all(list(~ sum(., na.rm = T))) %>% arrange(Gene)
#
#KosticData %>% select(-Gene) %>% summarize_all(list(~ ifelse(. > 0, 1,.))) %>% summarize_all(sum)
#
#ancientData <- ancientData %>% full_join(KosticData) %>%
#	group_by(Gene) %>% summarize_all(list(~ sum(., na.rm = T))) %>%
#	ungroup() %>%
#       	mutate(KaeroEcoli = ifelse(is.na(KaeroEcoli) | KaeroEcoli == 0, 0, 1),
#		Lib4_3_bin = ifelse(is.na(Lib4_3_bin)| Lib4_3_bin == 0, 0, 1),
#		Lib4_7_bin = ifelse(is.na(Lib4_7_bin)| Lib4_7_bin == 0, 0, 1),
#		Lib4_8_bin = ifelse(is.na(Lib4_8_bin)| Lib4_8_bin == 0, 0, 1))# %>%
##	select(-Gene) %>% summarize_all(list(~ ifelse(. > 0, 1,.))) %>% summarize_all(sum)

coreData <- coreDf %>% left_join(ancientData) %>%
       	mutate(KaeroEcoli = ifelse(is.na(KaeroEcoli) | KaeroEcoli == 0, 0, 1))

# Here, we're trying to plot the core P/A data as a heatmap and/or a hierachical cluster
#coreData <- depthDf %>% filter(Genome %in% strainListOlivier, MeanCoverage >= 1, PercentCoverage >= 0.9) %>% mutate(Status = ifelse(Gene %in% coreGenes, "Core", "Accessory")) %>% filter(Status == "Core") %>%
#	pivot_wider(names_from = c(Genome), values_from=MeanCoverage) %>% select(-c(sdCoverage, PercentCoverage, CV, Status))

#ancientData <- depthDf %>% filter(Genome %in% strainListFragmented, MeanCoverage >= 10, PercentCoverage >= 0.9) %>% mutate(Status = ifelse(Gene %in% coreGenes, "Core", "Accessory")) %>% filter(Status == "Core") %>%
#	pivot_wider(names_from = c(Genome), values_from=MeanCoverage) %>% select(-c(sdCoverage, PercentCoverage, CV, Status))
#coreData <- coreData %>% full_join(ancientData)
#
#coreData[,-1] <- sapply(coreData[,-1], function(x){ifelse(is.na(x), 0, x)})
#coreData <- coreData %>% group_by(Gene) %>% summarize_all(sum) %>% group_by(Gene) %>% mutate_all(as.logical)

# Preparing the data
coreData <- as.data.frame(coreData)
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
olivierTableHeatmap <- phylogroup[,c(4,2)] %>% as.data.frame()
olivierTableHeatmap$Pathovar <- sapply(olivierTableHeatmap$Pathovar, function(x){ifelse(x == "?", "Unknown",x)})
rownames(olivierTableHeatmap) <- phylogroup$Genome
olivierTableHeatmap <- olivierTableHeatmap[rownames(coreDataSub),]
#olivierTableHeatmap <- olivierTableHeatmap[complete.cases(olivierTableHeatmap),] %>%
#       	bind_rows(data.frame(row.names = c("Zape2.1","Zape2.2","Zape3"),ST = c("Kostic", "Kostic", "Kostic"), Pathovar = "Kostic"))
olivierTableHeatmap["AncientEcoli","ST"] <- "Ancient"

rownamesHeatmap <- ifelse(rownames(coreDataSub) == "AncientEcoli", "AncientEcoli","")
#rownames(olivierTableHeatmap) <- phylogroup$Genome
#rownames(olivierTableHeatmap) <- olivierTable$Genome

# Dendrogram and order
specDist <- dist(coreDataSub, method = "binary")
clusteredSpec <- hclust(specDist, method = "ward.D2") # 
rownameOrder <- clusteredSpec$order

treeColours <- as.list(ann_colors[[2]])[olivierTableHeatmap[rownameOrder,]$ST] %>% unlist()
treeColours <- gsub("#ffffff", "#000000", treeColours)
dend <- as.dendrogram(clusteredSpec)
labels_colors(dend) <- treeColours

pdf("~/TMPCoreTree.pdf", width = 24, height = 9)
#pdf("~/Documents/University/EcoliPaperV2/AdditionalFiles/PartsofFigures/SpeciesCoreTree.pdf", width = 6, height = 8)
plot(dend, horiz = F)
legend("topright", legend = names(ann_colors[[2]]), fill = unlist(ann_colors[[2]]), title = "ST")
dev.off()

# Now to cluster the trees
geneDist <- dist(coreDataSub)
geneClust <- hclust(geneDist, method = "ward.D2")

sil <- sapply(2:15, function(x){
	       tmp <- cutree(geneClust, k = x)
   	       summary(silhouette(tmp, geneDist))$avg.width
	 })

plot(y = sil, x = 2:15, type = "b")

small <- cutree(geneClust, k = 5)
clusterBorder <- colour[small] %>% unique()

pdf("~/KosticHierCoreClust.pdf", width = 12, height = 9)
plot(geneClust, hang = -1)
#rect.hclust(geneClust, k = 4, border = clusterBorder)
dev.off()

pheatmap(coreDataSub[rownameOrder,], annotation_row = olivierTableHeatmap,
	 color = c("#007dba", "#f8333c"), clustering_method = "ward.D2", legend = F, show_colnames = F,
	 gaps_row = c(which(rownames(coreDataSub)[rownameOrder] == "AncientEcoli"),which(rownames(coreDataSub)[rownameOrder] == "AncientEcoli") - 1),cluster_rows = F,
	 labels_row = rownamesHeatmap[rownameOrder], fontsize_row = 5, 
	 width = 12, height = 8, border_color = NA, annotation_colors = ann_colors)#, main = "Core Presence/Absence")
	 #filename = "~/Documents/University/EcoliPaperV2/AdditionalFiles/PartsofFigures/TestCore.png", width = 12, height = 8, border_color = NA, annotation_colors = ann_colors)#, main = "Core Presence/Absence")
#pheatmap(coreDataSub[rownameOrder,], annotation_row = olivierTableHeatmap,
#	 color = c("#007dba", "#f8333c"), clustering_method = "ward.D2", legend = F, show_colnames = F,
#	 gaps_row = c(which(rownames(coreDataSub)[rownameOrder] == "AncientEcoli"),which(rownames(coreDataSub)[rownameOrder] == "AncientEcoli") - 1),cluster_rows = F,
#	 cutree_cols = 4, labels_row = rownamesHeatmap[rownameOrder], fontsize_row = 5, 
#	 filename = "~/TestCore.png", width = 12, height = 8, border_color = NA, annotation_colors = ann_colors)#, main = "Core Presence/Absence")

# Core data counts for phylogeny
counts <- coreData %>% rowSums()
counts <- tibble(Genome = names(counts), Count = counts)
write.table(counts, sep = "\t", "../CoreGeneCounts.tab")

########## Accessory Genome Analysis ########
accessDf <- roaryOutput %>% filter(!(Gene %in% coreGenes))
tmp <- colnames(accessDf)
tmp <- sapply(tmp, function(x){
	       ifelse(x %in% etecTranslations$GenomeID, etecList[x], x)
		  }) %>% unlist() %>% gsub(pattern = "\\.result.*|_genome|_genomic|\\.scaffold", replacement = "") %>%
	gsub(pattern = "\\.",replacement = "-")
colnames(accessDf) <- tmp
#colnames(accessDf) %>% filter(!(Genome %in% c("IAI1", "ESC_NB8751AA_AS", "ATCC11231", "ESC_LB2165AA_AS")))
ancientData <- depthDf %>% filter(Genome == "KaeroEcoli", MeanCoverage >= 10) %>% mutate(Status = ifelse(Gene %in% coreGenes, "Core", "Accessory")) %>% filter(Status == "Accessory") %>%
	pivot_wider(names_from = c(Genome), values_from=MeanCoverage) %>% select(-c(sdCoverage, PercentCoverage, CV, Status))

KosticData <- depthDf %>% filter(Genome %in% c("Lib4_3_bin","Lib4_7_bin","Lib4_8_bin"), MeanCoverage >= 1, PercentCoverage >= 0.9) %>%
       	mutate(Status = ifelse(Gene %in% coreGenes, "Core", "Accessory")) %>%
       	filter(Status == "Accessory") %>%
	pivot_wider(names_from = c(Genome), values_from=MeanCoverage) %>%
       	select(-c(sdCoverage, PercentCoverage, CV, Status))

KosticData <- KosticData %>% group_by(Gene) %>% summarize_all(list(~ sum(., na.rm = T))) %>% arrange(Gene)

ancientData <- ancientData %>% full_join(KosticData) %>%
	group_by(Gene) %>% summarize_all(list(~ sum(., na.rm = T))) %>%
	ungroup() %>%
       	mutate(KaeroEcoli = ifelse(is.na(KaeroEcoli) | KaeroEcoli == 0, 0, 1),
		Lib4_3_bin = ifelse(is.na(Lib4_3_bin)| Lib4_3_bin == 0, 0, 1),
		Lib4_7_bin = ifelse(is.na(Lib4_7_bin)| Lib4_7_bin == 0, 0, 1),
		Lib4_8_bin = ifelse(is.na(Lib4_8_bin)| Lib4_8_bin == 0, 0, 1))# %>%

accessData <- accessDf %>% left_join(ancientData) %>%# mutate(KaeroEcoli = ifelse(is.na(KaeroEcoli), 0, 1))
       	mutate(KaeroEcoli = ifelse(is.na(KaeroEcoli) | KaeroEcoli == 0, 0, 1),
		Lib4_3_bin = ifelse(is.na(Lib4_3_bin)| Lib4_3_bin == 0, 0, 1),
		Lib4_7_bin = ifelse(is.na(Lib4_7_bin)| Lib4_7_bin == 0, 0, 1),
		Lib4_8_bin = ifelse(is.na(Lib4_8_bin)| Lib4_8_bin == 0, 0, 1))

# Now, we'll be doing the same thing for the accessory genome.  Firstly, however, is going to be the removal genes which are barely present.  We'll do a separate one with unique genes
#accessData <- depthDf %>% filter(Genome %in% strainListOlivier, MeanCoverage >= 1) %>% mutate(Status = ifelse(Gene %in% coreGenes, "Core", "Accessory")) %>% filter(Status == "Accessory") %>%
	#pivot_wider(names_from = c(Genome), values_from=MeanCoverage) %>% select(-c(sdCoverage, PercentCoverage, CV, Status))
#ancientData <- depthDf %>% filter(Genome %in% strainListFragmented, MeanCoverage >= 10) %>% mutate(Status = ifelse(Gene %in% coreGenes, "Core", "Accessory")) %>% filter(Status == "Accessory") %>%
#	pivot_wider(names_from = c(Genome), values_from=MeanCoverage) %>% select(-c(sdCoverage, PercentCoverage, CV, Status))
#
#accessData <- accessData %>% full_join(ancientData)

accessData[,-1] <- sapply(accessData[,-1], function(x){ifelse(is.na(x), 0, x)})
accessData <- accessData %>% group_by(Gene) %>% summarize_all(sum) %>% group_by(Gene) %>% mutate_all(as.logical)

# Preparing the data
accessData <- as.data.frame(accessData)
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
rownames(accessDataSub)[which(rownames(accessDataSub) %in% c("KaeroEcoli","Lib4_3_bin","Lib4_7_bin","Lib4_8_bin"))] <- c("AncientEcoli", "Zape2.1", "Zape2.2", "Zape3")
rownamesHeatmap <- ifelse(rownames(accessDataSub) == "AncientEcoli", "AncientEcoli","")
olivierTableHeatmap <- phylogroup[,c(4,2)] %>% as.data.frame()
olivierTableHeatmap$Pathovar <- sapply(olivierTableHeatmap$Pathovar, function(x){ifelse(x == "?", "Unknown",x)})
rownames(olivierTableHeatmap) <- phylogroup$Genome
olivierTableHeatmap <- olivierTableHeatmap[rownames(accessDataSub),]
olivierTableHeatmap <- olivierTableHeatmap[complete.cases(olivierTableHeatmap),] %>%
       	bind_rows(data.frame(row.names = c("Zape2.1","Zape2.2","Zape3"),ST = c("Kostic", "Kostic", "Kostic"), Pathovar = "Kostic"))
olivierTableHeatmap["AncientEcoli","ST"] <- "Ancient"

# Dendrogram and order
specDist <- dist(accessDataSub, method = "binary")
clusteredSpec <- hclust(specDist, method = "ward.D2") # 
rownameOrder <- clusteredSpec$order

treeColours <- as.list(ann_colors[[2]])[olivierTableHeatmap[rownameOrder,]$ST] %>% unlist()
treeColours <- gsub("#ffffff", "#000000", treeColours)
dend <- as.dendrogram(clusteredSpec)
labels_colors(dend) <- treeColours

pdf("~/KosticAccessTree.pdf", width = 24, height = 18)
#pdf("~/Documents/University/EcoliPaperV2/AdditionalFiles/PartsofFigures/SpeciesAccessTree.pdf", width = 6, height = 8)
plot(dend, horiz = F, main = "Accessory Genome Clustering")
legend("topright", legend = names(ann_colors[[2]]), fill = unlist(ann_colors[[2]]), title = "ST")
dev.off()

# Now to cluster the trees
#geneDist <- dist(t(accessDataSub))
geneDist <- dist(accessDataSub, method = "binary")
geneClust <- hclust(geneDist, method = "ward.D2")
######################################
# If I want to use a MDS plot
fit <- cmdscale(geneDist, eig = T, k = 4)
coord <- fit$points %>% as_tibble()
coord$Genome <- rownames(fit$points)

olivierTableHeatmap$Genome <- rownames(olivierTableHeatmap)

p1 <- coord %>% left_join(olivierTableHeatmap) %>%
       	ggplot(aes(x = V1, y = V2, label = Genome, color = ST)) +
	geom_hline(yintercept=0, lty = 2, colour = "grey90") + 
	geom_vline(xintercept=0, lty = 2, colour = "grey90") +
	geom_point() +
	geom_text_repel(show.legend = F) +
	theme_classic() + scale_colour_manual(values = ann_colors[["ST"]]) 

p2 <- coord %>% left_join(olivierTableHeatmap) %>%
       	ggplot(aes(x = V3, y = V4, label = Genome, color = ST)) +
	geom_hline(yintercept=0, lty = 2, colour = "grey90") + 
	geom_vline(xintercept=0, lty = 2, colour = "grey90") +
	geom_point() +
	geom_text_repel(show.legend = F) +
	theme_classic() + scale_colour_manual(values = ann_colors[["ST"]]) 

ggarrange(p1,p2, legend = "bottom", align = "hv", ncol = 1, common.legend = T)
#####################################

sil <- sapply(2:15, function(x){
	       tmp <- cutree(geneClust, k = x)
   	       summary(silhouette(tmp, geneDist))$avg.width
	 })

plot(y = sil, x = 2:15, type = "b")

small <- cutree(geneClust, k = 4)
clusterBorder <- colour[small] %>% unique()

pheatmap(accessDataSub[rownameOrder,], annotation_row = olivierTableHeatmap,
	 color = c("#007dba", "#f8333c"), clustering_method = "ward.D2", legend = F, show_colnames = F,
	 gaps_row = c(which(rownames(accessDataSub)[rownameOrder] == "AncientEcoli"),which(rownames(accessDataSub)[rownameOrder] == "AncientEcoli") - 1),cluster_rows = F,
	labels_row = rownamesHeatmap[rownameOrder], fontsize_row = 5, 
	 #border_color = NA, annotation_colors = ann_colors)#, main = "Accessory Presence/Absence")
	 filename = "~/AccessoryKosticHeatmap.png", width = 12, height = 8, border_color = NA, annotation_colors = ann_colors)#, main = "Accessory Presence/Absence")

################################################
### Now to perform the presence absence Test ###
################################################
#ancientGenes <- ancientGenesPan %>% filter(Genome == "PanGenomeMappingFeb", MeanCoverage > 0, PercentCoverage > 0.5)
ancientGenes <- depthDf %>% filter(Genome == "KaeroEcoli", MeanCoverage >= 10, PercentCoverage >= 0.9)
corePA <- coreGenes %in% ancientGenes$Gene
missingGenes <- coreGenes[!corePA]

status <- ifelse(corePA, "Core Found", "Core Missing")
status <- tibble(status) %>% count(status)
colnames(status)[1] <- "Status"
status
status %>% ggplot(aes(x = "", y = `n`, fill = Status)) + geom_col() + theme_classic() + scale_fill_manual(values = c("#f8333c", "#007dba")) +
	xlab("") + ylab("Gene Count")
ggsave("~/Documents/University/LabMeetings/Feb12BEAP/Figures/UpdatedCoreGenes.pdf", height = 6, width = 4)

tmp <- depthDf %>% filter(Gene %in% missingGenes, Genome == "KaeroEcoli")
tmp <- tmp %>% left_join(translatedNR[,c("Gene", "ID", "Name")])
tmp$Name <- sapply(1:nrow(tmp), function(x){ifelse(is.na(tmp$Name[x]),tmp$Gene[x] ,tmp$Name[x])})

tmp %>% select(Gene, MeanCoverage, PercentCoverage, CV) %>% as.data.frame() %>% xtable::xtable() %>% print(file = "~/Documents/University/EcoliPaperV2/MissingGenes.tex")
write.table(tmp,file = "~/Documents/University/EcoliPaperV2/AdditionalFiles/MissingCoreGenesFixed.tab", row.names = F, sep = "\t", quote = F)

depthDf %>% filter(Genome == "KaeroEcoli", Gene %in% coreGenes) %>% mutate(Presence = ifelse(MeanCoverage >= 10 & PercentCoverage >= 0.9, "Present", "Absent")) %>%
	mutate(Presence = factor(Presence, levels = c("Present", "Absent"))) %>%
	ggplot(aes(x = MeanCoverage, y = PercentCoverage, colour = Presence)) + geom_point(alpha = 0.5) + theme_bw() +
	geom_hline(yintercept = 0.9, col = "red", lty = 2) + geom_vline(xintercept = 10, col = "red", lty = 2) + xlab("Mean Read Coverage") +
	scale_y_continuous(breaks = pretty_breaks(n = 10)) + scale_x_continuous(breaks = pretty_breaks(n = 10)) +
       	ylab(expression("Nucleotide Coverage" >= 1)) + theme(legend.position = "bottom") + scale_colour_manual(values = list("Present" = "#3cb44b", "Absent" = "#e6194b"))

ggsave(file = "~/Documents/University/EcoliPaperV2/Figures/CoreGeneScatter.pdf", width = 8, height = 6)


#ancientGenes[corePA,]

#write.table(coreGenes[!corePA], file = "~/MissingCoreGenes.list", row.names = F, quote = F, col.names = F)
##modernPA <- coreGenes %in% gffNoCopy
##status <- ifelse(modernPA, "Core Found", "Core Missing")
##status <- tibble(status) %>% count(status)
#
##### Let's get those Density Plots going ####
#ancientGenesPan <- depthDf %>% filter(Genome == "MQ30ML30Depths", PercentCoverage >= 0.9)
ancientGenesOnly <- depthDf %>% filter(Genome == "KaeroEcoli", MeanCoverage >= 10, PercentCoverage >= 0.9)
ancientGenesOnly$Status <- ifelse(ancientGenesOnly$Gene %in% coreGenes, "Core", "Accessory")

# Quick summary Statistics
ancientGenesOnly %>% summarize(Mean = mean(MeanCoverage), SD = sd(MeanCoverage), confInt = qnorm(0.975)*SD/sqrt(length(MeanCoverage)), high = Mean + confInt, low = Mean - confInt) %>% as.data.frame()
ancientGenesOnly %>% group_by(Status)%>% summarize(Mean = mean(MeanCoverage), SD = sd(MeanCoverage), confInt = qnorm(0.975)*SD/sqrt(length(MeanCoverage)), high = Mean + confInt, low = Mean - confInt) %>% as.data.frame()

p1 <- ancientGenesOnly %>% ggplot(aes(x = MeanCoverage)) + geom_histogram(fill = "#2e294e", colour = "black", bins = 60) +
	scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
	ylab("Genes") + theme_bw() + xlab("") + theme(axis.text.x = element_blank()) + scale_y_continuous(breaks = scales::pretty_breaks(n = 10), limits = c(0, 550))

breaks <- pretty(range(ancientGenesOnly$MeanCoverage), n = nclass.FD(ancientGenesOnly$MeanCoverage), min.n = 1)
bwidth <- breaks[2]-breaks[1]

ancientGenesOnly <- ancientGenesOnly %>% mutate(Status = factor(Status, levels = c("Core", "Accessory")))
p2 <- ancientGenesOnly  %>% ggplot(aes(x = MeanCoverage, fill = Status)) + geom_histogram(position = "identity", alpha = 0.75, colour = "black", binwidth = 1) +
	#scale_fill_manual(values = c(Core = "#f8333c", Accessory = "#007dba", `K. aerogenes` = "#1b998b")) + theme_bw() +
	scale_fill_manual(values = c(Core = "#f8333c", Accessory = "#007dba")) + theme_bw() +
	xlab("Mean Read Depth") + ylab("Genes") + theme(legend.position = "bottom") + 
	scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) + scale_y_continuous(breaks = scales::pretty_breaks(n = 10), limits = c(0, 500))
p2
save(p2, file = "~/GeneHistogram.RData")
ggsave("~/UpdatedGeneDepthHistogram.pdf", width = 6, height = 4)

#ggarrange(p1, p2, nrow = 2, align = "hv", labels = "AUTO")
#ggsave("~/Documents/University/EcoliPaperV2/Figures/GeneDistributionUpdatedPercent90.pdf", width = 8, height = 6)

#####################################
#### Let's look at virulence now ####
#####################################

# First, we will want to replace the group names when possible
translatedNR <- as_tibble(read.delim("../ClermontST4995AmbigGenes/SmallPanGenomeGoodHits.tab", header = T))
colnames(translatedNR)[c(1,2,14)] <- c("Gene", "ID", "Name")
depthDfTrans <- depthDf %>% left_join(translatedNR[, c("Gene", "ID", "Name")])

depthDfTrans$Name <- sapply(1:nrow(depthDfTrans), function(x){ifelse(is.na(depthDfTrans$Name[x]),depthDfTrans$Gene[x] ,depthDfTrans$Name[x])})

depthDfTrans %>% filter(Genome == "KaeroEcoli", Gene %in% missingGenes) %>%
       	write.table(file = "~/Documents/University/EcoliPaperV2/AdditionalFiles/MissingGenesWithCov.tab", row.names = F, quote = F, sep = "\t")

virDf <- roaryOutput
tmp <- colnames(virDf)
tmp <- sapply(tmp, function(x){
	       ifelse(x %in% etecTranslations$GenomeID, etecList[x], x)
		  }) %>% unlist() %>% gsub(pattern = "\\.result.*|_genome|_genomic|\\.scaffold", replacement = "") %>%
	gsub(pattern = "\\.",replacement = "-")
colnames(virDf) <- tmp

virDf <- virDf %>% left_join(translatedNR[, c("Gene", "ID", "Name")])
virDf$Name <-  sapply(1:nrow(virDf), function(x){ifelse(is.na(virDf$Name[x]),virDf$Gene[x] ,virDf$Name[x])})

##########Getting some shared Gene Counts################################
ancientData <- depthDfTrans %>% filter(Genome == "KaeroEcoli", MeanCoverage >= 10, PercentCoverage >= 0.9) %>%
	pivot_wider(names_from = c(Genome), values_from=MeanCoverage) %>% select(-c(sdCoverage, PercentCoverage, CV, Name, ID))

tmp <- virDf %>% select(-c(Name, ID)) %>% left_join(ancientData) %>% mutate(KaeroEcoli = ifelse(is.na(KaeroEcoli),0,1))

tmp <- tmp %>% select(c(Gene,st4995, "KaeroEcoli"))
# Now to make the list
st4995Venn <- list()
st4995Venn$`Ancient Ecoli` <- tmp$Gene[as.logical(tmp$KaeroEcoli)]
st4995Venn$ESC_VA4573AA_AS <- tmp$Gene[as.logical(tmp$ESC_VA4573AA_AS)]
st4995Venn$ATCC11229 <- tmp$Gene[as.logical(tmp$ATCC11229)]
st4995Venn$ESC_ZA5349AA_AS <- tmp$Gene[as.logical(tmp$ESC_ZA5349AA_AS)]

ggvenn(st4995Venn)
ggsave(file = "~/Documents/University/EcoliPaperV2/AdditionalFiles/PartsofFigures/ST4995Venn.pdf", width = 9, height = 6)

# What are the missing genes?
tmp2 <- tmp %>% select(-Gene)
ind <- tmp2 %>% rowSums() %>% as.logical()
tmp2 <- tmp[ind,] %>% filter(KaeroEcoli == 0)
ind <- tmp2[,-1] %>% rowSums() > 1
missingST4995Genes <- tmp2[ind,"Gene"] %>% pull()
virDf %>% filter(Gene %in% missingST4995Genes) %>% select(Gene, ID, Name) %>% write.table(file = "~/Documents/University/EcoliPaperV2/AdditionalFiles/ST4995Missing.tab", row.names = F, quote = F, sep = "\t")
ancientOnly <- tmp[ind,] %>% filter(KaeroEcoli == 1, ATCC11229 == 0, ESC_ZA5349AA_AS == 0, ESC_VA4573AA_AS == 0) %>% pull(Gene)
virDf %>% filter(Gene %in% ancientOnly) %>% select(Gene, ID, Name) %>% write.table(file = "~/Documents/University/EcoliPaperV2/AdditionalFiles/ST4995AncientOnly.tab", row.names = F, quote = F, sep = "\t")
##########################################

#ancientData <- depthDf %>% filter(Genome == "KaeroEcoli", MeanCoverage >= 10, PercentCoverage >= 0.9) %>% mutate(Status = ifelse(Gene %in% coreGenes, "Core", "Accessory")) %>% filter(Status == "Core") %>%
#	pivot_wider(names_from = c(Genome), values_from=MeanCoverage) %>% select(-c(sdCoverage, PercentCoverage, CV, Status))
#
#coreData <- coreDf %>% left_join(ancientData) %>% mutate(KaeroEcoli = ifelse(is.na(KaeroEcoli), 0, 1))

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

## Accounting for the NR genes.  Not clear so we'll need to start with key words (secretion, invasion, virulence, toxin)
ind <- grep("invasion|secretion|enterotoxin", depthDfTrans$Name) %>% unique()
#

# The old way with all being mapped
testOut <- pbsapply(virGenes, cl = 10, function(x){grep(x, depthDfTrans$Gene)}) %>% unlist() %>% unique()
foundGenes <- depthDfTrans$Gene[testOut] %>% unique()
ancientData <- depthDfTrans %>% ungroup() %>% filter(Gene %in% foundGenes) %>% bind_rows(depthDfTrans[ind,]) %>%
       	filter(!is.na(CV), !grepl("group_",Name)) %>% distinct() %>% filter(Genome == "KaeroEcoli", MeanCoverage >= 10, PercentCoverage >= 0.9) %>%
	select(Gene, Name, MeanCoverage)
colnames(ancientData)[3] <- "KaeroEcoli"

# Now to get the new simpler method involved here as well
ind <- grep("invasion|secretion|enterotoxin", virDf$Name) %>% unique()
testOut <- sapply(virGenes, function(x){grep(x, virDf$Gene)}) %>% unlist() %>% unique()
virDf <- virDf[unique(c(testOut,ind)),]

virData <- virDf %>% left_join(ancientData, by = c("Gene","Name")) %>% mutate(KaeroEcoli = ifelse(is.na(KaeroEcoli),0,1))

##Getting them specifically now
#t3ss <- depthDfTrans %>% filter(grepl("^sep|^esc|^esp|^ces|^tir|^epr|^eiv|^epa|^etr|^mxi|^tss|^hcp|^vgr", Gene, ignore.case = T))
#upec <- depthDfTrans %>% filter(grepl("^hly|^sfa|^usp|^pap|^sat|^ipa|^yad|^foc|^hof|^fim", Gene, ignore.case = T)) #Morales-Espinosa et al 2016
#ehec <- depthDfTrans %>% filter(grepl("^iha|^eae|^stx|^nle|^bfp|^ehx|^esp|^etp|^kat|^ent", Gene, ignore.case = T)) #Bugarel et al 2011
#eaec <- depthDfTrans %>% filter(grepl("^agg|^astA|^tss|^aaf|^AAF|^pet", Gene, ignore.case = T)) #Croxen & Finlay 2010
#etec <- depthDfTrans %>% filter(grepl("^est|^elt|^lta|^cfa|^tib|^etp|^eat|^cyl|east", Gene, ignore.case = T))
#expec <- depthDfTrans %>% filter(grepl("^pap|^prf|^sfa|^gaf|^bma|^iha|^afa|^tsh|^ibea|^irea|^iuc|^ybt|^iro|^sita|^hlya|^cdt|^cnf|^hlyf|^clb|^sat|^pic|^trat|^ompt|^iss|^iss|^cva|^dsda|^malx", Gene, ignore.case =T))

#virDf <- virDf %>% bind_rows(list(t3ss, upec, ehec, eaec,etec,expec)) %>% distinct() 
#virDf <- virDf %>% filter(!grepl("hypothetical|conserved", Name))#, !grepl("elf|ecp|hcp|fim", Name))
#rm(ehec,etec,etecList,etecTranslations, expec, t3ss, testOut, upec, virGenes, ind, eaec)

#Now we'll be making the heatmap
#virData <- virDf %>% filter(Genome %in% strainListOlivier, MeanCoverage >= 1, PercentCoverage >= 0.9)  %>%
#	select(-c(sdCoverage, PercentCoverage, CV, Gene)) %>% pivot_wider(names_from = c(Genome), values_from=MeanCoverage) 
#
#ancientData <- virDf %>% filter(Genome %in% strainListFragmented, MeanCoverage >= 10, PercentCoverage >= 0.9) %>% 
#	select(-c(sdCoverage, PercentCoverage, CV, Gene)) %>% pivot_wider(names_from = c(Genome), values_from=MeanCoverage) 
#
FoundGenes <- ancientData %>% select(Name, KaeroEcoli) %>% filter(!(is.na(KaeroEcoli))) %>% pull(Name)
FoundGenes <- depthDfTrans %>% filter(Genome == "KaeroEcoli", Name %in% FoundGenes, !is.na(CV)) %>% select(-c(Genome))
write.table(FoundGenes, file ="~/SmallPanGenomeVirulenceList.tab", row.names = F, sep = "\t")

#virData <- virData %>% full_join(ancientData)
virData <- virData %>% mutate(Name = gsub(" \\[.*|MULTISPECIES: |, partial","",Name)) 
virData <- virData %>% select(-c(ID, Gene)) %>% group_by(Name) %>% summarize_all(sum) %>% group_by(Name) %>% mutate_all(as.logical)
virData[,-1] <- sapply(virData[,-1], function(x){ifelse(is.na(x),0,x)})
virData[,-1] <- sapply(virData[,-1], function(x){ifelse(x,1,0)})

virData <- as.data.frame(virData)
rownames(virData) <- virData[,1]
virData <- virData[,-1]

st4995 <- phylogroup[,1:2] %>% filter(ST == 4995) %>% pull(Genome)
st4995 <- st4995[st4995 %in% colnames(virData)]

tmp <- virData %>% select(c(st4995, "KaeroEcoli"))
tmp <- tmp[rowSums(tmp) > 0,]
rowSums(tmp) %>% table()



#st4995Genomes <- st4995Genomes[st4995Genomes %in% colnames(virData)]
#tmp <- virData[,c(st4995Genomes[-2], "KaeroEcoli")] %>% filter(KaeroEcoli == 1)
# Limiting the genes to only those which have a difference
maxCount <- virData %>% rowSums() %>% max()
ind <- virData %>% rowSums() == maxCount
virData <- virData[!ind,]
virData <- as.data.frame(t(virData))
rownames(virData)[which(rownames(virData) == "KaeroEcoli")] <- "AncientEcoli"

# Getting the Gene Virulence Clusters
virDist <- dist(t(virData))
clusteredVir <- hclust(virDist, method = "ward.D2")

sil <- sapply(2:15, function(x){
	       tmp <- cutree(clusteredVir, k = x)
   	       summary(silhouette(tmp, virDist))$avg.width
	 })

plot(y = sil, x = 2:15, type = "b")

small <- cutree(clusteredVir, k = 9)

clusterBorder <- colour[small + 8] %>% unique()
pdf("~/GeneVirulenceTree.pdf", width = 40, height = 12)
plot(clusteredVir, hang = -1, cex = 0.6)
rect.hclust(clusteredVir, k = 6, border = clusterBorder)
dev.off()

# Getting ready for the heatmap
specDist <- dist(virData)
clusteredSpec <- hclust(specDist, method = "ward.D2")
rownameOrder <- clusteredSpec$order
colnameOrder <- clusteredVir$order

ancientGenesFound <- as.logical(virData["AncientEcoli",])
ancientGenesFound <- data.frame(row.names = colnames(virData), "Ancient" = ifelse(ancientGenesFound, "Yes", "No"))

olivierTableHeatmap <- phylogroup[,c(4,2)] %>% as.data.frame()
olivierTableHeatmap$Pathovar <- sapply(olivierTableHeatmap$Pathovar, function(x){gsub("\\?", "Unknown",x)}) # 
rownamesHeatmap <- ifelse(rownames(virData) == "AncientEcoli", "AncientEcoli","")
rownames(olivierTableHeatmap) <- phylogroup$Genome
olivierTableHeatmap <- olivierTableHeatmap[rownames(virData),]

treeColours <- as.list(ann_colors[[1]])[olivierTableHeatmap[rownameOrder,]$Pathovar] %>% unlist()
treeColours <- gsub("#ffffff", "#000000", treeColours)
dend <- as.dendrogram(clusteredSpec)
labels_colors(dend) <- treeColours

pdf("~/Documents/University/EcoliPaperV2/AdditionalFiles/PartsofFigures/SpeciesVirulenceTreeA0.pdf", width = 6, height = 8)
plot(dend, horiz = T)
dev.off()

pheatmap(virData[rownameOrder,], annotation_row = olivierTableHeatmap, annotation_col = ancientGenesFound,
	 color = c("#007dba", "#f8333c"), clustering_method = "ward.D2", legend = F, show_colnames = F,
	 gaps_row = c(which(rownames(virData)[rownameOrder] == "AncientEcoli"),which(rownames(virData)[rownameOrder] == "AncientEcoli") - 1),cluster_rows = F,
	 labels_row = rownamesHeatmap[rownameOrder], fontsize_row = 5, 
	 filename = "~/Documents/University/EcoliPaperV2/AdditionalFiles/PartsofFigures/SmallVirHeatmapUpdateA0.png", width = 12, height = 8, border_color = NA, annotation_colors = ann_colors)#, main = "Core Presence/Absence")

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

RGIHomo %>% mutate(Drug.Class = gsub("; *", ";",Drug.Class)) %>% pull() %>%
	strsplit(";") %>% unlist() %>% table() %>% sort(decreasing = T) %>% length()

RGIHomo$Gene <- RGIHomo$Best_Hit_ARO %>% gsub(pattern = "Kleb.* |pneu.* |Escherichia |coli |beta-lactamase",replacement= "")

# Note, no coverage filtering on ancientGenesOnly because they all happen to be >= 10 in this case
#colnames(RGIHomo)[5] <- "Gene"

RGIHomo %>% filter(Gene %in% gsub("_.*","",ancientGenesOnly$Gene)) %>% inner_join(ancientGenesOnly) %>% select(Gene, MeanCoverage, sdCoverage, PercentCoverage, CV, Drug.Class, Resistance.Mechanism, Count) %>% arrange(-Count)
%>% write.csv("~/Documents/University/EcoliPaperV2/AdditionalFiles/AMRGenesFound.csv", row.names = F, quote = F)

tmp <- RGIHomo %>% filter(Gene %in%gsub("_.*","",ancientGenesOnly$Gene)) %>% mutate(Drug.Class = gsub("; *", ";",Drug.Class)) %>% pull(Drug.Class) %>%
	strsplit(";") %>% unlist() %>% tibble()
colnames(tmp) <- "Resistances"
withefflux <- tmp %>% count(Resistances, name = "With Efflux")

tmp <- RGIHomo %>% mutate(Resistance.Mechanism = gsub("antibiotic efflux","",Resistance.Mechanism),
		   Drug.Class = gsub("; *", ";",Drug.Class)) %>%
	filter(Gene %in% gsub("_.*","",ancientGenesOnly$Gene), Resistance.Mechanism != "") %>%
	pull(Drug.Class) %>%
	strsplit(";") %>% unlist() %>% tibble()
colnames(tmp) <- "Resistances"
noefflux <- tmp %>% count(Resistances, name = "No Efflux")
withefflux %>% full_join(noefflux) %>% arrange(-`With Efflux`, Resistances)%>% filter(!is.na(`No Efflux`)) %>%
       xtable:::xtable() %>% print(file = "~/Documents/University/EcoliPaperV2/AdditionalFiles/DrugClasses.tex")


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
