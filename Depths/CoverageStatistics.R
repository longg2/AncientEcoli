# Libraries
library(seqinr)
library(dplyr)
library(emmeans)

# Let's get the Plasmids here
depthOrig <- as_tibble(read.table("Genomes/CompPlasmids.tab.gz", header = F))
colnames(depthOrig) <- c("Chromosome", "Position", "Coverage")
depthCP012732 <- depthOrig %>% filter(Chromosome == "CP012732.1")
depthCP019906 <- depthOrig %>% filter(Chromosome == "CP019906.1")

# Now the T6SS
depthOrig <- as_tibble(read.table("Genomes/KaeroFinalDepths.tab.gz", header = F, col.names = c("Chromosome", "Position", "Coverage")))
deptht6ss <- depthOrig %>% group_by(Chromosome) %>% filter(`Position` >= 4130000, `Position` <= 4170000)# %>% summarize(Mean = mean(Coverage), SD = sd(Coverage), ConfInt = qnorm(0.975)*SD/sqrt(length(Coverage)), high = Mean + ConfInt, low = Mean - ConfInt, PercentCoverage10 = sum(Coverage >= 10)/length(Coverage)) %>% as.data.frame()

# Now to get the P/A data in
depthDfNotAvg <- lapply("../Virulence/WholePangenome/WholePanGeneome.tab.gz", function(f){
		tmp <- as_tibble(read.delim(f, header = F, col.names = c("Gene", "Pos", "Coverage")))
		tmp$Genome <- gsub(".*/", "", gsub(".tab.gz","",f)) 	
		return(tmp)
	}) %>% bind_rows()

identifiedPanGenes <- read.delim("../Virulence/WholePangenome/IdentifiedGenes.tab")
depthDfNotAvg <- depthDfNotAvg %>% filter(Gene %in% identifiedPanGenes$Gene)

# Tukey for Depths
depthDfNotAvgTuk <- depthDfNotAvg %>% rename(Position = Pos, Chromosome = Genome) # For compatibility
compData <- bind_rows(depthCP019906, depthCP012732) %>% bind_rows(depthDfNotAvgTuk[,-1]) %>% bind_rows(deptht6ss)
tmp <- lm(data = compData, Coverage ~ Chromosome) %>% emmeans("Chromosome") %>% pairs()
tmp
summary(tmp)$p.value
plot(tmp)
