# Libraries
library(roll)
library(ggplot2)
library(ggpubr)
library(circlize)
library(seqinr)
library(dplyr)
colour = c("#f8333c", "#007dba")

PlasmidMapping <- function(depthtmp, snpstmp,windowKProp = 0.1, Chromo, fasta, gff){
	fasta <- read.fasta(fasta)
	#bed <- read.delim(bed, header = F)[,2:3]
	gff <- read.delim(gff, header = F, comment.char = "#", col.names = c("Acc", "DB", "Type", "Start", "Stop","Other", "Strand", "Other1", "Details")) %>%
		as_tibble() %>% filter(DB == "Protein Homology") %>% select(c(Acc, DB, Type, Start, Stop, Strand, Details)) %>%
		mutate(Pseudo = grepl("pseudo=true", Details))

	if(sum(gff$Stop > nrow(depthtmp))){ # If we have something going over then end
		over <- gff %>% slice(tail(row_number(), 1))
		lengthGene <- abs(over$Start - over$Stop)
		overbefore <- overafter <- over
		overbefore$Stop <- nrow(depthtmp)
		overafter$Start <- 1
		overafter$Stop <- lengthGene - (ifelse(overbefore$Stop - overbefore$Start > 0,overbefore$Stop - overbefore$Start,1))
		gff <- gff %>% slice(1:(n() - 1)) %>% bind_rows(overbefore) %>% bind_rows(overafter)
	}

	# Want to get the orientation of the arrows correct

	depthtmpV2 <- depthtmp %>% filter(Chromosome == Chromo) %>% mutate(Base = as.character(fasta[[Chromo]]))
	snpstmpV2 <- snpstmp %>% filter(CHROM == Chromo) %>% select(CHROM, POS) %>% mutate(SNP = 1)
	colnames(snpstmpV2) <- c("Chromosome", "Position", "SNP")
	majorTicks <- nrow(depthtmpV2)/10
	depthtmpV3 <- depthtmpV2 %>% full_join(snpstmpV2) %>% mutate(windowK = floor(windowKProp * length(Chromosome)), Interval = floor(`Position`/windowK), SNP = ifelse(is.na(SNP), 0, 1)) %>% 
		group_by(Chromosome,Interval) %>%
		summarize(MeanCoverage = mean(Coverage), GCContent = sum(ifelse(grepl("g|c|G|C",Base), T, F))/length(Base), windowK = windowK, SNP = sum(SNP)) %>% distinct() %>%
		mutate(`Position` = Interval * windowK)
	#depthtmpV2 <- depthtmpV2 %>% mutate(Interval = floor(`Position`/windowK)) %>% group_by(Chromosome, Interval) %>%
	#	summarize(MeanCoverage = mean(Coverage), GCContent = sum(ifelse(grepl("g|c|G|C",Base), T, F))/length(Base)) %>%
	#	mutate(`Position` = Interval * windowK)

	# Getting GC ConfInt
	tmp <- depthtmpV2 %>% summarize(meanCoverage = mean(Coverage), sdCoverage = sd(Coverage), Error = qnorm(0.975) * sdCoverage/sqrt(length(Coverage)), low = meanCoverage - Error, hi = meanCoverage + Error)
	confIntCov <- c(tmp$meanCoverage, tmp$low, tmp$hi)
	#tmp <- depthtmpV2 %>% summarize(meanGC = mean(GCContent), sdGC = sd(GCContent), Error = qnorm(0.975) * sdGC/sqrt(length(GC)), low = meanGC - Error, hi = meanGC + Error)
	#confIntGC <- c(tmp$low, tmp$hi)
	
	#depthtmpV3 <- depthtmpV3 %>% mutate(SNP = ifelse(SNP != 0, 18, NA))
	# Running
	circos.par(cell.padding = c(0.00, 0, 0.00, 0), gap.after = 15, start.degree = -278)
	circos.initialize(sectors = depthtmpV3$Chromosome,x = depthtmpV3$Position)
	 # Making the Track
	circos.track(depthtmpV3$Chromosome, y = depthtmpV3$MeanCoverage,
		panel.fun = function(x,y){
			circos.text(max(depthtmpV3$Position)*1.02,
	            CELL_META$cell.ylim[2] + mm_y(4),
	            CELL_META$sector.index)
	        circos.genomicAxis(h = "top", major.by = majorTicks)
	})

	####### The Genome Coverage ######
	# Colouring in Zones
	circos.rect(xleft = gff$Start, xright = gff$Stop, ybottom =0, ytop =max(depthtmpV3$MeanCoverage),
		border = NA, col = ifelse(gff$Pseudo, "#22c3b0", "#1ab2ff"))
		#border = NA, col = ifelse(gff$Pseudo, "#22c3b0", "#007dba"))
#	circos.rect(xleft = 0, xright = max(depthtmpV3$Position), ybottom =0, ytop = 10,
#		 border = NULL, col = c("red"))
	# Axis
	circos.yaxis(side = "left", labels.cex = 0.5)
	# Lines
	circos.trackLines(sectors = depthtmpV3$Chromosome, x = depthtmpV3$Position, y = depthtmpV3$MeanCoverage, area = T, col = "#808080")
	circos.segments(x0 = 0, y0 = confIntCov[1],
			x1 = max(depthtmpV3$Position), y1 =confIntCov[1], col = "#f8333c", lty = 2)

	####### Gene Arrow ######
	circos.genomicTrack(gff[,c(1,4,5,8)], stack = T,
			    panel.fun = function(region, value,...){
				    circos.genomicLines(region, value, type = "segment", lwd = 2,
							col = ifelse(gff$Pseudo, "#22c3b0", "#1ab2ff"), y = getI(...)*0.25)
		 },
		 bg.border = NA, track.height = 0.01)
	####### SNP ########
	circos.track(depthtmpV3$Chromosome, y = depthtmpV3$SNP, ylim = c(0,3))
	circos.yaxis(side = "left", labels.cex = 0.5, at = 1:3)
	circos.trackLines(sectors = depthtmpV3$Chromosome, x = depthtmpV3$Position,,y = depthtmpV3$SNP, col = "#8279b9", area = T)

	####### The GC Content ######
	circos.track(depthtmpV3$Chromosome, y = depthtmpV3$GCContent, ylim = c(0,1))
	circos.yaxis(side = "left", labels.cex = 0.5)
	#circos.rect(xleft = 0, xright = max(depthtmpV3$Position), ybottom =confIntGC[1] , ytop = confIntGC[2],
	#	 border = NA, col = c("grey90")) # Making the confidence interval
	circos.trackLines(sectors = depthtmpV3$Chromosome, x = depthtmpV3$Position, y = depthtmpV3$GCContent, col = "orange")

	circos.clear()
}

GenomeMapping <- function(depthtmp, snpstmp,windowKProp = 0.1, Chromo, fasta){
	fasta <- read.fasta(fasta)
	depthtmpV2 <- depthtmp %>% filter(Chromosome == Chromo) %>% mutate(Base = as.character(fasta[[Chromo]]))
	snpstmpV2 <- snpstmp %>% filter(CHROM == Chromo) %>% select(CHROM, POS) %>% mutate(SNP = 1)
	colnames(snpstmpV2) <- c("Chromosome", "Position", "SNP")
	majorTicks <- nrow(depthtmpV2)/10
	depthtmpV3 <- depthtmpV2 %>% full_join(snpstmpV2) %>% mutate(windowK = floor(windowKProp * length(Chromosome)), Interval = floor(`Position`/windowK), SNP = ifelse(is.na(SNP), 0, 1)) %>% 
		group_by(Chromosome,Interval) %>%
		summarize(MeanCoverage = mean(Coverage), GCContent = sum(ifelse(grepl("g|c|G|C",Base), T, F))/length(Base), windowK = windowK, SNP = sum(SNP)) %>% distinct() %>%
		mutate(`Position` = Interval * windowK)
	#depthtmpV2 <- depthtmpV2 %>% mutate(Interval = floor(`Position`/windowK)) %>% group_by(Chromosome, Interval) %>%
	#	summarize(MeanCoverage = mean(Coverage), GCContent = sum(ifelse(grepl("g|c|G|C",Base), T, F))/length(Base)) %>%
	#	mutate(`Position` = Interval * windowK)

	# Getting the mean Coverage
	tmp <- depthtmpV2 %>% summarize(meanCoverage = mean(Coverage), sdCoverage = sd(Coverage), Error = qnorm(0.975) * sdCoverage/sqrt(length(Coverage)), low = meanCoverage - Error, hi = meanCoverage + Error)
	confIntCov <- c(tmp$meanCoverage, tmp$low, tmp$hi)
	## Getting GC ConfInt
	#tmp <- depthtmpV2 %>% summarize(meanGC = mean(GCContent), sdGC = sd(GCContent), Error = qnorm(0.975) * sdGC/sqrt(length(GC)), low = meanGC - Error, hi = meanGC + Error)
	#confIntGC <- c(tmp$low, tmp$hi)

	# Running
	circos.par(cell.padding = c(0.00, 0, 0.00, 0), gap.after = 15, start.degree = -278)
	circos.initialize(sectors = depthtmpV3$Chromosome,x = depthtmpV3$Position)
	 # Making the Track
	circos.track(depthtmpV3$Chromosome, y = depthtmpV3$MeanCoverage,
		panel.fun = function(x,y){
			circos.text(max(depthtmpV3$Position)*1.02,
	            CELL_META$cell.ylim[2] + mm_y(4),
	            CELL_META$sector.index)
	        circos.genomicAxis(h = "top", major.by = majorTicks)
	})

	####### The Genome Coverage ######
	# Colouring in Zones
#	circos.rect(xleft = 0, xright = max(depthtmpV3$Position), ybottom =0, ytop = 10,
#		 border = NULL, col = "red")
	# Axis
	circos.yaxis(side = "left", labels.cex = 0.5)
	# Lines
	circos.trackLines(sectors = depthtmpV3$Chromosome, x = depthtmpV3$Position, y = depthtmpV3$MeanCoverage, area = T, col = "#808080")
	circos.segments(x0 = 0, y0 = confIntCov[1],
			x1 = max(depthtmpV3$Position), y1 =confIntCov[1], col = "#f8333c", lty = 2)
	####### SNP ########
	circos.track(depthtmpV3$Chromosome, y = depthtmpV3$SNP)
	circos.yaxis(side = "left", labels.cex = 0.5)
	circos.trackLines(sectors = depthtmpV3$Chromosome, x = depthtmpV3$Position,,y = depthtmpV3$SNP, col = "#8279b9", area = T)

	####### The GC Content ######
	circos.track(depthtmpV3$Chromosome, y = depthtmpV3$GCContent, ylim = c(0,1))
	circos.yaxis(side = "left", labels.cex = 0.5)
#	circos.rect(xleft = 0, xright = max(depthtmpV3$Position), ybottom =confIntGC[1] , ytop = confIntGC[2],
#		 border = NA, col = c("grey90")) # Making the confidence interval
	circos.trackLines(sectors = depthtmpV3$Chromosome, x = depthtmpV3$Position, y = depthtmpV3$GCContent, col = "orange")
	circos.clear()
}

GenomeMappingAll <- function(depthtmp, snpstmp, windowKProp = 0.1, fasta){
	fasta <- read.fasta(fasta)
	depthtmpV2 <- split(depthtmp,f = depthtmp$Chromosome) %>%
	       	lapply(function(x){x %>% mutate(Base = as.character(fasta[[x$Chromosome[1]]]))}) %>% bind_rows()
	snpstmpV2 <- snpstmp %>% select(CHROM, POS) %>% mutate(SNP = 1)
	colnames(snpstmpV2) <- c("Chromosome", "Position", "SNP")
	#majorTicks <- nrow(depthtmpV2)/10
	depthtmpV3 <- depthtmpV2 %>% full_join(snpstmpV2) %>% group_by(Chromosome) %>%
	       	mutate(windowK = floor(windowKProp * length(Chromosome)), Interval = floor(`Position`/windowK), SNP = ifelse(is.na(SNP), 0, 1)) %>% 
		group_by(Chromosome,Interval) %>%
		summarize(MeanCoverage = mean(Coverage), GCContent = sum(ifelse(grepl("g|c|G|C",Base), T, F))/length(Base), windowK = windowK, SNP = sum(SNP)) %>% distinct() %>%
		mutate(`Position` = Interval * windowK)

	# Getting the mean Coverage
	tmp <- depthtmpV2 %>% group_by(Chromosome)%>% summarize(meanCoverage = mean(Coverage), sdCoverage = sd(Coverage), Error = qnorm(0.975) * sdCoverage/sqrt(length(Coverage)), low = meanCoverage - Error, hi = meanCoverage + Error)
	confIntCov <- c(tmp$meanCoverage, tmp$low, tmp$hi)
	depthtmpV3 <- depthtmpV3 %>% left_join(tmp) %>% ungroup() %>% distinct()
	## Getting GC ConfInt
	#tmp <- depthtmpV2 %>% summarize(meanGC = mean(GCContent), sdGC = sd(GCContent), Error = qnorm(0.975) * sdGC/sqrt(length(GC)), low = meanGC - Error, hi = meanGC + Error)
	#confIntGC <- c(tmp$low, tmp$hi)

	# Running
	#depthtmpV2$Chromosome <- as.factor(depthtmpV2$Chromosome)
	tmp2<- split(depthtmpV3, depthtmpV3$Chromosome)
	longerFragments <- sort(sapply(tmp2, function(x){max(x$Position)}), decreasing = T)
	tmp2 <- tmp2[names(longerFragments)[1:10]] %>% bind_rows()
	tmp2$Chromosome <- factor(tmp2$Chromosome, levels = names(longerFragments)[1:10])

	circos.par(cell.padding = c(0.02, 0, 0.02, 0), gap.degree = 5, start.degree = -292)
	#circos.par(cell.padding = c(0.02, 0, 0.02, 0))
	circos.initialize(sectors = tmp2$Chromosome,x = tmp2$Position)
	 # Making the Track
	circos.track(tmp2$Chromosome, y = tmp2$MeanCoverage,
		panel.fun = function(x,y){
	        circos.genomicAxis(h = "top")
	})

	####### The Genome Coverage ######
	# Colouring in Zones
#	circos.rect(xleft = 0, xright = max(tmp2$Position), ybottom =0, ytop = 10,
#		 border = NULL, col = "red")
	# Axis
	circos.yaxis(side = "left", labels.cex = 0.5)
	# Lines
	circos.trackLines(sectors = tmp2$Chromosome, x = tmp2$Position, y = tmp2$MeanCoverage, area = T, col = "#808080")
	circos.trackLines(sectors = tmp2$Chromosome, x = tmp2$Position, y = tmp2$meanCoverage, col = "#f8333c", lty = 2)

	####### SNP ########
	circos.track(tmp2$Chromosome, y = tmp2$SNP)
	circos.yaxis(side = "left", labels.cex = 0.5)
	circos.trackLines(sectors = tmp2$Chromosome, x = tmp2$Position,y = tmp2$SNP, col = "#8279b9", area = T)
	####### The GC Content ######
	circos.track(tmp2$Chromosome, y = tmp2$GCContent, ylim = c(0,1))
	circos.yaxis(side = "left", labels.cex = 0.25)
#	circos.rect(xleft = 0, xright = max(tmp2$Position), ybottom =confIntGC[1] , ytop = confIntGC[2],
#		 border = NA, col = c("grey90")) # Making the confidence interval
	circos.trackLines(sectors = tmp2$Chromosome, x = tmp2$Position, y = tmp2$GCContent, col = "orange")
	circos.clear()

	# Now for the ggplot version -- Because there are too many!!!
	depthtmpV4 <- depthtmpV3 %>% group_by(Chromosome) %>% mutate(Chromosome = gsub(".*(?=NODE)|_cov.*", "", Chromosome, perl = T))
	ordered <- depthtmpV4 %>% pull(Chromosome) %>% unique() %>% gsub(pattern = ".*_",replacement = "") %>% as.numeric() %>% order(decreasing = T)
	tmp <- depthtmpV4 %>% pull(Chromosome) %>% unique()
	depthtmpV4$Chromosome <- factor(depthtmpV4$Chromosome, level = tmp[ordered])
	depthtmpV4 <- depthtmpV4 %>% arrange(Chromosome)

	# Splitting the Genome into thirds
	ScaffLengths <- depthtmpV4 %>% pull(Chromosome) %>% unique() %>% gsub(pattern = ".*_",replacement = "") %>% as.numeric()
	AdditionalLengths <- sapply(1:length(ScaffLengths), function(x){
		       if(x == 1){
			       return(0)
		       }else{
			       return(sum(ScaffLengths[1:(x-1)])+1)
		       }
	})
	ScaffLengthsDf <- tibble("Chromosome" = depthtmpV4 %>% pull(Chromosome) %>% unique(),"Length" = ScaffLengths) %>% arrange(-Length) %>% 
		mutate(AdditionalLength = AdditionalLengths) 
	depthtmpV4 <- depthtmpV4 %>% left_join(ScaffLengthsDf) %>% mutate(Position = AdditionalLength + Position)

	thirds <- floor(max(depthtmpV4$Position)/3)
	depthtmpV4 <- depthtmpV4 %>% mutate(Thirds = ifelse(Position <= thirds, "First", ifelse(Position > thirds & Position <= 2*thirds, "Second", "Third"))) %>%
		mutate(Boundary = ifelse(Interval == 0, AdditionalLength,NA))

	ggCov <- depthtmpV4 %>% mutate(Position = ifelse(Thirds == "First", Position, ifelse(Thirds == "Second", Position - thirds, Position - 2*thirds))) %>%
		       mutate(Boundary = ifelse(Thirds == "First", Boundary, ifelse(Thirds == "Second", Boundary - thirds, Boundary - 2*thirds))) %>%
		ggplot(aes(x = Position, y = MeanCoverage)) + geom_vline(aes(xintercept = Boundary), lty = 2, col = "red") + geom_line() + theme_bw() +
		facet_grid(Thirds~.) + theme(axis.title.x = element_blank(), axis.text.x = element_blank())
	ggSNP <- depthtmpV4 %>% mutate(Position = ifelse(Thirds == "First", Position, ifelse(Thirds == "Second", Position - thirds, Position - 2*thirds))) %>%
		       mutate(Boundary = ifelse(Thirds == "First", Boundary, ifelse(Thirds == "Second", Boundary - thirds, Boundary - 2*thirds))) %>%
		ggplot(aes(x = Position, y = SNP)) + geom_vline(aes(xintercept = Boundary), lty = 2, col = "red") + geom_line() + theme_bw() +
		facet_grid(Thirds~.)
		

#	tmp <- depthtmpV4 %>% select(Chromosome, windowK) %>% distinct() %>% arrange(-windowK) %>% pull(Chromosome)
#	depthtmpV4$Chromosome <- factor(depthtmpV4$Chromosome, level = tmp)
#	ggCov <- depthtmpV4 %>% ggplot(aes(x = Position, y = MeanCoverage)) + geom_line() + theme_bw() + facet_wrap("Chromosome", scales = "free_x") +
#		geom_line(aes(y = meanCoverage), col = "#f8333c", lty = 2)
#	ggSNP <- depthtmpV4 %>% ggplot(aes(x = Position, y = SNP)) + geom_line() + theme_bw() + facet_wrap("Chromosome", scales = "free_x")
	return(list(ggCov, ggSNP))
}

T6SSKaeroPlot <- function(depthtmp, windowK = 1000, Chromo, fasta, gff){

	# Doing the basics here
	fasta <- read.fasta(fasta)
	depthtmpV2 <- depthtmp %>% filter(Chromosome == Chromo) %>% mutate(Base = as.character(fasta[[Chromo]]))
	depthtmpV2 <- depthtmpV2 %>% mutate(Interval = floor(`Position`/windowK)) %>%
		filter(`Position` >= 4130000, `Position` <= 4170000) %>%
		group_by(Chromosome, Interval) %>%
		summarize(MeanCoverage = mean(Coverage), GCContent = sum(ifelse(grepl("g|c|G|C",Base), T, F))/length(Base)) %>%
		mutate(`Position` = Interval * windowK)

	# Reading the GFF file
	gff <- read.delim(gff, header = F, comment.char = "#", col.names = c("Acc", "DB", "Type", "Start", "Stop","Other", "Strand", "Other1", "Details")) %>%
		as_tibble() %>% filter(DB == "Protein Homology", Start >= 4130000, Stop <= 4170000) %>% select(c(Acc, DB, Type, Start, Stop, Strand, Details))

	tmpgff <- gff$Details %>% strsplit(split = ";")
	
	gff$Names <- sapply(tmpgff, function(x){
		       y <- grep("gene=",x)
		       if(length(y) == 0){
				return(NA)
		       }
		       z <- x[y]
		       return(gsub(".*=","",z))
	})
	gff$Fill <- ifelse(is.na(gff$Names), "Unnamed", "Named")

	gff$LabelPos <- ifelse(1:nrow(gff) %% 2, 12.5,10)

	# Getting GC ConfInt
	tmp <- depthtmpV2 %>% summarize(meanCoverage = mean(MeanCoverage), sdCoverage = sd(MeanCoverage), Error = qnorm(0.975) * sdCoverage/sqrt(length(MeanCoverage)), low = meanCoverage - Error, hi = meanCoverage + Error)
	confIntCov <- c(tmp$meanCoverage, tmp$low, tmp$hi)

	depthtmpV2 %>% ggplot(aes(x = `Position`, y = MeanCoverage)) + theme_bw() +
		geom_rect(data = gff, mapping = aes(xmin = Start, xmax = Stop, ymin = 10, ymax = 30, fill = Fill),
			  inherit.aes = F, alpha = 0.75, colour = "white") + 
		scale_fill_manual(values = c("Named" = "#1ab2ff", "Unnamed" = "#22c3b0")) + 
		geom_line() +
		geom_hline(yintercept = tmp$meanCoverage, col = "#f8333c", lty = 2) +
		geom_label(data = gff, mapping = aes(x = (Start + Stop )/ 2, y = LabelPos, label = Names), fontface = "italic",inherit.aes = F, size = 2) +
		labs(fill = "Gene") +# ggtitle(subtitle = "Type VI Secretion System ", label = "Klebsiella aerogenes") +
		theme(legend.position = "bottom", plot.title = element_text(face = "italic"))

}


# Figure 2 of LID
pdf("~/Figure2.pdf", width = 12, height = 12)
par(mfrow = c(2,2))

## O104
#depthOrig <- as_tibble(read.table("Genomes/O104Depth.tab.gz", header = F))
#colnames(depthOrig) <- c("Chromosome", "Position", "Coverage")
#GenomeMapping(depthOrig, Chromo = "NC_018658.1", fasta = "FastaFiles/O104.fasta",windowKProp = 0.001)

# K12
depthOrig <- as_tibble(read.table("Genomes/K12Depth.tab.gz", header = F))
snpsOrig <- as_tibble(read.table("SNPS/K12Snps.vcf", comment.char = "#", header = T))
colnames(depthOrig) <- c("Chromosome", "Position", "Coverage")
GenomeMapping(depthOrig, snpsOrig,Chromo = "NC_000913.3", fasta = "FastaFiles/K12.fasta",windowKProp = 0.001)

depthOrig %>% summarize(Mean = mean(Coverage), SD = sd(Coverage), ConfInt = qnorm(0.975)*SD/sqrt(length(Coverage)), high = Mean + ConfInt, low = Mean - ConfInt) %>% as.data.frame()

# ESC_VA4573AA
depthOrig <- as_tibble(read.table("Genomes/ESC_VA4Depth.tab.gz", header =F))
snpsOrig <- as_tibble(read.table("SNPS/ESC_VA4Snps.vcf", comment.char = "#", header = T))
colnames(depthOrig) <- c("Chromosome", "Position", "Coverage")
pdf(file = "~/Documents/University/EcoliPaperV2/AdditionalFiles/ESCMapping.pdf", width = 6, height = 6)
allPlot <- GenomeMappingAll(depthOrig,snpsOrig, windowKProp = 0.01, fasta = "FastaFiles/ESC_VA4573AA_AS.fasta")
dev.off()

depthOrig %>% group_by(Chromosome) %>% summarize(Mean = mean(Coverage), SD = sd(Coverage),PercentCoverage1 = sum(Coverage >= 1)/length(Coverage)) %>% filter(Mean < 1)

tmp2 <- depthOrig %>% group_by(Chromosome) %>% summarize(Mean = mean(Coverage), SD = sd(Coverage),PercentCoverage1 = sum(Coverage >= 1)/length(Coverage)) %>% filter(Mean >= 1)

tmp2 %>% summarize(meanFinal = mean(Mean), sdFinal = sd(Mean)) %>% as.data.frame()

tmp2 %>% ungroup() %>% summarize(MeanOver = mean(Mean), SD = sd(Mean), ConfInt = qnorm(0.975)*SD/sqrt(length(Mean)), high = MeanOver + ConfInt, low = MeanOver - ConfInt) %>% as.data.frame()

# Klebsiella genome
pdf("~/Documents/University/EcoliPaperV2/Figures/KlebGenome.pdf", width = 6, height = 6)
depthOrig <- as_tibble(read.table("Genomes/KaeroFinalDepths.tab.gz", header = F))
colnames(depthOrig) <- c("Chromosome", "Position", "Coverage")
GenomeMappingKaero(depthOrig, Chromo = "NC_015663.1", windowK = 1e3, fasta = "FastaFiles/KaeroRef.fna")
dev.off()

#### The two Plasmids
png("~/Documents/University/EcoliPaperV2/Figures/PlasmidGenomes.png", width = 12, height = 9, res = 300, units = "in")
par(mfrow = c(1,2))
depthOrig <- as_tibble(read.table("Genomes/CP019906Depth.tab.gz", header = F))
colnames(depthOrig) <- c("Chromosome", "Position", "Coverage")
PlasmidMapping(depthOrig, Chromo = "CP019906.1", fasta = "FastaFiles/CP019906.fasta", windowK = 30, gff = "GeneFeatures/CP019906.gff3")

depthOrig %>% summarize(Mean = mean(Coverage), SD = sd(Coverage), ConfInt = qnorm(0.975)*SD/sqrt(length(Coverage)), high = Mean + ConfInt, low = Mean - ConfInt) %>% as.data.frame()

depthOrig <- as_tibble(read.table("Genomes/CP012732Depth.tab.gz", header = F))
colnames(depthOrig) <- c("Chromosome", "Position", "Coverage")
PlasmidMapping(depthOrig, Chromo = "CP012732.1", fasta = "FastaFiles/CP012732.fasta", windowK = 30, gff = "GeneFeatures/CP012732.gff3")
legend(x = 0.25, y = -1.2, legend = c("Protein-Coding", "Pseudo"), fill = c("#007dba","#22c3b0"))

depthOrig %>% summarize(Mean = mean(Coverage), SD = sd(Coverage), ConfInt = qnorm(0.975)*SD/sqrt(length(Coverage)), high = Mean + ConfInt, low = Mean - ConfInt) %>% as.data.frame()
dev.off()

# T6SS Testing
depthOrig <- as_tibble(read.table("Genomes/KaeroFinalDepths.tab.gz", header = F, col.names = c("Chromosome", "Position", "Coverage")))
depthOrig %>% group_by(Chromosome) %>% filter(`Position` >= 4130000, `Position` <= 4170000) %>% summarize(Mean = mean(Coverage), SD = sd(Coverage), ConfInt = qnorm(0.975)*SD/sqrt(length(Coverage)), high = Mean + ConfInt, low = Mean - ConfInt, PercentCoverage1 = sum(Coverage >= 1)/length(Coverage)) %>% as.data.frame()

fasta = "FastaFiles/KaeroRef.fna"
gff = "GeneFeatures/Kaero.gff3"
Chromo = "NC_015663.1"

t6ssPlot <- T6SSKaeroPlot(windowK = 100, gff = gff, fasta = fasta, Chromo = Chromo, depthtmp = depthOrig)
load("~/GeneHistogram.RData")
ggarrange(p2,plot.new(),plot.new(), t6ssPlot, nrow = 2, ncol = 2, labels = "AUTO")
ggsave("~/Figure2V2.pdf", width = 12, height = 8)
