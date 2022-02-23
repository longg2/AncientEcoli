library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)

vcf <- "~/Scratch/EcoliNodule/ML30MQ30/PanGenomeV2/Heterozygosity/test.vcf"
geneLengths <- read.delim("~/Scratch/EcoliNodule/ML30MQ30/PanGenomeV2/Heterozygosity/panGenomeGeneLengths.tab", header =F, col.names = c("Gene", "Length")) %>% as_tibble()

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

vcf <- VCFParsing(vcf)


vcfHet <- vcf %>% full_join(geneLengths, by = c("CHROM" = "Gene")) %>% mutate(Hetero = ifelse(grepl("0/1|1/0", GT), T, F))
vcfHet %>% summarize(Hetero = sum(Hetero)) / vcfHet %>% distinct(CHROM, Length) %>% summarize(sum(Length)) %>% pull()
vcfHet %>% filter(QUAL >= 100) %>% summarize(sum(Hetero)) %>% pull()/ vcfHet %>% distinct(CHROM, Length) %>% summarize(sum(Length)) %>% pull()
vcfHet %>% filter(REF != "C" & ALT != "T", REF != "G" & ALT != "A")%>% summarize(sum(Hetero)) %>% pull()/ vcfHet %>% distinct(CHROM, Length) %>% summarize(sum(Length)) %>% pull()
vcfHet %>% filter(QUAL >= 100,REF != "C" & ALT != "T", REF != "G" & ALT != "A")%>% summarize(sum(Hetero)) %>% pull()/ vcfHet %>% distinct(CHROM, Length) %>% summarize(sum(Length)) %>% pull()

# Getting those Histograms
vcfHet %>% filter(!is.na(FILTER)) %>% group_by(CHROM, Length) %>% summarize(Hetero = sum(Hetero)) %>% summarize(HeteroFrac = log10(Hetero/Length)) %>%
	mutate(HeteroFrac = replace(HeteroFrac, is.infinite(HeteroFrac), 0)) %>% 
	ggplot(aes(x = HeteroFrac)) + geom_histogram(fill = "white",colour = "black") + geom_vline(aes(xintercept = mean(HeteroFrac)), colour = "red", lty = 2) + theme_bw()

vcfPlot <- vcfHet %>% filter(!is.na(FILTER),QUAL >= 100) %>%
       	group_by(CHROM, Length) %>% summarize(Hetero = sum(Hetero)) %>% summarize(HeteroFrac = log10(Hetero/Length)) %>%
	mutate(HeteroFrac = replace(HeteroFrac, is.infinite(HeteroFrac), 0))
vcfPlot %>% filter(HeteroFrac < 0) %>% ggplot(aes(x = HeteroFrac)) + geom_histogram(fill = "white",colour = "black") +
       	geom_vline(xintercept = vcfPlot %>% filter(HeteroFrac < 0) %>% summarize(mean(HeteroFrac)) %>% pull(), colour = "red", lty = 2) + theme_bw() +
	annotate(geom = "text", x = -2, y = 100, label = paste0("Genes with no Heterozygous Variants: ",sum(vcfPlot$HeteroFrac == 0))) +
	ylab("Genes") + xlab(bquote(log[10]("F(Heterozygous)")))


vcfHet %>% filter(REF != "C" & ALT != "T", REF != "G" & ALT != "A") %>% group_by(CHROM, Length) %>%
       	summarize(Hetero = sum(Hetero)) %>% summarize(HeteroFrac = log10(Hetero/Length)) %>% pull(HeteroFrac) %>% hist()

vcfHet %>% filter(QUAL >= 100,REF != "C" & ALT != "T", REF != "G" & ALT != "A") %>%
       	group_by(CHROM, Length) %>% summarize(Hetero = sum(Hetero)) %>% summarize(HeteroFrac = log10(Hetero/Length)) %>% pull(HeteroFrac) %>% hist()


