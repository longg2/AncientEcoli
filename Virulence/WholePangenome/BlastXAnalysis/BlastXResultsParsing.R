library(dplyr)
library(parallel)
library(pbapply)

# Reading the results
tmp <- read.delim("groupPanGenome.tab.gz", header = F,
		  col.names = c("Query", "Match", "PIdent", "Length", "Mismatches", "Gapopen", "QStart", "QEnd", "SStart", "SSend", "Evalue", "Bitscore", "Taxa")) %>% as_tibble()

tmpList <- split(tmp, as.factor((tmp$Query)))

op <- pboptions(type = "timer")
ncores <- 10

FilteredAmbiguousGenes<- pblapply(cl = 10,tmpList, function(x){ 
	       	filtResults <- x %>% filter(Evalue == min(Evalue))

		if(nrow(filtResults) > 1){
			filtResults <- filtResults %>% filter(Bitscore == max(Bitscore))
			if(nrow(filtResults) > 1){
				filtResults <- filtResults %>% filter(PIdent == max(PIdent))
				}
				if(nrow(filtResults) > 1){
					filtResults <- filtResults %>% filter(Mismatches == min(Mismatches))
					}
			}
		return(filtResults[1,])
 }) %>% bind_rows() %>% filter(PIdent >= 90)

FilteredAmbiguousGenes %>% pull(Match) %>% unique() %>% write.table("AmbiguousHits.list", col.names = F, quote = F, row.names = F) # From here we run epost/esummary/xtract to get the titles

# Assuming its been done
classifiedHits <- read.delim("IdentifiedGenes.list", header = F, col.names = c("Match", "Name", "Organism"))
FilteredAmbiguousGenes <- FilteredAmbiguousGenes %>% mutate(Match = gsub("\\.[0-9]$","", Match))
tmp <- FilteredAmbiguousGenes %>% left_join(classifiedHits) %>% filter(grepl(" VI |T6SS|type VI", Name))
FilteredAmbiguousGenes %>% left_join(classifiedHits) %>% select(Query, Match, Name, Organism) %>% 
	write.table(file = "ClassifiedGenes.tab",row.name = F, quote = F, sep = "\t")

tmp %>% pull(Query) %>% write.table(file = "T6SSGenesFound.list", col.name = F, row.name = F, quote = F)

