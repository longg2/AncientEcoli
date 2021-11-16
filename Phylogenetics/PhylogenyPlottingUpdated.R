library(dplyr)
library(ggplot2)
library(ggtree)
library(phangorn)
library(treeio)
library(ggpubr)
library(ggstance)
library(ggnewscale)
library(phytools)
library(ape)

colour <- c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6',
'#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3',
'#808000', '#ffd8b1', '#000075', '#808080', '#f0f0ff', '#000000')

ann_colors = list(Pathovar = c("Ancient" = colour[1], "Commensal" = "#000000","EAEC" = colour[2], "EIEC" = "#22662A", "ETEC" = colour[16], "ExPEC" = colour[4],
			       "Hybrid ExPEC/InPEC" = colour[3], "STEC" = colour[9],"Unknown" = colour[20]),
		  #Host = c("Environment" = colour[10], "Food" = colour[11], "Human" = colour[12], "Livestock" = colour[13], "Porc" = colour[14], "Wild animal" = colour[15],"?" = colour[3]),
		  ST = c("Ancient" = colour[1],"1429" = colour[6], "325" = colour[7], "3630" = colour[8], "399" = colour[12], "4995" = colour[18], "Other" = colour[20]),
		  Ancient = c("Yes" = colour[22], "No" = colour[21]))
#ann_colors = list(Pathovar = c("Ancient" = colour[1], "Commensal" = colour[2],"EAEC" = colour[4], "EIEC" = colour[5], "ETEC" = colour[6], "ExPEC" = colour[7],
#			       "Hybrid ExPEC/InPEC" = colour[8], "STEC" = colour[9],"Unknown" = colour[3]),
#		  Host = c("Environment" = colour[10], "Food" = colour[11], "Human" = colour[12], "Livestock" = colour[13], "Porc" = colour[14], "Wild animal" = colour[15],"Unknown" = colour[3]),
#		  ST = c("1429" = colour[6], "325" = colour[7], "3630" = colour[8], "399" = colour[12], "4995" = colour[18], "Other" = colour[20]),
#		  Ancient = c("Yes" = colour[22],"No" = colour[21]))

# Basic Tree
tree <- read.newick(file = "OlivierST4995/GubbinsCleanedSNPsOnly.fasta.treefile")
tree <- phytools::reroot(tree, interactive = T)
Tiplabels <- tree$tip.label

##################################################
# Let's colour these based on their PhyloGroups ##
##################################################
phylogroupFrame <- read.csv("../Virulence/VirulenceHeatmapOlivier.csv")[,1:4]
rownames(phylogroupFrame) <- phylogroupFrame$Genome
phylogroupFrame["AncientEcoli",] <- list(Genome = "AncientEcoli", ST = "Ancient", Host = "Ancient", Pathovar = "Ancient")
phylogroupFrame <- phylogroupFrame %>% mutate(Pathovar = sapply(Pathovar, function(x){gsub(pattern = "\\?",replacement = "Unknown", x)}))
#phylogroupFrame <- read.delim("OlivierST4995/PhyloMetadata.tab", stringsAsFactors = F)

#phylogroupFrame$Genome <- sapply(phylogroupFrame$Genome, function(x){ifelse(grepl("^ETEC",x), gsub("ETEC", "ETEC ",x),x)})
#colnames(phylogroupFrame)[1] <- "Genome"
#
#phylogroupFrame$Pathovar <- sapply(phylogroupFrame$Source, function(x){ifelse(grepl("Blood", x), "ExPEC", ifelse(grepl("Faeces",x), "Commensal",
#													 ifelse(grepl("Mummy",x), "Ancient", "Unknown")))})
#phylogroupFrame$Pathovar <- sapply(1:length(phylogroupFrame$ST), function(x){ifelse(grepl("ETEC", phylogroupFrame$ST[x]), "ETEC",phylogroupFrame$Pathovar[x])})
#phylogroupFrame$Source <- sapply(phylogroupFrame$Source, function(x){ifelse(grepl("animal|pig|,|Swine|Bovine|scrofa", x), "Animal",
#									    ifelse(grepl("Homo|ND|Blood|Faeces",x), "Human",
#										   ifelse(grepl("Cur|Cor|Ground", x), "Food", 
#											  ifelse(grepl("^$",x),"Unknown",x))))})
#
#phylogroupFrame$Genome  <- sapply(phylogroupFrame$Genome, function(x){ifelse(grepl("^ESC_", x), paste0(x,"_AS"),x)})
phylogroupFrame$Genome  <- gsub("^ETEC ", "ETEC",phylogroupFrame$Genome)
phylogroupFrame$Genome  <- gsub("_AS| ", "",phylogroupFrame$Genome)
phylogroupFrame <- phylogroupFrame[match(Tiplabels, phylogroupFrame$Genome),] # Getting them in the Right Order
phylogroupFrame[nrow(phylogroupFrame),] <- list(Genome = "IAI1", ST = "Other", Host = "Human", Pathovar = "Commensal")

cophenetic(tree)[c("AncientEcoli", "ESC_VA4573AA", "STEC388", "035-003-ccl", "IAI1"),c("AncientEcoli", "ESC_VA4573AA")]

# Now to read in the Core Gene statuses
coreGenecounts <- read.delim("../CoreGeneCounts.tab",)
coreGenecounts$Genome  <- gsub("^ETEC ", "ETEC",coreGenecounts$Genome)
coreGenecounts$Genome[which(coreGenecounts$Genome == "KaeroEcoli")] <- "AncientEcoli"
#virGenecounts <- read.delim("../Virulence/VirulenceGeneCounts.tab",)
#virGenecounts$Genome[nrow(coreGenecounts)] <- "AncientEcoli"
#colnames(virGenecounts) <- c("Genome", "Virulence")
colnames(coreGenecounts) <- c("Genome", "Core")

coreGenecounts$Genome <- gsub("_AS| ","", coreGenecounts$Genome)
#virGenecounts$Genome <- gsub("_AS| ","", virGenecounts$Genome)

FinalcoreGene <- phylogroupFrame %>% left_join(coreGenecounts) %>% mutate(Frac = Core/2555) %>% as.data.frame()# %>% left_join(virGenecounts)
rownames(FinalcoreGene) <- FinalcoreGene$Genome
FinalcoreGene <- FinalcoreGene[tree$tip.label,]
phylogroupFrame <- as.data.frame(phylogroupFrame)

# Let's get the labels with bars on it

#ancestor(tree, "AncientEcoli")

rownames(phylogroupFrame)  <- phylogroupFrame$Genome
#FinalcoreGene <- FinalcoreGene %>% mutate(Source = gsub(" H", "H", Source), Source = gsub("n ", "n", Source))

#FinalcoreGene$Source <- sapply(FinalcoreGene$Source, function(x){ifelse(grepl("animal|pig|;", x), "Animal",
#									    ifelse(grepl("Homo|Blood|Faeces|ND",x), "Human",
#										   ifelse(grepl("Cur", x), "Plant", ifelse(grepl("^$",x),"Unknown",x))))})
#phylogroupFrame$ST <- sapply(phylogroupFrame$ST, function(x){ifelse(grepl("A0", x), "A",)})
# Going to try replacing the labels with the Phylogroups
p1 <- ggtree(tree, right = T) %<+% FinalcoreGene +
       	geom_nodepoint(size = 1.5, colour = ifelse(as.numeric(tree$node.label) >= 90, "black", ifelse(as.numeric(tree$node.label) >= 50, "grey",NA)),
		     shape = "square") +
	geom_tippoint(mapping = aes(colour = ST), size = 1.5) +
	#geom_tiplab(align = T) + 
	geom_treescale(linesize = 1, offset = 2) + theme_tree() + geom_rootedge(rootedge = 0.005) +
	scale_color_manual(values = ann_colors$ST, name = "Sequence Type") +
	theme(legend.position = "bottom")# + guides(colour = guide_legend(nrow = 1))
#p2 <-  ggtree(tree, right = T) %<+% phylogroupFrame +
#	geom_tippoint(mapping = aes(colour = Source), size = 1.5) + geom_nodelab(size = 1.5) +
#	geom_tiplab(align = T) + 
#	geom_treescale(linesize = 1, offset = 2) + theme_tree() + geom_rootedge(rootedge = 0.01) +
#	scale_color_manual(values = colour, name = "Source") +
#	theme(legend.position = "bottom") + guides(colour = guide_legend(nrow = 1))
#p2
	
#ggsave("~/Workshop/VacationFiles/Labmeetings/May25/Figures/PhyloTree.pdf", width = 8, height = 6)

#p2 <- ggtree(tree, right = T) %<+% FinalcoreGene +
#       	geom_nodepoint(size = 1.5, colour = ifelse(as.numeric(tree$node.label) >= 90, "black", ifelse(as.numeric(tree$node.label) >= 50, "grey",NA)),
#		     shape = "square") +
#	geom_tippoint(mapping = aes(colour = Pathovar), size = 1.5) +
#	#geom_tiplab(align = T) + 
#	geom_treescale(linesize = 1, offset = 2) + theme_tree() + geom_rootedge(rootedge = 0.01) +
#	scale_color_manual(values = ann_colors$Pathovar, name = "Pathovar") +
#	theme(legend.position = "bottom") + guides(colour = guide_legend(nrow = 2))
#p3 <- ggtree(tree, right = T) %<+% FinalcoreGene +
#	theme_tree() + geom_rootedge(rootedge = 0.01) +
#       	geom_nodepoint(size = 1.5, colour = ifelse(as.numeric(tree$node.label) >= 90, "black", ifelse(as.numeric(tree$node.label) >= 50, "grey",NA)),
#		     shape = "square") +
#	geom_tippoint(mapping = aes(colour = Frac)) + 
#	scale_color_continuous(name = "Core Genome Presence",oob=scales::squish, low = "#f8333c", high = "#007dba") +
#	theme(legend.position = "bottom")
ggsave(p1, file = "~/Documents/University/EcoliPaperV2/AdditionalFiles/PartsofFigures/STPhylo.pdf", width = 12, height = 9)
#ggsave(p2, file = "~/Documents/University/EcoliPaperV2/AdditionalFiles/PartsofFigures/PathPhylo.pdf", width = 12, height = 9)
load("FullPhylo.RData")
load("ST4995Only.RData")
p23 <- ggarrange(p1,p3, ncol = 2, labels = c("B", "C"), align = "hv")

p123 <- ggarrange(p4, p23, nrow =2, labels = c("A",""), align = "hv", heights = c(4,3), common.legend = F, legend = "bottom")
ggsave(plot = p123, "~/Documents/University/EcoliPaperV2/Figures/Figure3ML.pdf", width = 12, height = 8)

#write.tree(tree, file = "MidpointTreemerTree.nwk")

########################
### The Tempest Plot ###
########################

# Now to get the Tempest Plot
dat <- as_tibble(read.delim("OlivierST4995/TempestResultsML.tab")) %>% filter(date > 0)
# Getting the Phylogroups plotted

phylogroupFrame <- phylogroupFrame[match(dat$tip, phylogroupFrame$Genome),] # Getting them in the Right Order
dat$ST <- phylogroupFrame$ST

tempest <- dat %>% ggplot(aes(x = date, y  = distance)) + geom_smooth(method = "lm", show.legend = F, colour = "black") +
	geom_point(alpha = 0.75, aes(colour = ST), show.legend = T) + theme_bw() + ylab("Root to Tip Divergence") + xlab("Year") +
	scale_color_manual(values = ann_colors$ST, name = "Sequence Type")

model <- lm(distance ~date, data = dat)
summary(model)
r2 <- round(summary(model)$adj.r.squared,3)
tmp <- summary(model)$fstatistic
pval <- round(pf(tmp[1],tmp[2],tmp[3], lower.tail = F),3)

tempest <- tempest + annotate(geom = "text", x = 1700, y = 0.03,label = bquote(R[adj]^2 == .(r2))) +
	annotate(geom = "text", x = 1700, y = 0.025, label = bquote(P == .(pval))) +
	scale_y_continuous(breaks = scales:::pretty_breaks(n = 10)) + theme(legend.position = "bottom") + coord_cartesian(ylim = c(0,0.05))
load("TempestAll.RData")
load("TempestST4995.RData")

ggarrange(tempestAll, tempest,tempestST4995ONLY, ncol = 1, align = "hv", labels = "AUTO")
ggsave("~/Documents/University/EcoliPaperV2/Figures/TempestPlots.pdf", width = 9, height = 12)

# Now to make the distance matrix
d <- cophenetic.phylo(tree)

subsetKeep <- c("AncientEcoli","ATCC11229", "STEC388","035-003-ccl", "IAI1")
d[subsetKeep,subsetKeep] %>% xtable:::xtable() %>% print(file = "~/Documents/University/EcoliPaperV2/AdditionalFiles/DistMatrix.tex")


