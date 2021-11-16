library(dplyr)
library(ggplot2)
library(ggtree)
library(phangorn)
library(treeio)
library(ggpubr)
library(ggstance)
library(ggnewscale)

colour <- c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6',
'#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3',
'#808000', '#ffd8b1', '#000075', '#808080', '#ffffff', '#000000')

# Basic Tree
tree <- read.beast(file = "BeastLogRelaxedTree.nexus")
#tree <- read.beast(file = "BEAST/LogNormalDate/NoScaffoldsRecombFix.nwk")
tree@phylo <- midpoint(tree@phylo)
tree@phylo <- reorder(tree@phylo)

Tiplabels <- tree@phylo$tip.label

##################################################
# Let's colour these based on their PhyloGroups ##
##################################################
phylogroupFrame <- read.delim("~/Downloads/AllFinal_phylogroups.txt", header = F, stringsAsFactors = F)[,c(1,5)]
colnames(phylogroupFrame) <- c("Genome", "Phylogroup")
phylogroupFrame[456,] <- list("O104","B1")

Tiplabels <- gsub(".fna|.fasta|.fasta.ref","",Tiplabels) # Removing the last bit
phylogroupFrame$Genome <- gsub(".fna|.fasta|.fasta.ref","",phylogroupFrame$Genome) # Removing the last bit

phylogroupFrame <- phylogroupFrame[match(Tiplabels, phylogroupFrame$Genome),] # Getting them in the Right Order
notFound <- which(sapply(phylogroupFrame$Genome, is.na)) # Which were not present in the Phylogrouping?
phylogroupFrame$Genome[notFound] <- Tiplabels[notFound]
#phylogroupFrame$Phylogroup[notFound] <- c("Efergusonii", "Ancient")
phylogroupFrame$Phylogroup[notFound] <- c("Unknown")

# Removing Phylogroup Info from Shigella to make it stand out more
phylogroupFrame$Phylogroup[grep("Shigella", phylogroupFrame$Genome)] <- "Shigella"
phylogroupFrame$Phylogroup[grep("Ancient", phylogroupFrame$Genome)] <- "Ancient"
phylogroupFrame$Phylogroup[is.na(phylogroupFrame$Phylogroup)] <- "Unknown"

# Now to read in the Core Gene statuses
coreGenecounts <- read.delim("../Virulence/CoreGeneCounts.tab",)
coreGenecounts$Genome[108] <- "AncientChromosome"
virGenecounts <- read.delim("../Virulence/VirulenceGeneCounts.tab",)
virGenecounts$Genome[108] <- "AncientChromosome"
colnames(virGenecounts) <- c("Genome", "Virulence")
colnames(coreGenecounts) <- c("Genome", "Core")
FinalcoreGene <- phylogroupFrame %>% left_join(coreGenecounts) %>% mutate(Frac = Core/3144) %>% left_join(virGenecounts)
rownames(FinalcoreGene) <- FinalcoreGene$Genome
FinalcoreGene <- FinalcoreGene[tree@phylo$tip.label,]

# Going to try replacing the labels with the Phylogroups
p1 <- ggtree(tree, right = T, mrsd = "2019-08-02") %<+% FinalcoreGene +
	geom_tippoint(mapping = aes(colour = Phylogroup), size = 1.5, key_glyph = draw_key_point, align = T) +
	geom_range("height_0.95_HPD", color = "#f8333c", size = 2, alpha = 0.5) +
	#geom_text(aes(label=node))+
	geom_text(aes(label=round(as.numeric(posterior), 2), x = branch), vjust=0, size = 2.5) +
#	geom_treescale(linesize = 1, offset = 2) +
        theme_tree2() + geom_rootedge(rootedge = 0.01) +
	scale_color_manual(values = colour, name = "Phylogroup") +
	scale_x_continuous(breaks = scales:::pretty_breaks(n = 5), minor_breaks = scales:::pretty_breaks(n = 15)) +
	theme(legend.position = "bottom", panel.grid.major = element_line(color = "black", size = .2),
	      panel.grid.minor = element_line(color = "grey", size = .2),
	      panel.grid.major.y = element_blank(),
              panel.grid.minor.y = element_blank()) +
	guides(colour = guide_legend(nrow = 1)) 

p1 <- p1 + geom_cladelabel(node = c(197,193), label = "A", align = T) + geom_cladelabel(node = 149, label = "B1", align = T) + geom_cladelabel(node = 112, label = "B2", align = T) +
	geom_cladelabel(node = 152, "C", align = F) + geom_cladelabel(130, "D", align = T) + geom_cladelabel(142, "E", align = T) + geom_cladelabel(137, "F", align = T) + geom_cladelabel(127, "G", align = F)
p1
ggsave(plot = p1, "~/Documents/University/EcoliPaperV2/Figures/BeastTree.pdf", width = 8, height = 6)

p2 <- ggtree(tree, right = T, aes(color = Frac)) %<+% FinalcoreGene +
	theme_tree() + geom_rootedge(rootedge = 0.01) +
	scale_color_continuous(name = "Core Genome Presence", limits = c(0.8,1),oob=scales::squish, low = "#f8333c", high = "#007dba")
p2 <- p2 + new_scale_color() + geom_tippoint(mapping = aes(colour = Phylogroup), size = 1, show.legend = F) +
	scale_color_manual(values = colour, name = "Phylogroup")

p3 <- ggtree(tree, right = T, aes(color = Virulence)) %<+% FinalcoreGene +
	theme_tree() + geom_rootedge(rootedge = 0.01) +
	scale_color_continuous(name = "Virulence Genes", oob=scales::squish, low = "#f8333c", high = "#007dba")
p3 <- p3 + new_scale_color() + geom_tippoint(mapping = aes(colour = Phylogroup), size = 1, show.legend = F) +
	scale_color_manual(values = colour, name = "Phylogroup")
p23 <- ggarrange(p2,p3, ncol = 2, labels = c("B", "C"), align = "hv")

p123 <- ggarrange(p1, p23, nrow =2, labels = c("A",""), align = "hv", heights = c(4,2), common.legend = T, legend = "bottom")
p123
ggsave(plot = p123, "~/Documents/University/EcoliPaperV2/Figures/BeastTrees.pdf", width = 12, height = 8)

write.tree(tree@phylo, file = "MidpointTreemerBEASTTree.nwk")

########################
### The Tempest Plot ###
########################

# Now to get the Tempest Plot
dat <- as_tibble(read.delim("TreemerMidpointTempestBeast.tab"))
# Getting the Phylogroups plotted
phylogroupFrame <- read.delim("~/Downloads/AllFinal_phylogroups.txt", header = F, stringsAsFactors = F)[,c(1,5)]
colnames(phylogroupFrame) <- c("Genome", "Phylogroup")
phylogroupFrame[456,] <- list("O104","B1")

phylogroupFrame$Genome <- gsub(".fna|.fasta|.fasta.ref","",phylogroupFrame$Genome) # Removing the last bit

phylogroupFrame <- phylogroupFrame[match(dat$tip, phylogroupFrame$Genome),] # Getting them in the Right Order
notFound <- which(sapply(phylogroupFrame$Genome, is.na)) # Which were not present in the Phylogrouping?
phylogroupFrame$Genome[notFound] <- dat$tip[notFound]
phylogroupFrame$Phylogroup[notFound] <- c("Unknown")

# Removing Phylogroup Info from Shigella to make it stand out more
phylogroupFrame$Phylogroup[grep("Shigella", phylogroupFrame$Genome)] <- "Shigella"
phylogroupFrame$Phylogroup[grep("Ancient", phylogroupFrame$Genome)] <- "Ancient"
phylogroupFrame$Phylogroup[is.na(phylogroupFrame$Phylogroup)] <- "Unknown"

dat$Phylogroup <- phylogroupFrame$Phylogroup

tempest <- dat %>% ggplot(aes(x = date, y  = distance)) + geom_smooth(method = "lm", show.legend = F, colour = "black") +
	geom_point(alpha = 0.75, aes(colour = Phylogroup), show.legend = T) + theme_classic() + ylab("Root to Tip Divergence") + xlab("Year") +
	scale_color_manual(values = colour, name = "Phylogroups")

model <- lm(distance ~date, data = dat)
summary(model)
r2 <- round(summary(model)$adj.r.squared,3)
tmp <- summary(model)$fstatistic
pval <- round(pf(tmp[1],tmp[2],tmp[3], lower.tail = F),3)

tempest <- tempest + annotate(geom = "text", x = 1700, y = 415,label = bquote(R[adj]^2 == .(r2))) +
	annotate(geom = "text", x = 1700, y = 400, label = bquote(P == .(pval)))

ggsave(tempest, file = "~/Documents/University/EcoliPaperV2/Figures/MidpointTempestBeast.pdf", width = 8, height = 6)
