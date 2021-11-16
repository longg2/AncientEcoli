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

ann_colors = list(ST = c("Ancient" = colour[1],"1429" = colour[6], "325" = colour[7], "3630" = colour[8], "399" = colour[12], "4995" = colour[18], "Other" = colour[20]))

# Basic Tree
tree <- read.newick(file = "ST4995Only/clean.core.fasta.contree")
tree <- phytools::reroot(tree, interactive = T)
Tiplabels <- tree$tip.label

FinalcoreGene <- data.frame(ST = rep("4995", length(Tiplabels)))
rownames(FinalcoreGene) <- Tiplabels
FinalcoreGene$Genome <- Tiplabels
FinalcoreGene$ST[1] <- "Ancient"
FinalcoreGene$ST[22] <- "325"

##################################################
# Let's colour these based on their PhyloGroups ##
##################################################
# Going to try replacing the labels with the Phylogroups

test <- full_join(as_tibble(tree), FinalcoreGene, by = c("label" = "Genome"))
tree2 <- as.treedata(test)

p3 <- ggtree(tree2, right = T) %<+% FinalcoreGene +
       	geom_nodepoint(size = 1.5, colour = ifelse(as.numeric(tree$node.label) >= 90, "black", ifelse(as.numeric(tree$node.label) >= 50, "grey",NA)),
		     shape = "square") +
	geom_tippoint(mapping = aes(colour = ST), size = 1.5) +
	#geom_tiplab(align = T) + 
	geom_treescale(linesize = 1, offset = 2) + theme_tree() + geom_rootedge(rootedge = 0.0005) +
	scale_color_manual(values = ann_colors$ST, name = "Sequence Type") +
	theme(legend.position = "bottom")# + guides(colour = guide_legend(nrow = 1))

save(p3, file = "ST4995Only.RData")

#write.tree(tree, file = "MidpointTreemerTree.nwk")

########################
### The Tempest Plot ###
########################

# Now to get the Tempest Plot
dat <- as_tibble(read.delim("ST4995Only/TempestResultsML.tab")) %>% filter(date > 0)
# Getting the Phylogroups plotted

dat <- dat %>% left_join(FinalcoreGene, by = c("tip" = "Genome"))

tempest <- dat %>% ggplot(aes(x = date, y  = distance)) + geom_smooth(method = "lm", show.legend = F, colour = "black") +
	geom_point(alpha = 0.75, aes(colour = ST), show.legend = F) + theme_bw() + ylab("Root to Tip Divergence") + xlab("Year") +
	scale_color_manual(values = ann_colors$ST)

model <- lm(distance ~date, data = dat)
summary(model)
r2 <- round(summary(model)$adj.r.squared,3)
tmp <- summary(model)$fstatistic
pval <- round(pf(tmp[1],tmp[2],tmp[3], lower.tail = F),3)

tempestST4995ONLY <- tempest + annotate(geom = "text", x = 1700, y = 0.005,label = bquote(R[adj]^2 == .(r2))) +
	annotate(geom = "text", x = 1700, y = 0.00475, label = bquote(P == .(pval))) +
	scale_y_continuous(breaks = scales:::pretty_breaks(n = 10)) + theme(legend.position = "bottom")
save(tempestST4995ONLY, file = "TempestST4995.RData")
