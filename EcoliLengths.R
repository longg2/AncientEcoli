library(dplyr)
library(parallel)
library(pbapply)

op <- pboptions(type = "timer")
ncores = 10
depthDf <-pblapply(list.files("EcoliChromosomeLengths", full.names = T), cl = ncores,function(f){
		tmp <- as_tibble(read.table(f, header = F, col.names = c("Chrom", "Length")))
		if(nrow(tmp) == 1){
			return(tmp)
		}else{
			return(NA)
		}
	})

depthDf[!is.na(depthDf)] %>% bind_rows() %>% summarize(Mean = mean(Length), SD = sd(Length), ConfInt = qnorm(0.975)*SD/sqrt(length(Length)), hi = Mean + ConfInt, lo = Mean - ConfInt)
