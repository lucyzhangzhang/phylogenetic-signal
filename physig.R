#making trees

library(ape)
library(phangorn)
library(seqinr)

setwd("~/scratch/physig/alignment/")

#import files
genes <- c("Cyclophilin", "EF1a", "Hsp90",
           "PHO2", "PHR1", "SIZ1")
#read file into a vector
organisms <- scan("../config", character(), quote = "")
phylip.files <- (paste0(genes,".fa.phy"))

#testing models
seq <- read.phyDat(phylip.files[1], format = "interleaved")
dm <- dist.ml(seq)
UPGMA <- upgma(dm)
NJ <- NJ(dm)
fit <- pml(NJ, data = seq)

gtr <- update(fit, k = 4, inv = 0.2)
gtr <- optim.pml(gtr, model = "GTR", optInv = T, optGamma = T,
                 rearrangement = "stochastic", control = pml.control(trace = 0 ))
gtr

bs <- bootstrap.pml(gtr, bs = 1000, control = pml.control(trace = 0))

gtrM <- gtr
gtrM$tree[3] <- list(paste0(organisms,":Cyclophilin"))

plotBS(midpoint(gtrM$tree), bs, gtrM$edge.length, p = 50, type = "p")
