library(ape)
library(seqinr)
library(phangorn)

#making trees
setwd("~/scratch/physig/alignment/")

#import files
genes <- c("Cyclophilin", "EF1a", "Hsp90",
           "IPS1", "PHO2", "PHR1", "SIZ1")
#read file into a vector
organisms <- scan("../config", character(), quote = "")
phylip.files <- paste0(genes,".fa.phy")
fasta.files <- paste0(genes,".fa.aln")

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

gtr$tree$tip.label <- paste0(organisms,":Cyclophilin")
# bs <- bootstrap.pml(gtr, bs = 100, control = pml.control(trace = 0))

# plot(gtr)

# plotBS(midpoint(gtr$tree), bs, gtr$edge.length, p = 50, type = "p")
replace <- as.data.frame(cbind(tip=plot$tip.label[[1]], org=paste0(organisms,":Cyclophilin")))
# plot$tip.label[[1]] <- replace[[2]][match(plot$tip.label, replace[[1]])]
# plot$tip.label  <- sapply(plot$tip.label, function(x) parse(text=x))

#visualize trees
library(ggtree)

tree.bs <- read.newick("../raxml/raxml/RAxML_bipartitions.tres")
tree.plot <- { ggtree(tree.bs, aes(x, y), layout = "rectangular")  + geom_tree() +
  geom_tiplab() + ggplot2::xlim(0, 0.04) +
  geom_treescale(x=0, y=8.5) + geom_nodelab(vjust=-0.17, hjust=1.4) }
tree.plot
ggsave("IPS1.png", IPS1_tree, height=3, width=4, dpi = 200)


#phylogenetic signal
library(picante)
library(adephylo)
library(phylobase)
library(geiger)
library(phytools)
library(dplyr)

#traits generated from moleculatTrait.sh
traits <- read.table("../molecularTrait.out", header = T)
head(traits,7)

#trees generated from bootstrapping with RAxML
tree <- read.tree("../raxml/raxml/RAxML_bipartitionsBranchLabels.tres")

#get first 10 characters of ID because of stupid phylip :((((((
traits <- cbind(traits, phyID=substr(traits[ ,3], 1, 10)) 

Run <- {
Moran <- NULL
Lambda <- NULL
BlomK <- NULL

#test pipeline
for (g in 1:length(genes)) {
    dat <- dplyr::filter(traits, grepl(genes[g], traits$Gene))
    rownames(dat) <- dat[,1]
    dat <- dat[ ,c(1, 4:7)]
    phylotraits <- phylo4d(tree, dat[,2:ncol(dat)])
    
    geneName <- genes[g]
    
    #Moran's I
    #prox <- proxTips(phylotraits, method = "patristic", normalize = "none")
    moran <- abouheif.moran(phylotraits, method = "Abouheif")
    Moran <- rbind(Moran, cbind(geneName, moran$names, moran$obs, moran$pvalue))
    plot(moran)
    for (i in 2:ncol(dat)) {
        trait <- dat[ ,i]
        names(trait) <- dat$phyID 
        traitName <- colnames(dat)[i]
        #     print(traitName)
    
        #Pagel's lambda
        #     print("Lambda...")
        lambda <- phylosig(tree, trait, method = "lambda", test = T, nsim = 999)
        lambda$P <- ifelse(is.nan(lambda$P), NA, lambda$P)
        Lambda <- rbind(Lambda,c(geneName, traitName, lambda$lambda, lambda$P))
       
        #Blomberg's K
        #     print("K...")
        if (var(trait) != 0) {
            #         print("Var is not 0")
            K <- phylosig(tree, trait, method = "K", test = T, nsim = 999)
        }
        else {
            #0 variance fuck shit up
            #         print("Var is 0")
            K <- list(K=NaN, P=NaN)
        }
        K$P <- ifelse(is.nan(K$P), NA, K$P)
        K$K <- ifelse(is.nan(K$K), NA, K$K)
        BlomK  <- rbind(BlomK, c(geneName, traitName, K$K, K$P))
    }
}
    colnames(Moran) <- c("gene", "trait", "Obs", "Pval")
    colnames(BlomK) <- c("gene", "trait", "K", "Pval")
    colnames(Lambda) <- c("gene", "trait", "lambda", "Pval")
}

#get data frames
DF <- {
Moran2 <- data.frame(Moran, stringsAsFactors = F)
Moran2$Pval <- as.numeric(Moran2$Pval)
Moran2$Obs <- as.numeric(Moran2$Obs)
BlomK2 <- data.frame(BlomK, stringsAsFactors = F)
BlomK2$Pval <- as.numeric(BlomK2$Pval)
BlomK2$K <- as.numeric(BlomK2$K)
Lambda2 <- data.frame(Lambda, stringsAsFactors = F)
Lambda2$Pval <- as.numeric(Lambda2$Pval)
Lambda2$lambda <- as.numeric(Lambda2$lambda)
}

stats <- {
    traits <- c("exonNum", "ORFLen", "mRNAlen", "GC")
    StatsM <- data.frame()
    Statsl <- data.frame()
    StatsK <- data.frame()
    for (i in 1:length(traits)) {
        lambda <- dplyr::filter(Lambda2, grepl(traits[i], Lambda2$trait))
        Statsl <- rbind(Statsl, lambda)
        moran <- dplyr::filter(Moran2, grepl(traits[i], Moran2$trait))
        StatsM <- rbind(StatsM, moran)
        blomk <- dplyr::filter(BlomK2, grepl(traits[i], BlomK2$trait))
        StatsK <- rbind(StatsK, blomk)
    }
}
StatsM
Statsl
StatsK
Combined <- cbind(StatsM, Statsl[, 3:4], StatsK[, 3:4])
Combined[,3:8] <- format(Combined[,3:8], digits = 4, scientific = T)
colnames(Combined) <- c("Gene", "Feature", "I", "I:Pval", "l", "l:Pval", "K", "K:Pval")
head(Combined)
write.table(Combined, "../Physig.tab", quote = F, sep = " & ", row.names = F, eol = " \\\\ \n")
write.table(Combined, "../Phyvals.tab", quote = F, sep = " & ", row.names = F, eol = " \\\\ \n")

#Average values
averages <- {
    traits <- c("exonNum", "ORFLen", "mRNAlen", "GC")
Alltrait <- data.frame()
Comp <- cbind(StatsM, Statsl[, 3:4], StatsK[, 3:4]) 
colnames(Comp) <- c("Gene", "Feature", "I", "I:Pval", "l", "l:Pval", "K", "K:Pval")
for (i in 1:length(traits)) {
    #housekeeping average
    house <- c("Cyclophilin", "EF1a", "Hsp90") 
    H.read <- Comp %>% dplyr::filter(grepl(traits[i], Feature)) %>% filter(grepl(paste0(house, collapse = "|"), Gene))
    H.means <- format(colMeans(H.read[,3:8]), digits = 4, scientific = T)
    #phosphate genes average
    pho <- c("SIZ1", "PHR1", "PHO2")
    P.read <- Comp %>% dplyr::filter(grepl(traits[i], Feature)) %>% filter(grepl(paste0(pho, collapse = "|"), Gene))
    P.means <- format(colMeans(P.read[,3:8]), digits = 4, scientific = T)
    
    pre <- data.frame(Type=c("House", "Pho"), trait=traits[i], rbind(H.means, P.means), stringsAsFactors = F)
    Alltrait <- rbind(Alltrait, pre)

}
}
