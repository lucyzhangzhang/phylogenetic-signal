library(ape)
# library(seqinr)
# library(phangorn)
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
# seq <- read.phyDat(phylip.files[1], format = "interleaved")
# dm <- dist.ml(seq)
# UPGMA <- upgma(dm)
# NJ <- NJ(dm)
# fit <- pml(NJ, data = seq)
# 
# gtr <- update(fit, k = 4, inv = 0.2)
# gtr <- optim.pml(gtr, model = "GTR", optInv = T, optGamma = T,
#                  rearrangement = "stochastic", control = pml.control(trace = 0 ))
# gtr
# 
# gtr$tree$tip.label <- paste0(organisms,":Cyclophilin")
# bs <- bootstrap.pml(gtr, bs = 100, control = pml.control(trace = 0))
# 
# plot(gtr)
# 
#     plotBS(midpoint(gtr$tree), bs, gtr$edge.length, p = 50, type = "p")
# replace <- as.data.frame(cbind(tip=plot$tip.label[[1]], org=paste0(organisms,":Cyclophilin")))

library(phytools)
#useless

#best tree
bipartition <- read.tree("BLtrees")
av.tr <- averageTree(bipartition, method = "symmetric.difference")
plotTree(root(av.tr, outgroup = c("Zmays", "Osat"), resolve.root = T))  
bipart <- read.tree("bipart")
av.tr <- averageTree(bipart, method = "symmetric.difference")
write.tree(av.tr, file = "correctTopo")
plotTree(root(av.tr, outgroup = c("Zmays", "Osat"), resolve.root = T))  

root <- lapply(bipart, root, c("Zmays", "Osat"), resolve.root = T)
sp <- speciesTree(root, FUN = mean)
plotTree(speciesTree(bipart), resolve.root= T)

# bipart <- read.tree("RAxML_bipartitionsBranchLabels.bipartition")
# plotTree(bipart)

all <- read.tree("RAxML_bipartitionsBranchLabels.all")
rooted <- root(all, outgroup = c("Osat"), resolve.root = T)
plotTree(rooted)

chro <- chronos(rooted)
d <- c("M. truncatula", "Z. mays",  "O. sativa",  "B. napus", "B. rapa",  
        "E. salsugineum",  "A. thaliana",  "C. rubella",  "S. lycopersicum")
chro$tip.label <- d
pdf("chronotree.pdf", height = 4, width = 6)
plot(chro)
add.scale.bar(pos = 3)
dev.off()

#phylogenetic signal
library(picante)
library(adephylo)
library(phylobase)
library(geiger)
library(dplyr)
library(phylosignal)

#traits generated from moleculatTrait.sh
traits <- read.table("../molecularTrait.out", header = T)
head(traits,7)

#trees generated from bootstrapping with RAxML
tree <- all

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
    phylotraits <- phylo4d(tree, tip.data = dat[,2:ncol(dat)])
    
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
test <- phyloSignal(phylotraits, methods = c("Lambda", "K", "I"), reps = 1000)
round(test$stat,digits=3 )## REQUIRES BRANCH LENGTHS!!
round(test$pvalue, digits = 4)

gc.corr <- phyloCorrelogram(phylotraits, trait = "GC")
plot(gc.corr)
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

#columns 3 and 4 are the signals and P-values
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


write.table(Alltrait, "../Allt.tab", quote = F, sep = "\t", row.names = F)
write.table(Alltrait, "../Alltable.sand", quote = F, sep = " & ", row.names = F, eol = " \\\\ \n")

#moran's I auto correlation graphs

h0 <- -1/(length(organisms) - 1) #h0 is the value we expect without phylogenetic signal
#caitlin's function, modified for expected value
## extracting data from phyloCorrelogram for custom plotting
correlogram_cust <- function(gc.corr, orf.corr, length.corr, exon.corr){
  h0 <- -1/(length(organisms) - 1) #h0 is the value we expect without phylogenetic signal
#moran's I expected value, value significantly higher than this has positive spatial autocorrelation
  gc<- data.frame(gc.corr$res, trait=rep("GC", 1000)) %>% 
    `colnames<-`(c("dist", "min", "max", "cor", "trait"))  %>% 
    mutate(sig = case_when(min > h0 & max > h0 ~ "Significant",
                           min < h0 & max < h0 ~ "Significant", 
                           TRUE ~ "Insignificant"))
  orf <- data.frame(orf.corr$res, trait=rep("ORF", 1000))  %>% 
    `colnames<-`(c("dist", "min", "max", "cor", "trait"))  %>% 
    mutate(sig = case_when(min > h0 & max > h0 ~ "Significant",
                           min < h0 & max < h0 ~ "Significant", 
                           TRUE ~ "Insignificant"))
  length <- data.frame(length.corr$res, trait=rep("length", 1000)) %>% 
    `colnames<-`(c("dist", "min", "max", "cor", "trait"))  %>% 
    mutate(sig = case_when(min > h0 & max > h0 ~ "Significant",
                           min < h0 & max < h0 ~ "Significant", 
                           TRUE ~ "Insignificant"))
  exon <- data.frame(exon.corr$res, trait=rep("exon", 1000)) %>% 
    `colnames<-`(c("dist", "min", "max", "cor", "trait")) %>% 
    mutate(sig = case_when(min > h0 & max > h0 ~ "Significant",
                           min < h0 & max < h0 ~ "Significant", 
                           TRUE ~ "Insignificant"))
  plotdata <- rbind(gc, orf, length, exon)
  return(plotdata)
  }

#the features to test, extracted with moleculatTrair.sh
feat <- c("exonNum", "ORFLen", "mRNAlen", "GC")

#initate empty data frames
H.feat <- data.frame(org = organisms)
P.feat <- data.frame(org = organisms)
IPS1.feat <- data.frame(org = organisms)

for (i in 1:length(feat)) {
    #housekeeping average
    house <- c("Cyclophilin", "EF1a", "Hsp90") 
    H <- traits %>% filter(grepl(paste0(house, collapse = "|"), Gene))
    #plus 3 because of how the data is structured
    H.feat <- cbind(H.feat,  H[, i + 3])
    
    #phosphate genes average
    pho <- c("SIZ1", "PHR1", "PHO2")
    P <- traits %>% filter(grepl(paste0(pho, collapse = "|"), Gene))
    P.feat <- cbind(P.feat,  P[, i + 3])

    #IPS1
    IPS <- traits %>% filter(grepl("IPS1", Gene))
    IPS1.feat <- cbind(IPS1.feat, IPS[, i + 3])
}

    #make objects with tree and trait data
    H.mean <- aggregate(H.feat[,2:length(H.feat)], list(H.feat$org), mean)
    colnames(H.mean) <- c("org", feat)
    H.phylo <- phylo4d(tree, tip.data = H.mean[, 2:5])

    P.mean <- aggregate(P.feat[,2:length(P.feat)], list(P.feat$org), mean)
    colnames(P.mean) <- c("org", feat)
    P.phylo <- phylo4d(tree, tip.data = P.mean[, 2:5])

    IPS1.mean <- aggregate(IPS1.feat[,2:length(IPS1.feat)], list(IPS1.feat$org), mean)
    colnames(IPS1.mean) <- c("org", feat)
    IPS1.phylo <- phylo4d(tree, tip.data = IPS1.mean[, 2:5])

#phylogenetic signal testing
H.test <- phyloSignal(H.phylo, methods = c("Lambda", "K", "I"), reps = 1000)
round(H.test$stat,digits=3 )## REQUIRES BRANCH LENGTHS!!
round(H.test$pvalue, digits = 4)

P.test <- phyloSignal(P.phylo, methods = c("Lambda", "K", "I"), reps = 1000)
round(P.test$stat,digits=3 )## REQUIRES BRANCP.LENGTP.!!
round(P.test$pvalue, digits = 4)

IPS1.test <- phyloSignal(IPS1.phylo, methods = c("Lambda", "K", "I"), reps = 1000)
round(IPS1.test$stat,digits=3 )## REQUIRES BRANCIPS1.LENGTIPS1.!!
round(IPS1.test$pvalue, digits = 4)

#Correlogram of each feature (and together)
# H.gc.corr <- phyloCorrelogram(H.phylo, trait = "GC")
# H.orf.corr <- phyloCorrelogram(H.phylo, trait = "ORFLen")
# H.mrna.corr <- phyloCorrelogram(H.phylo, trait = "mRNAlen")
# H.exon.corr <- phyloCorrelogram(H.phylo, trait = "exonNum")
H.all.corr <- phyloCorrelogram(H.phylo, trait = feat)

H.plotDat <- correlogram_cust(H.gc.corr, H.orf.corr, H.mrna.corr, H.exon.corr)

# P.gc.corr <- phyloCorrelogram(P.phylo, trait = "GC")
# P.orf.corr <- phyloCorrelogram(P.phylo, trait = "ORFLen")
# P.mrna.corr <- phyloCorrelogram(P.phylo, trait = "mRNAlen")
# P.exon.corr <- phyloCorrelogram(P.phylo, trait = "exonNum")
P.all.corr <- phyloCorrelogram(P.phylo, trait = feat)

P.plotDat <- correlogram_cust(P.gc.corr, P.orf.corr, P.mrna.corr, P.exon.corr)

# IPS1.gc.corr <- phyloCorrelogram(IPS1.phylo, trait = "GC")
# IPS1.orf.corr <- phyloCorrelogram(IPS1.phylo, trait = "ORFLen")
# IPS1.mrna.corr <- phyloCorrelogram(IPS1.phylo, trait = "mRNAlen")
# IPS1.exon.corr <- phyloCorrelogram(IPS1.phylo, trait = "exonNum")
IPS1.all.corr <- phyloCorrelogram(IPS1.phylo, trait = feat)

IPS1.plotDat <- correlogram_cust(IPS1.gc.corr, IPS1.orf.corr, IPS1.mrna.corr, IPS1.exon.corr)

library(ggplot2)
{ ggplot(H.plotDat, aes(x = dist, y = cor, color = sig, group = trait)) +  
  geom_point(aes(shape=sig)) + geom_line() + 
  labs(x = "Distance", y = "Correlation") + 
  geom_ribbon(aes(ymin=min,ymax=max),alpha=0.3, color = NA) +
  facet_grid(~trait) +
  geom_hline(yintercept=(-1/H0), color="black")
}   
all.plotdata <- rbind(data.frame(H.plotDat, type=rep("house", 4000)), 
                      data.frame(P.plotDat, type=rep("pho", 4000)),
                      data.frame(IPS1.plotDat, type=rep("ips", 4000))) 

all.plotdata$trait <- factor(all.plotdata$trait, levels = c("ORF", "GC", "exon", "length"))
all.plotdata$sig <- factor(all.plotdata$sig, levels = c("Significant", "Insignificant"))

# legend_ord <- levels(with(all.plotdata, reorder(c("Significant", "Not Significant"))))
                           
feature_names <- as_labeller(c(
  'GC'="GC%",
  'length'="mRNA length",
  'exon'="Exon number",
  'ORF'="ORF length",
  'house'="Housekeeping genes",
  'pho'="PSR",
  'ips'="XLOC_008023"
))


correlogram <- { ggplot(all.plotdata, aes(x = dist, y = cor, color = sig, group = trait)) + 
  geom_point() + 
  geom_ribbon(aes(ymin=min,ymax=max),alpha=0.15, color = NA) +
  facet_grid(type ~ trait, labeller=feature_names) +
  geom_hline(yintercept=h0, color="black") + theme_bw() + ylab("Phylogenetic correlation") +
  xlab("Phylogenetic distance") + 
 #scale_color_manual(values=c("#ed7a53","#A7c6da")) +
  theme(legend.title=element_blank(), legend.position = "bottom")
}

# gridplot(H.phylo)
#really slow on remote X11 forwarding
 correlogram
ggsave("correlogram.pdf", correlogram, height = 7, width = 9)  

test.values <- rbind(H.test$stat, H.test$pvalue, 
P.test$stat, P.test$pvalue, 
IPS1.test$stat, IPS1.test$pvalue)
test.values <- format(test.values, digits = 3, scientific = F)
write.table(test.values, "physig.vals", sep = " & ", quote = F, eol = " \\\\\n")

#statistical analysis
