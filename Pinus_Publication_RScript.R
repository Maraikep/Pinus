#Loading packages
library(vegan)
library(multcomp)
library(compositions)
library(SpiecEasi)
library(igraph)
library(PMCMRplus)
library(PMCMR)
library(NST)
library(seqinr)
library(ALDEx2)
library(writexl)
make.italic <- function(x) as.expression(lapply(x, function(y) bquote(italic(.(y)))))

setwd("C:/Users/marai/Dropbox/UIBK/MICINSNOW/Labstuff/Analysis/Pinus")
#Loading data
load("seqtab_16S.nasrm.RData")
bac <- data.frame(seqtab.nochim.nasrm)
load("seqtab.nochim_ITS.RData")
fun <- data.frame(seqtab.nochim)
rm(seqtab.nochim, seqtab.nochim.nasrm)
load("taxa_16S.nasrm.RData")
bac.tax <- data.frame(taxa.nasrm)
load("taxa_ITS.RData")
fun.tax <- data.frame(taxa_ITS)
rm(taxa.nasrm, taxa_ITS)

#There are some fungal asvs that are unannotated at phylum level. If I drop these, we are loosing some samples, because those do not have reads anymore. We should look at this more carefully.
#summary(as.factor(fun.tax[,2]))
#summary(colSums(fun)[is.na(fun.tax[,2])])
#summary(colSums(fun))
#summary(rowSums(fun))
#summary(rowSums(fun[,is.na(fun.tax[,2])]))

# bac.seqs <- as.list(colnames(bac))
# for(i in 1:length(bac.seqs)){
#   bac.seqs[[i]] <- s2c(bac.seqs[[i]])
# }
#write.fasta(bac.seqs, paste("bacASVs", 1:ncol(bac), sep = "_"), "ROut_publication/bac_seqs.fasta", nbchar = 400)
names(bac) <- paste("bacASVs", 1:ncol(bac), sep = "_")
rownames(bac.tax) <- colnames(bac)


# fun.seqs <- as.list(colnames(fun))
# for(i in 1:length(fun.seqs)){
#   fun.seqs[[i]] <- s2c(fun.seqs[[i]])
# }
#write.fasta(fun.seqs, paste("funASVs", 1:ncol(fun), sep = "_"), "ROut_publication/fun_seqs.fasta", nbchar = 500)
names(fun) <- paste("funASVs", 1:ncol(fun), sep = "_")
rownames(fun.tax) <- names(fun)

meta <- read.delim("metadata.txt", header = T)
meta$Sample[1:84] <- paste(meta$Sample[1:84], meta$snowcover[1:84], sep = "")
#length(which(rownames(bac) %in% meta$bacids))
#length(which(rownames(fun) %in% meta$funids))
#names(meta)


#preparing the bacterial dataset for analysis
bac <- bac[which(rownames(bac) %in% meta$bacids),]
meta.bac <- meta[match(rownames(bac), meta$bacids),]
rownames(bac) <- meta.bac$Sample
#cbind(rownames(bac), meta.bac$bacids)

meta.bac <- meta.bac[which(rowSums(bac) != 0),]
bac <- bac[which(rowSums(bac) != 0),]

bac <- bac[-which(meta.bac$Sample %in% c("PRA4.1no", "PRA4.4no")),]
meta.bac <- meta.bac[-which(meta.bac$Sample %in% c("PRA4.1no", "PRA4.4no")),]

bac.tax <- bac.tax[colSums(bac) > 0,]
bac <- bac[,colSums(bac) > 0]


#Preparing the fungal dataset for analysis
fun <- fun[which(rownames(fun) %in% meta$funids),]
meta.fun <- meta[match(rownames(fun), meta$funids),]
rownames(fun) <- meta.fun$Sample
#cbind(rownames(fun), meta.fun$funids)

meta.fun <- meta.fun[which(rowSums(fun) != 0),]
fun <- fun[which(rowSums(fun) != 0),]

fun.tax <- fun.tax[colSums(fun) > 0,]
fun <- fun[,colSums(fun) > 0]

#################################################################################
#We decided to first focus on the soil samples and leave the meshbag samples aside for the time being. There are two reasons for this: First, their sequencing depth is very low. Second, their communities are very different. If both soil and meshbag samples are analysed jointly, I am afraid these differences will violate the assumptions of the tests applied and also inflate other analysis, especially the networks. Therefore, we decided it is better to analyse first the soil samples and then try to confirm the patterns found with the soil samples by the meshbag samples.

##Bacteria
bac.soil <- bac[which(meta.bac$snowcover != "bmb"),]
meta.bac.soil <- meta.bac[which(meta.bac$snowcover != "bmb"),]

bac.soil <- bac.soil[-which(meta.bac.soil$Sample %in% c("NEGbmb1", "NEGbmb2")),]
meta.bac.soil <- meta.bac.soil[-which(meta.bac.soil$Sample %in% c("NEGbmb1", "NEGbmb2")),]

bac.soil.tax <- bac.tax[colSums(bac.soil) > 0,]
bac.soil <- bac.soil[,colSums(bac.soil) > 0]
#summary(colSums(bac.soil))
#rownames(bac.soil)
#summary(rowSums(bac.soil))

##Fungi
fun.soil <- fun[which(meta.fun$snowcover != "bmb"),]
meta.fun.soil <- meta.fun[which(meta.fun$snowcover != "bmb"),]
#checked that the fungal dataset has the same order as the bacterial one: DONE.

fun.soil.tax <- fun.tax[colSums(fun.soil) > 0,]
fun.soil <- fun.soil[,colSums(fun.soil) > 0]
#summary(colSums(fun.soil))
#summary(rowSums(fun.soil))
rm(bac, bac.seqs, bac.tax, fun, fun.seqs, fun.tax, meta, meta.bac, meta.fun, i)

############################################################################
#Ordination analysis and the influence of experimental factors (location and snowcover) on the compositions.

#PCA
##Bacteria
#tried nmds on clr, rda, cca: They do not converge or give very wired patterns. They cannot be used to visualize these data. the stress for the nmds is rather high, however, it is probably the only choice.
#The PCA has very, very low variance explained. So, better not using it for visualization.
# bac.soil.pca <- prcomp(clr(bac.soil), retx = T)
# bac.pca <- round(bac.soil.pca$sdev^2/sum(bac.soil.pca$sdev^2), digits = 3)

bac.soil.nmds <- metaMDS(bac.soil)
plot(bac.soil.nmds$points, pch = 16, col = as.factor(meta.bac.soil$Factor))

adonis2(dist(clr(bac.soil)) ~ as.factor(meta.bac.soil$Location))
adonis2(dist(clr(bac.soil)) ~ as.factor(meta.bac.soil$Location) + as.factor(meta.bac.soil$snowcover))
#adonis2(dist(clr(bac.soil)) ~ as.factor(meta.bac.soil$Location) + as.factor(meta.bac.soil$snowcover) + as.factor(meta.bac.soil$Location)*as.factor(meta.bac.soil$snowcover)) #interaction is insignificant.
#adonis2(vegdist(bac.soil) ~ as.factor(meta.bac.soil$Location) + as.factor(meta.bac.soil$snowcover) + as.factor(meta.bac.soil$Location)*as.factor(meta.bac.soil$snowcover)) #Bray Curtis confirms the results for all factors.

##Fungi
# fun.soil.pca <- prcomp(clr(fun.soil), retx = T)
# fun.pca <- round(fun.soil.pca$sdev^2/sum(fun.soil.pca$sdev^2), digits = 3)

fun.soil.nmds <- metaMDS(fun.soil)
plot(fun.soil.nmds$points, pch = 16, col = as.factor(meta.fun.soil$Factor))

adonis2(dist(clr(fun.soil)) ~ as.factor(meta.fun.soil$Location))
adonis2(dist(clr(fun.soil)) ~ as.factor(meta.fun.soil$Location) + as.factor(meta.fun.soil$snowcover))


# png("ROut_publication/Ordination_NMDS_fun.png")
# plot(fun.soil.nmds$points, pch = c(16,8)[as.factor(meta.fun.soil$snowcover)], col = c("maroon", "steelblue", "darkgoldenrod")[as.factor(meta.fun.soil$Location)], cex.lab = 1.5, cex.axis = 1.5, cex = 2)
# abline(h = 0, col = "grey50", lty = 3)
# abline(v = 0, col = "grey50", lty = 3)
# dev.off()
# png("ROut_publication/Ordination_NMDS_bac.png")
# plot(bac.soil.nmds$points, pch = c(16,8)[as.factor(meta.bac.soil$snowcover)], col = c("maroon", "steelblue", "darkgoldenrod")[as.factor(meta.bac.soil$Location)], cex.lab = 1.5, cex.axis = 1.5, cex = 2)
# abline(h = 0, col = "grey50", lty = 3)
# abline(v = 0, col = "grey50", lty = 3)
# legend("topleft", legend = c("Kuhtai", "Patscherkofel", "Praxmar", "snow-covered", "snow-free"), col = c("maroon", "steelblue", "darkgoldenrod", "grey", "grey"), pch = c(15, 15, 15, 8, 16), bty = "n", cex = 1.3)
# dev.off()

#########################################################################
#Alpha diversity

bac.soil.rare <- bac.soil
rownames(bac.soil.rare) <- c(paste("K.snowfree.", 1:13, sep = ""), paste("K.snowcover.", 1:12, sep = ""), paste("PA.snowfree.", 1:19, sep = ""), paste("PA.snowcover", 1:14, sep = ""), paste("Pr.snowfree", 1:13, sep = ""), paste("Pr.snowcover", 1:8, sep = ""))
fun.soil.rare <- fun.soil
rownames(fun.soil.rare) <- c(paste("K.snowfree.", 1:13, sep = ""), paste("K.snowcover.", 1:12, sep = ""), paste("PA.snowfree.", 1:19, sep = ""), paste("PA.snowcover", 1:14, sep = ""), paste("Pr.snowfree", 1:13, sep = ""), paste("Pr.snowcover", 1:8, sep = ""))

png("Maraike_ROut/rarefy_bac.png")
rarecurve(bac.soil.rare, step = 100, ylab = "ASVs")
dev.off()
png("Maraike_ROut/rarefy_fun.png")
rarecurve(fun.soil.rare, step = 100, ylab = "ASVs")
dev.off()

##Bacteria
meta.bac.soil$richness <- specnumber(bac.soil)
meta.bac.soil$richness_rare <- specnumber(rrarefy(bac.soil, sample = min(rowSums(bac.soil))))
meta.bac.soil$diversity <- vegan::diversity(bac.soil)
meta.bac.soil$diversity_rare <- vegan::diversity(rrarefy(bac.soil, sample = min(rowSums(bac.soil))))
#boxplot(rowSums(bac.soil) ~ meta.bac.soil$Factor)
#boxplot(meta.bac.soil$richness ~ meta.bac.soil$Factor)
#boxplot(meta.bac.soil$richness ~ meta.bac.soil$Location)
#boxplot(meta.bac.soil$richness ~ meta.bac.soil$snowcover)

kruskal.test(meta.bac.soil$richness ~ meta.bac.soil$Factor)
kruskal.test(meta.bac.soil$richness_rare ~ meta.bac.soil$Factor)
kruskal.test(meta.bac.soil$richness ~ meta.bac.soil$Location)
kruskal.test(meta.bac.soil$richness_rare ~ meta.bac.soil$Location)
kruskal.test(meta.bac.soil$richness ~ meta.bac.soil$snowcover)
kruskal.test(meta.bac.soil$richness_rare ~ meta.bac.soil$snowcover)

aggregate(meta.bac.soil$richness_rare, list(meta.bac.soil$Location), FUN = mean)
aggregate(meta.bac.soil$richness_rare, list(meta.bac.soil$Location), FUN = sd)
aggregate(meta.bac.soil$richness_rare, list(meta.bac.soil$snowcover), FUN = mean)
aggregate(meta.bac.soil$richness_rare, list(meta.bac.soil$snowcover), FUN = sd)


posthoc.kruskal.dunn.test(meta.bac.soil$richness ~ as.factor(meta.bac.soil$Factor))
posthoc.kruskal.dunn.test(meta.bac.soil$richness_rare ~ as.factor(meta.bac.soil$Location))

#aggregate(meta.bac.soil$richness, list(meta.bac.soil$Location), FUN = mean)


kruskal.test(meta.bac.soil$diversity ~ meta.bac.soil$Factor)
kruskal.test(meta.bac.soil$diversity_rare ~ meta.bac.soil$Factor)
kruskal.test(meta.bac.soil$diversity ~ meta.bac.soil$Location)
kruskal.test(meta.bac.soil$diversity_rare ~ meta.bac.soil$Location)
kruskal.test(meta.bac.soil$diversity ~ meta.bac.soil$snowcover)
kruskal.test(meta.bac.soil$diversity_rare ~ meta.bac.soil$snowcover)

aggregate(meta.bac.soil$diversity_rare, list(meta.bac.soil$Location), FUN = mean)
aggregate(meta.bac.soil$diversity_rare, list(meta.bac.soil$Location), FUN = sd)

kruskal.test(rowSums(bac.soil) ~ meta.bac.soil$Factor)
kruskal.test(rowSums(bac.soil) ~ meta.bac.soil$Location)
kruskal.test(rowSums(bac.soil) ~ meta.bac.soil$snowcover)
x11(); plot(rowSums(bac.soil) ~ as.factor(meta.bac.soil$snowcover))

##Fungi
kruskal.test(rowSums(fun.soil) ~ meta.fun.soil$Location)
kruskal.test(rowSums(fun.soil) ~ meta.fun.soil$snowcover)
x11(); plot(rowSums(fun.soil) ~ as.factor(meta.fun.soil$Location))

meta.fun.soil$richness <- specnumber(fun.soil)
meta.fun.soil$richness_rare <- specnumber(rrarefy(fun.soil, 5000))
meta.fun.soil$richness_rare2 <- specnumber(rrarefy(fun.soil, min(rowSums(fun.soil))))
meta.fun.soil$diversity <- vegan::diversity(fun.soil)
meta.fun.soil$diversity_rare <- vegan::diversity(rrarefy(fun.soil, 5000))
meta.fun.soil$diversity_rare2 <- vegan::diversity(rrarefy(fun.soil, min(rowSums(fun.soil))))

#boxplot(meta.fun.soil$richness ~ meta.fun.soil$Factor)
#boxplot(meta.fun.soil$richness ~ meta.fun.soil$Location)
#boxplot(meta.fun.soil$richness ~ meta.fun.soil$snowcover)

kruskal.test(meta.fun.soil$richness ~ meta.fun.soil$Factor)
kruskal.test(meta.fun.soil$richness_rare ~ meta.fun.soil$Factor)
kruskal.test(meta.fun.soil$richness_rare2 ~ meta.fun.soil$Factor)
kruskal.test(meta.fun.soil$richness ~ meta.fun.soil$Location)
kruskal.test(meta.fun.soil$richness_rare ~ meta.fun.soil$Location)
kruskal.test(meta.fun.soil$richness_rare2 ~ meta.fun.soil$Location)
kruskal.test(meta.fun.soil$richness ~ meta.fun.soil$snowcover)
kruskal.test(meta.fun.soil$richness_rare ~ meta.fun.soil$snowcover)
kruskal.test(meta.fun.soil$richness_rare2 ~ meta.fun.soil$snowcover)
aggregate(meta.fun.soil$richness, list(meta.fun.soil$snowcover), FUN = mean)
aggregate(meta.fun.soil$richness_rare2, list(meta.fun.soil$snowcover), FUN = mean)
aggregate(meta.fun.soil$richness, list(meta.fun.soil$snowcover), FUN = sd)
aggregate(meta.fun.soil$richness_rare2, list(meta.fun.soil$snowcover), FUN = sd)

posthoc.kruskal.dunn.test(meta.fun.soil$richness ~ as.factor(meta.fun.soil$Factor))

p <- kruskal.test(rowSums(fun.soil) ~ meta.fun.soil$Factor)$p.value
p <- append(p, kruskal.test(meta.fun.soil$richness ~ meta.fun.soil$Factor)$p.value, length (p))
p <- append(p, kruskal.test(meta.fun.soil$diversity ~ meta.fun.soil$Factor)$p.value, length (p))
p <- append(p, kruskal.test(rowSums(bac.soil) ~ meta.bac.soil$Factor)$p.value, length (p))
p <- append(p, kruskal.test(meta.bac.soil$richness ~ meta.bac.soil$Factor)$p.value, length (p))
p <- append(p, kruskal.test(meta.bac.soil$diversity ~ meta.bac.soil$Factor)$p.value, length (p))
p <- round(p, digits = 4)


kruskal.test(meta.fun.soil$diversity ~ meta.fun.soil$Factor)
kruskal.test(meta.fun.soil$diversity_rare ~ meta.fun.soil$Factor)
kruskal.test(meta.fun.soil$diversity_rare2 ~ meta.fun.soil$Factor)
kruskal.test(meta.fun.soil$diversity ~ meta.fun.soil$Location)
kruskal.test(meta.fun.soil$diversity_rare ~ meta.fun.soil$Location)
kruskal.test(meta.fun.soil$diversity_rare2 ~ meta.fun.soil$Location)
kruskal.test(meta.fun.soil$diversity ~ meta.fun.soil$snowcover)
kruskal.test(meta.fun.soil$diversity_rare ~ meta.fun.soil$snowcover)
kruskal.test(meta.fun.soil$diversity_rare2 ~ meta.fun.soil$snowcover)

#png("ROut_publication/Alpha.png", height = 500)
#par(mfrow = c(2,3), mar = c(2,5,2,1))
#boxplot(rowSums(fun.soil) ~ meta.fun.soil$Factor, ylab = "#Reads", xlab = "", cex.lab = 2, cex.axis = 2, las = 3, ylim = c(0,15000), names = rep("", 6))
#boxplot(meta.fun.soil$richness ~ meta.fun.soil$Factor, ylab = "#ASVs", xlab = "", cex.lab = 2, cex.axis = 2, las = 3, ylim = c(0,100), names = rep("", 6))
#boxplot(meta.fun.soil$diversity ~ meta.fun.soil$Factor, ylab = "Shannon diversity", xlab = "", cex.lab = 2, cex.axis = 2, las = 3, ylim = c(0,4), names = rep("", 6))
#boxplot(rowSums(bac.soil) ~ meta.bac.soil$Factor, ylab = "#Reads", xlab = "", cex.lab = 2, cex.axis = 2, las = 3, ylim = c(0,120000), names = rep("", 6))
#boxplot(meta.bac.soil$richness ~ meta.bac.soil$Factor, ylab = "#ASVs", xlab = "", cex.lab = 2, cex.axis = 2, las = 3, ylim = c(0,2000), names = rep("", 6))
#boxplot(meta.bac.soil$diversity ~ meta.bac.soil$Factor, ylab = "Shannon diversity", xlab = "", cex.lab = 2, cex.axis = 2, las = 3, ylim = c(0,7), names = rep("", 6))
#dev.off()


#png("ROut_publication/Alpha_legend.png")
#plot.new()
#legend("topleft", legend = paste("p = ", p), bty = "n", cex = 2)
#dev.off()



p <- kruskal.test(rowSums(fun.soil) ~ meta.fun.soil$Location)$p.value
p <- append(p, kruskal.test(meta.fun.soil$richness ~ meta.fun.soil$Location)$p.value, length (p))
p <- append(p, kruskal.test(meta.fun.soil$diversity ~ meta.fun.soil$Location)$p.value, length (p))
p <- append(p, kruskal.test(rowSums(bac.soil) ~ meta.bac.soil$Location)$p.value, length (p))
p <- append(p, kruskal.test(meta.bac.soil$richness ~ meta.bac.soil$Location)$p.value, length (p))
p <- append(p, kruskal.test(meta.bac.soil$diversity ~ meta.bac.soil$Location)$p.value, length (p))
p <- round(p, digits = 4)

#png("ROut_publication/Alpha_Location.png", height = 500)
#par(mfrow = c(2,3), mar = c(2,5,2,1))
#boxplot(rowSums(fun.soil) ~ meta.fun.soil$Location, ylab = "#Reads", xlab = "", cex.lab = 2, cex.axis = 2, las = 3, ylim = c(0,15000), names = rep("", 3))
#boxplot(meta.fun.soil$richness ~ meta.fun.soil$Location, ylab = "#ASVs", xlab = "", cex.lab = 2, cex.axis = 2, las = 3, ylim = c(0,100), names = rep("", 3))
#boxplot(meta.fun.soil$diversity ~ meta.fun.soil$Location, ylab = "Shannon diversity", xlab = "", cex.lab = 2, cex.axis = 2, las = 3, ylim = c(0,4), names = rep("", 3))
#boxplot(rowSums(bac.soil) ~ meta.bac.soil$Location, ylab = "#Reads", xlab = "", cex.lab = 2, cex.axis = 2, las = 3, ylim = c(0,120000), names = rep("", 3))
#boxplot(meta.bac.soil$richness ~ meta.bac.soil$Location, ylab = "#ASVs", xlab = "", cex.lab = 2, cex.axis = 2, las = 3, ylim = c(0,2000), names = rep("", 3))
#boxplot(meta.bac.soil$diversity ~ meta.bac.soil$Location, ylab = "Shannon diversity", xlab = "", cex.lab = 2, cex.axis = 2, las = 3, ylim = c(0,7), names = rep("", 3))
#dev.off()


#png("ROut_publication/Alpha_Location_legend.png")
#plot.new()
#legend("topleft", legend = paste("p = ", p), bty = "n", cex = 2)
#dev.off()


p <- kruskal.test(rowSums(fun.soil) ~ meta.fun.soil$snowcover)$p.value
p <- append(p, kruskal.test(meta.fun.soil$richness ~ meta.fun.soil$snowcover)$p.value, length (p))
p <- append(p, kruskal.test(meta.fun.soil$diversity ~ meta.fun.soil$snowcover)$p.value, length (p))
p <- append(p, kruskal.test(rowSums(bac.soil) ~ meta.bac.soil$snowcover)$p.value, length (p))
p <- append(p, kruskal.test(meta.bac.soil$richness ~ meta.bac.soil$snowcover)$p.value, length (p))
p <- append(p, kruskal.test(meta.bac.soil$diversity ~ meta.bac.soil$snowcover)$p.value, length (p))
p <- round(p, digits = 4)

#png("ROut_publication/Alpha_Snowcover.png", height = 500)
#par(mfrow = c(2,3), mar = c(2,5,2,1))
#boxplot(rowSums(fun.soil) ~ meta.fun.soil$snowcover, ylab = "#Reads", xlab = "", cex.lab = 2, cex.axis = 2, las = 3, ylim = c(0,15000), names = rep("", 2))
#boxplot(meta.fun.soil$richness ~ meta.fun.soil$snowcover, ylab = "#ASVs", xlab = "", cex.lab = 2, cex.axis = 2, las = 3, ylim = c(0,100), names = rep("", 2))
#boxplot(meta.fun.soil$diversity ~ meta.fun.soil$snowcover, ylab = "Shannon diversity", xlab = "", cex.lab = 2, cex.axis = 2, las = 3, ylim = c(0,4), names = rep("", 2))
#boxplot(rowSums(bac.soil) ~ meta.bac.soil$snowcover, ylab = "#Reads", xlab = "", cex.lab = 2, cex.axis = 2, las = 3, ylim = c(0,120000), names = rep("", 2))
#boxplot(meta.bac.soil$richness ~ meta.bac.soil$snowcover, ylab = "#ASVs", xlab = "", cex.lab = 2, cex.axis = 2, las = 3, ylim = c(0,2000), names = rep("", 2))
#boxplot(meta.bac.soil$diversity ~ meta.bac.soil$snowcover, ylab = "Shannon diversity", xlab = "", cex.lab = 2, cex.axis = 2, las = 3, ylim = c(0,7), names = rep("", 2))
#dev.off()


#png("ROut_publication/Alpha_snowcover_legend.png")
#plot.new()
#legend("topleft", legend = paste("p = ", p), bty = "n", cex = 2)
#dev.off()

##########################################################################
#########Bacteria#############
#colors from the networks:
#Acidobacteriota   burlywood1
#Actinobacteriota    burlywood3
#Bacteroidota    steelblue4
#Chloroflexi   darkolivegreen3
#Myxococcota   darkseagreen
#Planctomycetota   lightgoldenrod1
#Proteobacteria    steelblue1
#Verrucomicrobiota   darkgoldenrod1

#missing: Firmicutes, WPS-2. Colring them lightgoldenrod4 and darkgoldenrod4.

#Barcharts for overview
##Bacteria
farbe <- data.frame(table(bac.soil.tax$Phylum))
farbe$abundance <- aggregate(colSums(bac.soil), list(bac.soil.tax$Phylum), FUN = sum)[,2]

sel <- c("Acidobacteriota", "Actinobacteriota", "Bacteroidota", "Chloroflexi", "Firmicutes", "Myxococcota", "Planctomycetota", "Proteobacteria", "Verrucomicrobiota", "WPS-2")
#col.sel <- c("burlywood1", "burlywood3", "steelblue4", "darkolivegreen3", "lightgoldenrod4", "darkseagreen", "lightgoldenrod1", "steelblue1", "darkgoldenrod1", "darkgoldenrod4")

bac.soil.tax$farbe <- "black"
bac.soil.tax$farbe[-which(bac.soil.tax$Phylum %in% sel)] <- "grey"

for(i in 1:length(sel)){
  bac.soil.tax$farbe[which(bac.soil.tax$Phylum == sel[i])] <- col.sel[i]
}

col.sel <- c("burlywood", "wheat", "lightsteelblue", "darkolivegreen", "lightgoldenrod", "darkseagreen", "khaki", "steelblue", "orange", "darkgoldenrod")
legendcolors <- c()
legendtext <- c()
bac.soil.tax$farbe_order <- "black"
bac.soil.tax$farbe_order[-which(bac.soil.tax$Phylum %in% sel)] <- "grey"
for(i in 1:length(sel)){
  rv <- levels(as.factor(bac.soil.tax$Order[which(bac.soil.tax$Phylum == sel[i])]))
  rvcol <- colorRampPalette(c(paste(col.sel[i], 1, sep = ""), paste(col.sel[i], 4, sep = "")))(length(rv))
  legendcolors <- append(legendcolors, rvcol, length(legendcolors))
  legendcolors <- append(legendcolors, col.sel[i], length(legendcolors))
  legendtext <- append(legendtext, rv, length(legendtext))
  legendtext <- append(legendtext, paste("other", sel[i], sep = " "), length(legendtext))
  bac.soil.tax$farbe_order[which(bac.soil.tax$Phylum == sel[i])] <- col.sel[i]
  for(j in 1:length(rv)){
    bac.soil.tax$farbe_order[which(bac.soil.tax$Order == rv[j])] <- rvcol[j]
  }
}

bac.rel <- t(bac.soil)
for(i in 1:ncol(bac.rel)){
  bac.rel[,i] <- bac.rel[,i]/rowSums(bac.soil)[i]
}
colSums(bac.rel)



#png("ROut_publication/Barplot_tax_bac_proportions.png", width = 1000)
#par(mar = c(12,5,2,2))
#barplot(as.matrix(bac.rel[order(bac.soil.tax$Phylum),]), las = 3, col = bac.soil.tax$farbe[order(bac.soil.tax$Phylum)], border = NA, ylab  = "Relative abundance", cex.axis = 1.3, cex.lab = 1.5, names.arg = meta.bac.soil$bacids)
#dev.off()

# png("ROut_publication/Barplot_bac_legend.png")
# plot.new()
# legend("topleft", legend = c(sel, "Others"), col = c(col.sel, "grey"), pch = 15, bty = "n", cex = 1.5, pt.cex = 3.3)
# dev.off()


bac.factor <- data.frame(matrix(data = NA, ncol = ncol(bac.soil), nrow = length(levels(as.factor(meta.bac.soil$Factor)))))
names(bac.factor) <- names(bac.soil)
rownames(bac.factor) <- levels(as.factor(meta.bac.soil$Factor))
for(i in 1:ncol(bac.factor)){
  bac.factor[,i] <- aggregate(bac.soil[,i], list(meta.bac.soil$Factor), FUN = sum)[,2]
}
summary(rowSums(bac.factor))
summary(rowSums(bac.soil))

#png("ROut_publication/Barplot_tax_bac_factor_reads.png")
#par(mar = c(10,5,2,2))
#barplot(t(bac.factor[,order(bac.soil.tax$Phylum)]), col = bac.soil.tax$farbe[order(bac.soil.tax$Phylum)], border = NA, ylab  = "Relative abundance", cex.axis = 1.3, cex.lab = 1.5, las = 3, names.arg = c("Kuhtai_free", "Kuhtai_covered", "Patscherkofel_free", "Patscherkofel_covered", "Praxmar_free", "Praxmar_covered"))
#dev.off()
#png("ROut_publication/Barplot_tax_bac_factor_subsample.png")
#par(mar = c(10,5,2,2))
#barplot(t(rrarefy(bac.factor[,order(bac.soil.tax$Phylum)], 55000)), col = bac.soil.tax$farbe[order(bac.soil.tax$Phylum)], border = NA, ylab  = "Relative abundance", cex.axis = 1.3, cex.lab = 1.5, cex.names = 1.5, las = 3, names.arg = c("Kuhtai", "Kuhtai", "Patscherkofel", "Patscherkofel", "Praxmar", "Praxmar"))
#dev.off()



#png("ROut_publication/Barplot_tax_bac_factor_reads_order.png")
#par(mar = c(10,5,2,2))
#barplot(t(bac.factor[,order(paste(bac.soil.tax$Phylum, bac.soil.tax$Order, sep = "_"))]), col = bac.soil.tax$farbe_order[order(paste(bac.soil.tax$Phylum, bac.soil.tax$Order, sep = "_"))], border = NA, ylab  = "Relative abundance", cex.axis = 1.3, cex.lab = 1.5, las = 3, names.arg = c("Kuhtai_free", "Kuhtai_covered", "Patscherkofel_free", "Patscherkofel_covered", "Praxmar_free", "Praxmar_covered"))
#dev.off()
#png("ROut_publication/Barplot_tax_bac_factor_subsample_order.png")
#par(mar = c(10,5,2,2))
#barplot(t(rrarefy(bac.factor[,order(paste(bac.soil.tax$Phylum, bac.soil.tax$Order, sep = "_"))], 55000)), col = bac.soil.tax$farbe_order[order(paste(bac.soil.tax$Phylum, bac.soil.tax$Order, sep = "_"))], border = NA, ylab  = "Relative abundance", cex.axis = 1.3, cex.lab = 1.5, cex.names = 1.5, las = 3, names.arg = c("Kuhtai", "Kuhtai", "Patscherkofel", "Patscherkofel", "Praxmar", "Praxmar"))
#dev.off()

#png("ROut_publication/Barplot_bac_legend_order.png")
#plot.new()
#legend("topleft", legend = c(sel, "Others"), col = c(col.sel, "grey"), pch = 15, bty = "n")
#dev.off()

#Tables for overview
bac.factor.rel <- t(bac.factor)
for(i in 1:ncol(bac.factor.rel)){
  bac.factor.rel[,i] <- bac.factor.rel[,i]/rowSums(bac.factor)[i]
}
colSums(bac.factor.rel)

bac.phylum <- data.frame(matrix(data = NA, ncol = length(levels(as.factor(bac.soil.tax$Phylum))), nrow = nrow(bac.factor)))
rownames(bac.phylum) <- rownames(bac.factor)
names(bac.phylum) <- levels(as.factor(bac.soil.tax$Phylum))
for(i in 1:nrow(bac.phylum)){
  bac.phylum[i,] <- aggregate(as.numeric(bac.factor[i,]), list(bac.soil.tax$Phylum), FUN = sum)[,2]
}

bac.phylum.rel <- t(bac.phylum)
for(i in 1:ncol(bac.phylum.rel)){
  bac.phylum.rel[,i] <- bac.phylum.rel[,i]/rowSums(bac.phylum)[i]
}
colSums(bac.phylum.rel)

bac.class <- data.frame(matrix(data = NA, ncol = length(levels(as.factor(bac.soil.tax$Class))), nrow = nrow(bac.factor)))
rownames(bac.class) <- rownames(bac.factor)
names(bac.class) <- levels(as.factor(bac.soil.tax$Class))
for(i in 1:nrow(bac.class)){
  bac.class[i,] <- aggregate(as.numeric(bac.factor[i,]), list(bac.soil.tax$Class), FUN = sum)[,2]
}

bac.class.rel <- t(bac.class)
for(i in 1:ncol(bac.class.rel)){
  bac.class.rel[,i] <- bac.class.rel[,i]/rowSums(bac.class)[i]
}
colSums(bac.class.rel)

bac.genus <- data.frame(matrix(data = NA, ncol = length(levels(as.factor(bac.soil.tax$Genus))), nrow = nrow(bac.factor)))
rownames(bac.genus) <- rownames(bac.factor)
names(bac.genus) <- levels(as.factor(bac.soil.tax$Genus))
for(i in 1:nrow(bac.genus)){
  bac.genus[i,] <- aggregate(as.numeric(bac.factor[i,]), list(bac.soil.tax$Genus), FUN = sum)[,2]
}

bac.genus.rel <- t(bac.genus)
for(i in 1:ncol(bac.genus.rel)){
  bac.genus.rel[,i] <- bac.genus.rel[,i]/rowSums(bac.genus)[i]
}
colSums(bac.genus.rel)

#write_xlsx(list(data.frame(rownames(bac.factor.rel),bac.factor.rel, bac.soil.tax[,2:6]), data.frame(rownames(bac.phylum.rel), bac.phylum.rel), data.frame(rownames(bac.class.rel), bac.class.rel), data.frame(rownames(bac.genus.rel), bac.genus.rel)),col_names = T, path = "ROut_publication/Bac_abundance_tables.xlsx")


#########Fungi#############
#colors from the networks:
#p__Ascomycota    pink
#p__Basidiomycota    maroon
#p__Mucoromycota   thistle
#p__Mortierellomycota    mediumpurple

farbe <- data.frame(table(fun.soil.tax$Phylum))
farbe$abundance <- aggregate(colSums(fun.soil), list(fun.soil.tax$Phylum), FUN = sum)[,2]

sel <- c("p__Ascomycota", "p__Basidiomycota", "p__Mucoromycota", "p__Mortierellomycota", "p__Rozellomycota")
col.sel <- c("pink", "maroon", "thistle", "mediumpurple", "plum4")

fun.soil.tax$farbe <- "black"
fun.soil.tax$farbe[-which(fun.soil.tax$Phylum %in% sel)] <- "grey"

for(i in 1:length(sel)){
  fun.soil.tax$farbe[which(fun.soil.tax$Phylum == sel[i])] <- col.sel[i]
}

col.sel <- c("pink", "maroon", "thistle", "mediumpurple", "plum")
legendcolors <- c()
legendtext <- c()
fun.soil.tax$farbe_order <- "black"
fun.soil.tax$farbe_order[-which(fun.soil.tax$Phylum %in% sel)] <- "grey"
for(i in 1:length(sel)){
  rv <- levels(as.factor(fun.soil.tax$Order[which(fun.soil.tax$Phylum == sel[i])]))
  rvcol <- colorRampPalette(c(paste(col.sel[i], 1, sep = ""), paste(col.sel[i], 4, sep = "")))(length(rv))
  legendcolors <- append(legendcolors, rvcol, length(legendcolors))
  legendcolors <- append(legendcolors, col.sel[i], length(legendcolors))
  legendtext <- append(legendtext, rv, length(legendtext))
  legendtext <- append(legendtext, paste("other", sel[i], sep = " "), length(legendtext))
  fun.soil.tax$farbe_order[which(fun.soil.tax$Phylum == sel[i])] <- col.sel[i]
  for(j in 1:length(rv)){
    fun.soil.tax$farbe_order[which(fun.soil.tax$Order == rv[j])] <- rvcol[j]
  }
}

fun.rel <- t(fun.soil)
for(i in 1:ncol(fun.rel)){
  fun.rel[,i] <- fun.rel[,i]/rowSums(fun.soil)[i]
}
colSums(fun.rel)


#png("ROut_publication/Barplot_tax_fun_proportions.png", width = 1000)
#par(mar = c(12,5,2,2))
#barplot(as.matrix(fun.rel[order(fun.soil.tax$Phylum),]), las = 3, col = fun.soil.tax$farbe[order(fun.soil.tax$Phylum)], border = NA, ylab  = "Relative abundance", cex.axis = 1.3, cex.lab = 1.5, names.arg = meta.fun.soil$funids)
#dev.off()

# png("ROut_publication/Barplot_fun_legend.png")
# plot.new()
# legend("topleft", legend = c("Ascomycota", "Basidiomycota", "Mucoromycota", "Mortierellomycota", "Rozellomycota", "Others"), col = c(col.sel, "grey"), pch = 15, bty = "n", cex = 1.5, pt.cex = 3.3)
# dev.off()


fun.factor <- data.frame(matrix(data = NA, ncol = ncol(fun.soil), nrow = length(levels(as.factor(meta.fun.soil$Factor)))))
names(fun.factor) <- names(fun.soil)
rownames(fun.factor) <- levels(as.factor(meta.fun.soil$Factor))
for(i in 1:ncol(fun.factor)){
  fun.factor[,i] <- aggregate(fun.soil[,i], list(meta.fun.soil$Factor), FUN = sum)[,2]
}
summary(rowSums(fun.factor))
summary(rowSums(fun.soil))

#png("ROut_publication/Barplot_tax_fun_factor_reads.png")
#par(mar = c(10,5,2,2))
#barplot(t(fun.factor[,order(fun.soil.tax$Phylum)]), col = fun.soil.tax$farbe[order(fun.soil.tax$Phylum)], border = NA, ylab  = "Relative abundance", cex.axis = 1.3, cex.lab = 1.5, las = 3, names.arg = c("Kuhtai_free", "Kuhtai_covered", "Patscherkofel_free", "Patscherkofel_covered", "Praxmar_free", "Praxmar_covered"))
#dev.off()
#png("ROut_publication/Barplot_tax_fun_factor_subsample.png")
#par(mar = c(10,5,2,2))
#barplot(t(rrarefy(fun.factor[,order(fun.soil.tax$Phylum)], 10000)), col = fun.soil.tax$farbe[order(fun.soil.tax$Phylum)], border = NA, ylab  = "Relative abundance", cex.axis = 1.3, cex.lab = 1.5, cex.names = 1.5, las = 3, names.arg = c("Kuhtai", "Kuhtai", "Patscherkofel", "Patscherkofel", "Praxmar", "Praxmar"))
#dev.off()


#png("ROut_publication/Barplot_taxorder_fun_factor_reads.png")
#par(mar = c(10,5,2,2))
#barplot(t(fun.factor[,order(paste(fun.soil.tax$Phylum, fun.soil.tax$Order, sep = "|"))]), col = fun.soil.tax$farbe_order[order(paste(fun.soil.tax$Phylum, fun.soil.tax$Order, sep = "|"))], border = NA, ylab  = "Relative abundance", cex.axis = 1.3, cex.lab = 1.5, las = 3, cex.names = 1.5, names.arg = c("Kuhtai", "Kuhtai", "Patscherkofel", "Patscherkofel", "Praxmar", "Praxmar"))
#dev.off()
#png("ROut_publication/Barplot_tax_funorder_factor_subsample.png")
#par(mar = c(10,5,2,2))
#barplot(t(rrarefy(fun.factor[,order(paste(fun.soil.tax$Phylum, fun.soil.tax$Order, sep = "|"))], 10000)), col = fun.soil.tax$farbe_order[order(paste(fun.soil.tax$Phylum, fun.soil.tax$Order, sep = "|"))], border = NA, ylab  = "Relative abundance", cex.axis = 1.3, cex.lab = 1.5, cex.names = 1.5, las = 3, names.arg = c("Kuhtai", "Kuhtai", "Patscherkofel", "Patscherkofel", "Praxmar", "Praxmar"))
#dev.off()

#png("ROut_publication/Barplot_fun_legend_order.png", height = 1500)
#plot.new()
#legend("topleft", legend = c(legendtext, "Others"), col = c(legendcolors, "grey"), pch = 15, bty = "n")
#dev.off()



#Tables for overview
fun.factor.rel <- t(fun.factor)
for(i in 1:ncol(fun.factor.rel)){
  fun.factor.rel[,i] <- fun.factor.rel[,i]/rowSums(fun.factor)[i]
}
colSums(fun.factor.rel)

fun.phylum <- data.frame(matrix(data = NA, ncol = length(levels(as.factor(fun.soil.tax$Phylum))), nrow = nrow(fun.factor)))
rownames(fun.phylum) <- rownames(fun.factor)
names(fun.phylum) <- levels(as.factor(fun.soil.tax$Phylum))
for(i in 1:nrow(fun.phylum)){
  fun.phylum[i,] <- aggregate(as.numeric(fun.factor[i,]), list(fun.soil.tax$Phylum), FUN = sum)[,2]
}

fun.phylum.rel <- t(fun.phylum)
for(i in 1:ncol(fun.phylum.rel)){
  fun.phylum.rel[,i] <- fun.phylum.rel[,i]/rowSums(fun.phylum)[i]
}
colSums(fun.phylum.rel)

fun.class <- data.frame(matrix(data = NA, ncol = length(levels(as.factor(fun.soil.tax$Class))), nrow = nrow(fun.factor)))
rownames(fun.class) <- rownames(fun.factor)
names(fun.class) <- levels(as.factor(fun.soil.tax$Class))
for(i in 1:nrow(fun.class)){
  fun.class[i,] <- aggregate(as.numeric(fun.factor[i,]), list(fun.soil.tax$Class), FUN = sum)[,2]
}

fun.class.rel <- t(fun.class)
for(i in 1:ncol(fun.class.rel)){
  fun.class.rel[,i] <- fun.class.rel[,i]/rowSums(fun.class)[i]
}
colSums(fun.class.rel)

fun.order <- data.frame(matrix(data = NA, ncol = length(levels(as.factor(fun.soil.tax$Order))), nrow = nrow(fun.factor)))
rownames(fun.order) <- rownames(fun.factor)
names(fun.order) <- levels(as.factor(fun.soil.tax$Order))
for(i in 1:nrow(fun.order)){
  fun.order[i,] <- aggregate(as.numeric(fun.factor[i,]), list(fun.soil.tax$Order), FUN = sum)[,2]
}

fun.order.rel <- t(fun.order)
for(i in 1:ncol(fun.order.rel)){
  fun.order.rel[,i] <- fun.order.rel[,i]/rowSums(fun.order)[i]
}
colSums(fun.order.rel)

fun.genus <- data.frame(matrix(data = NA, ncol = length(levels(as.factor(fun.soil.tax$Genus))), nrow = nrow(fun.factor)))
rownames(fun.genus) <- rownames(fun.factor)
names(fun.genus) <- levels(as.factor(fun.soil.tax$Genus))
for(i in 1:nrow(fun.genus)){
  fun.genus[i,] <- aggregate(as.numeric(fun.factor[i,]), list(fun.soil.tax$Genus), FUN = sum)[,2]
}

fun.genus.rel <- t(fun.genus)
for(i in 1:ncol(fun.genus.rel)){
  fun.genus.rel[,i] <- fun.genus.rel[,i]/rowSums(fun.genus)[i]
}
colSums(fun.genus.rel)

#write_xlsx(list(data.frame(rownames(fun.factor.rel),fun.factor.rel, fun.soil.tax[,2:6]), data.frame(rownames(fun.phylum.rel), fun.phylum.rel), data.frame(rownames(fun.class.rel), fun.class.rel), data.frame(rownames(fun.order.rel), fun.order.rel), data.frame(rownames(fun.genus.rel), fun.genus.rel)),col_names = T, path = "ROut_publication/Fun_abundance_tables.xlsx")


#TopASVs - were selected from the overview table based on high abundances of individual ASVs in each sample group.
#############Bacteria#####################
#bac.factor.sel <- data.frame(rrarefy(bac.factor, 55000))
#names(bac.factor.sel) <- names(bac.factor)
#bac.factor.sel <- bac.factor.sel[,which(names(bac.factor.sel) %in% c("bacASVs_5", "bacASVs_7", "bacASVs_1", "bacASVs_21", "bacASVs_3", "bacASVs_6", "bacASVs_12", "bacASVs_10", "bacASVs_4", "bacASVs_9", "bacASVs_15"))]
#On relative abundances summing to 1
bac.factor.nonrare.sel <- bac.factor.rel[which(rownames(bac.factor.rel) %in% c("bacASVs_5", "bacASVs_7", "bacASVs_1", "bacASVs_21", "bacASVs_3", "bacASVs_6", "bacASVs_12", "bacASVs_10", "bacASVs_4", "bacASVs_9", "bacASVs_15")),]

bac.soil.tax.sel <- bac.soil.tax[which(names(bac.factor) %in% c("bacASVs_5", "bacASVs_7", "bacASVs_1", "bacASVs_21", "bacASVs_3", "bacASVs_6", "bacASVs_12", "bacASVs_10", "bacASVs_4", "bacASVs_9", "bacASVs_15")),]

bac.soil.tax.sel <- bac.soil.tax.sel[order(paste(bac.soil.tax.sel$Phylum, bac.soil.tax.sel$Class, bac.soil.tax.sel$Order, bac.soil.tax.sel$Genus, sep = "|")),]
#bac.factor.sel <- bac.factor.sel[,match(rownames(bac.soil.tax.sel), names(bac.factor.sel))]
bac.factor.nonrare.sel <- bac.factor.nonrare.sel[match(rownames(bac.soil.tax.sel), rownames(bac.factor.nonrare.sel)),]

bac.soil.tax.sel$Genus[1] <- "unknown_Acidobacteriales"
bac.soil.tax.sel$Genus[2] <- "unknown_Acidobacteriae"
bac.soil.tax.sel$Genus[6] <- "unknown_Planococcaceae"
bac.soil.tax.sel$Genus[7] <- "unknown_Bacillales"
bac.soil.tax.sel$Genus[8] <- "unknown_Xanthobacteraceae"
bac.soil.tax.sel$Genus[9] <- "unknown_Xanthobacteraceae"
bac.soil.tax.sel$Genus[10] <- "unknown_Gammaproteobacteria"


bac.soil.tax.sel$farbe2 <- c("burlywood1", "burlywood4", "wheat", "wheat2", "wheat4", "lightgoldenrod", "lightgoldenrod3", "steelblue1", "steelblue2", "steelblue4", "orange")


#png("ROut_publication/Barplot_taxorder_bac_factor_rel_topASVs.png")
#par(mar = c(10,5,2,2))
#barplot(bac.factor.nonrare.sel, col = bac.soil.tax.sel$farbe2, ylab  = "Relative abundance", cex.axis = 1.3, cex.lab = 1.5, las = 3, cex.names = 1.5, names.arg = c("Kuhtai", "Kuhtai", "Patscherkofel", "Patscherkofel", "Praxmar", "Praxmar"))
#dev.off()
png("ROut_publication/Barplot_taxorder_bac_factor_reads_topASVs_legend.png")
plot.new()
legend("topleft", bty = "n", legend = c(c("unknown_Acidobacteriales"), c("unkown_Acidobacteriae"), make.italic(c("Acidothermus sp.", "Acidothermus sp.", "Conexibacter sp.")), c("unknown_Planococcaceae"), c("unknown_Bacillales"), c("unknown_Xanthobacteraceae"), c("unknown_Xanthobacteraceae"), c("unknown_Gammaproteobacteria"), make.italic(c("Candidatus Xiphinematobacter"))), pch = 15, col = bac.soil.tax.sel$farbe2, cex = 1.5, pt.cex = 3.3)
dev.off()


#png("ROut_publication/Barplot_taxorder_bac_factor_reads_topASVs.png")
#par(mar = c(10,5,2,2))
#barplot(t(bac.factor.sel)/55000, col = bac.soil.tax.sel$farbe2, border = NA, ylab  = "Relative abundance", cex.axis = 1.3, cex.lab = 1.5, las = 3, cex.names = 1.5, names.arg = c("Kuhtai", "Kuhtai", "Patscherkofel", "Patscherkofel", "Praxmar", "Praxmar"))
#dev.off()


#############Fungi#####################
#fun.factor.sel <- data.frame(rrarefy(fun.factor, 10000))
#names(fun.factor.sel) <- names(fun.factor)
#fun.factor.sel <- fun.factor.sel[,which(names(fun.factor.sel) %in% c("funASVs_2", "funASVs_1", "funASVs_7", "funASVs_12", "funASVs_14", "funASVs_1", "funASVs_4", "funASVs_18", "funASVs_3", "funASVs_8", "funASVs_11", "funASVs_17", "funASVs_5"))]
#On relative abundances summing to 1
fun.factor.nonrare.sel <- fun.factor.rel[which(rownames(fun.factor.rel) %in% c("funASVs_2", "funASVs_1", "funASVs_3", "funASVs_4", "funASVs_5", "funASVs_6", "funASVs_7", "funASVs_8", "funASVs_9", "funASVs_11", "funASVs_12", "funASVs_14", "funASVs_17", "funASVs_18")),]

fun.soil.tax.sel <- fun.soil.tax[which(names(fun.factor) %in% c("funASVs_2", "funASVs_1", "funASVs_3", "funASVs_4", "funASVs_5", "funASVs_6", "funASVs_7", "funASVs_8", "funASVs_9", "funASVs_11", "funASVs_12", "funASVs_14", "funASVs_17", "funASVs_18")),]

fun.soil.tax.sel <- fun.soil.tax.sel[order(paste(fun.soil.tax.sel$Phylum, fun.soil.tax.sel$Class, fun.soil.tax.sel$Order, fun.soil.tax.sel$Genus, sep = "|")),]
#fun.factor.sel <- fun.factor.sel[,match(rownames(fun.soil.tax.sel), names(fun.factor.sel))]
fun.factor.nonrare.sel <- fun.factor.nonrare.sel[match(rownames(fun.soil.tax.sel), rownames(fun.factor.nonrare.sel)),]

#fun.soil.tax.sel$Genus[5] <- "unknown_Tulsanellaceae"
#fun.soil.tax.sel$Genus[1] <- "unknown"
fun.soil.tax.sel$Species[1] <- "unknown"
fun.soil.tax.sel$Species[2] <- "Amanita submembranacea"
fun.soil.tax.sel$Species[3] <- "Hygrocybe conica"
fun.soil.tax.sel$Species[4] <- "Rhizopogon sp."
fun.soil.tax.sel$Species[5] <- "Suillus sp."
fun.soil.tax.sel$Species[6] <- "unknown Tulasnellaceae"
fun.soil.tax.sel$Species[7] <- "Russula decolorans"
fun.soil.tax.sel$Species[8] <- "Russula heterophylla"
fun.soil.tax.sel$Species[9] <- "Russula mustelina"
fun.soil.tax.sel$Species[10] <- "Russula decolorans"
fun.soil.tax.sel$Species[11] <- "Basidioascus sp."
fun.soil.tax.sel$Species[12] <- "Basidioascus sp."
fun.soil.tax.sel$Species[13] <- "Solicoccozyma terricola"
fun.soil.tax.sel$Species[14] <- "Mortierella macrocystis"


#f <- colorRampPalette(c("white", "maroon4"))(9)
#f <- f[-c(1:2)]
#fun.soil.tax.sel$farbe2 <- c("grey", "#FF34B3", "#EB30A5", "#D82C98", "#C4288A", rep("#B1247D", 4), "#9E206F", "#8B1C62", "mediumpurple")
#fun.soil.tax.sel$farbe2 <- c("grey", f[1:4], rep(f[5], 4), f[6:7], "mediumpurple")

#fun.soil.tax.sel$farbe2 <- c("grey", "#efd5df", "#dfacbf", "#cf829f", "#b7446f", "#b03060", "#9e2b56", "#8c264c", "#7b2143", "#691c39", "#581830", "mediumpurple")
fun.soil.tax.sel$farbe2 <- c("grey", "#efd5df", "#dfacbf", "#bf597f", "#8c264c", "#581830", "#ffc2e8", "#ff85d1", "#ff48ba", "#ffc2e8", "#cc298f", "#cc298f", "#7f1a59", "mediumpurple")


#png("ROut_publication/Barplot_taxorder_fun_factor_rel_topASVs.png")
#par(mar = c(10,5,2,2))
#barplot(fun.factor.nonrare.sel, col = fun.soil.tax.sel$farbe2, ylab  = "Relative abundance", cex.axis = 1.3, cex.lab = 1.5, las = 3, cex.names = 1.5, names.arg = c("Kuhtai", "Kuhtai", "Patscherkofel", "Patscherkofel", "Praxmar", "Praxmar"))
#dev.off()
# png("ROut_publication/Barplot_taxorder_fun_factor_reads_topASVs_legend.png")
# plot.new()
# legend("topleft", bty = "n", legend = c(c("unknown_fungus"), make.italic(c("Amanita submembranacea", "Hygrocybe conica", "Rhizopogon sp.", "Suillus sp.")), c("unknown_Tulasnellaceae"), make.italic(c("Russula decolorans", "Russula heterophylla", "Russula mustelina", "Basidioascus sp", "Solicoccozyma terricola", "Mortierella macrocystis"))), pch = 15, col = fun.soil.tax.sel$farbe2[-c(10,12)], cex = 1.5, pt.cex = 3.3)
# dev.off()


#png("ROut_publication/Barplot_taxorder_fun_factor_reads_topASVs.png")
#par(mar = c(10,5,2,2))
#barplot(t(fun.factor.sel)/10000, col = fun.soil.tax.sel$farbe2, border = NA, ylab  = "Relative abundance", cex.axis = 1.3, cex.lab = 1.5, las = 3, cex.names = 1.5, names.arg = c("Kuhtai", "Kuhtai", "Patscherkofel", "Patscherkofel", "Praxmar", "Praxmar"))
#dev.off()


######################################################################
#Venn diagramm - the Core

fun.bi <- fun.soil
fun.bi[fun.bi > 0] <- 1

bac.bi <- bac.soil
bac.bi[bac.bi > 0] <- 1

fun.ls <- vector("list", length = length(levels(as.factor(meta.fun.soil$Location))))
bac.ls <- vector("list", length = length(levels(as.factor(meta.bac.soil$Location))))
names(fun.ls) <- levels(as.factor(meta.fun.soil$Location))
names(bac.ls) <- names(fun.ls)

fungi <- c()
bacti <- c()
for(i in 1:length(fun.ls)){
  rdf.fun <- fun.bi[which(meta.fun.soil$Location == levels(as.factor(meta.fun.soil$Location))[i]),]
  fun.ls[[i]] <- colnames(rdf.fun)[colSums(rdf.fun) > 0]
  fungi <- append(fungi, fun.ls[[i]])
  fun.ls[[i]] <- data.frame(fun.ls[[i]])
  rdf.bac <- bac.bi[which(meta.bac.soil$Location == levels(as.factor(meta.bac.soil$Location))[i]),]
  bac.ls[[i]] <- colnames(rdf.bac)[colSums(rdf.bac) > 0]
  bacti <- append(bacti, bac.ls[[i]])
  bac.ls[[i]] <- data.frame(bac.ls[[i]])
}
#write_xlsx(fun.ls, path = "ROut_publication/Venn_fun.xlsx")
#write_xlsx(bac.ls, path = "ROut_publication/Venn_bac.xlsx")


fungi <- data.frame(table(fungi))
bacti <- data.frame(table(bacti))
# save(bacti, file = "ROut_publication/bacti.RData")
# save(fungi, file = "ROut_publication/fungi.RData")

table(fungi$Freq)
table(bacti$Freq)

aggregate(colSums(fun.rel[which(rownames(fun.rel) %in% fungi$fungi[which(fungi$Freq == 3)]),]), list(meta.fun.soil$Location), FUN = mean)
aggregate(colSums(fun.rel[which(rownames(fun.rel) %in% fungi$fungi[which(fungi$Freq == 3)]),]), list(meta.fun.soil$Location), FUN = sd)

kruskal.test(colSums(fun.rel[which(rownames(fun.rel) %in% fungi$fungi[which(fungi$Freq == 3)]),]), meta.fun.soil$Location)


aggregate(colSums(bac.rel[which(rownames(bac.rel) %in% bacti$bacti[which(bacti$Freq == 3)]),]), list(meta.bac.soil$Location), FUN = mean)
aggregate(colSums(bac.rel[which(rownames(bac.rel) %in% bacti$bacti[which(bacti$Freq == 3)]),]), list(meta.bac.soil$Location), FUN = sd)

kruskal.test(colSums(bac.rel[which(rownames(bac.rel) %in% bacti$bacti[which(bacti$Freq == 3)]),]), meta.bac.soil$Location)

###Barcharts to show the unique
#fun.factor.unique <- rrarefy(fun.factor, 10000)
#fun.factor.unique <- fun.factor.unique[,which(names(fun.factor) %in% fungi$fungi[which(fungi$Freq == 1)])]
fun.tax.unique <- fun.soil.tax[which(names(fun.factor) %in% fungi$fungi[which(fungi$Freq == 1)]),]
fun.factor.rel.unique <- fun.factor.rel[which(rownames(fun.factor.rel) %in% fungi$fungi[which(fungi$Freq == 1)]),]

#bac.factor.unique <- rrarefy(bac.factor, 55000)
#bac.factor.unique <- bac.factor.unique[,which(names(bac.factor) %in% bacti$bacti[which(bacti$Freq == 1)])]
bac.tax.unique <- bac.soil.tax[which(names(bac.factor) %in% bacti$bacti[which(bacti$Freq == 1)]),]
bac.factor.rel.unique <- bac.factor.rel[which(rownames(bac.factor.rel) %in% bacti$bacti[which(bacti$Freq == 1)]),]

#png("ROut_publication/Barplot_tax_funorder_factor_rel_unique.png")
#par(mar = c(10,5,2,2))
#barplot(fun.factor.rel.unique[order(paste(fun.tax.unique$Phylum, fun.tax.unique$Order, sep = "|")),], col = fun.tax.unique$farbe_order[order(paste(fun.tax.unique$Phylum, fun.tax.unique$Order, sep = "|"))], border = NA, ylab  = "Relative abundance", ylim = c(0,1), cex.axis = 1.3, cex.lab = 1.5, cex.names = 1.5, las = 3, names.arg = c("Kuhtai", "Kuhtai", "Patscherkofel", "Patscherkofel", "Praxmar", "Praxmar"))
#dev.off()

#png("ROut_publication/Barplot_tax_bacorder_factor_rel_unique.png")
#par(mar = c(10,5,2,2))
#barplot(bac.factor.rel.unique[order(paste(bac.tax.unique$Phylum, bac.tax.unique$Order, sep = "|")),], col = bac.tax.unique$farbe_order[order(paste(bac.tax.unique$Phylum, bac.tax.unique$Order, sep = "|"))], border = NA, ylab  = "Relative abundance", ylim = c(0,1), cex.axis = 1.3, cex.lab = 1.5, cex.names = 1.5, las = 3, names.arg = c("Kuhtai", "Kuhtai", "Patscherkofel", "Patscherkofel", "Praxmar", "Praxmar"))
#dev.off()


#png("ROut_publication/Barplot_tax_funorder_factor_subsample_unique.png")
#par(mar = c(10,5,2,2))
#barplot(t(fun.factor.unique[,order(paste(fun.tax.unique$Phylum, fun.tax.unique$Order, sep = "|"))]), col = fun.tax.unique$farbe_order[order(paste(fun.tax.unique$Phylum, fun.tax.unique$Order, sep = "|"))], border = NA, ylab  = "Relative abundance", ylim = c(0,10000), cex.axis = 1.3, cex.lab = 1.5, cex.names = 1.5, las = 3, names.arg = c("Kuhtai", "Kuhtai", "Patscherkofel", "Patscherkofel", "Praxmar", "Praxmar"))
#dev.off()



#png("ROut_publication/Barplot_tax_bac_factor_subsample_unique.png")
#par(mar = c(10,5,2,2))
#barplot(t(bac.factor.unique[,order(paste(bac.tax.unique$Phylum, bac.tax.unique$Order, sep = "|"))]), col = bac.tax.unique$farbe[order(paste(bac.tax.unique$Phylum, bac.tax.unique$Order, sep = "|"))], border = NA, ylab  = "Relative abundance", ylim = c(0,55000), cex.axis = 1.3, cex.lab = 1.5, cex.names = 1.5, las = 3, names.arg = c("Kuhtai", "Kuhtai", "Patscherkofel", "Patscherkofel", "Praxmar", "Praxmar"))
#dev.off()


###Barcharts to show the core
#fun.factor.core <- rrarefy(fun.factor, 10000)
#fun.factor.core <- fun.factor.core[,which(names(fun.factor) %in% fungi$fungi[which(fungi$Freq == 3)])]
fun.factor.rel.core <- fun.factor.rel[which(rownames(fun.factor.rel) %in% fungi$fungi[which(fungi$Freq == 3)]),]
fun.tax.core <- fun.soil.tax[which(names(fun.factor) %in% fungi$fungi[which(fungi$Freq == 3)]),]

#bac.factor.core <- rrarefy(bac.factor, 55000)
#bac.factor.core <- bac.factor.core[,which(names(bac.factor) %in% bacti$bacti[which(bacti$Freq == 3)])]
bac.factor.rel.core <- bac.factor.rel[which(rownames(bac.factor.rel) %in% bacti$bacti[which(bacti$Freq == 3)]),]
bac.tax.core <- bac.soil.tax[which(names(bac.factor) %in% bacti$bacti[which(bacti$Freq == 3)]),]


#png("ROut_publication/Barplot_tax_funorder_factor_rel_core.png")
#par(mar = c(10,5,2,2))
#barplot(fun.factor.rel.core[order(paste(fun.tax.core$Phylum, fun.tax.core$Order, sep = "|")),], col = fun.tax.core$farbe_order[order(paste(fun.tax.core$Phylum, fun.tax.core$Order, sep = "|"))], border = NA, ylab  = "Relative abundance", ylim = c(0,1), cex.axis = 1.3, cex.lab = 1.5, cex.names = 1.5, las = 3, names.arg = c("Kuhtai", "Kuhtai", "Patscherkofel", "Patscherkofel", "Praxmar", "Praxmar"))
#dev.off()

#png("ROut_publication/Barplot_tax_bacorder_factor_rel_core.png")
#par(mar = c(10,5,2,2))
#barplot(bac.factor.rel.core[order(paste(bac.tax.core$Phylum, bac.tax.core$Order, sep = "|")),], col = bac.tax.core$farbe_order[order(paste(bac.tax.core$Phylum, bac.tax.core$Order, sep = "|"))], border = NA, ylab  = "Relative abundance", ylim = c(0,1), cex.axis = 1.3, cex.lab = 1.5, cex.names = 1.5, las = 3, names.arg = c("Kuhtai", "Kuhtai", "Patscherkofel", "Patscherkofel", "Praxmar", "Praxmar"))
#dev.off()



#png("ROut_publication/Barplot_tax_funorder_factor_subsample_core.png")
#par(mar = c(10,5,2,2))
#barplot(t(fun.factor.core[,order(paste(fun.tax.core$Phylum, fun.tax.core$Order, sep = "|"))]), col = fun.tax.core$farbe_order[order(paste(fun.tax.core$Phylum, fun.tax.core$Order, sep = "|"))], border = NA, ylab  = "Relative abundance", ylim = c(0,10000), cex.axis = 1.3, cex.lab = 1.5, cex.names = 1.5, las = 3, names.arg = c("Kuhtai", "Kuhtai", "Patscherkofel", "Patscherkofel", "Praxmar", "Praxmar"))
#dev.off()

#png("ROut_publication/Barplot_tax_bac_factor_subsample_core.png")
#par(mar = c(10,5,2,2))
#barplot(t(bac.factor.core[,order(paste(bac.tax.core$Phylum, bac.tax.core$Order, sep = "|"))]), col = bac.tax.core$farbe[order(paste(bac.tax.core$Phylum, bac.tax.core$Order, sep = "|"))], border = NA, ylab  = "Relative abundance", ylim = c(0,55000), cex.axis = 1.3, cex.lab = 1.5, cex.names = 1.5, las = 3, names.arg = c("Kuhtai", "Kuhtai", "Patscherkofel", "Patscherkofel", "Praxmar", "Praxmar"))
#dev.off()

#write.csv(data.frame(fun.factor.rel.core, fun.tax.core), file = "ROut_publication/Core_fun_percentges.csv")
#write.csv(data.frame(bac.factor.rel.core, bac.tax.core), file = "ROut_publication/Core_bac_percentages.csv")


###Trying a heatmap approach to show the abundance change in ASVs
heat.fun.factor.rel.core <- fun.factor.rel.core
for(i in 1:ncol(heat.fun.factor.rel.core)){
  for(j in 1:nrow(heat.fun.factor.rel.core)){
    if(heat.fun.factor.rel.core[j,i] > 0 && heat.fun.factor.rel.core[j,i] < 0.025){
      heat.fun.factor.rel.core[j,i] <- 1
    }else if(heat.fun.factor.rel.core[j,i] >= 0.025 && heat.fun.factor.rel.core[j,i] < 0.05){
      heat.fun.factor.rel.core[j,i] <- 2
    }else if(heat.fun.factor.rel.core[j,i] >= 0.05 && heat.fun.factor.rel.core[j,i] < 0.075){
      heat.fun.factor.rel.core[j,i] <- 3
    }else if(heat.fun.factor.rel.core[j,i] >= 0.075 && heat.fun.factor.rel.core[j,i] < 0.1){
      heat.fun.factor.rel.core[j,i] <- 4
    }else if(heat.fun.factor.rel.core[j,i] >= 0.1 && heat.fun.factor.rel.core[j,i] < 0.125){
      heat.fun.factor.rel.core[j,i] <- 5
    }else if(heat.fun.factor.rel.core[j,i] >= 0.125){
      heat.fun.factor.rel.core[j,i] <- 6
    }
  }
}
remrow <- c()
for(i in 1:nrow(heat.fun.factor.rel.core)){
  r <- heat.fun.factor.rel.core[i,1]
  if(heat.fun.factor.rel.core[i,2] == r){
    if(heat.fun.factor.rel.core[i,3] == r){
      if(heat.fun.factor.rel.core[i,4] == r){
        if(heat.fun.factor.rel.core[i,5] == r){
          if(heat.fun.factor.rel.core[i,6] == r){
            remrow <- append(remrow, i, length(remrow))
          }
        }
      }
    }
  }
}
heat2.fun.factor.rel.core <- heat.fun.factor.rel.core[-remrow,]


heat.bac.factor.rel.core <- bac.factor.rel.core
for(i in 1:ncol(heat.bac.factor.rel.core)){
  for(j in 1:nrow(heat.bac.factor.rel.core)){
    if(heat.bac.factor.rel.core[j,i] > 0 && heat.bac.factor.rel.core[j,i] < 0.0025){
      heat.bac.factor.rel.core[j,i] <- 1
    }else if(heat.bac.factor.rel.core[j,i] >= 0.0025 && heat.bac.factor.rel.core[j,i] < 0.005){
      heat.bac.factor.rel.core[j,i] <- 2
    }else if(heat.bac.factor.rel.core[j,i] >= 0.005 && heat.bac.factor.rel.core[j,i] < 0.0075){
      heat.bac.factor.rel.core[j,i] <- 3
    }else if(heat.bac.factor.rel.core[j,i] >= 0.0075 && heat.bac.factor.rel.core[j,i] < 0.01){
      heat.bac.factor.rel.core[j,i] <- 4
    }else if(heat.bac.factor.rel.core[j,i] >= 0.01 && heat.bac.factor.rel.core[j,i] < 0.015){
      heat.bac.factor.rel.core[j,i] <- 5
    }else if(heat.bac.factor.rel.core[j,i] >= 0.015){
      heat.bac.factor.rel.core[j,i] <- 6
    }
  }
}
remrow <- c()
for(i in 1:nrow(heat.bac.factor.rel.core)){
  r <- heat.bac.factor.rel.core[i,1]
  if(heat.bac.factor.rel.core[i,2] == r){
    if(heat.bac.factor.rel.core[i,3] == r){
      if(heat.bac.factor.rel.core[i,4] == r){
        if(heat.bac.factor.rel.core[i,5] == r){
          if(heat.bac.factor.rel.core[i,6] == r){
            remrow <- append(remrow, i, length(remrow))
          }
        }
      }
    }
  }
}
heat2.bac.factor.rel.core <- heat.bac.factor.rel.core[-remrow,]

s <- c()
for(i in 1:nrow(heat2.bac.factor.rel.core)){
  if(heat.bac.factor.rel.core[i,1] <= 1){
    if(heat.bac.factor.rel.core[i,2] <= 1){
      if(heat.bac.factor.rel.core[i,3] <= 1){
        if(heat.bac.factor.rel.core[i,4] <= 1){
          if(heat.bac.factor.rel.core[i,5] <= 1){
            if(heat.bac.factor.rel.core[i,6] <= 1){
              s <- append(s, i, length(s))
            }
          }
        }
      }
    }
  }
}
heat3.bac.factor.rel.core <- heat2.bac.factor.rel.core[-s,]


pdf("ROut_publication/Heatmap_core_ASVs.pdf")
heatmap.2(as.matrix(heat2.bac.factor.rel.core), scale = "none", col = c("white", "grey90", "grey80", "grey60", "grey40", "grey20", "black"), trace = "none", RowSideColors = bac.tax.core$farbe[match(rownames(heat2.bac.factor.rel.core), rownames(bac.tax.core))], density.info = "none", key = F, cexRow = 0.2, labRow = "")
heatmap.2(as.matrix(heat2.fun.factor.rel.core), scale = "none", col = c("white", "grey90", "grey80", "grey60", "grey40", "grey20", "black"), trace = "none", RowSideColors = fun.tax.core$farbe[match(rownames(heat2.fun.factor.rel.core), rownames(fun.tax.core))], density.info = "none", key = F, cexRow = 0.4, labRow = paste(fun.tax.core$Phylum, fun.tax.core$Order, fun.tax.core$Genus, sep = "|"), margins = c(6,20))
heatmap.2(as.matrix(heat3.bac.factor.rel.core), scale = "none", col = c("white", "grey90", "grey80", "grey60", "grey40", "grey20", "black"), trace = "none", RowSideColors = bac.tax.core$farbe[match(rownames(heat3.bac.factor.rel.core), rownames(bac.tax.core))], density.info = "none", cexRow = 0.2, labRow = paste(bac.tax.core$Phylum, bac.tax.core$Order, bac.tax.core$Genus, sep = "|"), margins = c(6,20), key = F)
dev.off()

png("ROut_publication/Heatmap_core_ASVs_bac_legend.png")
plot.new()
legend("topleft", legend = c("< 0.25%", "< 0.50%", "< 0.75%", "< 1.00%", "< 1.50%", "> 1.50%"), pch = 15, col = c("grey90", "grey80", "grey60", "grey40", "grey20", "black"), bty = "n", pt.cex = 4.5, cex = 2)
dev.off()

png("ROut_publication/Heatmap_core_ASVs_fun_legend.png")
plot.new()
legend("topleft", legend = c("<   2.5%", "<   5.0%", "<   7.5%", "< 10.0%", "< 12.5%", "> 12.5%"), pch = 15, col = c("grey90", "grey80", "grey60", "grey40", "grey20", "black"), bty = "n", pt.cex = 4.5, cex = 2)
legend("topright", legend = c("Total read abundance"), bty = "n", cex = 2)
dev.off()


#Extracted the most abundant bacteria as there are so many bacteria across all locations.
ab <- c("bacASVs_1", "bacASVs_4", "bacASVs_3", "bacASVs_18", "bacASVs_14", "bacASVs_15", "bacASVs_11", "bacASVs_25", "bacASVs_26", "bacASVs_20", "bacASVs_36", "bacASVs_8", "bacASVs_5", "bacASVs_10", "bacASVs_9", "bacASVs_6", "bacASVs_28", "bacASVs_27", "bacASVs_12", "bacASVs_23", "bacASVs_19", "bacASVs_17", "bacASVs_7", "bacASVs_16", "bacASVs_13", "bacASVs_21", "bacASVs_53", "bacASVs_149")


aggregate(colSums(fun.rel[which(rownames(fun.rel) %in% fungi$fungi[which(fungi$Freq == 3)]),]), list(meta.fun.soil$Location), FUN = mean)
aggregate(colSums(fun.rel[which(rownames(fun.rel) %in% fungi$fungi[which(fungi$Freq == 3)]),]), list(meta.fun.soil$Location), FUN = sd)

kruskal.test(colSums(fun.rel[which(rownames(fun.rel) %in% fungi$fungi[which(fungi$Freq == 3)]),]), meta.fun.soil$Location)


aggregate(colSums(bac.rel[which(rownames(bac.rel) %in% bacti$bacti[which(bacti$Freq == 3)]),]), list(meta.bac.soil$Location), FUN = mean)
aggregate(colSums(bac.rel[which(rownames(bac.rel) %in% bacti$bacti[which(bacti$Freq == 3)]),]), list(meta.bac.soil$Location), FUN = sd)

kruskal.test(colSums(bac.rel[which(rownames(bac.rel) %in% bacti$bacti[which(bacti$Freq == 3)]),]), meta.bac.soil$Location)


#Are there top ASVs that are not core ASVs?
which(rownames(fun.soil.tax.sel) %in% rownames(fun.tax.core)) #this is 6 out of 14
which(rownames(bac.soil.tax.sel) %in% rownames(bac.tax.core)) #this is all

fun.tax.core[which(rownames(fun.soil.tax.sel) %in% rownames(fun.tax.core)),]

#################Snowcover
fun.snow.ls <- vector("list", length = length(levels(as.factor(meta.fun.soil$snowcover))))
bac.snow.ls <- vector("list", length = length(levels(as.factor(meta.bac.soil$snowcover))))
names(fun.snow.ls) <- levels(as.factor(meta.fun.soil$snowcover))
names(bac.snow.ls) <- names(fun.snow.ls)

fungi.snow <- c()
bacti.snow <- c()
for(i in 1:length(fun.snow.ls)){
  rdf.fun <- fun.bi[which(meta.fun.soil$snowcover == levels(as.factor(meta.fun.soil$snowcover))[i]),]
  fun.snow.ls[[i]] <- colnames(rdf.fun)[colSums(rdf.fun) > 0]
  fungi.snow <- append(fungi.snow, fun.snow.ls[[i]])
  fun.snow.ls[[i]] <- data.frame(fun.snow.ls[[i]])
  rdf.bac <- bac.bi[which(meta.bac.soil$snowcover == levels(as.factor(meta.bac.soil$snowcover))[i]),]
  bac.snow.ls[[i]] <- colnames(rdf.bac)[colSums(rdf.bac) > 0]
  bacti.snow <- append(bacti.snow, bac.snow.ls[[i]])
  bac.snow.ls[[i]] <- data.frame(bac.snow.ls[[i]])
}
#write_xlsx(fun.snow.ls, path = "ROut_publication/Venn_fun_snow.xlsx")
#write_xlsx(bac.snow.ls, path = "ROut_publication/Venn_bac_snow.xlsx")


fungi.snow <- data.frame(table(fungi.snow))
bacti.snow <- data.frame(table(bacti.snow))

table(fungi.snow$Freq)
table(bacti.snow$Freq)

snow <- as.character(fungi.snow$fungi.snow[which(fungi.snow$Freq == 1)])
snow.core <- snow[which(snow %in% fungi$fungi[which(fungi$Freq == 3)])]

fun.soil.tax[which(rownames(fun.soil.tax) %in% snow.core),]
fun.factor.rel[which(rownames(fun.factor.rel) %in% snow.core),]
colSums(fun.bi)[which(names(colSums(fun.bi)) %in% snow.core)]



snow.bac <- as.character(bacti.snow$bacti.snow[which(bacti.snow$Freq == 1)])
snow.bac.core <- snow.bac[which(snow.bac %in% bacti$bacti[which(bacti$Freq == 3)])]

bac.soil.tax[which(rownames(bac.soil.tax) %in% snow.bac.core),]
bac.factor.rel[which(rownames(bac.factor.rel) %in% snow.bac.core),]
colSums(fun.bi)[which(names(colSums(fun.bi)) %in% snow.core)]


#############Aldex2 approach - using location separate datasets, because I am afraid that the differences might mask the signal as both effects (location and snow-cover) are small on the entire composition###############################
bac.loc.ls <- list("Kuhtai" = bac.soil[which(meta.bac.soil$Location == "Kuhtai"),which(names(bac.soil) %in% bacti$bacti[which(bacti$Freq == 3)])], "Patscherkofel" = bac.soil[which(meta.bac.soil$Location == "Patscherkofel"),which(names(bac.soil) %in% bacti$bacti[which(bacti$Freq == 3)])], "Praxmar" = bac.soil[which(meta.bac.soil$Location == "Praxmar"),which(names(bac.soil) %in% bacti$bacti[which(bacti$Freq == 3)])])
fun.loc.ls <- list("Kuhtai" = fun.soil[which(meta.fun.soil$Location == "Kuhtai"),which(names(fun.soil) %in% fungi$fungi[which(fungi$Freq == 3)])], "Patscherkofel" = fun.soil[which(meta.fun.soil$Location == "Patscherkofel"),which(names(fun.soil) %in% fungi$fungi[which(fungi$Freq == 3)])], "Praxmar" = fun.soil[which(meta.fun.soil$Location == "Praxmar"),which(names(fun.soil) %in% fungi$fungi[which(fungi$Freq == 3)])])

bac.meta.loc.ls <- list("Kuhtai" = meta.bac.soil[which(meta.bac.soil$Location == "Kuhtai"),], "Patscherkofel" = meta.bac.soil[which(meta.bac.soil$Location == "Patscherkofel"),], "Praxmar" = meta.bac.soil[which(meta.bac.soil$Location == "Praxmar"),])
fun.meta.loc.ls <- list("Kuhtai" = meta.fun.soil[which(meta.fun.soil$Location == "Kuhtai"),], "Patscherkofel" = meta.fun.soil[which(meta.fun.soil$Location == "Patscherkofel"),], "Praxmar" = meta.fun.soil[which(meta.fun.soil$Location == "Praxmar"),])

aldex.bac.ls <- vector("list", length = 3)
names(aldex.bac.ls) <- names(bac.loc.ls)
aldex.fun.ls <- aldex.bac.ls

for(i in 1:length(aldex.bac.ls)){
  rdf.fun <- aldex.clr(t(fun.loc.ls[[i]]), model.matrix(~fun.meta.loc.ls[[i]]$snowcover), mc.samples = 128, verbose = T)
  aldex.fun.ls[[i]] <- aldex.glm(rdf.fun)
  aldex.fun.ls[[i]]$phylum <- fun.soil.tax[match(rownames(aldex.fun.ls[[i]]), rownames(fun.soil.tax)),2]
  aldex.fun.ls[[i]]$class <- fun.soil.tax[match(rownames(aldex.fun.ls[[i]]), rownames(fun.soil.tax)),3]
  aldex.fun.ls[[i]]$family <- fun.soil.tax[match(rownames(aldex.fun.ls[[i]]), rownames(fun.soil.tax)),4]
  aldex.fun.ls[[i]]$genus <- fun.soil.tax[match(rownames(aldex.fun.ls[[i]]), rownames(fun.soil.tax)),5]
  aldex.fun.ls[[i]]$species <- fun.soil.tax[match(rownames(aldex.fun.ls[[i]]), rownames(fun.soil.tax)),6]
  
  #rdf.bac <- aldex.clr(t(bac.loc.ls[[i]]), model.matrix(~bac.meta.loc.ls[[i]]$snowcover), mc.samples = 128, verbose = T)
  #aldex.bac.ls[[i]] <- aldex.glm(rdf.bac)
  #aldex.bac.ls[[i]]$phylum <- bac.soil.tax[match(rownames(aldex.bac.ls[[i]]), rownames(bac.soil.tax)),2]
  #aldex.bac.ls[[i]]$class <- bac.soil.tax[match(rownames(aldex.bac.ls[[i]]), rownames(bac.soil.tax)),3]
  #aldex.bac.ls[[i]]$family <- bac.soil.tax[match(rownames(aldex.bac.ls[[i]]), rownames(bac.soil.tax)),4]
  #aldex.bac.ls[[i]]$genus <- bac.soil.tax[match(rownames(aldex.bac.ls[[i]]), rownames(bac.soil.tax)),5]
  #aldex.bac.ls[[i]]$species <- bac.soil.tax[match(rownames(aldex.bac.ls[[i]]), rownames(bac.soil.tax)),6]
}
write_xlsx(aldex.bac.ls, path = "aldex.bac.ls.xlsx")
write_xlsx(aldex.fun.ls, path = "aldex.fun.ls.xlsx")

#There is nothing significant here. Fascinating.


###Doing Venn analysis for only the core ASVs testing for seasonal differences. 
###Maybe: Do Venn analysis separate for the locations, checking for those ASVs that are repeatedly (across locations) detected in only snow-free or snow-covered soils.

fun.snow.ls <- vector("list", length = length(levels(as.factor(meta.fun.soil$snowcover))))
bac.snow.ls <- vector("list", length = length(levels(as.factor(meta.bac.soil$snowcover))))
names(fun.snow.ls) <- levels(as.factor(meta.fun.soil$snowcover))
names(bac.snow.ls) <- names(fun.snow.ls)

bac.core.bi <- bac.bi[,which(names(bac.bi) %in% bacti$bacti[which(bacti$Freq == 3)])]
bac.core.bi <- bac.core.bi[,colSums(bac.core.bi) > 0]
fun.core.bi <- fun.bi[,which(names(fun.bi) %in% fungi$fungi[which(fungi$Freq == 3)])]
fun.core.bi <- fun.core.bi[,colSums(fun.core.bi) > 0]


fungi.snow <- c()
bacti.snow <- c()
for(i in 1:length(fun.snow.ls)){
  rdf.fun <- fun.core.bi[which(meta.fun.soil$snowcover == levels(as.factor(meta.fun.soil$snowcover))[i]),]
  fun.snow.ls[[i]] <- colnames(rdf.fun)[colSums(rdf.fun) > 0]
  fungi.snow <- append(fungi.snow, fun.snow.ls[[i]])
  fun.snow.ls[[i]] <- data.frame(fun.snow.ls[[i]])
  rdf.bac <- bac.core.bi[which(meta.bac.soil$snowcover == levels(as.factor(meta.bac.soil$snowcover))[i]),]
  bac.snow.ls[[i]] <- colnames(rdf.bac)[colSums(rdf.bac) > 0]
  bacti.snow <- append(bacti.snow, bac.snow.ls[[i]])
  bac.snow.ls[[i]] <- data.frame(bac.snow.ls[[i]])
}
#write_xlsx(fun.snow.ls, path = "ROut_publication/Venn_fun_core_snow.xlsx")
#write_xlsx(bac.snow.ls, path = "ROut_publication/Venn_bac_core_snow.xlsx")


fungi.snow <- data.frame(table(fungi.snow))
bacti.snow <- data.frame(table(bacti.snow))

table(fungi.snow$Freq)
table(bacti.snow$Freq)

fungi.snow$phylum <- fun.soil.tax$Phylum[match(fungi.snow$fungi.snow, rownames(fun.soil.tax))]
fungi.snow$class <- fun.soil.tax$Class[match(fungi.snow$fungi.snow, rownames(fun.soil.tax))]
fungi.snow$order <- fun.soil.tax$Order[match(fungi.snow$fungi.snow, rownames(fun.soil.tax))]
fungi.snow$family <- fun.soil.tax$Family[match(fungi.snow$fungi.snow, rownames(fun.soil.tax))]
fungi.snow$genus <- fun.soil.tax$Genus[match(fungi.snow$fungi.snow, rownames(fun.soil.tax))]
fungi.snow$species <- fun.soil.tax$Species[match(fungi.snow$fungi.snow, rownames(fun.soil.tax))]
fungi.snow$biab <- colSums(fun.bi)[match(fungi.snow$fungi.snow, names(colSums(fun.bi)))]
fungi.snow$season <- "no"
fungi.snow$season[match(fun.snow.ls$no$fun.snow.ls..i.., fungi.snow$fungi.snow)] <- "snowfree"
fungi.snow$season[match(fun.snow.ls$yes$fun.snow.ls..i.., fungi.snow$fungi.snow)] <- "snowcover"
fungi.snow$season[which(fungi.snow$Freq == 2)] <- "both"
fungi.snow$relab <- c(colSums(fun.soil)/sum(colSums(fun.soil)))[match(fungi.snow$fungi.snow, names(fun.soil))]

#write.csv(fungi.snow[which(fungi.snow$Freq == 1),], file = "ROut_publication/Fungi_snow_core.csv")


bacti.snow$phylum <- bac.soil.tax$Phylum[match(bacti.snow$bacti.snow, rownames(bac.soil.tax))]
bacti.snow$class <- bac.soil.tax$Class[match(bacti.snow$bacti.snow, rownames(bac.soil.tax))]
bacti.snow$order <- bac.soil.tax$Order[match(bacti.snow$bacti.snow, rownames(bac.soil.tax))]
bacti.snow$family <- bac.soil.tax$Family[match(bacti.snow$bacti.snow, rownames(bac.soil.tax))]
bacti.snow$genus <- bac.soil.tax$Genus[match(bacti.snow$bacti.snow, rownames(bac.soil.tax))]
bacti.snow$species <- bac.soil.tax$Species[match(bacti.snow$bacti.snow, rownames(bac.soil.tax))]
bacti.snow$biab <- colSums(bac.bi)[match(bacti.snow$bacti.snow, names(colSums(bac.bi)))]
bacti.snow$season <- "no"
bacti.snow$season[match(bac.snow.ls$no$bac.snow.ls..i.., bacti.snow$bacti.snow)] <- "snowfree"
bacti.snow$season[match(bac.snow.ls$yes$bac.snow.ls..i.., bacti.snow$bacti.snow)] <- "snowcover"
bacti.snow$season[which(bacti.snow$Freq == 2)] <- "both"
bacti.snow$relab <- c(colSums(bac.soil)/sum(colSums(bac.soil)))[match(bacti.snow$bacti.snow, names(bac.soil))]

bacti.snow[which(bacti.snow$Freq == 1),]
table(bacti.snow[which(bacti.snow$season == "snowfree"),]$phylum)/117
table(bacti.snow[which(bacti.snow$season == "snowcover"),]$phylum)/67

#write.table(table(bacti.snow[which(bacti.snow$season == "snowfree"),]$phylum), dec = ",", file = "ROut_publication/Bacti_snow_core_snowfree_phylum.csv")
#write.table(table(bacti.snow[which(bacti.snow$season == "snowcover"),]$phylum), dec = ",", file = "ROut_publication/Bacti_snow_core_snowcover_phylum.csv")
#write.csv(bacti.snow[which(bacti.snow$Freq == 1),], file = "ROut_publication/Bacti_snow_core.csv")







#Higher occurrences of certain bacterial phyla depending on snow-cover?
bacti.snow <- bacti.snow[match(names(bac.core.bi), bacti.snow$bacti.snow),]
bacti.snow$bacti.snow <- as.character(bacti.snow$bacti.snow)
bac.core.df <- data.frame(matrix(data = NA, ncol = 5, nrow = length(levels(as.factor(bacti.snow$phylum)))))
names(bac.core.df) <- c("median.snowcover", "mean.snowcover", "median.snowfree", "mean.snowfree", "p.value")
rownames(bac.core.df) <- levels(as.factor(bacti.snow$phylum))
bac.phy.rich.df <- data.frame(matrix(data = NA, ncol = length(levels(as.factor(bacti.snow$phylum))), nrow = nrow(bac.soil)))
rownames(bac.phy.rich.df) <- rownames(bac.soil)
names(bac.phy.rich.df) <- levels(as.factor(bacti.snow$phylum))
for(i in 1:nrow(bac.core.df)){
  if(length(dim(bac.core.bi[, which(bacti.snow$phylum == levels(as.factor(bacti.snow$phylum))[i])])) > 1){
    bac.phy.rich.df[,i] <- rowSums(bac.core.bi[, which(bacti.snow$phylum == levels(as.factor(bacti.snow$phylum))[i])])
  }else{
    bac.phy.rich.df[,i] <- bac.core.bi[, which(bacti.snow$phylum == levels(as.factor(bacti.snow$phylum))[i])]
  }
  bac.core.df$median.snowcover[i] <- median(bac.phy.rich.df[which(meta.bac.soil$snowcover == "yes"),i])
  bac.core.df$mean.snowcover[i] <- mean(bac.phy.rich.df[which(meta.bac.soil$snowcover == "yes"),i])
  
  bac.core.df$median.snowfree[i] <- median(bac.phy.rich.df[which(meta.bac.soil$snowcover == "no"),i])
  bac.core.df$mean.snowfree[i] <- mean(bac.phy.rich.df[which(meta.bac.soil$snowcover == "no"),i])
  
  bac.core.df$p.value[i] <- wilcox.test(bac.phy.rich.df[,i] ~ meta.bac.soil$snowcover)$p.value
}



bac.unicore.bi <- bac.bi[,which(names(bac.bi) %in% bacti.snow$bacti.snow[which(bacti.snow$Freq == 1)])]
bacti.snow.unicore <- bacti.snow[which(bacti.snow$Freq == 1),]
bac.unicore.df <- data.frame(matrix(data = NA, ncol = 5, nrow = length(levels(as.factor(bacti.snow.unicore$phylum)))))
names(bac.unicore.df) <- c("median.snowcover", "mean.snowcover", "median.snowfree", "mean.snowfree", "p.value")
rownames(bac.unicore.df) <- levels(as.factor(bacti.snow.unicore$phylum))

bac.unicore.phy.rich.df <- data.frame(matrix(data = NA, ncol = length(levels(as.factor(bacti.snow.unicore$phylum))), nrow = nrow(bac.soil)))
rownames(bac.unicore.phy.rich.df) <- rownames(bac.soil)
names(bac.unicore.phy.rich.df) <- levels(as.factor(bacti.snow.unicore$phylum))
for(i in 1:nrow(bac.unicore.df)){
  if(length(dim(bac.unicore.bi[, which(bacti.snow.unicore$phylum == levels(as.factor(bacti.snow.unicore$phylum))[i])])) > 1){
    bac.unicore.phy.rich.df[,i] <- rowSums(bac.unicore.bi[, which(bacti.snow.unicore$phylum == levels(as.factor(bacti.snow.unicore$phylum))[i])])
  }else{
    bac.unicore.phy.rich.df[,i] <- bac.unicore.bi[, which(bacti.snow.unicore$phylum == levels(as.factor(bacti.snow.unicore$phylum))[i])]
  }
  bac.unicore.df$median.snowcover[i] <- median(bac.unicore.phy.rich.df[which(meta.bac.soil$snowcover == "yes"),i])
  bac.unicore.df$mean.snowcover[i] <- mean(bac.unicore.phy.rich.df[which(meta.bac.soil$snowcover == "yes"),i])
  
  bac.unicore.df$median.snowfree[i] <- median(bac.unicore.phy.rich.df[which(meta.bac.soil$snowcover == "no"),i])
  bac.unicore.df$mean.snowfree[i] <- mean(bac.unicore.phy.rich.df[which(meta.bac.soil$snowcover == "no"),i])
  
  bac.unicore.df$p.value[i] <- wilcox.test(bac.unicore.phy.rich.df[,i] ~ meta.bac.soil$snowcover)$p.value
}

write.csv(bac.unicore.df, file= "ROut_publication/bac.unicore.df.csv")
write.csv(bac.core.df, file= "ROut_publication/bac.core.df.csv")



fungi.snow <- fungi.snow[match(names(fun.core.bi), fungi.snow$fungi.snow),]
fungi.snow$fungi.snow <- as.character(fungi.snow$fungi.snow)
fun.core.df <- data.frame(matrix(data = NA, ncol = 5, nrow = length(levels(as.factor(fungi.snow$phylum)))))
names(fun.core.df) <- c("median.snowcover", "mean.snowcover", "median.snowfree", "mean.snowfree", "p.value")
rownames(fun.core.df) <- levels(as.factor(fungi.snow$phylum))
fun.phy.rich.df <- data.frame(matrix(data = NA, ncol = length(levels(as.factor(fungi.snow$phylum))), nrow = nrow(fun.soil)))
rownames(fun.phy.rich.df) <- rownames(fun.soil)
names(fun.phy.rich.df) <- levels(as.factor(fungi.snow$phylum))
for(i in 1:nrow(fun.core.df)){
  if(length(dim(fun.core.bi[, which(fungi.snow$phylum == levels(as.factor(fungi.snow$phylum))[i])])) > 1){
    fun.phy.rich.df[,i] <- rowSums(fun.core.bi[, which(fungi.snow$phylum == levels(as.factor(fungi.snow$phylum))[i])])
  }else{
    fun.phy.rich.df[,i] <- fun.core.bi[, which(fungi.snow$phylum == levels(as.factor(fungi.snow$phylum))[i])]
  }
  fun.core.df$median.snowcover[i] <- median(fun.phy.rich.df[which(meta.fun.soil$snowcover == "yes"),i])
  fun.core.df$mean.snowcover[i] <- mean(fun.phy.rich.df[which(meta.fun.soil$snowcover == "yes"),i])
  
  fun.core.df$median.snowfree[i] <- median(fun.phy.rich.df[which(meta.fun.soil$snowcover == "no"),i])
  fun.core.df$mean.snowfree[i] <- mean(fun.phy.rich.df[which(meta.fun.soil$snowcover == "no"),i])
  
  fun.core.df$p.value[i] <- wilcox.test(fun.phy.rich.df[,i] ~ meta.fun.soil$snowcover)$p.value
}



fun.unicore.bi <- fun.bi[,which(names(fun.bi) %in% fungi.snow$fungi.snow[which(fungi.snow$Freq == 1)])]
fungi.snow.unicore <- fungi.snow[which(fungi.snow$Freq == 1),]
fun.unicore.df <- data.frame(matrix(data = NA, ncol = 5, nrow = length(levels(as.factor(fungi.snow.unicore$phylum)))))
names(fun.unicore.df) <- c("median.snowcover", "mean.snowcover", "median.snowfree", "mean.snowfree", "p.value")
rownames(fun.unicore.df) <- levels(as.factor(fungi.snow.unicore$phylum))

fun.unicore.phy.rich.df <- data.frame(matrix(data = NA, ncol = length(levels(as.factor(fungi.snow.unicore$phylum))), nrow = nrow(fun.soil)))
rownames(fun.unicore.phy.rich.df) <- rownames(fun.soil)
names(fun.unicore.phy.rich.df) <- levels(as.factor(fungi.snow.unicore$phylum))
for(i in 1:nrow(fun.unicore.df)){
  if(length(dim(fun.unicore.bi[, which(fungi.snow.unicore$phylum == levels(as.factor(fungi.snow.unicore$phylum))[i])])) > 1){
    fun.unicore.phy.rich.df[,i] <- rowSums(fun.unicore.bi[, which(fungi.snow.unicore$phylum == levels(as.factor(fungi.snow.unicore$phylum))[i])])
  }else{
    fun.unicore.phy.rich.df[,i] <- fun.unicore.bi[, which(fungi.snow.unicore$phylum == levels(as.factor(fungi.snow.unicore$phylum))[i])]
  }
  fun.unicore.df$median.snowcover[i] <- median(fun.unicore.phy.rich.df[which(meta.fun.soil$snowcover == "yes"),i])
  fun.unicore.df$mean.snowcover[i] <- mean(fun.unicore.phy.rich.df[which(meta.fun.soil$snowcover == "yes"),i])
  
  fun.unicore.df$median.snowfree[i] <- median(fun.unicore.phy.rich.df[which(meta.fun.soil$snowcover == "no"),i])
  fun.unicore.df$mean.snowfree[i] <- mean(fun.unicore.phy.rich.df[which(meta.fun.soil$snowcover == "no"),i])
  
  fun.unicore.df$p.value[i] <- wilcox.test(fun.unicore.phy.rich.df[,i] ~ meta.fun.soil$snowcover)$p.value
}

write.csv(fun.unicore.df, file= "ROut_publication/fun.unicore.df.csv")
write.csv(fun.core.df, file= "ROut_publication/fun.core.df.csv")


fun.phy.rich.df$Sample <- rownames(fun.phy.rich.df)
fun.phy.rich.df$Snowcover <- meta.fun.soil$snowcover

png("ROut_publication/Phylumrichness_basidiomycota.png", width = 225)
par(mar = c(3,7,1,1))
boxplot(fun.phy.rich.df$p__Basidiomycota ~ fun.phy.rich.df$Snowcover, col = "maroon", ylab = "ASV richness \n Basidiomycota", xlab = "", ylim = c(0,20), cex.axis = 2, cex.lab = 2, axes = FALSE, frame.plot = T)
Axis(side=2, labels=TRUE, cex.axis = 2)
dev.off()
png("ROut_publication/Phylumrichness_mortierellomycota.png", width = 200)
par(mar = c(3,5,1,1), mgp = c(2,0,0))
boxplot(fun.phy.rich.df$p__Mortierellomycota ~ fun.phy.rich.df$Snowcover, col = "mediumpurple", ylab = "Mortierellomycota", xlab = "", ylim = c(0,20), cex.axis = 2, cex.lab = 2, axes = FALSE, frame.plot = T)
Axis(side=2, labels=FALSE)
dev.off()

bac.phy.rich.df$Sample <- rownames(bac.phy.rich.df)
bac.phy.rich.df$Snowcover <- meta.bac.soil$snowcover

png("ROut_publication/Phylumrichness_acidobacteriota.png", width = 200)
par(mar = c(3,5,1,1))
boxplot(bac.phy.rich.df$Acidobacteriota ~ bac.phy.rich.df$Snowcover, col = "burlywood", ylab = "Acidobacteriota", xlab = "", ylim = c(0,200), cex.axis = 2, cex.lab = 2, axes = FALSE, frame.plot = T)
Axis(side=2, labels=TRUE, cex.axis = 2)
dev.off()
png("ROut_publication/Phylumrichness_bacteroidota.png", width = 200)
par(mar = c(3,5,1,1))
boxplot(bac.phy.rich.df$Bacteroidota ~ bac.phy.rich.df$Snowcover, col = "lightsteelblue", ylab = "Bacteroidota", xlab = "", ylim = c(0,200), cex.axis = 2, cex.lab = 2, axes = FALSE, frame.plot = T)
Axis(side=2, labels=FALSE)
dev.off()
png("ROut_publication/Phylumrichness_myxococcota.png", width = 200)
par(mar = c(3,5,1,1))
boxplot(bac.phy.rich.df$Myxococcota ~ bac.phy.rich.df$Snowcover, col = "darkseagreen", ylab = "Myxococcota", xlab = "", ylim = c(0,200), cex.axis = 2, cex.lab = 2, axes = FALSE, frame.plot = T)
Axis(side=2, labels=FALSE)
dev.off()
png("ROut_publication/Phylumrichness_planctomycetota.png", width = 200)
par(mar = c(3,5,1,1))
boxplot(bac.phy.rich.df$Planctomycetota ~ bac.phy.rich.df$Snowcover, col = "khaki", ylab = "Planctomycetota", xlab = "", ylim = c(0,200), cex.axis = 2, cex.lab = 2, axes = FALSE, frame.plot = T)
Axis(side=2, labels=FALSE)
dev.off()
png("ROut_publication/Phylumrichness_proteobacteria.png", width = 200)
par(mar = c(3,5,1,1))
boxplot(bac.phy.rich.df$Proteobacteria ~ bac.phy.rich.df$Snowcover, col = "steelblue", ylab = "Proteobacteria", xlab = "", ylim = c(0,200), cex.axis = 2, cex.lab = 2, axes = FALSE, frame.plot = T)
Axis(side=2, labels=FALSE)
dev.off()
png("ROut_publication/Phylumrichness_verrucomicrobia.png", width = 200)
par(mar = c(3,5,1,1))
boxplot(bac.phy.rich.df$Verrucomicrobiota ~ bac.phy.rich.df$Snowcover, col = "orange", ylab = "Verrucomicrobia", xlab = "", ylim = c(0,200), cex.axis = 2, cex.lab = 2, axes = FALSE, frame.plot = T)
Axis(side=2, labels=FALSE)
dev.off()

png("ROut_publication/Phylumrichness_pvalues.png")
plot.new()
legend("topright", bty = "n", legend = c("p = 0.031", "p = 0.015", "p = 0.027", "p = 0.050", "p = 0.027", "p = 0.008", "p = 0.075", "p = 0.005"), col = "white", pch = 15, cex = 2)
dev.off()





#####Rhizopogon
fun.soil.tax[which(fun.soil.tax$Genus == "g__Rhizopogon"),]
r <- fun.soil[,which(fun.soil.tax$Genus == "g__Rhizopogon")]
r.bi <- r
r.bi[r.bi > 0] <- 1
names(r.bi) <- c("funASVs_18_Rhizopogon sp.", "funASVs_66_Rhizopogon subbadius", "funASVs_72_Rhizopogon sp.", "funASVs_1012_Rhizopogon odoratus")
pdf("ROut_publication/Rhizopogon_Heatmap.pdf")
heatmap(as.matrix(r.bi), scale = "none", Rowv = NA, Colv = NA, margins = c(15,10), cexCol = 1)
dev.off()


r.factor <- fun.factor[,which(fun.soil.tax$Genus == "g__Rhizopogon")]
r.factor[r.factor > 0] <- 1           
names(r.factor) <- c("funASVs_18_Rhizopogon sp.", "funASVs_66_Rhizopogon subbadius", "funASVs_72_Rhizopogon sp.", "funASVs_1012_Rhizopogon odoratus")
pdf("ROut_publication/Rhizopogon_Heatmap.pdf")
heatmap(as.matrix(r.factor), scale = "none", Rowv = NA, Colv = NA, margins = c(15,10), cexCol = 1)
dev.off()
