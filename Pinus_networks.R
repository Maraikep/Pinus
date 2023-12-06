setwd("C:/Users/marai/Dropbox/UIBK/MICINSNOW/Labstuff/Analysis/Pinus")
##Loading packages
#library(vegan)
#library(multcomp)
#library(compositions)
library(SpiecEasi)
library(igraph)
#Edoardo had the idea of sankey plots. Trying this out.
library(alluvial)
library(ggalluvial)
library(writexl)
library(gridExtra)

##Loading data
load("seqtab_16S.nasrm.RData")
bac <- data.frame(seqtab.nochim.nasrm)
load("seqtab.nochim_ITS.RData")
fun <- data.frame(seqtab.nochim)
load("taxa_16S.nasrm.RData")
bac.tax <- data.frame(taxa.nasrm)
load("taxa_ITS.RData")
fun.tax <- data.frame(taxa_ITS)
rm(taxa.nasrm, taxa_ITS, seqtab.nochim, seqtab.nochim.nasrm)

names(bac) <- paste("bacASVs", 1:ncol(bac), sep = "_")
rownames(bac.tax) <- colnames(bac)
names(fun) <- paste("funASVs", 1:ncol(fun), sep = "_")
rownames(fun.tax) <- names(fun)

meta <- read.delim("metadata.txt", header = T)
meta$Sample[1:84] <- paste(meta$Sample[1:84], meta$snowcover[1:84], sep = "")

##Reducing samples to those needed for the networks.
#Dropping bmb samples.
#Removing those samples with too low sequencing depth.
#If necessary, limiting to those samples present in both, fungal and bacterial, dataset.
bac <- bac[which(rownames(bac) %in% meta$bacids),]
meta.bac <- meta[match(rownames(bac), meta$bacids),]
rownames(bac) <- meta.bac$Sample
bac <- bac[which(meta.bac$snowcover != "bmb"),]
meta.bac <- meta.bac[which(meta.bac$snowcover != "bmb"),]
bac <- bac[-which(meta.bac$Sample %in% c("NEGbmb1", "NEGbmb2", "PRA4.1no", "PRA4.4no")),]
meta.bac <- meta.bac[-which(meta.bac$Sample %in% c("NEGbmb1", "NEGbmb2", "PRA4.1no", "PRA4.4no")),]
meta.bac <- meta.bac[which(rowSums(bac) != 0),]
bac <- bac[which(rowSums(bac) != 0),]
bac.tax <- bac.tax[colSums(bac) > 0,]
bac <- bac[,colSums(bac) > 0]

fun <- fun[which(rownames(fun) %in% meta$funids),]
meta.fun <- meta[match(rownames(fun), meta$funids),]
rownames(fun) <- meta.fun$Sample
fun <- fun[which(meta.fun$snowcover != "bmb"),]
meta.fun <- meta.fun[which(meta.fun$snowcover != "bmb"),]
meta.fun <- meta.fun[which(rowSums(fun) != 0),]
fun <- fun[which(rowSums(fun) != 0),]
fun.tax <- fun.tax[colSums(fun) > 0,]
fun <- fun[,colSums(fun) > 0]
#cbind(rownames(fun), rownames(bac))
#checked that the fungal dataset has the same order as the bacterial one: DONE.
meta <- meta.bac
rm(meta.bac, meta.fun)

##Filtering the dataset
bac.bi <- bac
bac.bi[bac.bi > 0] <- 1
fun.bi <- fun
fun.bi[fun.bi > 0] <- 1

bac.table <- data.frame(table(colSums(bac.bi)))
bac.table$cumsum <- 0
bac.table$mindepth <- 0
bac.table$meddepth <- 0
bac.table$maxdepth <- 0
bac.table$avsize <- 0
bac.table$minsize <- 0
bac.table$maxsize <- 0
for(i in 1:nrow(bac.table)){
  bac.table$cumsum[i] <- sum(bac.table$Freq[i:nrow(bac.table)])
  rdf <- bac[,which(colSums(bac.bi) >= as.integer(bac.table$Var1[i]))]
  bac.table$mindepth[i] <- min(rowSums(rdf))
  bac.table$maxdepth[i] <- max(rowSums(rdf))
  bac.table$meddepth[i] <- median(rowSums(rdf))
  bac.table$avsize[i] <- mean(colSums(rdf)/length(which(rowSums(rdf) > 0)))
  bac.table$minsize[i] <- min(colSums(rdf)/length(which(rowSums(rdf) > 0)))
  bac.table$maxsize[i] <- max(colSums(rdf)/length(which(rowSums(rdf) > 0)))
}
bac.table$percent <- round(bac.table$cumsum/bac.table$cumsum[1]*100, digits = 3)

fun.table <- data.frame(table(colSums(fun.bi)))
fun.table$cumsum <- 0
fun.table$mindepth <- 0
fun.table$meddepth <- 0
fun.table$maxdepth <- 0
fun.table$avsize <- 0
fun.table$minsize <- 0
fun.table$maxsize <- 0
for(i in 1:nrow(fun.table)){
  fun.table$cumsum[i] <- sum(fun.table$Freq[i:nrow(fun.table)])
  rdf <- fun[,which(colSums(fun.bi) >= as.integer(fun.table$Var1[i]))]
  fun.table$mindepth[i] <- min(rowSums(rdf))
  fun.table$maxdepth[i] <- max(rowSums(rdf))
  fun.table$meddepth[i] <- median(rowSums(rdf))
  fun.table$avsize[i] <- mean(colSums(rdf)/length(which(rowSums(rdf) > 0)))
  fun.table$minsize[i] <- min(colSums(rdf)/length(which(rowSums(rdf) > 0)))
  fun.table$maxsize[i] <- max(colSums(rdf)/length(which(rowSums(rdf) > 0)))
}
fun.table$percent <- round(fun.table$cumsum/fun.table$cumsum[1]*100, digits = 3)

write_xlsx(list(fun.table, bac.table), path = "Maraike_ROut/funbac_tables.xlsx")


#Based on these tables and stats I decided to compare the association networks resulting from datasets filtered for the following numbers of samples (...an ASV needs to be present in at least x samples, in order to be considered relevant in interaction network.)
fun.sel <- c(3,5,9,11,19)
bac.sel <- c(5,10,20,31,39,46)



##Coloring approach for bacterial and fungal phyla
farben <- data.frame(phylum = c("p__Ascomycota", "p__Basidiomycota", "p__Mucoromycota", "p__Mortierellomycota", "Acidobacteriota", "Actinobacteriota", "Bacteroidota", "Chloroflexi", "Firmicutes", "Myxococcota", "Planctomycetota", "Proteobacteria", "Verrucomicrobiota", "WPS-2"), farbe = c("pink", "maroon", "thistle", "mediumpurple", "burlywood", "wheat", "lightsteelblue", "darkolivegreen", "lightgoldenrod", "darkseagreen", "khaki", "steelblue", "orange", "darkgoldenrod"))




###################Calculating networks###################################################
##Calculating the networks to individual RObjects
#Saving the networks, winter and summer, in a list. Saving both networks in one list-object; V()s attached.

dimensions <- data.frame(matrix(data = NA, ncol = 14, nrow = 30))
names(dimensions) <- c("funthres", "samples_total_fun", "ASVs_total_fun", "samples_summer_fun", "ASVs_summer_fun", "samples_winter_fun", "ASVs_winter_fun", "bacthres", "samples_total_bac", "ASVs_total_bac", "samples_summer_bac", "ASVs_summer_bac", "samples_winter_bac", "ASVs_winter_bac")
b <- 0
for(funthres in fun.sel){
  fun.ls <- vector("list", length = 2)
  names(fun.ls) <- c("summer", "winter")
  fun.rdf <- fun
  fun.rdf <- fun.rdf[,colSums(fun.bi) >= funthres]
  fun.ls[[1]] <- fun.rdf[which(meta$snowcover == "no"),]
  fun.ls[[1]] <- fun.ls[[1]][rowSums(fun.ls[[1]]) > 0, colSums(fun.ls[[1]]) > 0]
  fun.ls[[2]] <- fun.rdf[which(meta$snowcover == "yes"),]
  fun.ls[[2]] <- fun.ls[[2]][rowSums(fun.ls[[2]]) > 0, colSums(fun.ls[[2]]) > 0]
  for(bacthres in bac.sel){
    b <- b+1
    bac.ls <- vector("list", length = 2)
    names(bac.ls) <- c("summer", "winter")
    bac.rdf <- bac
    bac.rdf <- bac.rdf[,colSums(bac.bi) >= bacthres]
    bac.ls[[1]] <- bac.rdf[which(meta$snowcover == "no"),]
    bac.ls[[1]] <- bac.ls[[1]][rowSums(bac.ls[[1]]) > 0, colSums(bac.ls[[1]]) > 0]
    bac.ls[[2]] <- bac.rdf[which(meta$snowcover == "yes"),]
    bac.ls[[2]] <- bac.ls[[2]][rowSums(bac.ls[[2]]) > 0, colSums(bac.ls[[2]]) > 0]
    dimensions$funthres[b] <- funthres
    dimensions$samples_total_fun[b] <- dim(fun.rdf[rowSums(fun.rdf) > 0,])[1]
    dimensions$ASVs_total_fun[b] <- dim(fun.rdf[rowSums(fun.rdf) > 0,])[2]
    dimensions$samples_summer_fun[b] <- dim(fun.ls[[1]])[1]
    dimensions$ASVs_summer_fun[b] <- dim(fun.ls[[1]])[2]
    dimensions$samples_winter_fun[b] <- dim(fun.ls[[2]])[1]
    dimensions$ASVs_winter_fun[b] <- dim(fun.ls[[2]])[2]
    dimensions$bacthres[b] <- bacthres
    dimensions$samples_total_bac[b] <- dim(bac.rdf[rowSums(bac.rdf) > 0,])[1]
    dimensions$ASVs_total_bac[b] <- dim(bac.rdf[rowSums(bac.rdf) > 0,])[2]
    dimensions$samples_summer_bac[b] <- dim(bac.ls[[1]])[1]
    dimensions$ASVs_summer_bac[b] <- dim(bac.ls[[1]])[2]
    dimensions$samples_winter_bac[b] <- dim(bac.ls[[2]])[1]
    dimensions$ASVs_winter_bac[b] <- dim(bac.ls[[2]])[2]

    # net.ls <- vector("list", length = 2)
    # names(net.ls) <- c("summer", "winter")
    graph.ls <- vector("list", length = 2)
    names(graph.ls) <- c("summer", "winter")
    # 
    # net.ls[[1]] <- multi.spiec.easi(list(as.matrix(bac.ls[[1]]), as.matrix(fun.ls[[1]])), method = "mb", sel.criterion = "stars", verbose = T, pulsar.select = T)
    load(paste("Maraike_ROut/net.ls", b, funthres, bacthres, "RData", sep = "."))
    bm <- symBeta(getOptBeta(net.ls[[1]]), mode="maxabs")
    graph.ls[[1]] <- adj2igraph(bm*getRefit(net.ls[[1]]))
    #graph.ls[[1]] <- adj2igraph(getRefit(net.ls[[1]]))
    V(graph.ls[[1]])$name <- c(names(bac.ls[[1]]), names(fun.ls[[1]]))
    V(graph.ls[[1]])$pch <- c(rep("circle", ncol(bac.ls[[1]])), rep("square", ncol(fun.ls[[1]])))
    V(graph.ls[[1]])$phylum <- "dinosaur"
    V(graph.ls[[1]])$phylum[which(V(graph.ls[[1]])$pch == "circle")] <- bac.tax$Phylum[match(V(graph.ls[[1]])$name[which(V(graph.ls[[1]])$pch == "circle")], rownames(bac.tax))]
    V(graph.ls[[1]])$phylum[which(V(graph.ls[[1]])$pch == "square")] <- fun.tax$Phylum[match(V(graph.ls[[1]])$name[which(V(graph.ls[[1]])$pch == "square")], rownames(fun.tax))]
    V(graph.ls[[1]])$farbe <- "grey"
    V(graph.ls[[1]])$farbe[which(V(graph.ls[[1]])$phylum %in% farben$phylum)] <- farben$farbe[match(V(graph.ls[[1]])$phylum[which(V(graph.ls[[1]])$phylum %in% farben$phylum)], farben$phylum)]
    V(graph.ls[[1]])$abundance <- rep(1, length(V(graph.ls[[1]])$name))
    V(graph.ls[[1]])$abundance[which(V(graph.ls[[1]])$pch == "circle")] <- colSums(bac.ls[[1]])[match(V(graph.ls[[1]])$name[which(V(graph.ls[[1]])$pch == "circle")], names(bac.ls[[1]]))]/sum(colSums(bac[which(meta$snowcover == "no")]))
    V(graph.ls[[1]])$abundance[which(V(graph.ls[[1]])$pch == "square")] <- colSums(fun.ls[[1]])[match(V(graph.ls[[1]])$name[which(V(graph.ls[[1]])$pch == "square")], names(fun.ls[[1]]))]/sum(colSums(fun[which(meta$snowcover == "no")]))

    #net.ls[[2]] <- multi.spiec.easi(list(as.matrix(bac.ls[[2]]), as.matrix(fun.ls[[2]])), method = "mb", sel.criterion = "stars", verbose = T, pulsar.select = T)
    bm <- symBeta(getOptBeta(net.ls[[2]]), mode="maxabs")
    graph.ls[[2]] <- adj2igraph(bm*getRefit(net.ls[[2]]))
    #graph.ls[[2]] <- adj2igraph(getRefit(net.ls[[2]]))
    V(graph.ls[[2]])$name <- c(names(bac.ls[[2]]), names(fun.ls[[2]]))
    V(graph.ls[[2]])$pch <- c(rep("circle", ncol(bac.ls[[2]])), rep("square", ncol(fun.ls[[2]])))
    V(graph.ls[[2]])$phylum <- "dinosaur"
    V(graph.ls[[2]])$phylum[which(V(graph.ls[[2]])$pch == "circle")] <- bac.tax$Phylum[match(V(graph.ls[[2]])$name[which(V(graph.ls[[2]])$pch == "circle")], rownames(bac.tax))]
    V(graph.ls[[2]])$phylum[which(V(graph.ls[[2]])$pch == "square")] <- fun.tax$Phylum[match(V(graph.ls[[2]])$name[which(V(graph.ls[[2]])$pch == "square")], rownames(fun.tax))]
    V(graph.ls[[2]])$farbe <- "grey"
    V(graph.ls[[2]])$farbe[which(V(graph.ls[[2]])$phylum %in% farben$phylum)] <- farben$farbe[match(V(graph.ls[[2]])$phylum[which(V(graph.ls[[2]])$phylum %in% farben$phylum)], farben$phylum)]
    V(graph.ls[[2]])$abundance <- rep(1, length(V(graph.ls[[2]])$name))
    V(graph.ls[[2]])$abundance[which(V(graph.ls[[2]])$pch == "circle")] <- colSums(bac.ls[[2]])[match(V(graph.ls[[2]])$name[which(V(graph.ls[[2]])$pch == "circle")], names(bac.ls[[2]]))]/sum(colSums(bac[which(meta$snowcover == "yes")]))
    V(graph.ls[[2]])$abundance[which(V(graph.ls[[2]])$pch == "square")] <- colSums(fun.ls[[2]])[match(V(graph.ls[[2]])$name[which(V(graph.ls[[2]])$pch == "square")], names(fun.ls[[2]]))]/sum(colSums(fun[which(meta$snowcover == "yes")]))

    #save(net.ls, file = paste("Maraike_ROut/net.ls", b, funthres, bacthres, "RData", sep = "."))
    save(graph.ls, file = paste("Maraike_ROut/graph.ls", b, funthres, bacthres, "RData", sep = "."))
    # png(paste("Maraike_ROut/plot_network", b, "summer", funthres, bacthres, "png", sep = "."))
    # plot(graph.ls[[1]], vertex.size = log10(V(graph.ls[[1]])$abundance*2000000), vertex.label = "", vertex.shape = V(graph.ls[[1]])$pch, vertex.color = V(graph.ls[[1]])$farbe)
    # dev.off()
    # png(paste("Maraike_ROut/plot_network", b, "winter", funthres, bacthres, "png", sep = "."))
    # plot(graph.ls[[2]], vertex.size = log10(V(graph.ls[[2]])$abundance*2000000), vertex.label = "", vertex.shape = V(graph.ls[[2]])$pch, vertex.color = V(graph.ls[[2]])$farbe)
    # dev.off()
  }
}

load("dimensions.rData")
dimensions$sharedbac <- NA
dimensions$sbac <- NA
dimensions$wbac <- NA
dimensions$sharedfun <- NA
dimensions$sfun <- NA
dimensions$wfun <- NA
b <- 1
for(funthres in fun.sel){
  fun.rdf <- fun[,colSums(fun.bi) >= funthres]
  funw.rdf <- fun.rdf
  funs.rdf <- fun.rdf
  funw.rdf <- funw.rdf[which(meta$snowcover == "yes"),]
  funw.rdf <- funw.rdf[rowSums(funw.rdf) > 0, colSums(funw.rdf) > 0]
  funs.rdf <- funs.rdf[which(meta$snowcover == "no"),]
  funs.rdf <- funs.rdf[rowSums(funs.rdf) > 0, colSums(funs.rdf) > 0]
  for(bacthres in bac.sel){
    bac.rdf <- bac[,colSums(bac.bi) >= bacthres]
    bacw.rdf <- bac.rdf
    bacs.rdf <- bac.rdf
    bacw.rdf <- bacw.rdf[which(meta$snowcover == "yes"),]
    bacw.rdf <- bacw.rdf[rowSums(bacw.rdf) > 0, colSums(bacw.rdf) > 0]
    bacs.rdf <- bacs.rdf[which(meta$snowcover == "no"),]
    bacs.rdf <- bacs.rdf[rowSums(bacs.rdf) > 0, colSums(bacs.rdf) > 0]
    
    dimensions$sharedfun[b] <- length(intersect(names(funw.rdf), names(funs.rdf)))
    dimensions$sfun[b] <- length(setdiff(names(funs.rdf), names(funw.rdf)))
    dimensions$wfun[b] <- length(setdiff(names(funw.rdf), names(funs.rdf)))
    
    dimensions$sharedbac[b] <- length(intersect(names(bacw.rdf), names(bacs.rdf)))
    dimensions$sbac[b] <- length(setdiff(names(bacs.rdf), names(bacw.rdf)))
    dimensions$wbac[b] <- length(setdiff(names(bacw.rdf), names(bacs.rdf)))
    b <- b+1
  }
}
dimensions$sharedbac_wperc <- dimensions$sharedbac/dimensions$ASVs_winter_bac*100
dimensions$sharedbac_sperc <- dimensions$sharedbac/dimensions$ASVs_summer_bac*100
dimensions$sharedfun_wperc <- dimensions$sharedfun/dimensions$ASVs_winter_fun*100
dimensions$sharedfun_sperc <- dimensions$sharedfun/dimensions$ASVs_summer_fun*100

#save(dimensions, file = "dimensions.RData")
#write.csv(dimensions, file = "Maraike_ROut/network_dimensions.csv")


summary(dimensions$sharedbac_sperc)
summary(dimensions$sharedbac_wperc)
summary(dimensions$sharedfun_sperc)
summary(dimensions$sharedfun_wperc)

png("Maraike_ROut/Boxplot_shared_ASVs.png", width = 900, height = 150)
par(mfrow = c(1,2), mar = c(5,8,1,1)); boxplot(c(dimensions$sharedfun_sperc, dimensions$sharedfun_wperc) ~ as.factor(c(rep("Free", 30), rep("Covered", 30))), horizontal = T, xlab = "", ylab = "", ylim = c(0,100), las = 1, cex.names = 2, cex.axis = 2, cex.lab = 2); boxplot(c(dimensions$sharedbac_sperc, dimensions$sharedbac_wperc) ~ as.factor(c(rep("Free", 30), rep("Covered", 30))), horizontal = T, xlab = "", ylab = "", ylim = c(0,100), las = 1, cex.names = 2, cex.axis = 2, cex.lab = 2)
dev.off()

png("Maraike_ROut/Boxplot_shared_ASVs_legend.png", width = 900)
par(mar = c(1,1,1,1))
plot.new()
legend("topleft", legend = c("Fungal and bacterial ASVs in both snow free and -covered networks [%]"), bty = "n", col = "white", pch = 15, cex = 2)
dev.off()
####################Limiting the networks to fungal-bacterial associations.#########################
#Checking which associations are present - and which are present across networks, i.e. thresholds.

fun.sel <- c(3,5,9,11,19)
bac.sel <- c(5,10,20,31,39,46)
b <- 0

wintergraphs.ls <- vector("list", length = 30)
summergraphs.ls <- vector("list", length = 30)
winteredges.ls <- vector("list", length = 30)
summeredges.ls <- vector("list", length = 30)
summerwinter_pairs.ls <- vector("list", length = 30)
summerwinterpairs <- c()
phylumpairs.ls <- vector("list", length = 30)
for(funthres in fun.sel){
  for(bacthres in bac.sel){
    b <- b+1
    load(file = paste("Maraike_ROut/graph.ls", b, funthres, bacthres, "RData", sep = "."))
    names(wintergraphs.ls)[b] <- paste(b, funthres, bacthres, sep = "_")
    names(summergraphs.ls)[b] <- paste(b, funthres, bacthres, sep = "_")
    names(winteredges.ls)[b] <- paste(b, funthres, bacthres, sep = "_")
    names(summeredges.ls)[b] <- paste(b, funthres, bacthres, sep = "_")
    names(summerwinter_pairs.ls)[b] <- paste(b, funthres, bacthres, sep = "_")
    names(phylumpairs.ls)[b] <- paste(b, funthres, bacthres, sep = "_")
    
    summer.rdf <- as_edgelist(graph.ls[[1]], names = T)
    funbac.s <- c()
    if(nrow(summer.rdf) > 0){
      for(i in 1:nrow(summer.rdf)){
        if(summer.rdf[i,1] %in% rownames(bac.tax)){
          if(summer.rdf[i,2] %in% rownames(bac.tax)){
            funbac.s <- append(funbac.s, i, length(funbac.s))
          }
        }else{
          if(summer.rdf[i,2] %in% rownames(fun.tax)){
            funbac.s <- append(funbac.s, i, length(funbac.s))
          }
        }
      }
      if(nrow(summer.rdf[-funbac.s,]) > 0){
        summergraphs.ls[[b]] <- graph_from_edgelist(summer.rdf[-funbac.s,], directed = F)
        V(summergraphs.ls[[b]])$farbe <- V(graph.ls[[1]])$farbe[match(V(summergraphs.ls[[b]])$name, V(graph.ls[[1]])$name)]
        V(summergraphs.ls[[b]])$pch <- V(graph.ls[[1]])$pch[match(V(summergraphs.ls[[b]])$name, V(graph.ls[[1]])$name)]
        V(summergraphs.ls[[b]])$pyhlum <- V(graph.ls[[1]])$phylum[match(V(summergraphs.ls[[b]])$name, V(graph.ls[[1]])$name)]
        V(summergraphs.ls[[b]])$abundance <- V(graph.ls[[1]])$abundance[match(V(summergraphs.ls[[b]])$name, V(graph.ls[[1]])$name)]
        png(paste("Maraike_ROut/plot_fb_network", b, "summer", funthres, bacthres, "png", sep = "."))
        plot(summergraphs.ls[[b]], vertex.size = log10(V(summergraphs.ls[[b]])$abundance*2000000), vertex.label = "", vertex.shape = V(summergraphs.ls[[b]])$pch, vertex.color = V(summergraphs.ls[[b]])$farbe)
        dev.off()
        
        wg <- data.frame(as_ids(E(graph.ls[[1]])), E(graph.ls[[1]])$weight)
        names(wg) <- c("pair", "weight")
        #wg$p1 <- t(matrix(unlist(strsplit(wg$pair, "\\|")), nrow = 2))[,1]
        #wg$p2 <- t(matrix(unlist(strsplit(wg$pair, "\\|")), nrow = 2))[,2]
        wg <- wg[-funbac.s,]
        E(summergraphs.ls[[b]])$weight <- wg$weight
        summeredges.ls[[b]] <- data.frame(as_edgelist(summergraphs.ls[[b]], names = T))
        summeredges.ls[[b]]$weight <- wg$weight
        summeredges.ls[[b]]$weight.bi <- summeredges.ls[[b]]$weight
        summeredges.ls[[b]]$weight.bi[which(summeredges.ls[[b]]$weight < 0)] <- -1
        summeredges.ls[[b]]$weight.bi[which(summeredges.ls[[b]]$weight > 0)] <- 1
        summeredges.ls[[b]]$fungus <- "dinosaur"
        summeredges.ls[[b]]$bacterium <- "dinosaur"
        for(i in 1:nrow(summeredges.ls[[b]])){
          if(summeredges.ls[[b]]$X1[i] %in% rownames(bac.tax)){
            summeredges.ls[[b]]$bacterium[i] <- summeredges.ls[[b]]$X1[i]
            summeredges.ls[[b]]$fungus[i] <- summeredges.ls[[b]]$X2[i]
          }else{
            summeredges.ls[[b]]$fungus[i] <- summeredges.ls[[b]]$X1[i]
            summeredges.ls[[b]]$bacterium[i] <- summeredges.ls[[b]]$X2[i]
          }
        }
        summeredges.ls[[b]]$phylum_bac <- bac.tax$Phylum[match(summeredges.ls[[b]]$bacterium, rownames(bac.tax))]
        summeredges.ls[[b]]$phylum_fun <- fun.tax$Phylum[match(summeredges.ls[[b]]$fungus, rownames(fun.tax))]
        summeredges.ls[[b]]$class_bac <- bac.tax$Class[match(summeredges.ls[[b]]$bacterium, rownames(bac.tax))]
        summeredges.ls[[b]]$class_fun <- fun.tax$Class[match(summeredges.ls[[b]]$fungus, rownames(fun.tax))]
        summeredges.ls[[b]]$family_bac <- bac.tax$Family[match(summeredges.ls[[b]]$bacterium, rownames(bac.tax))]
        summeredges.ls[[b]]$family_fun <- fun.tax$Family[match(summeredges.ls[[b]]$fungus, rownames(fun.tax))]
        summeredges.ls[[b]]$genus_bac <- bac.tax$Genus[match(summeredges.ls[[b]]$bacterium, rownames(bac.tax))]
        summeredges.ls[[b]]$genus_fun <- fun.tax$Genus[match(summeredges.ls[[b]]$fungus, rownames(fun.tax))]
        summeredges.ls[[b]]$pair <- paste(summeredges.ls[[b]]$fungus, summeredges.ls[[b]]$bacterium, sep = "_")
        summeredges.ls[[b]]$phylumpair <- paste(summeredges.ls[[b]]$phylum_fun, summeredges.ls[[b]]$phylum_bac, sep = "_")
      }
    }
    
    winter.rdf <- as_edgelist(graph.ls[[2]], names = T)
    funbac.w <- c()
    if(nrow(winter.rdf) > 0){
      for(i in 1:nrow(winter.rdf)){
        if(winter.rdf[i,1] %in% rownames(bac.tax)){
          if(winter.rdf[i,2] %in% rownames(bac.tax)){
            funbac.w <- append(funbac.w, i, length(funbac.w))
          }
        }else{
          if(winter.rdf[i,2] %in% rownames(fun.tax)){
            funbac.w <- append(funbac.w, i, length(funbac.w))
          }
        }
      }
      if(nrow(winter.rdf[-funbac.w,]) > 0){
        wintergraphs.ls[[b]] <- graph_from_edgelist(winter.rdf[-funbac.w,], directed = F)
        V(wintergraphs.ls[[b]])$farbe <- V(graph.ls[[2]])$farbe[match(V(wintergraphs.ls[[b]])$name, V(graph.ls[[2]])$name)]
        V(wintergraphs.ls[[b]])$pch <- V(graph.ls[[2]])$pch[match(V(wintergraphs.ls[[b]])$name, V(graph.ls[[2]])$name)]
        V(wintergraphs.ls[[b]])$pyhlum <- V(graph.ls[[2]])$phylum[match(V(wintergraphs.ls[[b]])$name, V(graph.ls[[2]])$name)]
        V(wintergraphs.ls[[b]])$abundance <- V(graph.ls[[2]])$abundance[match(V(wintergraphs.ls[[b]])$name, V(graph.ls[[2]])$name)]
        png(paste("Maraike_ROut/plot_fb_network", b, "winter", funthres, bacthres, "png", sep = "."))
        plot(wintergraphs.ls[[b]], vertex.size = log10(V(wintergraphs.ls[[b]])$abundance*2000000), vertex.label = "", vertex.shape = V(wintergraphs.ls[[b]])$pch, vertex.color = V(wintergraphs.ls[[b]])$farbe)
        dev.off()
        
        wg <- data.frame(as_ids(E(graph.ls[[2]])), E(graph.ls[[2]])$weight)
        names(wg) <- c("pair", "weight")
        #wg$p1 <- t(matrix(unlist(strsplit(wg$pair, "\\|")), nrow = 2))[,1]
        #wg$p2 <- t(matrix(unlist(strsplit(wg$pair, "\\|")), nrow = 2))[,2]
        wg <- wg[-funbac.w,]
        E(wintergraphs.ls[[b]])$weight <- wg$weight
        
        winteredges.ls[[b]] <- data.frame(as_edgelist(wintergraphs.ls[[b]], names = T))
        winteredges.ls[[b]]$weight <- wg$weight
        winteredges.ls[[b]]$weight.bi <- winteredges.ls[[b]]$weight
        winteredges.ls[[b]]$weight.bi[which(winteredges.ls[[b]]$weight < 0)] <- -1
        winteredges.ls[[b]]$weight.bi[which(winteredges.ls[[b]]$weight > 0)] <- 1
        winteredges.ls[[b]]$fungus <- "dinosaur"
        winteredges.ls[[b]]$bacterium <- "dinosaur"
        for(i in 1:nrow(winteredges.ls[[b]])){
          if(winteredges.ls[[b]]$X1[i] %in% rownames(bac.tax)){
            winteredges.ls[[b]]$bacterium[i] <- winteredges.ls[[b]]$X1[i]
            winteredges.ls[[b]]$fungus[i] <- winteredges.ls[[b]]$X2[i]
          }else{
            winteredges.ls[[b]]$fungus[i] <- winteredges.ls[[b]]$X1[i]
            winteredges.ls[[b]]$bacterium[i] <- winteredges.ls[[b]]$X2[i]
          }
        }
        winteredges.ls[[b]]$phylum_bac <- bac.tax$Phylum[match(winteredges.ls[[b]]$bacterium, rownames(bac.tax))]
        winteredges.ls[[b]]$phylum_fun <- fun.tax$Phylum[match(winteredges.ls[[b]]$fungus, rownames(fun.tax))]
        winteredges.ls[[b]]$class_bac <- bac.tax$Class[match(winteredges.ls[[b]]$bacterium, rownames(bac.tax))]
        winteredges.ls[[b]]$class_fun <- fun.tax$Class[match(winteredges.ls[[b]]$fungus, rownames(fun.tax))]
        winteredges.ls[[b]]$family_bac <- bac.tax$Family[match(winteredges.ls[[b]]$bacterium, rownames(bac.tax))]
        winteredges.ls[[b]]$family_fun <- fun.tax$Family[match(winteredges.ls[[b]]$fungus, rownames(fun.tax))]
        winteredges.ls[[b]]$genus_bac <- bac.tax$Genus[match(winteredges.ls[[b]]$bacterium, rownames(bac.tax))]
        winteredges.ls[[b]]$genus_fun <- fun.tax$Genus[match(winteredges.ls[[b]]$fungus, rownames(fun.tax))]
        winteredges.ls[[b]]$pair <- paste(winteredges.ls[[b]]$fungus, winteredges.ls[[b]]$bacterium, sep = "_")
        winteredges.ls[[b]]$phylumpair <- paste(winteredges.ls[[b]]$phylum_fun, winteredges.ls[[b]]$phylum_bac, sep = "_")
      }
    }
  }
}
save(wintergraphs.ls, file = "wintergraphs.ls.RData")
save(summergraphs.ls, file = "summergraphs.ls.RData")
save(winteredges.ls, file = "winteredges.ls.RData")
save(summeredges.ls, file = "summeredges.ls.RData")

###########################Weights############################################
pdf("Maraike_ROut/Network_weights_summer.pdf")
par(mfrow = c(2,2))
for(i in c(1:length(summeredges.ls))[-c(1,7,13,19,25)]){
  plot(abs(summeredges.ls[[i]]$weight) ~ as.factor(summeredges.ls[[i]]$weight.bi), ylab = names(summeredges.ls)[i], xlab = "", ylim = c(0,1))
}
dev.off()

pdf("Maraike_ROut/Network_weights_winter.pdf")
par(mfrow = c(2,2))
for(i in c(1:length(winteredges.ls))[-c(1,7,13,19,25)]){
  plot(abs(winteredges.ls[[i]]$weight) ~ as.factor(winteredges.ls[[i]]$weight.bi), ylab = names(winteredges.ls)[i], xlab = "", ylim = c(0,1))
}
dev.off()


w.freq.s <- c()
w.freq.w <- c()
for(i in 1:length(summeredges.ls)){
  w.freq.s <- append(w.freq.s, length(which(summeredges.ls[[i]]$weight.bi == 1)), length(w.freq.s))
  w.freq.s <- append(w.freq.s, length(which(summeredges.ls[[i]]$weight.bi == -1)), length(w.freq.s))
  w.freq.w <- append(w.freq.w, length(which(winteredges.ls[[i]]$weight.bi == 1)), length(w.freq.w))
  w.freq.w <- append(w.freq.w, length(which(winteredges.ls[[i]]$weight.bi == -1)), length(w.freq.w))
}
w.freq.s <- t(matrix(w.freq.s, nrow = 2, dimnames = list(Summer = c("positive", "negative"))))
w.freq.w <- t(matrix(w.freq.w, nrow = 2, dimnames = list(Winter = c("positive", "negative"))))

w.freq <- data.frame(w.freq.s, w.freq.w)
names(w.freq) <- c(paste(c(rep("summer", 2), rep("winter", 2)), names(w.freq), sep = "|"))
w.freq <- w.freq[-c(1,7,13,19,25),]
chisq.test(w.freq$`summer|positive`, w.freq$`summer|negative`)
chisq.test(w.freq$`winter|positive`, w.freq$`winter|negative`)
#In both summer and winter, the numbers of positive edges is unrelated to the number of negative edges. 
wilcox.test(w.freq$`winter|positive`, w.freq$`winter|negative`, paired = T)
wilcox.test(w.freq$`summer|positive`, w.freq$`summer|negative`, paired = T)
#In both summer and winter, the number of positive edges is comparable to the number of negative edges.

w.freq.2 <- data.frame(rbind(w.freq.s, w.freq.w))
w.freq.2$season <- as.factor(c(rep("Summer", 30), rep("Winter", 30)))
gmbi <- glm(cbind(positive, negative) ~ season, data = w.freq.2, family = binomial)
summary(aov(gmbi))
summary(gmbi)
step(gmbi)
#The modell with season is the better one. The odds of positive interactions are higher in summer compared to winter. (Although the size effect is not that big.)

gmbi.pred <- predict(gmbi, type = "response", se.fit = T)
aggregate(gmbi$fitted.values, list(w.freq.2$season), mean)

png("Maraike_ROut/Network_binomial_posneg_Associations_season.png", width = 350)
par(mar = c(5,5,1,1))
plot(gmbi$fitted.values ~ w.freq.2$season, ylim = c(0.5,0.7), ylab = "Probability of positive association", xlab = "", cex.lab = 2, cex.axis = 2)
arrows(x0=1, y0=c(gmbi.pred$fit[which(w.freq.2$season == "Summer")]-gmbi.pred$se.fit[which(w.freq.2$season == "Summer")]), x1=1, y1=c(gmbi.pred$fit[which(w.freq.2$season == "Summer")]+gmbi.pred$se.fit[which(w.freq.2$season == "Summer")]), code=3, angle=90, length=0.1)
arrows(x0=2, y0=c(gmbi.pred$fit[which(w.freq.2$season == "Winter")]-gmbi.pred$se.fit[which(w.freq.2$season == "Winter")]), x1=2, y1=c(gmbi.pred$fit[which(w.freq.2$season == "Winter")]+gmbi.pred$se.fit[which(w.freq.2$season == "Winter")]), code=3, angle=90, length=0.1)
dev.off()

sink("Maraike_ROut/Network_binomial_posneg_associations_season.txt")
print(summary(aov(gmbi)))
print(summary(gmbi))
print(step(gmbi))
sink()


##############################################Comparing network stats across thresholds
load(file = "wintergraphs.ls.RData")
load(file = "summergraphs.ls.RData")
load(file = "winteredges.ls.RData")
load(file = "summeredges.ls.RData")

net.stats.df <- data.frame(matrix(data = NA, ncol = 14, nrow = 30))
names(net.stats.df) <- c("degree_summer", "degree_winter", "gsize_summer", "gsize_winter", "gorder_summer", "gorder_winter", "diameter_summer", "diameter_winter", "edge_density_summer", "edge_density_winter", "corabund_summer", "corabund_summer_p", "corabund_winter", "corabund_winter_p")
rownames(net.stats.df) <- names(winteredges.ls)
for(i in 1:nrow(net.stats.df)){
  if(class(summergraphs.ls[[i]]) == "igraph"){
    net.stats.df$degree_summer[i] <- median(degree(summergraphs.ls[[i]]))
    net.stats.df$gsize_summer[i] <- gsize(summergraphs.ls[[i]])
    net.stats.df$gorder_summer[i] <- gorder(summergraphs.ls[[i]])
    net.stats.df$degree_summer[i] <- median(degree(summergraphs.ls[[i]]))
    net.stats.df$diameter_summer[i] <- diameter(summergraphs.ls[[i]], directed = F, weights = NA)
    net.stats.df$edge_density_summer[i] <- edge_density(summergraphs.ls[[i]])
    net.stats.df$corabund_summer[i] <- cor.test(degree(summergraphs.ls[[i]]), V(summergraphs.ls[[i]])$abundance)$estimate
    net.stats.df$corabund_summer_p[i] <- cor.test(degree(summergraphs.ls[[i]]), V(summergraphs.ls[[i]])$abundance)$p.value
  }else{
    net.stats.df$degree_summer[i] <- "no graph"
    net.stats.df$gsize_summer[i] <- "no graph"
    net.stats.df$gorder_summer[i] <- "no graph"
    net.stats.df$degree_summer[i] <- "no graph"
    net.stats.df$diameter_summer[i] <- "no graph"
    net.stats.df$edge_density_summer[i] <- "no graph"
    net.stats.df$corabund_summer[i] <- "no graph"
    net.stats.df$corabund_summer_p[i] <- "no graph"
  }
  if(class(wintergraphs.ls[[i]]) == "igraph"){
    net.stats.df$degree_winter[i] <- median(degree(wintergraphs.ls[[i]]))
    net.stats.df$gsize_winter[i] <- gsize(wintergraphs.ls[[i]])
    net.stats.df$gorder_winter[i] <- gorder(wintergraphs.ls[[i]])
    net.stats.df$degree_winter[i] <- median(degree(wintergraphs.ls[[i]]))
    net.stats.df$diameter_winter[i] <- diameter(wintergraphs.ls[[i]], directed = F, weights = NA)
    net.stats.df$edge_density_winter[i] <- edge_density(wintergraphs.ls[[i]])
    net.stats.df$corabund_winter[i] <- cor.test(degree(wintergraphs.ls[[i]]), V(wintergraphs.ls[[i]])$abundance)$estimate
    net.stats.df$corabund_winter_p[i] <- cor.test(degree(wintergraphs.ls[[i]]), V(wintergraphs.ls[[i]])$abundance)$p.value
  }else{
    net.stats.df$degree_winter[i] <- "no graph"
    net.stats.df$gsize_winter[i] <- "no graph"
    net.stats.df$gorder_winter[i] <- "no graph"
    net.stats.df$degree_winter[i] <- "no graph"
    net.stats.df$diameter_winter[i] <- "no graph"
    net.stats.df$edge_density_winter[i] <- "no graph"
    net.stats.df$corabund_winter[i] <- "no graph"
    net.stats.df$corabund_winter_p[i] <- "no graph"
  }
}
#write.csv(net.stats.df, file = "Maraike_ROut/net.stats_across_thresholds.csv")


###Updating net.stats
##Updating dimensions
#load("dimensions.RData")


b <- 1
for(funthres in fun.sel){
  for(bacthres in bac.sel){
    load(paste("Maraike_ROut/graph.ls", b, funthres, bacthres, "RData", sep = "."))

    net.stats.df$orig_degree_summer[b] <- median(degree(graph.ls[[1]]))
    net.stats.df$orig_degree_winter[b] <- median(degree(graph.ls[[2]]))
    net.stats.df$orig_gsize_summer[b] <- gsize(graph.ls[[1]])
    net.stats.df$orig_gsize_winter[b] <- gsize(graph.ls[[2]])
    net.stats.df$orig_gorder_summer[b] <- gorder(graph.ls[[1]])
    net.stats.df$orig_gorder_winter[b] <- gorder(graph.ls[[2]])
    net.stats.df$orig_diameter_summer[b] <- diameter(graph.ls[[1]], directed = F, weights = NA)
    net.stats.df$orig_diameter_winter[b] <- diameter(graph.ls[[2]], directed = F, weights = NA)
    net.stats.df$orig_density_summer[b] <- edge_density(graph.ls[[1]])
    net.stats.df$orig_density_winter[b] <- edge_density(graph.ls[[2]])

    summer.rdf <- as_edgelist(graph.ls[[1]], names = T)
    sassos <- c(summer.rdf[,1], summer.rdf[,2])
    sassos <- unique(sassos)
    net.stats.df$orig_bacinnet_summer[b] <- length(which(sassos %in% rownames(bac.tax)))
    net.stats.df$orig_funinnet_summer[b] <- length(which(sassos %in% rownames(fun.tax)))
    sffbb1 <- summer.rdf[,1]
    sffbb1[which(sffbb1 %in% rownames(bac.tax))] <- "b"
    sffbb1[which(sffbb1 %in% rownames(fun.tax))] <- "f"
    sffbb2 <- summer.rdf[,2]
    sffbb2[which(sffbb2 %in% rownames(bac.tax))] <- "b"
    sffbb2[which(sffbb2 %in% rownames(fun.tax))] <- "f"
    sffbb <- paste(sffbb1, sffbb2, sep = "")
    
    winter.rdf <- as_edgelist(graph.ls[[2]], names = T)
    wassos <- c(winter.rdf[,1], winter.rdf[,2])
    wassos <- unique(wassos)
    net.stats.df$orig_bacinnet_winter[b] <- length(which(wassos %in% rownames(bac.tax)))
    net.stats.df$orig_funinnet_winter[b] <- length(which(wassos %in% rownames(fun.tax)))
    wffbb1 <- winter.rdf[,1]
    wffbb1[which(wffbb1 %in% rownames(bac.tax))] <- "b"
    wffbb1[which(wffbb1 %in% rownames(fun.tax))] <- "f"
    wffbb2 <- winter.rdf[,2]
    wffbb2[which(wffbb2 %in% rownames(bac.tax))] <- "b"
    wffbb2[which(wffbb2 %in% rownames(fun.tax))] <- "f"
    wffbb <- paste(wffbb1, wffbb2, sep = "")
    
    if(length(sffbb) > 0){
      net.stats.df$bb_summer[b] <- data.frame(table(sffbb))[1,2]
      net.stats.df$bf_summer[b] <- data.frame(table(sffbb))[2,2]
      net.stats.df$ff_summer[b] <- data.frame(table(sffbb))[3,2]
      net.stats.df$bb_winter[b] <- data.frame(table(wffbb))[1,2]
      net.stats.df$bf_winter[b] <- data.frame(table(wffbb))[2,2]
      net.stats.df$ff_winter[b] <- data.frame(table(wffbb))[3,2]
    }
    
    net.stats.df$orig_bacingraph_summer[b] <- length(which(V(graph.ls[[1]])$name %in% rownames(bac.tax)))
    net.stats.df$orig_funingraph_summer[b] <- length(which(V(graph.ls[[1]])$name %in% rownames(fun.tax)))
    net.stats.df$orig_bacingraph_winter[b] <- length(which(V(graph.ls[[2]])$name %in% rownames(bac.tax)))
    net.stats.df$orig_funingraph_winter[b] <- length(which(V(graph.ls[[2]])$name %in% rownames(fun.tax)))
    
    if(class(summeredges.ls[[b]]) == "data.frame"){
      net.stats.df$bacinnet_summer[b] <- length(unique(summeredges.ls[[b]]$bacterium))
      net.stats.df$funinnet_summer[b] <- length(unique(summeredges.ls[[b]]$fungus))
      net.stats.df$bacinnet_winter[b] <- length(unique(winteredges.ls[[b]]$bacterium))
      net.stats.df$funinnet_winter[b] <- length(unique(winteredges.ls[[b]]$fungus))
      net.stats.df$bacingraph_summer[b] <- length(which(V(summergraphs.ls[[b]])$name %in% rownames(bac.tax)))
      net.stats.df$funingraph_summer[b] <- length(which(V(summergraphs.ls[[b]])$name %in% rownames(fun.tax)))
      net.stats.df$bacingraph_winter[b] <- length(which(V(wintergraphs.ls[[b]])$name %in% rownames(bac.tax)))
      net.stats.df$funingraph_winter[b] <- length(which(V(wintergraphs.ls[[b]])$name %in% rownames(fun.tax)))
    }
    
    b <- b+1
  }
}

net.stats.df$orig_fb_density_summer <- net.stats.df$orig_gsize_summer/(net.stats.df$orig_bacingraph_summer*net.stats.df$orig_funingraph_summer)
net.stats.df$orig_fb_density_winter <- net.stats.df$orig_gsize_winter/(net.stats.df$orig_bacingraph_winter*net.stats.df$orig_funingraph_winter)
net.stats.df$fb_density_summer <- as.numeric(net.stats.df$gsize_summer)/(as.numeric(net.stats.df$bacinnet_summer)*as.numeric(net.stats.df$funinnet_summer))
net.stats.df$fb_density_winter <- as.numeric(net.stats.df$gsize_winter)/(as.numeric(net.stats.df$bacinnet_winter)*as.numeric(net.stats.df$funinnet_winter))

net.stats.df$bbper_summer <- net.stats.df$bb_summer/rowSums(net.stats.df[,c(33:35)], na.rm = T)*100
net.stats.df$bfper_summer <- net.stats.df$bf_summer/rowSums(net.stats.df[,c(33:35)], na.rm = T)*100
net.stats.df$ffper_summer <- net.stats.df$ff_summer/rowSums(net.stats.df[,c(33:35)], na.rm = T)*100
#rowSums(net.stats.df[,51:53], na.rm = T)
net.stats.df$bbper_winter <- net.stats.df$bb_winter/rowSums(net.stats.df[,c(36:38)], na.rm = T)*100
net.stats.df$bfper_winter <- net.stats.df$bf_winter/rowSums(net.stats.df[,c(36:38)], na.rm = T)*100
net.stats.df$ffper_winter <- net.stats.df$ff_winter/rowSums(net.stats.df[,c(36:38)], na.rm = T)*100
#rowSums(net.stats.df[,54:56], na.rm = T)

net.stats.df$orig_sharedpairs <- NA
b <- 1
origsharedpairs <- c()
for(funthres in fun.sel){
  for(bacthres in bac.sel){
    load(paste("Maraike_ROut/graph.ls", b, funthres, bacthres, "RData", sep = "."))
    rds.net <- as_edgelist(graph.ls[[1]])
    rdw.net <- as_edgelist(graph.ls[[2]])
    rds.net <- paste(rds.net[,1], rds.net[,2], sep = "_")
    rdw.net <- paste(rdw.net[,1], rdw.net[,2], sep = "_")
    net.stats.df$orig_sharedpairs[b] <- length(intersect(rds.net, rdw.net))
    origsharedpairs <- append(origsharedpairs, intersect(rds.net, rdw.net))
    b <- b+1
  }
}

net.stats.df$fb_sharedpairs <- NA
fbsharedpairs <- c()
for(i in 1:length(summeredges.ls)){
  if(class(summeredges.ls[[i]]) == "data.frame"){
    fbsharedpairs <- append(fbsharedpairs, intersect(summeredges.ls[[i]]$pair, winteredges.ls[[i]]$pair))
    net.stats.df$fb_sharedpairs[i] <- length(intersect(summeredges.ls[[i]]$pair, winteredges.ls[[i]]$pair))
  }
}
fbsharedpairs <- data.frame(table(fbsharedpairs))
table(fbsharedpairs$Freq)

net.stats.df$orig_sharedpairs_perc_s <- net.stats.df$orig_sharedpairs/net.stats.df$orig_gsize_summer*100
net.stats.df$orig_sharedpairs_perc_w <- net.stats.df$orig_sharedpairs/net.stats.df$orig_gsize_winter*100
net.stats.df$fb_sharedpairs_allassos_perc_s <- net.stats.df$fb_sharedpairs/net.stats.df$orig_gsize_summer*100
net.stats.df$fb_sharedpairs_allassos_perc_w <- net.stats.df$fb_sharedpairs/net.stats.df$orig_gsize_winter*100
net.stats.df$fb_sharedpairs_allfb_perc_s <- net.stats.df$fb_sharedpairs/as.numeric(net.stats.df$gsize_summer)*100
net.stats.df$fb_sharedpairs_allfb_perc_w <- net.stats.df$fb_sharedpairs/as.numeric(net.stats.df$gsize_winter)*100
#write.csv(net.stats.df, file = "Maraike_ROut/net.stats.df.2.csv")

summary(net.stats.df$orig_sharedpairs_perc_s)
summary(net.stats.df$orig_sharedpairs_perc_w)
wilcox.test(net.stats.df$orig_sharedpairs_perc_s, net.stats.df$orig_sharedpairs_perc_w, paired = T)

wilcox.test(as.numeric(net.stats.df$degree_summer), as.numeric(net.stats.df$degree_winter), paired = T)#ns p = 0.233
wilcox.test(as.numeric(net.stats.df$gsize_summer), as.numeric(net.stats.df$gsize_winter), paired = T)#ns p = 0.258
wilcox.test(as.numeric(net.stats.df$gorder_summer), as.numeric(net.stats.df$gorder_winter), paired = T)#ns p = 0.628
wilcox.test(as.numeric(net.stats.df$diameter_summer), as.numeric(net.stats.df$diameter_winter), paired = T)#ns p = 0.301
wilcox.test(as.numeric(net.stats.df$edge_density_summer), as.numeric(net.stats.df$edge_density_winter), paired = T)#p=0.0025
mean(as.numeric(net.stats.df$edge_density_summer), na.rm = T) # mean = 0.02336, median = 0.00858
mean(as.numeric(net.stats.df$edge_density_winter), na.rm = T) # mean = 0.02573, median = 0.00861
#these means differ. The difference is (0.00237)
#The networks are not particularly dense, however, the winter network is denser than the summer network. 
#I do not think that we can expect a very dense network... so, this actually makes sense. It is what we would expect.
wilcox.test(net.stats.df$orig_density_summer, net.stats.df$orig_density_winter, paired = T)
summary(c(net.stats.df$orig_density_summer, net.stats.df$orig_density_winter))
wilcox.test(net.stats.df$orig_fb_density_summer, net.stats.df$orig_fb_density_winter, paired = T)
wilcox.test(net.stats.df$fb_density_summer, net.stats.df$fb_density_winter, paired = T)#p = 0.002576; mean difference = 0.0035 (winter is denser than summer; mean summer = 0.0534, winter = 0.0569)

wilcox.test(net.stats.df$orig_gsize_summer, net.stats.df$orig_gsize_winter, paired = T)
wilcox.test(net.stats.df$orig_gorder_summer, net.stats.df$orig_gorder_winter, paired = T)
summary(net.stats.df$orig_gorder_summer)
summary(net.stats.df$orig_gorder_winter)
wilcox.test(net.stats.df$orig_diameter_summer, net.stats.df$orig_gorder_winter, paired = T)
summary(net.stats.df$orig_diameter_summer)
summary(net.stats.df$orig_diameter_winter)

wilcox.test(net.stats.df$bbper_summer, net.stats.df$bbper_winter, paired = T)
wilcox.test(net.stats.df$bfper_summer, net.stats.df$bfper_winter, paired = T)
wilcox.test(net.stats.df$ffper_summer, net.stats.df$ffper_winter, paired = T)
summary(c(net.stats.df$bbper_summer, net.stats.df$bbper_winter))
#sd(c(net.stats.df$bbper_summer, net.stats.df$bbper_winter), na.rm = T)
summary(c(net.stats.df$bfper_summer, net.stats.df$bfper_winter))
#sd(c(net.stats.df$bfper_summer, net.stats.df$bfper_winter), na.rm = T)
summary(c(net.stats.df$ffper_summer, net.stats.df$ffper_winter))
#sd(c(net.stats.df$ffper_summer, net.stats.df$ffper_winter), na.rm = T)


wilcox.test(net.stats.df$bacinnet_summer, net.stats.df$bacinnet_winter, paired = T)
wilcox.test(net.stats.df$funinnet_summer, net.stats.df$funinnet_winter, paired = T)#by trend: p = 0.063; mean diff = -2.44; there are more fungi in winter compared to summer.
#wilcox.test(net.stats.df$funingraph_summer, net.stats.df$funingraph_winter, paired = T)#by trend

wilcox.test(net.stats.df$orig_bacinnet_summer, net.stats.df$orig_bacinnet_winter, paired = T)
wilcox.test(net.stats.df$orig_funinnet_summer, net.stats.df$orig_funinnet_winter, paired = T)
wilcox.test(net.stats.df$orig_bacingraph_summer, net.stats.df$orig_bacingraph_winter, paired = T)
wilcox.test(net.stats.df$orig_funingraph_summer, net.stats.df$orig_funingraph_winter, paired = T)
wilcox.test(net.stats.df$bacingraph_summer, net.stats.df$bacingraph_winter, paired = T)
wilcox.test(net.stats.df$funingraph_summer, net.stats.df$funingraph_winter, paired = T)

wilcox.test(net.stats.df$orig_bacinnet_summer/net.stats.df$orig_bacingraph_summer, net.stats.df$orig_bacinnet_winter/net.stats.df$orig_bacingraph_winter, paired = T) #summer is lower than winter.
wilcox.test(net.stats.df$orig_bacingraph_summer-net.stats.df$orig_bacinnet_summer, net.stats.df$orig_bacingraph_winter-net.stats.df$orig_bacinnet_winter, paired = T)
summary(net.stats.df$orig_bacingraph_summer-net.stats.df$orig_bacinnet_summer)
summary(net.stats.df$orig_bacingraph_winter-net.stats.df$orig_bacinnet_winter)


wilcox.test(net.stats.df$orig_funinnet_summer/net.stats.df$orig_funingraph_summer, net.stats.df$orig_funinnet_winter/net.stats.df$orig_funingraph_winter, paired = T) #summer is lower than winter.
summary(net.stats.df$orig_funingraph_summer-net.stats.df$orig_funinnet_summer)
summary(net.stats.df$orig_funingraph_winter-net.stats.df$orig_funinnet_winter)

wilcox.test(as.numeric(net.stats.df$orig_sharedpairs_perc_s), as.numeric(net.stats.df$orig_sharedpairs_perc_w), paired = T)
wilcox.test(as.numeric(net.stats.df$fb_sharedpairs_allfb_perc_s), as.numeric(net.stats.df$fb_sharedpairs_allfb_perc_w), paired = T)
summary(c(net.stats.df$orig_fb_density_summer, net.stats.df$fb_sharedpairs_allfb_perc_w))



load("dimensions.RData")
nr <- c(2:6,8:12,14:18,20:24,26:30) #using only those networks that resulted with pairs
nr2 <- nr+30
test <- data.frame(Density = c(net.stats.df$orig_density_summer, net.stats.df$orig_density_winter), Season = as.factor(c(rep("Free", 30), rep("Covered", 30))), GSize = c(net.stats.df$orig_gsize_summer, net.stats.df$orig_gsize_winter), GOrder = c(net.stats.df$orig_gorder_summer, net.stats.df$orig_gorder_winter), Shared = c(net.stats.df$orig_sharedpairs_perc_s, net.stats.df$orig_sharedpairs_perc_w), sharedfun = c(dimensions$sharedfun_sperc, dimensions$sharedfun_wperc), sharedbac = c(dimensions$sharedbac_sperc, dimensions$sharedbac_wperc), group = as.factor(c(1:30,1:30)), fbshared = c(net.stats.df$fb_sharedpairs_allfb_perc_s, net.stats.df$fb_sharedpairs_allfb_perc_w))
test <- test[c(nr, nr2),]


png("Maraike_ROut/Orig_network_shared_assos.png", width = 900, height = 350)
print(grid.arrange(ggplot(data = test, aes(y = sharedfun, x = Season)) + geom_boxplot() + geom_point(aes(size = 3)) + geom_line(aes(group=group) , linetype = 2) + theme(legend.position = "none", axis.text = element_text(size = 20), axis.title = element_text(size = 22), panel.background = element_rect(fill = "white", color = "grey"), panel.grid.major = element_line(color = "grey", linetype = "dotted", size = 0.5), panel.grid.minor = element_line(color = "grey", size = 0.5, linetype = "dotted")) + xlab("") + ylab("Shared fungal ASVs [%]") + scale_y_continuous(limits = c(0, 100)), ggplot(data = test, aes(y = sharedbac, x = Season)) + geom_point(aes(size = 3)) + geom_boxplot() + geom_line(aes(group=group) , linetype = 2) + theme(legend.position = "none", axis.text = element_text(size = 20), axis.title = element_text(size = 22), panel.background = element_rect(fill = "white", color = "grey"), panel.grid.major = element_line(color = "grey", linetype = "dotted", size = 0.5), panel.grid.minor = element_line(color = "grey", size = 0.5, linetype = "dotted")) + xlab("") + ylab("Shared bacterial ASVs [%]") + scale_y_continuous(limits = c(0, 100)), ggplot(data = test, aes(y = Shared, x = Season)) + geom_boxplot() + geom_point(aes(size = 3)) + geom_line(aes(group=group) , linetype = 2) + theme(legend.position = "none", axis.text = element_text(size = 20), axis.title = element_text(size = 22), panel.background = element_rect(fill = "white", color = "grey"), panel.grid.major = element_line(color = "grey", linetype = "dotted", size = 0.5), panel.grid.minor = element_line(color = "grey", size = 0.5, linetype = "dotted")) + xlab("") + ylab("Shared associations [%]") + scale_y_continuous(limits = c(0, 25)), ggplot(data = test, aes(y = fbshared, x = Season)) + geom_boxplot() + geom_point(aes(size = 3)) + geom_line(aes(group=group) , linetype = 2) + theme(legend.position = "none", axis.text = element_text(size = 20), axis.title = element_text(size = 22), panel.background = element_rect(fill = "white", color = "grey"), panel.grid.major = element_line(color = "grey", linetype = "dotted", size = 0.5), panel.grid.minor = element_line(color = "grey", size = 0.5, linetype = "dotted")) + xlab("") + ylab("Shared fb associations [%]"), ncol = 4))
dev.off()


png("Maraike_ROut/Orig_network_stats.png", width = 900, height = 350)
print(grid.arrange(ggplot(data = test, aes(y = GSize, x = Season)) + geom_boxplot() + geom_point(aes(size = 3)) + geom_line(aes(group=group) , linetype = 2) + theme(legend.position = "none", axis.text = element_text(size = 20), axis.title = element_text(size = 22), panel.background = element_rect(fill = "white", color = "grey"), panel.grid.major = element_line(color = "grey", linetype = "dotted", size = 0.5), panel.grid.minor = element_line(color = "grey", size = 0.5, linetype = "dotted")) + xlab("") + ylab("#Nodes"), ggplot(data = test, aes(y = GOrder, x = Season)) + geom_point(aes(size = 3)) + geom_boxplot() + geom_line(aes(group=group) , linetype = 2) + theme(legend.position = "none", axis.text = element_text(size = 20), axis.title = element_text(size = 22), panel.background = element_rect(fill = "white", color = "grey"), panel.grid.major = element_line(color = "grey", linetype = "dotted", size = 0.5), panel.grid.minor = element_line(color = "grey", size = 0.5, linetype = "dotted")) + xlab("") + ylab("#Edges"), ggplot(data = test, aes(y = Density, x = Season)) + geom_boxplot() + geom_point(aes(size = 3)) + geom_line(aes(group=group) , linetype = 2) + theme(legend.position = "none", axis.text = element_text(size = 20), axis.title = element_text(size = 22), panel.background = element_rect(fill = "white", color = "grey"), panel.grid.major = element_line(color = "grey", linetype = "dotted", size = 0.5), panel.grid.minor = element_line(color = "grey", size = 0.5, linetype = "dotted")) + xlab("") + ylab("Density"), ncol = 3))
dev.off()


png("Maraike_ROut/Orig_network_shared_assos_2.png", width = 900, height = 350)
print(grid.arrange(ggplot(data = test, aes(y = sharedfun, x = Season)) + geom_boxplot() + geom_point(aes(size = 3)) + geom_line(aes(group=group) , linetype = 2) + theme(legend.position = "none", axis.text = element_text(size = 20), axis.title = element_text(size = 22), panel.background = element_rect(fill = "white", color = "grey"), panel.grid.major = element_line(color = "grey", linetype = "dotted", size = 0.5), panel.grid.minor = element_line(color = "grey", size = 0.5, linetype = "dotted")) + xlab("") + ylab("Shared fungal ASVs [%]") + scale_y_continuous(limits = c(0, 100)), ggplot(data = test, aes(y = sharedbac, x = Season)) + geom_point(aes(size = 3)) + geom_boxplot() + geom_line(aes(group=group) , linetype = 2) + theme(legend.position = "none", axis.text = element_text(size = 20), axis.title = element_text(size = 22), panel.background = element_rect(fill = "white", color = "grey"), panel.grid.major = element_line(color = "grey", linetype = "dotted", size = 0.5), panel.grid.minor = element_line(color = "grey", size = 0.5, linetype = "dotted")) + xlab("") + ylab("Shared bacterial ASVs [%]") + scale_y_continuous(limits = c(0, 100)), ggplot(data = test, aes(y = Shared, x = Season)) + geom_boxplot() + geom_point(aes(size = 3)) + geom_line(aes(group=group) , linetype = 2) + theme(legend.position = "none", axis.text = element_text(size = 20), axis.title = element_text(size = 22), panel.background = element_rect(fill = "white", color = "grey"), panel.grid.major = element_line(color = "grey", linetype = "dotted", size = 0.5), panel.grid.minor = element_line(color = "grey", size = 0.5, linetype = "dotted")) + xlab("") + ylab("Shared associations [%]") + scale_y_continuous(limits = c(0, 25)), ncol = 3))
dev.off()

png("Maraike_ROut/Orig_network_shared_assos_fbs_3.png", width = 300, height = 350)
print(ggplot(data = test, aes(y = fbshared, x = Season)) + geom_boxplot() + geom_point(aes(size = 3)) + geom_line(aes(group=group) , linetype = 2) + theme(legend.position = "none", axis.text = element_text(size = 20), axis.title = element_text(size = 22), panel.background = element_rect(fill = "white", color = "grey"), panel.grid.major = element_line(color = "grey", linetype = "dotted", size = 0.5), panel.grid.minor = element_line(color = "grey", size = 0.5, linetype = "dotted")) + xlab("") + ylab("Shared fb associations [%]"))
dev.off()



##Checking how many joint associations are expected to be shared, if distribution is random

nr <- c(2:6,8:12,14:18,20:24,26:30) #using only those networks that resulted with pairs
random_pairs.ls <- vector("list", length = 25)
for(i in 1:length(random_pairs.ls)){
  random_pairs.ls[[i]] <- data.frame(matrix(data = NA, ncol = 9, nrow = 999))
  names(random_pairs.ls[[i]]) <- c("orig_pairs", "ff_s", "fb_s", "bb_s", "ff_w", "fb_w", "bb_w", "orig_fb_pairs", "fb_pairs")
}
names(random_pairs.ls) <- names(summeredges.ls)[nr]


for(i in 1:length(nr)){
  origbac_s <- paste("bac", 1:net.stats.df$orig_bacinnet_summer[nr[i]], sep = "_")
  origbac_w <- paste("bac", 1:net.stats.df$orig_bacinnet_winter[nr[i]], sep = "_")
  origfun_s <- paste("fun", 1:net.stats.df$orig_funinnet_summer[nr[i]], sep = "_")
  origfun_w <- paste("fun", 1:net.stats.df$orig_funinnet_winter[nr[i]], sep = "_")
  
  origasvs_s <- c(origbac_s, origfun_s)
  origasvs_w <- c(origbac_w, origfun_w)
  
  fbbac_s <- paste("bac", 1:net.stats.df$bacinnet_summer[nr[i]], sep = "_")
  fbfun_s <- paste("fun", 1:net.stats.df$funinnet_summer[nr[i]], sep = "_")
  fbbac_w <- paste("bac", 1:net.stats.df$bacinnet_winter[nr[i]], sep = "_")
  fbfun_w <- paste("fun", 1:net.stats.df$funinnet_winter[nr[i]], sep = "_")
  
  for(n in 1:999){
    orig_s_p1 <- sample(origasvs_s, net.stats.df$orig_gsize_summer[nr[i]], replace = T)
    orig_s_p2 <- sample(origasvs_s, net.stats.df$orig_gsize_summer[nr[i]], replace = T)
    orig_s_p1fb <- substr(orig_s_p1, 1, 1)
    orig_s_p2fb <- substr(orig_s_p2, 1, 1)
    orig_edge_fb_s <- paste(orig_s_p1fb, orig_s_p2fb, sep = "")
    orig_edge_fb_s <- data.frame(table(orig_edge_fb_s))
    
    orig_w_p1 <- sample(origasvs_w, net.stats.df$orig_gsize_winter[nr[i]], replace = T)
    orig_w_p2 <- sample(origasvs_w, net.stats.df$orig_gsize_winter[nr[i]], replace = T)
    orig_w_p1fb <- substr(orig_w_p1, 1, 1)
    orig_w_p2fb <- substr(orig_w_p2, 1, 1)
    orig_edge_fb_w <- paste(orig_w_p1fb, orig_w_p2fb, sep = "")
    orig_edge_fb_w <- data.frame(table(orig_edge_fb_w))
    
    if("ff" %in% orig_edge_fb_s$orig_edge_fb_s){
      random_pairs.ls[[i]]$ff_s[n] <- orig_edge_fb_s$Freq[which(orig_edge_fb_s$orig_edge_fb_s == "ff")]
    }else{
      random_pairs.ls[[i]]$ff_s[n] <- 0
    }
    random_pairs.ls[[i]]$bb_s[n] <- orig_edge_fb_s$Freq[which(orig_edge_fb_s$orig_edge_fb_s == "bb")]
    if("ff" %in% orig_edge_fb_w$orig_edge_fb_w){
      random_pairs.ls[[i]]$ff_w[n] <- orig_edge_fb_w$Freq[which(orig_edge_fb_w$orig_edge_fb_w == "ff")]
    }else{
      random_pairs.ls[[i]]$ff_w[n] <- 0
    }
    random_pairs.ls[[i]]$bb_w[n] <- orig_edge_fb_w$Freq[which(orig_edge_fb_w$orig_edge_fb_w == "bb")]
    
    if("fb" %in% orig_edge_fb_s$orig_edge_fb_s){
      random_pairs.ls[[i]]$fb_s[n] <- orig_edge_fb_s$Freq[which(orig_edge_fb_s$orig_edge_fb_s == "fb")]
    }else{
      random_pairs.ls[[i]]$ff_s[n] <- 0
    }
    if("bf" %in% orig_edge_fb_s$orig_edge_fb_s){
      random_pairs.ls[[i]]$fb_s[n] <- random_pairs.ls[[i]]$fb_s[n] + orig_edge_fb_s$Freq[which(orig_edge_fb_s$orig_edge_fb_s == "bf")]
    }
    
    if("fb" %in% orig_edge_fb_w$orig_edge_fb_w){
      random_pairs.ls[[i]]$fb_w[n] <- orig_edge_fb_w$Freq[which(orig_edge_fb_w$orig_edge_fb_w == "fb")]
    }else{
      random_pairs.ls[[i]]$ff_w[n] <- 0
    }
    if("bf" %in% orig_edge_fb_w$orig_edge_fb_w){
      random_pairs.ls[[i]]$fb_w[n] <- random_pairs.ls[[i]]$fb_w[n] + orig_edge_fb_w$Freq[which(orig_edge_fb_w$orig_edge_fb_w == "bf")]
    }
    
    origfb_s <- data.frame(orig_s_p1, orig_s_p2, orig_s_p1fb, orig_s_p2fb, asso = paste(orig_s_p1fb, orig_s_p2fb, sep = ""))
    origfb_s <- origfb_s[-which(origfb_s$asso %in% c("ff", "bb")),]
    origfb_s$asso[which(origfb_s$asso == "fb")] <- paste(origfb_s$orig_s_p2[which(origfb_s$asso == "fb")], origfb_s$orig_s_p1[which(origfb_s$asso == "fb")], sep = "_")
    origfb_s$asso[which(origfb_s$asso == "bf")] <- paste(origfb_s$orig_s_p1[which(origfb_s$asso == "bf")], origfb_s$orig_s_p2[which(origfb_s$asso == "bf")], sep = "_")
    
    origfb_w <- data.frame(orig_w_p1, orig_w_p2, orig_w_p1fb, orig_w_p2fb, asso = paste(orig_w_p1fb, orig_w_p2fb, sep = ""))
    origfb_w <- origfb_w[-which(origfb_w$asso %in% c("ff", "bb")),]
    origfb_w$asso[which(origfb_w$asso == "fb")] <- paste(origfb_w$orig_w_p2[which(origfb_w$asso == "fb")], origfb_w$orig_w_p1[which(origfb_w$asso == "fb")], sep = "_")
    origfb_w$asso[which(origfb_w$asso == "bf")] <- paste(origfb_w$orig_w_p1[which(origfb_w$asso == "bf")], origfb_w$orig_w_p2[which(origfb_w$asso == "bf")], sep = "_")
    
    random_pairs.ls[[i]]$orig_fb_pairs[n] <- length(intersect(origfb_s$asso, origfb_w$asso))
    
    orig_edge_s <- paste(orig_s_p1, orig_s_p2, sep = "_")
    orig_edge_w <- c(paste(orig_w_p1, orig_w_p2, sep = "_"), paste(orig_w_p2, orig_w_p1, sep = "_"))
    random_pairs.ls[[i]]$orig_pairs[n] <- length(intersect(unique(orig_edge_s), unique(orig_edge_w)))
    
    
    fbedge_s <- paste(sample(fbbac_s, net.stats.df$gsize_summer[nr[i]], replace = T), sample(fbfun_s, net.stats.df$gsize_summer[nr[i]], replace = T), sep = "_")
    fbedge_w <- paste(sample(fbbac_w, net.stats.df$gsize_winter[nr[i]], replace = T), sample(fbfun_w, net.stats.df$gsize_winter[nr[i]], replace = T), sep = "_")
    
    random_pairs.ls[[i]]$fb_pairs[n] <- length(intersect(unique(fbedge_s), unique(fbedge_w)))
  }
}
save(random_pairs.ls, file = "random_pairs.ls.RData")

pdf("Maraike_ROut/Random_network_stats_comparison.pdf")
for(i in 1:length(nr)){
  print(grid.arrange(ggplot(random_pairs.ls[[i]], aes(x=orig_pairs)) + geom_density() + geom_vline(aes(xintercept=net.stats.df$orig_sharedpairs[nr[i]]), color="blue", linetype="dashed", size=1), ggplot(random_pairs.ls[[i]], aes(x=orig_fb_pairs)) + geom_density() + geom_vline(aes(xintercept=net.stats.df$fb_sharedpairs[nr[i]]), color="blue", linetype="dashed", size=1), ggplot(random_pairs.ls[[i]], aes(x=fb_pairs)) + geom_density() + geom_vline(aes(xintercept=net.stats.df$fb_sharedpairs[nr[i]]), color="blue", linetype="dashed", size=1), ggplot(random_pairs.ls[[i]], aes(x=ff_s)) + geom_density() + geom_vline(aes(xintercept=net.stats.df$ff_summer[nr[i]]), color="blue", linetype="dashed", size=1), ggplot(random_pairs.ls[[i]], aes(x=ff_w)) + geom_density() + geom_vline(aes(xintercept=net.stats.df$ff_winter[nr[i]]), color="blue", linetype="dashed", size=1), ggplot(random_pairs.ls[[i]], aes(x=bb_s)) + geom_density() + geom_vline(aes(xintercept=net.stats.df$bb_summer[nr[i]]), color="blue", linetype="dashed", size=1), ggplot(random_pairs.ls[[i]], aes(x=bb_w)) + geom_density() + geom_vline(aes(xintercept=net.stats.df$bb_winter[nr[i]]), color="blue", linetype="dashed", size=1), ggplot(random_pairs.ls[[i]], aes(x=fb_s)) + geom_density() + geom_vline(aes(xintercept=net.stats.df$bf_summer[nr[i]]), color="blue", linetype="dashed", size=1), ggplot(random_pairs.ls[[i]], aes(x=fb_w)) + geom_density() + geom_vline(aes(xintercept=net.stats.df$bf_winter[nr[i]]), color="blue", linetype="dashed", size=1), ncol = 3))
}
dev.off()


load("random_pairs.ls.RData")
value.pos <- data.frame(matrix(data = NA, ncol = ncol(random_pairs.ls[[1]]), nrow = length(random_pairs.ls)))
names(value.pos) <- names(random_pairs.ls[[1]])
rownames(value.pos) <- names(random_pairs.ls)
n <- c(57,35,34,33,38,37,36,58,58)
for(i in 1:length(random_pairs.ls)){
  for(j in 1:ncol(random_pairs.ls[[i]])){
    rv <- random_pairs.ls[[i]][,j]
    rv <- c(rv, net.stats.df[nr[i],n[j]])
    rv <- data.frame(rv, rank(rv))
    value.pos[i,j] <- rv[1000,2]
  }
}
write.csv(value.pos, file = "Maraike_ROut/Ranks_random_value.pos.csv")

x11()
plot(0, xlim = c(1,1000), axes=F, type = "n", xlab = "Ranks", ylab = "")
axis(1)
abline(v = value.pos[,2], col = rainbow(25), lwd = 2)

pdf("Maraike_ROut/Random_net.stats_ranks.pdf")
par(mfrow = c(4,2), mar = c(5,1,2,2))
for(i in c(1,9,4,7,3,6,2,5)){
  plot(0, xlim = c(1,1000), axes=F, type = "n", xlab = "Ranks", ylab = "", cex.lab = 2)
  axis(1, cex.axis = 2)
  abline(v = jitter(value.pos[,i], amount = 10), col = c("slateblue", "slateblue1", "slateblue2", "slateblue3", "slateblue4", "orchid", "orchid1", "orchid2", "orchid3", "orchid4", "orange", "orange1", "orange2", "orange3", "orange4", "orangered", "orangered1", "orangered2", "orangered3", "orangered4", "skyblue", "skyblue1", "skyblue2", "skyblue3", "skyblue4"), lwd = 2)
}
dev.off()
png("Maraike_ROut/Random_net.stats_ranks_legend.png", height = 500)
plot.new()
legend("topleft", legend= c("f3_b10", "f3_b20", "f3_b31", "f3_b39", "f3_b46", "f5_b10", "f5_b20", "f5_b31", "f5_b39", "f5_b46", "f9_b10", "f9_b20", "f9_b31", "f9_b39", "f9_b46", "f11_b10", "f11_b20", "f11_b31", "f11_b39", "f11_b46", "f19_b10", "f19_b20", "f19_b31", "f19_b39", "f19_b46"), pch = 15, pt.cex = 2, col = c("slateblue", "slateblue1", "slateblue2", "slateblue3", "slateblue4", "orchid", "orchid1", "orchid2", "orchid3", "orchid4", "orange", "orange1", "orange2", "orange3", "orange4", "orangered", "orangered1", "orangered2", "orangered3", "orangered4", "skyblue", "skyblue1", "skyblue2", "skyblue3", "skyblue4"), bty = "n")
dev.off()

#Comparing association distribution of pairs across thresholds and between summer and winter on phylum level and ASVs pairs
#Are there ASV-pairs detected across thresholds?
######################################PAIRS############################################
#load(file = "wintergraphs.ls.RData")
#load(file = "summergraphs.ls.RData")
load(file = "winteredges.ls.RData")
load(file = "summeredges.ls.RData")

su_asvpairs <- c()
su_class <- 0
wi_asvpairs <- c()
wi_class <- 0
for(i in 1:length(winteredges.ls)){
  if(class(winteredges.ls[[i]]) == "data.frame"){
    wi_asvpairs <- append(wi_asvpairs, winteredges.ls[[i]]$pair)
    wi_class <- wi_class+1
  }
  if(class(summeredges.ls[[i]]) == "data.frame"){
    su_asvpairs <- append(su_asvpairs, summeredges.ls[[i]]$pair)
    su_class <- su_class+1
  }
}

su_asvpairs <- data.frame(table(su_asvpairs))
wi_asvpairs <- data.frame(table(wi_asvpairs))

png("Maraike_ROut/Histogramm_pairs_across_networks.png", width = 850)
par(mfrow = c(1,2), mar = c(5,5,1,1)); hist(su_asvpairs$Freq, breaks = 25, main = "", xlab = "#Snow free networks", ylab = "#ASV pairs", cex.lab = 2, cex.axis = 2); hist(wi_asvpairs$Freq, breaks = 25, main = "", xlab = "#Snow-covered networks", ylab = "#ASV pairs", cex.lab = 2, cex.axis = 2)
dev.off()

#checking if frequent interactions can be primarily found in those networks with significant correlation between size of ASV and connections.
#--> Clearly, the frequent ones are independent of the correlations. Those networks with insignificant correlations are those that were filtered more strict. In those, the percentage of frequently detected pairs is much higher.
#-->So, I would trust all these networks. Taken together they should provide a very good picture of what is happening. Quite robust. We might miss a lot of associations, but a pattern stays a pattern...


###As I was afraid that the associations found will be location specific rather than pinus specific, I checked, if the associations of the core asvs correspond to frequent associations. This is the case. Above a threshold of 10, there is no difference between frequent associations and the core associations.
#Looking at the distribution of numbers of pairs for thresholds, the threshold to set is either 5 or 10. Checking both.

# #How many of the frequent assos are in which network:
# su_freqassos <- su_asvpairs$su_asvpairs[which(su_asvpairs$Freq > 10)]
# su_freqassos_netmatch <- c()
# su_freqassos_netmatch_per <- c()
# wi_freqassos <- wi_asvpairs$wi_asvpairs[which(wi_asvpairs$Freq > 10)]
# wi_freqassos_netmatch <- c()
# wi_freqassos_netmatch_per <- c()
# for(i in 1:length(summeredges.ls)){
#   if(class(summeredges.ls[[i]]) == "data.frame"){
#     su_freqassos_netmatch <- append(su_freqassos_netmatch, length(intersect(su_freqassos, summeredges.ls[[i]]$pair)), length(su_freqassos_netmatch))
#     su_freqassos_netmatch_per <- append(su_freqassos_netmatch_per, length(intersect(su_freqassos, summeredges.ls[[i]]$pair))/length(summeredges.ls[[i]]$pair), length(su_freqassos_netmatch_per))
#   }
#   if(class(winteredges.ls[[i]]) == "data.frame"){
#     wi_freqassos_netmatch <- append(wi_freqassos_netmatch, length(intersect(wi_freqassos, winteredges.ls[[i]]$pair)), length(wi_freqassos_netmatch))
#     wi_freqassos_netmatch_per <- append(wi_freqassos_netmatch_per, length(intersect(wi_freqassos, winteredges.ls[[i]]$pair))/length(winteredges.ls[[i]]$pair), length(wi_freqassos_netmatch_per))
#   }
# }
# 
# net_pairs <- data.frame(matrix(data = NA, ncol = 2, nrow = 25))
# names(net_pairs) <- c("summer", "winter")
# for(i in 1:nrow(net_pairs)){
#   net_pairs$summer[i] <- length(su_asvpairs$Freq[which(su_asvpairs$Freq >= i)])
#   net_pairs$winter[i] <- length(wi_asvpairs$Freq[which(wi_asvpairs$Freq >= i)])
# }
# 
# write.csv(data.frame(winter=wi_freqassos_netmatch, wi_per=round(wi_freqassos_netmatch_per*100, digits = 0), summer=su_freqassos_netmatch, su_per=round(su_freqassos_netmatch_per*100, digits = 0), net_pairs), file = "Maraike_ROut/networks_assos_thres_freq_bigger11.csv")
# 
# 
##Closer look on pairs
wi_asvpairs$fungus <- matrix(unlist(strsplit(as.character(wi_asvpairs$wi_asvpairs), "_b")), nrow = 2)[1,]
wi_asvpairs$bacterium <- paste("b",matrix(unlist(strsplit(as.character(wi_asvpairs$wi_asvpairs), "_b")), nrow = 2)[2,], sep = "")
wi_asvpairs$phylum_bac <- bac.tax$Phylum[match(wi_asvpairs$bacterium, rownames(bac.tax))]
wi_asvpairs$class_bac <- bac.tax$Class[match(wi_asvpairs$bacterium, rownames(bac.tax))]
wi_asvpairs$order_bac <- bac.tax$Order[match(wi_asvpairs$bacterium, rownames(bac.tax))]
wi_asvpairs$family_bac <- bac.tax$Family[match(wi_asvpairs$bacterium, rownames(bac.tax))]
wi_asvpairs$genus_bac <- bac.tax$Genus[match(wi_asvpairs$bacterium, rownames(bac.tax))]
wi_asvpairs$phylum_fun <- fun.tax$Phylum[match(wi_asvpairs$fungus, rownames(fun.tax))]
wi_asvpairs$class_fun <- fun.tax$Class[match(wi_asvpairs$fungus, rownames(fun.tax))]
wi_asvpairs$order_fun <- fun.tax$Order[match(wi_asvpairs$fungus, rownames(fun.tax))]
wi_asvpairs$family_fun <- fun.tax$Family[match(wi_asvpairs$fungus, rownames(fun.tax))]
wi_asvpairs$genus_fun <- fun.tax$Genus[match(wi_asvpairs$fungus, rownames(fun.tax))]
wi_asvpairs$species_fun <- paste(wi_asvpairs$genus_fun, fun.tax$Species[match(wi_asvpairs$fungus, rownames(fun.tax))], sep = "__")
wi_asvpairs$phylumpair <- paste(wi_asvpairs$phylum_fun, wi_asvpairs$phylum_bac, sep = "_")
wi_asvpairs$speciesphylumpair <- paste(wi_asvpairs$species_fun, wi_asvpairs$phylum_bac, sep = "_")
wi_asvpairs$speciesgenuspair <- paste(wi_asvpairs$species_fun, wi_asvpairs$genus_bac, sep = "_")
wi_asvpairs$genuspair <- paste(wi_asvpairs$genus_fun, wi_asvpairs$genus_bac, sep = "_")

su_asvpairs$fungus <- matrix(unlist(strsplit(as.character(su_asvpairs$su_asvpairs), "_b")), nrow = 2)[1,]
su_asvpairs$bacterium <- paste("b",matrix(unlist(strsplit(as.character(su_asvpairs$su_asvpairs), "_b")), nrow = 2)[2,], sep = "")
su_asvpairs$phylum_bac <- bac.tax$Phylum[match(su_asvpairs$bacterium, rownames(bac.tax))]
su_asvpairs$class_bac <- bac.tax$Class[match(su_asvpairs$bacterium, rownames(bac.tax))]
su_asvpairs$order_bac <- bac.tax$Order[match(su_asvpairs$bacterium, rownames(bac.tax))]
su_asvpairs$family_bac <- bac.tax$Family[match(su_asvpairs$bacterium, rownames(bac.tax))]
su_asvpairs$genus_bac <- bac.tax$Genus[match(su_asvpairs$bacterium, rownames(bac.tax))]
su_asvpairs$phylum_fun <- fun.tax$Phylum[match(su_asvpairs$fungus, rownames(fun.tax))]
su_asvpairs$class_fun <- fun.tax$Class[match(su_asvpairs$fungus, rownames(fun.tax))]
su_asvpairs$order_fun <- fun.tax$Order[match(su_asvpairs$fungus, rownames(fun.tax))]
su_asvpairs$family_fun <- fun.tax$Family[match(su_asvpairs$fungus, rownames(fun.tax))]
su_asvpairs$genus_fun <- fun.tax$Genus[match(su_asvpairs$fungus, rownames(fun.tax))]
su_asvpairs$species_fun <- paste(su_asvpairs$genus_fun, fun.tax$Species[match(su_asvpairs$fungus, rownames(fun.tax))], sep = "__")
su_asvpairs$phylumpair <- paste(su_asvpairs$phylum_fun, su_asvpairs$phylum_bac, sep = "_")
su_asvpairs$speciesphylumpair <- paste(su_asvpairs$species_fun, su_asvpairs$phylum_bac, sep = "_")
su_asvpairs$speciesgenuspair <- paste(su_asvpairs$species_fun, su_asvpairs$genus_bac, sep = "_")
su_asvpairs$genuspair <- paste(su_asvpairs$genus_fun, su_asvpairs$genus_bac, sep = "_")

######Checking for Suillus
#suillus <- rownames(fun.tax)[which(fun.tax$Genus == "g__Suillus")]
# write.csv(wi_asvpairs[which(wi_asvpairs$fungus %in% suillus),], file = "Maraike_ROut/Pairs_Suillus_winter.csv")
# write.csv(su_asvpairs[which(su_asvpairs$fungus %in% suillus),], file = "Maraike_ROut/Pairs_Suillus_summer.csv")

###Extracting those pairs that appear frequently across thresholds
# write.csv(su_asvpairs[which(su_asvpairs$Freq > 5),], file = "Maraike_ROut/su_asvpairs_5.csv")
# write.csv(wi_asvpairs[which(su_asvpairs$Freq > 5),], file = "Maraike_ROut/wi_asvpairs_5.csv")
# write.csv(su_asvpairs[which(su_asvpairs$Freq > 10),], file = "Maraike_ROut/su_asvpairs_10.csv")
# write.csv(wi_asvpairs[which(su_asvpairs$Freq > 10),], file = "Maraike_ROut/wi_asvpairs_10.csv")

# nrow(su_asvpairs[which(su_asvpairs$su_asvpairs %in% intersect(su_asvpairs$su_asvpairs[which(su_asvpairs$Freq > 5)], wi_asvpairs$wi_asvpairs[which(wi_asvpairs$Freq > 5)])),]) #this is 19 shared pairs for a threshold of 5 and 5 shared pairs for a threshold of 10.
# write.csv(su_asvpairs[which(su_asvpairs$su_asvpairs %in% intersect(su_asvpairs$su_asvpairs[which(su_asvpairs$Freq > 5)], wi_asvpairs$wi_asvpairs[which(wi_asvpairs$Freq > 5)])),], file = "Maraike_ROut/shared_pairs_threshold5.csv")


##Referring back to season-specific ASVs mentioned:
wi_asvpairs[intersect(which(wi_asvpairs$fungus %in% c("funASVs_2", "funASVs_17", "funASVs_89", "funASVs_150", "funASVs_168")), which(wi_asvpairs$Freq > 4)),]
write.csv(wi_asvpairs[intersect(which(wi_asvpairs$fungus %in% c("funASVs_2", "funASVs_17", "funASVs_89", "funASVs_150", "funASVs_168")), which(wi_asvpairs$Freq > 4)),], file = "Maraike_ROut/season_specific_ASVs_discussed_refound_network.csv")

sps <- c()
spscore <- c()
for(i in 1:25){
  sps <- append(sps, length(intersect(wi_asvpairs$wi_asvpairs[which(wi_asvpairs$Freq > i)], su_asvpairs$su_asvpairs[which(su_asvpairs$Freq > i)])), length(sps))
  spscore <- append(spscore, length(intersect(wi_asv_core_pairs$wi_asv_core_pairs[which(wi_asv_core_pairs$Freq > i)], su_asv_core_pairs$su_asv_core_pairs[which(su_asv_core_pairs$Freq > i)])), length(spscore))
}


#There is little overlap in pairs. Therefore, it does not make sense to visualize the common and unique pairs. Going with separate visualization for winter and summer.
##On species/genus level: Do we have frequent associations? Are they specific for a season?
genspec_freq.ls <- vector("list", length = 4)
names(genspec_freq.ls) <- c("summer_5", "summer_10", "winter_5", "winter_10")
genspec_freq.ls$summer_5 <- data.frame(table(su_asvpairs$genuspair[which(su_asvpairs$Freq > 5)]))
genspec_freq.ls$winter_5 <- data.frame(table(wi_asvpairs$genuspair[which(wi_asvpairs$Freq > 5)]))
genspec_freq.ls$summer_10 <- data.frame(table(su_asvpairs$genuspair[which(su_asvpairs$Freq > 10)]))
genspec_freq.ls$winter_10 <- data.frame(table(wi_asvpairs$genuspair[which(wi_asvpairs$Freq > 10)]))


farben3 <- c()
farben_genus <- c()
for(i in 1:length(genspec_freq.ls)){
  genspec_freq.ls[[i]]$pair <- genspec_freq.ls[[i]]$Var1
  genspec_freq.ls[[i]]$Var1 <- gsub("g__", "", genspec_freq.ls[[i]]$Var1)
  genspec_freq.ls[[i]]$fun.genus <- matrix(unlist(strsplit(genspec_freq.ls[[i]]$Var1, split = "_")), nrow = 2)[1,]
  genspec_freq.ls[[i]]$bac.genus <- matrix(unlist(strsplit(genspec_freq.ls[[i]]$Var1, split = "_")), nrow = 2)[2,]
  genspec_freq.ls[[i]]$fun.phylum <- fun.tax[!duplicated(fun.tax$Genus),]$Phylum[match(genspec_freq.ls[[i]]$fun.genus, gsub("g__", "", fun.tax[!duplicated(fun.tax$Genus),]$Genus))]
  genspec_freq.ls[[i]]$bac.phylum <- bac.tax[!duplicated(bac.tax$Genus),]$Phylum[match(genspec_freq.ls[[i]]$bac.genus, bac.tax[!duplicated(bac.tax$Genus),]$Genus)]
  genspec_freq.ls[[i]]$fun.genus[which(genspec_freq.ls[[i]]$fun.genus == "NA")] <- "unknown_fungus"
  genspec_freq.ls[[i]]$bac.genus[which(genspec_freq.ls[[i]]$bac.genus == "NA")] <- "unknown_bacterium"
  farben3 <- append(farben3, c(genspec_freq.ls[[i]]$fun.genus, genspec_freq.ls[[i]]$bac.genus), length(farben3))
  farben_genus <- append(farben_genus, c(genspec_freq.ls[[i]]$fun.phylum, genspec_freq.ls[[i]]$bac.phylum), length(farben_genus))
}

farben3 <- data.frame(genus = farben3, phylum = farben_genus)
farben3 <- farben3[!duplicated(farben3$genus),]
farben3$farbe <- "grey"
farben3$farbe[which(farben3$phylum %in% farben$phylum)] <- farben$farbe[match(farben3$phylum[which(farben3$phylum %in% farben$phylum)], farben$phylum)]
farben3 <- farben3[order(farben3$genus),]

farben_genus.ls <- vector("list", length = 4)
names(farben_genus.ls) <- names(genspec_freq.ls)
for(i in 1:length(farben_genus.ls)){
  farben_genus.ls[[i]] <- farben3
  rv <- unique(c(genspec_freq.ls[[i]]$fun.genus, genspec_freq.ls[[i]]$bac.genus))
  farben_genus.ls[[i]] <- farben_genus.ls[[i]][which(farben_genus.ls[[i]]$genus %in% rv),]
  png(paste("Maraike_ROut/sakey_freqassos_genera_", names(farben_genus.ls)[i], ".png", sep = ""))
  print(ggplot(genspec_freq.ls[[i]], aes(axis1=fun.genus, axis2=bac.genus, y=Freq)) + geom_alluvium(aes(fill = fun.genus)) +  geom_stratum(width = 1/3) + theme_void() + geom_stratum(aes(fill = bac.genus)) + geom_text(stat = "stratum", aes(label = after_stat(stratum))) + scale_fill_manual(values = farben_genus.ls[[i]]$farbe, guide = "none"))
  dev.off()
}


##Producing plots which exclude unknown fungi and associations that have been observed only once.
#for nicer figure, editing the labelling.
for(i in 1:length(genspec_freq.ls)){
  genspec_freq.ls[[i]]$bac.genus <- gsub("Burkholderia-Caballeronia-Paraburkholderia", "BCP", genspec_freq.ls[[i]]$bac.genus)
  genspec_freq.ls[[i]]$bac.genus <- gsub("Candidatus", "Cand.", genspec_freq.ls[[i]]$bac.genus)
}
farben3$genus <- gsub("Burkholderia-Caballeronia-Paraburkholderia", "BCP", farben3$genus)
farben3$genus <- gsub("Candidatus", "Cand.", farben3$genus)

farben_genus_sel.ls <- vector("list", length = 4)
names(farben_genus_sel.ls) <- names(genspec_freq.ls)
genspec_freq_sel.ls <- genspec_freq.ls
for(i in 1:length(farben_genus_sel.ls)){
  genspec_freq_sel.ls[[i]] <- genspec_freq_sel.ls[[i]][which(genspec_freq_sel.ls[[i]]$Freq > 1),]
  genspec_freq_sel.ls[[i]] <- genspec_freq_sel.ls[[i]][which(genspec_freq_sel.ls[[i]]$fun.genus != "unknown_fungus"),]
  genspec_freq_sel.ls[[i]] <- genspec_freq_sel.ls[[i]][which(genspec_freq_sel.ls[[i]]$bac.genus != "unknown_bacterium"),]
  
  farben_genus_sel.ls[[i]] <- farben3
  rv <- unique(c(genspec_freq_sel.ls[[i]]$fun.genus, genspec_freq_sel.ls[[i]]$bac.genus))
  farben_genus_sel.ls[[i]] <- farben_genus_sel.ls[[i]][which(farben_genus_sel.ls[[i]]$genus %in% rv),]
  farben_genus_sel.ls[[i]] <- farben_genus_sel.ls[[i]][order(farben_genus_sel.ls[[i]]$genus),]
  pdf(paste("Maraike_ROut/sakey_freqassos_genera_selected_", names(farben_genus_sel.ls)[i], ".pdf", sep = ""))
  print(ggplot(genspec_freq_sel.ls[[i]], aes(axis1=fun.genus, axis2=bac.genus, y=Freq)) + geom_alluvium(aes(fill = fun.genus)) +  geom_stratum(width = 1/3) + theme_void() + geom_stratum(aes(fill = bac.genus)) + geom_text(stat = "stratum", aes(label = after_stat(stratum))) + scale_fill_manual(values = farben_genus_sel.ls[[i]]$farbe, guide = "none"))
  dev.off()
}

write_xlsx(genspec_freq_sel.ls, path = "Maraike_ROut/genspec_freq_sel.ls.xlsx")
write_xlsx(genspec_freq.ls, path = "Maraike_ROut/genspec_freq.ls.xlsx")

length(su_asvpairs$fungus[which(su_asvpairs$Freq > 4)] %in% fungi$fungi[which(fungi$Freq == 3)])
nrow(su_asvpairs[which(su_asvpairs$Freq > 4),])
length(su_asvpairs$bacterium[which(su_asvpairs$Freq > 4)] %in% bacti$bacti[which(bacti$Freq == 3)])

length(wi_asvpairs$fungus[which(wi_asvpairs$Freq > 4)] %in% fungi$fungi[which(fungi$Freq == 3)])
nrow(wi_asvpairs[which(wi_asvpairs$Freq > 4),])
length(wi_asvpairs$bacterium[which(wi_asvpairs$Freq > 4)] %in% bacti$bacti[which(bacti$Freq == 3)])
#These ASVs in the associations > 4 are all core asvs.


#############Network figures for overview. Including all frequent associations > 5

wgraph <- graph_from_edgelist(as.matrix(wi_asvpairs[which(wi_asvpairs$Freq > 5),3:4]), directed = F)
sgraph <- graph_from_edgelist(as.matrix(su_asvpairs[which(su_asvpairs$Freq > 5),3:4]), directed = F)

V(wgraph)$phylum <- "dinosaur"
V(wgraph)$phylum[which(V(wgraph)$name %in% rownames(fun.tax))] <- fun.tax$Phylum[match(V(wgraph)$name[which(V(wgraph)$name %in% rownames(fun.tax))], rownames(fun.tax))]
V(wgraph)$phylum[which(V(wgraph)$name %in% rownames(bac.tax))] <- bac.tax$Phylum[match(V(wgraph)$name[which(V(wgraph)$name %in% rownames(bac.tax))], rownames(bac.tax))]

V(wgraph)$pch <- NA
V(wgraph)$pch[which(V(wgraph)$name %in% rownames(fun.tax))] <- "square"
V(wgraph)$pch[which(V(wgraph)$name %in% rownames(bac.tax))] <- "circle"

V(wgraph)$farbe <- "grey"
V(wgraph)$farbe[which(V(wgraph)$phylum %in% farben$phylum)] <- farben$farbe[match(V(wgraph)$phylum[which(V(wgraph)$phylum %in% farben$phylum)], farben$phylum)]

V(sgraph)$phylum <- "dinosaur"
V(sgraph)$phylum[which(V(sgraph)$name %in% rownames(fun.tax))] <- fun.tax$Phylum[match(V(sgraph)$name[which(V(sgraph)$name %in% rownames(fun.tax))], rownames(fun.tax))]
V(sgraph)$phylum[which(V(sgraph)$name %in% rownames(bac.tax))] <- bac.tax$Phylum[match(V(sgraph)$name[which(V(sgraph)$name %in% rownames(bac.tax))], rownames(bac.tax))]

V(sgraph)$pch <- NA
V(sgraph)$pch[which(V(sgraph)$name %in% rownames(fun.tax))] <- "square"
V(sgraph)$pch[which(V(sgraph)$name %in% rownames(bac.tax))] <- "circle"

V(sgraph)$farbe <- "grey"
V(sgraph)$farbe[which(V(sgraph)$phylum %in% farben$phylum)] <- farben$farbe[match(V(sgraph)$phylum[which(V(sgraph)$phylum %in% farben$phylum)], farben$phylum)]



pdf("Maraike_ROut/network_fig4_winter.pdf")
plot(wgraph, vertex.size = 4, vertex.label = "", vertex.shape = V(wgraph)$pch, vertex.color = V(wgraph)$farbe)
dev.off()
pdf("Maraike_ROut/network_fig4_summer.pdf")
plot(sgraph, vertex.size = 4, vertex.label = "", vertex.shape = V(sgraph)$pch, vertex.color = V(sgraph)$farbe)
dev.off()


pdf("Maraike_ROut/Legend.pdf", paper = "a4r")
par(mar = c(0,0,0,0))
plot.new()
legend("topleft", legend = c(gsub("p__", "", farben$phylum), "Other fungi", "Other bacteria"), col = c(farben$farbe, "grey", "grey"), pch = c(rep(15, 4), rep(16, 10), 15,16), cex = 1, bty = "n", pt.cex = 2, ncol = 1)
dev.off()





##########################Core ASVs and seasonal ASVs inside the networks###############################
load(file = "winteredges.ls.RData")
load(file = "summeredges.ls.RData")
load(file = "wintergraphs.ls.RData")
load(file = "summergraphs.ls.RData")
load(file = "ROut_publication/bacti.RData")
load(file = "ROut_publication/fungi.RData")

bactop <- c("bacASVs_5", "bacASVs_7", "bacASVs_1", "bacASVs_21", "bacASVs_3", "bacASVs_6", "bacASVs_12", "bacASVs_10", "bacASVs_4", "bacASVs_9", "bacASVs_15")
funtop <- c("funASVs_2", "funASVs_1", "funASVs_3", "funASVs_4", "funASVs_5", "funASVs_6", "funASVs_7", "funASVs_8", "funASVs_9", "funASVs_11", "funASVs_12", "funASVs_14", "funASVs_17", "funASVs_18")



##########CORE##########################
# 
# winteredges.core.ls <- winteredges.ls
# wintergraphs.core.ls <- wintergraphs.ls
# summeredges.core.ls <- summeredges.ls
# summergraphs.core.ls <- summergraphs.ls
# for(i in 1:length(winteredges.core.ls)){
#   if(class(winteredges.core.ls[[i]]) == "data.frame"){
#     print(dim(winteredges.core.ls[[i]]))
#     s <- which(winteredges.core.ls[[i]]$bacterium %in% bacti$bacti[which(bacti$Freq == 3)])
#     s <- append(s, which(winteredges.core.ls[[i]]$fungus %in% fungi$fungi[which(fungi$Freq == 3)]))
#     s <- unique(s)
#     winteredges.core.ls[[i]] <- winteredges.core.ls[[i]][s,]
#     print(dim(winteredges.core.ls[[i]]))
#     
#     wintergraphs.core.ls[[i]] <- graph_from_edgelist(as.matrix(winteredges.core.ls[[i]][,5:6]), directed = F)
#     V(wintergraphs.core.ls[[i]])$farbe <- V(wintergraphs.ls[[i]])$farbe[match(V(wintergraphs.core.ls[[i]])$name, V(wintergraphs.ls[[i]])$name)]
#     V(wintergraphs.core.ls[[i]])$pch <- V(wintergraphs.ls[[i]])$pch[match(V(wintergraphs.core.ls[[i]])$name, V(wintergraphs.ls[[i]])$name)]
#     V(wintergraphs.core.ls[[i]])$abundance <- V(wintergraphs.ls[[i]])$abundance[match(V(wintergraphs.core.ls[[i]])$name, V(wintergraphs.ls[[i]])$name)]
#     if(gorder(wintergraphs.core.ls[[i]]) > 0){
#       png(paste("Maraike_ROut/plot_fbcore_network", "winter", names(wintergraphs.core.ls)[i], "png", sep = "."))
#       plot(wintergraphs.core.ls[[i]], vertex.size = log10(V(wintergraphs.core.ls[[i]])$abundance*2000000), vertex.label = "", vertex.shape = V(wintergraphs.core.ls[[i]])$pch, vertex.color = V(wintergraphs.core.ls[[i]])$farbe)
#       dev.off()
#     }
#   }
#   if(class(summeredges.core.ls[[i]]) == "data.frame"){
#     print(dim(summeredges.core.ls[[i]]))
#     s <- which(summeredges.core.ls[[i]]$bacterium %in% bacti$bacti[which(bacti$Freq == 3)])
#     s <- append(s, which(summeredges.core.ls[[i]]$fungus %in% fungi$fungi[which(fungi$Freq == 3)]))
#     s <- unique(s)
#     summeredges.core.ls[[i]] <- summeredges.core.ls[[i]][s,]
#     print(dim(summeredges.core.ls[[i]]))
#     
#     summergraphs.core.ls[[i]] <- graph_from_edgelist(as.matrix(summeredges.core.ls[[i]][,5:6]), directed = F)
#     V(summergraphs.core.ls[[i]])$farbe <- V(summergraphs.ls[[i]])$farbe[match(V(summergraphs.core.ls[[i]])$name, V(summergraphs.ls[[i]])$name)]
#     V(summergraphs.core.ls[[i]])$pch <- V(summergraphs.ls[[i]])$pch[match(V(summergraphs.core.ls[[i]])$name, V(summergraphs.ls[[i]])$name)]
#     V(summergraphs.core.ls[[i]])$abundance <- V(summergraphs.ls[[i]])$abundance[match(V(summergraphs.core.ls[[i]])$name, V(summergraphs.ls[[i]])$name)]
#     if(gorder(summergraphs.core.ls[[i]]) > 0){
#       png(paste("Maraike_ROut/plot_fbcore_network", "summer", names(summergraphs.core.ls)[i], "png", sep = "."))
#       plot(summergraphs.core.ls[[i]], vertex.size = log10(V(summergraphs.core.ls[[i]])$abundance*2000000), vertex.label = "", vertex.shape = V(summergraphs.core.ls[[i]])$pch, vertex.color = V(summergraphs.core.ls[[i]])$farbe)
#       dev.off()
#     }
#   }
# }

# save(summergraphs.core.ls, file = "summergraphs.core.ls.RData")
# save(wintergraphs.core.ls, file = "wintergraphs.core.ls.RData")
# save(summeredges.core.ls, file = "summeredges.core.ls.RData")
# save(winteredges.core.ls, file = "winteredges.core.ls.RData")

#Checking if this distribution is very different from the overall taxonomic distribtion
#load(file = "summergraphs.core.ls.RData")
#load(file = "wintergraphs.core.ls.RData")
load(file = "summeredges.core.ls.RData")
load(file = "winteredges.core.ls.RData")


su_asv_core_pairs <- c()
su_class <- 0
wi_asv_core_pairs <- c()
wi_class <- 0
for(i in 1:length(winteredges.core.ls)){
  if(class(winteredges.core.ls[[i]]) == "data.frame"){
    wi_asv_core_pairs <- append(wi_asv_core_pairs, winteredges.core.ls[[i]]$pair)
    wi_class <- wi_class+1
  }
  if(class(summeredges.core.ls[[i]]) == "data.frame"){
    su_asv_core_pairs <- append(su_asv_core_pairs, summeredges.core.ls[[i]]$pair)
    su_class <- su_class+1
  }
}

su_asv_core_pairs <- data.frame(table(su_asv_core_pairs))
wi_asv_core_pairs <- data.frame(table(wi_asv_core_pairs))

length(which(su_asv_core_pairs$Freq > 5))
length(which(wi_asv_core_pairs$Freq > 5))
length(which(su_asvpairs$Freq > 5))
length(which(wi_asvpairs$Freq > 5))


ov.freq <- data.frame(matrix(data = NA, ncol = 4, nrow = 25))
names(ov.freq) <- c("summer", "winter", "summer_core", "winter_core")
for(i in 1:25){
  ov.freq[i,1] <- length(su_asvpairs$Freq[which(su_asvpairs$Freq > i)])
  ov.freq[i,2] <- length(wi_asvpairs$Freq[which(wi_asvpairs$Freq > i)])
  ov.freq[i,3] <- length(su_asv_core_pairs$Freq[which(su_asv_core_pairs$Freq > i)])
  ov.freq[i,4] <- length(wi_asv_core_pairs$Freq[which(wi_asv_core_pairs$Freq > i)])
}

ov.freq2 <- ov.freq
ov.freq2$summer <- cumsum(ov.freq$summer[25:1])
ov.freq2$summer_core <- cumsum(ov.freq$summer_core[25:1])
ov.freq2$winter <- cumsum(ov.freq$winter[25:1])
ov.freq2$winter_core <- cumsum(ov.freq$winter_core[25:1])
write.csv(ov.freq2, file = "Maraike_ROut/ov.freq.csv")


su_asvpairs$su_asvpairs[which(su_asvpairs$Freq > 10)][-which(su_asvpairs$su_asvpairs[which(su_asvpairs$Freq > 10)] %in% su_asv_core_pairs$su_asv_core_pairs)]
wi_asvpairs$wi_asvpairs[which(wi_asvpairs$Freq > 10)][-which(wi_asvpairs$wi_asvpairs[which(wi_asvpairs$Freq > 10)] %in% wi_asv_core_pairs$wi_asv_core_pairs)]

#At a threshold of 10, there is no difference in pairs comparing core pairs and pairs. At a level of 5, there are three pairs (2 summer and 2 winter, one of them shared) that are not in the core. Perfect. So, all the Sakeys from above can be used as are also for the core.



d <- c()
for(i in 1:length(summeredges.core.ls)){
  if(class(summeredges.core.ls[[i]]) == "data.frame"){
    d <- append(d, i, length(d))
  }
}

#write_xlsx(summeredges.core.ls[d], path = "Maraike_ROut/summeredges.core.ls.xlsx")
#write_xlsx(winteredges.core.ls[d], path = "Maraike_ROut/winteredges.core.ls.xlsx")

# 
# sakeys.summer.ls <- vector("list", length = length(d))
# sakeys.neg.summer.ls <- vector("list", length = length(d))
# sakeys.pos.summer.ls <- vector("list", length = length(d))
# names(sakeys.summer.ls) <- names(summeredges.core.ls)[d]
# names(sakeys.neg.summer.ls) <- names(summeredges.core.ls)[d]
# names(sakeys.pos.summer.ls) <- names(summeredges.core.ls)[d]
# 
# sakeys.winter.ls <- vector("list", length = length(d))
# sakeys.neg.winter.ls <- vector("list", length = length(d))
# sakeys.pos.winter.ls <- vector("list", length = length(d))
# names(sakeys.winter.ls) <- names(winteredges.core.ls)[d]
# names(sakeys.neg.winter.ls) <- names(winteredges.core.ls)[d]
# names(sakeys.pos.winter.ls) <- names(winteredges.core.ls)[d]
# 
# for(i in 1:length(d)){
#   sakeys.summer.ls[[i]] <- data.frame(table(summeredges.core.ls[[d[i]]]$phylumpair))
#   sakeys.neg.summer.ls[[i]] <- data.frame(table(summeredges.core.ls[[d[i]]][which(summeredges.core.ls[[d[i]]]$weight < 0),]$phylumpair))
#   sakeys.pos.summer.ls[[i]] <- data.frame(table(summeredges.core.ls[[d[i]]][which(summeredges.core.ls[[d[i]]]$weight > 0),]$phylumpair))
#   
#     sakeys.summer.ls[[i]]$Var1 <- gsub("p__", "", sakeys.summer.ls[[i]]$Var1)
#   sakeys.summer.ls[[i]]$fungus <- matrix(unlist(strsplit(as.character(sakeys.summer.ls[[i]]$Var1), "_")), nrow = 2)[1,]
#   sakeys.summer.ls[[i]]$bacterium <- matrix(unlist(strsplit(as.character(sakeys.summer.ls[[i]]$Var1), "_")), nrow = 2)[2,]
#   
#   if(nrow(sakeys.neg.summer.ls[[i]]) > 0){
#     sakeys.neg.summer.ls[[i]]$Var1 <- gsub("p__", "", sakeys.neg.summer.ls[[i]]$Var1)
#     sakeys.neg.summer.ls[[i]]$fungus <- matrix(unlist(strsplit(as.character(sakeys.neg.summer.ls[[i]]$Var1), "_")), nrow = 2)[1,]
#     sakeys.neg.summer.ls[[i]]$bacterium <- matrix(unlist(strsplit(as.character(sakeys.neg.summer.ls[[i]]$Var1), "_")), nrow = 2)[2,]
#     
#     if(length(which(sakeys.neg.summer.ls[[i]]$Freq > 1)) > 0){
#       rvfarbe <- unique(as.character(sakeys.neg.summer.ls[[i]][which(sakeys.neg.summer.ls[[i]]$Freq > 1),]$bacterium))
#       rvfarbe <- rvfarbe[order(rvfarbe)]
#       rvfarbe <- farben$farbe[match(rvfarbe, farben$phylum)]
#       png(paste("Maraike_ROut/sakey_summer_neg_", names(sakeys.neg.summer.ls)[i], ".png", sep = ""))
#       print(ggplot(sakeys.neg.summer.ls[[i]][which(sakeys.neg.summer.ls[[i]]$Freq > 1),], aes(axis1=fungus, axis2=bacterium, y=Freq)) + geom_alluvium(aes(fill = bacterium)) +  geom_stratum(width = 1/3) + theme_void() + geom_stratum(aes(fill = bacterium)) + geom_text(stat = "stratum", aes(label = after_stat(stratum))) + theme(legend.position="none") + scale_fill_manual(values = rvfarbe, guide = "none"))
#       dev.off()
#     }
#   }
#   
#   if(nrow(sakeys.pos.summer.ls[[i]]) > 0){
#     sakeys.pos.summer.ls[[i]]$Var1 <- gsub("p__", "", sakeys.pos.summer.ls[[i]]$Var1)
#     sakeys.pos.summer.ls[[i]]$fungus <- matrix(unlist(strsplit(as.character(sakeys.pos.summer.ls[[i]]$Var1), "_")), nrow = 2)[1,]
#     sakeys.pos.summer.ls[[i]]$bacterium <- matrix(unlist(strsplit(as.character(sakeys.pos.summer.ls[[i]]$Var1), "_")), nrow = 2)[2,]
#     
#     if(length(which(sakeys.pos.summer.ls[[i]]$Freq > 1)) > 0){
#       rvfarbe <- unique(as.character(sakeys.pos.summer.ls[[i]][which(sakeys.pos.summer.ls[[i]]$Freq > 1),]$bacterium))
#       rvfarbe <- rvfarbe[order(rvfarbe)]
#       rvfarbe <- farben$farbe[match(rvfarbe, farben$phylum)]
#       png(paste("Maraike_ROut/sakey_summer_pos_", names(sakeys.pos.summer.ls)[i], ".png", sep = ""))
#       print(ggplot(sakeys.pos.summer.ls[[i]][which(sakeys.pos.summer.ls[[i]]$Freq > 1),], aes(axis1=fungus, axis2=bacterium, y=Freq)) + geom_alluvium(aes(fill = bacterium)) +  geom_stratum(width = 1/3) + theme_void() + geom_stratum(aes(fill = bacterium)) + geom_text(stat = "stratum", aes(label = after_stat(stratum))) + theme(legend.position="none") + scale_fill_manual(values = rvfarbe, guide = "none"))
#       dev.off()
#     }
#   }
#   
#   if(length(which(sakeys.summer.ls[[i]]$Freq > 1)) > 0){
#     rvfarbe <- unique(as.character(sakeys.summer.ls[[i]][which(sakeys.summer.ls[[i]]$Freq > 1),]$bacterium))
#     rvfarbe <- rvfarbe[order(rvfarbe)]
#     rvfarbe <- farben$farbe[match(rvfarbe, farben$phylum)]
#     png(paste("Maraike_ROut/sakey_summer_", names(sakeys.summer.ls)[i], ".png", sep = ""))
#     print(ggplot(sakeys.summer.ls[[i]][which(sakeys.summer.ls[[i]]$Freq > 1),], aes(axis1=fungus, axis2=bacterium, y=Freq)) + geom_alluvium(aes(fill = bacterium)) +  geom_stratum(width = 1/3) + theme_void() + geom_stratum(aes(fill = bacterium)) + geom_text(stat = "stratum", aes(label = after_stat(stratum))) + theme(legend.position="none") + scale_fill_manual(values = rvfarbe, guide = "none"))
#     dev.off()
#   }
#   
#   
#   sakeys.winter.ls[[i]] <- data.frame(table(winteredges.core.ls[[d[i]]]$phylumpair))
#   sakeys.neg.winter.ls[[i]] <- data.frame(table(winteredges.core.ls[[d[i]]][which(winteredges.core.ls[[d[i]]]$weight < 0),]$phylumpair))
#   sakeys.pos.winter.ls[[i]] <- data.frame(table(winteredges.core.ls[[d[i]]][which(winteredges.core.ls[[d[i]]]$weight > 0),]$phylumpair))
#   
#     sakeys.winter.ls[[i]]$Var1 <- gsub("p__", "", sakeys.winter.ls[[i]]$Var1)
#   sakeys.winter.ls[[i]]$fungus <- matrix(unlist(strsplit(as.character(sakeys.winter.ls[[i]]$Var1), "_")), nrow = 2)[1,]
#   sakeys.winter.ls[[i]]$bacterium <- matrix(unlist(strsplit(as.character(sakeys.winter.ls[[i]]$Var1), "_")), nrow = 2)[2,]
#   
#   if(nrow(sakeys.neg.winter.ls[[i]]) > 0){
#     sakeys.neg.winter.ls[[i]]$Var1 <- gsub("p__", "", sakeys.neg.winter.ls[[i]]$Var1)
#     sakeys.neg.winter.ls[[i]]$fungus <- matrix(unlist(strsplit(as.character(sakeys.neg.winter.ls[[i]]$Var1), "_")), nrow = 2)[1,]
#     sakeys.neg.winter.ls[[i]]$bacterium <- matrix(unlist(strsplit(as.character(sakeys.neg.winter.ls[[i]]$Var1), "_")), nrow = 2)[2,]
#     
#     if(length(which(sakeys.neg.winter.ls[[i]]$Freq > 1)) > 0){
#       rvfarbe <- unique(as.character(sakeys.neg.winter.ls[[i]][which(sakeys.neg.winter.ls[[i]]$Freq > 1),]$bacterium))
#       rvfarbe <- rvfarbe[order(rvfarbe)]
#       rvfarbe <- farben$farbe[match(rvfarbe, farben$phylum)]
#       png(paste("Maraike_ROut/sakey_winter_neg_", names(sakeys.neg.winter.ls)[i], ".png", sep = ""))
#       print(ggplot(sakeys.neg.winter.ls[[i]][which(sakeys.neg.winter.ls[[i]]$Freq > 1),], aes(axis1=fungus, axis2=bacterium, y=Freq)) + geom_alluvium(aes(fill = bacterium)) +  geom_stratum(width = 1/3) + theme_void() + geom_stratum(aes(fill = bacterium)) + geom_text(stat = "stratum", aes(label = after_stat(stratum))) + theme(legend.position="none") + scale_fill_manual(values = rvfarbe, guide = "none"))
#       dev.off()
#     }
#   }
#   
#   if(nrow(sakeys.pos.winter.ls[[i]]) > 0){
#     sakeys.pos.winter.ls[[i]]$Var1 <- gsub("p__", "", sakeys.pos.winter.ls[[i]]$Var1)
#     sakeys.pos.winter.ls[[i]]$fungus <- matrix(unlist(strsplit(as.character(sakeys.pos.winter.ls[[i]]$Var1), "_")), nrow = 2)[1,]
#     sakeys.pos.winter.ls[[i]]$bacterium <- matrix(unlist(strsplit(as.character(sakeys.pos.winter.ls[[i]]$Var1), "_")), nrow = 2)[2,]
#     
#     if(length(which(sakeys.pos.winter.ls[[i]]$Freq > 1)) > 0){
#       rvfarbe <- unique(as.character(sakeys.pos.winter.ls[[i]][which(sakeys.pos.winter.ls[[i]]$Freq > 1),]$bacterium))
#       rvfarbe <- rvfarbe[order(rvfarbe)]
#       rvfarbe <- farben$farbe[match(rvfarbe, farben$phylum)]
#       png(paste("Maraike_ROut/sakey_winter_pos_", names(sakeys.pos.winter.ls)[i], ".png", sep = ""))
#       print(ggplot(sakeys.pos.winter.ls[[i]][which(sakeys.pos.winter.ls[[i]]$Freq > 1),], aes(axis1=fungus, axis2=bacterium, y=Freq)) + geom_alluvium(aes(fill = bacterium)) +  geom_stratum(width = 1/3) + theme_void() + geom_stratum(aes(fill = bacterium)) + geom_text(stat = "stratum", aes(label = after_stat(stratum))) + theme(legend.position="none") + scale_fill_manual(values = rvfarbe, guide = "none"))
#       dev.off()
#     }
#   }
#   
#   if(length(which(sakeys.winter.ls[[i]]$Freq > 1)) > 0){
#     rvfarbe <- unique(as.character(sakeys.winter.ls[[i]][which(sakeys.winter.ls[[i]]$Freq > 1),]$bacterium))
#     rvfarbe <- rvfarbe[order(rvfarbe)]
#     rvfarbe <- farben$farbe[match(rvfarbe, farben$phylum)]
#     png(paste("Maraike_ROut/sakey_winter_", names(sakeys.winter.ls)[i], ".png", sep = ""))
#     print(ggplot(sakeys.winter.ls[[i]][which(sakeys.winter.ls[[i]]$Freq > 1),], aes(axis1=fungus, axis2=bacterium, y=Freq)) + geom_alluvium(aes(fill = bacterium)) +  geom_stratum(width = 1/3) + theme_void() + geom_stratum(aes(fill = bacterium)) + geom_text(stat = "stratum", aes(label = after_stat(stratum))) + theme(legend.position="none") + scale_fill_manual(values = rvfarbe, guide = "none"))
#     dev.off()
#   }
# }
# 
# save(sakeys.summer.ls, file = "sakeys.summer.ls.RData")
# save(sakeys.winter.ls, file = "sakeys.winter.ls.RData")
# save(sakeys.neg.summer.ls, file = "sakeys.neg.summer.ls.RData")
# save(sakeys.pos.summer.ls, file = "sakeys.pos.summer.ls.RData")
# save(sakeys.neg.winter.ls, file = "sakeys.neg.winter.ls.RData")
# save(sakeys.pos.winter.ls, file = "sakeys.pos.winter.ls.RData")
# 

#Modelling the numbers of positive and negative and total edges between different phyla
load(file = "summergraphs.core.ls.RData")
load(file = "wintergraphs.core.ls.RData")
load(file = "summeredges.core.ls.RData")
load(file = "winteredges.core.ls.RData")
load(file = "sakeys.summer.ls.RData")
load(file = "sakeys.winter.ls.RData")
#load(file = "sakeys.neg.summer.ls.RData")
#load(file = "sakeys.pos.summer.ls.RData")
#load(file = "sakeys.neg.winter.ls.RData")
#load(file = "sakeys.pos.winter.ls.RData")


sakeys.ls <- vector("list", length = length(sakeys.summer.ls))
names(sakeys.ls) <- names(sakeys.summer.ls)
for(i in 1:length(sakeys.ls)){
  sakeys.ls[[i]] <- data.frame(rbind(sakeys.summer.ls[[i]], sakeys.winter.ls[[i]]))
  sakeys.ls[[i]]$season <- c(rep("summer", nrow(sakeys.summer.ls[[i]])), rep("winter", nrow(sakeys.winter.ls[[i]])))
  if(length(which(sakeys.ls[[i]]$fungus == "NA")) > 0){
    sakeys.ls[[i]] <- sakeys.ls[[i]][-which(sakeys.ls[[i]]$fungus == "NA"),]
  }
  if(length(which(sakeys.ls[[i]]$bacterium == "NA")) > 0){
    sakeys.ls[[i]] <- sakeys.ls[[i]][-which(sakeys.ls[[i]]$bacterium == "NA"),]
  }
  sakeys.ls[[i]]$pair <- as.factor(paste(sakeys.ls[[i]]$fungus, sakeys.ls[[i]]$bacterium, sep = "_"))
  sakeys.ls[[i]]$fungus <- as.factor(sakeys.ls[[i]]$fungus)
  sakeys.ls[[i]]$bacterium <- as.factor(sakeys.ls[[i]]$bacterium)
  sakeys.ls[[i]]$season <- as.factor(sakeys.ls[[i]]$season)
}


#Poisson model to check if there are differences in numbers of associations among bacterial and fungal phyla depending on season.

pglm.int.ls <- vector("list", length = length(sakeys.ls))
names(pglm.int.ls) <- names(sakeys.ls)
pglm.sf.ls <- vector("list", length = length(sakeys.ls))
names(pglm.sf.ls) <- names(sakeys.ls)
pglm.s.ls <- vector("list", length = length(sakeys.ls))
names(pglm.s.ls) <- names(sakeys.ls)
pglm.f.ls <- vector("list", length = length(sakeys.ls))
names(pglm.f.ls) <- names(sakeys.ls)
for(i in 1:length(sakeys.ls)){
  pglm.int.ls[[i]] <- glm(Freq ~ season*fungus, family = poisson, data = sakeys.ls[[i]])
  pglm.sf.ls[[i]] <- glm(Freq ~ season + fungus, family = poisson, data = sakeys.ls[[i]])
  pglm.f.ls[[i]] <- glm(Freq ~ fungus, family = poisson, data = sakeys.ls[[i]])
  pglm.s.ls[[i]] <- glm(Freq ~ season, family = poisson, data = sakeys.ls[[i]])
}

poi.ov.ls <- vector("list", length = length(pglm.int.ls))
names(poi.ov.ls) <- names(pglm.int.ls)
for(j in 1:length(poi.ov.ls)){
  poi.ov.ls[[j]] <- data.frame(matrix(data = NA, ncol = 6, nrow = 10))
  names(poi.ov.ls[[j]]) <- c("interaction", "season_fungus", "season", "fungus", "null", "variable")
  raov <- summary(aov(pglm.int.ls[[j]]))
  rs <- summary(pglm.int.ls[[j]])
  poi.ov.ls[[j]][5,1] <- raov[[1]][1,5] #pvalue season
  poi.ov.ls[[j]][7,1] <- raov[[1]][2,5] #pvalue fungus
  poi.ov.ls[[j]][9,1] <- raov[[1]][3,5] #pvalue fungus*season
  poi.ov.ls[[j]][4,1] <- raov[[1]][1,2] #ss season
  poi.ov.ls[[j]][6,1] <- raov[[1]][2,2] #ss fungus
  poi.ov.ls[[j]][8,1] <- raov[[1]][3,2] #ss fungus*season
  poi.ov.ls[[j]][10,1] <- raov[[1]][4,2] #ss res
  poi.ov.ls[[j]][2,1] <- rs$deviance
  poi.ov.ls[[j]][1,1] <- rs$null.deviance
  poi.ov.ls[[j]][3,1] <- rs$aic
  
  raov <- summary(aov(pglm.sf.ls[[j]]))
  rs <- summary(pglm.sf.ls[[j]])
  poi.ov.ls[[j]][5,2] <- raov[[1]][1,5] #pvalue season
  poi.ov.ls[[j]][7,2] <- raov[[1]][2,5] #pvalue fungus
  poi.ov.ls[[j]][4,2] <- raov[[1]][1,2] #ss season
  poi.ov.ls[[j]][6,2] <- raov[[1]][2,2] #ss fungus
  poi.ov.ls[[j]][10,2] <- raov[[1]][3,2] #ss res
  poi.ov.ls[[j]][2,2] <- rs$deviance
  poi.ov.ls[[j]][1,2] <- rs$null.deviance
  poi.ov.ls[[j]][3,2] <- rs$aic
  
  raov <- summary(aov(pglm.s.ls[[j]]))
  rs <- summary(pglm.s.ls[[j]])
  poi.ov.ls[[j]][5,3] <- raov[[1]][1,5] #pvalue season
  poi.ov.ls[[j]][4,3] <- raov[[1]][1,2] #ss season
  poi.ov.ls[[j]][10,3] <- raov[[1]][2,2] #ss res
  poi.ov.ls[[j]][2,3] <- rs$deviance
  poi.ov.ls[[j]][1,3] <- rs$null.deviance
  poi.ov.ls[[j]][3,3] <- rs$aic
  
  raov <- summary(aov(pglm.f.ls[[j]]))
  rs <- summary(pglm.f.ls[[j]])
  poi.ov.ls[[j]][7,4] <- raov[[1]][1,5] #pvalue fungus
  poi.ov.ls[[j]][6,4] <- raov[[1]][1,2] #ss fungus
  poi.ov.ls[[j]][10,4] <- raov[[1]][2,2] #ss res
  poi.ov.ls[[j]][2,4] <- rs$deviance
  poi.ov.ls[[j]][1,4] <- rs$null.deviance
  poi.ov.ls[[j]][3,4] <- rs$aic
  
  rmod <- glm(Freq ~ 1, family = poisson, data = sakeys.ls[[j]])
  raov <- summary(aov(rmod))
  rs <- summary(rmod)
  poi.ov.ls[[j]][10,5] <- raov[[1]][1,2] #ss res
  poi.ov.ls[[j]][2,5] <- rs$deviance
  poi.ov.ls[[j]][1,5] <- rs$null.deviance
  poi.ov.ls[[j]][3,5] <- rs$aic
  
  poi.ov.ls[[j]]$variable <- c("null.deviance", "deviance", "aic", "ss.season", "p.season", "ss.fungus", "p.fungus", "ss.seasonfungus", "p.seasonfungus", "ss.resid")
}

write_xlsx(poi.ov.ls, path = "Maraike_ROut/poi.ov.ls.xlsx")

bm <- c()
for(i in 1:length(poi.ov.ls)){
  bm <- append(bm, which(poi.ov.ls[[i]][3,1:5] == min(as.numeric(poi.ov.ls[[i]][3,1:5]))), length(bm))
}
table(bm)

sink("Maraike_ROut/Poisson_model_coefficients.txt")
for(i in c(1,2,6,7,11,12,16,17,21)){
  print(names(pglm.int.ls)[i])
  print(summary(pglm.int.ls[[i]]))
}
sink()


####Preparing a boxplot for the visualization of the poisson model results across thresholds
pglm.int.estimates.df <- data.frame(sakeys.ls[[1]], predict(pglm.int.ls[[1]]))
pglm.int.estimates.df$fungusseason <- paste(pglm.int.estimates.df$fungus, pglm.int.estimates.df$season, sep = "_")
names(pglm.int.estimates.df)[7] <- "prediction"
pglm.int.estimates.df <- pglm.int.estimates.df[!duplicated(pglm.int.estimates.df$fungusseason),]
pglm.int.estimates.df$estimate <- exp(pglm.int.estimates.df$prediction)
pglm.int.estimates.df <- pglm.int.estimates.df[,c(3,5,8,9)]
pglm.int.estimates.df$group <- names(pglm.int.ls)[1]
f <- pglm.int.estimates.df[which(pglm.int.estimates.df$fungus %in% as.character(data.frame(table(pglm.int.estimates.df$fungus))$Var1[which(data.frame(table(pglm.int.estimates.df$fungus))$Freq == 2)])),]
fs <- f$estimate[1:(nrow(f)/2)]
fw <- f$estimate[(nrow(f)-(nrow(f)/2)+1):nrow(f)]
pglm.int.estimates.df$fsfw <- NA
pglm.int.estimates.df$fsfw[which(pglm.int.estimates.df$fungus %in% as.character(data.frame(table(pglm.int.estimates.df$fungus))$Var1[which(data.frame(table(pglm.int.estimates.df$fungus))$Freq == 2)]))] <- c(fs/fw, fs/fw)


for(i in c(2,5,6,7,8,11,12,13,16,17,21)){
  rdf <- data.frame(sakeys.ls[[i]], predict(pglm.int.ls[[i]]))
  rdf$fungusseason <- paste(rdf$fungus, rdf$season, sep = "_")
  names(rdf)[7] <- "prediction"
  rdf <- rdf[!duplicated(rdf$fungusseason),]
  rdf$estimate <- exp(rdf$prediction)
  rdf <- rdf[,c(3,5,8,9)]
  rdf$group <- names(pglm.int.ls)[i]
  
  f <- rdf[which(rdf$fungus %in% as.character(data.frame(table(rdf$fungus))$Var1[which(data.frame(table(rdf$fungus))$Freq == 2)])),]
  fs <- f$estimate[1:(nrow(f)/2)]
  fw <- f$estimate[(nrow(f)-(nrow(f)/2)+1):nrow(f)]
  rdf$fsfw <- NA
  rdf$fsfw[which(rdf$fungus %in% as.character(data.frame(table(rdf$fungus))$Var1[which(data.frame(table(rdf$fungus))$Freq == 2)]))] <- c(fs/fw, fs/fw)
  
  pglm.int.estimates.df <- data.frame(rbind(pglm.int.estimates.df, rdf))
}


for(i in c(3,4,9,10,14,18,19,22:24)){
  rdf <- data.frame(sakeys.ls[[i]], predict(pglm.f.ls[[i]]))
  rdf$fungusseason <- paste(rdf$fungus, rdf$season, sep = "_")
  names(rdf)[7] <- "prediction"
  rdf <- rdf[!duplicated(rdf$fungusseason),]
  rdf$estimate <- exp(rdf$prediction)
  rdf <- rdf[,c(3,5,8,9)]
  rdf$group <- names(pglm.int.ls)[i]
  
  f <- rdf[which(rdf$fungus %in% as.character(data.frame(table(rdf$fungus))$Var1[which(data.frame(table(rdf$fungus))$Freq == 2)])),]
  fs <- f$estimate[1:(nrow(f)/2)]
  fw <- f$estimate[(nrow(f)-(nrow(f)/2)+1):nrow(f)]
  rdf$fsfw <- NA
  rdf$fsfw[which(rdf$fungus %in% as.character(data.frame(table(rdf$fungus))$Var1[which(data.frame(table(rdf$fungus))$Freq == 2)]))] <- c(fs/fw, fs/fw)
  
  pglm.int.estimates.df <- data.frame(rbind(pglm.int.estimates.df, rdf))
}
#write.csv(pglm.int.estimates.df, file = "Maraike_ROut/pglm_int_estimates_df.csv")

table(pglm.int.estimates.df$fungus) #this has to be divided by 2, because of season.

png("Maraike_ROut/Poisson_model_funphylum_season_phylumassos.png", width = 200)
print(ggplot(data = pglm.int.estimates.df[which(pglm.int.estimates.df$fungus %in% c("Ascomycota", "Basidiomycota", "Enthorrhizomycota", "Mortierellomycota", "Mucoromycota")),], aes(y = fsfw, x = fungus)) + geom_boxplot(aes(fill = fungus)) + scale_fill_manual(values = c("pink", "maroon", "mediumpurple", "thistle")) + geom_point(size = 2) + theme(legend.position = "none", axis.text = element_text(size = 20), axis.title = element_text(size = 22), panel.background = element_rect(fill = "white", color = "grey"), panel.grid.major = element_line(color = "grey", linetype = "dotted", size = 0.5), panel.grid.minor = element_line(color = "grey", size = 0.5, linetype = "dotted"), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + xlab("") + ylab("Snow free/snow-covered \n associations") + scale_y_continuous(limits = c(0.5, 1.5)))
dev.off()


png("Maraike_ROut/Poisson_model_funphylum_season_phylumassos_all.png", width = 550, height = 450)
print(ggplot(data = pglm.int.estimates.df[which(pglm.int.estimates.df$fungus != c("Basidiobolomycota")),], aes(y = fsfw, x = fungus)) + geom_boxplot(aes(fill = fungus)) + scale_fill_manual(values = c("pink", "maroon", rep("grey", 4), "mediumpurple", "thistle", rep("grey", 3))) + geom_point(size = 2) + theme(legend.position = "none", axis.text = element_text(size = 20), axis.title = element_text(size = 22), panel.background = element_rect(fill = "white", color = "grey"), panel.grid.major = element_line(color = "grey", linetype = "dotted", size = 0.5), panel.grid.minor = element_line(color = "grey", size = 0.5, linetype = "dotted"), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + xlab("") + ylab("Snow free/snow-covered \n associations") + scale_y_continuous(limits = c(0.5, 1.5)))
dev.off()


#####Binomial model for considering weights:
sakeys.bino.ls <- vector("list", length = length(sakeys.ls))
names(sakeys.bino.ls) <- names(sakeys.ls)
for(i in 1:length(sakeys.ls)){
  rdfcore.s <- summeredges.core.ls[[d[i]]]
  rdfcore.s <- rdfcore.s[!is.na(rdfcore.s$phylum_fun),]
  rdfcore.s <- rdfcore.s[!is.na(rdfcore.s$phylum_bac),]
  rdfcore.s$phylumpair <- as.factor(gsub("p__", "", rdfcore.s$phylumpair))
  rdf.s <- data.frame(table(rdfcore.s[which(rdfcore.s$weight < 0),]$phylumpair))
  names(rdf.s) <- c("pair", "neg")
  rdf.s$pos <- data.frame(table(rdfcore.s[which(rdfcore.s$weight > 0),]$phylumpair))[,2]
  rdfcore.w <- winteredges.core.ls[[d[i]]]
  rdfcore.w <- rdfcore.w[!is.na(rdfcore.w$phylum_fun),]
  rdfcore.w <- rdfcore.w[!is.na(rdfcore.w$phylum_bac),]
  rdfcore.w$phylumpair <- as.factor(gsub("p__", "", rdfcore.w$phylumpair))
  rdf.w <- data.frame(table(rdfcore.w[which(rdfcore.w$weight < 0),]$phylumpair))
  names(rdf.w) <- c("pair", "neg")
  rdf.w$pos <- data.frame(table(rdfcore.w[which(rdfcore.w$weight > 0),]$phylumpair))[,2]
  
  rdf <- data.frame(matrix(data = 0, ncol = 6, nrow = (length(unique(c(levels(rdfcore.s$phylumpair), levels(rdfcore.w$phylumpair))))*2)))
  names(rdf) <- c("pair", "fungus", "bacterium", "season", "freq.pos", "freq.neg")
  rdf$pair <- rep(unique(c(levels(rdfcore.s$phylumpair), levels(rdfcore.w$phylumpair))),2)
  rdf$season <- c(rep("summer", (nrow(rdf)/2)), rep("winter", (nrow(rdf)/2)))
  rdf$fungus <- matrix(unlist(strsplit(rdf$pair, "_")), nrow = 2)[1,]
  rdf$bacterium <- matrix(unlist(strsplit(rdf$pair, "_")), nrow = 2)[2,]
  rdf$freq.pos[which(rdf$season == "summer")] <- rdf.s$pos[match(rdf$pair[which(rdf$season == "summer")], rdf.s$pair)]
  rdf$freq.neg[which(rdf$season == "summer")] <- rdf.s$neg[match(rdf$pair[which(rdf$season == "summer")], rdf.s$pair)]
  rdf$freq.pos[which(rdf$season == "winter")] <- rdf.w$pos[match(rdf$pair[which(rdf$season == "winter")], rdf.w$pair)]
  rdf$freq.neg[which(rdf$season == "winter")] <- rdf.w$neg[match(rdf$pair[which(rdf$season == "winter")], rdf.w$pair)]
  sakeys.bino.ls[[i]] <- rdf
}


bglm.ls <- vector("list", length = length(sakeys.ls))
names(bglm.ls) <- names(sakeys.ls)
for(i in 1:length(bglm.ls)){
  rdf <- sakeys.bino.ls[[i]]
  rdf <- rdf[which(rdf$fungus == "Basidiomycota"),]
  
  bglm.ls[[i]] <- glm(cbind(freq.pos, freq.neg) ~ season + pair, data = rdf, family = binomial)
}




bino.set <- sakeys.bino.ls[[1]]
for(i in 2:length(sakeys.bino.ls)){
  bino.set <- rbind(bino.set, sakeys.bino.ls[[i]])
}
bino.set[is.na(bino.set)] <- 0

#This overall modell turns out insignificant. So, this cannot be used for the publication.
bglm <- glm(cbind(freq.pos, freq.neg) ~ season * pair, data = bino.set, family = binomial)
summary(bglm)
summary(aov(bglm))
s <- step(bglm)

bglm.pred <- predict(bglm, type = "response", se.fit = T)
bino.set.2 <- data.frame(bino.set, bglm.pred)
bino.set.2$factor <- paste(bino.set.2$pair, bino.set.2$season, sep = "|")

bino.mort <- bino.set.2[which(bino.set.2$fungus == "Mortierellomycota"),]
y0 <- c(bino.mort$fit-bino.mort$se.fit)[!duplicated(bino.mort$factor)]
y0 <- y0[order(unique(bino.mort$factor))]
y1 <- c(bino.mort$fit+bino.mort$se.fit)[!duplicated(bino.mort$factor)]
y1 <- y1[order(unique(bino.mort$factor))]

x11(); par(mar = c(15,5,2,2)); plot(fit ~ as.factor(factor), data = bino.mort, las = 3, xlab = "", ylab = "Probability of positive association", cex.axis = 0.7, ylim = c(0,1))
arrows(x0=1:length(unique(bino.mort$factor)), y0=y0, x1=1:length(unique(bino.mort$factor)), y1=y1, code=3, angle=90, length=0.1)


# sink("Maraike_ROut/Network_binomial_posneg_associations_season.txt")
# print(summary(aov(gmbi)))
# print(summary(gmbi))
# print(step(gmbi))
# sink()

bino.set.mort <- bino.set[which(bino.set$fungus == "Mortierellomycota"),]
bino.set.mort <- bino.set.mort[which((bino.set.mort$freq.pos + bino.set.mort$freq.neg) > 2),]
bglm.mort <- glm(cbind(freq.pos, freq.neg) ~ season + pair + season*pair, data = bino.set.mort, family = binomial)
summary(bglm)
summary(aov(bglm))
s <- step(bglm.mort, direction = "both")

bino.set.mort <- data.frame(bino.set.mort, predict(bglm.mort, type = "response", se.fit = T))
bino.set.mort$factor <- paste(bino.set.mort$pair, bino.set.mort$season, sep = "|")

y0 <- c(bino.set.mort$fit-bino.set.mort$se.fit)[!duplicated(bino.set.mort$factor)]
y0 <- y0[order(unique(bino.set.mort$factor))]
y1 <- c(bino.set.mort$fit+bino.set.mort$se.fit)[!duplicated(bino.set.mort$factor)]
y1 <- y1[order(unique(bino.set.mort$factor))]

x11(); par(mar = c(15,5,2,2)); plot(fit ~ as.factor(factor), data = bino.set.mort, las = 3, xlab = "", ylab = "Probability of positive association", cex.axis = 0.7, ylim = c(0,1))
arrows(x0=1:length(unique(bino.set.mort$factor)), y0=y0, x1=1:length(unique(bino.set.mort$factor)), y1=y1, code=3, angle=90, length=0.1)






####More network stats- focussing on network closeness, betweenness and clustering
x11(); hist(closeness(wintergraphs.ls$'2_3_10'))
x11(); hist(betweenness(wintergraphs.ls$'2_3_10'))
transitivity(wintergraphs.ls$'2_3_10')

for(i in 1:length(wintergraphs.ls)){
  if(class(wintergraphs.ls[[i]]) == "igraph"){
    print("winter")
    print(transitivity(wintergraphs.ls[[i]]))
    print("summer")
    print(transitivity(summergraphs.ls[[i]]))
  }
}


for(i in 1:length(wintergraphs.ls)){
  if(class(wintergraphs.ls[[i]]) == "igraph"){
    print("winter")
    print(summary(closeness(wintergraphs.ls[[i]])))
    print("summer")
    print(summary(closeness(summergraphs.ls[[i]])))
  }
}

for(i in 1:length(wintergraphs.ls)){
  if(class(wintergraphs.ls[[i]]) == "igraph"){
    print("winter")
    print(summary(betweenness(wintergraphs.ls[[i]])))
    print("summer")
    print(summary(betweenness(summergraphs.ls[[i]])))
  }
}

