# Libraries----

library(pander)
library(ape)
library(phytools)
library(phylobase)
library(vegan)
library(ggplot2)
library(dplyr)
library(reshape2)
library(pheatmap)
library(dnar)
library(Rtsne)
library(lattice)
library(grid)
library(gridExtra)
library(ggimage)
library(phyloseq)
library(eclectic)
library(ade4)

# Setup----

work_dir <- "/home/kevin/projects/bilateria/"

mapping_file_fp <- file.path(work_dir, "map/map.tsv")
otu_table_fp <- file.path(work_dir, "otu/otu_table.txt")
uu_fp <- file.path(work_dir, "beta_diversity/unweighted_unifrac_otu_table.original.txt")
wu_fp <- file.path(work_dir, "beta_diversity/weighted_unifrac_otu_table.original.txt")

# Functions----

resetValidationData <- function(E, status = "all") {
  library(qiimer)
  E$o <- read_qiime_otu_table(E$otuTableFile)
  E$s <- read_qiime_mapping_file(E$mapFile)
  if(status=="toUse") {
    E$s <- E$s[E$s$toKeep=="keep" & !is.na(E$s$weight),]
  }
  
  
  alignQiimeData(E)
  if(!status=="all") {
    E$o$metadata <- E$o$metadata[!grepl("*Chloroplast*",E$o$metadata)]
  }
  E$o$counts <- E$o$counts[names(E$o$metadata),]
  E$o$otu_ids <- names(E$o$metadata)
  
  E$sampleOrdering <- with(E$s,order(ave(1:nrow(E$s),E$s$phylum,FUN=length),phylum,
                                         ave(1:nrow(E$s),E$s$class,FUN=length),class,
                                         ave(1:nrow(E$s),E$s$order,FUN=length),order,
                                         ave(1:nrow(E$s),E$s$family,FUN=length),family,
                                         ave(1:nrow(E$s),E$s$genus,FUN=length),genus,
                                         ave(1:nrow(E$s),E$s$species,FUN=length),species,
                                         ave(1:nrow(E$s),E$s$common,FUN=length),common))
  
  E$speciesTax <- unique(E$s[,c("phylum","class","order","family","genus","common")])
  E$speciesTax <- unique(E$speciesTax)
  E$speciesOrdering <- with(E$speciesTax,
                              order(ave(1:nrow(E$speciesTax),E$speciesTax$phylum,FUN=length),phylum,
                                    ave(1:nrow(E$speciesTax),E$speciesTax$class,FUN=length),class,
                                    ave(1:nrow(E$speciesTax),E$speciesTax$order,FUN=length),order,
                                    ave(1:nrow(E$speciesTax),E$speciesTax$family,FUN=length),family,
                                    ave(1:nrow(E$speciesTax),E$speciesTax$genus,FUN=length),genus))
  
  #  E$figure1SpeciesOrdering <- match(E$,E$speciesTax$common)
  E$s$weight <- as.numeric(E$s$weight)
  E$s$log.weight <- log10(E$s$weight)
  E$s$numOtus <- colSums(E$o$counts>0)
  E$s$filteredReadCount <- colSums(E$o$counts)
  
  md <- sub("(; [kpcofgs]__)+$", "", E$o$metadata, perl=T)
  adf <- split_assignments(md)
  E$a <- simplify_assignments(adf)
  rm(md); rm(adf)
}

resetData <- function(status = "all") {
  library(qiimer)
  library(phylobase)
  library(adephylo)
  
  G <<- new.env()
  G$s <- read_qiime_mapping_file(mapping_file_fp)
  if(status=="toUse") {
    G$s <- G$s[G$s$toKeep=="keep" & !is.na(G$s$Weight_to_use),]
  }
  
  if(status=="toUsePlusControls") {
    G$s <- G$s[G$s$toKeep=="keep" | G$s$sample_type %in% c("GeneBlock","Extraction_blank",
                                                           "Dissection_blank", "PCR_H2O"),]
  }

  # Code used to create om.csv (om = OTU melted)
  # om <- melt(G$o$counts)
  # om <- om[om$value>0,]
  # om$metadata <- G$o$metadata[om$Var1]
  # colnames(om) <- c("OTU","SampleID","count","metadata")
  # write.csv(om, paste0(work_dir,"om.csv"), row.names = F)

  # Code to read in om.csv and create G$o object  
  # om <- read.csv(paste0(work_dir,"om.csv"))
  # o <- dcast(data = om,formula = OTU~SampleID,fun.aggregate = sum,value.var = "count")
  # rownames(o) <- o$OTU
  # G$o$counts <- as.matrix(o)
  # uniq_metadata <- unique(om[,c("OTU","metadata")])
  # G$o$sample_ids <- colnames(G$o$counts)
  # G$o$otu_ids <- rownames(G$o$counts)
  # G$o$metadata <- as.character(uniq_metadata$metadata)
  # names(G$o$metadata) <- as.character(uniq_metadata$OTU)
  # rm(uniq_metadata)
  # rm(o)
  # rm(om)
  
  # Use this code to create OTU table to save time on future loads
  # write.table(cbind(G$o$counts,G$o$metadata),
  #            paste0(work_dir,"otu/otu_table.txt"), quote=FALSE, sep = "\t")
  # system(paste0("sed -i '1s/^/# Constructed from biom file\\n#OTU ID\t/' ",work_dir,"otu/otu_table.txt"))
  # system(paste0("sed -i '2s/$/Consensus Lineage/' ",work_dir,"otu/otu_table.txt"))
  
  # Use this line in place of the above code block defining G$o if you have already created OTU table
  G$o <- read_qiime_otu_table(otu_table_fp)
  
  alignQiimeData(G)
  rownames(G$s) <- G$s$SampleID
  colnames(G$o$counts) <- G$s$SampleID
  G$o$sample_ids <- G$s$SampleID
  if(!status=="all") {
    G$o$metadata <- G$o$metadata[!grepl("*Chloroplast*",G$o$metadata)]
  }
  G$o$counts <- G$o$counts[names(G$o$metadata),]
  G$o$otu_ids <- names(G$o$metadata)
  
  G$sampleOrdering <- with(G$s,order(ave(1:nrow(G$s),G$s$phylum,FUN=length),phylum,
                                     ave(1:nrow(G$s),G$s$supersuperclass,FUN=length),supersuperclass,
                                     ave(1:nrow(G$s),G$s$superclass,FUN=length),superclass,
                                     ave(1:nrow(G$s),G$s$class,FUN=length),class,
                                     ave(1:nrow(G$s),G$s$superorder,FUN=length),superorder,
                                     ave(1:nrow(G$s),G$s$clade,FUN=length),clade,
                                     ave(1:nrow(G$s),G$s$order,FUN=length),order,
                                     ave(1:nrow(G$s),G$s$family,FUN=length),family,
                                     ave(1:nrow(G$s),G$s$genus,FUN=length),genus,
                                     ave(1:nrow(G$s),G$s$species,FUN=length),species,
                                     ave(1:nrow(G$s),G$s$common,FUN=length),common))
  
  G$speciesTax <- unique(G$s[,c("phylum","supersuperclass","superclass","class","superorder","clade",
                                "order","family","genus","common","timetree_proxy")])
  G$speciesTax$genus <- as.character(G$speciesTax$genus)
  G$speciesTax$genus[G$speciesTax$common=="Bee"] <- "Bee genus"
  G$speciesTax <- unique(G$speciesTax)
  G$speciesOrdering <- with(G$speciesTax,
                            order(ave(1:nrow(G$speciesTax),G$speciesTax$phylum,FUN=length),phylum,
                                  ave(1:nrow(G$speciesTax),G$speciesTax$supersuperclass,FUN=length),supersuperclass,
                                  ave(1:nrow(G$speciesTax),G$speciesTax$superclass,FUN=length),superclass,
                                  ave(1:nrow(G$speciesTax),G$speciesTax$class,FUN=length),class,
                                  ave(1:nrow(G$speciesTax),G$speciesTax$superorder,FUN=length),superorder,
                                  ave(1:nrow(G$speciesTax),G$speciesTax$clade,FUN=length),clade,
                                  ave(1:nrow(G$speciesTax),G$speciesTax$order,FUN=length),order,
                                  ave(1:nrow(G$speciesTax),G$speciesTax$family,FUN=length),family,
                                  ave(1:nrow(G$speciesTax),G$speciesTax$genus,FUN=length),genus))

  G$s$speciesOrder <- match(G$s$common,G$speciesTax$common[G$figure1SpeciesOrdering])
    
  md <- sub("(; [kpcofgs]__)+$", "", G$o$metadata, perl=T)
  adf <- split_assignments(md)
  G$a <- simplify_assignments(adf)
  rm(md); rm(adf)
  
  G$newick <- read.tree(paste0(work_dir, "speciesPhylogeny.nwk"))

  G$timetreeLabels <-
    tolower(gsub("_"," ",
                 G$newick$tip.label[G$newick$edge[,2][G$newick$edge[,2]<=length(G$newick$tip.label)]]))
  
  G$figure1SpeciesOrdering <- match(G$timetreeLabels,G$speciesTax$timetree_proxy)
  G$s$Weight_to_use <- as.numeric(G$s$Weight_to_use)
  G$s$log.weight <- log10(G$s$Weight_to_use)
  
  G$species <-
    as.data.frame(unique(G$s[match(G$speciesTax$common[G$figure1SpeciesOrdering],G$s$common),
                             c("common","phylum","class","order","family","diet","specificDiet","timetree_proxy")]))
  G$species <- G$species[G$species$common %in% G$speciesTax$common[G$figure1SpeciesOrdering],]
  rownames(G$species) <- G$species$common
}

resetUnifracData <- function(){
  library(qiimer)
  G$uur_fp <- file.path(work_dir, "bd_uur", "unweighted_unifrac_otu_rare.txt")
  G$uur <- read_qiime_distmat(G$uur_fp)
  G$uur <- dist_subset(G$uur, G$s$SampleID)

  G$wnu_fp <- file.path(work_dir, "bd_wnu", "weighted_unifrac_otu_prop.txt")
  G$wnu <- read_qiime_distmat(G$wnu_fp)
  G$wnu <- dist_subset(G$wnu, G$s$SampleID)
}

alignQiimeData <- function(E) {
  o <- E$o
  s <- E$s
  if(!is.list(o)) stop("Object o is not of type list")
  if(!is.data.frame(s)) stop("Object s is not of type data frame")
  if(!"sample_ids" %in% names(o)) stop("Object o needs element named 'sample_ids'")
  if(!"otu_ids" %in% names(o)) stop("Object o needs element named 'otu_ids'")
  if(!"metadata" %in% names(o)) stop("Object o needs element named 'metadata'")
  if(!"counts" %in% names(o)) stop("Object o needs element named 'counts'")
  if(!"SampleID" %in% colnames(s)) stop("Object s needs column named 'SampleID'")
  if(is.null(colnames(o$counts))) stop("o$counts needs named columns")
  if(is.null(rownames(o$counts))) stop("o$counts needs named rows")
  if(is.null(names(o$metadata))) stop("o$metadata must be a named vector")
  
  s <- s[s$SampleID %in% colnames(o$counts),]
  s <- droplevels(s)
  
  o$counts <- o$counts[,colnames(o$counts) %in% s$SampleID]
  o$counts <- o$counts[rowSums(o$counts)>0,]
  o$counts <- o$counts[,match(s$SampleID, colnames(o$counts))]
  o$sample_ids <- colnames(o$counts)
  o$metadata <- o$metadata[rownames(o$counts)]
  o$otu_ids <- o$otu_ids[o$otu_ids %in% rownames(o$counts)]
  E$s <- s
  E$o <- o
}

expandTree <- function(tree,replacingLabel,newLabel,newEdgeLength){
  replacingTip <- which(tree$tip.label==replacingLabel)
  replacingEdge <- which(tree$edge[,2]==replacingTip)
  
  numNodes <- max(tree$edge)
  newNode <- numNodes + 2
  newTip <- replacingTip + 1
  numTips <- length(tree$tip.label)
  
  tree$edge[tree$edge > replacingTip] <- tree$edge[tree$edge > replacingTip] + 1
  tree$edge[tree$edge==replacingTip] <- newNode
  tree$edge <- cbind(append(append(tree$edge[,1],newNode,replacingEdge),newNode,replacingEdge),
                     append(append(tree$edge[,2],newTip,replacingEdge),replacingTip,replacingEdge))
  
  tree$edge.length[replacingEdge] <- tree$edge.length[replacingEdge] - newEdgeLength
  tree$edge.length <- append(append(tree$edge.length,newEdgeLength,replacingEdge),newEdgeLength,replacingEdge)
  
  tree$tip.label <- append(tree$tip.label,newLabel,replacingTip)
  tree$Nnode <- tree$Nnode + 1
  tree$node.label <- c(tree$node.label, paste0("'",max(as.numeric(gsub("'","",tree$node.label[-1])))+1,"'"))
  
  return(tree)
}

expandTreeInternal <- function(tree,replacingNode,newLabel,newTipEdgeLength){
  replacingEdge <- which(tree$edge[,2]==replacingNode)
  replacingNodeLeaves <- as(subset(as(tree,"phylo4"),node.subtree=replacingNode),"phylo")$tip.label
  newTip <- max(match(replacingNodeLeaves,tree$tip.label)) + 1
  
  numNodes <- max(tree$edge)
  newNode <- numNodes + 2
  numTips <- length(tree$tip.label)
  replacingNodeHeight <- max(node.depth.edgelength(as(subset(as(tree,"phylo4"),
                                                             node.subtree=replacingNode),"phylo")))
  newInternalEdgeLength <- newTipEdgeLength - replacingNodeHeight
  newReplacingEdgeLength <- tree$edge.length[replacingEdge] - newInternalEdgeLength
  
  
  tree$edge[tree$edge >= newTip] <- tree$edge[tree$edge >= newTip] + 1
  replacingNode <- replacingNode + 1
  tree$edge[tree$edge[,1]==replacingNode,1] <- newNode
  tree$edge <- cbind(append(append(tree$edge[,1],replacingNode,replacingEdge),replacingNode,replacingEdge),
                     append(append(tree$edge[,2],newNode,replacingEdge),newTip,replacingEdge))
  
  tree$edge.length <- append(append(tree$edge.length,newInternalEdgeLength,replacingEdge),
                             newTipEdgeLength,replacingEdge)
  tree$edge.length[replacingEdge] <- newReplacingEdgeLength
  
  tree$tip.label <- append(tree$tip.label,newLabel,newTip - 1)
  tree$Nnode <- tree$Nnode + 1
  tree$node.label <- c(tree$node.label, paste0("'",max(as.numeric(gsub("'","",tree$node.label[-1])))+1,"'"))
  
  return(tree)
}

plotFigure1 <- function(E, study) {
  E$aPhylum <- sub("; c__(.)*$", "", E$o$metadata, perl=T)
  E$oByPhylum <- rowsum(E$o$counts,E$aPhylum)
  
  E$obpProp <- apply(E$oByPhylum,2,function(x) ifelse (x,x / sum(x),0))
  E$obpProp <- data.frame(t(E$obpProp[rowSums(E$obpProp)>0,]))
  
  E$obpProp$common <- as.character(E$s$common[match(rownames(E$obpProp),E$s$SampleID)])
  E$obpPropBySpecies <- aggregate(.~common,data=E$obpProp,mean)
  
  E$obpMaxProps <- apply(E$obpPropBySpecies[,-1],2,max)
  
  numOtuPhylaToUse <- 10
  E$otuPhylaToUse <- G$otuPhylaToUse
  
  E$obpPropBySpeciesToUse <- E$obpPropBySpecies[,names(E$obpPropBySpecies) == 'common' |
                                                  names(E$obpPropBySpecies) %in% E$otuPhylaToUse]
  
  colnames(E$obpPropBySpeciesToUse)[-1] <-
    colnames(E$obpPropBySpeciesToUse)[rev(order(colSums(E$obpPropBySpeciesToUse[-1])))+1]
  E$obpPropBySpeciesToUse[,-1] <-
    E$obpPropBySpeciesToUse[,rev(order(colSums(E$obpPropBySpeciesToUse[-1])))+1]
  
  E$obpPropBySpeciesToUse$Other <- 1 - rowSums(E$obpPropBySpeciesToUse[-1])
  
  E$obpHeatmapData <- melt(cbind(E$obpPropBySpeciesToUse), id.vars=c('common'))
  
  plotHeatmap <- ggplot(E$obpHeatmapData, aes(x=variable, y=common, fill=value)) +
    scale_y_discrete(limits=E$speciesTax$common[E$speciesOrdering]) +
    scale_x_discrete(limits=c(G$otuPhylaToUse,"Other")) +
    geom_tile(color="grey80", size=0.4) + theme_grey() +
    theme(
      legend.position = "none",
      strip.text.y = element_text(angle=0, vjust=0),
      strip.text.x = element_text(angle=90, vjust=0),
      panel.border = element_blank(),
      axis.text.y = element_text(size=12),
      axis.text.x = element_text(angle = 45, hjust = 1, size=10),
      axis.title.x = element_blank(),
      axis.title.y = element_blank()) +
    eclectic::saturated_rainbow_pct()
  
  
  # Bar chart of max OTU
  
  E$o$prop <- apply(E$o$counts, 2, function(x) ifelse (x, x / sum(x),0))
  
  E$o$maxProp <- apply(E$o$prop, 2, max)
  E$o$maxProp <- E$o$maxProp[match(E$s$SampleID[E$sampleOrdering],
                                   names(E$o$maxProp))]
  E$maxPropBySpecies <-
    tapply(E$o$maxProp[match(E$s$SampleID, names(E$o$maxProp))],
           E$s$common,mean)
  E$maxPropBySpecies <-
    E$maxPropBySpecies[match(E$speciesTax$common[E$speciesOrdering],
                             names(E$maxPropBySpecies))]
  
  E$mpbsDf <- data.frame(common=names(E$maxPropBySpecies),value=unname(E$maxPropBySpecies))
  E$mpDf <- data.frame(common=E$s$common[match(names(E$o$maxProp),E$s$SampleID)],
                       value=E$o$maxProp)
  plotMax <- ggplot(E$mpbsDf,aes(x=common,y=value,fill="cyan")) +
    geom_bar(stat="identity",color="cyan") +
    scale_fill_manual(values=c("cyan")) + guides(fill=FALSE) +
    scale_x_discrete(limits=E$speciesTax$common[E$speciesOrdering]) +
    coord_flip() + geom_col(show.legend=FALSE) +
    theme_bw() + labs(y="Most dominant OTU") +
    theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank(),
          panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
    geom_point(data=E$mpDf,aes(x=common,y=value),color="black",shape="+",size=4,show.legend = FALSE)
  
  
  # Bar chart of shannon diversity
  
  E$otuShannon <- shannon(E$o$counts)
  E$otuShannonDf <- data.frame(SampleID=as.character(E$otuShannon[,1]),
                               shannon=as.numeric(E$otuShannon[,2]))
  E$otuShannonDf$common <- E$s$common[match(E$s$SampleID,E$otuShannonDf$SampleID)]
  E$otuShannonBySpecies <- aggregate(shannon ~ common, data=E$otuShannonDf, mean)
  E$otuShannonBySpecies <-
    E$otuShannonBySpecies[match(E$speciesTax$common[E$speciesOrdering],
                                E$otuShannonBySpecies$common),]
  
  plotShannon <- ggplot(E$otuShannonBySpecies,aes(x=common,y=shannon, fill="purple")) +
    geom_bar(stat="identity",color="purple") +
    scale_fill_manual(values=c("purple")) + guides(fill=FALSE) +
    scale_x_discrete(limits=E$speciesTax$common[E$speciesOrdering]) +
    coord_flip() + geom_col(show.legend=FALSE) +
    theme_bw() + labs(y="Shannon Diversity") +
    theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank(),
          panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
    geom_point(data=E$otuShannonDf,aes(x=common,y=shannon),color="black",shape="+",size=4,show.legend = FALSE)
  
  
  # Bar chart of percent reads unassigned
  
  E$unassignedOtuIds <- names(E$a)[E$a=="Unassigned"]
  E$percentTotalUnassignedBySample <-
    apply(E$o$counts,2,function(X) sum(X[names(X) %in% E$unassignedOtuIds]) / sum(X))
  E$percentTotalUnassignedBySample <-
    E$percentTotalUnassignedBySample[match(E$s$SampleID[E$sampleOrdering],
                                           names(E$percentTotalUnassignedBySample))]
  E$s$percentTotalUnassigned <-
    E$percentTotalUnassignedBySample[match(E$s$SampleID, names(E$percentTotalUnassignedBySample))]
  
  
  E$percentTotalUnassignedBySpecies <-
    tapply(E$percentTotalUnassignedBySample[match(E$s$SampleID, names(E$percentTotalUnassignedBySample))],
           E$s$common,mean)
  E$percentTotalUnassignedBySpecies <-
    E$percentTotalUnassignedBySpecies[match(E$speciesTax$common[E$speciesOrdering],
                                            names(E$percentTotalUnassignedBySpecies))]
  
  E$ptubsDf <- data.frame(common=names(E$percentTotalUnassignedBySpecies),
                          value=unname(E$percentTotalUnassignedBySpecies))
  E$ptuDf <- data.frame(common=E$s$common[match(names(E$percentTotalUnassignedBySample),E$s$SampleID)],
                        value=E$percentTotalUnassignedBySample)
  
  plotUnassigned <- ggplot(E$ptubsDf,aes(x=common,y=value,fill="red")) +
    geom_bar(stat="identity",color="red") +
    scale_fill_manual(values=c("red")) + guides(fill=FALSE) +
    scale_x_discrete(limits=E$speciesTax$common[E$speciesOrdering]) +
    coord_flip() + geom_col(show.legend=FALSE) +
    theme_bw() + labs(y="% OTUs unassigned") +
    theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank(),
          panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
    geom_point(data=E$ptuDf,aes(x=common,y=value),color="black",shape="+",size=4,show.legend = FALSE)
  
  
  
  # Barchart of Wolbachia OTUs
  
  E$wolbachiaOtus <- names(E$a[grep("Wolbachia",E$a)])
  E$percentWolbachiaBySample <-
    apply(E$o$counts,2,function(X) sum(X[names(X) %in% E$wolbachiaOtus]) / sum(X))
  E$percentWolbachiaBySample <-
    E$percentWolbachiaBySample[match(E$s$SampleID[E$sampleOrdering],
                                     names(E$percentWolbachiaBySample))]
  
  E$percentWolbachiaBySpecies <-
    tapply(E$percentWolbachiaBySample[match(E$s$SampleID, names(E$percentWolbachiaBySample))],
           E$s$common,mean)
  E$percentWolbachiaBySpecies <-
    E$percentWolbachiaBySpecies[match(E$speciesTax$common[E$speciesOrdering],
                                      names(E$percentWolbachiaBySpecies))]
  
  E$pwbsDf <- data.frame(common=names(E$percentWolbachiaBySpecies),
                         value=unname(E$percentWolbachiaBySpecies))
  E$pwDf <- data.frame(common=E$s$common[match(names(E$percentWolbachiaBySample),E$s$SampleID)],
                       value=E$percentWolbachiaBySample)
  plotWolbachia <- ggplot(E$pwbsDf,aes(x=common,y=value, fill="orange")) +
    geom_bar(stat="identity",color="orange") +
    scale_fill_manual(values=c("orange")) + guides(fill=FALSE) +
    scale_x_discrete(limits=E$speciesTax$common[E$speciesOrdering]) +
    coord_flip() + geom_col(show.legend=FALSE) +
    theme_bw() + labs(y="% OTUs in Wolbachia") +
    theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank(),
          panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
    geom_point(data=E$pwDf,aes(x=common,y=value),color="black",shape="+",size=4,show.legend = FALSE)
  
  
  
  # Barcharts of % reads unique to each tax level
  
  E$oByHostSpecies <- t(rowsum(t(E$o$counts),E$s$common[match(colnames(E$o$counts),E$s$SampleID)]))
  E$otuNumSpecies <- rowSums(E$oByHostSpecies>0)
  E$oUniqueToSpecies <- E$o$counts[E$otuNumSpecies==1,]
  
  E$uniqueOtusBySpecies <- E$oByHostSpecies[E$otuNumSpecies==1,]
  E$s$uniqueOtusToSpecies <- colSums(E$oUniqueToSpecies > 0)[match(colnames(E$oUniqueToSpecies),
                                                                   E$s$SampleID)]
  E$s$propUniqueOtusToSpecies <- E$s$uniqueOtusToSpecies / E$s$numOtus
  E$meanPropUniqueOtusToSpecies <- tapply(E$s$propUniqueOtusToSpecies,E$s$common,mean)
  E$meanPropUniqueOtusToSpecies <-
    E$meanPropUniqueOtusToSpecies[match(E$speciesTax$common[E$speciesOrdering],
                                        names(E$meanPropUniqueOtusToSpecies))]
  E$s$uniqueReadsToSpecies <- colSums(E$oUniqueToSpecies)[match(colnames(E$oUniqueToSpecies),
                                                                E$s$SampleID)]
  E$s$propUniqueReadsToSpecies <- E$s$uniqueReadsToSpecies / E$s$filteredReadCount
  E$meanPropUniqueReadsToSpecies <- tapply(E$s$propUniqueReadsToSpecies,E$s$common,mean)
  E$meanPropUniqueReadsToSpecies <-
    E$meanPropUniqueReadsToSpecies[match(E$speciesTax$common[E$speciesOrdering],
                                         names(E$meanPropUniqueReadsToSpecies))]
  
  E$pubsDf <- data.frame(common=names(E$meanPropUniqueReadsToSpecies),
                         value=unname(E$meanPropUniqueReadsToSpecies))
  E$puDf <- data.frame(common=E$s$common,
                       value=E$s$propUniqueReadsToSpecies)
  
  plotUniqueToSpecies <- ggplot(E$pubsDf,aes(x=common,y=value,fill="green")) +
    geom_bar(stat="identity",color="green") +
    scale_fill_manual(values=c("green")) + guides(fill=FALSE) +
    scale_x_discrete(limits=E$speciesTax$common[E$speciesOrdering]) +
    coord_flip() + geom_col(show.legend=FALSE) +
    theme_bw() + labs(y="% reads in OTUs unique to species") +
    theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank(),
          panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
    geom_point(data=E$puDf,aes(x=common,y=value),color="black",shape="+",size=4,show.legend = FALSE)
  
  
  
  vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
  
  pdf(paste0(paste0(work_dir, "R_plots/figure1.",study,".pdf")),
      height=15, width=20, onefile = FALSE)
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(1000, 13)))
  print(plotHeatmap, vp = vplayout(35:1000, 5:8))
  print(plotMax, vp = vplayout(33:919, 9))
  print(plotShannon, vp = vplayout(33:921, 10))
  print(plotUnassigned, vp = vplayout(33:921, 11))
  print(plotUniqueToSpecies, vp = vplayout(33:921, 12))
  print(plotWolbachia, vp = vplayout(33:921, 13))
  
  dev.off()
  
}

plotTsneLegends <- function() {
  legend(
    par('usr')[2]+.01*diff(par('usr')[1:2]), # x pos a bit to right of right border
    par('usr')[4], #y pos in middle of plot
    names(G$classColors), col=G$classColors, pch=1, # labels and annotations here
    xjust=0,xpd=NA #left justify and plot outside region
  )
  legend(
    par('usr')[2]+.01*diff(par('usr')[1:2]), # x pos a bit to right of right border
    par('usr')[4]-0.4*diff(par('usr')[3:4]), #y pos in middle of plot
    names(G$mammalOrderShapeLegend), pch=as.numeric(G$mammalOrderShapeLegend),
    col=G$classColors['Mammalia'],
    xjust=0,xpd=NA #left justify and plot outside region
  )
  legend(
    par('usr')[2]+.01*diff(par('usr')[1:2]), # x pos a bit to right of right border
    par('usr')[4]-0.7*diff(par('usr')[3:4]), #y pos in middle of plot
    names(G$insectOrderShapeLegend), pch=as.numeric(G$insectOrderShapeLegend),
    col=G$classColors['Insecta'],
    xjust=0,xpd=NA #left justify and plot outside region
  )
}

countsToProp <- function(counts) {
  return(apply(counts, 2, function(x)   ifelse (x, x / sum(x),0)))
}



# Export OTU tables----
# NEEDS UPDATING----

# All samples
resetData()
write.table(cbind(G$o$metadata[rownames(G$o$counts)],G$o$counts),
            paste0(work_dir, "toShare/otuCount.all.csv"),
            row.names=FALSE,sep=",",quote = FALSE)

# Remove singleton OTUs
G$o$counts2 <- G$o$counts[rowSums(G$o$counts)>1,]
write.table(cbind(G$o$metadata[rownames(G$o$counts2)],G$o$counts2),
            paste0(work_dir, "toShare/otuCount.all.noSingletons.csv"),
            row.names=FALSE,sep=",",quote = FALSE)

# Remove controls and failed samples
resetData("toUse")
write.table(cbind(G$o$metadata[rownames(G$o$counts)],G$o$counts),
            paste0(work_dir, "toShare/otuCount.final.csv"),
            row.names=FALSE,sep=",",quote = FALSE)

# Remove OTUs with < 10 reads
G$o$counts2 <- G$o$counts[rowSums(G$o$counts)>=10,]
write.table(cbind(G$o$metadata[rownames(G$o$counts2)],G$o$counts2),
            paste0(work_dir, "toShare/otuCount.final.min10.csv"),
            sep=",",quote = FALSE)

# Aggregate on taxonomy
G$o$countsWithTaxa <- cbind(G$o$metadata[rownames(G$o$counts)],G$o$counts)
colnames(G$o$countsWithTaxa)[1] <- "taxonomy"
G$o$countsByTaxa <- aggregate(G$o$counts,by=list(G$o$metadata[rownames(G$o$counts)]),FUN=sum)
colnames(G$o$countsByTaxa)[1] <- "Taxonomy"
G$o$countsByTaxa$numOtus <-sapply(G$o$countsByTaxa$Taxonomy,
                                  function(X){ sum(G$o$metadata==X) })
G$o$countsByTaxa$otuPhylum <- sub("; c__(.)*$", "", G$o$countsByTaxa$Taxonomy, perl=T)
G$o$countsByTaxa$numReads <- sapply(G$o$countsByTaxa$Taxonomy,
                                    function(X) { sum(G$o$counts[G$o$metadata==X,]) })
G$o$phylaData <- aggregate(G$o$countsByTaxa[,c("numOtus","numReads")],
                           by=list(G$o$countsByTaxa$otuPhylum),sum)

resetData("toUsePlusControls")
G$o$prop <- apply(G$o$counts, 2, function(x) ifelse (x, x / sum(x),0))
G$maxOtuPropAll <- apply(G$o$prop, 1, max)
G$finalOtuTable <- G$o$counts[names(sort(G$maxOtuPropAll, decreasing=T))[1:10000],]
write.table(cbind(G$o$metadata[rownames(G$finalOtuTable)],G$finalOtuTable),
            paste0(work_dir, "toShare/otuTable.final.csv"),
            row.names=T,sep=",",quote = FALSE)




# Run Unifrac----
resetData("toUse")
G$o$prop <- countsToProp(G$o$counts)
write.table(cbind(G$o$prop,G$o$metadata),paste0(work_dir,"otu/otu_prop.txt"), quote = F, sep = "\t")
system(paste0("sed -i '1s/^/# Constructed from biom file\\n#OTU ID\t/' ",work_dir,"otu/otu_prop.txt"))
system(paste0("sed -i '2s/$/Consensus Lineage/' ",work_dir,"otu/otu_prop.txt"))
system(paste0("biom convert -i ",work_dir,"otu/otu_prop.txt -o ",work_dir,
              "otu/otu_prop.biom --table-type=\"OTU table\" --to-hdf5"))
system(paste0("beta_diversity.py -i ",work_dir,"otu/otu_prop.biom -o ",work_dir,
              "bd_wnu -t ",work_dir,"otu/rep_set.tre -m weighted_unifrac"))

G$o$rare <- apply(G$o$counts,2,function(x) rarefyCounts(x,1000))
rownames(G$o$rare) <- rownames(G$o$counts)
write.table(cbind(G$o$rare,G$o$metadata),paste0(work_dir,"otu/otu_rare.txt"), quote = F, sep = "\t")
system(paste0("sed -i '1s/^/# Constructed from biom file\\n#OTU ID\t/' ",work_dir,"otu/otu_rare.txt"))
system(paste0("sed -i '2s/$/Consensus Lineage/' ",work_dir,"otu/otu_rare.txt"))
system(paste0("biom convert -i ",work_dir,"otu/otu_rare.txt -o ",work_dir,
              "otu/otu_rare.biom --table-type=\"OTU table\" --to-hdf5"))
system(paste0("beta_diversity.py -i ",work_dir,"otu/otu_rare.biom -o ",work_dir,
              "bd_uur -t ",work_dir,"otu/rep_set.tre -m weighted_unifrac"))


# Figure 1----

resetData("toUse")

library(reshape2)
#plot(newick)
#plotTree(newick,node.numbers=TRUE)
#plot(newick,show.node.label=TRUE)
#edgelabels(round(newick$edge.length,0), font=1)
#plot_tree(newick,"treeonly",nodeplotblank,base.spacing=0.55)

newickPlot <- plot_tree(G$newick,"treeonly",nodeplotblank,base.spacing=0.55)

# Heatmap of OTUs by phylum grouped by species

G$aPhylum <- sub("; c__(.)*$", "", G$o$metadata, perl=T)
G$oByPhylum <- rowsum(G$o$counts,G$aPhylum)

G$obpProp <- apply(G$oByPhylum,2,function(x) ifelse (x,x / sum(x),0))
G$obpProp <- data.frame(t(G$obpProp[rowSums(G$obpProp)>0,]))

G$obpProp$common <- as.character(G$s$common[match(rownames(G$obpProp),G$s$SampleID)])
G$obpPropBySpecies <- aggregate(.~common,data=G$obpProp,mean)

G$obpMaxProps <- apply(G$obpPropBySpecies[,-1],2,max)

numOtuPhylaToUse <- 10
G$otuPhylaToUse <- names(sort(G$obpMaxProps,decreasing=TRUE))[1:(numOtuPhylaToUse-1)]

G$obpPropBySpeciesToUse <- G$obpPropBySpecies[,names(G$obpPropBySpecies) == 'common' |
                                                names(G$obpPropBySpecies) %in% G$otuPhylaToUse]

colnames(G$obpPropBySpeciesToUse)[-1] <-
  colnames(G$obpPropBySpeciesToUse)[rev(order(colSums(G$obpPropBySpeciesToUse[-1])))+1]
G$obpPropBySpeciesToUse[,-1] <-
  G$obpPropBySpeciesToUse[,rev(order(colSums(G$obpPropBySpeciesToUse[-1])))+1]

G$obpPropBySpeciesToUse$Other <- 1 - rowSums(G$obpPropBySpeciesToUse[-1])

G$obpHeatmapData <- melt(cbind(G$obpPropBySpeciesToUse), id.vars=c('common'))

plotHeatmap <- ggplot(G$obpHeatmapData, aes(x=variable, y=common, fill=value)) +
  scale_y_discrete(limits=G$speciesTax$common[G$figure1SpeciesOrdering]) +
  geom_tile(color="grey80", size=0.4) + theme_grey() +
  theme(
    legend.position = "none",
    strip.text.y = element_text(angle=0, vjust=0),
    strip.text.x = element_text(angle=90, vjust=0),
    panel.border = element_blank(),
    axis.text.y = element_text(size=12),
    axis.text.x = element_text(angle = 45, hjust = 1, size=10),
    axis.title.x = element_blank(),
    axis.title.y = element_blank()) +
  eclectic::saturated_rainbow_pct()


# Bar chart of max OTU

G$o$prop <- apply(G$o$counts, 2, function(x) ifelse (x, x / sum(x),0))

G$o$maxProp <- apply(G$o$prop, 2, max)
G$o$maxProp <- G$o$maxProp[match(G$s$SampleID[G$sampleOrdering],
                                 names(G$o$maxProp))]
G$maxPropBySpecies <-
  tapply(G$o$maxProp[match(G$s$SampleID, names(G$o$maxProp))],
         G$s$common,mean)
G$maxPropBySpecies <-
  G$maxPropBySpecies[match(G$speciesTax$common[G$figure1SpeciesOrdering],
                           names(G$maxPropBySpecies))]

G$mpbsDf <- data.frame(common=names(G$maxPropBySpecies),value=unname(G$maxPropBySpecies))
G$mpDf <- data.frame(common=G$s$common[match(names(G$o$maxProp),G$s$SampleID)],
                     value=G$o$maxProp)
plotMax <- ggplot(G$mpbsDf,aes(x=common,y=value,fill="cyan")) +
  geom_bar(stat="identity",color="cyan") +
  scale_fill_manual(values=c("cyan")) + guides(fill=FALSE) +
  scale_x_discrete(limits=G$speciesTax$common[G$figure1SpeciesOrdering]) +
  coord_flip() + geom_col(show.legend=FALSE) +
  theme_bw() + labs(y="Most dominant OTU") +
  theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank(),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  geom_point(data=G$mpDf,aes(x=common,y=value),color="black",shape="+",size=4,show.legend = FALSE)


# Bar chart of shannon diversity

G$otuShannon <- shannon(G$o$counts)
G$otuShannonDf <- data.frame(SampleID=as.character(G$otuShannon[,1]),
                             shannon=as.numeric(G$otuShannon[,2]))
G$otuShannonDf$common <- G$s$common[match(G$s$SampleID,G$otuShannonDf$SampleID)]
G$otuShannonBySpecies <- aggregate(shannon ~ common, data=G$otuShannonDf, mean)
G$otuShannonBySpecies <-
  G$otuShannonBySpecies[match(G$speciesTax$common[G$figure1SpeciesOrdering],
                              G$otuShannonBySpecies$common),]

plotShannon <- ggplot(G$otuShannonBySpecies,aes(x=common,y=shannon, fill="purple")) +
  geom_bar(stat="identity",color="purple") +
  scale_fill_manual(values=c("purple")) + guides(fill=FALSE) +
  scale_x_discrete(limits=G$speciesTax$common[G$figure1SpeciesOrdering]) +
  coord_flip() + geom_col(show.legend=FALSE) +
  theme_bw() + labs(y="Shannon Diversity") +
  theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank(),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  geom_point(data=G$otuShannonDf,aes(x=common,y=shannon),color="black",shape="+",size=4,show.legend = FALSE)


# Bar chart of percent reads unassigned

G$unassignedOtuIds <- names(G$a)[G$a=="Unassigned"]
G$percentTotalUnassignedBySample <-
  apply(G$o$counts,2,function(X) sum(X[names(X) %in% G$unassignedOtuIds]) / sum(X))
G$percentTotalUnassignedBySample <-
  G$percentTotalUnassignedBySample[match(G$s$SampleID[G$sampleOrdering],
                                         names(G$percentTotalUnassignedBySample))]
G$s$percentTotalUnassigned <-
  G$percentTotalUnassignedBySample[match(G$s$SampleID, names(G$percentTotalUnassignedBySample))]


G$percentTotalUnassignedBySpecies <-
  tapply(G$percentTotalUnassignedBySample[match(G$s$SampleID, names(G$percentTotalUnassignedBySample))],
         G$s$common,mean)
G$percentTotalUnassignedBySpecies <-
  G$percentTotalUnassignedBySpecies[match(G$speciesTax$common[G$figure1SpeciesOrdering],
                                          names(G$percentTotalUnassignedBySpecies))]

G$ptubsDf <- data.frame(common=names(G$percentTotalUnassignedBySpecies),
                        value=unname(G$percentTotalUnassignedBySpecies))
G$ptuDf <- data.frame(common=G$s$common[match(names(G$percentTotalUnassignedBySample),G$s$SampleID)],
                      value=G$percentTotalUnassignedBySample)

plotUnassigned <- ggplot(G$ptubsDf,aes(x=common,y=value,fill="red")) +
  geom_bar(stat="identity",color="red") +
  scale_fill_manual(values=c("red")) + guides(fill=FALSE) +
  scale_x_discrete(limits=G$speciesTax$common[G$figure1SpeciesOrdering]) +
  coord_flip() + geom_col(show.legend=FALSE) +
  theme_bw() + labs(y="% OTUs unassigned") +
  theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank(),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  geom_point(data=G$ptuDf,aes(x=common,y=value),color="black",shape="+",size=4,show.legend = FALSE)


# Barchart of Wolbachia OTUs

G$wolbachiaOtus <- names(G$a[grep("Wolbachia",G$a)])
G$percentWolbachiaBySample <-
  apply(G$o$counts,2,function(X) sum(X[names(X) %in% G$wolbachiaOtus]) / sum(X))
G$percentWolbachiaBySample <-
  G$percentWolbachiaBySample[match(G$s$SampleID[G$sampleOrdering],
                                   names(G$percentWolbachiaBySample))]

G$percentWolbachiaBySpecies <-
  tapply(G$percentWolbachiaBySample[match(G$s$SampleID, names(G$percentWolbachiaBySample))],
         G$s$common,mean)
G$percentWolbachiaBySpecies <-
  G$percentWolbachiaBySpecies[match(G$speciesTax$common[G$figure1SpeciesOrdering],
                                          names(G$percentWolbachiaBySpecies))]

G$pwbsDf <- data.frame(common=names(G$percentWolbachiaBySpecies),
                       value=unname(G$percentWolbachiaBySpecies))
G$pwDf <- data.frame(common=G$s$common[match(names(G$percentWolbachiaBySample),G$s$SampleID)],
                     value=G$percentWolbachiaBySample)
plotWolbachia <- ggplot(G$pwbsDf,aes(x=common,y=value, fill="orange")) +
  geom_bar(stat="identity",color="orange") +
  scale_fill_manual(values=c("orange")) + guides(fill=FALSE) +
  scale_x_discrete(limits=G$speciesTax$common[G$figure1SpeciesOrdering]) +
  coord_flip() + geom_col(show.legend=FALSE) +
  theme_bw() + labs(y="% OTUs in Wolbachia") +
  theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank(),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  geom_point(data=G$pwDf,aes(x=common,y=value),color="black",shape="+",size=4,show.legend = FALSE)


# Barcharts of % reads unique to each tax level

G$oByHostSpecies <- t(rowsum(t(G$o$counts),G$s$common[match(colnames(G$o$counts),G$s$SampleID)]))
G$otuNumSpecies <- rowSums(G$oByHostSpecies>0)
G$oUniqueToSpecies <- G$o$counts[G$otuNumSpecies==1,]

G$uniqueOtusBySpecies <- G$oByHostSpecies[G$otuNumSpecies==1,]
G$s$uniqueOtusToSpecies <- colSums(G$oUniqueToSpecies > 0)[match(colnames(G$oUniqueToSpecies),
                                                                 G$s$SampleID)]
G$s$propUniqueOtusToSpecies <- G$s$uniqueOtusToSpecies / G$s$numOtus
G$meanPropUniqueOtusToSpecies <- tapply(G$s$propUniqueOtusToSpecies,G$s$common,mean)
G$meanPropUniqueOtusToSpecies <-
  G$meanPropUniqueOtusToSpecies[match(G$speciesTax$common[G$figure1SpeciesOrdering],
                                       names(G$meanPropUniqueOtusToSpecies))]
G$s$uniqueReadsToSpecies <- colSums(G$oUniqueToSpecies)[match(colnames(G$oUniqueToSpecies),
                                                                      G$s$SampleID)]
G$s$propUniqueReadsToSpecies <- G$s$uniqueReadsToSpecies / G$s$filteredReadCount
G$meanPropUniqueReadsToSpecies <- tapply(G$s$propUniqueReadsToSpecies,G$s$common,mean)
G$meanPropUniqueReadsToSpecies <-
  G$meanPropUniqueReadsToSpecies[match(G$speciesTax$common[G$figure1SpeciesOrdering],
                                  names(G$meanPropUniqueReadsToSpecies))]

G$pubsDf <- data.frame(common=names(G$meanPropUniqueReadsToSpecies),
                       value=unname(G$meanPropUniqueReadsToSpecies))
G$puDf <- data.frame(common=G$s$common,
                     value=G$s$propUniqueReadsToSpecies)

plotUniqueToSpecies <- ggplot(G$pubsDf,aes(x=common,y=value,fill="green")) +
  geom_bar(stat="identity",color="green") +
  scale_fill_manual(values=c("green")) + guides(fill=FALSE) +
  scale_x_discrete(limits=G$speciesTax$common[G$figure1SpeciesOrdering]) +
  coord_flip() + geom_col(show.legend=FALSE) +
  theme_bw() + labs(y="% reads in OTUs unique to species") +
  theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank(),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  geom_point(data=G$puDf,aes(x=common,y=value),color="black",shape="+",size=4,show.legend = FALSE)


# Barchart of obligate anaerobes

G$anaerobeTable <-read.csv(paste0(work_dir,"bilaterian_aerotolerance.csv"),header=T,
                           na.strings="<NA>",stringsAsFactors = F)
G$anaerobeTable$SampleID[G$anaerobeTable$SampleID=="Bonobo.BI0331.Feces.1.1"] <- "Pig.BI0331.Feces.1.1"
G$anaerobeTable$speciesOrder <- G$s[G$anaerobeTable$SampleID,"speciesOrder"]
G$anaerobeTable <-
  G$anaerobeTable[!G$anaerobeTable$aerobic_status_corrected=="NA",]
G$anaerobeTable <-
  G$anaerobeTable[!G$anaerobeTable$aerobic_status_corrected=="no taxonomic assignment",]
G$anaerobeTable <-
  G$anaerobeTable[!G$anaerobeTable$aerobic_status_corrected=="Streptophyta",]
G$anaerobeTable$Proportion <-
  ave(G$anaerobeTable$Counts,G$anaerobeTable$SampleID,
      FUN=function(X) X / sum(X))


G$obligateAnaerobeTable <-
  G$anaerobeTable[G$anaerobeTable$aerobic_status_corrected=="obligate anaerobe",]
G$s$obligateAnaerobeProp <- 0
G$s$obligateAnaerobeProp[match(G$obligateAnaerobeTable$SampleID,G$s$SampleID)] <-
  G$obligateAnaerobeTable$Proportion
G$obligateAnaerobeSpecies <- tapply(G$s$obligateAnaerobeProp,G$s$common,mean)
G$obligateAnaerobeSpecies <-
  G$obligateAnaerobeSpecies[match(G$speciesTax$common[G$figure1SpeciesOrdering],
                                  names(G$obligateAnaerobeSpecies))]

G$oapbsDf <- data.frame(common=names(G$obligateAnaerobeSpecies),
                        value=unname(G$obligateAnaerobeSpecies))
G$oapDf <- data.frame(common=G$s$common,
                      value=G$s$obligateAnaerobeProp)

plotObligateAnaerobes <- ggplot(G$oapbsDf,aes(x=common,y=value,fill="pink")) +
  geom_bar(stat="identity",color="pink") +
  scale_fill_manual(values=c("pink")) + guides(fill=FALSE) +
  scale_x_discrete(limits=G$speciesTax$common[G$figure1SpeciesOrdering]) +
  coord_flip() + geom_col(show.legend=FALSE) +
  theme_bw() + labs(y="% reads assigned to obligate anaerobes") +
  theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank(),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  geom_point(data=G$oapDf,aes(x=common,y=value),color="black",shape="+",size=4,show.legend = FALSE)


# Print out figure

#save.image(paste0(work_dir, "R//good.20171016.RData")
#load(paste0(work_dir, "R//good.20171016.RData")

vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)

pdf(paste0(work_dir, "R_plots/figure1.final.pdf"),
    height=15, width=20, onefile = FALSE)
library(grid)
grid.newpage()
pushViewport(viewport(layout = grid.layout(1000, 14)))
print(newickPlot, vp = vplayout(1:935, 1:2))
print(plotHeatmap, vp = vplayout(35:1000, 5:8))
print(plotMax, vp = vplayout(33:919, 9))
print(plotShannon, vp = vplayout(33:921, 10))
print(plotUnassigned, vp = vplayout(33:921, 11))
print(plotUniqueToSpecies, vp = vplayout(33:921, 12))
print(plotWolbachia, vp = vplayout(33:921, 13))
print(plotObligateAnaerobes, vp = vplayout(33:921, 14))
dev.off()


# Figure 2 (tSNE)----
# PCoA, tsne, heatmaps for final figures

resetData("toUse")
resetUnifracData()

set.seed(1001)
library(Rtsne)
G$uur_tsne <- Rtsne(G$uur,is_distance=TRUE,verbose=TRUE,perplexity=10,max_iter=3000)

G$classColors <- makeLegendColors(G$s$class)
G$s$classColor <- G$classColors[G$s$class]

G$classOrderNumber <-
  unlist(sapply(unique(G$s$class),function(X) makeLegendShapes(G$s$order[G$s$class==X])))

G$s$orderShape <- as.numeric(G$classOrderNumber[as.character(G$s$order)])

G$mammalOrderShapeLegend <-
  G$classOrderNumber[names(G$classOrderNumber) %in% G$s$order[G$s$class=='Mammalia']]
G$mammalOrderShapeLegend <- G$mammalOrderShapeLegend[order(G$mammalOrderShapeLegend)]

G$insectOrderShapeLegend <-
  G$classOrderNumber[names(G$classOrderNumber) %in% G$s$order[G$s$class=='Insecta']]
G$insectOrderShapeLegend <- G$insectOrderShapeLegend[order(G$insectOrderShapeLegend)]

pdf(paste0(work_dir, "R_plots/figure2.tsne.uur.pdf"),width=11,height=8)
par(mar=c(4,4,1,12)) # need to leave margin room
plot(G$uur_tsne$Y, col=G$s$classColor, pch=as.numeric(G$s$orderShape),
     xlab="TSNE1",ylab="TSNE2", main="t-SNE plot for 266 samples")
plotTsneLegends()
dev.off()




# Figure 1G/H----
# plotDNA for unassigned vs assigned reads, both 5' and 3'

library(dnaplotr)
library(dnar)
library(ape)

GH <- new.env()

GH$otuTree <- read.tree(paste0(work_dir, "otu_good/rep_set.tre"))

GH$repSetAligned <-
  read.fa(paste0(work_dir, "otu_good/pynast_aligned_seqs/goodSeqs_rep_set_aligned.fasta"))
rownames(GH$repSetAligned) <- sapply(GH$repSetAligned$name,function(X) strsplit(X," ")[[1]][1])

GH$repSetAA <- GH$repSetAligned[!rownames(GH$repSetAligned) %in% names(GH$a)[GH$a=="Unassigned"],c("name","seq")]

GH$repSetAAordered <-
  GH$repSetAA[match(GH$otuTree$tip.label[GH$otuTree$tip.label %in% rownames(GH$repSetAA)],
                         rownames(GH$repSetAA)),]

set.seed(1001)
GH$repSetAAsample <- GH$repSetAA[sample(nrow(GH$repSetAA),10000),]

GH$repSetAAsample1 <- GH$repSetAA[sample(nrow(GH$repSetAA),50),]

# code to identify alignment positions with < 10% gaps
# For each position (1-~7000), count # of gaps
# Return True if < 10% of positions are gaps, False otherwise

maxGaps <- 0.9 * nrow(GH$repSetAAsample1)
positionsToKeep <-
  sapply(1:nchar(GH$repSetAAsample1$seq[1]),
         function(X) { return(
           sum(sapply(GH$repSetAAsample1$seq,
                      function(Y) { unlist(strsplit(Y, ""))[X] }) == "-") < maxGaps) })

GH$repSetAAordered$seq2 <-
  sapply(GH$repSetAAordered$seq,
         function(X) { paste0(unlist(strsplit(X,""))[positionsToKeep],collapse="") })
#plotDNA(GH$repSetAAordered$seq2, res=4000)


GH$repSetAAsample$seq2 <-
  sapply(GH$repSetAAsample$seq,
         function(X) { paste0(unlist(strsplit(X,""))[positionsToKeep],collapse="") })

GH$repSetAAsampleOrdered <-
  GH$repSetAAsample[match(GH$otuTree$tip.label[GH$otuTree$tip.label %in% rownames(GH$repSetAAsample)],
                         rownames(GH$repSetAAsample)),]

#plotDNA(GH$repSetAAsample$seq2, res=4000)
plotDNA(GH$repSetAAsampleOrdered$seq2, res=4000)
#plotDNA(GH$repSetAAsampleOrdered$seq2[c(700:1100,5050:5600)], res=4000)

readIds <- sapply(GH$repSetAAsampleOrdered$name,function(X) paste(strsplit(X," ")[[1]][1]))
head(sort(table(GH$a[readIds]), decreasing = T))

head(sort(table(GH$a[readIds[700:1100]]), decreasing = T))
head(sort(table(GH$a[readIds[4900:5400]]), decreasing = T))

#sampleSeqNames <- sapply(GH$repSetAAsample$name, function(X) {paste0(">",paste(strsplit(X," ")[[1]][1:2],collapse=" "))})
#write(sampleSeqNames,file = paste0(work_dir, "repSetSample.txt")

GH$repSetAU <- GH$repSetAligned[rownames(GH$repSetAligned) %in% names(GH$a)[GH$a=="Unassigned"],c("name","seq")]
GH$repSetAUsample <- GH$repSetAU[sample(nrow(GH$repSetAU),10000),]

GH$repSetAUsample$seq2 <-
  sapply(GH$repSetAUsample$seq,
         function(X) { paste0(unlist(strsplit(X,""))[positionsToKeep],collapse="") })

GH$repSetAUsampleOrdered <-
  GH$repSetAUsample[match(GH$otuTree$tip.label[GH$otuTree$tip.label %in% rownames(GH$repSetAUsample)],
                         rownames(GH$repSetAUsample)),]


#save.image(paste0(work_dir, "R//good.20170830.RData")
load(paste0(work_dir, "R//good.20170830.RData"))

pdf(paste0(work_dir, "R_plots/figure1.panelG.pdf"),
    height=4, width=4, onefile = FALSE)
plotDNA(GH$repSetAAsample$seq2, res=4000)
plotDNA(GH$repSetAAsampleOrdered$seq2, res=4000)

dev.off()
pdf(paste0(work_dir, "R_plots/figure1.panelH.pdf"),
    height=4, width=4, onefile = FALSE)
plotDNA(GH$repSetAUsample$seq, res=4000)
dev.off()

GH$repSetAA <-
  GH$repSetAA[match(GH$otuTree$tip.label[GH$otuTree$tip.label %in% rownames(GH$repSetAA)],
                         rownames(GH$repSetAA)),]
GH$repSetAU <-
  GH$repSetAU[match(GH$otuTree$tip.label[GH$otuTree$tip.label %in% rownames(GH$repSetAU)],
                         rownames(GH$repSetAU)),]

pdf(paste0(work_dir, "R_plots/figure1.panelG.full.pdf"),
    height=4, width=4, onefile = FALSE)
plotDNA(GH$repSetAA$seq, res=4000)
dev.off()
pdf(paste0(work_dir, "R_plots/figure1.panelH.full.pdf"),
    height=4, width=4, onefile = FALSE)
plotDNA(GH$repSetAU$seq, res=4000)
dev.off()




# Figure A1----

resetData("toUse")

G$sampleMetadataDf <- do.call(rbind,
                              list(data.frame(level="Phylum",
                                              count=table(G$s$phylum[G$sampleOrdering])),
                                   data.frame(level="Class",
                                              count=table(G$s$class[G$sampleOrdering])),
                                   data.frame(level="Order",
                                              count=table(G$s$order[G$sampleOrdering])),
                                   data.frame(level="Family",
                                              count=table(G$s$family[G$sampleOrdering])),
                                   data.frame(level="Species",
                                              count=table(G$s$common[G$sampleOrdering]))))

G$sampleMetadataDf$count.Var1 <-
  factor(G$sampleMetadataDf$count.Var1,
         levels = c(unique(as.character(G$speciesTax$phylum[G$figure1SpeciesOrdering])),
                    unique(as.character(G$speciesTax$class[G$figure1SpeciesOrdering])),
                    unique(as.character(G$speciesTax$order[G$figure1SpeciesOrdering])),
                    unique(as.character(G$speciesTax$family[G$figure1SpeciesOrdering])),
                    unique(as.character(G$speciesTax$common[G$figure1SpeciesOrdering]))))
G$sampleMetadataDf <- G$sampleMetadataDf[match(levels(G$sampleMetadataDf$count.Var1),
                                               G$sampleMetadataDf$count.Var1),]
G$sampleMetadataDf <- G$sampleMetadataDf %>% group_by(level) %>%
  mutate(cum.freq = cumsum(count.Freq) - 0.5*count.Freq)
pal <- sapply(table(G$sampleMetadataDf$level),
              function(X) colorRampPalette(c('red','blue','dark green'))(X))
pal <- unname(unlist(pal))

palSet <- colorRampPalette(c('blue','dark green'))(nrow(G$s))
pal <- sapply(unique(G$sampleMetadataDf$level),
              function(X) palSet[G$sampleMetadataDf$cum.freq[G$sampleMetadataDf$level==X]+0.5])
pal <- unname(unlist(pal))

G$sampleMetadataDf$labelSize <- sapply(G$sampleMetadataDf$count.Freq, function(X) 2 * min(X,3))

pdf(paste0(work_dir, "R_plots/A1.final.pdf"),
    height=30, width=15, onefile = FALSE)
ggplot(G$sampleMetadataDf, aes(x = level, y = count.Freq, fill = count.Var1)) + 
  geom_bar(stat = "identity", show.legend=FALSE) + 
  scale_x_discrete(limits = c("Phylum","Class","Order","Family","Species"))  + 
  scale_fill_manual(values=pal) +
  theme_bw() + labs(y="Percent of OTUs in Wolbachia") +
  theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank(),
        panel.border = element_blank(), panel.grid.major = element_blank(), axis.title.x=element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_blank()) +
  geom_text(aes(label=paste0(count.Var1,": ", count.Freq),
                y=max(cum.freq)-cum.freq+0.5), colour="white",size=G$sampleMetadataDf$labelSize) +
  scale_y_continuous(labels=waiver(), trans="reverse")
dev.off()


# Figure A2----

resetData("all")

pdf(paste0(work_dir, "R_plots/A2.boxplot.readCounts.pdf"),width=11,height=8)
par(mar=c(5,7,6,6))
par(cex.lab=2)
par(cex.main=2)
boxplot(readCount ~ sample_type2, data = G$s, log = "y",
        xlab="Sample type", ylab="Number of reads",
        main="Number of reads per library by sample type")
dev.off()


# Figure A3----
# Reads by OTU Phylum

resetData("toUse")

#load(paste0(work_dir, "R//good.20170830.RData")

G$aPhylum <- sub("; c__(.)*$", "", G$o$metadata, perl=T)
G$oByPhylum <- rowsum(G$o$counts,G$aPhylum)

dataToPlot <- sort(rowSums(G$oByPhylum), decreasing = T)[1:9]
dataToPlot <- data.frame(value=dataToPlot,phylum=names(dataToPlot))

dataToPlot$phylum <- with(dataToPlot, factor(phylum, levels=c("Other",as.character(phylum[order(value)]))))
dataToPlot[10,] <- cbind(sum(sort(rowSums(G$oByPhylum), decreasing = T)[-c(1:9)]), "Other")

dataToPlot$value <- as.numeric(dataToPlot$value) / 1e6

pdf(paste0(work_dir, "R_plots/supp3.pdf"), width = 10, height = 7)
ggplot(dataToPlot, aes(x=phylum, y=value)) + geom_bar(stat="identity") +
  coord_flip() + theme_bw() +
  labs(y = "Number of reads (in millions)", x = "OTU Phylum") +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor = element_blank())
dev.off()





### Figures A4-A13 (Val.)----

VA <- new.env()
VA$otuTableFile <- paste0(work_dir, "validation/anteater/otu/otu_table.txt")
VA$mapFile <- paste0(work_dir, "validation/anteater/map/map.tsv")
resetValidationData(VA, "toUse")
plotFigure1(VA, "anteater")

VBi <- new.env()
VBi$otuTableFile <- paste0(work_dir, "validation/bird/otu/otu_table.txt")
VBi$mapFile <- paste0(work_dir, "validation/bird/map/map.tsv")
resetValidationData(VBi, "toUse")
plotFigure1(VBi, "bird")

VBu <- new.env()
VBu$otuTableFile <- paste0(work_dir, "validation/bug/otu/otu_table.txt")
VBu$mapFile <- paste0(work_dir, "validation/bug/map/map.tsv")
resetValidationData(VBu, "toUse")
plotFigure1(VBu, "bug")

VBR <- new.env()
VBR$otuTableFile <- paste0(work_dir, "validation/bugRev/otu/otu_table.txt")
VBR$mapFile <- paste0(work_dir, "validation/bugRev/map/map.tsv")
resetValidationData(VBR, "toUse")
plotFigure1(VBR, "bugRev")

VF <- new.env()
VF$otuTableFile <- paste0(work_dir, "validation/fish/otu/otu_table.txt")
VF$mapFile <- paste0(work_dir, "validation/fish/map/map.tsv")
resetValidationData(VF, "toUse")
plotFigure1(VF, "fish")

VM <- new.env()
VM$otuTableFile <- paste0(work_dir, "validation/muegge/otu/otu_table.txt")
VM$mapFile <- paste0(work_dir, "validation/muegge/map/map.tsv")
resetValidationData(VM, "toUse")
plotFigure1(VM, "muegge454")

VMI <- new.env()
VMI$otuTableFile <- paste0(work_dir, "validation/mueggeIllumina/otu/otu_table.txt")
VMI$mapFile <- paste0(work_dir, "validation/mueggeIllumina/map/map.tsv")
resetValidationData(VMI, "toUse")
plotFigure1(VMI, "mueggeIllumina")

VP <- new.env()
VP$otuTableFile <- paste0(work_dir, "validation/primates/otu/otu_table.txt")
VP$mapFile <- paste0(work_dir, "validation/primates/map/map.tsv")
resetValidationData(VP, "toUse")
plotFigure1(VP, "primates")

VW <- new.env()
VW$otuTableFile <- paste0(work_dir, "validation/whale454/otu/otu_table.txt")
VW$mapFile <- paste0(work_dir, "validation/whale454/map/map.tsv")
resetValidationData(VW, "toUse")
plotFigure1(VW, "whale454")

VB <- new.env()
VB$otuTableFile <- paste0(work_dir, "validation/whaleIllumina/otu/otu_table.txt")
VB$mapFile <- paste0(work_dir, "validation/whaleIllumina/map/map.tsv")
resetValidationData(VB, "toUse")
plotFigure1(VB, "whaleIllumina")







# Figure A14----
# Aerobic status of bacteria

G$anaerobeTable <-read.csv(paste0(work_dir,"bilaterian_aerotolerance.csv"),header=T,
                           na.strings="<NA>",stringsAsFactors = F)
G$anaerobeTable$speciesOrder <- G$s[G$anaerobeTable$SampleID,"speciesOrder"]
G$anaerobeTable <-
  G$anaerobeTable[!G$anaerobeTable$aerobic_status_corrected=="NA",]
G$anaerobeTable <-
  G$anaerobeTable[!G$anaerobeTable$aerobic_status_corrected=="no taxonomic assignment",]
G$anaerobeTable <-
  G$anaerobeTable[!G$anaerobeTable$aerobic_status_corrected=="Streptophyta",]
G$anaerobeTable$Proportion <-
  ave(G$anaerobeTable$Counts,G$anaerobeTable$SampleID,
      FUN=function(X) X / sum(X))

G$aerobeTable <-
  G$anaerobeTable[G$anaerobeTable$aerobic_status_corrected=="aerobe",]
G$s$aerobeProp <- 0
G$s$aerobeProp[match(G$aerobeTable$SampleID,G$s$SampleID)] <-
  G$aerobeTable$Proportion

G$obligateAnaerobeTable <-
  G$anaerobeTable[G$anaerobeTable$aerobic_status_corrected=="obligate anaerobe",]
G$s$obligateAnaerobeProp <- 0
G$s$obligateAnaerobeProp[match(G$obligateAnaerobeTable$SampleID,G$s$SampleID)] <-
  G$obligateAnaerobeTable$Proportion

G$s$aerobeToObligateAnaerobeRatio <-
  log10((G$s$aerobeProp + 0.00001) / (G$s$obligateAnaerobeProp + 0.00001))

G$anaerobeTableSpecies <- aggregate(Counts ~ common + aerobic_status_corrected + speciesOrder,
                                    G$anaerobeTable,sum)
G$anaerobeTableSpecies$Proportion <-
  ave(G$anaerobeTableSpecies$Counts,G$anaerobeTableSpecies$common,
      FUN=function(X) X / sum(X))

pdf(paste0(work_dir, "R_plots/aerobic_status_species.pdf"),
    height=11, width=8, onefile = FALSE)
G$anaerobeTableSpecies %>%
  mutate(ProportionNA = ifelse(is.na(aerobic_status_corrected), Proportion, 0)) %>%
  droplevels() %>%
  mutate(common = reorder(common, speciesOrder)) %>%
  ggplot() +
  geom_bar(aes(x=common, y=Proportion, fill=aerobic_status_corrected), stat="identity", position = "stack") +
  coord_flip() +
  scale_fill_brewer(palette="Paired", na.value="#CCCCCC") +
  theme_bw()
dev.off()


pdf(paste0(work_dir, "R_plots/aerobe_ratio.pdf"),
    height=8, width=11, onefile = FALSE)
plot(G$s$log.weight,G$s$aerobeToObligateAnaerobeRatio,
     col=G$s$classColor,xlab="Weight (grams, log10)",
     ylab="Ratio of aerobic reads to anaerobic reads",
     main="Ratio of aerobic/anaerobic reads vs sample weight.\nColored by host class.")
dev.off()

G$species$aerobeToObligateAnaerobeRatio <-
  with(aggregate(aerobeToObligateAnaerobeRatio ~ common,G$s, mean),
       aerobeToObligateAnaerobeRatio[match(G$species$common,common)])
G$species$log.weight <-
  with(aggregate(log.weight ~ common,G$s, mean),
       log.weight[match(G$species$common,common)])
G$classColors <- makeLegendColors(G$species$class)
G$species$classColor <- G$classColors[G$species$class]

pdf(paste0(work_dir, "R_plots/aerobe_ratio_species.pdf"),
    height=8, width=11, onefile = FALSE)
plot(G$species$log.weight,G$species$aerobeToObligateAnaerobeRatio,
     col=G$species$classColor,xlab="Weight (grams, log10)",
     ylab="Ratio of aerobic reads to anaerobic reads",
     main="Avg. ratio of aerobic/anaerobic reads vs avg. sample weight.\nAverages are by host species.\nColored by host class.")
dev.off()


### Marine bacteria----

resetData("toUse")

marineBacteriaStrings <- c("Synechococcales","Pelagibac","hodobacter")

isoTax <- function(metadata, tax, counts) {
  taxOtus <- names(metadata[grep(tax,metadata)])
  percentBySample <- apply(counts, 2, function(X) sum(X[names(X) %in% taxOtus]) / sum(X))
  return(percentBySample)
}

G$percentSynechocoBySample <- isoTax(G$o$metadata,"Synechococcales",G$o$counts)
G$percentSynechocoBySample <-
  G$percentSynechocoBySample[match(G$s$SampleID[G$sampleOrdering],
                                   names(G$percentSynechocoBySample))]


G$percentSynechocoBySpecies <-
  tapply(G$percentSynechocoBySample[match(G$s$SampleID, names(G$percentSynechocoBySample))],
         G$s$common,mean)
G$percentSynechocoBySpecies <-
  G$percentSynechocoBySpecies[match(G$speciesTax$common[G$figure1SpeciesOrdering],
                                    names(G$percentSynechocoBySpecies))]

G$percentPelagibacBySample <- isoTax(G$o$metadata,"Pelagibac",G$o$counts)
G$percentPelagibacBySample <-
  G$percentPelagibacBySample[match(G$s$SampleID[G$sampleOrdering],
                                   names(G$percentPelagibacBySample))]
G$percentPelagibacBySpecies <-
  tapply(G$percentPelagibacBySample[match(G$s$SampleID, names(G$percentPelagibacBySample))],
         G$s$common,mean)
G$percentPelagibacBySpecies <-
  G$percentPelagibacBySpecies[match(G$speciesTax$common[G$figure1SpeciesOrdering],
                                    names(G$percentPelagibacBySpecies))]

G$percentRoseobacterBySample <- isoTax(G$o$metadata,"hodobacter",G$o$counts)
G$percentRoseobacterBySample <-
  G$percentRoseobacterBySample[match(G$s$SampleID[G$sampleOrdering],
                                     names(G$percentRoseobacterBySample))]
G$percentRoseobacterBySpecies <-
  tapply(G$percentRoseobacterBySample[match(G$s$SampleID, names(G$percentRoseobacterBySample))],
         G$s$common,mean)
G$percentRoseobacterBySpecies <-
  G$percentRoseobacterBySpecies[match(G$speciesTax$common[G$figure1SpeciesOrdering],
                                      names(G$percentRoseobacterBySpecies))]






























# Adonis Unifrac stats----
# NEEDS UPDATING----

G$wnuMatrix <- as.matrix(G$wnu)
pheatmap(G$wnuMatrix, fontsize = 6)

G$wnuDist <- G$phyMatrix
G$wnuDist[,] <- 0
for(i in rownames(G$wnuDist)) {
  for(j in colnames(G$wnuDist)) {
    if(i==j) {
      G$wnuDist[i,j] <- 0
    } else {
      G$wnuDist[i,j] <- cd(G$wnuMatrix,
                           rownames(G$wnuMatrix)[rownames(G$wnuMatrix) %in% G$s$SampleID[G$s$common==i]],
                           rownames(G$wnuMatrix)[rownames(G$wnuMatrix) %in% G$s$SampleID[G$s$common==j]])
    }
  }
}

set.seed(1017)
G$wnuDist_tsne <- Rtsne(G$wnuDist,is_distance=TRUE,verbose=TRUE,perplexity=20,max_iter=3000)
pdf(paste0(work_dir, "R_plots/A16.wnu.speciesCentroid.tsne.pdf"),
    height=8, width=11, onefile = FALSE)
par(mar=c(4,4,1,12)) # need to leave margin room
plot(G$wnuDist_tsne$Y,col=G$speciesClassColors,xlab="",ylab="", cex=0.1,
     main="t-SNE clustering of species centroid of weighted Unifrac distances between normalized OTU counts")
text(G$wnuDist_tsne$Y[,1],G$wnuDist_tsne$Y[,2],rownames(G$wnuDist),col = G$speciesClassColors,
     cex=0.5)
dev.off()

pdf(paste0(work_dir, "R_plots/A16.wnu.speciesCentroid.tsne.points.pdf"),
    height=8, width=11, onefile = FALSE)
par(mar=c(4,4,1,12)) # need to leave margin room
plot(G$wnuDist_tsne$Y,col=G$speciesClassColors,xlab="",ylab="", cex=1,
     main="t-SNE clustering of species centroid of weighted Unifrac distances between normalized OTU counts")
dev.off()

#mantel.test(G$uurDistM,G$phyDist,nperm=1e6)

G$o$prop <- apply(G$o$counts, 2, function(x) ifelse (x, x / sum(x),0))
G$tree <- ape::multi2di(phyloseq::read_tree(paste0(work_dir, "otu_good/rep_set.tre")))
G$qiimeDataProp <- phyloseq(otu_table=otu_table(G$o$prop, taxa_are_rows = TRUE),phy_tree=G$tree)

G$o$rare <- apply(G$o$counts,2,function(x) rrarefy(x,1000))
rownames(G$o$rare) <- rownames(G$o$counts)
G$qiimeDataRare <- phyloseq(otu_table=otu_table(G$o$rare, taxa_are_rows = TRUE),phy_tree=G$tree)

G$uur2 <- UniFrac(G$qiimeDataRare,weighted=FALSE)
G$bcr2 <- distance(G$qiimeDataRare,'bray',binary=TRUE)
G$wnu2 <- UniFrac(G$qiimeDataProp,weighted=TRUE)
G$bcp2 <- distance(G$qiimeDataProp,'bray',binary=FALSE)

adonis(G$uur2)

G$mantels<-list(
  'uniW'=ade4::mantel.rtest(G$uniDist,G$uniDistW,nrepet=1e4),
  'brayW'=ade4::mantel.rtest(G$uniDist,G$brayDistW,nrepet=1e4),
  'brayUW'=ade4::mantel.rtest(G$uniDist,G$brayDist,nrepet=1e4)
)

adonis(G$uur2 ~ common,
       G$s[match(G$s$SampleID,attr(G$uur2,"Labels")),],
       permutations = nperm)

adonis(G$bcr2 ~ common,
       G$s[match(G$s$SampleID,attr(G$bcr2,"Labels")),],
       permutations = nperm)

adonis(G$wnu2 ~ common,
       G$s[match(G$s$SampleID,attr(G$wnu2,"Labels")),],
       permutations = nperm)

adonis(G$bcp2 ~ common,
       G$s[match(G$s$SampleID,attr(G$bcp2,"Labels")),],
       permutations = nperm)

G$s$classColor <- G$classColors[G$s$class]

set.seed(1001)
G$uur2_tsne <- Rtsne(G$uur2,is_distance=TRUE,verbose=TRUE,perplexity=20,max_iter=3000)

#pdf(paste0(work_dir, "R_plots/pcoa/wnu.speciesCentroid.tsne.pdf",
#    height=8, width=11, onefile = FALSE)
par(mar=c(4,4,1,12)) # need to leave margin room
plot(G$uur2_tsne$Y,col=G$s$classColor,xlab="",ylab="", cex=0.1,
     main="t-SNE clustering of species centroid of weighted Unifrac distances between normalized OTU counts")
text(G$uur2_tsne$Y[,1],G$uur2_tsne$Y[,2],labels(G$uur2),col = G$s$classColor,
     cex=0.5)
#dev.off()


nperm <- 1000000
nperm <- 1000

adonis(G$uur ~ phylum + class + order + family + genus + common,
       G$s[match(G$s$SampleID,attr(G$uur,"Labels")),],
       permutations = nperm)

adonis(G$uur ~ log.weight,
       G$s[match(G$s$SampleID,attr(G$uur,"Labels")),],
       permutations = nperm)

adonis(G$uur ~ Lab.Wild.Domestic,
       G$s[match(G$s$SampleID,attr(G$uur,"Labels")),],
       permutations = nperm)

adonis(G$uur ~ Lab.Wild.Domestic + phylum,
       G$s[match(G$s$SampleID,attr(G$uur,"Labels")),],
       permutations = nperm)

adonis(G$uur ~ phylum + Lab.Wild.Domestic,
       G$s[match(G$s$SampleID,attr(G$uur,"Labels")),],
       permutations = nperm)

adonis(G$uur ~ diet,
       G$s[match(G$s$SampleID,attr(G$uur,"Labels")),],
       permutations = nperm)

adonis(G$uur ~ phylum + diet,
       G$s[match(G$s$SampleID,attr(G$uur,"Labels")),],
       permutations = nperm)

adonis(G$uur ~ log.weight + diet,
       G$s[match(G$s$SampleID,attr(G$uur,"Labels")),],
       permutations = nperm)

adonis(G$uur ~ phylum + log.weight,
       G$s[match(G$s$SampleID,attr(G$uur,"Labels")),],
       permutations = nperm)

adonis(G$uur ~ class + log.weight,
       G$s[match(G$s$SampleID,attr(G$uur,"Labels")),],
       permutations = nperm)

adonis(G$uur ~ order + log.weight,
       G$s[match(G$s$SampleID,attr(G$uur,"Labels")),],
       permutations = nperm)

adonis(G$uur ~ Lab.Wild.Domestic + diet,
       G$s[match(G$s$SampleID,attr(G$uur,"Labels")),],
       permutations = nperm)

adonis(G$uur ~ diet + phylum,
       G$s[match(G$s$SampleID,attr(G$uur,"Labels")),],
       permutations = nperm)

adonis(G$uur ~ diet + log.weight,
       G$s[match(G$s$SampleID,attr(G$uur,"Labels")),],
       permutations = nperm)

adonis(G$uur ~ diet + Lab.Wild.Domestic,
       G$s[match(G$s$SampleID,attr(G$uur,"Labels")),],
       permutations = nperm)

adonis(G$uur ~ common,
       G$s[match(G$s$SampleID,attr(G$uur,"Labels")),],
       permutations = nperm)


G$uurMammal <- as.dist(as.matrix(G$uur)[match(G$s$SampleID[G$s$class=="Mammalia"],attr(G$uur,"Labels")),
                                        match(G$s$SampleID[G$s$class=="Mammalia"],attr(G$uur,"Labels"))])

adonis(G$uurMammal ~ order == "Carnivora",G$s[match(attr(G$uurMammal,"Labels"),G$s$SampleID),],
       permutations = nperm)

adonis(G$uur ~ class + I(order == "Carnivora"),G$s[match(attr(G$uur,"Labels"),G$s$SampleID),],
       permutations = nperm)

adonis(G$uur ~ I(class == "Mammalia") + I(order == "Carnivora"),G$s[match(attr(G$uur,"Labels"),G$s$SampleID),],
       permutations = nperm)

adonis(G$uur ~ order == "Carnivora",G$s[match(attr(G$uur,"Labels"),G$s$SampleID),],
       permutations = nperm)

adonis(G$uur ~ order == "Cetacea",G$s[match(attr(G$uur,"Labels"),G$s$SampleID),],
       permutations = nperm)

adonis(G$uur ~ class == "Chondrichthyes",G$s[match(attr(G$uur,"Labels"),G$s$SampleID),],
       permutations = nperm)

adonis(G$uur ~ order == "Artiodactyla" | order == "Lagomorpha",G$s[match(attr(G$uur,"Labels"),G$s$SampleID),],
       permutations = nperm)

adonis(G$uur ~ order == "Primates" | order == "Rodentia",G$s[match(attr(G$uur,"Labels"),G$s$SampleID),],
       permutations = nperm)

adonis(G$uur ~ class =="Insecta", G$s[match(attr(G$uur,"Labels"),G$s$SampleID),],
       permutations = nperm)



G$uurMatrix <- as.matrix(G$uur)
G$bcpMatrix <- as.matrix(G$bcp)






adonis(G$uur ~ common, G$s[match(attr(G$uur,"Labels"),G$s$SampleID),],
       permutations = nperm)

adonis(G$bcr ~ common, G$s[match(attr(G$uur,"Labels"),G$s$SampleID),],
       permutations = nperm)

adonis(G$wu ~ common, G$s[match(attr(G$uur,"Labels"),G$s$SampleID),],
       permutations = nperm)

adonis(G$bcp ~ common, G$s[match(attr(G$uur,"Labels"),G$s$SampleID),],
       permutations = nperm)






### Adonis on species phylogeny

resetData("toUse")
resetUnifracData()

G$phyDist <- as.dist(G$phyMatrix)

adonis(G$phyDist ~ diet, G$species, permutations = nperm)
adonis(G$phyDist ~ specificDiet, G$species, permutations = nperm)
adonis(G$phyDist ~ phylum, G$species, permutations = nperm)
adonis(G$phyDist ~ class, G$species, permutations = nperm)




###
### Adonis on species centroids
###

adonis(G$uurDistM ~ diet, G$species, permutations = nperm)
adonis(G$uurDistM ~ specificDiet, G$species, permutations = nperm)
adonis(G$uurDistM ~ diet + specificDiet, G$species, permutations = nperm)
adonis(G$uurDistM ~ phylum, G$species, permutations = nperm)
adonis(G$uurDistM ~ phylum + diet, G$species, permutations = nperm)
adonis(G$uurDistM ~ class + diet, G$species, permutations = nperm)
adonis(G$uurDistM ~ order + diet, G$species, permutations = nperm)
adonis(G$uurDistM ~ diet + order, G$species, permutations = nperm)










# Fisher's tests----
# Fisher's test for unassigned OTUs vs usearch chimeras
fisher.test(as.matrix(rbind(c(145842,8814225),c(630689,22588000))))
# Fisher's test for singletons vs usearch chimeras
fisher.test(as.matrix(rbind(c(326406,16058371),c(450125,15343854))))
# Fisher's test for unassigned OTUs vs chimera_slayer chimeras
fisher.test(as.matrix(rbind(c(241,150944),c(1024,627523))))
# Fisher's test for singletons vs chimera_slayer chimeras
fisher.test(as.matrix(rbind(c(525,334908),c(740,443559))))
