#packages
library(ape)
library(phytools)
library(phangorn)


#Import Data - trees
#################################################################################
#Make original consensus tree

#all_tree <- read.tree("species_trees_renamed.tre")
#all_tree <- all_tree[sapply(all_tree, inherits, "phylo")]
#tr1 <- consensus(all_tree, p = 0.5, check.labels = TRUE, rooted = FALSE)
#tr1<-consensus.edges(all_tree,if.absent="ignore")
#trimmed_tree <- drop.tip(tr1, c("0168", "New_sequence"))
#root_node <- getMRCA(trimmed_tree, c("L._rupicola_1", "L._rupicola_2"))
#tr1reroot <- reroot(trimmed_tree, node.number=root_node, position=NULL, interactive=FALSE)
#trimmed_tree$tip.label
#is.rooted(tr1reroot)
#plot(tr1reroot, show.node.label = TRUE, use.edge.length = T)
#write.tree(tr1reroot, file="original_consensus.tre")

#open original consensus tree
tr1reroot <- read.tree("original_consensus.tre")
##################################################################################
#Large consensus tree
#all_tree_L <- read.tree("large_trees_consensus_renamed.tre")
#all_tree_L <- all_tree_L[sapply(all_tree_L, inherits, "phylo")]
#tr2 <- consensus(all_tree_L, p = 0.5, check.labels = TRUE, rooted = FALSE)
#tr2<-consensus.edges(all_tree_L,if.absent="ignore")
#plot(tr2)

#root_node <- getMRCA(tr2, c("D._album", "R._elegans", "R._augusta", "I._mirabilis", "I._polymorpha", "I._thalangensis", "I._namsabiensis" ,"L._rupicola"))
#tr2reroot <- reroot(tr2, node.number=root_node, position=NULL, interactive=FALSE)
plot(tr2reroot, show.node.label = TRUE, use.edge.length = T, cex=0.4)
#title("Large consensus")
#is.rooted(tr2reroot)
#write.tree(tr2reroot, file="large_consensus.tre")

tr2reroot <- read.tree("large_consensus.tre")

################################################################################
#Large ML wAstral tree all genes (freerate)
tr3 <- read.tree("large_ML_renamed.tre")

#rerooting the tree
root_node <- getMRCA(tr3, c("D._album", "R._elegans", "R._augusta", "I._mirabilis", "I._polymorpha", "I._thalangensis", "I._namsabiensis" ,"L._rupicola"))
tr3reroot <- reroot(tr3, node.number=root_node, position=NULL, interactive=FALSE)
################################################################################
#Large Bayes wAstral tree
tr4 <- read.tree("large_Bayes_renamed.tre")
#rerooting the tree
root_node <- getMRCA(tr4, c("D._album", "R._elegans", "R._augusta", "I._mirabilis", "I._polymorpha", "I._thalangensis", "I._namsabiensis" ,"L._rupicola"))
tr4reroot <- reroot(tr4, node.number=root_node, interactive=FALSE)

################################################################################
#Large ML wAstral tree converged genes + mrb modesl only
tr5 <- read.tree("large_ML_converged_renamed.tre")

#rerooting the tree
root_node <- getMRCA(tr5, c("D._album", "R._elegans", "R._augusta", "I._mirabilis", "I._polymorpha", "I._thalangensis", "I._namsabiensis" ,"L._rupicola"))
tr5reroot <- reroot(tr5, node.number=root_node, position=NULL, interactive=FALSE)

################################################################################
#Original ML wASTRAL (all genes only mrb models) 
tr6 <- read.tree("Heyduk_astral_ML_renamed.tre")

#rerooting the tree
root_node <- getMRCA(tr6, c("L._rupicola_1", "L._rupicola_2"))
tr6reroot <- reroot(tr6, node.number=root_node, position=NULL, interactive=FALSE)

################################################################################
# Original bayesian wASTRAL tree
tr7<- read.tree("Original_astral_MrB_renamed.tre")
trimmed_tree <- drop.tip(tr7, c("0168", "New_sequence"))
root_node <- getMRCA(trimmed_tree, c("L._rupicola_1", "L._rupicola_2"))
tr7reroot <- reroot(trimmed_tree, node.number=root_node, interactive=FALSE)
is.rooted(tr7reroot)

################################################################################
#Heyduk IQ-tree wASTRAL (only genes that converged in MrBayes+ only mrb models)
tr8 <- read.tree("original_converged_ML_renamed.tre")

#rerooting the tree
root_node <- getMRCA(tr8, c("L._rupicola_1", "L._rupicola_2"))
tr8reroot <- reroot(tr8, node.number=root_node, position=NULL, interactive=FALSE)
tr7$tip.label
################################################################################
#Heyduk IQ-tree wASTRAL freerate
tr9 <- read.tree("original_astral_MLfreerate_renamed.tre")


#rerooting the tree
root_node <- getMRCA(tr9, c("L._rupicola_1", "L._rupicola_2"))
tr9reroot <- reroot(tr9, node.number=root_node, position=NULL, interactive=FALSE)

################################################################################
#Large ML wAstral tree all genes + only models in mrb
tr10 <- read.tree("large_astral_ML_mrb_model_renamed.tre")

#rerooting the tree
root_node <- getMRCA(tr10, c("D._album", "R._elegans", "R._augusta", "I._mirabilis", "I._polymorpha", "I._thalangensis", "I._namsabiensis" ,"L._rupicola"))
tr10reroot <- reroot(tr10, node.number=root_node, position=NULL, interactive=FALSE)




################################################################################
#Importing scored trees
################################################################################
tr3_info <- read.tree("large_scored_ML_freerate_scored.tre")
root_node <- getMRCA(tr3_info, c("D._album", "R._elegans", "R._augusta", "I._mirabilis", "I._polymorpha", "I._thalangensis", "I._namsabiensis" ,"L._rupicola"))
tr3_info<- reroot(tr3_info, node.number=root_node, interactive=FALSE)
EN_values <- sub(".*EN=(.*?)\\].*", "\\1",tr3_info$node.label)
length(EN_values) == Nnode(tr3reroot)
EN_values<- format(round(as.numeric(EN_values), 2), nsmall = 2)

EN_values_max=0
EN_values_min=0
EN_values_mean=0
EN_values_max <- max(EN_values)
EN_values_min <- min(as.numeric(EN_values), na.rm = T)
EN_values_mean <- mean(as.numeric(EN_values), na.rm = T)

EN_values_tr3 <- cbind(EN_values_min, EN_values_max, EN_values_mean)
################################################################################

tr4_info<- read.tree("large_scored_wastral_Bayes_renamed.tre")
root_node <- getMRCA(tr4_info, c("D._album", "R._elegans", "R._augusta", "I._mirabilis", "I._polymorpha", "I._thalangensis", "I._namsabiensis" ,"L._rupicola"))
tr4_info<- reroot(tr4_info, node.number=root_node, interactive=FALSE)
EN_values <- sub(".*EN=(.*?)\\].*", "\\1",tr4_info$node.label)
length(EN_values) == Nnode(tr4reroot)
EN_values<- format(round(as.numeric(EN_values), 2), nsmall = 2)

EN_values_max=0
EN_values_min=0
EN_values_mean=0
EN_values_max <- max(EN_values)
EN_values_min <- min(as.numeric(EN_values), na.rm = T)
EN_values_mean <- mean(as.numeric(EN_values), na.rm = T)

EN_values_tr4 <- cbind(EN_values_min, EN_values_max, EN_values_mean)

#extra for first test
#PP1_values <- gsub(".*pp1=(.*?)\\;.*", "\\1",tr4_info$node.label)
#PP1_values <- gsub(".*pp1=(.*?)\\-.*", "\\1",PP1_values)
#M <- cbind(PP1_values, tr4reroot$node.label)
# Match nodes from tr4_info (rerooted tree with labels) to tr4 (original tree)
matched_nodes <- matchNodes(tr4_info, tr4reroot)



################################################################################
tr5_info <- read.tree("large_scored_ML_converged.tre")
root_node <- getMRCA(tr5_info, c("D._album", "R._elegans", "R._augusta", "I._mirabilis", "I._polymorpha", "I._thalangensis", "I._namsabiensis" ,"L._rupicola"))
tr5_info<- reroot(tr5_info, node.number=root_node, interactive=FALSE)
EN_values <- sub(".*EN=(.*?)\\].*", "\\1",tr5_info$node.label)
length(EN_values) == Nnode(tr5reroot)
EN_values<- format(round(as.numeric(EN_values), 2), nsmall = 2)

EN_values_max=0
EN_values_min=0
EN_values_mean=0
EN_values_max <- max(EN_values)
EN_values_min <- min(as.numeric(EN_values), na.rm = T)
EN_values_mean <- mean(as.numeric(EN_values), na.rm = T)

EN_values_tr5 <- cbind(EN_values_min, EN_values_max, EN_values_mean)

################################################################################
tr6_info<- read.tree("org_scored_wastral_ML_renamed.tre")
tr6_info<- drop.tip(tr6_info, c("0168", "New_sequence"))
root_node <- getMRCA(tr6_info, c("L._rupicola_1", "L._rupicola_2"))
tr6_info<- reroot(tr6_info, node.number=root_node, interactive=FALSE)
EN_values <- as.numeric(gsub(".*EN=(.*?)\\].*", "\\1",tr6_info$node.label))
length(EN_values) == Nnode(tr6reroot)

EN_values<- format(round(as.numeric(EN_values), 2), nsmall = 2)

EN_values_max=0
EN_values_min=0
EN_values_mean=0
EN_values_max <- max(EN_values)
EN_values_min <- min(as.numeric(EN_values), na.rm = T)
EN_values_mean <- mean(as.numeric(EN_values), na.rm = T)

EN_values_tr6 <- cbind(EN_values_min, EN_values_max, EN_values_mean)

################################################################################
tr7_info<- read.tree("Original_scored_astral_MrB_renamed.tre")
tr7_info<- drop.tip(tr7_info, c("0168", "New_sequence"))
root_node <- getMRCA(tr7_info, c("L._rupicola_1", "L._rupicola_2"))
tr7_info<- reroot(tr7_info, node.number=root_node, interactive=FALSE)
EN_values <- as.numeric(gsub(".*EN=(.*?)\\].*", "\\1",tr7_info$node.label))
length(EN_values) == Nnode(tr7reroot)
EN_values<- format(round(as.numeric(EN_values), 2), nsmall = 2)

EN_values_max=0
EN_values_min=0
EN_values_mean=0
EN_values_max <- max(EN_values)
EN_values_min <- min(as.numeric(EN_values), na.rm = T)
EN_values_mean <- mean(as.numeric(EN_values), na.rm = T)

EN_values_tr7 <- cbind(EN_values_min, EN_values_max, EN_values_mean)

#PP1_values <- as.numeric(gsub(".*pp1=(.*?)\\;.*", "\\1",info_nodelabels))
#M <- cbind(PP1_values, tr7reroot$node.label)
################################################################################
tr8_info<- read.tree("org_scored_ML_converged_renamed.tre")
tr8_info<- drop.tip(tr8_info, c("0168", "New_sequence"))
root_node <- getMRCA(tr8_info, c("L._rupicola_1", "L._rupicola_2"))
tr8_info<- reroot(tr8_info, node.number=root_node, interactive=FALSE)
EN_values <- as.numeric(gsub(".*EN=(.*?)\\].*", "\\1",tr8_info$node.label))
length(EN_values) == Nnode(tr8reroot)

EN_values<- format(round(as.numeric(EN_values), 2), nsmall = 2)

EN_values_max=0
EN_values_min=0
EN_values_mean=0
EN_values_max <- max(EN_values)
EN_values_min <- min(as.numeric(EN_values), na.rm = T)
EN_values_mean <- mean(as.numeric(EN_values), na.rm = T)

EN_values_tr8 <- cbind(EN_values_min, EN_values_max, EN_values_mean)

################################################################################
tr9_info <- read.tree("org_scored__ML_freerate_renamed.tre")
tr9_info<- drop.tip(tr9_info, c("0168", "New_sequence"))
root_node <- getMRCA(tr9_info, c("L._rupicola_1", "L._rupicola_2"))
tr9_info<- reroot(tr9_info, node.number=root_node, interactive=FALSE)
EN_values <- as.numeric(gsub(".*EN=(.*?)\\].*", "\\1",tr9_info$node.label))
length(EN_values) == Nnode(tr9reroot)

EN_values<- format(round(as.numeric(EN_values), 2), nsmall = 2)

EN_values_max=0
EN_values_min=0
EN_values_mean=0
EN_values_max <- max(EN_values)
EN_values_min <- min(as.numeric(EN_values), na.rm = T)
EN_values_mean <- mean(as.numeric(EN_values), na.rm = T)

EN_values_tr9 <- cbind(EN_values_min, EN_values_max, EN_values_mean)


################################################################################
tr10_info<- read.tree("large_scored_ML_renamed.tre")
root_node <- getMRCA(tr10_info, c("D._album", "R._elegans", "R._augusta", "I._mirabilis", "I._polymorpha", "I._thalangensis", "I._namsabiensis" ,"L._rupicola"))
tr10_info<- reroot(tr10_info, node.number=root_node, interactive=FALSE)
EN_values <- sub(".*EN=(.*?)\\].*", "\\1",tr10_info$node.label)
length(EN_values) == Nnode(tr10reroot)
EN_values<- format(round(as.numeric(EN_values), 2), nsmall = 2)

EN_values_max=0
EN_values_min=0
EN_values_mean=0
EN_values_max <- max(EN_values)
EN_values_min <- min(as.numeric(EN_values), na.rm = T)
EN_values_mean <- mean(as.numeric(EN_values), na.rm = T)

EN_values_tr10 <- cbind(EN_values_min, EN_values_max, EN_values_mean)

################################################################################

EN_values_all <- rbind(EN_values_tr3, EN_values_tr4,EN_values_tr5,EN_values_tr6,EN_values_tr7,EN_values_tr8,EN_values_tr9 ,EN_values_tr10)




################################################################################
#Plotting the trees in PDF
#################################################################################
#from liam revell phytools blog
split.plotTree<-function(tree,splits=NULL,file=NULL,fn,...){
  ef<-0.037037037037
  if(!is.null(file)) pdf(file, height=15, width=10)
  if(is.null(splits)) splits<-(floor(0.5*Ntip(tree))+0.5)/Ntip(tree)
  S<-matrix(c(0,splits,splits,1+1/Ntip(tree)),length(splits)+1,2)
  S<-cbind(S[,1]+ef*(S[,2]-S[,1]),S[,2]-ef*(S[,2]-S[,1]))
  for(i in nrow(S):1){
    if(is.null(file)&&i<nrow(S)) par(ask=TRUE)
    plotTree(tree,ylim=Ntip(tree)*S[i,],use.edge.length=FALSE,...)
    fn()
  }
  if(!is.null(file)) oo<-dev.off()
}

################################################################################
#plotting tr1
#removing decimals from LPP values for plotting
for (i in seq_along(tr1reroot$node.label)) {
  if(tr1reroot$node.label[i] == "Root" ) {next
  } else {
    tr1reroot$node.label[i] <- format(round(as.numeric(tr1reroot$node.label[i]), 2), nsmall = 2)
  }
}
tr1reroot$node.label<- ifelse(tr1reroot$node.label == "1.00", "1", tr1reroot$node.label)
tr1reroot$tip.label <-gsub("_"," ",tr1reroot$tip.label)
pdf("original_consensus_done.pdf", width = 11, height = 27)
plot.phylo(tr1reroot,use.edge.length = TRUE, show.tip.label = FALSE, x.lim = c(0, max(branching.times(tr1reroot)) + 0.06))
obj<-get("last_plot.phylo",envir=.PlotPhyloEnv)
text(rep(max(obj$xx[1:Ntip(tr1reroot)]),Ntip(tr1reroot)),obj$yy[1:Ntip(tr1reroot)],
     labels=tr1reroot$tip.label,pos=4)
for(i in 1:Ntip(tr1reroot)) lines(c(obj$xx[i],max(obj$xx[1:Ntip(tr1reroot)])),
                                  rep(obj$yy[i],2),lty="dotted")
nodelabels(tr1reroot$node.label,adj = c(0,0.4), bg="transparent", frame="none", cex=0.7)
title("Original 50 majority rule consensus tree")
dev.off()



#################################################################################
#plotting tr 2
#moving decimals from PP values for plotting
for (i in seq_along(tr2reroot$node.label)) {
  if(tr2reroot$node.label[i] == "Root" ) {next
  } else {
    tr2reroot$node.label[i] <- format(round(as.numeric(tr2reroot$node.label[i]), 2), nsmall = 2)
  }
}

tr2reroot$node.label<- ifelse(tr2reroot$node.label == "1.00", "1", tr2reroot$node.label)

pdf("consensus_large.pdf", height=35, width=20)
plot.phylo(tr2reroot,use.edge.length = TRUE, label.offset=0.0005)
obj<-get("last_plot.phylo",envir=.PlotPhyloEnv)
text(rep(max(obj$xx[1:Ntip(tr2reroot)]),Ntip(tr2reroot)),obj$yy[1:Ntip(tr2reroot)],
     labels=tr2reroot$tip.label,pos=4)
for(i in 1:Ntip(tr1reroot)) lines(c(obj$xx[i],max(obj$xx[1:Ntip(tr1reroot)])),
                                  rep(obj$yy[i],2),lty="dotted")
nodelabels(tr2reroot$node.label,adj = c(1,- 0.2), bg="transparent", frame="none")
title("Large majority rule consensus")
dev.off()




#################################################################################
#plotting tr 3

#Adding terminal branchlegths
tr3reroot$edge.length[is.na(tr3reroot$edge.length)] <- 0.2
#moving decimals from LPP values for plotting
for (i in seq_along(tr3reroot$node.label)) {
  if(tr3reroot$node.label[i] == "Root" ) {next
  } else {
    tr3reroot$node.label[i] <- format(round(as.numeric(tr3reroot$node.label[i]), 2), nsmall = 2)
  }
}
tr3reroot$node.label<- ifelse(tr3reroot$node.label == "1.00", "1", tr3reroot$node.label)
#plotting the tree
EN_tr3 <- sub(".*EN=(.*?)\\].*", "\\1",tr3_info$node.label)
EN_tr3 <- as.numeric(EN_tr3)
EN_tr3[is.na(EN_tr3)] <- 0
length(EN_tr3) == tr3reroot$Nnode
EN_tr3 <- round(EN_tr3, 0)
matched_nodes_3 <- matchNodes(tr3_info, tr3reroot)
x <- cbind(matched_nodes_3, EN_tr3)

add_support_labels<-function(node=NULL,support,
                             cols=c("white","black"),cex=1.5){
  scale<-range(support)
  support<-(support-scale[1])/diff(scale)
  pp<-get("last_plot.phylo",envir=.PlotPhyloEnv)
  if(is.null(node)) node<-1:pp$Nnode+pp$Ntip
  root_node <- pp$Ntip + 1 
  colfunc<-colorRamp(cols)
  node_cols<-colfunc(support)/255
  x<-pp$xx[node] 
  y<-pp$yy[node]
  for(i in 1:nrow(node_cols)){
    if(!is.na(support[i]) && node[i] != root_node)
      points(x[i],y[i],pch=21,,cex=0.9,
             bg=rgb(node_cols[i,1],node_cols[i,2],node_cols[i,3]))
  }
  invisible(scale)
}



tree <- tr3reroot

tree$node.label<- ifelse(tree$node.label == "1.00", "1", tree$node.label)

foo<-function(){
  par(fg="black")
  nodelabels(text=nodelabels_for_plot,  bg= "white" , frame="none", width = 0.2, adj=c(1.05, -0.1)) 
  minmax<-add_support_labels(support=EN_tr3)
  add.color.bar(leg=0.3*max(nodeHeights(tree)),
                cols=colorRampPalette(c("white","black"))(100),
                lims=minmax,title="EN",
                subtitle="",prompt=FALSE,x=0.01*par()$usr[2],
                y=0.02*par()$usr[4])
  
}
splits<- 0.543
tree <- compute.brlen(tree, method = "Grafen", power = 0.3) 
terminal_edges <- which(tree$edge[,2] <= length(tree$tip.label))
tree$edge.length[terminal_edges] <- tree$edge.length[terminal_edges]-0.18
nodelabels_for_plot <- paste(tree$node.label,"/", EN_tr3)

split.plotTree(tree,splits,file= "large_ML_freerate_9.pdf", ftype="i",fn=foo,lwd=1, cex=0.7, oma = c(1, 1, 1, 1), mar = c(1, 1, 1, 1))


tr3_info$node.label


#pdf("large_freerate_ML.pdf", height=35, width=20)
#plotBranchbyTrait(tr3reroot, EN_tr3, mode = "node", type = "phylogram", edge.width = 2, palette = gray_pal)
#nodelabels(tr3reroot$node.label,adj = c(1,- 0.2), bg="transparent", frame="none")
#title("Large ML all genes (free rate models)")
#dev.off()

#plotBranchbyTrait(tr3reroot, EN_tr3, mode = "node", type = "phylogram", edge.width = 2, palette = "gray")


#################################################################################
#plotting tr 4
#Adding terminal branchlegths
tr4reroot$edge.length[is.na(tr4reroot$edge.length)] <- 0.2
#moving decimals from LPP values
for (i in seq_along(tr4reroot$node.label)) {
  if(tr4reroot$node.label[i] == "Root" ) {next
  } else {
    tr4reroot$node.label[i] <- format(round(as.numeric(tr4reroot$node.label[i]), 2), nsmall = 2)
  }
}

tr4reroot$node.label<- ifelse(tr4reroot$node.label == "1.00", "1", tr4reroot$node.label)
EN_tr4 <-  as.numeric(sub(".*EN=(.*?)\\].*", "\\1",tr4_info$node.label))
EN_tr4[is.na(EN_tr4)] <- 0
length(EN_tr4) == tr4reroot$Nnode
matched_nodes_4 <- matchNodes(tr4_info, tr4reroot)
EN_tr4 <- round(EN_tr4, 0)


#from liam revell phytools blog
add_support_labels<-function(node=NULL,support,
                             cols=c("white","black"),cex=1.5){
  scale<-range(support)
  support<-(support-scale[1])/diff(scale)
  pp<-get("last_plot.phylo",envir=.PlotPhyloEnv)
  if(is.null(node)) node<-1:pp$Nnode+pp$Ntip
  root_node <- pp$Ntip + 1 
  colfunc<-colorRamp(cols)
  node_cols<-colfunc(support)/255
  x<-pp$xx[node] 
  y<-pp$yy[node]
  for(i in 1:nrow(node_cols)){
    if(!is.na(support[i]) && node[i] != root_node)
      points(x[i],y[i],pch=21,,cex=0.9,
             bg=rgb(node_cols[i,1],node_cols[i,2],node_cols[i,3]))
  }
  invisible(scale)
}


tree <- tr4reroot

tree$node.label<- ifelse(tree$node.label == "1.00", "1", tree$node.label)

foo<-function(){
  par(fg="black")
  nodelabels(text=nodelabels_for_plot,  bg= "white" , frame="none", width = 0.2, adj=c(1.05, -0.1)) 
  minmax<-add_support_labels(support=EN_tr4)
  add.color.bar(leg=0.3*max(nodeHeights(tree)),
                cols=colorRampPalette(c("white","black"))(100),
                lims=minmax,title="EN",
                subtitle="",prompt=FALSE,x=0.01*par()$usr[2],
                y=0.02*par()$usr[4])
  
}
splits<- 0.533
tree <- compute.brlen(tree, method = "Grafen", power = 0.3) 
terminal_edges <- which(tree$edge[,2] <= length(tree$tip.label))
tree$edge.length[terminal_edges] <- tree$edge.length[terminal_edges]-0.18
nodelabels_for_plot <- paste(tree$node.label)

split.plotTree(tree,splits,file= "large_ML_Bayesian_genes_5.pdf", ftype="i",fn=foo,lwd=1, cex=0.7, oma = c(1, 1, 1, 1), mar = c(1, 1, 1, 1))






pdf("large_Bayesian.pdf", height=35, width=20)
plot.phylo(tr4reroot,use.edge.length = TRUE, label.offset=0.0005)
nodelabels(tr4reroot$node.label,adj = c(1,- 0.2), bg="transparent", frame="none")
title("Large dataset Bayesian tree")
dev.off()


#################################################################################
#plotting tr 5
#Adding terminal branch legths
tr5reroot$edge.length[is.na(tr5reroot$edge.length)] <- 0.2

#moving decimals from LPP values
for (i in seq_along(tr5reroot$node.label)) {
  if(tr5reroot$node.label[i] == "Root" ) {next
  } else {
    tr5reroot$node.label[i] <- format(round(as.numeric(tr5reroot$node.label[i]),2), nsmall = 2)
  }
}
tr5reroot$node.label<- ifelse(tr5reroot$node.label == "1.00", "1", tr5reroot$node.label)


EN_tr5 <-  as.numeric(sub(".*EN=(.*?)\\].*", "\\1",tr5_info$node.label))
EN_tr5[is.na(EN_tr5)] <- 0
length(EN_tr5) == tr5reroot$Nnode
matched_nodes_5 <- matchNodes(tr5_info, tr5reroot)
EN_tr5 <- round(EN_tr5, 0)


#from liam revell phytools blog
add_support_labels<-function(node=NULL,support,
                             cols=c("white","black"),cex=1.5){
  scale<-range(support)
  support<-(support-scale[1])/diff(scale)
  pp<-get("last_plot.phylo",envir=.PlotPhyloEnv)
  if(is.null(node)) node<-1:pp$Nnode+pp$Ntip
  root_node <- pp$Ntip + 1 
  colfunc<-colorRamp(cols)
  node_cols<-colfunc(support)/255
  x<-pp$xx[node] 
  y<-pp$yy[node]
  for(i in 1:nrow(node_cols)){
    if(!is.na(support[i]) && node[i] != root_node)
      points(x[i],y[i],pch=21,,cex=0.9,
             bg=rgb(node_cols[i,1],node_cols[i,2],node_cols[i,3]))
  }
  invisible(scale)
}


tree <- tr5reroot

tree$node.label<- ifelse(tree$node.label == "1.00", "1", tree$node.label)

foo<-function(){
  par(fg="black")
  nodelabels(text=nodelabels_for_plot,  bg= "white" , frame="none", width = 0.2, adj=c(1.05, -0.1)) 
  minmax<-add_support_labels(support=EN_tr5)
  add.color.bar(leg=0.3*max(nodeHeights(tree)),
                cols=colorRampPalette(c("white","black"))(100),
                lims=minmax,title="EN",
                subtitle="",prompt=FALSE,x=0.01*par()$usr[2],
                y=0.02*par()$usr[4])
  
}
splits<- 0.533
tree <- compute.brlen(tree, method = "Grafen", power = 0.3) 
terminal_edges <- which(tree$edge[,2] <= length(tree$tip.label))
tree$edge.length[terminal_edges] <- tree$edge.length[terminal_edges]-0.18
nodelabels_for_plot <- paste(tree$node.label)

split.plotTree(tree,splits,file= "large_ML_converged_Done.pdf", ftype="i",fn=foo,lwd=1, cex=0.7, oma = c(1, 1, 1, 1), mar = c(1, 1, 1, 1))





pdf("large_ML_mrb_models_converged.pdf", height=35, width=20)
plot.phylo(tr5reroot,use.edge.length = TRUE, label.offset=0.0005)
nodelabels(tr5reroot$node.label, bg="transparent", frame="none")
title("Large dataset ML converged genes")
dev.off()



#################################################################################
#plotting tr 6
#Adding terminal branch legths
tr6reroot$edge.length[is.na(tr6reroot$edge.length)] <- 0.2


#moving decimals from LPP values
for (i in seq_along(tr6reroot$node.label)) {
  if(tr6reroot$node.label[i] == "Root" ) {next
  } else {
    tr6reroot$node.label[i] <- format(round(as.numeric(tr6reroot$node.label[i]), 2), nsmall = 2)
    if (i %% 1 == 0) as.integer(i) else i
  }
}
tr6reroot$node.label<- ifelse(tr6reroot$node.label == "1.00", "1", tr6reroot$node.label)

EN_tr6 <-  as.numeric(sub(".*EN=(.*?)\\].*", "\\1",tr6_info$node.label))
EN_tr6[is.na(EN_tr6)] <- 0
length(EN_tr6) == tr6reroot$Nnode
matched_nodes_6 <- matchNodes(tr6_info, tr6reroot)
EN_tr6 <- round(EN_tr6, 0)


#from liam revell phytools blog
add_support_labels<-function(node=NULL,support,
                             cols=c("white","black"),cex=1.5){
  scale<-range(support)
  support<-(support-scale[1])/diff(scale)
  pp<-get("last_plot.phylo",envir=.PlotPhyloEnv)
  if(is.null(node)) node<-1:pp$Nnode+pp$Ntip
  root_node <- pp$Ntip + 1 
  colfunc<-colorRamp(cols)
  node_cols<-colfunc(support)/255
  x<-pp$xx[node] 
  y<-pp$yy[node]
  for(i in 1:nrow(node_cols)){
    if(!is.na(support[i]) && node[i] != root_node)
      points(x[i],y[i],pch=21,,cex=0.9,
             bg=rgb(node_cols[i,1],node_cols[i,2],node_cols[i,3]))
  }
  invisible(scale)
}

tree <- tr6reroot

tree$node.label<- ifelse(tree$node.label == "1.00", "1", tree$node.label)

foo<-function(){
  par(fg="black")
  nodelabels(text=nodelabels_for_plot,  bg= "white" , frame="none", width = 0.2, adj=c(1.05, -0.1)) 
  minmax<-add_support_labels(support=EN_tr6)
  add.color.bar(leg=0.3*max(nodeHeights(tree)),
                cols=colorRampPalette(c("white","black"))(100),
                lims=minmax,title="EN",
                subtitle="",prompt=FALSE,x=0.01*par()$usr[2],
                y=0.02*par()$usr[4])
  
}
splits<- 0.586
tree <- compute.brlen(tree, method = "Grafen", power = 0.3) 
terminal_edges <- which(tree$edge[,2] <= length(tree$tip.label))
tree$edge.length[terminal_edges] <- tree$edge.length[terminal_edges]-0.18
nodelabels_for_plot <- paste(tree$node.label)

split.plotTree(tree,splits,file= "Original_ML_all_genes_done.pdf", ftype="i",fn=foo,lwd=1, cex=0.7, oma = c(1, 1, 1, 1), mar = c(1, 1, 1, 1))





pdf("original_ML_all.pdf", height=35, width=18)
plot.phylo(tr6reroot,use.edge.length = TRUE, cex=1, label.offset=0.0005)
nodelabels(tr6reroot$node.label,adj = c(1,- 0.2), bg="transparent", frame="none")
title("Original dataset ML all genes")
dev.off()


#################################################################################
#plotting tr 7

#Adding terminal branch legths
tr7reroot$edge.length[is.na(tr7reroot$edge.length)] <- 0.1


#moving decimals from LPP values
for (i in seq_along(tr7reroot$node.label)) {
  if(tr7reroot$node.label[i] == "Root" ) {next
  } else {
    tr7reroot$node.label[i] <- format(round(as.numeric(tr7reroot$node.label[i]), 2), nsmall = 2)
  }
}
tr7reroot$node.label<- ifelse(tr7reroot$node.label == "1.00", "1", tr7reroot$node.label)


EN_tr7 <-  as.numeric(sub(".*EN=(.*?)\\].*", "\\1",tr7_info$node.label))
EN_tr7[is.na(EN_tr7)] <- 0
length(EN_tr7) == tr7reroot$Nnode
matched_nodes_7 <- matchNodes(tr7_info, tr7reroot)
EN_tr7 <- round(EN_tr7, 0)


#from liam revell phytools blog
add_support_labels<-function(node=NULL,support,
                             cols=c("white","black"),cex=1.5){
  scale<-range(support)
  support<-(support-scale[1])/diff(scale)
  pp<-get("last_plot.phylo",envir=.PlotPhyloEnv)
  if(is.null(node)) node<-1:pp$Nnode+pp$Ntip
  root_node <- pp$Ntip + 1 
  colfunc<-colorRamp(cols)
  node_cols<-colfunc(support)/255
  x<-pp$xx[node] 
  y<-pp$yy[node]
  for(i in 1:nrow(node_cols)){
    if(!is.na(support[i]) && node[i] != root_node)
      points(x[i],y[i],pch=21,,cex=0.9,
             bg=rgb(node_cols[i,1],node_cols[i,2],node_cols[i,3]))
  }
  invisible(scale)
}

tree <- tr7reroot

tree$node.label<- ifelse(tree$node.label == "1.00", "1", tree$node.label)

foo<-function(){
  par(fg="black")
  nodelabels(text=nodelabels_for_plot,  bg= "white" , frame="none", width = 0.2, adj=c(1.05, -0.1)) 
  minmax<-add_support_labels(support=EN_tr7)
  add.color.bar(leg=0.3*max(nodeHeights(tree)),
                cols=colorRampPalette(c("white","black"))(100),
                lims=minmax,title="EN",
                subtitle="",prompt=FALSE,x=0.01*par()$usr[2],
                y=0.02*par()$usr[4])
  
}
splits<- 0.564
tree <- compute.brlen(tree, method = "Grafen", power = 0.3) 
terminal_edges <- which(tree$edge[,2] <= length(tree$tip.label))
tree$edge.length[terminal_edges] <- tree$edge.length[terminal_edges]-0.18
nodelabels_for_plot <- paste(tree$node.label)

split.plotTree(tree,splits,file= "Original_Bayesian_genes_done.pdf", ftype="i",fn=foo,lwd=1, cex=0.7, oma = c(1, 1, 1, 1), mar = c(1, 1, 1, 1))




plot(tr7reroot, show.node.label = TRUE, use.edge.length =F, cex=0.5)
title("Original Bayes tree")


pdf("original_Bayesian_1.pdf", height=35, width=20)
plot.phylo(tr7reroot,use.edge.length = FALSE, label.offset=0.0005)
plotBranchbyTrait(tr7reroot, EN_values, mode="nodes", palette="grey")
nodelabels(EN_values,adj = c(1,- 0.2), bg="transparent", frame="none")
title("original dataset Bayesian")
dev.off()




#################################################################################
#plotting tr 8

tr8reroot$edge.length[is.na(tr8reroot$edge.length)] <- 0.2

#removing decimals from PP values
for (i in seq_along(tr8reroot$node.label)) {
  if(tr8reroot$node.label[i] == "Root" ) {next
  } else {
    tr8reroot$node.label[i] <- format(round(as.numeric(tr8reroot$node.label[i]), 2), nsmall = 2)
  }
}
tr8reroot$node.label<- ifelse(tr8reroot$node.label == "1.00", "1", tr8reroot$node.label)

EN_tr8 <-  as.numeric(sub(".*EN=(.*?)\\].*", "\\1",tr8_info$node.label))
EN_tr8[is.na(EN_tr8)] <- 0
length(EN_tr8) == tr8reroot$Nnode
matched_nodes_8 <- matchNodes(tr8_info, tr8reroot)
EN_tr8 <- round(EN_tr8, 0)


#from liam revell phytools blog
add_support_labels<-function(node=NULL,support,
                             cols=c("white","black"),cex=1.5){
  scale<-range(support)
  support<-(support-scale[1])/diff(scale)
  pp<-get("last_plot.phylo",envir=.PlotPhyloEnv)
  if(is.null(node)) node<-1:pp$Nnode+pp$Ntip
  root_node <- pp$Ntip + 1 
  colfunc<-colorRamp(cols)
  node_cols<-colfunc(support)/255
  x<-pp$xx[node] 
  y<-pp$yy[node]
  for(i in 1:nrow(node_cols)){
    if(!is.na(support[i]) && node[i] != root_node)
      points(x[i],y[i],pch=21,,cex=0.9,
             bg=rgb(node_cols[i,1],node_cols[i,2],node_cols[i,3]))
  }
  invisible(scale)
}

tree <- tr8reroot

tree$node.label<- ifelse(tree$node.label == "1.00", "1", tree$node.label)

foo<-function(){
  par(fg="black")
  nodelabels(text=nodelabels_for_plot,  bg= "white" , frame="none", width = 0.2, adj=c(1.05, -0.1)) 
  minmax<-add_support_labels(support=EN_tr8)
  add.color.bar(leg=0.3*max(nodeHeights(tree)),
                cols=colorRampPalette(c("white","black"))(100),
                lims=minmax,title="EN",
                subtitle="",prompt=FALSE,x=0.01*par()$usr[2],
                y=0.02*par()$usr[4])
  
}
splits<- 0.576
tree <- compute.brlen(tree, method = "Grafen", power = 0.3) 
terminal_edges <- which(tree$edge[,2] <= length(tree$tip.label))
tree$edge.length[terminal_edges] <- tree$edge.length[terminal_edges]-0.18
nodelabels_for_plot <- paste(tree$node.label)

split.plotTree(tree,splits,file= "Original_ML_coverged_done.pdf", ftype="i",fn=foo,lwd=1, cex=0.7, oma = c(1, 1, 1, 1), mar = c(1, 1, 1, 1))





pdf("original_ML_mrb_models_converged_1.pdf", height=35, width=18)
plot.phylo(tr8reroot,use.edge.length = T, cex=1, label.offset=0.05)
nodelabels(tr8reroot$node.label,adj = c(0,0.4), bg="transparent", frame="none")
title("Original dataset ML converged genes")
dev.off()


#################################################################################
#plotting tr 9
tr9reroot$edge.length[is.na(tr9reroot$edge.length)] <- 0.2

#removing decimals from PP values
for (i in seq_along(tr9reroot$node.label)) {
  if(tr9reroot$node.label[i] == "Root" ) {next
  } else {
    tr9reroot$node.label[i] <- format(round(as.numeric(tr9reroot$node.label[i]), 2), nsmall = 2)
  }
}
tr9reroot$node.label<- ifelse(tr9reroot$node.label == "1.00", "1", tr9reroot$node.label)

EN_tr9 <-  as.numeric(sub(".*EN=(.*?)\\].*", "\\1",tr9_info$node.label))
EN_tr9[is.na(EN_tr9)] <- 0
length(EN_tr9) == tr9reroot$Nnode
matched_nodes_9 <- matchNodes(tr9_info, tr9reroot)
EN_tr9 <- round(EN_tr9, 0)


#from liam revell phytools blog
add_support_labels<-function(node=NULL,support,
                             cols=c("white","black"),cex=1.5){
  scale<-range(support)
  support<-(support-scale[1])/diff(scale)
  pp<-get("last_plot.phylo",envir=.PlotPhyloEnv)
  if(is.null(node)) node<-1:pp$Nnode+pp$Ntip
  root_node <- pp$Ntip + 1 
  colfunc<-colorRamp(cols)
  node_cols<-colfunc(support)/255
  x<-pp$xx[node] 
  y<-pp$yy[node]
  for(i in 1:nrow(node_cols)){
    if(!is.na(support[i]) && node[i] != root_node)
      points(x[i],y[i],pch=21,,cex=0.9,
             bg=rgb(node_cols[i,1],node_cols[i,2],node_cols[i,3]))
  }
  invisible(scale)
}

tree <- tr9reroot

tree$node.label<- ifelse(tree$node.label == "1.00", "1", tree$node.label)

foo<-function(){
  par(fg="black")
  nodelabels(text=nodelabels_for_plot,  bg= "white" , frame="none", width = 0.2, adj=c(1.05, -0.1)) 
  minmax<-add_support_labels(support=EN_tr9)
  add.color.bar(leg=0.3*max(nodeHeights(tree)),
                cols=colorRampPalette(c("white","black"))(100),
                lims=minmax,title="EN",
                subtitle="",prompt=FALSE,x=0.01*par()$usr[2],
                y=0.02*par()$usr[4])
  
}
splits<- 0.565
tree <- compute.brlen(tree, method = "Grafen", power = 0.3) 
terminal_edges <- which(tree$edge[,2] <= length(tree$tip.label))
tree$edge.length[terminal_edges] <- tree$edge.length[terminal_edges]-0.18
nodelabels_for_plot <- paste(tree$node.label)

split.plotTree(tree,splits,file= "Original_ML_freerate_done.pdf", ftype="i",fn=foo,lwd=1, cex=0.7, oma = c(1, 1, 1, 1), mar = c(1, 1, 1, 1))





pdf("original_ML_all_sub_models.pdf", height=35, width=18)
plot.phylo(tr9reroot,use.edge.length = TRUE, cex=1, label.offset=0.0005)
nodelabels(tr9reroot$node.label,adj = c(1,- 0.2), bg="transparent", frame="none")
title("Original dataset ML all substitution models")
dev.off()



#################################################################################
#plotting tr 10
#Adding terminal branchlegths
tr10reroot$edge.length[is.na(tr10reroot$edge.length)] <- 0.2


#moving decimals from LPP values for plotting
for (i in seq_along(tr10reroot$node.label)) {
  if(tr10reroot$node.label[i] == "Root" ) {next
  } else {
    tr10reroot$node.label[i] <- format(round(as.numeric(tr10reroot$node.label[i]), 2), nsmall = 2)
  }
}
tr10reroot$node.label<- ifelse(tr10reroot$node.label == "1.00", "1", tr10reroot$node.label)
#plotting the tree

EN_tr10 <-  as.numeric(sub(".*EN=(.*?)\\].*", "\\1",tr10_info$node.label))
EN_tr10[is.na(EN_tr10)] <- 0
length(EN_tr10) == tr10reroot$Nnode
matched_nodes_10 <- matchNodes(tr10_info, tr10reroot)
EN_tr10 <- round(EN_tr10, 0)


#from liam revell phytools blog
add_support_labels<-function(node=NULL,support,
                             cols=c("white","black"),cex=1.5){
  scale<-range(support)
  support<-(support-scale[1])/diff(scale)
  pp<-get("last_plot.phylo",envir=.PlotPhyloEnv)
  if(is.null(node)) node<-1:pp$Nnode+pp$Ntip
  root_node <- pp$Ntip + 1 
  colfunc<-colorRamp(cols)
  node_cols<-colfunc(support)/255
  x<-pp$xx[node] 
  y<-pp$yy[node]
  for(i in 1:nrow(node_cols)){
    if(!is.na(support[i]) && node[i] != root_node)
      points(x[i],y[i],pch=21,,cex=0.9,
             bg=rgb(node_cols[i,1],node_cols[i,2],node_cols[i,3]))
  }
  invisible(scale)
}


tree <- tr10reroot

tree$node.label<- ifelse(tree$node.label == "1.00", "1", tree$node.label)

foo<-function(){
  par(fg="black")
  nodelabels(text=nodelabels_for_plot,  bg= "white" , frame="none", width = 0.2, adj=c(1.05, -0.1)) 
  minmax<-add_support_labels(support=EN_tr10)
  add.color.bar(leg=0.3*max(nodeHeights(tree)),
                cols=colorRampPalette(c("white","black"))(100),
                lims=minmax,title="EN",
                subtitle="",prompt=FALSE,x=0.01*par()$usr[2],
                y=0.02*par()$usr[4])
  
}
splits<- 0.533
tree <- compute.brlen(tree, method = "Grafen", power = 0.3) 
terminal_edges <- which(tree$edge[,2] <= length(tree$tip.label))
tree$edge.length[terminal_edges] <- tree$edge.length[terminal_edges]-0.18
nodelabels_for_plot <- paste(tree$node.label)

split.plotTree(tree,splits,file= "large_ML_all_genes_done.pdf", ftype="i",fn=foo,lwd=1, cex=0.7, oma = c(1, 1, 1, 1), mar = c(1, 1, 1, 1))





pdf("large_ML_mrb_models_all_genes.pdf", height=35, width=20)
plot.phylo(tr10reroot,use.edge.length = TRUE, label.offset=0.0005)
nodelabels(tr10reroot$node.label,adj = c(1,- 0.2), bg="transparent", frame="none")
title("Large ML MrBayes substitution models with all genes")
dev.off()







################################################################################
#############################Move PP ###########################################
################################################################################

#move PP Large Bayesian consensus and Bayesian tree
#Check the tip names
tr4plot <- drop.tip(tr4reroot, c("D._album", "R._elegans", "R._augusta", "I._mirabilis", "I._polymorpha", "I._thalangensis", "I._namsabiensis" ,"L._rupicola"))
tr2plot <- drop.tip(tr2reroot, c("D._album", "R._elegans", "R._augusta", "I._mirabilis", "I._polymorpha", "I._thalangensis", "I._namsabiensis" ,"L._rupicola"))

tr4plot$tip.label <-gsub("_"," ",tr4plot$tip.label)
tr2plot$tip.label <-gsub("_"," ",tr2plot$tip.label)
match(tr4plot$tip.label, tr2plot$tip.label)

tr4_labels <- tr4plot$tip.label
tr2_labels <- tr2plot$tip.label
missing_in_tr4plot <- setdiff(tr2_labels, tr4_labels)
missing_in_tr4plot
missing_in_tr2plot <- setdiff(tr4_labels, tr2_labels)
missing_in_tr2plot  


biparts <- prop.clades(tr4plot, tr2plot)
biparts 

large_M<-matchNodes(tr4plot,tr2plot)
#in M tr4=tr1 and tr2=tr2
large_M

new_number_of_node_different <- colSums(is.na(large_M))

# list of nodes from tr2rerooted that match tr1rerooted
NN_tr2_NL=c()
for(i in 1:nrow(large_M)){
  if(!is.na(large_M[i,2])) NN_tr2_NL[length(NN_tr2_NL)+1]= large_M[i,2]
}
length(NN_tr2_NL)

#matrix where the nodes that match tr1 are combined with the correct PP value
tr2_NL=c()
for(i in 1:nrow(large_M)){
  if(!is.na(large_M[i,2])) { tr2_NL[length(tr2_NL)+1]= tr2plot$node.label[large_M[i,2] - Ntip(tr2plot)]
  }else{ tr2_NL[length(tr2_NL)+1]= "*"
  }
}
tr2_NL
large_M <- cbind(large_M, tr2_NL)
large_M <- cbind(large_M, tr4plot$node.label)



#difference between LPP and PP
#(difference = PP - LPP)

l_Difference <- data.frame(Nodetr4 = character(), Nodetr2 = character(), Diff = numeric(), stringsAsFactors = FALSE)

for (i in 1:nrow(large_M)) {
  if (!is.na(as.numeric(large_M[i,3]) - as.numeric(large_M[i,4]))) {
    l_Difference <- rbind(l_Difference, data.frame(
      Nodetr4 = large_M[i,1],
      Nodetr2 = large_M[i,2],
      Diff = as.numeric(large_M[i,3]) - as.numeric(large_M[i,4])
    ))
  }
}

l_Difference

PP_higher <- l_Difference$Diff[l_Difference$Diff > 0]
difference_max <- max(PP_higher)
difference_min <- min(PP_higher)


large_pp_bigger_than_LPP=0
large_LPP_equals_PP=0
large_LPP_bigger_than_PP=0
node_does_not_exist =0

for (i in 1:nrow(l_Difference)){
  if (l_Difference[i,3]< 0) {large_LPP_bigger_than_PP <- large_LPP_bigger_than_PP+1
  }else if (l_Difference[i,3] == 0){large_LPP_equals_PP <- large_LPP_equals_PP+1
  }else if (l_Difference[i,3] > 0) {large_pp_bigger_than_LPP <- large_pp_bigger_than_LPP+1
  }else {node_does_not_exist <- node_does_not_exist+1}
}
large_LPP_bigger_than_PP
large_LPP_equals_PP
large_pp_bigger_than_LPP
node_does_not_exist

EN_tr4 <-  as.numeric(sub(".*EN=(.*?)\\].*", "\\1",tr4_info$node.label))
EN_tr4[is.na(EN_tr4)] <- 0
length(EN_tr4) == tr4reroot$Nnode
matched_nodes_4 <- matchNodes(tr4_info, tr4reroot)
EN_tr4 <- round(EN_tr4, 0)


#from liam revell phytools blog
add_support_labels<-function(node=NULL,support,
                             cols=c("white","black"),cex=1.5){
  scale<-range(support)
  support<-(support-scale[1])/diff(scale)
  pp<-get("last_plot.phylo",envir=.PlotPhyloEnv)
  if(is.null(node)) node<-1:pp$Nnode+pp$Ntip
  root_node <- pp$Ntip + 1 
  colfunc<-colorRamp(cols)
  node_cols<-colfunc(support)/255
  x<-pp$xx[node] 
  y<-pp$yy[node]
  for(i in 1:nrow(node_cols)){
    if(!is.na(support[i]))
      points(x[i],y[i],pch=21,,cex=0.9,
             bg=rgb(node_cols[i,1],node_cols[i,2],node_cols[i,3]))
  }
  invisible(scale)
}

#from liam revell phytools blog
split.plotTree<-function(tree,splits=NULL,file=NULL,fn,...){
  ef<-0.037037037037
  if(!is.null(file)) pdf(file, height=15, width=10)
  if(is.null(splits)) splits<-(floor(0.5*Ntip(tree))+0.5)/Ntip(tree)
  S<-matrix(c(0,splits,splits,1+1/Ntip(tree)),length(splits)+1,2)
  S<-cbind(S[,1]+ef*(S[,2]-S[,1]),S[,2]-ef*(S[,2]-S[,1]))
  for(i in nrow(S):1){
    if(is.null(file)&&i<nrow(S)) par(ask=TRUE)
    plotTree(tree,ylim=Ntip(tree)*S[i,],use.edge.length=FALSE,...)
    fn()
  }
  if(!is.null(file)) oo<-dev.off()
}

foo<-function(){
  par(fg="black")
  nodelabels(text=nodelabels_for_plot,  bg= "white" , frame="none", width = 0.2, adj=c(1.05, -0.1), cex=0.8)
  minmax<-add_support_labels(support=EN_tr4)
  add.color.bar(leg=0.3*max(nodeHeights(tree)),
                cols=colorRampPalette(c("white","black"))(100),
                lims=minmax,title="EN",
                subtitle="",prompt=FALSE,x=0.01*par()$usr[2],
                y=0.02*par()$usr[4])
}

splits<- 0.556
tree <- tr4plot

#removing decimals from PP values
for (i in seq_along(tree$node.label)) {
  if(tree$node.label[i] == "Root" ) {next
  } else {
    tree$node.label[i] <- format(round(as.numeric(tree$node.label[i]), 2), nsmall = 2)
  }
}
tree$node.label<- ifelse(tree$node.label == "1.00", "1", tree$node.label)

tr2_NL <- format(round(as.numeric(tr2_NL), 2), nsmall = 2)

tr2_NL <- ifelse(tr2_NL == "1.00", "1", tr2_NL)
tr2_NL <- ifelse(tr2_NL =="  NA", "*", tr2_NL)

nodelabels_for_plot <- paste(tree$node.label,"/", tr2_NL)
nodelabels_for_plot <- ifelse(nodelabels_for_plot =="1 / 1", "", nodelabels_for_plot)
tree <- compute.brlen(tree, method = "Grafen", power = 0.3) 
terminal_edges <- which(tree$edge[,2] <= length(tree$tip.label))
tree$edge.length[terminal_edges] <- tree$edge.length[terminal_edges]-0.15

split.plotTree(tree,splits,file= "large_PP_move_MrB_51.pdf", ftype="i",fn=foo,lwd=1, cex=0.7, oma = c(1, 1, 1, 1), mar = c(1, 1, 1, 1), use.edge.lenght=F)


################################################################################
#move PP Original MrB consensus and MrB tree

tr7plot <-drop.tip(tr7reroot, c("L._rupicola_1", "L._rupicola_2"))
tr1plot <-drop.tip(tr1reroot, c("L._rupicola_1", "L._rupicola_2"))

#Check the tip names
tr7plot$tip.label <-gsub("_"," ",tr7plot$tip.label)
tr1plot$tip.label <-gsub("_"," ",tr1plot$tip.label)
match(tr7plot$tip.label, tr1plot$tip.label)

tr7_labels <- tr7plot$tip.label
tr1_labels <- tr1plot$tip.label
missing_in_tr7 <- setdiff(tr1_labels, tr7_labels)
missing_in_tr7
missing_in_tr1 <- setdiff(tr7_labels, tr1_labels)
missing_in_tr1  

M<-matchNodes(tr7plot,tr1plot)
# in M tr1=tr7 and tr2=tr1
M

number_of_node_different <- colSums(is.na(M))

# list of nodes from tr1rerooted that match tr7rerooted
NN_tr1_NL=c()
for(i in 1:nrow(M)){
  if(!is.na(M[i,2])) NN_tr1_NL[length(NN_tr1_NL)+1]= M[i,2]
}
length(NN_tr1_NL)

#matrix where the nodes that match tr1 are combined with the correct PP value
tr1_NL=c()
for(i in 1:nrow(M)){
  if(!is.na(M[i,2])) { tr1_NL[length(tr1_NL)+1]= tr1plot$node.label[M[i,2] - Ntip(tr1plot)]
  }else{ tr1_NL[length(tr1_NL)+1]= "*"
  }
}
tr1_NL
M <- cbind(M, tr1_NL)


tr7_NL=c()
tr7_NL <- tr7plot$node.label
M <- cbind(M, tr7_NL)
M                    

#difference between LPP and PP
#Difference= PP - LPP 
Difference <- data.frame(Nodetr7 = character(), Nodetr1 = character(), Diff = numeric(), stringsAsFactors = FALSE)
for (i in 1:nrow(M)) {
  if (!is.na(as.numeric(M[i,3]) - as.numeric(M[i,4]))) {
    Difference <- rbind(Difference, data.frame(
      Nodetr7 = M[i,1],
      Nodetr1 = M[i,2],
      Diff = as.numeric(M[i,3]) - as.numeric(M[i,4])
    ))
  }
}

Difference

pp_bigger_than_LPP=0
LPP_equals_PP=0
LPP_bigger_than_PP=0
does_not_exist =0

for (i in 1:nrow(Difference)){
  if (Difference[i,3]< 0) {LPP_bigger_than_PP <- LPP_bigger_than_PP+1
  }else if (Difference[i,3] == 0){LPP_equals_PP <- LPP_equals_PP+1
  }else if (Difference[i,3] > 0) {pp_bigger_than_LPP <- pp_bigger_than_LPP+1
  }else {node_does_not_exist <- node_does_not_exist+1}
}
LPP_bigger_than_PP
LPP_equals_PP
pp_bigger_than_LPP
does_not_exist

PP_higher <- Difference$Diff[Difference$Diff > 0]
difference_max <- max(PP_higher)
difference_min <- min(PP_higher)

EN_tr7 <-  as.numeric(sub(".*EN=(.*?)\\].*", "\\1",tr7_info$node.label))
EN_tr7[is.na(EN_tr7)] <- 0
length(EN_tr7) == tr7reroot$Nnode
matched_nodes_7 <- matchNodes(tr7_info, tr7reroot)
EN_tr7 <- round(EN_tr7, 0)


#from liam revell phytools blog
add_support_labels<-function(node=NULL,support,
                             cols=c("white","black"),cex=1.5){
  scale<-range(support)
  support<-(support-scale[1])/diff(scale)
  pp<-get("last_plot.phylo",envir=.PlotPhyloEnv)
  if(is.null(node)) node<-1:pp$Nnode+pp$Ntip
  root_node <- pp$Ntip + 1 
  colfunc<-colorRamp(cols)
  node_cols<-colfunc(support)/255
  x<-pp$xx[node] 
  y<-pp$yy[node]
  for(i in 1:nrow(node_cols)){
    if(!is.na(support[i]))
      points(x[i],y[i],pch=21,,cex=0.9,
             bg=rgb(node_cols[i,1],node_cols[i,2],node_cols[i,3]))
  }
  invisible(scale)
}
#Making split plot
#from liam revell phytools blog
split.plotTree<-function(tree,splits=NULL,file=NULL,fn,...){
  ef<-0.037037037037
  if(!is.null(file)) pdf(file, height=15, width=14)
  if(is.null(splits)) splits<-(floor(0.5*Ntip(tree))+0.5)/Ntip(tree)
  S<-matrix(c(0,splits,splits,1+1/Ntip(tree)),length(splits)+1,2)
  S<-cbind(S[,1]+ef*(S[,2]-S[,1]),S[,2]-ef*(S[,2]-S[,1]))
  for(i in nrow(S):1){
    if(is.null(file)&&i<nrow(S)) par(ask=TRUE)
    plotTree(tree,ylim=Ntip(tree)*S[i,],...)
    fn()
  }
  if(!is.null(file)) oo<-dev.off()
}
foo<-function(){
  par(fg="black")
  nodelabels(text=nodelabels_for_plot,  bg="transparent", frame="none", width = 0.2, adj=c(1.05, -0.1), cex=0.8) 
  minmax<-add_support_labels(support=EN_tr7)
  add.color.bar(leg=0.3*max(nodeHeights(tree)),
                cols=colorRampPalette(c("white","black"))(100),
                lims=minmax,title="EN",
                subtitle="",prompt=FALSE,x=0.01*par()$usr[2],
                y=0.02*par()$usr[4])
}

splits<- 0.57
tree <- tr7plot
#removing decimals from PP values
for (i in seq_along(tree$node.label)) {
  if(tree$node.label[i] == "Root" ) {next
  } else {
    tree$node.label[i] <- format(round(as.numeric(tree$node.label[i]), 2), nsmall = 2)
  }
}
tree$node.label<- ifelse(tree$node.label == "1.00", "1", tree$node.label)

tr1_NL <- format(round(as.numeric(tr1_NL), 2), nsmall = 2)

tr1_NL <- ifelse(tr1_NL == "1.00", "1", tr1_NL)
tr1_NL <- ifelse(tr1_NL =="  NA", "*", tr1_NL)

nodelabels_for_plot <- paste(tree$node.label,"/", tr1_NL)
nodelabels_for_plot <- ifelse(nodelabels_for_plot =="1 / 1", "", nodelabels_for_plot)
tree <- compute.brlen(tree, method = "Grafen", power = 0.3) 
terminal_edges <- which(tree$edge[,2] <= length(tree$tip.label))
tree$edge.length[terminal_edges] <- tree$edge.length[terminal_edges]-0.18
#old
#tree$edge.length[is.na(tree$edge.length)] <- 99
#tree$edge.length <- ifelse(tree$edge.length < 1, tree$edge.length+1.3, tree$edge.length)
#tree$edge.length <- ifelse(tree$edge.length ==99, 0.1, tree$edge.length)
#tree$edge.length <- ifelse(tree$edge.length > 2, tree$edge.length/2, tree$edge.length)

split.plotTree(tree,splits,file= "PP_move_MrB_Original_51.pdf", ftype="i",fn=foo,lwd=1, cex=0.5, x.lim = c(0, max(branching.times(tree)) + 0.6))



#making regular pdf plot
pdf("Heyduk_MrB_with_PP_from_cons.pdf", height=30, width=25)
plot(tr7reroot, use.edge.length = T)
nodelabels(text=paste(tr7reroot$node.label, "/", M[,3], sep=" "),  bg="transparent", frame="none", width = 0.2, adj=c(1.05, -0.5), cex=0.8) 
title("Original dataset Bayesian tree with PP from Bayesian consensus tre")
dev.off()


################################################################################
#################### Making split plot #########################################
################################################################################

# Split plot Large MrB tree
#Making split plot
#from liam revell phytools blog
split.plotTree<-function(tree,splits=NULL,file=NULL,fn,...){
  ef<-0.037037037037
  if(!is.null(file)) pdf(file,width=30,height=30)
  if(is.null(splits)) splits<-(floor(0.5*Ntip(tree))+0.5)/Ntip(tree)
  S<-matrix(c(0,splits,splits,1+1/Ntip(tree)),length(splits)+1,2)
  S<-cbind(S[,1]+ef*(S[,2]-S[,1]),S[,2]-ef*(S[,2]-S[,1]))
  for(i in nrow(S):1){
    if(is.null(file)&&i<nrow(S)) par(ask=TRUE)
    plotTree(tree,ylim=Ntip(tree)*S[i,],...)
    fn()
  }
  if(!is.null(file)) oo<-dev.off()
}
foo<-function(){
  nodelabels(text=paste(tr4reroot$node.label, "/", M[,3], sep=" "),  bg="transparent", frame="none", width = 0.2, adj=c(1.05, -0.5)) 
  
}

splits<- 0.515
tree <- tr4reroot
for (i in seq_along(tree$node.label)) {
  if(tree$node.label[i] == "Root" ) {next
  } else {
    tree$node.label[i] <- format(round(as.numeric(tree$node.label[i]), 2), nsmall = 2)
  }
}
tree$node.label<- ifelse(tree$node.label == "1.00", "1", tree$node.label)
split.plotTree(tree,splits,file= "split_PP_move_MrB_15.pdf",ftype="i",mar=rep(1.1,4),fn=foo,lwd=1)


dev.off()





################################################################################
################# Make plot with Baitkit information ###########################
################################################################################

import data - metadata
m_data<-read.csv(file="My_dataset_info.csv",sep=";")
head(m_data)
m_data$Heyduk <- toupper(m_data$Heyduk)
m_data$PopcornPalms <- toupper(m_data$PopcornPalms)
m_data$PhyloPALMS <- toupper(m_data$PhyloPALMS)
# Convert the column to UTF-8 encoding (this should handle invalid characters)
m_data$Harmonization.11.May.2023 <- iconv(m_data$Harmonization.11.May.2023, from = "latin1", to = "UTF-8")
# Replace left and right single quotes with a regular apostrophe
m_data$Harmonization.11.May.2023 <- gsub("\u0091", "‘", m_data$Harmonization.11.May.2023)
m_data$Harmonization.11.May.2023 <- gsub("\u0092", "’", m_data$Harmonization.11.May.2023)


for (i in 1:nrow(m_data)){
  if (m_data$Heyduk[i] == "X"){ m_data$Baitkit[i] <- "HE"
  }else if (m_data$PopcornPalms[i] == "X"){m_data$Baitkit[i] <- "PO"
  }else if (m_data$PhyloPALMS[i] == "X") {m_data$Baitkit[i] <- "PP"
  } else {m_data$Baitkit[i] <- NA
  }
}
m_data$Baitkit
m_data$Harmonization.11.May.2023

species_names <- tr4reroot$tip.label

baitkit= c()

for (i in species_names) {
  m_data_indx <- grep(i, m_data[,10], fixed = TRUE) 
  if (length(m_data_indx) > 0) {
    baitkit[[i]] <- m_data[m_data_indx[1], 29]
  } else {
    baitkit[[i]] <- NA
  }
}

print(baitkit)
#from phytools blog:
## custom function I'm going to use for the box labels
boxlabel<-function(x,y,text,cex=1,bg="transparent",offset=0){
  w<-strwidth(text)*cex*1.1
  h<-strheight(text)*cex*1.4
  os<-offset*strwidth("W")*cex
  rect(x+os,y-0.5*h,x+w+os,y+0.5*h,col=bg,border=0)
  text(x,y,text,pos=4,offset=offset,font=3)
}

baitkit <- unlist(baitkit)  # Convert list to character vector

cols<-setNames(RColorBrewer::brewer.pal(length(unique(baitkit)),"Accent"),sort(unique(baitkit)))
cols

Gene_recov= c()
for (i in species_names) {
  m_data_indx <- grep(i, m_data[,10], fixed = TRUE) 
  if (length(m_data_indx) > 0) {
    Gene_recov[[i]] <- as.numeric(m_data[m_data_indx[1], 14])
  } else {
    Gene_recov[[i]] <- NA
  }
}
Gene_recov_df <- data.frame(Species = names(Gene_recov), recov_after_trim = unlist(Gene_recov), stringsAsFactors = FALSE)

tr4reroot$tip.label <- ifelse(Gene_recov_df$recov_after_trim < 500, 
                              paste0(tr4reroot$tip.label, " *"), 
                              tr4reroot$tip.label)

pdf("MrB_with_PP_baitkit_7.pdf", height=70, width=35)
plot(tr4reroot, use.edge.length = T, tip.label=F)
new_nodelabels <- nodelabels(text=paste(tr4reroot$node.label, "/PP", M[,3], sep=" "), cex=0.5, horiz =0.6)
pp<-get("last_plot.phylo",envir=.PlotPhyloEnv)
N<-Ntip(tr4reroot)
par(fg="black")
for(i in 1:Ntip(tr4reroot)) boxlabel(pp$xx[i],pp$yy[i],tr4reroot$tip.label[i],bg=cols[baitkit[i]])
legend(x=5.7,y=26,legend = c("HE", "PO", "PP"), col=cols, fill = cols)

dev.off()


################################################################################
############# plot coloured clades #############################################
################################################################################
#Making clades coloured 
bl_tr3reroot <-tr3reroot
bl_tr3reroot$edge.length <-bl_tr3reroot$edge.length[is.na(bl_tr3reroot$edge.length)] <- 1

root_node <- getMRCA(bl_tr1reroot, c("Loxococcus_rupicola_1", "Loxococcus_rupicola_2"))

painted<-paintSubTree(bl_tr3reroot,184,"outgroup","0")
painted<-paintSubTree(painted,185,"Masoala")
painted<-paintSubTree(painted,186,"Lemurophoenix")
painted<-paintSubTree(painted,187,"Vonitra")
painted<-paintSubTree(painted,188,"Chrysadliocarpus")
painted<-paintSubTree(painted,189,"Marojejya")

painted$maps

pdf("ML_with_PP_colors.pdf", height=60, width=20)
plotSimmap(painted, ftype="off", lwd=3, ylim=c(0, 200), xlim=c(0, 20))
cladelabels(tr3reroot,c("Outgroup","Masoala","Lemurophoenix", "Vonitra", "Chrysadliocarpus", "Marojejya"),
            c(359,185,348,349,282, 189),wing.length=0,offset=0.5)
tiplabels(painted$tip.label, frame = "n", adj = 1, cex = 0.4)
nodelabels(text=paste(painted$node.label, "/", M[,3], sep=" "), cex=0.4, adj = c(0.8,- 0.2), bg="transparent")
dev.off()



################################################################################
#calculating genera support avarages for the original dataset

original_dypsis <- c("D._lantzeana_1", "D._lantzeana_2",  "D._forficifolia" , "D._andapae"  ,       
"D._curtisii" ,        "D._occidentalis",     "D._humilis"   ,       "D._spicata"   ,      
"D._pachyramea",       "D._eriostachys"  ,    "D._mahia"      ,      "D._singularis" ,     
"D._laevis"     ,      "D._angusta"       ,   "D._thermarum"   ,     "D._hildebrandtii" ,  
"D._catatiana"   ,     "D._louvelii"       ,  "D._turkii"       ,    "D._angustifolia"   , 
"D._ambilaensis"  ,    "D._lanuginosa"      , "D._bosseri"       ,   "D._linearis"  ,      
"D._delicatula"    ,   "D._beentjei"    ,     "D._glabrescens"    ,  "D._viridis"    ,     
"D._ramentacea"     ,  "D._mocquerysiana",    "D._humbertii"       , "D._cookei"      ,    
 "D._gautieri"       ,  "D._rivularis"  ,
 "D._pustulata"       , "D._subacaulis"  ,     "D._brevicaulis" ,     "D._integra"    ,     
"D._intermedia" ,      "D._culminis"      ,   "D._aurantiaca"    ,   "D._simianensis"  ,   
"D._corniculata" ,     "D._gronophyllum"   ,  "D._fasciculata"    ,  "D._minuta"  ,        
"D._thiryana"     ,    "D._lutea"           , "D._heterophylla"    , "D._coriacea" ,       
"D._bonsai"        ,   "D._schatzii"    ,     "D._poivreana"  ,      "D._concinna"  ,      
"D._procumbens"     ,  "D._dracaenoides" ,    "D._bernieriana" ,     "D._digitata"   ,     
"D._tenuissima"      , "D._elegans"       ,   "D._henrici"      ,    "D._commersoniana",   
"D._scottiana"  ,      "D._lilacina"       ,  "D._marojejyi"     ,   "D._coursii"       ,  
"D._aff._coursii",     "D._aquatilis"       , "D._montana"        ,  "D._brittiana"      , 
"D._rakotonasoloi",    "D._betsimisarakae",   "D._nodifera" ,        "D._sp."             ,
"D._paludosa"      ,   "D._lokohoensis"    ,  "D._hiarakae"  ,       "D._fanjana"         ,
"D._faneva"         ,  "D._interrupta"      , "D._jeremiei"   ,      "D._pinnatifrons"    ,
"D._rosea"           , "D._sancta"    ,       "D._scandens"    ,     "D._mirabilis"       ,
"D._jumelleana",       "D._remotiflora",      "D._boiviniana"   ,    "D._pervillei"       ,
"D._makirae"    ,      "D._procera"     ,     "D._caudata"       ,   "D._metallica"       ,
"D._confusa"     ,     "D._reflexa"      ,    "D._acaulis")
 
original_vonitra <- c("V._utilis", "V._crinita", "V._fibrosa","V._vonitrandambo", "V._pusilla", "V._antanambensis",
"V._moorei","V._nossibensis","V._perrieri", "V._dransfieldii")

original_Masoala <- c("M._madagascariensis","M._kona")

original_Lemurophoenix <- c("L._halleuxii", "L._laevis")

original_Marojejya <- c("M._insignis_2","M._insignis_1", "M._sp.", "M._darianii")

 original_chrysadliocarpus <- c("C._mijoroanus","C._lastellianus","C._nauseosus","C._leptocheilos",   
"C._titan",  "C._bejofo", "C._canaliculatus", "C._andrianatonga",
"C._ambositrae_2","C._ambositrae_1","C._ambanjae","C._onilahensis",
"C._pumilus","C._acuminum",  "C._sp._4","C._ovobontsira",
"C._leucomallus","C._oropedionis","C._sp._2","C._pilulifer",
"C._mananjarensis","C._malcomberi" ,"C._sp._5","C._ovojavavy",
"C._sp._3","C._sanctaemariae", "C._psammophilus" ,"C._albofarinosus",
"C._gracilis" ,"C._ankirindro", "C._serpentinus", "C._tsaravoasira",
"C._tanalensis" ,"C._sp._1","C._baronii", "C._decipiens",
"C._sp._6","C._carlsmithii", "C._arenarum","C._rabepierrei",
"C._lutescens","C._decaryi_2","C._decaryi_1", "C._heteromorphus",
"C._plumosus","C._tokoravina","C._robustus","C._prestonianus",
"C._madagascariensis", "C._pembanus","C._humblotianus","C._cabadae",
"C._lanceolatus","C._basilongus","C._ceraceus","C._saintelucei")

################################################################################
##################### calculating support averages##############################
################################################################################
 #calculating tr1
 vonitra_clade_node=getMRCA(tr1reroot, original_vonitra)
 vonitra <- Descendants(tr1reroot, node=vonitra_clade_node, "all")
 vonitra_tips <- Descendants(tr1reroot, node=vonitra_clade_node, "tips")
 vonitra_nodes <- c(vonitra[!(vonitra %in% vonitra_tips[[1]])], vonitra_clade_node)
 LPP_vonitra_1 <- sum(as.numeric(tr1reroot$node.label[vonitra_nodes-Ntip(tr1reroot)]))/length(vonitra_nodes)
 LPP_vonitra_1
 
 chrys_clade_node=getMRCA(tr1reroot, original_chrysadliocarpus)
 chrys <- Descendants(tr1reroot, node=chrys_clade_node, "all")
 chrys_tips <- Descendants(tr1reroot, node=chrys_clade_node, "tips")
 chrys_nodes <- c(chrys[!(chrys %in% chrys_tips[[1]])], chrys_clade_node)
 LPP_chrys_1 <- sum(as.numeric(tr1reroot$node.label[chrys_nodes-Ntip(tr1reroot)]))/length(chrys_nodes)
 LPP_chrys_1
 
 dypsis_clade_node=getMRCA(tr1reroot, original_dypsis)
 dypsis <- Descendants(tr1reroot, node=dypsis_clade_node, "all")
 dypsis 
 dypsis_tips <- Descendants(tr1reroot, node=dypsis_clade_node, "tips")
 dypsis_tips
 dypsis_nodes <- c(dypsis[!(dypsis %in% dypsis_tips[[1]])], dypsis_clade_node)
 dypsis_nodes
 LPP_dypsis_1 <- sum(as.numeric(tr1reroot$node.label[dypsis_nodes-Ntip(tr1reroot)]),na.rm = TRUE)/length(dypsis_nodes)
 LPP_dypsis_1
 
#calculating t6
vonitra_clade_node=getMRCA(tr6reroot, original_vonitra)
vonitra <- Descendants(tr6reroot, node=vonitra_clade_node, "all")
vonitra 
vonitra_tips <- Descendants(tr6reroot, node=vonitra_clade_node, "tips")
vonitra_tips
vonitra_nodes <- c(vonitra[!(vonitra %in% vonitra_tips[[1]])], vonitra_clade_node)
vonitra_nodes
LPP_vonitra_6 <- sum(as.numeric(tr6reroot$node.label[vonitra_nodes-Ntip(tr6reroot)]))/length(vonitra_nodes)
LPP_vonitra_6

chrys_clade_node=getMRCA(tr6reroot, original_chrysadliocarpus)
  chrys <- Descendants(tr6reroot, node=chrys_clade_node, "all")
chrys 
chrys_tips <- Descendants(tr6reroot, node=chrys_clade_node, "tips")
chrys_tips
chrys_nodes <- c(chrys[!(chrys %in% chrys_tips[[1]])], chrys_clade_node)
chrys_nodes
LPP_chrys_6 <- sum(as.numeric(tr6reroot$node.label[chrys_nodes-Ntip(tr6reroot)]))/length(chrys_nodes)
LPP_chrys_6

dypsis_clade_node=getMRCA(tr6reroot, original_dypsis)
  dypsis <- Descendants(tr6reroot, node=dypsis_clade_node, "all")
dypsis 
dypsis_tips <- Descendants(tr6reroot, node=dypsis_clade_node, "tips")
dypsis_tips
dypsis_nodes <- c(dypsis[!(dypsis %in% dypsis_tips[[1]])], dypsis_clade_node)
dypsis_nodes
LPP_dypsis_6 <- sum(as.numeric(tr6reroot$node.label[dypsis_nodes-Ntip(tr6reroot)]),na.rm = TRUE)/length(dypsis_nodes)
LPP_dypsis_6

#calculating tr7
vonitra_clade_node=getMRCA(tr7reroot, original_vonitra)
vonitra <- Descendants(tr7reroot, node=vonitra_clade_node, "all")
vonitra 
vonitra_tips <- Descendants(tr7reroot, node=vonitra_clade_node, "tips")
vonitra_tips
vonitra_nodes <- c(vonitra[!(vonitra %in% vonitra_tips[[1]])], vonitra_clade_node)
vonitra_nodes
LPP_vonitra_7 <- sum(as.numeric(tr7reroot$node.label[vonitra_nodes-Ntip(tr7reroot)]))/length(vonitra_nodes)
PP_vonitra_7

chrys_clade_node=getMRCA(tr7reroot, original_chrysadliocarpus)
chrys <- Descendants(tr7reroot, node=chrys_clade_node, "all")
chrys 
chrys_tips <- Descendants(tr7reroot, node=chrys_clade_node, "tips")
chrys_tips
chrys_nodes <- c(chrys[!(chrys %in% chrys_tips[[1]])], chrys_clade_node)
chrys_nodes
LPP_chrys_7 <- sum(as.numeric(tr7reroot$node.label[chrys_nodes-Ntip(tr7reroot)]))/length(chrys_nodes)
LPP_chrys_7

dypsis_clade_node=getMRCA(tr7reroot, original_dypsis)
dypsis <- Descendants(tr7reroot, node=dypsis_clade_node, "all")
dypsis 
dypsis_tips <- Descendants(tr7reroot, node=dypsis_clade_node, "tips")
dypsis_tips
dypsis_nodes <- c(dypsis[!(dypsis %in% dypsis_tips[[1]])], dypsis_clade_node)
dypsis_nodes
LPP_dypsis_7 <- sum(as.numeric(tr7reroot$node.label[dypsis_nodes-Ntip(tr7reroot)]),na.rm = TRUE)/length(dypsis_nodes)
LPP_dypsis_7


#calculating tr8
vonitra_clade_node=getMRCA(tr8reroot, original_vonitra)
vonitra <- Descendants(tr8reroot, node=vonitra_clade_node, "all")
vonitra 
vonitra_tips <- Descendants(tr8reroot, node=vonitra_clade_node, "tips")
vonitra_tips
vonitra_nodes <- c(vonitra[!(vonitra %in% vonitra_tips[[1]])], vonitra_clade_node)
vonitra_nodes
LPP_vonitra_8 <- sum(as.numeric(tr8reroot$node.label[vonitra_nodes-Ntip(tr8reroot)]))/length(vonitra_nodes)
LPP_vonitra_8

chrys_clade_node=getMRCA(tr8reroot, original_chrysadliocarpus)
chrys <- Descendants(tr8reroot, node=chrys_clade_node, "all")
chrys 
chrys_tips <- Descendants(tr8reroot, node=chrys_clade_node, "tips")
chrys_tips
chrys_nodes <- c(chrys[!(chrys %in% chrys_tips[[1]])], chrys_clade_node)
chrys_nodes
LPP_chrys_8 <- sum(as.numeric(tr8reroot$node.label[chrys_nodes-Ntip(tr8reroot)]))/length(chrys_nodes)
LPP_chrys_8

dypsis_clade_node=getMRCA(tr8reroot, original_dypsis)
dypsis <- Descendants(tr8reroot, node=dypsis_clade_node, "all")
dypsis 
dypsis_tips <- Descendants(tr8reroot, node=dypsis_clade_node, "tips")
dypsis_tips
dypsis_nodes <- c(dypsis[!(dypsis %in% dypsis_tips[[1]])], dypsis_clade_node)
dypsis_nodes
LPP_dypsis_8 <- sum(as.numeric(tr8reroot$node.label[dypsis_nodes-Ntip(tr8reroot)]),na.rm = TRUE)/length(dypsis_nodes)
LPP_dypsis_8


#calculating tr9
vonitra_clade_node=getMRCA(tr9reroot, original_vonitra)
vonitra <- Descendants(tr9reroot, node=vonitra_clade_node, "all")
vonitra 
vonitra_tips <- Descendants(tr9reroot, node=vonitra_clade_node, "tips")
vonitra_tips
vonitra_nodes <- c(vonitra[!(vonitra %in% vonitra_tips[[1]])], vonitra_clade_node)
vonitra_nodes
LPP_vonitra_9 <- sum(as.numeric(tr9reroot$node.label[vonitra_nodes-Ntip(tr9reroot)]))/length(vonitra_nodes)
LPP_vonitra_9

chrys_clade_node=getMRCA(tr9reroot, original_chrysadliocarpus)
chrys <- Descendants(tr9reroot, node=chrys_clade_node, "all")
chrys 
chrys_tips <- Descendants(tr9reroot, node=chrys_clade_node, "tips")
chrys_tips
chrys_nodes <- c(chrys[!(chrys %in% chrys_tips[[1]])], chrys_clade_node)
chrys_nodes
LPP_chrys_9 <- sum(as.numeric(tr9reroot$node.label[chrys_nodes-Ntip(tr9reroot)]))/length(chrys_nodes)
LPP_chrys_9

dypsis_clade_node=getMRCA(tr9reroot, original_dypsis)
dypsis <- Descendants(tr9reroot, node=dypsis_clade_node, "all")
dypsis 
dypsis_tips <- Descendants(tr9reroot, node=dypsis_clade_node, "tips")
dypsis_tips
dypsis_nodes <- c(dypsis[!(dypsis %in% dypsis_tips[[1]])], dypsis_clade_node)
dypsis_nodes
LPP_dypsis_9 <- sum(as.numeric(tr9reroot$node.label[dypsis_nodes-Ntip(tr9reroot)]),na.rm = TRUE)/length(dypsis_nodes)
LPP_dypsis_9


################################################################################
# calculating support averages new

new_dypsis <- c("D._brevicaulis"  , "D._subacaulis" ,    "D._integra"  , "D._intermedia", "D._culminis", "D._aurantiaca"  ,   "D._simianensis", "D._gronophyllum"   , "D._fasciculata" ,   "D._schatzii"  ,   "D._sp._2", "D._sp._3"           
, "D._poivreana"  ,      "D._sp._4"  ,          "D._heterophylla"  ,   "D._coriacea",        
 "D._bonsai" ,          "D._minuta"  ,         "D._thiryana"    ,     "D._lutea" ,          
 "D._procumbens"     ,  "D._concinna"       ,  "D._bernieriana"   ,   "D._dracaenoides" ,   
 "D._digitata"   ,      "D._tenuissima"   ,    "D._elegans"    ,      "D._scottiana"  ,     
 "D._lilacina"  ,       "D._henrici"  ,        "D._commersoniana"  ,  "D._corniculata" ,     
 "D._mcdonaldiana"  ,   "D._thermarum"   ,     "D._angusta"   ,       "D._hildebrandtii",
"D._catatiana"  ,      "D._ambilaensis"  ,    "D._laevis"   ,        "D._singularis"  ,    
 "D._mahia"   ,         "D._eriostachys"  ,    "D._louvelii"  ,       "D._turkii",          
 "D._angustifolia" ,    "D._glabrescens" ,     "D._beentjei"     ,    "D._delicatula"  ,    
"D._viridis"        ,  "D._bosseri"       ,   "D._lanuginosa"     ,  "D._linearis"  ,      
"D._mocquerysiana"   , "D._ramentacea"     ,  "D._spicata"         , "D._pachyramea" ,     
"D._humbertii",        "D._forficifolia"    , "D._humilis"          ,"D._acaulis"     ,    
"D._sahanofensis",     "D._andapae",          "D._occidentalis" ,    "D._curtisii"     ,   
"D._gautieri"     ,    "D._cookei"  ,         "D._rivularis"     ,   "D._confusa"       ,  
"D._metallica"     ,   "D._caudata"  ,        "D._reflexa"        ,  "D._boiviniana"     , 
 "D._pervillei"     ,   "D._makirae"  ,        "D._remotiflora"    ,  "D._hiarakae"        ,
 "D._fanjana"        ,  "D._faneva"    ,       "D._lokohoensis"     , "D._paludosa"        ,
"D._sp.", "D._nodifera"   ,"D._betsimisarakae", "D._interrupta",       "D._scandens"    , "D._mirabilis" , "D._jumelleana",
"D._pinnatifrons",     "D._rosea"         ,   "D._sancta"       ,    "D._rakotonasoloi"   ,
"D._jeremiei"     ,    "D._brittiana"      ,  "D._andilamenensis",   "D._aquatilis"       ,
"D._coursii"       ,   "D._marojejyi"       , "D._pustulata"  )

new_Mayojejya <- c(  "M._darianii"  , "M._insignis" )

new_Chrysadliocarpus <- c( "C._mananjarensis",    "C._malcomberi" ,      "C._sp._2"       ,    
"C._oropedionis" ,     "C._ovojavavy"      ,  "C._pilulifer"    ,    "C._hovomantsina"  ,  
"C._sp._1"        ,    "C._sanctaemariae"   , "C._ifanadianae"   ,   "C._sp._3"           , 
"C._psammophilus"  ,   "C._albofarinosus"    ,"C._ovobontsira"    ,  "C._loucoubensis"    ,
"C._sp._14"         ,  "C._sp._9"  ,          "C._tanalensis"      , "C._tsaravoasira"    ,
"C._ankirindro"      , "C._gracilis",         "C._serpentinus",      "C._hankona"    ,     
"C._decipiens",        "C._sp._7"    ,        "C._sp._6"       ,     "C._onilahensis" ,    
"C._ambanjae"  ,       "C._ambositrae",       "C._andrianatonga",    "C._pumilus",         
"C._acuminum"   ,      "C.sp._4"       ,"C._sp._8"    ,        
"C._arenarum"    ,     "C._baronii"     ,     "C._lutescens"      ,  "C._carlsmithii",     
"C._cabadae"      ,    "C._sp._13"       ,    "C._lanceolatus"     , "C._pembanus"    ,    
"C._humblotianus"  ,   "C._madagascariensis", "C._sp._12"    ,       "C._prestonianus" ,   
"C._robustus"       ,  "C._tokoravina",       "C._rufescens"  ,      "C._decaryi"  ,       
"C._plumosus"   ,      "C._sp._10"     ,      "C._ceraceus"    ,     "C._basilongus",      
"C._saintelucei" ,     "C._mijoroanus"  ,     "C._sp._11"       ,    "C._lastellianus",    
"C._nauseosus"    ,    "C._leptocheilos" ,    "C._sp._13"         ,   "C._canaliculatus",   
"C._bejofo"        ,   "C._titan"  )
                       
new_Vonitra <- c("V._vonitrandambo","V._pusilla","V._antanambensis","V._crinita","V._utilis","V._fibrosa", "V._nossibensis","V._perrieri", "V._dransfieldii")
                       
new_Lemurophoenix <- c( "L._halleuxii " , "L._laevis")
                       
new_Masoala <- c( "M._kona" , "M._madagascariensis")
Outgroup <- c("I._thalangensis","I._polymorpha","I._mirabilis","I._namsabiensis","R._elegans","R._augusta","D._album","L._rupicola")

#calculating tr2
vonitra_clade_node=getMRCA(tr2reroot, new_Vonitra)
vonitra <- Descendants(tr2reroot, node=vonitra_clade_node, "all")
vonitra 
vonitra_tips <- Descendants(tr2reroot, node=vonitra_clade_node, "tips")
vonitra_tips
vonitra_nodes <- c(vonitra[!(vonitra %in% vonitra_tips[[1]])], vonitra_clade_node)
vonitra_nodes
PP_vonitra_2 <- sum(as.numeric(tr2reroot$node.label[vonitra_nodes-Ntip(tr2reroot)]))/length(vonitra_nodes)
PP_vonitra_2

chrys_clade_node=getMRCA(tr2reroot, new_Chrysadliocarpus)
chrys <- Descendants(tr2reroot, node=chrys_clade_node, "all")
chrys 
chrys_tips <- Descendants(tr2reroot, node=chrys_clade_node, "tips")
chrys_tips
chrys_nodes <- c(chrys[!(chrys %in% chrys_tips[[1]])], chrys_clade_node)
chrys_nodes
PP_chrys_2 <- sum(as.numeric(tr2reroot$node.label[chrys_nodes-Ntip(tr2reroot)]))/length(chrys_nodes)
PP_chrys_2

dypsis_clade_node=getMRCA(tr2reroot, new_dypsis)
dypsis <- Descendants(tr2reroot, node=dypsis_clade_node, "all")
dypsis 
dypsis_tips <- Descendants(tr2reroot, node=dypsis_clade_node, "tips")
dypsis_tips
dypsis_nodes <- c(dypsis[!(dypsis %in% dypsis_tips[[1]])], dypsis_clade_node)
dypsis_nodes
PP_dypsis_2 <- sum(as.numeric(tr2reroot$node.label[dypsis_nodes-Ntip(tr2reroot)]),na.rm = TRUE)/length(dypsis_nodes)
PP_dypsis_2

#calculating tr3
vonitra_clade_node=getMRCA(tr3reroot, new_Vonitra)
vonitra <- Descendants(tr3reroot, node=vonitra_clade_node, "all")
vonitra 
vonitra_tips <- Descendants(tr3reroot, node=vonitra_clade_node, "tips")
vonitra_tips
vonitra_nodes <- c(vonitra[!(vonitra %in% vonitra_tips[[1]])], vonitra_clade_node)
vonitra_nodes
LPP_vonitra_3 <- sum(as.numeric(tr3reroot$node.label[vonitra_nodes-Ntip(tr3reroot)]))/length(vonitra_nodes)
LPP_vonitra_3

chrys_clade_node=getMRCA(tr3reroot, new_Chrysadliocarpus)
chrys <- Descendants(tr3reroot, node=chrys_clade_node, "all")
chrys 
chrys_tips <- Descendants(tr3reroot, node=chrys_clade_node, "tips")
chrys_tips
chrys_nodes <- c(chrys[!(chrys %in% chrys_tips[[1]])], chrys_clade_node)
chrys_nodes
LPP_chrys_3 <- sum(as.numeric(tr3reroot$node.label[chrys_nodes-Ntip(tr3reroot)]))/length(chrys_nodes)
LPP_chrys_3

dypsis_clade_node=getMRCA(tr3reroot, new_dypsis)
dypsis <- Descendants(tr3reroot, node=dypsis_clade_node, "all")
dypsis 
dypsis_tips <- Descendants(tr3reroot, node=dypsis_clade_node, "tips")
dypsis_tips
dypsis_nodes <- c(dypsis[!(dypsis %in% dypsis_tips[[1]])], dypsis_clade_node)
dypsis_nodes
LPP_dypsis_3 <- sum(as.numeric(tr3reroot$node.label[dypsis_nodes-Ntip(tr3reroot)]),na.rm = TRUE)/length(dypsis_nodes)
LPP_dypsis_3

#Calculating tr4
vonitra_clade_node= getMRCA(tr4reroot, new_Vonitra)
vonitra <- Descendants(tr4reroot, node=vonitra_clade_node, "all")
vonitra 
vonitra_tips <- Descendants(tr4reroot, node=vonitra_clade_node, "tips")
vonitra_tips
vonitra_nodes <- c(vonitra[!(vonitra %in% vonitra_tips[[1]])], vonitra_clade_node)
vonitra_nodes
LPP_vonitra_4 <- sum(as.numeric(tr4reroot$node.label[vonitra_nodes-Ntip(tr4reroot)]))/length(vonitra_nodes)
LPP_vonitra_4

chrys_clade_node <- getMRCA(tr4reroot, new_Chrysadliocarpus)
chrys <- Descendants(tr4reroot, node=chrys_clade_node, "all")
chrys 
chrys_tips <- Descendants(tr4reroot, node=chrys_clade_node, "tips")
chrys_tips
chrys_nodes <- c(chrys[!(chrys %in% chrys_tips[[1]])], chrys_clade_node)
chrys_nodes
LPP_chrys_4 <- sum(as.numeric(tr4reroot$node.label[chrys_nodes-Ntip(tr4reroot)]))/length(chrys_nodes)
LPP_chrys_4

dypsis_clade_node <- getMRCA(tr4reroot, new_dypsis)
dypsis <- Descendants(tr4reroot, node=dypsis_clade_node, "all")
dypsis 
dypsis_tips <- Descendants(tr4reroot, node=dypsis_clade_node, "tips")
dypsis_tips
dypsis_nodes <- c(dypsis[!(dypsis %in% dypsis_tips[[1]])], dypsis_clade_node)
dypsis_nodes
LPP_dypsis_4 <- sum(as.numeric(tr4reroot$node.label[dypsis_nodes-Ntip(tr4reroot)]),na.rm = TRUE)/length(dypsis_nodes)
LPP_dypsis_4


#Calculating tr5
vonitra_clade_node= getMRCA(tr5reroot, new_Vonitra)
vonitra <- Descendants(tr5reroot, node=vonitra_clade_node, "all")
vonitra 
vonitra_tips <- Descendants(tr5reroot, node=vonitra_clade_node, "tips")
vonitra_tips
vonitra_nodes <- c(vonitra[!(vonitra %in% vonitra_tips[[1]])], vonitra_clade_node)
vonitra_nodes
LPP_vonitra_5 <- sum(as.numeric(tr5reroot$node.label[vonitra_nodes-Ntip(tr5reroot)]))/length(vonitra_nodes)
LPP_vonitra_5

chrys_clade_node <- getMRCA(tr5reroot, new_Chrysadliocarpus)
chrys <- Descendants(tr5reroot, node=chrys_clade_node, "all")
chrys 
chrys_tips <- Descendants(tr5reroot, node=chrys_clade_node, "tips")
chrys_tips
chrys_nodes <- c(chrys[!(chrys %in% chrys_tips[[1]])], chrys_clade_node)
chrys_nodes
LPP_chrys_5 <- sum(as.numeric(tr5reroot$node.label[chrys_nodes-Ntip(tr5reroot)]))/length(chrys_nodes)
LPP_chrys_5

dypsis_clade_node <- getMRCA(tr5reroot, new_dypsis)
dypsis <- Descendants(tr5reroot, node=dypsis_clade_node, "all")
dypsis 
dypsis_tips <- Descendants(tr5reroot, node=dypsis_clade_node, "tips")
dypsis_tips
dypsis_nodes <- c(dypsis[!(dypsis %in% dypsis_tips[[1]])], dypsis_clade_node)
dypsis_nodes
LPP_dypsis_5 <- sum(as.numeric(tr5reroot$node.label[dypsis_nodes-Ntip(tr5reroot)]),na.rm = TRUE)/length(dypsis_nodes)
LPP_dypsis_5

#Calculating tr10
vonitra_clade_node= getMRCA(tr10reroot, new_Vonitra)
vonitra <- Descendants(tr10reroot, node=vonitra_clade_node, "all")
vonitra 
vonitra_tips <- Descendants(tr10reroot, node=vonitra_clade_node, "tips")
vonitra_tips
vonitra_nodes <- c(vonitra[!(vonitra %in% vonitra_tips[[1]])], vonitra_clade_node)
vonitra_nodes
LPP_vonitra_10 <- sum(as.numeric(tr10reroot$node.label[vonitra_nodes-Ntip(tr10reroot)]))/length(vonitra_nodes)
LPP_vonitra_10

chrys_clade_node <- getMRCA(tr10reroot, new_Chrysadliocarpus)
chrys <- Descendants(tr10reroot, node=chrys_clade_node, "all")
chrys 
chrys_tips <- Descendants(tr10reroot, node=chrys_clade_node, "tips")
chrys_tips
chrys_nodes <- c(chrys[!(chrys %in% chrys_tips[[1]])], chrys_clade_node)
chrys_nodes
LPP_chrys_10 <- sum(as.numeric(tr10reroot$node.label[chrys_nodes-Ntip(tr10reroot)]))/length(chrys_nodes)
LPP_chrys_10

dypsis_clade_node <- getMRCA(tr10reroot, new_dypsis)
dypsis <- Descendants(tr10reroot, node=dypsis_clade_node, "all")
dypsis 
dypsis_tips <- Descendants(tr10reroot, node=dypsis_clade_node, "tips")
dypsis_tips
dypsis_nodes <- c(dypsis[!(dypsis %in% dypsis_tips[[1]])], dypsis_clade_node)
dypsis_nodes
LPP_dypsis_10 <- sum(as.numeric(tr10reroot$node.label[dypsis_nodes-Ntip(tr10reroot)]),na.rm = TRUE)/length(dypsis_nodes)
LPP_dypsis_10



org_chrys_mean <- mean(LPP_chrys_6, LPP_chrys_7, LPP_chrys_8, LPP_chrys_9)
org_dypsis_mean <- mean(LPP_dypsis_6, LPP_dypsis_7, LPP_dypsis_8, LPP_dypsis_9)
LPP_chrys_1
LPP_chrys_6
LPP_chrys_7
LPP_chrys_8
LPP_chrys_9
LPP_dypsis_1
LPP_dypsis_6
LPP_dypsis_7
LPP_dypsis_8
LPP_dypsis_9

new_chrys_mean <- mean(LPP_chrys_3,LPP_chrys_4, LPP_chrys_5, LPP_chrys_10)
new_dypsis_mean <- mean(LPP_dypsis_3, LPP_dypsis_4, LPP_dypsis_5, LPP_dypsis_10)
PP_dypsis_2
LPP_dypsis_3
LPP_dypsis_4
LPP_dypsis_5
LPP_dypsis_10
PP_chrys_2
LPP_chrys_3
LPP_chrys_4
LPP_chrys_5
LPP_chrys_10


################################################################################
######################## mean support value across tree ########################
################################################################################
#original cosensus
#large ML
n_val_tr1reroot <- as.numeric(tr1reroot$node.label)
node_val_tr1reroot <- mean(n_val_tr1reroot, na.rm = TRUE)
node_val_tr1reroot

#large consensus
n_val_tr2reroot <- as.numeric(tr2reroot$node.label)
node_val_tr2reroot <- mean(n_val_tr2reroot, na.rm = TRUE)
node_val_tr2reroot

#large ML
n_val_tr3reroot <- as.numeric(tr3reroot$node.label)
node_val_tr3reroot <- mean(n_val_tr3reroot, na.rm = TRUE)

print(node_val_tr3reroot)

#large Bayes
n_val_tr4reroot <- as.numeric(tr4reroot$node.label)
node_val_tr4reroot <- mean(n_val_tr4reroot, na.rm = TRUE)
node_val_tr4reroot

#large Ml converged
n_val_tr5reroot <- as.numeric(tr5reroot$node.label)
node_val_tr5reroot <- mean(n_val_tr5reroot, na.rm = TRUE)
node_val_tr5reroot

#original
n_val_tr6reroot <- as.numeric(tr6reroot$node.label)
node_val_tr6reroot <- mean(n_val_tr6reroot, na.rm = TRUE)
node_val_tr6reroot

#original
n_val_tr7reroot <- as.numeric(tr7reroot$node.label)
node_val_tr7reroot <- mean(n_val_tr7reroot, na.rm = TRUE)
node_val_tr7reroot

#original
n_val_tr8reroot <- as.numeric(tr8reroot$node.label)
node_val_tr8reroot <- mean(n_val_tr8reroot, na.rm = TRUE)
node_val_tr8reroot

#original
n_val_tr9reroot <- as.numeric(tr9reroot$node.label)
node_val_tr9reroot <- mean(n_val_tr9reroot, na.rm = TRUE)
node_val_tr9reroot

#original
n_val_tr10reroot <- as.numeric(tr10reroot$node.label)
node_val_tr10reroot <- mean(n_val_tr10reroot, na.rm = TRUE)
node_val_tr10reroot

node_val_tr1reroot
node_val_tr2reroot
node_val_tr3reroot
node_val_tr4reroot
node_val_tr5reroot
node_val_tr6reroot
node_val_tr7reroot
node_val_tr8reroot
node_val_tr9reroot
node_val_tr10reroot

################################################################################
###################### RF distances#############################################
################################################################################

all_trees <- c(tr4reroot,tr5reroot,tr10reroot,tr3reroot, tr7reroot,tr8reroot,tr6reroot, tr9reroot)
all_original <- c(tr6reroot, tr7reroot, tr8reroot, tr9reroot)
all_new <- c(tr3reroot, tr4reroot, tr5reroot, tr10reroot)
treefs_RF <- multiRF(all_trees, multi2di = FALSE )

tre_RF_new <- multiRF(all_new, multi2di = FALSE )
tre_RF_original <- multiRF(all_original, multi2di = FALSE )






################################################################################
#####################Make co-phylo plot ########################################
################################################################################
# All the Co-phylo plot codes are edited versions of a co-phyloplot script from Paola De Lima Ferreira


##MAKING CO-PHYLO Large and small dataset soncensus trees
obj<-cophylo(tr3reroot,tr3_info,space= 7, use.edge.length	= FALSE)
pdf("tr3_vs_tr3_info.pdf", height=70, width=45)
plot(obj,link.type="curved",link.lwd=0.8, link.lty="solid",cex=0.1, link.col=make.transparent("black", 0.8), fsize =2, label.offset = 0.5, use.edge.lenth = FALSE)
nodelabels.cophylo(obj$trees[[1]]$node.label,1:obj$trees[[1]]$Nnode+
                     Ntip(obj$trees[[1]]), cex=1.2)
title("tr3", adj=0.1, line=-5, cex.main=4)
nodelabels.cophylo(obj$trees[[2]]$node.label,1:obj$trees[[2]]$Nnode+
                     Ntip(obj$trees[[2]]),which="right", cex=1.2)
title("tr3_info", adj=0.85, line=-5, cex.main=4)



dev.off()




