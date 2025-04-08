#packages
library(phytools)
library(phangorn)
library(ape)


setwd("C:/Users/kristine/documents/")
#Import Data - trees
#################################################################################
#Heyduk consensus tree
tr1 <- read.tree("Heyduk_consensus.tre")
trimmed_tree <- drop.tip(tr1, c("0168", "New_sequence"))
plot(trimmed_tree, show.node.label = TRUE, use.edge.length = T)
title("Heyduk")
nodelabels()

root_node <- getMRCA(trimmed_tree, c("Loxococcus_rupicola_1", "Loxococcus_rupicola_2"))
tr1reroot <- reroot(trimmed_tree, node.number=root_node, interactive=FALSE)
plot(tr1reroot, show.node.label = TRUE, use.edge.length = T, cex=0.5)
is.rooted(tr1reroot)
nodelabels()

pdf("Heyduk_consensus.pdf", height=35, width=18)
plot.phylo(tr1reroot,use.edge.length = TRUE, cex=1, label.offset=0.0005)
nodelabels(tr1reroot$node.label,adj = c(1,- 0.2), bg="transparent", frame="none")
title("Heyduk")
dev.off()

##################################################################################
#Phylopalm consensus tree
tr2 <- read.tree("Phylopalm_consensus.tre")
plot(tr2, show.node.label = TRUE, cex=0.5)
nodelabels()

tr2reroot <- reroot(tr2, node.number=265,interactive=FALSE)
plot(tr2reroot, show.node.label = TRUE, use.edge.length = T, cex=0.4)
title("consensus")
nodelabels()
is.rooted(tr2reroot)

pdf("consensus_phylopalm_2.pdf", height=35, width=20)
plot.phylo(tr2reroot,use.edge.length = TRUE, label.offset=0.0005)
nodelabels(tr2reroot$node.label,adj = c(1,- 0.2), bg="transparent", frame="none")
title("Phylopalm majority rule consensus")
dev.off()

################################################################################
#IQ-TREE wAstral tree
tr3 <- read.tree("renamed_astral_ML_jan25.tre")
plot(tr3, show.node.label = TRUE, use.edge.length = T)
title("MrBayes")
nodelabels()

#rerooting the tree
root_node <- getMRCA(tr3, c("Dictyosperma_album", "Rhopaloblaste_elegans", "Rhopaloblaste_augusta", "Iguanura_mirabilis", "Iguanura_polymorpha", "Iguanura_thalangensis", "Iguanura_namsabiensis" ,"Loxococcus_rupicola"))
tr3reroot <- reroot(tr3, node.number=root_node, position=NULL, interactive=FALSE)
plot(tr3reroot, show.node.label = TRUE, use.edge.length = T, cex=0.4)


#moving decimals from PP values
for (i in seq_along(tr3reroot$node.label)) {
  if(tr3reroot$node.label[i] == "Root" ) {next
  } else {
    tr3reroot$node.label[i] <- format(round(as.numeric(tr3reroot$node.label[i]), 3), nsmall = 3)
  }
}


tr3reroot$edge.length[is.na(tr3reroot$edge.length)] <- 0.2



################################################################################
#MrBayes wAstral tree
tr4 <- read.tree("renamed_astral_Bayes_jan25.tre")
plot(tr4, show.node.label = TRUE, use.edge.length = F)
title("MrBayes")
nodelabels()

tr4reroot <- reroot(tr4, node.number=278, interactive=FALSE)
tr4reroot$edge.length[is.na(tr4reroot$edge.length)] <- 0.02
plot(tr4reroot, show.node.label = TRUE, use.edge.length = T, cex=0.4)

for (i in seq_along(tr4reroot$node.label)) {
  if(tr4reroot$node.label[i] == "Root" ) {next
  } else {
    tr4reroot$node.label[i] <- format(round(as.numeric(tr4reroot$node.label[i]), 3), nsmall = 3)
  }
}

################################################################################
#IQ-TREE wAstral tree converged genes
tr5 <- read.tree("renamed_astral_ML_converged_feb25.tre")
plot(tr5, show.node.label = TRUE, use.edge.length = F)
title("converged_IQtree")
nodelabels()

#rerooting the tree
root_node <- getMRCA(tr5, c("Dictyosperma_album", "Rhopaloblaste_elegans", "Rhopaloblaste_augusta", "Iguanura_mirabilis", "Iguanura_polymorpha", "Iguanura_thalangensis", "Iguanura_namsabiensis" ,"Loxococcus_rupicola"))
tr5reroot <- reroot(tr5, node.number=root_node, position=NULL, interactive=FALSE)
plot(tr5reroot, show.node.label = TRUE, use.edge.length = T, cex=0.4)


#moving decimals from PP values
for (i in seq_along(tr5reroot$node.label)) {
  if(tr5reroot$node.label[i] == "Root" ) {next
  } else {
    tr5reroot$node.label[i] <- format(round(as.numeric(tr5reroot$node.label[i]), 3), nsmall = 3)
  }
}


tr5reroot$edge.length[is.na(tr5reroot$edge.length)] <- 0.2


################################################################################


#Heyduk IQ-tree wASTRAL
tr_6 <- read.tree("Heyduk_astral_ML_renamed.tre")
plot(tr_extra, show.node.label = TRUE, use.edge.length = F)
nodelabels()

#rerooting the tree
root_node <- getMRCA(tr_extra, c("Loxococcus_rupicola1", "Loxococcus_rupicola2"))
tr_extrareroot <- reroot(tr_extra, node.number=root_node, position=NULL, interactive=FALSE)
plot(tr_extrareroot, show.node.label = TRUE, use.edge.length = T, cex=0.4)



################################################################################
# Heyduk MrB tree
tr7<- read.tree("Heyduk_astral_MrB_renamed.tre")
trimmed_tree <- drop.tip(tr7, c("0168", "New_sequence"))
plot(trimmed_tree, show.node.label = TRUE, use.edge.length = T)
title("Heyduk_MrB")
nodelabels()


tr7$tip.label

root_node <- getMRCA(trimmed_tree, c("Loxococcus_rupicola_1", "Loxococcus_rupicola_2"))
tr7reroot <- reroot(trimmed_tree, node.number=root_node, interactive=FALSE)
is.rooted(tr7reroot)
tr7reroot$edge.length[is.na(tr7reroot$edge.length)] <- 0.2
plot(tr7reroot, show.node.label = TRUE, use.edge.length = T, cex=0.5)

nodelabels()

#moving decimals from PP values
for (i in seq_along(tr7reroot$node.label)) {
  if(tr7reroot$node.label[i] == "Root" ) {next
  } else {
    tr7reroot$node.label[i] <- format(round(as.numeric(tr7reroot$node.label[i]), 3), nsmall = 3)
  }
}


pdf("MrB_Heyduk.pdf", height=35, width=20)
plot.phylo(tr7reroot,use.edge.length = TRUE, label.offset=0.0005)
nodelabels(tr7reroot$node.label,adj = c(1,- 0.2), bg="transparent", frame="none")
title("Heyduk MrB")
dev.off()


################################################################################
#####################Make co-phylo plot ########################################
################################################################################
# All Co-phylo plot codes is an edited version of a co-phyloplot script from Paola De Lima Ferreira


##MAKING CO-PHYLO Large and small dataset soncensus trees
obj<-cophylo(tr1reroot,tr2reroot,space= 7, use.edge.length	= FALSE)
pdf("heyduk_vs_phylopalm_1.pdf", height=70, width=45)
plot(obj,link.type="curved",link.lwd=0.8, link.lty="solid",cex=0.1, link.col=make.transparent("black", 0.8), fsize =2, label.offset = 0.5, use.edge.lenth = FALSE)
nodelabels.cophylo(obj$trees[[1]]$node.label,1:obj$trees[[1]]$Nnode+
                     Ntip(obj$trees[[1]]), cex=1.2)
title("Heyduk", adj=0.1, line=-5, cex.main=4)
nodelabels.cophylo(obj$trees[[2]]$node.label,1:obj$trees[[2]]$Nnode+
                     Ntip(obj$trees[[2]]),which="right", cex=1.2)
title("PhyloPalm", adj=0.85, line=-5, cex.main=4)





################################################################################
#############################Move PP ###########################################
################################################################################
#Move PP large dataset consensus and IQtree ML 
#Check the tip names
tr3reroot$tip.label <-gsub("_"," ",tr3reroot$tip.label)
tr2reroot$tip.label <-gsub("_"," ",tr2reroot$tip.label)
match(tr3reroot$tip.label, tr2reroot$tip.label)

tr3_labels <- tr3reroot$tip.label
tr2_labels <- tr2reroot$tip.label
missing_in_tr3 <- setdiff(tr2_labels, tr3_labels)
missing_in_tr3
missing_in_tr2 <- setdiff(tr3_labels, tr2_labels)
missing_in_tr2  

prop.part(tr3reroot, tr2reroot)

biparts <- prop.clades(tr3reroot, tr2reroot)
biparts 

M<-matchNodes(tr3reroot,tr2reroot)
M

# list of nodes from tr2rerooted that match tr1rerooted
NN_tr2_NL=c()
for(i in 1:nrow(M)){
  if(!is.na(M[i,2])) NN_tr2_NL[length(NN_tr2_NL)+1]= M[i,2]
}
length(NN_tr2_NL)

#matrix where the nodes that match tr1 are combined with the correct PP value
tr2_NL=c()
for(i in 1:nrow(M)){
  if(!is.na(M[i,2])) { tr2_NL[length(tr2_NL)+1]= tr2reroot$node.label[M[i,2] - Ntip(tr2reroot)]
  }else{ tr2_NL[length(tr2_NL)+1]= NA
  }
}
tr2_NL
M <- cbind(M, tr2_NL)
M

#making regular pdf plot
pdf("ML_with_PP_from_cons_4.pdf", height=30, width=25)
plot(tr3reroot, use.edge.length = T)
nodelabels(text=paste(tr3reroot$node.label, "/", M[,3], sep=" "),  bg="transparent", frame="none", width = 0.2, adj=c(1.05, -0.5), cex=0.8) 
title("ML with PP from Majority rule consensus")
dev.off()

#Making split plot
#from liam revell phytools blog
split.plotTree<-function(tree,splits=NULL,file=NULL,fn,...){
  ef<-0.037037037037
  if(!is.null(file)) pdf(file,width=30,height=40)
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
  nodelabels(text=paste(tr3reroot$node.label, "/", M[,3], sep=" "),  bg="transparent", frame="none", width = 0.2, adj=c(1.05, -0.5)) 
  
}

splits<- 0.522
tree <- tr3reroot
tree$edge.length[is.na(tree$edge.length)] <- 0.2
split.plotTree(tree,splits,file= "split_PP_move_ML_feb_test.pdf",ftype="i",mar=rep(1.1,4),fn=foo,lwd=1)
dev.off()


################################################################################
#move PP phylopalm MrB consensus and MrB tree
#Check the tip names
tr4reroot$tip.label <-gsub("_"," ",tr4reroot$tip.label)
tr2reroot$tip.label <-gsub("_"," ",tr2reroot$tip.label)
match(tr4reroot$tip.label, tr2reroot$tip.label)

tr4_labels <- tr4reroot$tip.label
tr2_labels <- tr2reroot$tip.label
missing_in_tr4 <- setdiff(tr2_labels, tr4_labels)
missing_in_tr4
missing_in_tr2 <- setdiff(tr4_labels, tr2_labels)
missing_in_tr2  


biparts <- prop.clades(tr4reroot, tr2reroot)
biparts 

M<-matchNodes(tr4reroot,tr2reroot)
M

# list of nodes from tr2rerooted that match tr1rerooted
NN_tr2_NL=c()
for(i in 1:nrow(M)){
  if(!is.na(M[i,2])) NN_tr2_NL[length(NN_tr2_NL)+1]= M[i,2]
}
length(NN_tr2_NL)

#matrix where the nodes that match tr1 are combined with the correct PP value
tr2_NL=c()
for(i in 1:nrow(M)){
  if(!is.na(M[i,2])) { tr2_NL[length(tr2_NL)+1]= tr2reroot$node.label[M[i,2] - Ntip(tr2reroot)]
  }else{ tr2_NL[length(tr2_NL)+1]= NA
  }
}
tr2_NL
M <- cbind(M, tr2_NL)
M <- cbind(M, tr4reroot$node.label)


tr4_NL=c()
tr4_NL <- tr4reroot$node.label
M <- cbind(M, tr4_NL)
M                    

#difference between LPP and PP
Difference= c()
for (i in 1:nrow(M)){
  if (!is.na(as.numeric(M[i,3])- as.numeric(M[i,4]))){
    Difference[[length(Difference)+1]] <- c(M[i,1],M[i,2],as.numeric(M[i,3])- as.numeric(M[i,4]))
  }
}
Difference


################################################################################
# move PP from large dataset pp to large dataset converged genes IQtree ML tree
#Check the tip names
tr5reroot$tip.label <-gsub("_"," ",tr5reroot$tip.label)
tr2reroot$tip.label <-gsub("_"," ",tr2reroot$tip.label)
match(tr5reroot$tip.label, tr2reroot$tip.label)

tr5_labels <- tr5reroot$tip.label
tr2_labels <- tr2reroot$tip.label
missing_in_tr5 <- setdiff(tr2_labels, tr5_labels)
missing_in_tr5
missing_in_tr2 <- setdiff(tr5_labels, tr2_labels)
missing_in_tr2  

prop.part(tr5reroot, tr2reroot)

biparts <- prop.clades(tr5reroot, tr2reroot)
biparts 

M<-matchNodes(tr5reroot,tr2reroot)
M

# list of nodes from tr2rerooted that match tr1rerooted
NN_tr2_NL=c()
for(i in 1:nrow(M)){
  if(!is.na(M[i,2])) NN_tr2_NL[length(NN_tr2_NL)+1]= M[i,2]
}
length(NN_tr2_NL)

#matrix where the nodes that match tr1 are combined with the correct PP value
tr2_NL=c()
for(i in 1:nrow(M)){
  if(!is.na(M[i,2])) { tr2_NL[length(tr2_NL)+1]= tr2reroot$node.label[M[i,2] - Ntip(tr2reroot)]
  }else{ tr2_NL[length(tr2_NL)+1]= NA
  }
}
tr2_NL
M <- cbind(M, tr2_NL)
M

#making regular pdf plot
pdf("ML_with_PP_from_cons_4.pdf", height=30, width=25)
plot(tr3reroot, use.edge.length = T)
nodelabels(text=paste(tr3reroot$node.label, "/", M[,3], sep=" "),  bg="transparent", frame="none", width = 0.2, adj=c(1.05, -0.5), cex=0.8) 
title("ML with PP from Majority rule consensus")
dev.off()

#Making split plot
#from liam revell phytools blog
split.plotTree<-function(tree,splits=NULL,file=NULL,fn,...){
  ef<-0.037037037037
  if(!is.null(file)) pdf(file,width=30,height=40)
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
  nodelabels(text=paste(tr5reroot$node.label, "/", M[,3], sep=" "),  bg="transparent", frame="none", width = 0.2, adj=c(1.05, -0.5)) 
  
}

splits<- 0.522
tree <- tr5reroot
tree$edge.length[is.na(tree$edge.length)] <- 0.2
split.plotTree(tree,splits,file= "split_PP_move_ML_converged.pdf",ftype="i",mar=rep(1.1,4),fn=foo,lwd=1)
dev.off()



################################################################################
#move PP Heyduk MrB consensus and MrB tree
#Check the tip names
tr7reroot$tip.label <-gsub("_"," ",tr7reroot$tip.label)
tr1reroot$tip.label <-gsub("_"," ",tr1reroot$tip.label)
match(tr7reroot$tip.label, tr1reroot$tip.label)

tr7_labels <- tr7reroot$tip.label
tr1_labels <- tr1reroot$tip.label
missing_in_tr7 <- setdiff(tr1_labels, tr7_labels)
missing_in_tr7
missing_in_tr1 <- setdiff(tr7_labels, tr1_labels)
missing_in_tr1  


biparts <- prop.clades(tr7reroot, tr1reroot)
biparts 

M<-matchNodes(tr7reroot,tr1reroot)
M

# list of nodes from tr1rerooted that match tr7rerooted
NN_tr1_NL=c()
for(i in 1:nrow(M)){
  if(!is.na(M[i,2])) NN_tr1_NL[length(NN_tr1_NL)+1]= M[i,2]
}
length(NN_tr1_NL)

#matrix where the nodes that match tr1 are combined with the correct PP value
tr1_NL=c()
for(i in 1:nrow(M)){
  if(!is.na(M[i,2])) { tr1_NL[length(tr1_NL)+1]= tr1reroot$node.label[M[i,2] - Ntip(tr1reroot)]
  }else{ tr1_NL[length(tr1_NL)+1]= NA
  }
}
tr1_NL
M <- cbind(M, tr1_NL)
M <- cbind(M, tr7reroot$node.label)


tr7_NL=c()
tr7_NL <- tr7reroot$node.label
M <- cbind(M, tr7_NL)
M                    

#difference between LPP and PP
Difference= c()
for (i in 1:nrow(M)){
  if (!is.na(as.numeric(M[i,3])- as.numeric(M[i,4]))){
    Difference[[length(Difference)+1]] <- c(M[i,1],M[i,2],as.numeric(M[i,3]),as.numeric(M[i,4]), as.numeric(M[i,3])- as.numeric(M[i,4]))
  }
}
Difference

#making regular pdf plot
pdf("Heyduk_MrB_with_PP_from_cons.pdf", height=30, width=25)
plot(tr7reroot, use.edge.length = T)
nodelabels(text=paste(tr7reroot$node.label, "/", M[,3], sep=" "),  bg="transparent", frame="none", width = 0.2, adj=c(1.05, -0.5), cex=0.8) 
title("Original dataset Bayesian tree with PP from Bayesian consensus tre")
dev.off()

#Making split plot
#from liam revell phytools blog
split.plotTree<-function(tree,splits=NULL,file=NULL,fn,...){
  ef<-0.037037037037
  if(!is.null(file)) pdf(file,width=30,height=40)
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
  nodelabels(text=paste(tr7reroot$node.label,"/", M[,3], sep=" "),  bg="transparent", frame="none", width = 0.2, adj=c(1.05, -0.5)) 
  
}

splits<- 0.563
tree <- tr7reroot
tree$edge.length[is.na(tree$edge.length)] <- 0.2
split.plotTree(tree,splits,file= "split_PP_move_MrB_Heyduk_1.pdf",ftype="i",mar=rep(1.1,4),fn=foo,lwd=1)
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
tree$edge.length[is.na(tree$edge.length)] <- 0.2
split.plotTree(tree,splits,file= "split_PP_move_MrB_15.pdf",ftype="i",mar=rep(1.1,4),fn=foo,lwd=1)


dev.off()

################################################################################
#Making split plot 





################################################################################
################# Make plot with Baitkit infroamtion ###########################
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
##################### calculating support averages##############################
################################################################################
#calculating tr1
vonitra_clade_node=348
vonitra <- Descendants(tr1reroot, node=vonitra_clade_node, "all")
vonitra 
vonitra_tips <- Descendants(tr1reroot, node=vonitra_clade_node, "tips")
vonitra_tips
vonitra_nodes <- c(vonitra[!(vonitra %in% vonitra_tips[[1]])], vonitra_clade_node)
vonitra_nodes
LPP_vonitra <- sum(as.numeric(tr1reroot$node.label[vonitra_nodes-Ntip(tr1reroot)]))/length(vonitra_nodes)
LPP_vonitra

chrys_clade_node=282
chrys <- Descendants(tr1reroot, node=chrys_clade_node, "all")
chrys 
chrys_tips <- Descendants(tr1reroot, node=chrys_clade_node, "tips")
chrys_tips
chrys_nodes <- c(chrys[!(chrys %in% chrys_tips[[1]])], chrys_clade_node)
chrys_nodes
LPP_chrys <- sum(as.numeric(tr1reroot$node.label[chrys_nodes-Ntip(tr1reroot)]))/length(chrys_nodes)
LPP_chrys

dypsis_clade_node=189
dypsis <- Descendants(tr1reroot, node=dypsis_clade_node, "all")
dypsis 
dypsis_tips <- Descendants(tr1reroot, node=dypsis_clade_node, "tips")
dypsis_tips
dypsis_nodes <- c(dypsis[!(dypsis %in% dypsis_tips[[1]])], dypsis_clade_node)
dypsis_nodes
dypsis_nl <- na.omit(as.numeric(tr1reroot$node.label[dypsis_nodes-Ntip(tr1reroot)]))
LPP_dypsis <- sum(dypsis_nl)/length(dypsis_nl)
LPP_dypsis

nl_chrys <- as.numeric(tr1reroot$node.label[chrys_nodes-Ntip(tr1reroot)])

median(nl_chrys)

################################################################################
#calculating tr2
vonitra_clade_node=344
vonitra <- Descendants(tr2reroot, node=vonitra_clade_node, "all")
vonitra 
vonitra_tips <- Descendants(tr2reroot, node=vonitra_clade_node, "tips")
vonitra_tips
vonitra_nodes <- c(vonitra[!(vonitra %in% vonitra_tips[[1]])], vonitra_clade_node)
vonitra_nodes
PP_vonitra <- sum(as.numeric(tr2reroot$node.label[vonitra_nodes-Ntip(tr2reroot)]))/length(vonitra_nodes)
PP_vonitra

chrys_clade_node=278
chrys <- Descendants(tr2reroot, node=chrys_clade_node, "all")
chrys 
chrys_tips <- Descendants(tr2reroot, node=chrys_clade_node, "tips")
chrys_tips
chrys_nodes <- c(chrys[!(chrys %in% chrys_tips[[1]])], chrys_clade_node)
chrys_nodes
PP_chrys <- sum(as.numeric(tr2reroot$node.label[chrys_nodes-Ntip(tr2reroot)]))/length(chrys_nodes)
PP_chrys

dypsis_clade_node=190
dypsis <- Descendants(tr2reroot, node=dypsis_clade_node, "all")
dypsis 
dypsis_tips <- Descendants(tr2reroot, node=dypsis_clade_node, "tips")
dypsis_tips
dypsis_nodes <- c(dypsis[!(dypsis %in% dypsis_tips[[1]])], dypsis_clade_node)
dypsis_nodes
PP_dypsis <- sum(as.numeric(tr2reroot$node.label[dypsis_nodes-Ntip(tr2reroot)]))/length(dypsis_nodes)
PP_dypsis




################################################################################
# calculating support averages

vonitra_clade_node=348
vonitra <- Descendants(tr3reroot, node=vonitra_clade_node, "all")
vonitra 
vonitra_tips <- Descendants(tr3reroot, node=vonitra_clade_node, "tips")
vonitra_tips
vonitra_nodes <- c(vonitra[!(vonitra %in% vonitra_tips[[1]])], vonitra_clade_node)
vonitra_nodes
LPP_vonitra <- sum(as.numeric(tr1reroot$node.label[vonitra_nodes-Ntip(tr1reroot)]))/length(vonitra_nodes)
LPP_vonitra

chrys_clade_node=282
chrys <- Descendants(tr3reroot, node=chrys_clade_node, "all")
chrys 
chrys_tips <- Descendants(tr3reroot, node=chrys_clade_node, "tips")
chrys_tips
chrys_nodes <- c(chrys[!(chrys %in% chrys_tips[[1]])], chrys_clade_node)
chrys_nodes
LPP_chrys <- sum(as.numeric(tr3reroot$node.label[chrys_nodes-Ntip(tr3reroot)]))/length(chrys_nodes)
LPP_chrys

dypsis_clade_node=189
dypsis <- Descendants(tr3reroot, node=dypsis_clade_node, "all")
dypsis 
dypsis_tips <- Descendants(tr3reroot, node=dypsis_clade_node, "tips")
dypsis_tips
dypsis_nodes <- c(dypsis[!(dypsis %in% dypsis_tips[[1]])], dypsis_clade_node)
dypsis_nodes
dypsis_nl <- na.omit(as.numeric(tr3reroot$node.label[dypsis_nodes-Ntip(tr3reroot)]))
LPP_dypsis <- sum(dypsis_nl)/length(dypsis_nl)
LPP_dypsis


################################################################################
######################## mean support value across tree ########################
################################################################################

#phylopalm ML
n_val_tr3 <- as.numeric(tr3$node.label)
node_val_tr3 <- mean(n_val_tr3, na.rm = TRUE)

print(node_val_tr3)

#heyduk ML

n_val_tr_ex <- as.numeric(tr_extrareroot$node.label)
node_val_tr_ex <- mean(n_val_tr_ex, na.rm = TRUE)  # Computes the mean while ignoring NA values

print(node_val_tr_ex)


################################################################################
###################### RF distances#############################################
################################################################################

all_trees <- c(tr1, tr2, tr3, tr5)
treefs_RF <- multiRF(all_trees, multi2di = TRUE )














