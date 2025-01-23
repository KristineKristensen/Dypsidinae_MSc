
library(phytools)
setwd("C:/Users/kristine/documents/")

#Import Data - trees

#IQ-TREE wAstral tree
tr1 <- read.tree("renamed_astral_ML_converged.tre")
plot(tr1, show.node.label = TRUE, use.edge.length = F)
title("MrBayes")
nodelabels()

tr1reroot <- reroot(tr1, node.number=326, position=NULL, interactive=FALSE)
plot(tr1reroot, show.node.label = TRUE, use.edge.length = F, cex=0.4)


#consensus tree
tr2 <- readNexus("MrBayes_consensus_renamed_extra.tre")
plot(tr2, show.node.label = TRUE, cex=0.5, label.offset = 1.2)
nodelabels()

tr2reroot <- reroot(tr2, node.number=265,interactive=FALSE)
plot(tr2reroot, show.node.label = TRUE, use.edge.length = F, cex=0.4)
title("consensus")
nodelabels()


#Check the tip names
tr1reroot$tip.label <-gsub("_"," ",tr1reroot$tip.label)
tr2reroot$tip.label <-gsub("_"," ",tr2reroot$tip.label)
match(tr1reroot$tip.label, tr2reroot$tip.label)

tr1_labels <- tr1reroot$tip.label
tr2_labels <- tr2reroot$tip.label
missing_in_tr1 <- setdiff(tr2_labels, tr1_labels)
missing_in_tr1
missing_in_tr2 <- setdiff(tr1_labels, tr2_labels)
missing_in_tr2  

#check if nodes are the same
biparts <- prop.clades(tr1reroot, tr2reroot)
biparts 

M<-matchNodes(tr1reroot,tr2reroot)
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


#Making the figure
pdf("ML_with_PP.pdf", height=70, width=35)
plot(tr1reroot, use.edge.length = F)
new_nodelabels <- nodelabels(text=paste(tr1reroot$node.label, "/PP", M[,3], sep=" "), cex=0.5, horiz =0.6)

dev.off()

