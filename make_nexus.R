library(ape)
setwd("C:/Users/kristine/documents")
tree <- read.tree("astral_for_mrbayes.tre")
write.nexus(tree, file = "mrbayes.run1.t", translate = TRUE)
