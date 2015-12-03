library(gplots)
library(RColorBrewer)

x=as.matrix(read.table("fungi-projects/20120606/kegg_pathway_projection_0.02", header=T, sep="\t"))

pal = brewer.pal(10, "Purples")
fungi = read.table("fungi-data/tree/FungalSpeciesClass.clf", header=T)
o = read.table("fungi-data/tree/order")
col = as.character(fungi[as.vector(o[,2]), "Color"])

draw <- function(group, x){
  png(paste("fungi-results/kegg-pathways/kegg_pathway_projection_002_", group, ".png", sep = ""), width=2048, height=2048, pointsize=28)
  par(oma=c(0,0,0,14))
  z = x[x[,1] == group, 2:dim(x)[2]]
  class(z) = "numeric"
  heatmap.2(z[,o[,2]], scale="none", Rowv=NA, Colv=NA, dendrogram="none", trace="none", col=pal, key=TRUE, keysize=0.9, density.info=c("histogram"), densadj=0.25, ColSideColors=col)
  dev.off()
}

draw("Carbohydrate", x)
draw("Amino Acid", x)
draw("Carbohydrate", x)
draw("Cofactors and Vitamins", x)
draw("Energy", x)
draw("Glycan", x)
draw("Lipid", x)
draw("Metabolic pathways", x)
draw("Nucleotide", x)
draw("Other Amino Acids", x)
draw("Other Secondary Metabolites", x)
draw("Terpenoids and Polyketides", x)
draw("Translation", x)
draw("Xenobiotics Biodegradation", x)
