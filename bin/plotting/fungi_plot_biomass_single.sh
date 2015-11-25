#!/bin/bash
#
# This is called from fungi_plot_biomass.sh to plot biomass component yield heatmaps
#
# - Plots only nonancestral species
# - Orders species by their phylogenetic order
# - Plots both scaled original and binary data


if [[ -z "$2" ]]; then
    echo "usage: $0 <directory with biomass component yields> <target file>"
    exit 2;
fi

DDIR=$1
TFN=$2

python fluxes_to_matrix.py $DDIR $TFN

# Use only nonancestral species data, remove KEGG
grep -v "N[[:digit:]]\+" $TFN | grep -v "KEGG" >$TFN.noanc

echo "
library(RColorBrewer)
library(gplots)
#pal=brewer.pal(10, \"RdBu\")
pal=brewer.pal(10, \"Spectral\")
#binpal=c(\"#FE0000\", \"#FEFEFE\")
# white, greenish
binpal=c(\"#FEFEFE\", \"#339900\")
#scalepal=c(\"#FE0000\", \"#00FE00\", \"#0000FE\", \"#FEFEFE\")
#scalepal=c(\"#FE0000\", \"#DD6000\", \"#0060FE\", \"#FEFEFE\")
#scalepal=brewer.pal(5,\"Purples\")

#scalepal <- c('red',brewer.pal(3,'Purples'))
#scalepal <- c('white',brewer.pal(3,'Greens'))

scalepal = c('#AAAAAA', 'white', '#D5E5C0', '#A1D99B', '#31A354')
scalepal = c('white', '#AAAAAA', '#D5E5C0', '#A1D99B', '#31A354')

fungi = read.table(\"fungi-data/tree/FungalSpeciesClass.clf\", header=T)

o = read.table(\"fungi-data/tree/order\")
x = as.matrix(read.table(\"$TFN.noanc\", comment.char=\"%\"), type=\"numeric\")

# scale values of y columnwise between 0 and 1
y = t(x)
#yy = t((y - apply(y, 1, min)) / (apply(y, 1, max) - apply(y, 1, min)))
yy = t(y / apply(y, 1, max))
#yy = x

yy[yy >= 0.9] = 4
yy[yy > 0.1 & yy < 0.9] = 3
yy[yy > 0 & yy <= 0.1] = 2
yy[yy == \"NaN\"] = 1

#yy[yy >= 1.0] = 5
#yy[yy >= 0.75 & yy < 1.0] = 4
#yy[yy >= 0.5 & yy < 0.75] = 3
#yy[yy >= 0.25 & yy < 0.5] = 2
#yy[yy > 0 & yy < 0.25] = 1
#yy[yy == 0] = 0 
#yy[yy == \"NaN\"] = 3

col = as.character(fungi[as.vector(o[,2]), \"Color\"])

pdf(\"$TFN.pdf\", pointsize=10)
par(omi=c(1.3,0,0,1))
heatmap.2(x[o[,2],], Rowv=NA, Colv=NA, dendrogram=\"none\", scale=\"col\", col=pal, trace=\"none\",
   labRow=fungi[as.character(o[,2]),4],
   key=TRUE,
   keysize = 0.9,
   density.info=c(\"histogram\"),
   densadj = 0.25,
   RowSideColors=col)
dev.off()

png(\"$TFN.png\", pointsize=36, width=2048, height=2048)
par(omi=c(1.3,0,0,1))
heatmap.2(x[o[,2],], Rowv=NA, Colv=NA, dendrogram=\"none\", scale=\"col\", col=pal, trace=\"none\",
   labRow=fungi[as.character(o[,2]),4],
   key=TRUE,
   keysize = 0.9,
   density.info=c(\"histogram\"),
   densadj = 0.25,
   RowSideColors=col)
dev.off()

png(\"$TFN-dendro.png\", pointsize=36, width=2048, height=2048)
par(omi=c(1.3,0,0,1))
heatmap.2(x[o[,2],], Rowv=TRUE, Colv=TRUE, scale=\"col\", col=pal, trace=\"none\",
   labRow=fungi[as.character(o[,2]),4],
   key=TRUE,
   keysize = 0.9,
   density.info=c(\"histogram\"),
   densadj = 0.25,
   RowSideColors=col)
dev.off()

pdf(\"$TFN-bin.pdf\", pointsize=10)
par(omi=c(1.3,0,0,1))
x[x > 0] = 1
heatmap.2(x[o[,2],], Rowv=NA, Colv=NA, dendrogram=\"none\", scale=\"none\", col=binpal, trace=\"none\",
   labRow=fungi[as.character(o[,2]),4],
    key=FALSE,
#   key=TRUE,
#   keysize = 0.9,
#   density.info=c(\"histogram\"),
#   densadj = 0.25,
   RowSideColors=col)
dev.off()

png(\"$TFN-bin.png\", pointsize=36, width=2048, height=2048)
par(omi=c(1.3,0,0,1))
heatmap.2(x[o[,2],], Rowv=NA, Colv=NA, dendrogram=\"none\", scale=\"none\", col=binpal, trace=\"none\",
   labRow=fungi[as.character(o[,2]),4],
    key=FALSE,
#   key=TRUE,
#   keysize = 0.9,
#   density.info=c(\"histogram\"),
#   densadj = 0.25,
   RowSideColors=col)
dev.off()

png(\"$TFN-bin-dendro.png\", pointsize=36, width=2048, height=2048)
par(omi=c(1.3,0,0,1))
heatmap.2(x[o[,2],], Rowv=TRUE, Colv=TRUE, scale=\"none\", col=binpal, trace=\"none\",
   labRow=fungi[as.character(o[,2]),4],
    key=FALSE,
#   key=TRUE,
#   keysize = 0.9,
#   density.info=c(\"histogram\"),
#   densadj = 0.25,
   RowSideColors=col)
dev.off()

pdf(\"$TFN-scale.pdf\", pointsize=10)
par(omi=c(1.3,0,0,1))
heatmap.2(yy[o[,2],], Rowv=NA, Colv=NA, dendrogram=\"none\", scale=\"none\", col=scalepal, trace=\"none\",
   labRow=fungi[as.character(o[,2]),4],
    key=FALSE,
#   key=TRUE,
#   keysize = 0.9,
#   density.info=c(\"histogram\"),
#   densadj = 0.25,
   RowSideColors=col)
dev.off()

png(\"$TFN-scale.png\", pointsize=36, width=2048, height=2048)
par(omi=c(1.3,0,0,1))
heatmap.2(yy[o[,2],], Rowv=NA, Colv=NA, dendrogram=\"none\", scale=\"none\", col=scalepal, trace=\"none\",
   labRow=fungi[as.character(o[,2]),4],
    key=FALSE,
#   key=TRUE,
#   keysize = 0.9,
#   density.info=c(\"histogram\"),
#   densadj = 0.25,
   RowSideColors=col)
dev.off()

png(\"$TFN-scale-dendro.png\", pointsize=36, width=2048, height=2048)
par(omi=c(1.3,0,0,1))
heatmap.2(yy[o[,2],], Rowv=TRUE, Colv=TRUE, scale=\"none\", col=scalepal, trace=\"none\", 
   labRow=fungi[as.character(o[,2]),4],
    key=FALSE,
#   key=TRUE,
#   keysize = 0.9,
#   density.info=c(\"histogram\"),
#   densadj = 0.25,
   RowSideColors=col)
dev.off()

scalepal
" >fluxrtmp.R

R --no-save <fluxrtmp.R

pdfcrop $TFN.pdf
pdfcrop $TFN-bin.pdf
pdfcrop $TFN-scale.pdf

rm fluxrtmp.R
