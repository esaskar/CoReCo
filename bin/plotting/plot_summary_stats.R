library(lattice)
library(MASS)
library(ellipse)

o = read.table("fungi-data/tree/order")

x=as.matrix(read.table("fungi-projects/20120724-species/reconstruction_stats_noanc.txt", row.names=1, header=T))
x=x[o[,2], ]
y=x[!rownames(x) %in% "Ecun",]
yy=y[,c("NumReactions","Gapfills","GapfillRatio")]

fungi = read.table("fungi-data/tree/FungalSpeciesClass.clf", header=T)
colors = as.character(fungi[as.vector(o[,2]), "Color"])

smyco = rownames(fungi[fungi["Subphylum"] == "Saccharomycotina", ])
smix = rownames(x) %in% smyco

emyco = rownames(fungi[fungi["Class"] == "Eurotiomycetes", ])
emix = rownames(x) %in% emyco

ix1 = rownames(x) %in% rownames(fungi[fungi["Class"] == "Dothideomycetes", ])
ix2 = rownames(x) %in% rownames(fungi[fungi["Class"] == "Leotiomycetes", ])
ix3 = rownames(x) %in% rownames(fungi[fungi["Class"] == "Mucoromycotina", ])
ix4 = rownames(x) %in% rownames(fungi[fungi["Class"] == "Schizosaccharomycetes", ])
ix5 = rownames(x) %in% rownames(fungi[fungi["Class"] == "Agaricomycetes", ])
ix6 = rownames(x) %in% rownames(fungi[fungi["Class"] == "Sordariomycetes", ])

#"Eurotiomycetes"
#"Saccharomycetes"

# add color here

pdf("reco_summary.pdf")
splom(yy, pch=19, col=colors, xlab = "",
    upper.panel=function(x, y, ...) {
        sx = x[smix]
	sy = y[smix]
	sel = ellipse(cor(sx, sy), scale=c(sd(sx), sd(sy)), centre=c(mean(sx), mean(sy)))
	panel.lines(sel[,1], sel[,2], col='#FCCDE5', lty=2)
        ex = x[emix]
	ey = y[emix]
	eel = ellipse(cor(ex, ey), scale=c(sd(ex), sd(ey)), centre=c(mean(ex), mean(ey)))
	panel.lines(eel[,1], eel[,2], col='#BEBADA', lty=2)
	panel.splom(x, y, ...)
    },
    lower.panel=function(x, y, ...) {
#        x1 = x[ix1]; y1 = y[ix1]; el = ellipse(cor(x1, y1), scale=c(sd(x1), sd(y1)), centre=c(mean(x1), mean(y1)))
#	panel.lines(el[,1], el[,2], col='#FFFFB3', lty=2)
#        x1 = x[ix2]; y1 = y[ix2]; el = ellipse(cor(x1, y1), scale=c(sd(x1), sd(y1)), centre=c(mean(x1), mean(y1)))
#	panel.lines(el[,1], el[,2], col='#FB8072', lty=2)
#        x1 = x[ix3]; y1 = y[ix3]; el = ellipse(cor(x1, y1), scale=c(sd(x1), sd(y1)), centre=c(mean(x1), mean(y1)))
#	panel.lines(el[,1], el[,2], col='#FDB462', lty=2)
#        x1 = x[ix4]; y1 = y[ix4]; el = ellipse(cor(x1, y1), scale=c(sd(x1), sd(y1)), centre=c(mean(x1), mean(y1)))
#	panel.lines(el[,1], el[,2], col='#D9D9D9', lty=2)
        x1 = x[ix5]; y1 = y[ix5]; el = ellipse(cor(x1, y1), scale=c(sd(x1), sd(y1)), centre=c(mean(x1), mean(y1)))
	panel.lines(el[,1], el[,2], col='#8DD3C7', lty=2)
        x1 = x[ix6]; y1 = y[ix6]; el = ellipse(cor(x1, y1), scale=c(sd(x1), sd(y1)), centre=c(mean(x1), mean(y1)))
	panel.lines(el[,1], el[,2], col='#BC80BD', lty=2)
	panel.splom(x, y, ...)
    },
    diag.panel=function(x, ...) {
        yrng=current.panel.limits()$ylim;
 	d=density(x)
	d$y=with(d, yrng[1] + 0.95 * diff(yrng) * y / max(y))
	panel.lines(d)
	diag.panel.splom(x, ...)
    }
)
dev.off()