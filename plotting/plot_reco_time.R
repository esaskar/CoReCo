

r = as.matrix(read.table("reco_times.species", row.names=1))

v = read.csv("visuals/FungalSpeciesClass_20110516.csv")

rownames(r) = v[,1]
col = as.character(v[,8])

pdf("plots/reco_times.pdf")

plot(r[,2], r[,1], xlab="Number of reactions", ylab="Reconstruction time (seconds)", bty="n", col=col, bg=col, pch=21, cex=2.5, xlim=c(2600,4300), ylim=c(23000,48000))

text(r[,2], r[,1], rownames(r), cex=0.7)

dev.off()