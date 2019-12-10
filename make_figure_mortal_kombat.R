mem = log10(as.matrix(read.table("results_mem.txt")))
skip = log10(as.matrix(read.table("results_skip.txt")))

k = 20:160
ylim = c(-10,0)

pdf("mortal_kombat.pdf", height=6, width=12, useDingbats=FALSE)

par(mfrow=c(1,2))
par(mar=c(3.5,3.5,1.5,1.5))

plot(k, mem[1,], type="n",
   panel.first=grid(), ylim=ylim,
   bty="n", xaxt="n", yaxt="n", xlab="", ylab="")
# Plot first so that it is below.
legend(x="topright", inset=.02, bg="white", box.col="white",
   legend="MEM seeds", text.col="grey25")
lines(k, mem[1,], col="grey50", lwd=1.05)
lines(k, mem[2,], col="grey45", lwd=1.10)
lines(k, mem[3,], col="grey40", lwd=1.15)
lines(k, mem[4,], col="grey35", lwd=1.20)
lines(k, mem[5,], col="grey30", lwd=1.25)
lines(k, mem[6,], col="grey25", lwd=1.30)
lines(k, mem[7,], col="grey20", lwd=1.35)
lines(k, mem[8,], col="grey15", lwd=1.40)
lines(k, mem[9,], col="grey10", lwd=1.45)
lines(k, mem[10,], col="grey5", lwd=1.50)

title(xlab="Read size (nt)",
   ylab="Probability that seeding is off target (log10)",
   col.lab="grey25", line=2.2, cex.lab=1.2)
axis(side=1, cex.axis=.8, col="grey50", col.axis="grey25")
axis(side=2, cex.axis=.8, col="grey50", col.axis="grey25")

legend(x="bottomleft", lwd=c(0, 1.05, 1.10, 1.15, 1.20, 1.25, 1.30,
   1.35, 1.40, 1.45, 1.50), col=c("white", "grey50", "grey45", "grey40",
   "grey35", "grey30", "grey25", "grey20", "grey15", "grey10", "grey5"),
   legend=c("", "1","5","10","20","50","100","200","400", "800","1600"),
   inset=.02, bg="white", box.col="white", text.col="grey25")
text(x=20, y=-5.4, "Duplicates", col="grey25", adj=c(0,0))

plot(k, skip[1,], type="n",
   panel.first=grid(), ylim=ylim,
   bty="n", xaxt="n", yaxt="n", xlab="", ylab="")
# Plot first so that it is below.
legend(x="topright", inset=.02, bg="white", box.col="white",
   legend="Skip seeds", text.col="grey25")
lines(k, skip[1,], col="grey50", lwd=1.05)
lines(k, skip[2,], col="grey45", lwd=1.10)
lines(k, skip[3,], col="grey40", lwd=1.15)
lines(k, skip[4,], col="grey35", lwd=1.20)
lines(k, skip[5,], col="grey30", lwd=1.25)
lines(k, skip[6,], col="grey25", lwd=1.30)
lines(k, skip[7,], col="grey20", lwd=1.35)
lines(k, skip[8,], col="grey15", lwd=1.40)
lines(k, skip[9,], col="grey10", lwd=1.45)
lines(k, skip[10,], col="grey5", lwd=1.50)

title(xlab="Read size (nt)",
   ylab="Probablity that seeding is off target (log10)",
   col.lab="grey25", line=2.2, cex.lab=1.2)
axis(side=1, cex.axis=.8, col="grey50", col.axis="grey25")
axis(side=2, cex.axis=.8, col="grey50", col.axis="grey25")

dev.off()
