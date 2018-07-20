d = 17
k = 30:100
k3 = rep(30:100, each=3)

L = list()
M = list()

M17_2 = as.matrix(read.table("../False_Positives/M17.txt"))
z3_2 = M17_2[151:200,c(1,3,5)]
C3_2 = M17_2[151:200,c(2,4,6)]

M17_3 = as.matrix(read.table("M17.txt"))
z3_3 = M17_3[151:200,c(1,3,5)]
C3_3 = M17_3[151:200,c(2,4,6)]

for (i in 1:50) {
   # All the terms
   #L[[i]] = Re(tapply(X = C3_3[i,]/z3_3[i,]^(k3+1) - C3_2[i,]/z3_2[i,]^(k3+1), INDEX=k3, sum))
   # P2 simplified
   #L[[i]] = Re(tapply(X = C3_3[i,]/z3_3[i,]^(k3+1), INDEX=k3, sum)) - Re(C3_2[i,1]/z3_2[i,1]^(k+1))
   # P3 simplified
   #L[[i]] = Re(C3_3[i,1]/z3_3[i,1]^(k+1)) - Re(tapply(X = C3_2[i,]/z3_2[i,]^(k3+1), INDEX=k3, sum)) 
   # Both simplified
   L[[i]] = Re(C3_3[i,1]/z3_3[i,1]^(k+1)) - Re(C3_2[i,1]/z3_2[i,1]^(k+1))
}

pdf("sht.pdf")
ymin = 0
ymax = max(unlist(L))
plot(L[[1]], type='l', ylim=c(ymin,ymax))
for (i in 33:50) {
   lines(L[[i]])
}
dev.off()
