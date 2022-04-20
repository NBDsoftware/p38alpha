library(bio3d)
setwd("/Users/mesguerra/projects/mapk14/7pv3/cleanLIGCNA")
pdbfile <- "../repeat1/gppmin_ref2.pdb"
pdb <- read.pdb(pdbfile)

# 1 frame is 0.2 ns
for (val in 1: 5) {
  assign(paste("dcd",val , sep=""), read.dcd(sprintf("../repeat%s/lig%s.dcd",val ,val)))
}


fit.inds <- atom.select(pdb, elety="CA")
xyz1 <- fit.xyz(fixed=pdb$xyz, 
               mobile=dcd1,
               fixed.inds=fit.inds$xyz,
               mobile.inds=fit.inds$xyz)
dim(xyz1) == dim(dcd1)

xyz2 <- fit.xyz(fixed=pdb$xyz, 
                mobile=dcd2,
                fixed.inds=fit.inds$xyz,
                mobile.inds=fit.inds$xyz)

xyz3 <- fit.xyz(fixed=pdb$xyz, 
                mobile=dcd3,
                fixed.inds=fit.inds$xyz,
                mobile.inds=fit.inds$xyz)

xyz4 <- fit.xyz(fixed=pdb$xyz, 
                mobile=dcd4,
                fixed.inds=fit.inds$xyz,
                mobile.inds=fit.inds$xyz)

xyz5 <- fit.xyz(fixed=pdb$xyz, 
                mobile=dcd5,
                fixed.inds=fit.inds$xyz,
                mobile.inds=fit.inds$xyz)
system("rm -rf ncdfs")
system("mkdir -p ncdfs")
write.ncdf(xyz1, trjfile="ncdfs/lig1.ncdf", cell = NULL)
write.ncdf(xyz2, trjfile="ncdfs/lig2.ncdf", cell = NULL)
write.ncdf(xyz3, trjfile="ncdfs/lig3.ncdf", cell = NULL)
write.ncdf(xyz4, trjfile="ncdfs/lig4.ncdf", cell = NULL)
write.ncdf(xyz5, trjfile="ncdfs/lig5.ncdf", cell = NULL)

files <- list.files('ncdfs', full.names=TRUE)

chunks <- lapply(files, read.ncdf, first=51, time=FALSE, at.sel=atom.select(pdb, elety="CA"))
traj <- do.call(rbind, chunks)
cijs <- lapply(chunks, dccm)
cm <- cmap.xyz(traj, dcut=10, scut=0, pcut=0.75)
# calculate consensus cij
cij <- filter.dccm(cijs, cutoff.cij=0.60, cmap=cm)

# View the correlations in pymol.
pymol(cij, pdb, type = "script")

# build the network (Note the cutoff of cij here is 0 because we have already removed weak couplings from above calling of the **filter.dccm()**
net <- cna(cij, cutoff.cij=0, cluster.method="btwn")
pruned <- prune.cna(net, size.min=10)


plot.dccm(cij, margin.segments = pruned$communities$membership,
          at=c(0,0.40, 0.41, 0.42, 0.44, 0.46, 0.50, 0.54, 0.56, 0.60),
          main="Dynamical cross-correlation map (DCCM) ",
          sub="Pair-wise residue CC using the MD simulation of 7pv3 with compound 37",
          sse=dssp(pdb, exefile='/opt/local/bin/mkdssp'),
          helix.col = "red",
          sheet.col = "blue",
          contour = FALSE,
          col.regions = hcl.colors(16, "RdPu", rev = TRUE))

# Check modularity of network
tree <- community.tree(net, rescale=TRUE)
plot( tree$num.of.comms, tree$modularity, xlab="Communities", ylab="Modularity", typ="o", pch=20, panel.first = grid() ) 
plot( tree$num.of.comms, tree$modularity, xlab="Communities", ylab="Modularity", typ="o", pch=20, panel.first = grid(), xlim=c(0,50) ) 
max.mod.ind <- which.max(tree$modularity)
tree$num.of.comms[ max.mod.ind ]
all( net$communities$membership == tree$tree[max.mod.ind,] )
# Inspect the clustering dendrogram
h <- as.hclust(net$communities)
numclus=9
hclustplot(h, k=numclus)
memb.new <- tree$tree[ tree$num.of.comms == numclus, ]
net.new <- network.amendment(net, memb.new)
prune.new <- prune.cna(net.new, size.min=15)
cent.pdbxy = layout.cna(prune.new, pdb, k=3)[,1:2]
cent.pdbxy2 = layout.cna(prune.new, pdb, full=TRUE, k=3)[,1:2]
plot(prune.new, pdb, layout=cent.pdbxy,full=F, vertex.label=NA, interactive=F,
     col=c("blue","red","gray30","yellow","chartreuse","white"))
plot.cna(prune.new, layout=cent.pdbxy2, full=TRUE, vertex.label=NA, vertex.size=5, weights=0.2)
edges <- E(prune.new$community.network)
vertices <- V(prune.new$community.network)
plot(prune.new, pdb, layout=cent.pdbxy,full=F, , interactive=F,
     col=c("blue","red","gray30","yellow","chartreuse","white"), weights = (edges$weight*20))



cent.full = layout.cna(pruned, pdb=pdb, full=TRUE, k=3)[,1:2]
cent.net = layout.cna(net, pdb=pdb, k=3)[,1:2]
cent = layout.cna(pruned, pdb=pdb, k=3)[,1:2]
layout3D <- layout.cna(net, pdb, k=3)[,c(2,3)]
#plot.cna(pruned, interactive=F) 
plot.cna(pruned, layout=cent.full, full=TRUE, vertex.label=NA, vertex.size=5, weights=0.2)
plot.cna(pruned, layout=cent, interactive=F)
plot.cna(pruned, pdb, layout=igraph::layout.reingold.tilford(pruned$community.network),
          interactive=F, vertex.size=30) 

plot(net, pdb, interactive=F, col=c("blue","red","gray30","yellow","chartreuse","white",
                                       "khaki4","orange","pink","purple","violet"), layout=cent)

node.betweenness <- betweenness(net$network)
a <- node.betweenness
diffmat <- t(outer(a, a, '-'))
tobinary <- as.matrix((diffmat<0)+0)
sumall <- rowSums(tobinary)
normsum <- sumall/length(a)
gt985x <- which(normsum >= 0.985)+(unique(pdb$atom$resno)[1]-1)
gt985y <- normsum[normsum >= 0.985]
gt985x
gt985y
pdf("centrality.pdf")
plot(unique(pdb$atom$resno),normsum*100, xlab="Residue", ylab="Centrality", 
     panel.first = grid(), type="h", col="blue",  las=1, ylim=c(0,105),
     xaxs = "i", yaxs = "i")
abline(h=98.5, col="orange", lwd=2, lty=2)
mtext("Betweenness Centrality Index for empty 7PVU WT", side=3, line=0.8, adj=0.5, cex=1)
text(c(50,100,150,200,250), 102, paste(gt985x), cex=0.8)
dev.off()

# Output a PDB file with normalized betweenness mapped to B-factor
ca.pdb <- trim.pdb(pdb, "calpha")
write.pdb(ca.pdb, b=normalize.vector(node.betweenness), file="BCI.pdb")

vmd(prune.new, pdb, launch=F, 
    radius = "4",
    alpha=0.5,
    exefile="/Applications/VMD1.9.4.app/Contents/MacOS/startup.command")


