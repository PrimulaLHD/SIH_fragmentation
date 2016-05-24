#plot network examples####
require(igraph)
source("./Functions/maketriangles.r")

numCom<-100
eAMP<-1
ePeriod<-40000


shapeV<-rep("circle",numCom)

add.vertex.shape("uptriangle",plot=uptriangle)
add.vertex.shape("downtriangle",plot=downtriangle)

success<-F
while(!success){
  landscape<-round(data.frame(x = runif(numCom, min = 1, max = 1000), y = runif(numCom, min = 1, max = 1000)))
  distance_mat1<-as.matrix(dist(landscape,method = "euclidean",diag = T,upper=T))
  
  distance_mat<-1*(distance_mat1<200)
  diag(distance_mat)<-0
  connections<-distance_mat
  distance_mat[upper.tri(distance_mat)]<-0
  
  graph<-as.undirected(graph.adjacency(distance_mat))
  graph<-set.vertex.attribute(graph,"x coordinate",value=landscape$x)
  graph<-set.vertex.attribute(graph,"y coordinate",value=landscape$y)
  graph<-set.edge.attribute(graph,"weight",value=distance_mat1[cbind(as.numeric(get.edgelist(graph)[,1]),  as.numeric(get.edgelist(graph)[,2]))])
  
  
  if(components(graph)$no == 1){success<-T}}

envt.v<-0.5*eAMP*(sin((2*pi/ePeriod)*1+1+(landscape$y)*2*pi/1000)+1)


n.removed<-35

ColV2<-c(1,"red","dodgerblue")

pdf(file="./Figures/1. Network figure.pdf",width=7,height=7)
par(mfrow=c(2,2), mar=c(3,3,3,3),pty='s',oma=c(1,1,1,1))

ColV<-heat.colors(100)[1+(envt.v*99)]

#ColV[order(betweenness(graph),decreasing=T)[1]]<-ColV2[2]
#ColV[order(betweenness(graph),decreasing=F)[1]]<-ColV2[3]

shapeVr<-shapeV
shapeVr[order(betweenness(graph),decreasing=T)[1]]<-"uptriangle"
shapeVr[order(betweenness(graph),decreasing=F)[1]]<-"downtriangle"
plot.igraph(graph,layout=as.matrix(landscape), vertex.color=ColV, vertex.size=5000,vertex.label=NA, rescale=F, ylim=c(0,1000),xlim=c(0,1000),vertex.shape=shapeVr)
title("Intact\nnetwork",line = 1)

weightedgraph<-graph
landscape2<-landscape
for(d in 1:n.removed){
  patch.delete<-order(betweenness(weightedgraph),decreasing=F)[1]
  landscape2<-landscape2[-patch.delete,]
  weightedgraph<-delete.vertices(weightedgraph,patch.delete)
}
GrayColVr<-ColV[as.numeric(V(weightedgraph)$name)]
#GrayColVr[order(betweenness(weightedgraph),decreasing=F)[1]]<-"dodgerblue"
shapeVr<-rep("circle",length(V(weightedgraph)$name))
shapeVr[order(betweenness(weightedgraph),decreasing=F)[1]]<-"downtriangle"
plot.igraph(weightedgraph,layout=as.matrix(landscape2), vertex.color=GrayColVr, vertex.size=5000,vertex.label=NA, rescale=F, ylim=c(0,1000),xlim=c(0,1000),vertex.shape=shapeVr)

title("Remove\nmin betweenness",line = 1)

weightedgraph<-graph
landscape2<-landscape
for(d in 1:n.removed){
  patch.delete<-sample(vcount(weightedgraph),1)
  landscape2<-landscape2[-patch.delete,]
  weightedgraph<-delete.vertices(weightedgraph,patch.delete)
}

GrayColVr<-ColV[as.numeric(V(weightedgraph)$name)]
#GrayColVr[order(betweenness(weightedgraph),decreasing=F)[1]]<-"dodgerblue"
shapeVr<-rep("circle",length(V(weightedgraph)$name))
plot.igraph(weightedgraph,layout=as.matrix(landscape2), vertex.color=GrayColVr, vertex.size=5000,vertex.label=NA, rescale=F, ylim=c(0,1000),xlim=c(0,1000),vertex.shape=shapeVr)

title("Remove\nmin random",line = 1)

weightedgraph<-graph
landscape2<-landscape
for(d in 1:n.removed){
  patch.delete<-order(betweenness(weightedgraph),decreasing=T)[1]
  landscape2<-landscape2[-patch.delete,]
  weightedgraph<-delete.vertices(weightedgraph,patch.delete)
}
GrayColVr<-ColV[as.numeric(V(weightedgraph)$name)]
#GrayColVr[order(betweenness(weightedgraph),decreasing=F)[1]]<-"dodgerblue"
shapeVr<-rep("circle",length(V(weightedgraph)$name))
shapeVr[order(betweenness(weightedgraph),decreasing=T)[1]]<-"uptriangle"
plot.igraph(weightedgraph,layout=as.matrix(landscape2), vertex.color=GrayColVr, vertex.size=5000,vertex.label=NA, rescale=F, ylim=c(0,1000),xlim=c(0,1000),vertex.shape=shapeVr)

title("Remove\nmax betweenness",line = 1)
dev.off()
