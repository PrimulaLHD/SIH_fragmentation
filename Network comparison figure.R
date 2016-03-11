#plot network examples####
require(igraph)
source("./Functions/rewire.R")
source("./Functions/create_random_net.r")
source("./Functions/addweights.r")
source("./Functions/maketriangles.r")

numCom<-30
numLinks<-numCom*2

colfunc <- colorRampPalette(c("white", "green4"))
GrayColV <- colfunc(16)
GrayColV<-c(GrayColV,rev(GrayColV[2:15]))
GrayColV<-GrayColV[c(9:30,1:8)]
shapeV<-rep("circle",numCom)

add.vertex.shape("uptriangle",plot=uptriangle)
add.vertex.shape("downtriangle",plot=downtriangle)

rand<-50
numEdgesRewired<-rand/100*(numCom*2) 
success<-FALSE
while(!success){unweightedgraph<- if(rand==100) create_random_net(numCom, numLinks) else rewire(numCom,numLinks,numEdgesRewired)
success<-length(V(unweightedgraph))==30}
weightedgraph<-addweights(unweightedgraph,numLinks,numCom)
holdgraph50<-weightedgraph

n.removed<-14

ColV<-c(1,"red","dodgerblue")

pdf(file="./Figures/1. Network figure.pdf",width=7,height=7)
par(mfrow=c(2,2), mar=c(3,3,3,3),pty='s',oma=c(2,2,2,2))

  GrayColVr<-GrayColV
  GrayColVr[order(betweenness(holdgraph50),decreasing=T)[1]]<-ColV[2]
  GrayColVr[order(betweenness(holdgraph50),decreasing=F)[1]]<-ColV[3]
  shapeVr<-shapeV
  shapeVr[order(betweenness(holdgraph50),decreasing=T)[1]]<-"uptriangle"
  shapeVr[order(betweenness(holdgraph50),decreasing=F)[1]]<-"downtriangle"
  plot(holdgraph50, ylim=c(-1,1),xlim=c(-1,1), vertex.color= GrayColVr, edge.width=2, 
       edge.curved=E(holdgraph50)$weight<0.42, vertex.size=10, vertex.label=NA,vertex.shape=shapeVr)
  title("Intact\nnetwork",line = 1)
  
  weightedgraph<-holdgraph50
  for(d in 1:n.removed){
    patch.delete<-order(betweenness(weightedgraph),decreasing=F)[1]
    weightedgraph<-delete.vertices(weightedgraph,patch.delete)
  }
  GrayColVr<-GrayColV[as.numeric(V(weightedgraph)$name)]
  GrayColVr[order(betweenness(weightedgraph),decreasing=F)[1]]<-"dodgerblue"
  shapeVr<-shapeV[as.numeric(V(weightedgraph)$name)]
  shapeVr[order(betweenness(weightedgraph),decreasing=F)[1]]<-"downtriangle"
  plot.igraph(weightedgraph,edge.width=2,layout=layout.circle(holdgraph50)[as.numeric(V(weightedgraph)$name),],
              ylim=c(-1,1),xlim=c(-1,1), vertex.color= GrayColVr, rescale=F, 
              edge.curved=E(weightedgraph)$weight<0.42, vertex.size=10, vertex.label=NA,vertex.shape=shapeVr)
  title("Remove\nmin betweenness",line = 1)

  weightedgraph<-holdgraph50
  for(d in 1:n.removed){
    patch.delete<-sample(vcount(weightedgraph),1)
    weightedgraph<-delete.vertices(weightedgraph,patch.delete)
  }
  plot.igraph(weightedgraph, layout=layout.circle(holdgraph50)[as.numeric(V(weightedgraph)$name),],
              ylim=c(-1,1),xlim=c(-1,1), vertex.color= GrayColV[as.numeric(V(weightedgraph)$name)],edge.width=2 , rescale=F, 
              edge.curved=E(weightedgraph)$weight<0.42, vertex.size=10, vertex.label=NA)
  title("Remove\nmin random",line = 1)
  
  weightedgraph<-holdgraph50
  for(d in 1:n.removed){
    btw<-betweenness(weightedgraph)
    if(sum(btw==0)){
      patch.delete<-order(degree(weightedgraph),decreasing = T)[1]
    } else{patch.delete<-order(btw,decreasing = T)[1] }
    weightedgraph<-delete.vertices(weightedgraph,patch.delete)
  }
  GrayColVr<-GrayColV[as.numeric(V(weightedgraph)$name)]
  shapeVr<-shapeV[as.numeric(V(weightedgraph)$name)]
  
  if(sum(btw==0)){
    patch.delete<-order(degree(weightedgraph),decreasing = T)[1]
    GrayColVr[order(degree(weightedgraph),decreasing = T)[1]]<-"#e41a1c" 
    shapeVr[order(degree(weightedgraph),decreasing=T)[1]]<-"uptriangle"
  } else{GrayColVr[order(betweenness(weightedgraph),decreasing=T)[1]]<-"#e41a1c" 
  shapeVr[order(betweenness(weightedgraph),decreasing=T)[1]]<-"uptriangle"
  }
  
  plot.igraph(weightedgraph, layout=layout.circle(holdgraph50)[as.numeric(V(weightedgraph)$name),],
              ylim=c(-1,1),xlim=c(-1,1), vertex.color= GrayColVr,edge.width=2 , rescale=F, 
              edge.curved=E(weightedgraph)$weight<0.42, vertex.size=10, vertex.label=NA,vertex.shape=shapeVr)
  title("Remove\nmax betweenness",line = 1)
dev.off()
