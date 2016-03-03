#Function to add edge and node spatial weights to a graph
addweights<-function(unweightedgraph,numLinks,numNodes){
  weightedgraph<-unweightedgraph
   #determine node coordinates
   coordinates<-matrix(NA,2,numNodes)
   coordinates[,1]<-c(1,0)
   for (i in 2:numNodes){
     coordinates[,i]<-c(cos((i-1)*pi/(numNodes/2)),sin((i-1)*pi/(numNodes/2)))
     }

   #plot(coordinates[1,],coordinates[2,], xlim=c(-1.2,1.2), ylim=c(-1.2,1.2))
   
   #assign node coordinates as node weights
   weightedgraph<-set.vertex.attribute(weightedgraph,"x coordinate",value=coordinates[1,])
   weightedgraph<-set.vertex.attribute(weightedgraph,"y coordinate",value=coordinates[2,])

   #determine edge distances (topological and geodesic)
   numEdges<-numLinks
   topolength<-c(0)[rep(c(1), times=numEdges-1)]
   relativeposition<-c(0)[rep(c(1), times=numEdges-1)]
   geolength<-c(0)[rep(c(1), times=numEdges-1)]
   
   distances<-c(rep(NA,numNodes-1))
   for (i in 1:(numNodes/2)){
     distances[i]<-sqrt((coordinates[1,1] - coordinates[1,i+1])^2+(coordinates[2,1]-coordinates[2,i+1])^2)
   }
   for (i in (numNodes/2+1):(numNodes-1)){
     dif_from_med<-i-numNodes/2
     distances[i]<-distances[numNodes/2-dif_from_med]
     }
   
   #distance1 <- sqrt((coordinate1[1] - coordinate2[1])^2 + (coordinate1[2] - coordinate2[2])^2)
   #distance2 <- sqrt((coordinate1[1] - coordinate3[1])^2 + (coordinate1[2] - coordinate3[2])^2)
   #distance3 <- sqrt((coordinate1[1] - coordinate4[1])^2 + (coordinate1[2] - coordinate4[2])^2)
   #distance4 <- sqrt((coordinate1[1] - coordinate5[1])^2 + (coordinate1[2] - coordinate5[2])^2)
   #distance5 <- sqrt((coordinate1[1] - coordinate6[1])^2 + (coordinate1[2] - coordinate6[2])^2)
   
   graph_edges<-get.edgelist(weightedgraph, names=FALSE)
#   graph_edges_10<-get.edgelist(weightedgraph, names=FALSE)
#   graph_edges_10[graph_edges_10 == 0]<-10
   shortest_paths<-shortest.paths(weightedgraph)
#   for(i in 1:(numEdges-1)){
#     topolength[i]<-shortest_paths[graph_edges_10[i,1],graph_edges_10[i,2]]
#   }

   for(i in 1:numEdges){
     relativeposition[i]<-abs(graph_edges[i,1]-graph_edges[i,2])
   }
   
   for (i in 1:numEdges){
     geolength[i]<-distances[relativeposition[i]]
   }

   #geolength[topolength==1]<-distances[1]
   #geolength[topolength==2]<-distances[2]
   #geolength[topolength==3]<-distances[3]
   #geolength[topolength==4]<-distances[4]
   #geolength[topolength==5]<-distances[5]
   #geolength[topolength==6]<-distances[4]
   #geolength[topolength==7]<-distances[3]
   #geolength[topolength==8]<-distances[2]
   #geolength[topolength==9]<-distances[1]
   
#   weightedgraph<-set.edge.attribute(weightedgraph,"topodistance",value=topolength)
   weightedgraph<-set.edge.attribute(weightedgraph,"weight",value=geolength)
  
  return(weightedgraph)
  }