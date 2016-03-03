rewire<-function(numPatches,numLinks,numEdgesRewired){
   numRewiring<-numEdgesRewired
   numNodes<-numPatches
   numEdges<-numLinks
   lattice <- graph.lattice(numNodes, nei = 2, circular = TRUE)
   latticerewired<-lattice                               
   lattice_edges<-get.edgelist(lattice)
   latticerewired_edges<-lattice_edges
   donedone<-c(0)[rep(c(1), times=numEdges)]  # keep track of which edges have already been rewired so as not to rewire them again
   if(numRewiring>0)
   {
      for(i in 1:numRewiring)
      {
          repeat{                                            ##prevent rewiring the same link multiple times
             rewired_edge_temp<-as.integer(runif(1,1,numEdges))
             if(donedone[rewired_edge_temp]==0){
                node1_temp<-latticerewired_edges[rewired_edge_temp,1]
                repeat{                                            ##prevent self-link loops
                   random_node<-as.integer(runif(1,0,numNodes))
                   node2_temp<-random_node
                   continue<-FALSE
                   if(node2_temp!=node1_temp){
                      #four indices to check for 4 cases of whether link already exists (in each direction (2) and in each lattice (2))
                      indx1<-1
                      indx2<-1
                      indx3<-1
                      indx4<-1
                      indx1<-which(latticerewired_edges[,1]==node2_temp & latticerewired_edges[,2]==node1_temp)  #check if link already exists in one direction in rewired lattice
                      if(length(indx1)==0) indx2<-which(lattice_edges[,1]==node2_temp & lattice_edges[,2]==node1_temp) #check if link already exists in one direction in original lattice
                      if(length(indx2)==0) indx3<-which(latticerewired_edges[,1]==node1_temp & latticerewired_edges[,2]==node2_temp) #check if link already exists in other direction in rewired lattice
                      if(length(indx3)==0) indx4<-which(lattice_edges[,1]==node1_temp & lattice_edges[,2]==node2_temp) #check if link already exists in other direction in original lattice
                      if(length(indx4)==0) continue<-TRUE
                   }
                  if(continue==TRUE) break
                }
          break
          }     
      }
      rewired_edge<-rewired_edge_temp
      donedone[rewired_edge]<-1
      node1<-node1_temp
      node2<-node2_temp
      latticerewired_edges[rewired_edge,1]<-node1
      latticerewired_edges[rewired_edge,2]<-node2
      }
      latticerewired <- graph.data.frame(as.data.frame(latticerewired_edges), directed=FALSE)
   }
   latticerewired$layout<-layout.circle
   #plot(latticerewired)
   return(latticerewired)
   }