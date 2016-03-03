create_random_net<-function(numPatches,numLinks){
   repeat{
      randomgraph<-erdos.renyi.game(numPatches,numLinks,type="gnm",directed=FALSE,loops=FALSE)
      if(is.connected(randomgraph)==TRUE) break
   }
   randomgraph$layout<-layout.circle
#   plot(randomgraph)
#   ecount(randomgraph)
   return(randomgraph)
   }