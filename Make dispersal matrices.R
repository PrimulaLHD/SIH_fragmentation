library(igraph)
reps<-100
numCom<-100
CS_exe<-'C:/"Program Files"/Circuitscape/cs_run.exe'

dispersal_mats<-list()
net_list<-list()
landscape_list<-list()


for(r in 26:reps){
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
    
    
    if(components(graph)$no == 1 & sum(duplicated(landscape))==0){success<-T}}
  name<-paste("Network",r,sep="_")
  net_list[[name]]<-graph
  
  name<-paste("Landscape",r,sep="_")
  landscape_list[[name]]<-landscape
  
  for(j in 1:3){
    landscape_run<-landscape
    graph_run<-graph
    for(p in 1:100){
      if(max(degree(graph_run))<2){
        dispersal_matrix<-get.adjacency(graph_run,sparse=F)
      } else{
        graph_circuit<-data.frame(get.edgelist(graph_run), E(graph_run)$weight)
        write.table(graph_circuit, paste("./Circuits/network_graph",r,".txt", sep="_"), row.names = FALSE, col.names = FALSE, quote = FALSE)
        write.table(V(graph_run)$name, paste("./Circuits/network_focal_nodes",r,".txt", sep="_"), row.names = FALSE, col.names = FALSE, quote = FALSE)
        CS_ini <- c("[circuitscape options]",            
                    "data_type = network",
                    "scenario = pairwise",
                    paste(c("point_file =",
                            "habitat_file =",
                            "output_file ="),
                          c(paste("./Circuits/network_focal_nodes",r,".txt", sep="_"),
                            paste("./Circuits/network_graph",r,".txt", sep="_"),
                            paste("./Circuits/CS",r,".out",sep="_"))))
        
        writeLines(CS_ini, paste("./Circuits/my_ini",r,".ini", sep="_"))
        CS_run <- paste(CS_exe, paste("./Circuits/my_ini",r, ".ini", sep="_")) # Make the cmd
        system(CS_run)
        d <- read.table(paste("./Circuits/CS",r,"_resistances.out",sep="_"), row.names=1,header=TRUE)
        d[d<0]<-Inf
        #d<-shortest.paths(graph_run, mode="all", weights=NULL, algorithm="automatic")
        #d_exp<-exp(-0.002*d) - diag(nrow(d))  #dispersal kernel function of the d matrix
        d_exp<-exp(-0.05*d) - diag(nrow(d))
        d_exp[d_exp==1]<-0
        dispersal_matrix <- apply(d_exp, 1, function(x) x/sum(x)) #divides the d_exp matrix by the column sums to make it a conservative dispersal matrix
        dispersal_matrix[is.nan(dispersal_matrix)]<-0
        rownames(dispersal_matrix)<-rownames(landscape_run)
      }
      name<-paste("DM",r,j,101-p,sep="_")
      dispersal_mats[[name]]<-dispersal_matrix
      
      if(j==1){btw<-betweenness(graph_run)
      if(sum(btw)==0){
        patch.delete<-order(degree(graph_run),decreasing = T)[1]
      } else{patch.delete<-order(btw,decreasing = T)[1] }
      } else{
        if(j==2){patch.delete<-order(betweenness(graph_run),decreasing=F)[1]} else{
          patch.delete<-sample(vcount(graph_run),1)}}    
      graph_run<-delete.vertices(graph_run,patch.delete)
      landscape_run<-landscape_run[-patch.delete,]
    }
  }
}

save(dispersal_mats, net_list,landscape_list,file="Dispersal matrices.RData")
