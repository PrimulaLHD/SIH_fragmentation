require(igraph)
require(dplyr)
require(ggplot2)
require(tidyr)
require(doParallel)
require(foreach)

SIH_function<-function(dispersal=0.005,species=30,numCom=100){
  #Constants####
  #N<- matrix(10,ncol=species,nrow=numCom) # Community x Species abundance matrix
  N<-matrix(10,numCom,species)
  R<-rep(10*(species/10),numCom) #Initial resources
  N0<-N
  R0<-R
  
  rInput<-150 #resource input
  rLoss<-10 #resource loss 
  eff<-0.2 #conversion efficiency
  mort<-0.2 #mortality
  Ext<- 0.1 #extinction Threshold
  
  ePeriod<-40000 #period of env sinusoidal fluctuations
  eAMP<-1 #amplitude of envrionment sinusoidal fluctuations
  
  Tmax<-140000 #number of time steps in Sim
  DT<- 0.08 # % size of discrete "time steps" - this is the Euler value
  
  #vectors####
  eOptimum<-1-seq(0,eAMP, by=eAMP/(species-1)) #species environmental optima
  
  #network####
  success<-F
  while(!success){
    landscape<-round(data.frame(x = runif(numCom, min = 1, max = 1000), y = runif(numCom, min = 1, max = 1000)))
    distance_mat1<-as.matrix(dist(landscape,method = "euclidean",diag = T,upper=T))
    
    dist.lower<-dist(landscape)
    weights<-dist.lower[dist.lower<200 & dist.lower>0]
    distance_mat<-1*(distance_mat1<200)
    diag(distance_mat)<-0
    connections<-distance_mat
    distance_mat[upper.tri(distance_mat)]<-0
    
    graph<-as.undirected(graph.adjacency(distance_mat))
    graph<-set.vertex.attribute(graph,"x coordinate",value=landscape$x)
    graph<-set.vertex.attribute(graph,"y coordinate",value=landscape$y)
    graph<-set.edge.attribute(graph,"weight",value=weights)
    
    if(components(graph)$no == 1){success<-T}}
  
  colV<-heat.colors(100)[round(landscape$y/10)]
  plot.igraph(graph,layout=as.matrix(landscape), vertex.color=colV, vertex.size=10)
  
  #dispersal conditions####
  dispersal_matrix <- apply(connections, 1, function(x) x/sum(x)) #divides the d_exp matrix by the column sums to make it a conservative dispersal matrix
  
  calc.immigration <- function(N,a,dispersal_matrix) dispersal_matrix%*%N*rep(a,each=numCom)
  
  Prod<-matrix(NA,species*numCom,40000)
  Abund<-Prod
  
  Meta_dyn<-data.frame(Species_sorting=rep(NA,40000),Mass_effects=NA,Base_growth=NA)
  Species_data<-array(NA,dim=c(40000,species,2),dimnames = list(1:40000,1:species,c("Abundance","Occupancy")))
  
  for(TS in 1:Tmax){
    envt.v<-0.5*eAMP*(sin((2*pi/ePeriod)*TS+1+(landscape$y)*2*pi/1000)+1)
    consume <- 0.1*(1.5-abs(sapply(eOptimum,'-',envt.v)))
    Immigrants<-calc.immigration(N,dispersal,dispersal_matrix)
    Nt <- N*(1+DT*(eff*R*consume - dispersal - mort)) + DT*Immigrants
    
    Immigrants0<-calc.immigration(N0,0,dispersal_matrix)
    Nt0 <- N0*(1+DT*(eff*R0*consume -0 - mort)) + DT*Immigrants0
    
    Rt <- DT*rInput+R*(1-DT*(rLoss + rowSums(consume*N))) #resource step   
    Rt0 <- DT*rInput+R0*(1-DT*(rLoss + rowSums(consume*N0))) #resource step  
    
    if(TS>=100000){
      Prod[,(TS-100000)] <- c(t(eff*consume*R*N))
      Abund[,(TS-100000)] <- c(t(N))
      
      fitness<-((N*(1+DT*(eff*R*consume - dispersal - mort)))-N)*(Nt>Ext)
      fitness_w_disp<-((N*(1+DT*(eff*R*consume - dispersal - mort)) + DT*Immigrants)-N)*(Nt>Ext)
      fitness0<-(N0*(1+DT*(eff*R0*consume - mort))-N0)*(Nt0>Ext)
      home_prod<-mean(rowSums(fitness_w_disp*(fitness>0)))
      disp_prod_ME<-mean(rowSums(fitness_w_disp*(fitness<0 & fitness_w_disp>=0)))
      
      base_prod<-mean(rowSums(fitness0*(fitness0>0)))
      total_prod<-home_prod+disp_prod_ME
      
      home_prod_prop<-home_prod/total_prod
      SS_prod<-home_prod-base_prod
      SS_prod[SS_prod<0]<-0
      if(mean(rowSums(N>0))<=1){SS_prod<-0}
      SS<-(SS_prod/home_prod)*home_prod_prop
      SS[is.nan(SS)]<-0
      if(total_prod==0){SS<-NA}
      Meta_dyn$Species_sorting[(TS-100000)]<-SS
      
      ME<-(disp_prod_ME)/total_prod
      ME[is.nan(ME)]<-0
      if(total_prod==0){ME<-NA}
      Meta_dyn$Mass_effects[(TS-100000)]<-ME
      
      BP<-home_prod_prop*(1-(SS_prod/home_prod))
      BP[is.nan(BP)]<-0
      if(total_prod==0){BP<-NA}
      Meta_dyn$Base_growth[(TS-100000)]<-BP
      
      Species_data[(TS-100000),,1]<-colSums(N)
      Species_data[(TS-100000),,2]<-colSums(N>0)
    }
    N <- Nt * (Nt>Ext) # set to 0 if below extinction threshold
    R <- Rt
    
    N0 <- Nt0 * (Nt0>Ext) # set to 0 if below extinction threshold
    R0 <- Rt0
  } 
  
  Prod<-array(t(Prod),dim=c(40000,species,numCom))
  Prod<-Prod[seq(1,40000,100),,]
  
  Abund<-array(t(Abund),dim=c(40000,species,numCom))
  Abund<-Abund[seq(1,40000,100),,]
  #matplot(Abund[,,1], type ='l', lty=1)
  return(list(Prod=Prod,Abund=Abund,Meta_dyn=Meta_dyn,Spec_data=apply(Species_data,3,colMeans)))
}

vect<-c(0.0001,0.00015,0.00025,0.0005,0.00075)
dispV<-c(vect,vect*10,vect*100,vect*1000,1)
dispV<-dispV[-c(17:20)]

dispV<-c(0.0001,0.0005,0.001,0.005,0.01,0.05,0.1,0.5,1)

species<-50
reps<-2
Meta_dyn.df<-data.frame(Species_sorting=NA,Mass_effects=NA,Base_growth=NA,Dispersal=dispV,Rep=rep(1:reps,each=length(dispV)))

#make parallel####
cl<-makeCluster(2)
registerDoParallel(cl)
getDoParWorkers()

Sim_data_parallel<-foreach(r = 1:reps,.packages=c("igraph","dplyr","tidyr")) %dopar% sapply(dispV,SIH_function,species=species)

for(r in 1:reps){
  print(r)
  for(i in 1:length(dispV)){
    Meta_dyn.df[Meta_dyn.df$Dispersal==dispV[i] & Meta_dyn.df$Rep==r,1:3]<-colMeans(Sim_data_parallel[[r]]["Meta_dyn",i]$Meta_dyn,na.rm=T)
  }}


stopCluster(cl)

save(Meta_dyn.df,file="Meta_dynamics.RData")