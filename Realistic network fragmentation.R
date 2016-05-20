require(igraph)
require(dplyr)
require(ggplot2)
require(tidyr)
require(doParallel)
require(foreach)

#set up parallel####
cl<-makeCluster(detectCores())
registerDoParallel(cl)
getDoParWorkers()

#simulation code####
reps<-5
print.plots<-T # set this to true if you want to see the network as the sim runs - it makes it slower

nSpecies<-30
numCom<-100
dispV<-c(0.0005,0.005,0.015)

rInput<-150 #resource input
rLoss<-10 #resource loss 
eff<-0.2 #conversion efficiency
mort<-0.2 #mortality
Ext<- 0.1 #extinction Threshold

ePeriod<-40000#40000 #period of env sinusoidal fluctuations
eAMP<-1 #amplitude of envrionment sinusoidal fluctuations

drop_length<-ePeriod/4

Tmax<-100000+drop_length*(numCom-0) #number of time steps in Sim
Tdata<- seq(1, Tmax)
DT<- 0.08 # % size of discrete "time steps"
sampleV<-seq(102000,Tmax,by=2000)
removeV<-c("Max betweenness","Min betweenness","Random")

Meta_dyn_reps<-data.frame(Rep=rep(1:reps,each=(numCom-0)*3),Dispersal=rep(dispV,each=reps*(numCom-0)*3),Patch_remove=rep(factor(removeV,levels = c("Min betweenness","Random","Max betweenness"),ordered = T),each=length(dispV)*reps*(numCom-0)*3),Patches=NA,Dynamic=rep(factor(c("Species sorting", "Mass effects", "Base growth"),levels = c("Base growth","Species sorting","Mass effects")),each=numCom-0),Proportion=NA)
SIH_data_reps<-data.frame(Rep=rep(1:reps,each=(numCom-0)),Dispersal=rep(dispV,each=reps*(numCom-0)),Patch_remove=rep(factor(removeV,levels = c("Min betweenness","Random","Max betweenness"),ordered = T),each=length(dispV)*reps*(numCom-0)),Patches=NA,Regional_SR=NA,Local_SR=NA,Biomass=NA,Occupancy=NA,Regional_CV=NA,Local_CV=NA)
Component_data_reps<-data.frame(Rep=rep(1:reps,each=(numCom-0)),Dispersal=rep(dispV,each=reps*(numCom-0)),Patch_remove=rep(factor(removeV,levels = c("Min betweenness","Random","Max betweenness"),ordered = T),each=length(dispV)*reps*(numCom-0)),Patches=NA,Component_num=NA,Component_size=NA, Component_range=NA)

#the simulation function####  
SIH_frag<-function(){
  Meta_dyn_r1<-data.frame(Rep=r,Dispersal=rep(dispV,each=(numCom-0)*3),Patch_remove=rep(factor(removeV,levels = c("Min betweenness","Random","Max betweenness"),ordered = T),each=length(dispV)*(numCom-0)*3),Patches=NA,Dynamic=rep(factor(c("Species sorting", "Mass effects", "Base growth"),levels = c("Base growth","Species sorting","Mass effects")),each=numCom-0),Proportion=NA)
  SIH_data_r1<-data.frame(Rep=r,Dispersal=rep(dispV,each=(numCom-0)),Patch_remove=rep(factor(removeV,levels = c("Min betweenness","Random","Max betweenness"),ordered = T),each=length(dispV)*(numCom-0)),Patches=NA,Regional_SR=NA,Local_SR=NA,Biomass=NA,Occupancy=NA,Regional_CV=NA,Local_CV=NA)
  Component_data_r1<-data.frame(Rep=r,Dispersal=rep(dispV,each=(numCom-0)),Patch_remove=rep(factor(removeV,levels = c("Min betweenness","Random","Max betweenness"),ordered = T),each=length(dispV)*(numCom-0)),Patches=NA,Component_num=NA,Component_size=NA, Component_range=NA)
  
  for(i in 1:length(dispV)){
    disp<-dispV[i]
    
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
    
    
    #vectors####
    eOptimum<-1-seq(0,eAMP, by=eAMP/(nSpecies-1)) #species environmental optima
    
    calc.immigration <- function(N,a,dispersal_matrix) dispersal_matrix%*%N*rep(a,each=length(R))
    
    for(j in 1:3){
      landscape_run<-landscape
      graph_run<-graph
      dispersal_matrix <- apply(connections, 1, function(x) x/sum(x)) #divides the d_exp matrix by the column sums to make it a conservative dispersal matrix
      
      Prod<-array(NA,dim=c(numCom,nSpecies,length(sampleV)))
      Abund<-Prod
      
      N0<-N<- matrix(10,ncol=nSpecies,nrow=numCom) # Community x Species abundance matrix
      R0<-R<-rep(10*(nSpecies/10),numCom)
      
      Meta_dyn<-data.frame(Species_sorting=rep(NA,length(sampleV)),Mass_effects=NA,Base_growth=NA,Patches=NA)
      Species_data<-array(NA,dim=c(length(sampleV),nSpecies,2),dimnames = list(sampleV,1:nSpecies,c("Abundance","Occupancy")))
      Components<-data.frame(Number_components=rep(NA, length(sampleV)),Component_size=NA,Component_envt_range=NA)
      
      for(TS in 1:Tmax){
        #print(TS)
        Immigrants<-calc.immigration(N,disp,dispersal_matrix)
        envt.v<-0.5*eAMP*(sin((2*pi/ePeriod)*TS+1+(landscape_run$y)*2*pi/1000)+1)
        #plot.igraph(graph_run,layout=as.matrix(landscape), vertex.color=heat.colors(100)[1+(envt.v*99)], vertex.size=10)
        if(is.null(rownames(dispersal_matrix))){
          consume <- 0.1*(1.5-abs(sapply(eOptimum,'-',envt.v)))
          Nt <- N*(1+DT*(eff*R*consume - disp - mort)) + DT*Immigrants #abundance step
          Rt <- DT*rInput+R*(1-DT*(rLoss + sum(consume*N))) #resource step    
          
          Nt0 <- N0*(1+DT*(eff*R0*consume - mort))
          Rt0 <- DT*rInput+R0*(1-DT*(rLoss + sum(consume*N0)))} else { #resource step  
            consume <- 0.1*(1.5-abs(sapply(eOptimum,'-',envt.v)))
            Nt <- N*(1+DT*(eff*R*consume - disp - mort)) + DT*Immigrants #abundance step
            Rt <- DT*rInput+R*(1-DT*(rLoss + rowSums(consume*N))) #resource step    
            
            Nt0 <- N0*(1+DT*(eff*R0*consume - mort))
            Rt0 <- DT*rInput+R0*(1-DT*(rLoss + rowSums(consume*N0)))
          }
        
        if(sum(TS==sampleV)==1){
          sample_id<-which(TS==sampleV)
          Components$Number_components[sample_id]<-components(graph_run)$no
          Components$Component_size[sample_id]<-mean(components(graph_run)$csize)
          members<-components(graph_run)$membership
          envt.ranges<-sapply(unique(members),function(x){range(envt.v[members==x])})
          Components$Component_envt_range[sample_id]<-mean(envt.ranges[2,]-envt.ranges[1,])
          
          if(is.null(rownames(N))){
            Prod[as.numeric(names(dispersal_matrix)),,sample_id] <- eff*consume*R*N
            Abund[as.numeric(names(dispersal_matrix)),,sample_id] <- N
            
            fitness<-((N*(1+DT*(eff*R*consume - disp - mort)))-N)*(Nt>Ext)
            fitness_w_disp<-((N*(1+DT*(eff*R*consume - disp - mort)) + DT*Immigrants)-N)*(Nt>Ext)
            fitness0<-(N0*(1+DT*(eff*R0*consume - mort))-N0)*(Nt0>Ext)
            home_prod<-mean(rowSums(fitness_w_disp*(fitness>0)))
            disp_prod_ME<-mean(rowSums(fitness_w_disp*(fitness<0 & fitness_w_disp>=0)))
            
            base_prod<-mean(sum(fitness0*(fitness0>0)))
            total_prod<-home_prod+disp_prod_ME
            
            home_prod_prop<-home_prod/total_prod
            SS_prod<-home_prod-base_prod
            SS_prod[SS_prod<0]<-0
            if(mean(sum(N>0))<=1){SS_prod<-0}
            SS<-(SS_prod/home_prod)*home_prod_prop
            SS[is.nan(SS)]<-0
            if(total_prod==0){SS<-NA}
            Meta_dyn$Species_sorting[sample_id]<-SS
            
            ME<-(disp_prod_ME)/total_prod
            ME[is.nan(ME)]<-0
            if(total_prod==0){ME<-NA}
            Meta_dyn$Mass_effects[sample_id]<-ME
            
            BP<-home_prod_prop*(1-(SS_prod/home_prod))
            BP[is.nan(BP)]<-0
            if(total_prod==0){BP<-NA}
            Meta_dyn$Base_growth[sample_id]<-BP
            
            Meta_dyn$Patches[sample_id]<-1
            
            Species_data[sample_id,,1]<-N
            Species_data[sample_id,,2]<-N>0
          } else{
            Prod[as.numeric(rownames(N)),,sample_id] <- eff*consume*R*N
            Abund[as.numeric(rownames(N)),,sample_id] <- N
            
            fitness<-((N*(1+DT*(eff*R*consume - disp - mort)))-N)*(Nt>Ext)
            fitness_w_disp<-((N*(1+DT*(eff*R*consume - disp - mort)) + DT*Immigrants)-N)*(Nt>Ext)
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
            Meta_dyn$Species_sorting[sample_id]<-SS
            
            ME<-(disp_prod_ME)/total_prod
            ME[is.nan(ME)]<-0
            if(total_prod==0){ME<-NA}
            Meta_dyn$Mass_effects[sample_id]<-ME
            
            BP<-home_prod_prop*(1-(SS_prod/home_prod))
            BP[is.nan(BP)]<-0
            if(total_prod==0){BP<-NA}
            Meta_dyn$Base_growth[sample_id]<-BP
            
            Meta_dyn$Patches[sample_id]<-nrow(N)
            
            Species_data[sample_id,,1]<-colSums(N)
            Species_data[sample_id,,2]<-colSums(N>0)
          }
        }
        
        N <- Nt * (Nt>Ext) # set to 0 if below extinction threshold
        R <- Rt
        
        N0 <- Nt0 * (Nt0>Ext) # set to 0 if below extinction threshold
        R0 <- Rt0
        
        if(max(TS==seq(100000+drop_length,Tmax-1,by=drop_length))){
          if(j==1){btw<-betweenness(graph_run)
          if(sum(btw)==0){
            patch.delete<-order(degree(graph_run),decreasing = T)[1]
          } else{patch.delete<-order(btw,decreasing = T)[1] }
          } else{
            if(j==2){patch.delete<-order(betweenness(graph_run),decreasing=F)[1]} else{
              patch.delete<-sample(nrow(N),1)}}    
          graph_run<-delete.vertices(graph_run,patch.delete)
          landscape_run<-landscape_run[-patch.delete,]
          
          dispersal_matrix <- apply(get.adjacency(graph_run,sparse=F), 1, function(x) x/sum(x)) #divides the d_exp matrix by the column sums to make it a conservative dispersal matrix
          dispersal_matrix[is.nan(dispersal_matrix)]<-0
          
          envt.v<-0.5*eAMP*(sin((2*pi/ePeriod)*TS+1+(landscape_run$y)*2*pi/1000)+1)
          if(print.plots==T){
            plot.igraph(graph_run,layout=as.matrix(landscape_run), vertex.color=heat.colors(100)[1+(envt.v*99)], vertex.size=5000,vertex.label=NA, rescale=F, ylim=c(0,1000),xlim=c(0,1000))
            }
            N<-N[-patch.delete,]
            R<-R[-patch.delete]
            N0<-N0[-patch.delete,]
            R0<-R0[-patch.delete]
          }  
        } 
        
        L_Bmass<-colMeans(apply(Abund,3,rowSums),na.rm=T)
        L_Bmass_sep<-data.frame(t(apply(Abund,3,rowSums)))
        R_Bmass<-apply(Abund,3,sum,na.rm=T)
        R_SR<-colSums(apply(Abund,3,colSums, na.rm=T)>0)
        L_SR<-colMeans(apply((Abund>0),3,rowSums),na.rm=T)
        L_Occ<-colMeans(apply((Abund>0),3,rowSums)>0,na.rm=T)
        
        cv<-function(x){sd(x,na.rm=T)/mean(x,na.rm=T)}
        
        CVdf<-cbind(L_Bmass_sep,data.frame(R_Bmass=R_Bmass,Patches=rep(numCom:1,each=drop_length/2000))) %>%
          group_by(Patches) %>%
          summarise_each(funs(cv))
        L_CV<-rowMeans(CVdf[,2:31],na.rm=T)
        R_CV<-CVdf$R_Bmass
        
        SIH_data_means<-data.frame(R_SR=R_SR,L_SR=L_SR,L_Bmass=L_Bmass,L_Occ=L_Occ,Patches=numCom-colMeans(apply(is.na(Abund),3,colSums))) %>%
          group_by(Patches) %>%
          summarise_each(funs(mean))
        SIH_data_means$R_CV<-R_CV
        SIH_data_means$L_CV<-L_CV
        
        Component_data_means<-data.frame(Patches=numCom-colMeans(apply(is.na(Abund),3,colSums)),Component_num=Components$Number_components,Component_size=Components$Component_size,Component_range=Components$Component_envt_range)%>%
          group_by(Patches) %>%
          summarise_each(funs(mean))
        
        SIH_data_r1[SIH_data_r1$Dispersal==dispV[i] & SIH_data_r1$Patch_remove==removeV[j],-c(1:3)]<-SIH_data_means
        Component_data_r1[SIH_data_r1$Dispersal==dispV[i] & SIH_data_r1$Patch_remove==removeV[j],-c(1:3)]<-Component_data_means
        
        mean.df<-summarise(group_by(Meta_dyn,Patches),Species_sorting=mean(Species_sorting,na.rm=T),Mass_effects=mean(Mass_effects,na.rm=T),Base_growth=mean(Base_growth,na.rm=T))
        Meta.dyn.long<-gather(mean.df,key = Dynamic,value=Proportion,-Patches)
        
        Meta_dyn_r1[Meta_dyn_r1$Dispersal==dispV[i] & Meta_dyn_r1$Patch_remove==removeV[j],c(4,6)]<-Meta.dyn.long[,-2]
      }}
    return(list(Meta_dyn_r1,SIH_data_r1,Component_data_r1))
  }
  
  #run simulation function in parallel
  Sim_data_parallel<-foreach(r = 1:reps,.packages=c("igraph","dplyr","tidyr")) %dopar% SIH_frag()
  for(r in 1:reps){
    Sim_data<-Sim_data_parallel[[r]]
    Meta_dyn_reps[Meta_dyn_reps$Rep==r,]<-Sim_data[[1]]
    SIH_data_reps[SIH_data_reps$Rep==r,]<-Sim_data[[2]]
    Component_data_reps[Component_data_reps$Rep==r,]<-Sim_data[[3]]
  }  
  
  save(Meta_dyn_reps,Component_data_reps,SIH_data_reps,file="Fragmentation.RData")