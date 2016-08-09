library(igraph)
library(dplyr)
library(ggplot2)
library(tidyr)
library(ggExtra)
library(RColorBrewer)
library(doParallel)
library(foreach)

load("Dispersal matrices.RData")

#set up parallel####
detectCores()
cl<-makeCluster(32,outfile="")
registerDoParallel(cl)
getDoParWorkers()

#simulation code####
reps<-32

#the simulation function####  
SIH_frag<-function(numCom=100,nSpecies=10){
  dispV<-c(0.0005,0.005,0.015)
  print.plots<-F# set this to true if you want to see the network as the sim runs - it makes it slower
  
  
  rInput<-150 #resource input
  rLoss<-10 #resource loss 
  eff<-0.2 #conversion efficiency
  mort<-0.2 #mortality
  Ext<- 0.1 #extinction Threshold
  
  ePeriod<-40000#40000 #period of env sinusoidal fluctuations
  eAMP<-1 #amplitude of envrionment sinusoidal fluctuations
  
  drop_length<-ePeriod
  
  Tmax<-100000+drop_length*(numCom-0) #number of time steps in Sim
  Tdata<- seq(1, Tmax)
  DT<- 0.08 # % size of discrete "time steps"
  sampleV<-seq(102000,Tmax,by=2000)
  removeV<-c("Max betweenness","Min betweenness","Random")
  
  eOptimum<-1-seq(0,eAMP, by=eAMP/(nSpecies-1)) #species environmental optima
  
  #network####
  name<-paste("Network",r,sep="_")
  graph<-net_list[[name]]
  
  name<-paste("Landscape",r,sep="_")
  landscape<-landscape_list[[name]]
  
  envt.v<-0.5*eAMP*(sin((2*pi/ePeriod)*1+1+(landscape$y)*2*pi/1000)+1)
  
  #plot.igraph(graph,layout=as.matrix(landscape), vertex.color=heat.colors(100)[1+(envt.v*99)], vertex.size=5000,vertex.label=NA, rescale=F, ylim=c(0,1000),xlim=c(0,1000))
  
  for(i in 1:length(dispV)){
    disp<-dispV[i]
    
    
    #dispersal conditions####
    calc.immigration <- function(N,a,dispersal_matrix) dispersal_matrix%*%N*rep(a,each=length(R))
    
    for(j in 1:3){
      landscape_run<-landscape
      graph_run<-graph
      
      name<-paste("DM",r,j,100,sep="_")
      dispersal_matrix<-dispersal_mats[[name]]
      
      Abund<-array(NA,dim=c(numCom,nSpecies,length(sampleV)))
      
      N0<-N<- matrix(10,ncol=nSpecies,nrow=numCom,dimnames = list(1:numCom,1:nSpecies)) # Community x Species abundance matrix
      R0<-R<-matrix(rep(10*(nSpecies/10),numCom),nrow=numCom,ncol = 1,dimnames=list(1:numCom,"R"))
      
      Meta_dyn<-data.frame(Species_sorting=rep(NA,length(sampleV)),Mass_effects=NA,Base_growth=NA,Patches=NA)
      Species_data<-array(NA,dim=c(length(sampleV),nSpecies,2),dimnames = list(sampleV,1:nSpecies,c("Abundance","Occupancy")))
      Components<-data.frame(Number_components=rep(NA, length(sampleV)),Component_size=NA,Component_envt_range=NA)
      
      
      for(TS in 1:Tmax){
        #print(TS)
        Immigrants<-calc.immigration(N,disp,dispersal_matrix)
        envt.v<-0.5*eAMP*(sin((2*pi/ePeriod)*TS+1+(landscape_run$y)*2*pi/1000)+1)
        if(vcount(graph_run)==1){
          consume <- 0.1*(1.5-abs(sapply(eOptimum,'-',envt.v)))
          Nt <- N*(1+DT*(eff*c(R)*consume - disp - mort)) + DT*Immigrants #abundance step
          Rt <- DT*rInput+R*(1-DT*(rLoss + sum(consume*N))) #resource step    
          
          Nt0 <- N0*(1+DT*(eff*c(R0)*consume - mort))
          Rt0 <- DT*rInput+R0*(1-DT*(rLoss + sum(consume*N0)))} else { #resource step  
            consume <- 0.1*(1.5-abs(sapply(eOptimum,'-',envt.v)))
            Nt <- N*(1+DT*(eff*c(R)*consume - disp - mort)) + DT*Immigrants #abundance step
            Rt <- DT*rInput+R*(1-DT*(rLoss + rowSums(consume*N))) #resource step    
            
            Nt0 <- N0*(1+DT*(eff*c(R0)*consume - mort))
            Rt0 <- DT*rInput+R0*(1-DT*(rLoss + rowSums(consume*N0)))
          }
        
        if(sum(TS==sampleV)==1){
          sample_id<-which(TS==sampleV)
          Components$Number_components[sample_id]<-components(graph_run)$no
          Components$Component_size[sample_id]<-mean(components(graph_run)$csize)
          members<-components(graph_run)$membership
          envt.ranges<-sapply(unique(members),function(x){range(envt.v[members==x])})
          Components$Component_envt_range[sample_id]<-mean(envt.ranges[2,]-envt.ranges[1,])
          
          if(vcount(graph_run)==1){
            Abund[as.numeric(rownames(landscape_run)),,sample_id] <- N
            
            fitness<-((N*(1+DT*(eff*c(R)*consume - disp - mort)))-N)*(Nt>Ext)
            fitness_w_disp<-((N*(1+DT*(eff*c(R)*consume - disp - mort)) + DT*Immigrants)-N)*(Nt>Ext)
            fitness0<-(N0*(1+DT*(eff*c(R0)*consume - mort))-N0)*(Nt0>Ext)
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
            Abund[as.numeric(rownames(N)),,sample_id] <- N
            
            fitness<-((N*(1+DT*(eff*c(R)*consume - disp - mort)))-N)*(Nt>Ext)
            fitness_w_disp<-((N*(1+DT*(eff*c(R)*consume - disp - mort)) + DT*Immigrants)-N)*(Nt>Ext)
            fitness0<-(N0*(1+DT*(eff*c(R0)*consume - mort))-N0)*(Nt0>Ext)
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
        
        if(max(TS==seq(100000+drop_length,Tmax-1,by=drop_length))==1){
          dropV<-seq(100000+drop_length,Tmax-1,by=drop_length)
          
          name<-paste("DM",r,j,100-which(dropV==TS),sep="_")
          dispersal_matrix<-dispersal_mats[[name]]
          
          landscape_run<-landscape[rownames(dispersal_matrix),]
          graph_run<-induced_subgraph(graph,rownames(dispersal_matrix))
          
          envt.v<-0.5*eAMP*(sin((2*pi/ePeriod)*TS+1+(landscape_run$y)*2*pi/1000)+1)
          if(print.plots==T){
            plot.igraph(graph_run,layout=as.matrix(landscape_run), vertex.color=heat.colors(100)[1+(envt.v*99)], vertex.size=5000,vertex.label=NA, rescale=F, ylim=c(0,1000),xlim=c(0,1000))
          }
          
          N<-N[rownames(dispersal_matrix),]
          R<-matrix(R[rownames(dispersal_matrix),],ncol = 1,dimnames = list(rownames(dispersal_matrix),"R"))
          N0<-N0[rownames(dispersal_matrix),]
          R0<-matrix(R0[rownames(dispersal_matrix),],ncol = 1,dimnames = list(rownames(dispersal_matrix),"R"))
        }  
      } 
      
      L_Bmass<-colMeans(apply(Abund,3,rowSums),na.rm=T)
      L_Bmass_sep<-data.frame(t(apply(Abund,3,rowSums)))
      R_Bmass<-apply(Abund,3,sum,na.rm=T)
      R_SR<-colSums(apply(Abund,3,colSums, na.rm=T)>0)
      L_SR<-colMeans(apply((Abund>0),3,rowSums),na.rm=T)
      L_Occ<-colMeans(apply((Abund>0),3,rowSums)>0,na.rm=T)
      
      cv<-function(x){sd(x,na.rm=T)/mean(x,na.rm=T)}
      
      L_Bmass_sep[L_Bmass_sep==0]<-NA
      
      CVdf<-cbind(L_Bmass_sep,data.frame(R_Bmass=R_Bmass,Patches=rep(numCom:1,each=drop_length/2000))) %>%
        group_by(Patches) %>%
        summarise_each(funs(cv))
      L_CV<-rowMeans(CVdf[,2:101],na.rm=T)
      R_CV<-CVdf$R_Bmass
      
      SIH_data_means<-data.frame(R_SR=R_SR,L_SR=L_SR,L_Bmass=L_Bmass,L_Occ=L_Occ,Patches=numCom-colMeans(apply(is.na(Abund),3,colSums))) %>%
        group_by(Patches) %>%
        summarise_each(funs(mean(.,na.rm=T)))
      SIH_data_means$R_CV<-R_CV
      SIH_data_means$L_CV<-L_CV
      
      Component_data_means<-data.frame(Patches=numCom-colMeans(apply(is.na(Abund),3,colSums)),Component_num=Components$Number_components,Component_size=Components$Component_size,Component_range=Components$Component_envt_range)%>%
        group_by(Patches) %>%
        summarise_each(funs(mean(.,na.rm=T)))
      
      mean.df<-summarise(group_by(Meta_dyn,Patches),Species_sorting=mean(Species_sorting,na.rm=T),Mass_effects=mean(Mass_effects,na.rm=T),Base_growth=mean(Base_growth,na.rm=T))
      
      SIH.df_temp<-cbind(SIH_data_means,Component_data_means[,-1],mean.df[,-1])
      SIH.df_temp$Dispersal<-disp
      SIH.df_temp$Scenario<-removeV[j]
      if(i==1 & j ==1){
        SIH.df<-SIH.df_temp
      } else {
        SIH.df<-rbind(SIH.df,SIH.df_temp)
      }
    }}
  return(SIH.df)
}

#run simulation function in parallel
Sim_data_parallel<-foreach(r = 37:100,.packages=c("igraph","dplyr","tidyr")) %dopar% SIH_frag()

stopCluster(cl)

save(Sim_data_parallel,file = "Savedata2.RData")

Sim_data<-do.call("rbind",Sim_data_parallel)

Sim_data_long<-gather(Sim_data,key = Response,value = Value,R_SR:Base_growth)

SIH_means<-Sim_data_long%>%
  group_by(Dispersal,Response,Scenario,Patches)%>%
  summarise(Mean=mean(Value,na.rm=T),Lower=quantile(Value,probs = c(0.025),na.rm=T),Upper=quantile(Value,probs=c(0.975),na.rm=T),LQ=quantile(Value,probs=0.25,na.rm=T),UQ=quantile(Value,probs=0.75,na.rm=T),Median=median(Value,na.rm=T))

save(SIH_means,file="Fragmentation.RData")

SIH_means$Scenario<-factor(SIH_means$Scenario,levels=c("Min betweenness", "Random", "Max betweenness"),ordered = T)

#Figures####
options(scipen=9)
#Figure 3####
Components.df<-filter(SIH_means,Response == "Component_num"|
                        Response == "Component_range"|
                        Response == "Component_size",
                      Dispersal==0.0005)

Components.df$Response<-replace(Components.df$Response,Components.df$Response=="Component_num", "Component number")
Components.df$Response<-replace(Components.df$Response,Components.df$Response=="Component_size", "Component size")
Components.df$Response<-replace(Components.df$Response,Components.df$Response=="Component_range", "Component range")

ggplot(Components.df,aes(x=Patches,y=Mean, color=Scenario, fill=Scenario))+
  geom_ribbon(aes(ymin = LQ, ymax = UQ),alpha=0.3, color=NA)+
  geom_line(size=1.2)+
  facet_grid(Response~.,scale='free')+
  scale_color_manual(values = c("dodgerblue1","black","red"),name="")+
  scale_fill_manual(values = c("dodgerblue1","black","red"),name="")+
  xlim(100,0)+
  theme_bw(base_size = 16)+
  ylab("Mean value")+
  theme(legend.position="top")+
  removeGrid()
ggsave("./Figures/3. Network components with fragmentation.pdf",width = 6,height = 8.5)

#Figure 4####
SIH_trad.df<-filter(SIH_means,Response == "R_SR"|
                      Response == "L_SR"|
                      Response == "L_Occ"|
                      Response == "L_Bmass"|
                      Response == "R_CV"|
                      Response == "L_CV")

SIH_trad.df$Response<-factor(SIH_trad.df$Response,levels=c("R_SR","L_SR","L_Occ","L_Bmass","R_CV","L_CV"),ordered=T)
SIH_trad.df$cleanNames<-factor(SIH_trad.df$Response)
levels(SIH_trad.df$cleanNames)<-c("Regional\nspecies\nrichness","Local\nspecies\nrichness","Local\noccupancy","Local\nbiomass","Regional\nbiomass\nvariability","Local\nbiomass\nvariability")

SIH_trad.df$Dispersal_text<-paste("Dispersal =",SIH_trad.df$Dispersal)

ggplot(SIH_trad.df,aes(x=Patches,y=Mean, color=Scenario, fill=Scenario))+
  geom_ribbon(aes(ymin = LQ, ymax = UQ),alpha=0.2)+
  geom_line(size=1.2)+
  facet_grid(cleanNames~Dispersal_text,scale='free_y')+
  scale_color_manual(values = c("dodgerblue1","black","red"),name="")+
  scale_fill_manual(values = c("dodgerblue1","black","red"),name="")+
  xlim(100,0)+
  theme_bw(base_size = 16)+
  theme(legend.position="top")+
  ylab("Mean value")+
  removeGrid()
ggsave("./Figures/4. Diversity and biomass with fragmentation.pdf",width = 11,height = 8.5)

#Figure 5####
BE_Fcurve.df<-ungroup(filter(SIH_trad.df,Response == "L_Bmass" | Response == "L_SR"))
BEF_curve.df<-spread(BE_Fcurve.df[,1:5],key = Response,value = Mean)


ggplot(BEF_curve.df,aes(x=L_SR,y=L_Bmass, color=Scenario, fill=Scenario))+
  #geom_ribbon(aes(ymin = SD_min, ymax = SD_max),alpha=0.2)+
  geom_path(size=1.2)+
  facet_grid(.~Dispersal)+
  scale_shape_manual(values = c(25,19,24),name="")+
  scale_color_manual(values = c("dodgerblue1","black","red"),name="")+
  scale_fill_manual(values = c("dodgerblue1","black","red"),name="")+
  theme_bw(base_size = 16)+
  theme(legend.position="top")+
  geom_point(data=filter(BEF_curve.df,Patches==70),aes(x=L_SR,y=L_Bmass, color=Scenario,shape=Scenario,fill=Scenario),size=3)+
  xlab("Local species richness")+
  ylab("Local biomass")+
  removeGrid()
ggsave("./Figures/5. BEF curves.pdf", width = 11,height = 5)


#Figure 6####
Dynamics<-filter(SIH_means,Response == "Base_growth"|
                   Response == "Species_sorting"|
                   Response == "Mass_effects")

Dynamics$Response<-factor(Dynamics$Response,levels=c("Base_growth","Species_sorting","Mass_effects"),ordered = T)
levels(Dynamics$Response)<-c("Base growth","Species sorting","Mass effects")
Dynamics$Dispersal<-paste("Dispersal = ",Dynamics$Dispersal,sep="")


ggplot(Dynamics,aes(x=Patches,y=Mean, color=Response, fill=Response))+
  geom_ribbon(aes(ymin = LQ, ymax = UQ),alpha=0.2,color=NA)+
  geom_line(size=1.2)+
  facet_grid(Scenario~Dispersal)+
  scale_color_manual(values = brewer.pal(3,"Dark2")[c(2,1,3)],name="")+
  scale_fill_manual(values = brewer.pal(3,"Dark2")[c(2,1,3)],name="")+
  xlim(100,0)+
  theme_bw(base_size = 16)+
  scale_y_continuous(breaks=seq(0,1,length=3))+
  theme(legend.position="top")+
  ylab("Proportion of biomass production")+
  removeGrid()
ggsave("./Figures/6. SIH dynamics with fragmentation.pdf",width = 11,height = 8.5)
