source("./Functions/rewire.R")
source("./Functions/create_random_net.r")
source("./Functions/addweights.r")
require(igraph)
require(dplyr)
require(ggplot2)
require(tidyr)

SIH_function<-function(dispersal=0.005,species=9,patches=30,rand=50){
  #Constants####
  #N<- matrix(10,ncol=species,nrow=patches) # Community x Species abundance matrix
  N<-matrix(10,patches,species)
  R<-rep(10*(species/10),patches) #Initial resources
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
  numEdgesRewired<-rand/100*(patches*2) 
  numLinks<-patches*2
  success<-FALSE
  while(!success){unweightedgraph<- if(rand==100) create_random_net(patches, numLinks) else rewire(patches,numLinks,numEdgesRewired)
  success<-length(V(unweightedgraph))==30}
  weightedgraph<-addweights(unweightedgraph,numLinks,patches)
  
  #dispersal conditions####
  d<-shortest.paths(weightedgraph, mode="all", weights=NULL, algorithm="automatic")
  d_exp<-exp(-1*d) - diag(nrow(d))  #dispersal kernel function of the d matrix
  dispersal_matrix <- apply(d_exp, 1, function(x) x/sum(x)) #divides the d_exp matrix by the column sums to make it a conservative dispersal matrix
  
  calc.immigration <- function(N,a,dispersal_matrix) dispersal_matrix%*%N*rep(a,each=patches)
  
  Prod<-matrix(NA,species*patches,40000)
  Abund<-Prod
  
  Meta_dyn<-data.frame(Species_sorting=rep(NA,40000),Mass_effects=NA,Base_growth=NA)
  Species_data<-array(NA,dim=c(40000,species,2),dimnames = list(1:40000,1:species,c("Abundance","Occupancy")))
  
  for(TS in 1:Tmax){
    envt.v<-0.5*eAMP*(sin((2*pi/ePeriod)*TS+1+(1:patches)*2*pi/patches)+1)
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
  
  Prod<-array(t(Prod),dim=c(40000,species,patches))
  Prod<-Prod[seq(1,40000,100),,]
  
  Abund<-array(t(Abund),dim=c(40000,species,patches))
  Abund<-Abund[seq(1,40000,100),,]
  matplot(Abund[,,1], type ='l', lty=1)
  return(list(Prod=Prod,Abund=Abund,Meta_dyn=Meta_dyn,Spec_data=apply(Species_data,3,colMeans)))
}

dispV<-c(0.0001,0.0005,0.001,0.005,0.01,0.05,0.1,0.5,1)
vect<-c(0.0001,0.00025,0.0005,0.00075)
dispV<-c(vect,vect*10,vect*100,vect*1000,1)
#dispV<-dispV[-12]
dispV<-c(0.0001,0.0003,0.0004,0.0005,0.001,0.005,0.01,0.015,0.02,0.05,0.1,0.5,1)
dispV<-c(0.0005,0.005,0.015)

species<-9
reps<-20
SIH_results<-data.frame(SR=NA,Diversity=NA,Biomass=NA,Biomass_CV=NA,Species_sorting=NA,Mass_effects=NA,Base_growth=NA,Dispersal=dispV,Scale=rep(c("Local","Regional"),each=length(dispV)),Rep=rep(1:reps,each=length(dispV)*2))
Species_data<-data.frame(Abundance=NA,Occupancy=NA,Species=1:species,Dispersal=rep(dispV,each=species),Rep=rep(1:reps,each=length(dispV)*species))
for(r in 1:reps){
  print(r)
  SIH_data<-sapply(dispV,SIH_function,species=species,patches=30,rand=50)
  require(vegan)
  for(i in 1:length(dispV)){
    SIH_results$SR[SIH_results$Dispersal==dispV[i] & SIH_results$Scale == "Local" & SIH_results$Rep==r]<-mean(apply(SIH_data[["Abund",i]]>0,3,rowSums))
    SIH_results$SR[SIH_results$Dispersal==dispV[i] & SIH_results$Scale == "Regional" & SIH_results$Rep==r]<-mean(rowSums(apply(SIH_data[["Abund",i]],2,rowSums)>0))
    SIH_results$Diversity[SIH_results$Dispersal==dispV[i] & SIH_results$Scale == "Local" & SIH_results$Rep==r]<-mean(exp(apply(SIH_data[["Abund",i]],3,diversity)))
    SIH_results$Diversity[SIH_results$Dispersal==dispV[i] & SIH_results$Scale == "Regional" & SIH_results$Rep==r]<-mean(exp(diversity(apply(SIH_data[["Abund",i]],2,rowSums))))
    SIH_results$Biomass[SIH_results$Dispersal==dispV[i] & SIH_results$Scale == "Local" & SIH_results$Rep==r]<-mean(apply(SIH_data[["Abund",i]],3,rowSums))
    SIH_results$Biomass[SIH_results$Dispersal==dispV[i] & SIH_results$Scale == "Regional" & SIH_results$Rep==r]<-mean(rowSums(apply(SIH_data[["Abund",i]],2,rowSums)))
    SIH_results$Biomass_CV[SIH_results$Dispersal==dispV[i] & SIH_results$Scale == "Local" & SIH_results$Rep==r]<-sd(apply(SIH_data[["Abund",i]],3,rowSums))/mean(apply(SIH_data[["Abund",i]],3,rowSums))
    SIH_results$Biomass_CV[SIH_results$Dispersal==dispV[i] & SIH_results$Scale == "Regional" & SIH_results$Rep==r]<-sd(rowSums(apply(SIH_data[["Abund",i]],2,rowSums)))/mean(rowSums(apply(SIH_data[["Abund",i]],2,rowSums)))
    SIH_results[SIH_results$Dispersal==dispV[i] & SIH_results$Rep==r,5:7 ]<-rep(colMeans(SIH_data[["Meta_dyn",i]],na.rm=T),each=2)
    Species_data[Species_data$Dispersal==dispV[i],1:2]<-SIH_data[["Spec_data",i]]
  }}

Meta_dynamics<-gather(Meta_dyn.df,key = Dynamic,value=Proportion_of_production,Species_sorting:Base_growth)
Meta_dynamics$Dynamic<-factor(Meta_dynamics$Dynamic,levels = c("Base_growth","Species_sorting","Mass_effects"),ordered = T)
Meta_dynamics$Dynamic_clean<-factor(Meta_dynamics$Dynamic,levels = c("Base growth","Species sorting","Mass effects"),ordered=T)
Meta_dynamics$Dynamic_clean<-factor(Meta_dynamics$Dynamic,labels = c("Base growth","Species sorting","Mass effects"),ordered=T)


Meta_dynamics_means<-Meta_dynamics%>%
  group_by(Dispersal,Dynamic_clean)%>%
  summarise(Mean=mean(Proportion_of_production), SD=sd(Proportion_of_production))%>%
  mutate(Max_sd=Mean+SD,Min_sd=Mean-SD)%>%
  mutate(Max_sd=replace(Max_sd,Max_sd>1,1),Min_sd=replace(Min_sd,Min_sd<0,0))

require(ggplot2)
require(dplyr)
require(RColorBrewer)

pdf("./Figures/2. Metacommunity dynamics.pdf",width = 8,height = 6)
ggplot(Meta_dynamics_means,aes(x=Dispersal,y=Mean,group=Dynamic_clean,color=Dynamic_clean, fill=Dynamic_clean))+
  geom_ribbon(aes(ymin = Min_sd, ymax = Max_sd),alpha=0.3)+
  geom_line(size=1.5)+
  theme_bw(base_size = 15)+
  scale_color_manual(values = brewer.pal(3,"Set1")[c(1,3,2)],name="")+
  scale_fill_manual(values = brewer.pal(3,"Set1")[c(1,3,2)],name="")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_x_log10(breaks=c(0.0001,0.001,0.01,0.1,1),labels=c("0.0001","0.001","0.01","0.1","1"))+
  scale_y_continuous(breaks=seq(0,1,length=5))+
  ylab("Propotion of biomass production")+
  theme(legend.justification=c(1,0),legend.position=c(1,0.5))+
  geom_vline(xintercept = c(0.0005,0.005,0.015), linetype=2)
dev.off()

ggplot(SIH_results,aes(x=Dispersal,y=SR,group=Scale,color=Scale))+
  geom_line(size=1.5)+
  scale_y_continuous(breaks=0:9)+
  theme_bw(base_size = 15)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_x_log10(breaks=c(0.0001,0.001,0.01,0.1,1),labels=c("0.0001","0.001","0.01","0.1","1"))

ggplot(SIH_results,aes(x=Dispersal,y=Diversity,group=Scale,color=Scale))+
  geom_line(size=1.5)+
  scale_y_continuous(breaks=0:9)+
  theme_bw(base_size = 15)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_x_log10(breaks=c(0.0001,0.001,0.01,0.1,1),labels=c("0.0001","0.001","0.01","0.1","1"))

ggplot(SIH_results,aes(x=Dispersal,y=Biomass,group=Scale,color=Scale))+
  geom_line(size=1.5)+
  scale_y_continuous()+
  theme_bw(base_size = 15)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_x_log10(breaks=c(0.0001,0.001,0.01,0.1,1),labels=c("0.0001","0.001","0.01","0.1","1"))

ggplot(SIH_results,aes(x=Dispersal,y=Biomass_CV,group=Scale,color=Scale))+
  geom_line(size=1.5)+
  scale_y_continuous()+
  theme_bw(base_size = 15)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_x_log10(breaks=c(0.0001,0.001,0.01,0.1,1),labels=c("0.0001","0.001","0.01","0.1","1"))

pdf(file = "2.Diversity - ME.pdf",width = 8,height = 8)
ggplot(SIH_results,aes(y=Diversity,x=Mass_effects,group=Scale))+
  geom_path(size=1.5,aes(color=Dispersal))+
  scale_color_gradient(low = "blue",high="red",trans = "log10")+
  facet_grid(Scale~.)+
  scale_y_continuous()+
  theme_bw(base_size = 15)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()

pdf(file = "3.Diversity - SS.pdf",width = 8,height = 8)
ggplot(SIH_results,aes(y=Diversity,x=Species_sorting,group=Scale))+
  geom_path(size=1.5,aes(color=Dispersal))+
  scale_color_gradient(low = "blue",high="red",trans = "log10")+
  facet_grid(Scale~.)+
  scale_y_continuous()+
  theme_bw(base_size = 15)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()

pdf(file = "5.Biomass - SS.pdf",width = 8,height = 5)
ggplot(SIH_results[SIH_results$Scale=="Local",],aes(y=Biomass,x=Species_sorting))+
  geom_path(size=1.5,aes(color=Dispersal))+
  scale_color_gradient(low = "blue",high="red",trans = "log10")+
  scale_y_continuous()+
  theme_bw(base_size = 15)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()

pdf(file = "4.Biomass - ME.pdf",width = 8,height = 5)
ggplot(SIH_results[SIH_results$Scale=="Local",],aes(y=Biomass,x=Mass_effects))+
  geom_path(size=1.5,aes(color=Dispersal))+
  scale_color_gradient(low = "blue",high="red",trans = "log10")+
  scale_y_continuous()+
  theme_bw(base_size = 15)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()

pdf(file = "7.Biomass CV - SS.pdf",width = 8,height = 8)
ggplot(SIH_results,aes(y=Biomass_CV,x=Species_sorting,group=Scale))+
  geom_path(size=1.5,aes(color=Dispersal))+
  scale_color_gradient(low = "blue",high="red",trans = "log10")+
  facet_grid(Scale~.)+
  scale_y_continuous()+
  theme_bw(base_size = 15)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()

pdf(file = "6.Biomass CV - ME.pdf",width = 8,height = 8)
ggplot(SIH_results,aes(y=Biomass_CV,x=Mass_effects,group=Scale))+
  geom_path(size=1.5,aes(color=Dispersal))+
  scale_color_gradient(low = "blue",high="red",trans = "log10")+
  facet_grid(Scale~.)+
  scale_y_continuous()+
  theme_bw(base_size = 15)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()

require(scatterplot3d)
require(dplyr)
colF<-colorRampPalette(colors = c("blue","red"))
Local_data<-subset(SIH_results,Scale=="Local")
pdf(file = "8.Diversity 3D.pdf",width = 8,height = 8)
scatterplot3d(y=Local_data$Species_sorting,x=Local_data$Mass_effects,z=Local_data$Diversity,color = colF(length(dispV)),pch=19, type='h', xlab="Mass effects",ylab="Species sorting",zlab="Diversity")
dev.off()

pdf(file = "9.Biomass 3D.pdf",width = 8,height = 8)
scatterplot3d(y=Local_data$Species_sorting,x=Local_data$Mass_effects,z=Local_data$Biomass,color = colF(length(dispV)),pch=19, type='h', xlab="Mass effects",ylab="Species sorting",zlab="Biomass")
dev.off()

pdf(file = "10.Biomass 3D.pdf",width = 8,height = 8)
scatterplot3d(y=Local_data$Species_sorting,x=Local_data$Mass_effects,z=Local_data$Biomass_CV,color = colF(length(dispV)),pch=19, type='h', xlab="Mass effects",ylab="Species sorting",zlab="Biomass CV")
dev.off()

#species plots####
pdf(file = "11.Abundance by dispersal.pdf",width = 8,height = 8)
ggplot(subset(Species_data,Species<6),aes(x=Dispersal,y=Abundance,group=Species,color=factor(Species),linetype=factor(Species)))+
  geom_point(size=1.5)+
  geom_line(size=1.5)+
  theme_bw(base_size = 15)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_x_log10(breaks=c(0.0001,0.001,0.01,0.1,1),labels=c("0.0001","0.001","0.01","0.1","1"))
dev.off()

pdf(file = "12.Occupancy by dispersal.pdf",width = 8,height = 8)
ggplot(subset(Species_data,Species<6),aes(x=Dispersal,y=Occupancy,group=Species,color=factor(Species),linetype=factor(Species)))+
  geom_point(size=1.5)+
  geom_line(size=1.5)+
  theme_bw(base_size = 15)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_x_log10(breaks=c(0.0001,0.001,0.01,0.1,1),labels=c("0.0001","0.001","0.01","0.1","1"))
dev.off()

pdf(file = "13.Abundance by Species sorting.pdf",width = 8,height = 8)
ggplot(subset(Meta_species_data,Species<6),aes(x=Species_sorting,y=Abundance,group=Species,color=factor(Species),linetype=factor(Species)))+
  geom_point(size=1.5)+
  geom_line(size=1.5)+
  theme_bw(base_size = 15)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()

pdf(file = "14.Abundance by mass effects.pdf",width = 8,height = 8)
ggplot(subset(Meta_species_data,Species<6),aes(x=Mass_effects,y=Abundance,group=Species,color=factor(Species),linetype=factor(Species)))+
  geom_point(size=1.5)+
  geom_line(size=1.5)+
  theme_bw(base_size = 15)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()

pdf(file = "15.Abundance by base growth.pdf",width = 8,height = 8)
ggplot(subset(Meta_species_data,Species<6),aes(x=Base_growth,y=Abundance,group=Species,color=factor(Species),linetype=factor(Species)))+
  geom_point(size=1.5)+
  geom_line(size=1.5)+
  theme_bw(base_size = 15)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()

pdf(file = "16.Occupancy by species sorting.pdf",width = 8,height = 8)
ggplot(subset(Meta_species_data,Species<6),aes(x=Species_sorting,y=Occupancy,group=Species,color=factor(Species),linetype=factor(Species)))+
  geom_point(size=1.5)+
  geom_line(size=1.5)+
  theme_bw(base_size = 15)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()

pdf(file = "17.Occupancy by mass effects.pdf",width = 8,height = 8)
ggplot(subset(Meta_species_data,Species<6),aes(x=Mass_effects,y=Occupancy,group=Species,color=factor(Species),linetype=factor(Species)))+
  geom_point(size=1.5)+
  geom_line(size=1.5)+
  theme_bw(base_size = 15)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()

pdf(file = "18.Occupancy by base growth.pdf",width = 8,height = 8)
ggplot(subset(Meta_species_data,Species<6),aes(x=Base_growth,y=Occupancy,group=Species,color=factor(Species),linetype=factor(Species)))+
  geom_point(size=1.5)+
  geom_line(size=1.5)+
  theme_bw(base_size = 15)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()
