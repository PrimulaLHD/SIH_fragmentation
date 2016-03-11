source("./Functions/rewire.R")
source("./Functions/create_random_net.r")
source("./Functions/addweights.r")
require(igraph)
require(dplyr)
require(ggplot2)
require(tidyr)

options(scipen=9)

reps<-10
print.plots<-F # set this to true if you want to see the network as the sim runs - it makes it slower

nSpecies<-9
numCom<-30
randV<-50#seq(10,90,by=20)
dispV<- c(0.0005,0.005,0.015)
dd<-1
numLinks<-numCom*2


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

Meta_dyn_reps<-data.frame(Rep=rep(1:reps,each=(numCom-0)*3),Dispersal=rep(dispV,each=reps*(numCom-0)*3),Patch_remove=rep(factor(removeV,levels = c("Min betweenness","Random","Max betweenness"),ordered = T),each=length(dispV)*reps*(numCom-0)*3),Patches=NA,Dynamic=rep(factor(c("Species sorting", "Mass effects", "Base growth"),levels = c("Base growth","Species sorting","Mass effects")),each=numCom-0),Proportion=NA)
SIH_data_reps<-data.frame(Rep=rep(1:reps,each=(numCom-0)),Dispersal=rep(dispV,each=reps*(numCom-0)),Patch_remove=rep(factor(removeV,levels = c("Min betweenness","Random","Max betweenness"),ordered = T),each=length(dispV)*reps*(numCom-0)),Patches=NA,Regional_SR=NA,Local_SR=NA,Biomass=NA,Occupancy=NA,Regional_CV=NA,Local_CV=NA)
Component_data_reps<-data.frame(Rep=rep(1:reps,each=(numCom-0)),Dispersal=rep(dispV,each=reps*(numCom-0)),Patch_remove=rep(factor(removeV,levels = c("Min betweenness","Random","Max betweenness"),ordered = T),each=length(dispV)*reps*(numCom-0)),Patches=NA,Component_num=NA,Component_size=NA, Component_range=NA)



#initialize community network use rewire for lattice or small world - use random for random
pb <- txtProgressBar(min = 0, max = reps, style = 3)
for(r in 1:reps){
  for(i in 1:length(dispV)){
    disp<-dispV[i]
    rand<-randV[1]
    numEdgesRewired<-rand/100*(numCom*2) 
    success<-FALSE
    while(!success){unweightedgraph<- if(rand==100) create_random_net(numCom, numLinks) else rewire(numCom,numLinks,numEdgesRewired)
    success<-length(V(unweightedgraph))==30}
    for(j in 1:3){
      weightedgraph<-addweights(unweightedgraph,numLinks,numCom)
      holdgraph<-weightedgraph
      if(print.plots==T){plot(holdgraph, ylim=c(-1,1),xlim=c(-1,1))}
      d<-shortest.paths(weightedgraph, mode="all", weights=NULL, algorithm="automatic")
      d_exp<-exp(-dd*d) - diag(nrow(d))  #dispersal kernel function of the d matrix
      dispersal_matrix <- apply(d_exp, 1, function(x) x/sum(x)) #divides the d_exp matrix by the column sums to make it a conservative dispersal matrix
      
      #vectors####
      eOptimum<-1-seq(0,eAMP, by=eAMP/(nSpecies-1)) #species environmental optima
      
      calc.immigration <- function(N,a,dispersal_matrix) dispersal_matrix%*%N*rep(a,each=length(R))
      
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
        envt.v<-0.5*eAMP*(sin((2*pi/ePeriod)*TS+1+(1:numCom)*2*pi/numCom)+1)
        if(is.null(rownames(dispersal_matrix))){
          envt.v<-envt.v[as.numeric(names(dispersal_matrix))]
          consume <- 0.1*(1.5-abs(sapply(eOptimum,'-',envt.v)))
          Nt <- N*(1+DT*(eff*R*consume - disp - mort)) + DT*Immigrants #abundance step
          Rt <- DT*rInput+R*(1-DT*(rLoss + sum(consume*N))) #resource step    
          
          Nt0 <- N0*(1+DT*(eff*R0*consume - mort))
          Rt0 <- DT*rInput+R0*(1-DT*(rLoss + sum(consume*N0)))} else { #resource step  
            envt.v<-envt.v[as.numeric(rownames(dispersal_matrix))]
            consume <- 0.1*(1.5-abs(sapply(eOptimum,'-',envt.v)))
            Nt <- N*(1+DT*(eff*R*consume - disp - mort)) + DT*Immigrants #abundance step
            Rt <- DT*rInput+R*(1-DT*(rLoss + rowSums(consume*N))) #resource step    
            
            Nt0 <- N0*(1+DT*(eff*R0*consume - mort))
            Rt0 <- DT*rInput+R0*(1-DT*(rLoss + rowSums(consume*N0)))
          }
        
        if(sum(TS==sampleV)==1){
          sample_id<-which(TS==sampleV)
          Components$Number_components[sample_id]<-components(weightedgraph)$no
          Components$Component_size[sample_id]<-mean(components(weightedgraph)$csize)
          members<-components(weightedgraph)$membership
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
          if(j==1){btw<-betweenness(weightedgraph)
          if(sum(btw==0)){
            patch.delete<-order(degree(weightedgraph),decreasing = T)[1]
          } else{patch.delete<-order(btw,decreasing = T)[1] }
          } else{
            if(j==2){patch.delete<-order(betweenness(weightedgraph),decreasing=F)[1]} else{
              patch.delete<-sample(nrow(N),1)}}    
          weightedgraph<-delete.vertices(weightedgraph,patch.delete)
          d<-shortest.paths(weightedgraph, mode="all", weights=NULL, algorithm="automatic")
          d_exp<-exp(-dd*d) - diag(nrow(d))  #dispersal kernel function of the d matrix
          dispersal_matrix <- apply(d_exp, 1, function(x) x/sum(x)) #divides the d_exp matrix by the column sums to make it a conservative dispersal matrix
          dispersal_matrix[is.nan(dispersal_matrix)]<-0
          if(print.plots==T){
            if(length(V(weightedgraph))>1){plot(weightedgraph,layout=layout.circle(holdgraph)[as.numeric(colnames(dispersal_matrix)),],ylim=c(-1,1),xlim=c(-1,1))} else{plot(weightedgraph)}}
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
      
      CVdf<-cbind(L_Bmass_sep,data.frame(R_Bmass=R_Bmass,Patches=rep(30:1,each=drop_length/2000))) %>%
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

      SIH_data_reps[SIH_data_reps$Rep==r & SIH_data_reps$Dispersal==dispV[i] & SIH_data_reps$Patch_remove==removeV[j],-c(1:3)]<-SIH_data_means
      Component_data_reps[SIH_data_reps$Rep==r & SIH_data_reps$Dispersal==dispV[i] & SIH_data_reps$Patch_remove==removeV[j],-c(1:3)]<-Component_data_means
      
      mean.df<-summarise(group_by(Meta_dyn,Patches),Species_sorting=mean(Species_sorting,na.rm=T),Mass_effects=mean(Mass_effects,na.rm=T),Base_growth=mean(Base_growth,na.rm=T))
      Meta.dyn.long<-gather(mean.df,key = Dynamic,value=Proportion,-Patches)
      
      Meta_dyn_reps[Meta_dyn_reps$Rep==r & Meta_dyn_reps$Dispersal==dispV[i] & Meta_dyn_reps$Patch_remove==removeV[j],-c(1:3,5)]<-Meta.dyn.long[,-2]
    }}
  Sys.sleep(0.1)
  setTxtProgressBar(pb, r)
}

save(Meta_dyn_reps,Component_data_reps,SIH_data_reps,file="Fragmentation.RData")



hold<-Meta_dyn_reps%>%
  group_by(Dispersal,Patch_remove,Patches,Dynamic)%>%
  summarise(Mean=mean(Proportion,na.rm=T),Proportion_sd=sd(Proportion,na.rm = T))%>%
  mutate(Max_sd=Mean+Proportion_sd,Min_sd=Mean-Proportion_sd)%>%
  mutate(Max_sd=replace(Max_sd,Max_sd>1,1),Min_sd=replace(Min_sd,Min_sd<0,0))
  

hold$Dispersal_text<-paste("Dispersal =",hold$Dispersal)


require(RColorBrewer)
pdf("./Figures/4. SIH dynamics with fragmentation.pdf",width = 11,height = 8.5)
ggplot(hold,aes(x=Patches,y=Mean, group=interaction(Dynamic,Patch_remove), color=Dynamic, fill=Dynamic))+
  geom_ribbon(aes(ymin = Min_sd, ymax = Max_sd),alpha=0.2)+
  geom_line(size=1.2)+
  facet_grid(Patch_remove~Dispersal_text)+
  scale_color_manual(values = brewer.pal(3,"Set1")[c(1,3,2)],name="")+
  scale_fill_manual(values = brewer.pal(3,"Set1")[c(1,3,2)],name="")+
  xlim(30,0)+
  theme_bw(base_size = 16)+
  scale_y_continuous(breaks=seq(0,1,length=3))+
  theme(legend.position="top")+
  ylab("Proportion of biomass production")
dev.off()



SIH_long<-gather(SIH_data_reps,key = SIH_attribute,value = Value,Regional_SR:Local_CV)

pdf("./Figures/6. BEF curves.pdf", width = 11,height = 4)
ggplot(summarise(group_by(SIH_data_reps,Patch_remove,Dispersal,Patches),Local_SR=mean(Local_SR),Biomass=mean(Biomass)),aes(x=Local_SR,y=Biomass, color=Patch_remove))+
  geom_path(size=1.2)+
  facet_grid(.~Dispersal)+
  scale_color_manual(values = c("dodgerblue1","black","red"),name="")+
  theme(legend.position="top")+
  theme_bw(base_size = 16)
dev.off()


SIH_means<-SIH_long %>%
  select(Rep,Dispersal,Patch_remove,Patches,SIH_attribute,Value) %>%
  group_by(Dispersal,Patch_remove,Patches,SIH_attribute) %>%
  summarise(Mean=mean(Value,na.rm=T),value_sd=sd(Value,na.rm=T),Max=max(Value,na.rm=T),Min=min(Value,na.rm=T))%>%
  group_by(Dispersal,Patch_remove,Patches,SIH_attribute) %>%
  mutate(Max_sd=Mean+value_sd,Min_sd=Mean-value_sd)%>%
  mutate(Max_sd=replace(Max_sd,Max_sd>Max,Max),Min_sd=replace(Min_sd,Min_sd<Min,Min))

SIH_means$SIH_attribute<-factor(SIH_means$SIH_attribute,levels=c("Regional_SR","Local_SR","Occupancy","Biomass","Regional_CV","Local_CV"),ordered=T)
SIH_means$cleanNames<-factor(SIH_means$SIH_attribute,levels=c("Regional\nspecies\nrichness","Local\nspecies\nrichness","Local\noccupancy","Local\nbiomass","Regional\nbiomass\nvariability","Local\nbiomass\nvariability"),ordered=T)
SIH_means$cleanNames<-factor(SIH_means$SIH_attribute,labels=c("Regional\nspecies\nrichness","Local\nspecies\nrichness","Local\noccupancy","Local\nbiomass","Regional\nbiomass\nvariability","Local\nbiomass\nvariability"),ordered=T)

SIH_means$Dispersal_text<-paste("Dispersal =",SIH_means$Dispersal)

mt<-ggplot(SIH_means,aes(x=Patches,y=Mean, group=interaction(Dispersal_text,Patch_remove), color=Patch_remove, fill=Patch_remove))+
  geom_ribbon(aes(ymin = Min_sd, ymax = Max_sd),alpha=0.2)+
  geom_line(size=1.2)+
  facet_grid(cleanNames~Dispersal_text,scale='free_y')+
  scale_color_manual(values = c("dodgerblue1","black","red"),name="")+
  scale_fill_manual(values = c("dodgerblue1","black","red"),name="")+
  xlim(30,0)+
  theme_bw(base_size = 16)+
  theme(legend.position="top")+
  ylab("Mean value")

p1 <- mt
g1 <- ggplotGrob(p1)

p2 <- mt + coord_cartesian(ylim = c(0,1))
g2 <- ggplotGrob(p2)

# ----------------------------------------------------------
# Replace the upper panels and upper axis of p1 with that of p2
# Tweak panels of second plot - the upper panels
library(gtable)
g1[["grobs"]][[15]] <- g2[["grobs"]][[15]]
g1[["grobs"]][[16]] <- g2[["grobs"]][[16]]
g1[["grobs"]][[21]] <- g2[["grobs"]][[21]]
g1[["grobs"]][[22]] <- g2[["grobs"]][[22]]
g1[["grobs"]][[27]] <- g2[["grobs"]][[27]]
g1[["grobs"]][[28]] <- g2[["grobs"]][[28]]

#Tweak axis
g1[["grobs"]][[9]] <- g2[["grobs"]][[9]]
g1[["grobs"]][[10]] <- g2[["grobs"]][[10]]


grid.newpage()

pdf("./Figures/5. Diversity and biomass with fragmentation.pdf",width = 11,height = 8.5)
grid.draw(g1)
dev.off()




Network_components<-gather(Component_data_reps,key = Component_attribute,value = Value,Component_num:Component_range)

Component_means<-Network_components %>%
  select(Rep,Patch_remove,Patches,Component_attribute,Value) %>%
  group_by(Patch_remove,Patches,Component_attribute) %>%
  summarise(Mean=mean(Value,na.rm=T),value_sd=sd(Value,na.rm=T))

Component_means$Min_sd<-Component_means$Mean-Component_means$value_sd
Component_means$Min_sd[Component_means$Min_sd<0]<-0
Component_means$Max_sd<-Component_means$Mean+Component_means$value_sd
Component_means$Max_sd[Component_means$Max_sd>30]<-30
Component_means$Max_sd[Component_means$Component_attribute=="Component_range" & Component_means$Max_sd>1]<-1

Component_means$cleanNames<-factor(Component_means$Component_attribute,levels=c("Component number","Component size","Component range"))
Component_means$cleanNames<-factor(Component_means$Component_attribute,labels=c("Component number","Component size","Component range"))


pdf("./Figures/3. Network components with fragmentation.pdf",width = 6,height = 8.5)
ggplot(Component_means,aes(x=Patches,y=Mean, group=Patch_remove, color=Patch_remove, fill=Patch_remove))+
  geom_ribbon(aes(ymin = Min_sd, ymax = Max_sd),alpha=0.3)+
  geom_line(size=1.2)+
  facet_grid(cleanNames~.,scale='free')+
  scale_color_manual(values = c("dodgerblue1","black","red"),name="")+
  scale_fill_manual(values = c("dodgerblue1","black","red"),name="")+
  xlim(30,0)+
  theme_bw(base_size = 16)+
  ylab("Mean value")+
  theme(legend.position="top")
dev.off()

