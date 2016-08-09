require(dplyr)
require(ggplot2)
require(tidyr)
require(ggExtra)
require(RColorBrewer)

options(scipen=9)
#Figure 1####
#plot network examples####
require(igraph)
source("./Functions/maketriangles.r")

numCom<-100
eAMP<-1
ePeriod<-40000


shapeV<-rep("circle",numCom)

add.vertex.shape("uptriangle",plot=uptriangle)
add.vertex.shape("downtriangle",plot=downtriangle)

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
  
  
  if(components(graph)$no == 1){success<-T}}

envt.v<-0.5*eAMP*(sin((2*pi/ePeriod)*1+1+(landscape$y)*2*pi/1000)+1)


n.removed<-35

ColV2<-c(1,"red","dodgerblue")

pdf(file="./Figures/1. Network figure.pdf",width=7,height=7)
par(mfrow=c(2,2), mar=c(3,3,3,3),pty='s',oma=c(1,1,1,1))

ColV<-heat.colors(100)[1+(envt.v*99)]

#ColV[order(betweenness(graph),decreasing=T)[1]]<-ColV2[2]
#ColV[order(betweenness(graph),decreasing=F)[1]]<-ColV2[3]

shapeVr<-shapeV
shapeVr[order(betweenness(graph),decreasing=T)[1]]<-"uptriangle"
shapeVr[order(betweenness(graph),decreasing=F)[1]]<-"downtriangle"
plot.igraph(graph,layout=as.matrix(landscape), vertex.color=ColV, vertex.size=5000,vertex.label=NA, rescale=F, ylim=c(0,1000),xlim=c(0,1000),vertex.shape=shapeVr)
title("Intact\nnetwork",line = 1)

weightedgraph<-graph
landscape2<-landscape
for(d in 1:n.removed){
  patch.delete<-order(betweenness(weightedgraph),decreasing=F)[1]
  landscape2<-landscape2[-patch.delete,]
  weightedgraph<-delete.vertices(weightedgraph,patch.delete)
}
GrayColVr<-ColV[as.numeric(V(weightedgraph)$name)]
#GrayColVr[order(betweenness(weightedgraph),decreasing=F)[1]]<-"dodgerblue"
shapeVr<-rep("circle",length(V(weightedgraph)$name))
shapeVr[order(betweenness(weightedgraph),decreasing=F)[1]]<-"downtriangle"
plot.igraph(weightedgraph,layout=as.matrix(landscape2), vertex.color=GrayColVr, vertex.size=5000,vertex.label=NA, rescale=F, ylim=c(0,1000),xlim=c(0,1000),vertex.shape=shapeVr)

title("Remove\nmin betweenness",line = 1)

weightedgraph<-graph
landscape2<-landscape
for(d in 1:n.removed){
  patch.delete<-sample(vcount(weightedgraph),1)
  landscape2<-landscape2[-patch.delete,]
  weightedgraph<-delete.vertices(weightedgraph,patch.delete)
}

GrayColVr<-ColV[as.numeric(V(weightedgraph)$name)]
#GrayColVr[order(betweenness(weightedgraph),decreasing=F)[1]]<-"dodgerblue"
shapeVr<-rep("circle",length(V(weightedgraph)$name))
plot.igraph(weightedgraph,layout=as.matrix(landscape2), vertex.color=GrayColVr, vertex.size=5000,vertex.label=NA, rescale=F, ylim=c(0,1000),xlim=c(0,1000),vertex.shape=shapeVr)

title("Remove\nmin random",line = 1)

weightedgraph<-graph
landscape2<-landscape
for(d in 1:n.removed){
  patch.delete<-order(betweenness(weightedgraph),decreasing=T)[1]
  landscape2<-landscape2[-patch.delete,]
  weightedgraph<-delete.vertices(weightedgraph,patch.delete)
}
GrayColVr<-ColV[as.numeric(V(weightedgraph)$name)]
#GrayColVr[order(betweenness(weightedgraph),decreasing=F)[1]]<-"dodgerblue"
shapeVr<-rep("circle",length(V(weightedgraph)$name))
shapeVr[order(betweenness(weightedgraph),decreasing=T)[1]]<-"uptriangle"
plot.igraph(weightedgraph,layout=as.matrix(landscape2), vertex.color=GrayColVr, vertex.size=5000,vertex.label=NA, rescale=F, ylim=c(0,1000),xlim=c(0,1000),vertex.shape=shapeVr)

title("Remove\nmax betweenness",line = 1)
dev.off()

#Figure 2####
load("Meta_dynamics.RData")

ggplot(Meta.dyn.df,aes(x=Dispersal,y=Mean,color=Response, fill=Response))+
  geom_ribbon(aes(ymin = Lower, ymax = Upper),alpha=0.3,color=NA)+
  geom_line(size=1.5)+
  theme_bw(base_size = 15)+
  scale_color_manual(values = brewer.pal(3,"Dark2")[c(2,1,3)],name="",labels=c("Base growth","Species sorting","Mass effects"))+
  scale_fill_manual(values = brewer.pal(3,"Dark2")[c(2,1,3)],name="",labels=c("Base growth","Species sorting","Mass effects"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_x_log10(breaks=c(0.0001,0.001,0.01,0.1,1),labels=c("0.0001","0.001","0.01","0.1","1"))+
  scale_y_continuous(breaks=seq(0,1,length=5))+
  ylab("Propotion of biomass production")+
  theme(legend.justification=c(1,0),legend.position=c(1,0.5))+
  geom_vline(xintercept = c(0.0005,0.005,0.015), linetype=2)
ggsave("./Figures/2. Metacommunity dynamics.pdf",width = 8.5,height = 8.5)


#Figure 3####
load("Fragmentation.RData")

Components.df<-SIH_means%>%
  filter(Response=="Component_num" | Response=="Component_range" | Response=="Component_size" | Response == "Metacom_range"& Dispersal == 0.01)

Components.df$Response<-replace(Components.df$Response,Components.df$Response=="Component_num", "Component number")
Components.df$Response<-replace(Components.df$Response,Components.df$Response=="Component_size", "Component size")
Components.df$Response<-replace(Components.df$Response,Components.df$Response=="Component_range", "Component range")
Components.df$Response<-replace(Components.df$Response,Components.df$Response=="Metacom_range", "Network range")


ggplot(Components.df,aes(x=Patches,y=Mean, color=Scenario, fill=Scenario))+
  geom_ribbon(aes(ymin = Lower, ymax = Upper),alpha=0.3, color=NA)+
  geom_line(size=1.2)+
  facet_wrap(~Response,scale='free')+
  scale_color_manual(values = c("dodgerblue1","black","red"),name="")+
  scale_fill_manual(values = c("dodgerblue1","black","red"),name="")+
  xlim(100,0)+
  theme_bw(base_size = 16)+
  ylab("Mean value")+
  theme(legend.position="top")+
  removeGrid()
ggsave("./Figures/3. Network components with fragmentation.pdf",width = 8.5,height = 8.5)

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
  geom_ribbon(aes(ymin = Lower, ymax = Upper),alpha=0.3,color=NA)+
  geom_line(size=1.2)+
  facet_grid(cleanNames~Dispersal_text,scale='free_y')+
  scale_color_manual(values = c("dodgerblue1","black","red"),name="")+
  scale_fill_manual(values = c("dodgerblue1","black","red"),name="")+
  xlim(100,0)+
  theme_bw(base_size = 16)+
  theme(legend.position="top")+
  ylab("Mean value")+
  removeGrid()
ggsave("./Figures/4. Diversity and biomass with fragmentation - 95.pdf",width = 11,height = 8.5)

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
  geom_ribbon(aes(ymin = Lower, ymax = Upper),alpha=0.3,color=NA)+
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
ggsave("./Figures/6. SIH dynamics with fragmentation 95.pdf",width = 11,height = 8.5