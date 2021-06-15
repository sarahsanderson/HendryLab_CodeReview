
#
##  PCA PLOT FUNCTION for lake means (creates a pca biplot that is nicer looking than )
# pca is a summary of a pca object (i.e., summary(pca))
# xlab and ylab are character strings to be used for axis labels (e.g., "PC1 - XX.X%")
plot.Lmean<-function(pca,xlab,ylab){
  ggplot() + 
    geom_segment(aes(x=0, y=0, xend=PC1, yend=PC2), arrow=arrow(length = unit(0.15,"cm")), 
                 data=pca$species, color="coral1", alpha=0.5) + theme_classic() + coord_fixed()+
    geom_text_repel(data=pca$species, aes(x=PC1,y=PC2,label = row.names(pca$species)),box.padding=.1,point.padding = .01,
                    colour="red4",fontface="bold", size =2.5, show.legend = F)+
    geom_text(data=pca$sites, aes(x=PC1,y=PC2,label = row.names(pca$sites)), size = 3, colour="steelblue4", show.legend = F)+
    xlab(xlab)+ylab(ylab)
}


#
##########
#         ECOLOGICAL DATA
##########
#
#



###### Depth profiles 
# load depth profiles
depth<-read.csv("depth_profiles.csv")[-(90:92),c(1:7)]
# rename cols
colnames(depth)<-c("Lake","Depth.ft","Depth.m","image.area","acres","hectares","scale.fac")
depth$Lake<-droplevels(depth$Lake) #Drops lakes w/o depth data
# the following fits each depth profile with a monotonic spline and interpolates
#   area at each .5 m depth increment
depth.by.lake<-by(depth,depth$Lake, FUN = function(x) unique(data.frame(
  area=splinefun(x=x$Depth.m,y=x$hectares, ties = mean,method="hyman")(
    c((seq(from =0,to= min(x$Depth.m),by=-.5)),min(x$Depth.m))),
  depth=c((seq(from =0,to= min(x$Depth.m),by=-.5)),min(x$Depth.m)),
  Lake=rep(x$Lake, length(seq(from =0,to= min(x$Depth.m),by=-.5))+1))))
depth.splines<-do.call(rbind,depth.by.lake)

# Uses splines calculated above to estimate proportion of lake shallower than 3 meters
litt.area<-c(by(depth,depth$Lake, FUN = function(x) litt.area=1-splinefun(
  x=x$Depth.m,y=x$hectares, ties = mean,method="hyman")(-3.0)/max(x$hectares)))

# Creates figure showing area:depth relationships using interpolated splines
depth.map.splines<-ggplot(data = depth.splines, aes(x=area,y=depth))+
  geom_rect(xmin=0,xmax=Inf,ymax=0,ymin=-Inf,fill = "honeydew2")+
  labs(x="Area (ha)", y= "Depth (m)")+
  geom_ribbon(data=depth.splines, aes(ymin=depth,ymax=0), fill="skyblue1")+
  geom_line(aes(y=-3),linetype="dashed", color="skyblue3")+
  facet_wrap(~Lake)+theme_classic()

# Creates figure showing area:depth relationships using only measured depths
depth.map.linear<-ggplot(data = depth, aes(x=hectares,y=Depth.m))+
  geom_rect(xmin=0,xmax=Inf,ymax=0,ymin=-Inf,fill = "honeydew2")+
  geom_ribbon(data=depth, aes(ymin=Depth.m,ymax=0), fill="skyblue1")+
  geom_line(aes(y=-3),linetype="dashed", color="skyblue3")+
  facet_wrap(~Lake)+theme_classic()

#
#
#
#
#
eco.dat<-read.csv("AK_ENV_zoo_invert_2018-2019_v2.csv")
eco.dat<-data.frame(eco.dat[1:14,c("Lake","LATITUDE","LONGITUDE","Region","Elevation..m.","area..ha.",
                                   "approx.max.z..m.","DOC..mg.L.","TP..ug.L.","TN..ug.L.","ChlA..ug.L.",
                                   "Sp.Cond..uS.cm.","pH","CALCIUM..mg.L.","Total.zoo.not.nauplii..invid.L.",
                                   "Total..Daphnia..indiv.L.","Total.macroinvert..m2","gamaridae",
                                   "chironomidae")],row.names = "Lake")
colnames(eco.dat)<-c("Lat.","Lon.","Region","Elev.","Area","Max_Depth","DOC","TP","TN",
                     "ChlA","Cond.","pH","Ca","T.Zoop","Daph.","T.macro",
                     "Gamm.","Chiro.")
eco.dat<-merge(eco.dat,litt.area,by ="row.names",all=T)%>%
  rename(Lake=Row.names,littoral=y)
rownames(eco.dat)<-eco.dat$Lake
(lake.pca<-rda(na.omit(eco.dat[,c("Area","Max_Depth","DOC","TP","TN",
                                  "ChlA","Cond.","pH","Ca")]),scale = T))%>%biplot()
(lake.pca.w.prey<-rda(na.omit(eco.dat[,c("Area","Max_Depth","DOC","TP","TN",
                                         "ChlA","Cond.","pH","Ca","Daph.",
                                         "Gamm.","littoral")]),scale = T))%>%biplot()
(lake.pca.trophic<-rda(na.omit(eco.dat[,c("Daph.","Gamm.","littoral")]),scale = T))%>%biplot()

# dimensionality and default biplots
biplot(lake.pca)
estimate.ED.eig(summary(lake.pca)$cont[["importance"]][2,])
biplot(lake.pca.w.prey)
estimate.ED.eig(summary(lake.pca.w.prey)$cont[["importance"]][2,])
biplot(lake.pca.trophic)
estimate.ED.eig(summary(lake.pca.trophic)$cont[["importance"]][2,])

Env.pca<-plot.Lmean(summary(lake.pca),"PC1 - 50.3%","PC2 - 17.4%")
Env.pca_w.prey<-plot.Lmean(summary(lake.pca.w.prey), "PC1 - 46.2%","PC2 - 20.9%")
Env.pca_trophic<-plot.Lmean(summary(lake.pca.trophic),"PC1 - 54.6%","PC2 - 33.5%")
Env.pca/Env.pca_w.prey.cond # 12.5 x 6 inches

(env.pca.vectors<-round(rbind("Eigenvalues"=lake.pca$CA$eig,"Singular Values"=sqrt(lake.pca$CA$eig),
                              "Proportion Explained"=lake.pca$CA$eig/sum(lake.pca$CA$eig),
                              lake.pca$CA$u,lake.pca$CA$v),3))%>%write.csv("env.pca.vectors.csv")
(env.pca.vectors_w.prey<-round(rbind("Eigenvalues"=lake.pca.w.prey$CA$eig,"Singular Values"=sqrt(lake.pca.w.prey$CA$eig),
                                     "Proportion Explained"=lake.pca.w.prey$CA$eig/sum(lake.pca.w.prey$CA$eig),
                                     lake.pca.w.prey$CA$u,lake.pca.w.prey$CA$v),3))%>%write.csv("env.pca.w.prey.vectors.csv")
(env.pca.vectors_trophic<-round(rbind("Eigenvalues"=lake.pca.trophic$CA$eig,"Singular Values"=sqrt(lake.pca.trophic$CA$eig),
                                      "Proportion Explained"=lake.pca.trophic$CA$eig/sum(lake.pca.trophic$CA$eig),
                                      lake.pca.trophic$CA$u,lake.pca.trophic$CA$v),3))%>%write.csv("env.pca.trophic.vectors.csv")


# Diet plot showing Daphnia and gammarid abundance w/ % littoral area as the color scale

benthlim.plot<-ggplot(data = eco.dat,aes(x=Gamm.,y=Daph., color=littoral*100)) + theme_classic()+
  labs(x="Gammarids",y="Daphnia",color="% littoral")+
  geom_point(size=2.5)+ scale_color_gradient(low="steelblue4",high = "seagreen2")+
  geom_text_repel(data = eco.dat,aes(label=Lake),color="black",box.padding=.5,point.padding = .2,size=2.5)




#
# MAP
library(rnaturalearth)
library(rnaturalearthdata)
library(ggrepel)

#map
#loads map data
world <- ne_countries(scale=50,returnclass = 'sf')
usa <- subset(world, admin == "United States of America")
can <- subset(world, admin == "Canada")
rus <- subset(world, admin == "Russia")
usacanrus<-rbind(usa,can,rus)
eco.dat$Lake

#inset map of entire state of AK
alaska<-ggplot(data = usacanrus) +
  geom_sf(fill = "honeydew2") +
  panel_border(color = "grey50")+ theme_grey()+
  geom_rect(xmin = -152.5, xmax = -147.5, ymin = 59, ymax = 62, 
            fill = NA, colour = "black", size = .45)+
  coord_sf(xlim = c(-176, -130), 
           ylim = c(53, 73), expand = FALSE, datum = NA)
# Kenai Peninsula local map
alaskalocal <- ggplot(data = usa) + theme_cowplot(font_size = 9)+
  geom_sf(fill = "honeydew2") +coord_sf(xlim = c(-152.5, -147.5), ylim = c(59, 62), expand = F)
# combines inset and local Kenai map, and adds points and labels
AKlocal_inset<-alaskalocal + annotation_custom(
  grob = ggplotGrob(alaska),
  xmin = -149.4,
  xmax = -147.5,
  ymin = 57.9,
  ymax = 59.8+1.2) +
  geom_point(data=eco.dat,aes(x=Lon.,y=Lat.),size = 2, color = "skyblue4")+
  geom_text_repel(data = eco.dat,aes(x=Lon.,y=Lat.,label = Lake),
                  box.padding=.5,point.padding = .2,size = 3, fontface = "bold")
ggsave("AK lake map.pdf", width = 5, height = 6, dpi = "screen")

#
#
#  Groups of environmental traits 
prey<-c("Daph.","Gamm.","littoral")
physical<-c("Area","Max_Depth","DOC")
chemical<-c("TP","TN","ChlA","Cond.","pH","Ca")
env<-c(prey,physical,chemical)


get.lake.means<-function(traits,traitnames,eco,eco.var){
  traits.eco=transform(merge(traits,eco[,eco.var],by.x="Lake",by.y = "row.names", all = T),row.names = ID)
  lake.means=group_by(na.omit(as.data.table(traits.eco[,c(-(2))]),cols = traitnames),
                      Lake)%>%summarise_all(mean,na.rm=T)%>%as.data.frame() 
  rownames(lake.means)<-lake.means$Lake
  lake.means
}

def.lake.means<-get.lake.means(Defense.df,def.trait.list,eco.dat,env)
swi.lake.means<-get.lake.means(Swimming.df,swi.trait.list,eco.dat,env)
tro.lake.means<-get.lake.means(Trophic.df,tro.trait.list,eco.dat,env)



#
eco.dat<-data.frame(eco.dat[,])# exclude Echo, b/c partial eco data
# POPULATION MEANS FOR TRAITS
adj.data.unscaled<-merge(Defense.df[c("ID","Lake",def.trait.list)],
                         merge(Swimming.df[c("ID",swi.trait.list)],merge(Trophic.df[,c("ID",tro.trait.list)],two.d.coords, 
                             by.x = "ID",by.y = "row.names", all = T), by= "ID",all= T), by = "ID", all = T)


lake.means.adj.data<-group_by(na.omit(adj.data.unscaled[,c(-1,-(19:60))]),Lake)%>%summarise_all(mean,na.rm=T)%>%as.data.frame()# non-landmark traits only
lake.means.data.eco<-merge(lake.means.adj.data,eco.dat[,env],by.x = "Lake", by.y = "row.names")
rownames(lake.means.data.eco)<-lake.means.data.eco$Lake

LD.means<-group_by(merge(merge(LD.def.df,LD.swim.df[-2],by="FishID"),LD.tro.df[-2],by="FishID"),Lake)%>%
  summarise_all(mean,na.rm=T)%>%as.data.frame()%>%merge(eco.dat[,env],by.x = "Lake", by.y = "row.names")



lake.means.4tile<-as.matrix(lake.means.data.eco[,-1])
colnames(lake.means.4tile)<-c("LP","DS1l","DS2l","PSl","LPm","SL",
                              "BD","CP","BC","GW","PW","GR","JL","SN","ED","HL",env)

trait.corr.tile<-rcorr(lake.means.4tile,type="spearman")[[1]][c("LP","DS1l","DS2l",
              "PSl","LPm","SL","BD","CP","BC","GW","PW","GR","JL","SN","ED","HL"),env]%>%
  corrplot(method ="color",addgrid.col=NA,sig.level=c(.001,.01,.05), insig = "label_sig",tl.col = "black",tl.cex = .5,pch="",
           pch.cex = .75,p.mat = rcorr(lake.means.4tile,type="spearman")[[3]][c("LP","DS1l",
               "DS2l","PSl","LPm","SL","BD","CP","BC","GW","PW","GR","JL","SN","ED","HL"),env])

LD.corr.tile<-rcorr(as.matrix(LD.means[,-c(1:2)]),type="spearman")[[1]][
  c("defense.LD1","swimming.LD1", "trophic.LD1","defense.LD2","swimming.LD2","trophic.LD2","defense.LD3","swimming.LD3","trophic.LD3"),env]%>%
  corrplot(method ="color",addgrid.col=NA,sig.level=c(.001,.01,.05), tl.col = "black",tl.cex = .5,pch="",pch.cex = .75,
           p.mat = rcorr(as.matrix(LD.means[,-c(1:2)]),type="spearman")[[1]][
             c("defense.LD1","swimming.LD1", "trophic.LD1","defense.LD2","swimming.LD2","trophic.LD2","defense.LD3","swimming.LD3","trophic.LD3"),env])


# TRAIT X ENV PLS

# z-transforms trait and eco data 
lake.means.data.eco_adj<-data.frame(row.names=(lake.means.data.eco)$Lake,
                                    sapply((lake.means.data.eco[,-1]), function(x) (x-mean(na.omit(x)))/sd(na.omit(x))))

pls.vector.df<-function(pls){
  round(data.frame(rbind("singular values"=pls$svd$d,"SV sq"=(pls$svd$d)^2,
                         pls$left.pls.vectors,pls$right.pls.vectors)),4)
}

# 2-B PLS for Defense, Swimming, and Trophic traits
(eco.def.pls<-two.b.pls(na.omit(lake.means.data.eco_adj)[,env],
                        na.omit(lake.means.data.eco_adj)[,def.trait.list]))%>%plot()
eco.def.pls.bestfit<-odregress(eco.def.pls$XScores[,1],eco.def.pls$YScores[,1])$coeff[c(2,1)]
(eco.swi.pls<-two.b.pls(na.omit(lake.means.data.eco_adj)[,env],
                        na.omit(lake.means.data.eco_adj)[,swi.trait.list]))%>%plot()
eco.swi.pls.bestfit<-odregress(eco.swi.pls$XScores[,1],eco.swi.pls$YScores[,1])$coeff[c(2,1)]
(eco.tro.pls<-two.b.pls(na.omit(lake.means.data.eco_adj)[,env],
                        na.omit(lake.means.data.eco_adj)[,tro.trait.list]))%>%plot()
eco.tro.pls.bestfit<-odregress(eco.tro.pls$XScores[,1],eco.tro.pls$YScores[,1])$coeff[c(2,1)]

# full (except shape)
(eco.full.pls<-two.b.pls(na.omit(lake.means.data.eco_adj)[,env],
                         na.omit(lake.means.data.eco_adj)[,c(def.trait.list,swi.trait.list,tro.trait.list)]))%>%plot()
eco.full.pls_bestfit<-odregress(eco.full.pls$XScores[,1],eco.full.pls$YScores[,1])$coeff[c(2,1)]

(eco.def.vectors<-round(pls.vector.df(eco.def.pls),3))%>%write.csv("eco-def_PLSvec.csv")
(eco.swi.vectors<-round(pls.vector.df(eco.swi.pls),3))%>%write.csv("eco-swi_PLSvec.csv")
(eco.tro.vectors<-round(pls.vector.df(eco.tro.pls),3))%>%write.csv("eco-tro_PLSvec.csv")
(eco.full.vectors<-round(pls.vector.df(eco.full.pls),3))%>%write.csv("eco-full_PLSvec.csv")

PhenoEco.stats<-data.frame(rbind(c(eco.def.pls$r.pls,eco.def.pls$P.value,eco.def.pls$Z),
                                 c(eco.swi.pls$r.pls,eco.swi.pls$P.value,eco.swi.pls$Z),
                                 c(eco.tro.pls$r.pls,eco.tro.pls$P.value,eco.tro.pls$Z),
                                 c(eco.full.pls$r.pls,eco.full.pls$P.value,eco.full.pls$Z)),
                           row.names=c("Defense x Eco","Swimming x Eco","Trophic x Eco",
                                       "Def-Swi-Tro x Eco"))
colnames(PhenoEco.stats)<-c("r-PLS","p","Z")
PhenoEco.stats[1:3,"adj.p"]<-p.adjust(PhenoEco.stats[1:3,"p"], method = "fdr")
PhenoEco.stats<-round(PhenoEco.stats,3)
write.csv(PhenoEco.stats,"Pheno-Eco_PLS-stats.csv")





# 2-B PLS for Defense, Swimming, and Trophic traits W/O feeding-related env. traits
(eco.physchem.def.pls<-two.b.pls(na.omit(lake.means.data.eco_adj[,env[-c(1:3)]]),
                                 na.omit(lake.means.data.eco_adj[,c(env[-c(1:3)],def.trait.list)])[,def.trait.list]))%>%plot()
eco.physchem.def.pls.bestfit<-odregress(eco.physchem.def.pls$XScores[,1],eco.physchem.def.pls$YScores[,1])$coeff[c(2,1)]
(eco.physchem.swi.pls<-two.b.pls(na.omit(lake.means.data.eco_adj[,env[-c(1:3)]]),
                                 na.omit(lake.means.data.eco_adj[,c(env[-c(1:3)],swi.trait.list)])[,swi.trait.list]))%>%plot()
eco.physchem.swi.pls.bestfit<-odregress(eco.physchem.swi.pls$XScores[,1],eco.physchem.swi.pls$YScores[,1])$coeff[c(2,1)]
(eco.physchem.tro.pls<-two.b.pls(na.omit(lake.means.data.eco_adj[,env[-c(1:3)]]),
                                 na.omit(lake.means.data.eco_adj[,c(env[-c(1:3)],tro.trait.list)])[,tro.trait.list]))%>%plot()
eco.physchem.tro.pls.bestfit<-odregress(eco.physchem.tro.pls$XScores[,1],eco.physchem.tro.pls$YScores[,1])$coeff[c(2,1)]

# full (except shape)
(eco.physchem.full.pls<-two.b.pls(na.omit(lake.means.data.eco_adj[,env[-c(1:3)]]),
                                  na.omit(lake.means.data.eco_adj[,c(env[-c(1:3)],
                                        def.trait.list,swi.trait.list,tro.trait.list)])
                                  [,c(def.trait.list,swi.trait.list,tro.trait.list)]))%>%plot()
eco.physchem.full.pls.bestfit<-odregress(eco.physchem.full.pls$XScores[,1],eco.physchem.full.pls$YScores[,1])$coeff[c(2,1)]

(eco.physchem.def.vectors<-round(pls.vector.df(eco.physchem.def.pls),3))%>%write.csv("eco.physchem-def_PLSvec.csv")
(eco.physchem.swi.vectors<-round(pls.vector.df(eco.physchem.swi.pls),3))%>%write.csv("eco.physchem-swi_PLSvec.csv")
(eco.physchem.tro.vectors<-round(pls.vector.df(eco.physchem.tro.pls),3))%>%write.csv("eco.physchem-tro_PLSvec.csv")
(eco.physchem.full.vectors<-round(pls.vector.df(eco.physchem.full.pls),3))%>%write.csv("eco.physchem-full_PLSvec.csv")

PhenoEco.physchem.stats<-data.frame(rbind(c(eco.physchem.def.pls$r.pls,eco.physchem.def.pls$P.value,eco.physchem.def.pls$Z),
                                          c(eco.physchem.swi.pls$r.pls,eco.physchem.swi.pls$P.value,eco.physchem.swi.pls$Z),
                                          c(eco.physchem.tro.pls$r.pls,eco.physchem.tro.pls$P.value,eco.physchem.tro.pls$Z),
                                          c(eco.physchem.full.pls$r.pls,eco.physchem.full.pls$P.value,eco.physchem.full.pls$Z)),
                                    row.names=c("Defense x Eco","Swimming x Eco","Trophic x Eco",
                                                "Def-Swi-Tro x Eco"))
colnames(PhenoEco.physchem.stats)<-c("r-PLS","p","Z")
PhenoEco.physchem.stats[1:3,"adj.p"]<-p.adjust(PhenoEco.physchem.stats[1:3,"p"], method = "fdr")
PhenoEco.physchem.stats<-round(PhenoEco.physchem.stats,3)
write.csv(PhenoEco.physchem.stats,"Pheno-Eco.physchem_PLS-stats.csv")



####### THIS EXTRACTS ORTHOGONAL REGRESSION LINE!!!!
eco.pls.plot<-function(pls,bestfit,xlab,ylab){
  pls.dat=data.frame(x=pls$XScores[,1],y=pls$YScores[,1])
  ggplot(data=pls.dat,aes(x=x,y=y))+ theme_classic()+labs(x=xlab,y=ylab)+
    geom_abline(intercept=bestfit[1],slope = bestfit[2], color = "grey80")+
    geom_point(color="seagreen3", size=2)+
    geom_text_repel(aes(label=rownames(pls.dat)),color="black",box.padding=.3,point.padding = .03,size=2.5)
}

eco.def.pls.plot<-eco.pls.plot(eco.def.pls,eco.def.pls.bestfit,"L Block - Eco","R Block - Defense")
eco.swi.pls.plot<-eco.pls.plot(eco.swi.pls,eco.swi.pls.bestfit,"L Block - Eco","R Block - Swimming")
eco.tro.pls.plot<-eco.pls.plot(eco.tro.pls,eco.tro.pls.bestfit,"L Block - Eco","R Block - Trophic")
eco.full.pls.plot<-eco.pls.plot(eco.full.pls,eco.full.pls_bestfit,"L Block - Eco","R Block - Pheno")
# compound pls plot for full eco variables including littoral, daph., and gamm. 
eco.def.pls.plot+eco.tro.pls.plot+eco.swi.pls.plot+eco.full.pls.plot

def.pls.eig<-ED((eco.def.pls$svd$d)^2)
swi.pls.eig<-ED((eco.swi.pls$svd$d)^2)
tro.pls.eig<-ED((eco.tro.pls$svd$d)^2)
full.pls.eig<-ED((eco.full.pls$svd$d)^2)



eco.physchem.def.pls.plot<-eco.pls.plot(eco.physchem.def.pls,eco.physchem.def.pls.bestfit,"L Block - Eco","R Block - Defense")
eco.physchem.swi.pls.plot<-eco.pls.plot(eco.physchem.swi.pls,eco.physchem.swi.pls.bestfit,"L Block - Eco","R Block - Swimming")
eco.physchem.tro.pls.plot<-eco.pls.plot(eco.physchem.tro.pls,eco.physchem.tro.pls.bestfit,"L Block - Eco","R Block - Trophic")
eco.physchem.full.pls.plot<-eco.pls.plot(eco.physchem.full.pls,eco.physchem.full.pls.bestfit,"L Block - Eco","R Block - Pheno")
eco.physchem.def.pls.plot+eco.physchem.tro.pls.plot+
  eco.physchem.swi.pls.plot+eco.physchem.full.pls.plot

def.physchem.pls.eig<-ED((eco.physchem.def.pls$svd$d)^2)
swi.physchem.pls.eig<-ED((eco.physchem.swi.pls$svd$d)^2)
tro.physchem.pls.eig<-ED((eco.physchem.tro.pls$svd$d)^2)
full.physchem.pls.eig<-ED((eco.physchem.full.pls$svd$d)^2)




ecoblock.theta<-function(vec1,vec2){
  theta.180=vec.theta(vec1,vec2)
  theta=if(theta.180>90) 180-theta.180 else theta.180
  theta
}
eco.pls.theta<-data.frame(matrix(as.numeric(c(ecoblock.theta(eco.def.pls$left.pls.vectors[,1],eco.swi.pls$left.pls.vectors[,1]),"",
                                              ecoblock.theta(eco.def.pls$left.pls.vectors[,1],eco.tro.pls$left.pls.vectors[,1]),
                                              ecoblock.theta(eco.swi.pls$left.pls.vectors[,1],eco.tro.pls$left.pls.vectors[,1]))),nrow = 2,byrow = T,),
                          row.names  = c("Swimming","Trophic"))
colnames(eco.pls.theta)<-c("Defense","Swimming")
eco.pls.theta<-round(eco.pls.theta,2)
write.csv(eco.pls.theta,"pw.pls.eco_AngleMatrix.csv")

ecoblock.theta(eco.full.pls$left.pls.vectors[,1],eco.def.pls$left.pls.vectors[,1])%>%round(2)
ecoblock.theta(eco.full.pls$left.pls.vectors[,1],eco.swi.pls$left.pls.vectors[,1])%>%round(2)
ecoblock.theta(eco.full.pls$left.pls.vectors[,1],eco.tro.pls$left.pls.vectors[,1])%>%round(2)


eco.physchem.pls.theta<-data.frame(matrix(as.numeric(c(ecoblock.theta(eco.physchem.def.pls$left.pls.vectors[,1],eco.physchem.swi.pls$left.pls.vectors[,1]),"",
                                                       ecoblock.theta(eco.physchem.def.pls$left.pls.vectors[,1],eco.physchem.tro.pls$left.pls.vectors[,1]),
                                                       ecoblock.theta(eco.physchem.swi.pls$left.pls.vectors[,1],eco.physchem.tro.pls$left.pls.vectors[,1]))),nrow = 2,byrow = T,),
                                   row.names  = c("Swimming","Trophic"))
colnames(eco.physchem.pls.theta)<-c("Defense","Swimming")
eco.physchem.pls.theta<-round(eco.physchem.pls.theta,2)
write.csv(eco.physchem.pls.theta,"pw.pls.eco.physchem_AngleMatrix.csv")

ecoblock.theta(eco.physchem.full.pls$left.pls.vectors[,1],eco.physchem.def.pls$left.pls.vectors[,1])%>%round(2)
ecoblock.theta(eco.physchem.full.pls$left.pls.vectors[,1],eco.physchem.swi.pls$left.pls.vectors[,1])%>%round(2)
ecoblock.theta(eco.physchem.full.pls$left.pls.vectors[,1],eco.physchem.tro.pls$left.pls.vectors[,1])%>%round(2)






