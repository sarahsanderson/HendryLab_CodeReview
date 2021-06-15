# LD frequency plot function - plots freq. distribution of LD1 for each lake
# data= dataframe containing Lake vector and LD1 scores
# X = LD1 vector (e.g., data$LD1)
# Lake = Lake vector (e.g., data$Lake)
# Lakelab = title for lake legend ("Lake")
# title = plot title (as character string)
# LD.prop = % of trace explained by LD1 (as character string)

LD.freq<-function(data,X,Lake,Lakelab,title,LD.prop){data %>%
    mutate(Lake = fct_reorder(Lake, X, .fun='mean' )) %>%ggplot( aes(x = X, y = Lake, fill = ..x..)) + 
    scale_fill_viridis(option = "C") +
    xlab(paste("LD1",LD.prop))+ylab(Lakelab)+
    geom_density_ridges_gradient(scale = 1.85, rel_min_height = 0.01,jittered_points = TRUE, 
                                 point_shape = "|", point_size = 2, size = 0.25,
                                 position = position_points_jitter(height = 0)) +
    labs(title = title)+
    theme_minimal_hgrid()+theme(legend.position="none")
}

# function threeD.scatter makes a 3d scatter plot
# data = a dataframe including groupinf factor (e.g., "Lake") 
#   and at least 3 vectors of PC scores, LD scores, etc.
# anchor.axis = axis droplines should run along
# anchor.val = value along anchor.axis at which to start droplines
# col.vec = a vector of color names to apply to groups
# axis.index = vector of column indices of "data" that contain X, Y, and Z values
# axis.labs = X, Y, and Z axis labels
threeD.scatter<-function(data,anchor.axis,anchor.val,col.vec,axis.index,axis.labs){
  dat.2 <- replicate(2, data, simplify = F)
  dat.2[[2]][anchor.axis] <- anchor.val
  dat.3 <- dat.2 %>%bind_rows() %>%group2NA("FishID", "Lake")
  
  plot_ly(color = ~factor(Lake), showlegend = T, 
          colors = col.vec) %>%
    add_markers(data = data, x = ~data[,axis.index[1]], y = ~data[,axis.index[2]], z = ~data[,axis.index[3]],size = 100) %>%
    add_paths(data = dat.3,  x = ~dat.3[,axis.index[1]], y = ~dat.3[,axis.index[2]], z = ~dat.3[,axis.index[3]], opacity = 0.1,showlegend =F)%>%
    layout(scene=list(xaxis=list(title = axis.labs[[1]]),yaxis=list(title = axis.labs[[2]]), zaxis=list(title= axis.labs[[3]])))
  
}

# data = a dataframe including groupinf factor (e.g., "Lake") 
#   and at least 3 vectors of PC scores, LD scores, etc.
# axis.index = vector of column indices of "data" that contain X, Y, and Z values
# col.vec = a vector of color names to apply to groups
# vector.sum.df = dataframe summarizing magnitude and loadings of LDA vectors
#   row #2 is proportion of trace
# axis.labs = vector of X and Y axis names, pasted to proportion of trace for axis label,
#   (e.g., c("LD1 - ","LD2 - "))
# which.vec = indices of columns representing vectors on X and Y axes in vector.sum.df
#   (e.g., for X = LD1 and Y = LD2, c(1,2), for X = LD1 and Y = LD3, c(1,3))
twoD.scatter<-function(data,axis.index,col.vec,vector.sum.df,axis.labs,which.vec){
  ggplot(data = data, aes(x=data[,axis.index[1]],y=data[,axis.index[2]],color=Lake))+
    geom_point(alpha=.7,stroke=0,size=2)+theme_classic()+
    scale_color_manual(values=col.vec)+coord_fixed()+
  guides(color = guide_legend(override.aes = list(size = 4, alpha=1)))+
  labs(x=paste(axis.labs[1],round(vector.sum.df[2,which.vec[1]]/sum(vector.sum.df[2,])*100,1),"%"),
       y=paste(axis.labs[2],round(vector.sum.df[2,which.vec[2]]/sum(vector.sum.df[2,])*100,1),"%"))
}


### Effective dimension matrix
dim.stats<-function(lda,euc.mat,lake.means){
  data.frame(dim.SVsq=estimate.ED.eig(prop.trace(lda$svd)),
             dim.matrix=unlist(estimate.ED(na.omit(euc.mat))),
             pop.mean.matrix=unlist(estimate.ED(subset(lake.means,select=-Lake))))
}

#LDAs#
#########
######    
#########

#### Defense ####
#Defense data frame#
Defense.df <- cbind.data.frame(ID=data$FishID,Lake= data$Lake, adj.PS_mean = data$adj.PS.mean, 
                               LatPlNum_mean = rowMeans(subset(data, select = c("LatPlNum_R", "LatPlNum_L")), na.rm = TRUE), 
                               MidPl_mean=rowMeans(subset(data, select = c("FirstPl_L", "LastPl_L","FirstPl_R","LastPl_R")), na.rm = TRUE), 
                               adj.DS1L= data$adj.DS1L,adj.DS2L = data$adj.DS2L)
rownames(Defense.df)<-Defense.df$ID


# LDAprep.A) This chunk splits Defense.df by Lake, then interpolates missing trait values
#   within lake, unlists the interpolated vector, and replaces the 
#   old vector with the interpolated one (Defense.df is already sorted by lake)
Defense.df.split<-split(Defense.df,Defense.df$Lake, drop = F)
NEW.adj.DS1L<-unlist(lapply(Defense.df.split,function(x){
  interp(data=x,y=x[,"adj.DS1L"],a=x[,"adj.PS_mean"],b=x[,"adj.DS2L"])$SecondaryNew.Vec}))
NEW.adj.DS2L<-unlist(lapply(Defense.df.split,function(x){
  interp(data=x,b=x[,"adj.DS1L"],a=x[,"adj.PS_mean"],y=x[,"adj.DS2L"])$SecondaryNew.Vec}))
NEW.adj.PS_mean<-unlist(lapply(Defense.df.split,function(x){
  if(var(na.omit(x$adj.PS_mean))==0)x$adj.PS_mean else(interp(data=x,b=x[,"adj.DS1L"],y=x[,"adj.PS_mean"],a=x[,"adj.DS2L"])$SecondaryNew.Vec)}))
Defense.df<-bind_rows(Defense.df.split)
Defense.df[,"adj.DS1L"]<-NEW.adj.DS1L
Defense.df[,"adj.DS2L"]<-NEW.adj.DS2L
Defense.df[,"adj.PS_mean"]<-NEW.adj.PS_mean
rownames(Defense.df)<-Defense.df$ID

# LDAprep.B) scaling to sd and mean-centering of traits
Defense.df.euc<-data.frame(Defense.df["ID"],Defense.df["Lake"],
                           LatPlNum_mean=((Defense.df[,"LatPlNum_mean"]*1-mean(na.omit(Defense.df[,"LatPlNum_mean"])))/sd(na.omit(Defense.df[,"LatPlNum_mean"]))),
                           adj.DS1L=((Defense.df[,"adj.DS1L"]*1-mean(na.omit(Defense.df[,"adj.DS1L"])))/sd(na.omit(Defense.df[,"adj.DS1L"]))),
                           adj.DS2L=((Defense.df[,"adj.DS2L"]*1-mean(na.omit(Defense.df[,"adj.DS2L"])))/sd(na.omit(Defense.df[,"adj.DS2L"]))),
                           adj.PS_mean=((Defense.df[,"adj.PS_mean"]*1-mean(na.omit(Defense.df[,"adj.PS_mean"])))/sd(na.omit(Defense.df[,"adj.PS_mean"]))),
                           MidPl_mean=((Defense.df[,"MidPl_mean"]*1-mean(na.omit(Defense.df[,"MidPl_mean"])))/sd(na.omit(Defense.df[,"MidPl_mean"]))))
# makes numeric matrix for dimensionality based on eigendecomp. (by individuals, then means)
Defense.df.euc.mat<-Defense.df.euc[,c("adj.PS_mean" ,"LatPlNum_mean","MidPl_mean","adj.DS1L", "adj.DS2L")]
def.means<-aggregate(data=na.omit(Defense.df.euc[,c("Lake","LatPlNum_mean","adj.DS1L","adj.DS2L", "adj.PS_mean", "MidPl_mean")]), .~Lake, FUN=mean)

########## LDA
Defense.lda<-lda(Lake ~ adj.PS_mean + LatPlNum_mean+
                   MidPl_mean+  adj.DS1L + adj.DS2L , data = Defense.df.euc, CV = F)

predictions.Defense <- predict(Defense.lda)
LD.def.df<-data.frame(FishID=na.omit(Defense.df.euc)$ID,Lake=na.omit(Defense.df.euc)$Lake, defense.LD1=predictions.Defense$x[,1],
                      defense.LD2=predictions.Defense$x[,2], defense.LD3=predictions.Defense$x[,3]) 
# LDA dimension summary
(lda.def.df<-round(lda.vector.df(Defense.lda),3))%>%write.csv("lda-def_df.csv")
# Effective dimensions (Defense)
defense.dim.matrix<-dim.stats(Defense.lda,Defense.df.euc.mat,def.means)


### Defense LDA dim. Plots
# Defense LD1 frequency plot
freq.def<-LD.freq(LD.def.df,LD.def.df$defense.LD1,LD.def.df$Lake,
        "Lake",'Discriminant Score of Defensive Traits by Lake',"92.1%")
# 3D scatter
threeD.scatter(LD.def.df,"defense.LD3",-4,c("steelblue4","firebrick3","palegreen4","palegreen2","steelblue2",
                                            "firebrick1","turquoise4","turquoise1","tan3","tan1","darkgreen",
                                            "sienna4","grey80","grey50"),c(3,4,5),c('LD1-92.1%', 'LD2-4.9%','LD3-1.6%'))
# 2D scatter 

def.plot.12<-twoD.scatter(LD.def.df,c(3,4),c("steelblue4","firebrick3","palegreen4","palegreen2","steelblue2",
                                "firebrick1","turquoise4","turquoise1","tan3","tan1","darkgreen",
                                "sienna4","grey80","grey50"),lda.def.df,c("LD1 - ","LD2 - "),c(1,2))
def.plot.13<-twoD.scatter(LD.def.df,c(3,5),c("steelblue4","firebrick3","palegreen4","palegreen2","steelblue2",
                                             "firebrick1","turquoise4","turquoise1","tan3","tan1","darkgreen",
                                             "sienna4","grey80","grey50"),lda.def.df,c("LD1 - ","LD3 - "),c(1,3))
def.plot.12/def.plot.13+ plot_layout(guides = 'collect')


##################
#### Swimming ####

####LDA Swimming Trait####
#Swimming data frame#
Swimming.df <- cbind.data.frame(ID=data$FishID,Lake= data$Lake, 
                                SL_mm = log(data$SL_mm), 
                                adj.BD= data$adj.BD, 
                                adj.CP= data$adj.CP)
rownames(Swimming.df)<-Swimming.df$ID

# LDAprep.A) # SL data not interpolated because all 
#   individuals have SL data (and traits size corrected based on SL.
#   Body depth not interpolated because all 
#   individulas with SL data also have BD data
Swimming.df.split<-split(Swimming.df,Swimming.df$Lake, drop = F)
NEW.adj.CP<-unlist(lapply(Swimming.df.split,function(x){
  interp(data=x,y=x[,"adj.CP"],a=x[,"SL_mm"],b=x[,"adj.BD"])$SecondaryNew.Vec}))
Swimming.df<-bind_rows(Swimming.df.split)
Swimming.df[,"adj.CP"]<-NEW.adj.CP 
rownames(Swimming.df)<-Swimming.df$ID

# LDAprep.B) scaling to sd and mean-centering of traits
Swimming.df.euc<-data.frame(ID=Swimming.df["ID"],Lake=Swimming.df["Lake"],
                            SL_mm=((Swimming.df[,"SL_mm"]*1-mean(na.omit(Swimming.df[,"SL_mm"])))/sd(na.omit(Swimming.df[,"SL_mm"]))),
                            adj.CP=((Swimming.df[,"adj.CP"]*1-mean(na.omit(Swimming.df[,"adj.CP"])))/sd(na.omit(Swimming.df[,"adj.CP"]))),
                            adj.BD=((Swimming.df[,"adj.BD"]*1-mean(na.omit(Swimming.df[,"adj.BD"])))/sd(na.omit(Swimming.df[,"adj.BD"]))))
# makes numeric matrix for dimensionality based on eigendecomp. (by individuals, then means)
Swimming.df.euc.mat<-Swimming.df.euc[,c("SL_mm","adj.CP","adj.BD")]
swi.means<-aggregate(data=na.omit(Swimming.df.euc[,c("SL_mm","adj.CP","adj.BD")]), .~Lake, FUN=mean)


########## LDA Swimming
Swimming.lda <- lda(Lake ~ SL_mm + adj.BD + adj.CP, data = Swimming.df.euc)
predictions.Swimming <- predict(Swimming.lda)
LD.swim.df<-data.frame(FishID=na.omit(Swimming.df.euc)$ID,Lake=na.omit(Swimming.df.euc)$Lake, swimming.LD1=predictions.Swimming$x[,1],
                       swimming.LD2=predictions.Swimming$x[,2], swimming.LD3=predictions.Swimming$x[,3]) 
# LDA dimension summary
(lda.swi.df<-round(lda.vector.df(Swimming.lda),3))%>%write.csv("lda-swi_df.csv")
# Effective dimensions (Swimming)
swimming.dim.matrix<-dim.stats(Swimming.lda,Swimming.df.euc.mat,swi.means)


### Swimming LDA dim. Plots
# Swimming LD1 frequency plot
freq.swi<-LD.freq(LD.swim.df,LD.swim.df$swimming.LD1,LD.swim.df$Lake,
        "Lake","Discriminant Score of Swimming Traits by Lake","67.5%")
# 3D scatter
threeD.scatter(LD.swim.df,"swimming.LD3",-4,c("steelblue4","firebrick3","palegreen4","palegreen2","steelblue2",
                                              "firebrick1","turquoise4","turquoise1","tan3","tan1","darkgreen",
                                              "sienna4","grey80","grey50"),c(3,4,5),c('LD1-67.5%', 'LD2-25.6%','LD3-6.9%'))

# 2D scatter 
swi.plot.12<-twoD.scatter(LD.swim.df,c(3,4),c("steelblue4","firebrick3","palegreen4","palegreen2","steelblue2",
                                             "firebrick1","turquoise4","turquoise1","tan3","tan1","darkgreen",
                                             "sienna4","grey80","grey50"),lda.swi.df,c("LD1 - ","LD2 - "),c(1,2))



swi.plot.13<-twoD.scatter(LD.swim.df,c(3,5),c("steelblue4","firebrick3","palegreen4","palegreen2","steelblue2",
                                             "firebrick1","turquoise4","turquoise1","tan3","tan1","darkgreen",
                                             "sienna4","grey80","grey50"),lda.swi.df,c("LD1 - ","LD3 - "),c(1,3))
swi.plot.12/swi.plot.13+ plot_layout(guides = 'collect')



##################
#### Trophic ####

####LDA Trophic Traits####
#Trophic data frame#
Trophic.df <- cbind.data.frame(ID=data$FishID,Lake= data$Lake, 
                               adj.BCL= data$adj.BCL, 
                               adj.GW= data$adj.GW, 
                               adj.PW= data$adj.PW, 
                               RakerNum= data$RakerNum, 
                               adj.JL= data$adj.JL,
                               adj.SnL= data$adj.SnL, 
                               adj.ED= data$adj.ED, 
                               adj.HL= data$adj.HL)
rownames(Trophic.df)<-Trophic.df$ID

# LDAprep.A)  RakerNum not interpolated because meristic (not continuous)
Trophic.df.split<-split(Trophic.df,Trophic.df$Lake, drop = F)
NEW.adj.GW<-unlist(lapply(Trophic.df.split,function(x){
  interp(data=x,y=x[,"adj.GW"],a=x[,"adj.BCL"],b=x[,"adj.JL"], c=x[,"adj.PW"],d=x[,"adj.SnL"],e=x[,"adj.ED"], f=x[,"adj.HL"])$SecondaryNew.Vec}))
NEW.adj.PW<-unlist(lapply(Trophic.df.split,function(x){
  interp(data=x,y=x[,"adj.PW"],a=x[,"adj.BCL"],b=x[,"adj.JL"],c=x[,"adj.GW"],d=x[,"adj.SnL"],e=x[,"adj.ED"], f=x[,"adj.HL"])$SecondaryNew.Vec}))
NEW.adj.JL<-unlist(lapply(Trophic.df.split,function(x){
  interp(data=x,y=x[,"adj.JL"],a=x[,"adj.BCL"],b=x[,"adj.PW"],c=x[,"adj.GW"],d=x[,"adj.SnL"],e=x[,"adj.ED"], f=x[,"adj.HL"])$SecondaryNew.Vec}))
NEW.adj.BCL<-unlist(lapply(Trophic.df.split,function(x){
  interp(data=x,y=x[,"adj.BCL"],a=x[,"adj.JL"],b=x[,"adj.PW"],c=x[,"adj.GW"],d=x[,"adj.SnL"],e=x[,"adj.ED"], f=x[,"adj.HL"])$SecondaryNew.Vec}))
NEW.adj.SnL<-unlist(lapply(Trophic.df.split,function(x){
  interp(data=x,y=x[,"adj.SnL"],a=x[,"adj.JL"],b=x[,"adj.PW"],c=x[,"adj.GW"],d=x[,"adj.BCL"],e=x[,"adj.ED"], f=x[,"adj.HL"])$SecondaryNew.Vec}))
NEW.adj.ED<-unlist(lapply(Trophic.df.split,function(x){
  interp(data=x,y=x[,"adj.ED"],a=x[,"adj.JL"],b=x[,"adj.PW"],c=x[,"adj.GW"],d=x[,"adj.SnL"],e=x[,"adj.BCL"], f=x[,"adj.HL"])$SecondaryNew.Vec}))
NEW.adj.HL<-unlist(lapply(Trophic.df.split,function(x){
  interp(data=x,y=x[,"adj.HL"],a=x[,"adj.JL"],b=x[,"adj.PW"],c=x[,"adj.GW"],d=x[,"adj.SnL"],e=x[,"adj.ED"], f=x[,"adj.BCL"])$SecondaryNew.Vec}))
Trophic.df<-bind_rows(Trophic.df.split)
Trophic.df[,"adj.GW"]<-NEW.adj.GW
Trophic.df[,"adj.PW"]<-NEW.adj.PW
Trophic.df[,"adj.JL"]<-NEW.adj.JL
Trophic.df[,"adj.BCL"]<-NEW.adj.BCL
Trophic.df[,"adj.SnL"]<-NEW.adj.SnL
Trophic.df[,"adj.ED"]<-NEW.adj.ED
Trophic.df[,"adj.HL"]<-NEW.adj.HL

# LDAprep.B)
Trophic.df.euc<-data.frame(ID=Trophic.df["ID"],Lake=Trophic.df["Lake"],
                           adj.GW=((Trophic.df[,"adj.GW"]*1-mean(na.omit(Trophic.df[,"adj.GW"])))/sd(na.omit(Trophic.df[,"adj.GW"]))),
                           adj.PW=((Trophic.df[,"adj.PW"]*1-mean(na.omit(Trophic.df[,"adj.PW"])))/sd(na.omit(Trophic.df[,"adj.PW"]))),
                           RakerNum=((Trophic.df[,"RakerNum"]*1-mean(na.omit(Trophic.df[,"RakerNum"])))/sd(na.omit(Trophic.df[,"RakerNum"]))),
                           adj.BCL=((Trophic.df[,"adj.BCL"]*1-mean(na.omit(Trophic.df[,"adj.BCL"])))/sd(na.omit(Trophic.df[,"adj.BCL"]))),
                           adj.JL=((Trophic.df[,"adj.JL"]*1-mean(na.omit(Trophic.df[,"adj.JL"])))/sd(na.omit(Trophic.df[,"adj.JL"]))),
                           adj.SnL=((Trophic.df[,"adj.SnL"]*1-mean(na.omit(Trophic.df[,"adj.SnL"])))/sd(na.omit(Trophic.df[,"adj.SnL"]))),
                           adj.ED=((Trophic.df[,"adj.ED"]*1-mean(na.omit(Trophic.df[,"adj.ED"])))/sd(na.omit(Trophic.df[,"adj.ED"]))),
                           adj.HL=((Trophic.df[,"adj.HL"]*1-mean(na.omit(Trophic.df[,"adj.HL"])))/sd(na.omit(Trophic.df[,"adj.HL"]))))
# makes numeric matrix for dimensionality based on eigendecomp. (by individuals, then means)
Trophic.df.euc.mat<-Trophic.df.euc[,c("adj.GW","adj.PW","RakerNum","adj.BCL",
                                      "adj.JL","adj.SnL","adj.ED","adj.HL")]
Trophic.df.euc.mat.noGR<-Trophic.df.euc[,c("adj.GW","adj.PW","adj.BCL",
                                           "adj.JL","adj.SnL","adj.ED","adj.HL")]
tro.means<-aggregate(data=na.omit(Trophic.df.euc[,c("Lake","adj.GW","adj.PW","RakerNum","adj.BCL",
                                                         "adj.JL","adj.SnL","adj.ED","adj.HL")]), .~Lake, FUN=mean)
tro.means.noGR<-aggregate(data=na.omit(Trophic.df.euc[,c("Lake","adj.GW","adj.PW","adj.BCL",
                                                              "adj.JL","adj.SnL","adj.ED","adj.HL")]), .~Lake, FUN=mean)

#LDA
Trophic.lda <- lda(Lake ~ adj.BCL + adj.GW + adj.PW + RakerNum + adj.JL +
                     adj.SnL + adj.ED + adj.HL, data = Trophic.df.euc)
Trophic.lda.noGR <- lda(Lake ~ adj.BCL + adj.GW + adj.PW  + adj.JL +
                          adj.SnL + adj.ED + adj.HL, data = Trophic.df.euc)
#         # note the following two lines are not replicated 
#           for noGR because these are not plotted 
predictions.Trophic <- predict(Trophic.lda)
LD.tro.df<-data.frame(FishID=na.omit(Trophic.df.euc)$ID,Lake=na.omit(Trophic.df.euc)$Lake, trophic.LD1=predictions.Trophic$x[,1],
                      trophic.LD2=predictions.Trophic$x[,2], trophic.LD3=predictions.Trophic$x[,3]) 

# LDA dim summary (Trophic)
(lda.tro.df<-round(lda.vector.df(Trophic.lda),3))%>%write.csv("lda-tro_df.csv")
(lda.tro.df.noGR<-lda.vector.df(Trophic.lda.noGR))%>%write.csv("lda-tro_dfnoGR.csv")
# Effective Dimensions (Trophic)
trophic.dim.matrix<-dim.stats(Trophic.lda,Trophic.df.euc.mat,tro.means)
trophic.dim.matrix.noGR<-dim.stats(Trophic.lda.noGR,Trophic.df.euc.mat.noGR,tro.means.noGR)


# Trophic LD1 frequency plot
freq.tro<-LD.freq(LD.tro.df,LD.tro.df$trophic.LD1,LD.tro.df$Lake,
        "Lake","Discriminant Score of Trophic Traits by Lake","44.5%")
# 3D scatter
threeD.scatter(LD.tro.df,"trophic.LD3",-4,c("steelblue4","firebrick3","palegreen4","palegreen2","steelblue2",
                                              "firebrick1","turquoise4","turquoise1","tan3","tan1","darkgreen",
                                              "sienna4","grey80","grey50"),c(3,4,5),c('LD1-44.5%', 'LD2-25.5%','LD3-13.2%'))
# 2D scatter 
tro.plot.12<-twoD.scatter(LD.tro.df,c(3,4),c("steelblue4","firebrick3","palegreen4","palegreen2","steelblue2",
                                              "firebrick1","turquoise4","turquoise1","tan3","tan1","darkgreen",
                                              "sienna4","grey80","grey50"),lda.tro.df,c("LD1 - ","LD2 - "),c(1,2))
tro.plot.13<-twoD.scatter(LD.tro.df,c(3,5),c("steelblue4","firebrick3","palegreen4","palegreen2","steelblue2",
                                              "firebrick1","turquoise4","turquoise1","tan3","tan1","darkgreen",
                                              "sienna4","grey80","grey50"),lda.tro.df,c("LD1 - ","LD3 - "),c(1,3))
tro.plot.12/tro.plot.13+ plot_layout(guides = 'collect')



#######  GEOMORPH  ############
######
#Read raw landmarks coordinates
AKmorph <- readland.tps("Edited_8_semi.TPS", specID="imageID")
ID <- as.factor(sapply(strsplit(dimnames(AKmorph)[[3]], "_"),"[[",3))
dimnames(AKmorph)<-list(NULL,NULL,as.vector(ID))
#slider file
slider <- as.matrix(read.table("slider_8_semi.txt"))

#Using bending energy for semi-landmark superimposition and General Procrustes Analysis (GPA)
AKmorph.sup <- gpagen(AKmorph, curves= slider, ProcD=FALSE)
#Creating geomorph dataframe#
AKmorph.df <- geomorph.data.frame(coords = AKmorph.sup$coords, Lake= data[match(ID,data$FishID),"Lake"], ID=ID)

# LDA & CVA (The same thing, but CVA results make visualizing shape change easier)
AKmorph.cva<-CVA(AKmorph.df$coords,AKmorph.df$Lake, weighting = T, plot = T, cv=T)
Geomorph.df<-data.frame(FishID=AKmorph.df$ID,Lake=AKmorph.df$Lake,AKmorph.cva$CVscores)
AKmorph.lda<-lda(AKmorph.df$Lake~two.d.array(AKmorph.df$coords))

LD.geo.df<-data.frame(FishID=Geomorph.df$FishID,Lake=Geomorph.df$Lake,CV1=Geomorph.df$CV.1,
                      CV2=Geomorph.df$CV.2, CV3=Geomorph.df$CV.3) 

# LDA dimension summary
(lda.geo.df<-lda.vector.df(AKmorph.lda))%>%write.csv("lda-geo_df.csv")
geo.means<-aggregate(data=na.omit(data.frame(two.d.array(AKmorph.df$coords),Lake=AKmorph.df$Lake)), .~Lake, FUN=mean)

# Effective Dimensions (Trophic)
geo.dim.matrix<-dim.stats(AKmorph.lda,data.frame(two.d.array(AKmorph.df$coords)),geo.means)


#Visualize individuals
plotAllSpecimens(AKmorph.sup$coords)
plotOutliers(AKmorph.sup$coords)
plotTangentSpace(AKmorph.sup$coords, groups = AKmorph.df$Lake, legend = T)
CVposfive<-5*AKmorph.cva$CVvis[,1]+AKmorph.cva$Grandm
CVnegfive<-(-5)*AKmorph.cva$CVvis[,1]+AKmorph.cva$Grandm
deformGrid2d(CVposfive,CVnegfive,ngrid = 10,wireframe = c(1,2,7:21,1),
             cex1=.75,cex2=.75,col1 = "firebrick4",col2 = "steelblue2")
CV2posfive<-5*AKmorph.cva$CVvis[,2]+AKmorph.cva$Grandm
CV2negfive<-(-5)*AKmorph.cva$CVvis[,2]+AKmorph.cva$Grandm
deformGrid2d(CV2posfive,CV2negfive,ngrid = 10,wireframe = c(1,2,7:21,1),
             cex1=.5,cex2=.5,col1 = "firebrick4",col2 = "steelblue2")

# Geomorph CV1 frequency plot
LD.freq(Geomorph.df,Geomorph.df$CV.1,Geomorph.df$Lake,
        "Lake","Shape CV 1 by Lake", "35.2%")
# 3D scatter
threeD.scatter(LD.geo.df,"CV3",-4.5,c("steelblue4","firebrick3","palegreen4","steelblue2",
                                              "firebrick1","turquoise4","turquoise1","tan3","tan1","darkgreen",
                                              "grey80","grey50"),c(3,4,5),c('LD1-35.2%', 'LD2-17.0%','LD3-11.9%'))
# 2D scatter 
geo.plot.12<-twoD.scatter(LD.geo.df,c(3,4),c("steelblue4","firebrick3","palegreen4","steelblue2",
                                              "firebrick1","turquoise4","turquoise1","tan3","tan1","darkgreen",
                                              "grey80","grey50"),lda.geo.df,c("LD1 - ","LD2 - "),c(1,2))
geo.plot.13<-twoD.scatter(LD.geo.df,c(3,5),c("steelblue4","firebrick3","palegreen4","steelblue2",
                                              "firebrick1","turquoise4","turquoise1","tan3","tan1","darkgreen",
                                              "grey80","grey50"),lda.geo.df,c("LD1 - ","LD3 - "),c(1,3))
geo.plot.12/geo.plot.13+ plot_layout(guides = 'collect')



#
#

# adj.data uses adjusted trait scaled to SD, adj.data.unscaled has unscaled adjusted data for later use in RDAs

##### CHECK TEXT FOR WHY "Pelvic_Total" excluded!!!!
adj.data<-merge(Defense.df.euc[c("ID","Lake","LatPlNum_mean","adj.DS1L","adj.DS2L","adj.PS_mean","MidPl_mean")],
                merge(Swimming.df.euc[c("ID","SL_mm","adj.BD","adj.CP")],
                      Trophic.df.euc[c("ID","adj.BCL","adj.GW","adj.EP","RakerNum","adj.JL","adj.SnL","adj.ED","adj.HL")],
                              by= "ID",all= T), by = "ID", all = T)
rownames(adj.data)<-adj.data$ID

adj.data.mat<-subset(adj.data,select = -c(ID,Lake))
adj.data.mat.noGR<-subset(adj.data,select = -c(ID,Lake,RakerNum))

full.means<-aggregate(data=na.omit(subset(adj.data,select = -ID)), .~Lake, FUN=mean)
full.means.noGR<-aggregate(data=na.omit(subset(adj.data,select = -c(ID,RakerNum))), .~Lake, FUN=mean)


#LDA
total.lda <- lda(Lake ~ LatPlNum_mean+adj.DS1L+adj.DS2L+adj.PS_mean+MidPl_mean +SL_mm+adj.BD+adj.CP+
                   adj.BCL + adj.GW + adj.PW + RakerNum + adj.JL +adj.SnL + adj.ED + adj.HL, data = adj.data)
total.lda.noGR <- lda(Lake ~ LatPlNum_mean+adj.DS1L+adj.DS2L+adj.PS_mean+MidPl_mean +SL_mm+adj.BD+adj.CP+
                        adj.BCL + adj.GW + adj.PW +adj.JL +adj.SnL + adj.ED + adj.HL, data = adj.data)
predictions.total <- predict(total.lda)
LD.total.df<-data.frame(FishID=na.omit(adj.data)$ID,Lake=na.omit(adj.data)$Lake, total.LD1=predictions.total$x[,1],
                        total.LD2=predictions.total$x[,2], total.LD3=predictions.total$x[,3]) 
predictions.noGR.total <- predict(total.lda.noGR)
LD.total.noGR.df<-data.frame(ID=na.omit(subset(adj.data,select=-RakerNum))$ID,Lake=na.omit(subset(adj.data,select=-RakerNum))$Lake, 
                             total.LD1=predictions.noGR.total$x[,1],total.LD2=predictions.noGR.total$x[,2], 
                             total.LD3=predictions.noGR.total$x[,3])

# LDA dim summary (full/fullnoGR)
(lda.full.df<-round(lda.vector.df(total.lda),3))%>%write.csv("lda-full_df.csv")
(lda.full.df.noGR<-round(lda.vector.df(total.lda.woGR),3))%>%write.csv("lda-full_dfnoGR.csv")

# Effective Dimensions (full/fullnoGR)
full.dim.matrix<-dim.stats(total.lda,adj.data.mat,full.means)
full.dim.matrix.noGR<-dim.stats(total.lda.woGR,adj.data.mat.noGR,full.means.noGR)
# Table of effective dimensions for all trait suites
do.call(rbind,list(defense=round(defense.dim.matrix,2),trophic=round(trophic.dim.matrix,2),
                   swimming=round(swimming.dim.matrix,2),Morph=round(geo.dim.matrix,2),
                   total.noGeo=round(full.dim.matrix,2), total.noGeo.noGR=round(full.dim.matrix.noGR,2)))%>%
  write.csv("dimMatrices.csv")



# Total LD1 frequency plot
freq.total<-LD.freq(LD.total.df,LD.total.df$total.LD1,LD.total.df$Lake,
        "Lake","Discriminant Score of D/Sw/T Traits by Lake","69.0%")
# Total LD 3D scatter plot
# 3D scatter
threeD.scatter(LD.total.df,"total.LD3",-5,c("steelblue4","firebrick3","palegreen4","palegreen2","steelblue2",
                                               "firebrick1","turquoise4","turquoise1","tan3","tan1","darkgreen",
                                               "sienna4","grey80","grey50"),c(3,4,5),c('LD1-69.0%', 'LD2-10.4%','LD3-7.8%'))

total.plot.12<-twoD.scatter(LD.total.df,c(3,4),c("steelblue4","firebrick3","palegreen4","palegreen2","steelblue2",
                                                "firebrick1","turquoise4","turquoise1","tan3","tan1","darkgreen",
                                                "sienna4","grey80","grey50"),lda.full.df,c("LD1 - ","LD2 - "),c(1,2))
total.plot.13<-twoD.scatter(LD.total.df,c(3,5),c("steelblue4","firebrick3","palegreen4","palegreen2","steelblue2",
                                                "firebrick1","turquoise4","turquoise1","tan3","tan1","darkgreen",
                                                "sienna4","grey80","grey50"),lda.full.df,c("LD1 - ","LD3 - "),c(1,3))
total.plot.12/total.plot.13+ plot_layout(guides = 'collect')




(def.plot.12+swi.plot.12+tro.plot.12+
    def.plot.13+swi.plot.13+tro.plot.13)+ 
  plot_layout(guides = 'collect',nrow=2,widths=c(1,1,1))
geo.plot.12/geo.plot.13+ plot_layout(nrow=2,guides = 'collect')

total.plot.12/total.plot.13+ plot_layout(guides = 'collect')



# Total LD1 frequency plot w/o GR
freq.total.noGR<-LD.freq(LD.total.woGR.df,LD.total.woGR.df$total.LD1,LD.total.woGR.df$Lake,
        "Lake","Discriminant Score of D/S/T Traits by Lake","71.4%")


all.LD<-(merge(merge(LD.swim.df[c("FishID","Lake","swimming.LD1")],LD.def.df[c("FishID","defense.LD1")],by="FishID",all = TRUE, no.dups = T),
              LD.tro.df[c("FishID","trophic.LD1")],by.x ="FishID",by.y="FishID", all = TRUE, no.dups = T)%>%
           merge(LD.geo.df[c("CV1","FishID")],by.x="FishID",by.y = "FishID",all = TRUE))
rownames(all.LD)<-all.LD$FishID


######
(freq.def+freq.swi+freq.tro)/freq.total
