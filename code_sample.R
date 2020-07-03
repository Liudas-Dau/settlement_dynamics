
#First round of variable importance evaluation, FULL DATA, Neolithic
#Importing,  Extracting and Formatting relevant data
{
  REG_DATA <- read.table("C://REG_DATA_TRU.txt", 
                         header = TRUE, sep = ";",dec = ",")
  #Information about variables contained in REG_DATA:
  # column 1, "eolas" - a logical indicator, whether a point is overlying a soil of eolian origin (1: TRU, 0: FALSE)
  # column 2, "litologija" - lithology of soil (coded)
  # column 3, "Lok" - classification of settlement localization accuracy (coded PROLIGIS dataset values); 22 #- non-site point
  # column 4, "Lok_adjust" - adjusted classification of settlement localization accuracy according to the #Cultural Heritage Department data; 11 - adjusted localization
  # column 5-6, "POINT_X", "POINT_Y" - coordinates of points. The projection is LKS94 / Lithuania TM; crs: #3346
  # column 7-30, "PAL"..."VIVEGA_spec" - logical indicator, whether site is present in a given dataset #(SPEC/FULL), during a given period. 
  #PAL - Paleolithic, MEZO - Mesolithic, NEO - Neolithic, BRO - Bronze Age, ASEGA - Early and Old Iron Age #(EOIA), VIVEGA - Middle and Late Iron Age (MLIA)
  #variable without _spec at the end belongs to FULL dataset, with _spec - to SPEC dataset.
  # column 31-51, "pc7" ... "reg_wetnes" - regional environmental variables used in the study.
  #pc1-pc7 - principal components of lithology, reg_avht - EACN, topo_vid_d - WAT_DIST(t), vand_top_p - #WAT_AREA(t), r_wdi_geo_ - WAT_DIST(g),
  #prop_van_g - WAT_AREA(g), dist_krant - SEA_DIST, s_litho - SOIL_SP, litho_div - SOIL_DIV, slait_r10k - #SLOPE, titnago_at - FLINT,
  #reg_topex_ - TOPEX, tri_reg_r1 - TRI, reg_wetnes - WET, reg_ls_eu - LS.
  mode(attributes(REG_DATA)$row.names)
  rownames(REG_DATA)<-as.character(seq(nrow(REG_DATA)))
  #settlements extraction
  neo_dat_full<-REG_DATA[,c(names(REG_DATA)[1:6],"neo","neo_spec",names(REG_DATA)[31:51])]
  neo_set<-neo_dat_full[which(neo_dat_full[,3]!=22&neo_dat_full[,"neo"]==1),] # Neolithic sites only
  NET<-neo_dat_full[which(neo_dat_full[,3]==22),] # non-sites
  dat_neo<-rbind(neo_set,NET) # sites & non-sites combained
  # training data
  train_neo_full<-dat_neo[,-c(1:4,8)] # leaving only coordinates, dependent variable (status of Neolithic sites) and predictor variables 
  #
  rownames(train_neo_full)<-as.character(seq(nrow(train_neo_full)))
  #
  attach(train_neo_full)
  train_neo_full[,3]<-factor(neo)
}
## Simulating spatially autocorrelated data
{
  data_sp<-train_neo_full
  library(gstat)
  library(sp)
  coordinates(data_sp)=~POINT_X+POINT_Y
  xy<-data.frame(POINT_X,POINT_Y)
  names(xy) <- c('x','y')
  #1
  a<-variogram(slait_r10k~POINT_X+POINT_Y, data_sp)
  b<-fit.variogram(a, vgm("Sph"))
  g.dummy <- gstat(formula=z~x+y, locations=~x+y, dummy=T, beta=1, model=vgm(psill=b[2,2], range=b[2,3], model='Exp'), nmax=20)
  sim_slope <- predict(g.dummy, newdata=xy, nsim=30)[,c(3:32)]
  #2
  a<-variogram(tri_reg_r1~POINT_X+POINT_Y, data_sp)
  b<-fit.variogram(a, vgm("Sph"))
  g.dummy <- gstat(formula=z~x+y, locations=~x+y, dummy=T, beta=1, model=vgm(psill=b[2,2], range=b[2,3], model='Exp'), nmax=20)
  sim_TRI <- predict(g.dummy, newdata=xy, nsim=30)[,c(3:32)]
  #3
  a<-variogram(dist_krant~POINT_X+POINT_Y, data_sp)
  b<-fit.variogram(a, vgm("Sph"))
  g.dummy <- gstat(formula=z~x+y, locations=~x+y, dummy=T, beta=1, model=vgm(psill=b[2,2], range=b[2,3], model='Exp'), nmax=20)
  sim_dist_krant <- predict(g.dummy, newdata=xy, nsim=30)[,c(3:32)]
  #4
  a<-variogram(prop_van_g~POINT_X+POINT_Y, data_sp)
  b<-fit.variogram(a, vgm("Sph"))
  g.dummy <- gstat(formula=z~x+y, locations=~x+y, dummy=T, beta=1, model=vgm(psill=b[2,2], range=b[2,3], model='Exp'), nmax=20)
  sim_vandpro_Geo <- predict(g.dummy, newdata=xy, nsim=30)[,c(3:32)]
  #5
  a<-variogram(reg_avht~POINT_X+POINT_Y, data_sp)
  b<-fit.variogram(a, vgm("Sph"))
  g.dummy <- gstat(formula=z~x+y, locations=~x+y, dummy=T, beta=1, model=vgm(psill=b[2,2], range=b[2,3], model='Exp'), nmax=20)
  sim_AVHT <- predict(g.dummy, newdata=xy, nsim=30)[,c(3:32)]
  #6
  a<-variogram(vand_top_p~POINT_X+POINT_Y, data_sp)
  b<-fit.variogram(a, vgm("Sph"))
  g.dummy <- gstat(formula=z~x+y, locations=~x+y, dummy=T, beta=1, model=vgm(psill=b[2,2], range=b[2,3], model='Exp'), nmax=20)
  sim_vandpro_Topo <- predict(g.dummy, newdata=xy, nsim=30)[,c(3:32)]
  #7
  a<-variogram(topo_vid_d~POINT_X+POINT_Y, data_sp)
  b<-fit.variogram(a, vgm("Sph"))
  g.dummy <- gstat(formula=z~x+y, locations=~x+y, dummy=T, beta=1, model=vgm(psill=b[2,2], range=b[2,3], model='Exp'), nmax=20)
  sim_vid_v_Topo <- predict(g.dummy, newdata=xy, nsim=30)[,c(3:32)]
  #8
  a<-variogram(litho_div~POINT_X+POINT_Y, data_sp)
  b<-fit.variogram(a, vgm("Sph"))
  g.dummy <- gstat(formula=z~x+y, locations=~x+y, dummy=T, beta=1, model=vgm(psill=b[2,2], range=b[2,3], model='Exp'), nmax=20)
  sim_litho_div <- predict(g.dummy, newdata=xy, nsim=30)[,c(3:32)]
  #9
  a<-variogram(pc1~POINT_X+POINT_Y, data_sp)
  b<-fit.variogram(a, vgm("Sph"))
  g.dummy <- gstat(formula=z~x+y, locations=~x+y, dummy=T, beta=1, model=vgm(psill=b[2,2], range=b[2,3], model='Exp'), nmax=20)
  sim_pc1 <- predict(g.dummy, newdata=xy, nsim=30)[,c(3:32)]
  #10
  a<-variogram(pc2~POINT_X+POINT_Y, data_sp)
  b<-fit.variogram(a, vgm("Sph"))
  g.dummy <- gstat(formula=z~x+y, locations=~x+y, dummy=T, beta=1, model=vgm(psill=b[2,2], range=b[2,3], model='Exp'), nmax=20)
  sim_pc2 <- predict(g.dummy, newdata=xy, nsim=30)[,c(3:32)]
  #11
  a<-variogram(pc3~POINT_X+POINT_Y, data_sp)
  b<-fit.variogram(a, vgm("Sph"))
  g.dummy <- gstat(formula=z~x+y, locations=~x+y, dummy=T, beta=1, model=vgm(psill=b[2,2], range=b[2,3], model='Exp'), nmax=20)
  sim_pc3 <- predict(g.dummy, newdata=xy, nsim=30)[,c(3:32)]
  #12
  a<-variogram(pc4~POINT_X+POINT_Y, data_sp)
  b<-fit.variogram(a, vgm("Sph"))
  g.dummy <- gstat(formula=z~x+y, locations=~x+y, dummy=T, beta=1, model=vgm(psill=b[2,2], range=b[2,3], model='Exp'), nmax=20)
  sim_pc4 <- predict(g.dummy, newdata=xy, nsim=30)[,c(3:32)]
  #13
  a<-variogram(pc5~POINT_X+POINT_Y, data_sp)
  b<-fit.variogram(a, vgm("Sph"))
  g.dummy <- gstat(formula=z~x+y, locations=~x+y, dummy=T, beta=1, model=vgm(psill=b[2,2], range=b[2,3], model='Exp'), nmax=20)
  sim_pc5 <- predict(g.dummy, newdata=xy, nsim=30)[,c(3:32)]
  #14
  a<-variogram(pc6~POINT_X+POINT_Y, data_sp)
  b<-fit.variogram(a, vgm("Sph"))
  g.dummy <- gstat(formula=z~x+y, locations=~x+y, dummy=T, beta=1, model=vgm(psill=b[2,2], range=b[2,3], model='Exp'), nmax=20)
  sim_pc6 <- predict(g.dummy, newdata=xy, nsim=30)[,c(3:32)]
  #15
  a<-variogram(pc7~POINT_X+POINT_Y, data_sp)
  b<-fit.variogram(a, vgm("Sph"))
  g.dummy <- gstat(formula=z~x+y, locations=~x+y, dummy=T, beta=1, model=vgm(psill=b[2,2], range=b[2,3], model='Exp'), nmax=20)
  sim_pc7 <- predict(g.dummy, newdata=xy, nsim=30)[,c(3:32)]
  #16
  a<-variogram(r_wdi_geo_~POINT_X+POINT_Y, data_sp)
  b<-fit.variogram(a, vgm("Sph"))
  g.dummy <- gstat(formula=z~x+y, locations=~x+y, dummy=T, beta=1, model=vgm(psill=b[2,2], range=b[2,3], model='Exp'), nmax=20)
  sim_vid_v_GEO <- predict(g.dummy, newdata=xy, nsim=30)[,c(3:32)]
  #17
  a<-variogram(s_litho~POINT_X+POINT_Y, data_sp)
  b<-fit.variogram(a, vgm("Sph"))
  g.dummy <- gstat(formula=z~x+y, locations=~x+y, dummy=T, beta=1, model=vgm(psill=b[2,2], range=b[2,3], model='Exp'), nmax=20)
  sim_litho_spR <- predict(g.dummy, newdata=xy, nsim=30)[,c(3:32)]
  for (i in seq(30)){
    sim_litho_spR[,i]<-17/(max(sim_litho_spR[,i])-min(sim_litho_spR[,i]))*(sim_litho_spR[,i]-min(sim_litho_spR[,i]))+1
  }
  sim_litho_spR<-round(sim_litho_spR,0)
  #18
  a<-variogram(reg_topex_~POINT_X+POINT_Y, data_sp)
  b<-fit.variogram(a, vgm("Sph"))
  g.dummy <- gstat(formula=z~x+y, locations=~x+y, dummy=T, beta=1, model=vgm(psill=b[2,2], range=b[2,3], model='Exp'), nmax=20)
  sim_topex<- predict(g.dummy, newdata=xy, nsim=30)[,c(3:32)]
  #19
  a<-variogram(reg_ls_eu~POINT_X+POINT_Y, data_sp)
  b<-fit.variogram(a, vgm("Sph"))
  g.dummy <- gstat(formula=z~x+y, locations=~x+y, dummy=T, beta=1, model=vgm(psill=b[2,2], range=b[2,3], model='Exp'), nmax=20)
  sim_ls <- predict(g.dummy, newdata=xy, nsim=30)[,c(3:32)]
  #20
  a<-variogram(reg_wetnes~POINT_X+POINT_Y, data_sp)
  b<-fit.variogram(a, vgm("Sph"))
  g.dummy <- gstat(formula=z~x+y, locations=~x+y, dummy=T, beta=1, model=vgm(psill=b[2,2], range=b[2,3], model='Exp'), nmax=20)
  sim_Wet <- predict(g.dummy, newdata=xy, nsim=30)[,c(3:32)]
  #
  kasyklos<-xy[round(runif(30,1,nrow(dat_neo)),0),]
  kasyklos1<-xy[round(runif(30,1,nrow(dat_neo)),0),]
  kasyklos2<-xy[round(runif(30,1,nrow(dat_neo)),0),]
  
  di1<-list()
  di2<-list()
  di3<-list()
  
  a<-rep(NA,nrow(dat_neo))
  sim_titnagas<-data.frame(a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a)
  
  for (i in 1:30) {
    di1[[i]]<-sqrt((xy[,1]-kasyklos[i,1])^2+(xy[,2]-kasyklos[i,2])^2)
    di2[[i]]<-sqrt((xy[,1]-kasyklos1[i,1])^2+(xy[,2]-kasyklos1[i,2])^2)
    di3[[i]]<-sqrt((xy[,1]-kasyklos2[i,1])^2+(xy[,2]-kasyklos2[i,2])^2)
    for (b in 1:nrow(dat_neo)) {
      sim_titnagas[,i][b]<-min(data.frame(di1[[i]],di2[[i]],di3[[i]])[b,])
    }
    print(i)
  }
}
detach(train_neo_full)
#loading neccesery packages
library("randomForest")
library(caret)
library(ROCR)
####
data<-train_neo_full[,-c(1:2)] #removing coordinates

lcp<-42 # number of variables + number of simulations
repeats<-30 
# creating empty data structures
train<-vector(mode="list", repeats)

Gini1<-vector(mode="list", repeats)
Gini2<-vector(mode="list", repeats)
Gini3<-vector(mode="list", repeats)
Gini4<-vector(mode="list", repeats)
Gini5<-vector(mode="list", repeats)
Gini6<-vector(mode="list", repeats)
Gini7<-vector(mode="list", repeats)
Gini8<-vector(mode="list", repeats)
Gini9<-vector(mode="list", repeats)
Gini10<-vector(mode="list", repeats)

Acc1<-vector(mode="list", repeats)
Acc2<-vector(mode="list", repeats)
Acc3<-vector(mode="list", repeats)
Acc4<-vector(mode="list", repeats)
Acc5<-vector(mode="list", repeats)
Acc6<-vector(mode="list", repeats)
Acc7<-vector(mode="list", repeats)
Acc8<-vector(mode="list", repeats)
Acc9<-vector(mode="list", repeats)
Acc10<-vector(mode="list", repeats)

AUC1<-vector(mode="numeric", repeats)
AUC2<-vector(mode="numeric", repeats)
AUC3<-vector(mode="numeric", repeats)
AUC4<-vector(mode="numeric", repeats)
AUC5<-vector(mode="numeric", repeats)
AUC6<-vector(mode="numeric", repeats)
AUC7<-vector(mode="numeric", repeats)
AUC8<-vector(mode="numeric", repeats)
AUC9<-vector(mode="numeric", repeats)
AUC10<-vector(mode="numeric", repeats)
#performing repeated (x30) stratified 10-fold cross-validation
for (a in seq(repeats)){
  Acc1[[a]] <- NaN*seq(lcp)
  Acc2[[a]] <- NaN*seq(lcp)
  Acc3[[a]] <- NaN*seq(lcp)
  Acc4[[a]] <- NaN*seq(lcp)
  Acc5[[a]] <- NaN*seq(lcp)
  Acc6[[a]] <- NaN*seq(lcp)
  Acc7[[a]] <- NaN*seq(lcp)
  Acc8[[a]] <- NaN*seq(lcp)
  Acc9[[a]] <- NaN*seq(lcp)
  Acc10[[a]] <- NaN*seq(lcp)
  
  Gini1[[a]] <- NaN*seq(lcp)
  Gini2[[a]] <- NaN*seq(lcp)
  Gini3[[a]] <- NaN*seq(lcp)
  Gini4[[a]] <- NaN*seq(lcp)
  Gini5[[a]] <- NaN*seq(lcp)
  Gini6[[a]] <- NaN*seq(lcp)
  Gini7[[a]] <- NaN*seq(lcp)
  Gini8[[a]] <- NaN*seq(lcp)
  Gini9[[a]] <- NaN*seq(lcp)
  Gini10[[a]] <- NaN*seq(lcp)
  
  ### combining data with a set of simulations
  train[[a]]<-data.frame(data,sim_pc2[,a],sim_AVHT[,a],sim_pc3[,a],sim_dist_krant[,a],
                         sim_litho_div[,a],sim_litho_spR[,a],sim_ls[,a],
                         sim_pc4[,a],sim_pc5[,a],sim_pc6[,a],sim_pc7[,a],
                         sim_slope[,a],sim_pc1[,a],sim_titnagas[,a],sim_topex[,a],
                         sim_TRI[,a],sim_vandpro_Geo[,a],sim_vandpro_Topo[,a],
                         sim_vid_v_GEO[,a],sim_vid_v_Topo[,a],sim_Wet[,a])
  # creating stratified folds & sampling data
  folds <- createFolds(data$neo, k = 10, list = FALSE,returnTrain = FALSE)
  #data used for testing
  fold1<-train[[a]][which(folds==1),]
  fold2<-train[[a]][which(folds==2),]
  fold3<-train[[a]][which(folds==3),]
  fold4<-train[[a]][which(folds==4),]
  fold5<-train[[a]][which(folds==5),]
  fold6<-train[[a]][which(folds==6),]
  fold7<-train[[a]][which(folds==7),]
  fold8<-train[[a]][which(folds==8),]
  fold9<-train[[a]][which(folds==9),]
  fold10<-train[[a]][which(folds==10),]
  #data used to train models
  train1<-train[[a]][-which(folds==1),]
  train2<-train[[a]][-which(folds==2),]
  train3<-train[[a]][-which(folds==3),]
  train4<-train[[a]][-which(folds==4),]
  train5<-train[[a]][-which(folds==5),]
  train6<-train[[a]][-which(folds==6),]
  train7<-train[[a]][-which(folds==7),]
  train8<-train[[a]][-which(folds==8),]
  train9<-train[[a]][-which(folds==9),]
  train10<-train[[a]][-which(folds==10),]
  
  ### kryzminimas
  
  ### 1 
  miskas<-randomForest(x=train1[,-1],y=train1[,1],ntree=500,
                       mtry=1,replace=TRUE, strata= train1[,1],
                       sampsize = c(round(length(train1[,1][which(train1[,1]==1)])*0.75,0),
                                    round(length(train1[,1][which(train1[,1]==1)])*0.75,0)),nodesize=2,
                       importance = TRUE, norm.votes = TRUE, keep.forest = TRUE)
  # Importance
  Gini1[[a]]<-importance(miskas)[,4]
  Acc1[[a]]<-(importance(miskas,scale = F)[,1]+importance(miskas,scale = F)[,2])/2
  # AUROC calculation
  pred1<-predict(miskas,fold1,type = "prob")
  predi1 <- prediction(pred1[,2],fold1[,1])
  perf1 <- performance(predi1, measure = "auc")
  AUC1[a]<-as.numeric(perf1@y.values)
  
  print((((a-1)*10+1)/300)*100)
  
  ### 2  
  miskas<-randomForest(x=train2[,-1],y=train2[,1],ntree=500,
                       mtry=1,replace=TRUE, strata= train2[,1],
                       sampsize = c(round(length(train2[,1][which(train2[,1]==1)])*0.75,0),
                                    round(length(train2[,1][which(train2[,1]==1)])*0.75,0)),nodesize=2,
                       importance = TRUE, norm.votes = TRUE, keep.forest = TRUE)
  # Importance
  Gini2[[a]]<-importance(miskas)[,4]
  Acc2[[a]]<-(importance(miskas,scale = F)[,1]+importance(miskas,scale = F)[,2])/2
  # AUROC calculation
  pred1<-predict(miskas,fold2,type = "prob")
  predi1 <- prediction(pred1[,2],fold2[,1])
  perf1 <- performance(predi1, measure = "auc")
  AUC2[a]<-as.numeric(perf1@y.values)
  
  print((((a-1)*10+2)/300)*100)
  
  ### 3  
  miskas<-randomForest(x=train3[,-1],y=train3[,1],ntree=500,
                       mtry=1,replace=TRUE, strata= train3[,1],
                       sampsize = c(round(length(train3[,1][which(train3[,1]==1)])*0.75,0),
                                    round(length(train3[,1][which(train3[,1]==1)])*0.75,0)),nodesize=2,
                       importance = TRUE, norm.votes = TRUE, keep.forest = TRUE)
  # Importance
  Gini3[[a]]<-importance(miskas)[,4]
  Acc3[[a]]<-(importance(miskas,scale = F)[,1]+importance(miskas,scale = F)[,2])/2
  # AUROC calculation
  pred1<-predict(miskas,fold3,type = "prob")
  predi1 <- prediction(pred1[,2],fold3[,1])
  perf1 <- performance(predi1, measure = "auc")
  AUC3[a]<-as.numeric(perf1@y.values)
  
  print((((a-1)*10+3)/300)*100)
  
  ### 4  
  miskas<-randomForest(x=train4[,-1],y=train4[,1],ntree=500,
                       mtry=1,replace=TRUE, strata= train4[,1],
                       sampsize = c(round(length(train4[,1][which(train4[,1]==1)])*0.75,0),
                                    round(length(train4[,1][which(train4[,1]==1)])*0.75,0)),nodesize=2,
                       importance = TRUE, norm.votes = TRUE, keep.forest = TRUE)
  # Importance
  Gini4[[a]]<-importance(miskas)[,4]
  Acc4[[a]]<-(importance(miskas,scale = F)[,1]+importance(miskas,scale = F)[,2])/2
  # AUROC calculation
  pred1<-predict(miskas,fold4,type = "prob")
  predi1 <- prediction(pred1[,2],fold4[,1])
  perf1 <- performance(predi1, measure = "auc")
  AUC4[a]<-as.numeric(perf1@y.values)
  
  print((((a-1)*10+4)/300)*100)
  
  ### 5  
  miskas<-randomForest(x=train5[,-1],y=train5[,1],ntree=500,
                       mtry=1,replace=TRUE, strata= train5[,1],
                       sampsize = c(round(length(train5[,1][which(train5[,1]==1)])*0.75,0),
                                    round(length(train5[,1][which(train5[,1]==1)])*0.75,0)),nodesize=2,
                       importance = TRUE, norm.votes = TRUE, keep.forest = TRUE)
  # Importance
  Gini5[[a]]<-importance(miskas)[,4]
  Acc5[[a]]<-(importance(miskas,scale = F)[,1]+importance(miskas,scale = F)[,2])/2
  # AUROC calculation
  pred1<-predict(miskas,fold5,type = "prob")
  predi1 <- prediction(pred1[,2],fold5[,1])
  perf1 <- performance(predi1, measure = "auc")
  AUC5[a]<-as.numeric(perf1@y.values)
  
  print((((a-1)*10+5)/300)*100)
  
  ### 6  
  miskas<-randomForest(x=train6[,-1],y=train6[,1],ntree=500,
                       mtry=1,replace=TRUE, strata= train6[,1],
                       sampsize = c(round(length(train6[,1][which(train6[,1]==1)])*0.75,0),
                                    round(length(train6[,1][which(train6[,1]==1)])*0.75,0)),nodesize=2,
                       importance = TRUE, norm.votes = TRUE, keep.forest = TRUE)
  # Importance
  Gini6[[a]]<-importance(miskas)[,4]
  Acc6[[a]]<-(importance(miskas,scale = F)[,1]+importance(miskas,scale = F)[,2])/2
  # AUROC calculation
  pred1<-predict(miskas,fold6,type = "prob")
  predi1 <- prediction(pred1[,2],fold6[,1])
  perf1 <- performance(predi1, measure = "auc")
  AUC6[a]<-as.numeric(perf1@y.values)
  
  print((((a-1)*10+6)/300)*100)
  
  ### 7  
  miskas<-randomForest(x=train7[,-1],y=train7[,1],ntree=500,
                       mtry=1,replace=TRUE, strata= train7[,1],
                       sampsize = c(round(length(train7[,1][which(train7[,1]==1)])*0.75,0),
                                    round(length(train7[,1][which(train7[,1]==1)])*0.75,0)),nodesize=2,
                       importance = TRUE, norm.votes = TRUE, keep.forest = TRUE)
  # Importance
  Gini7[[a]]<-importance(miskas)[,4]
  Acc7[[a]]<-(importance(miskas,scale = F)[,1]+importance(miskas,scale = F)[,2])/2
  # AUROC calculation
  pred1<-predict(miskas,fold7,type = "prob")
  predi1 <- prediction(pred1[,2],fold7[,1])
  perf1 <- performance(predi1, measure = "auc")
  AUC7[a]<-as.numeric(perf1@y.values)
  
  print((((a-1)*10+7)/300)*100)
  
  ### 8  
  miskas<-randomForest(x=train8[,-1],y=train8[,1],ntree=500,
                       mtry=1,replace=TRUE, strata= train8[,1],
                       sampsize = c(round(length(train8[,1][which(train8[,1]==1)])*0.75,0),
                                    round(length(train8[,1][which(train8[,1]==1)])*0.75,0)),nodesize=2,
                       importance = TRUE, norm.votes = TRUE, keep.forest = TRUE)
  # Importance
  Gini8[[a]]<-importance(miskas)[,4]
  Acc8[[a]]<-(importance(miskas,scale = F)[,1]+importance(miskas,scale = F)[,2])/2
  # AUROC calculation
  pred1<-predict(miskas,fold8,type = "prob")
  predi1 <- prediction(pred1[,2],fold8[,1])
  perf1 <- performance(predi1, measure = "auc")
  AUC8[a]<-as.numeric(perf1@y.values)
  
  print((((a-1)*10+8)/300)*100)
  
  ### 9  
  miskas<-randomForest(x=train9[,-1],y=train9[,1],ntree=500,
                       mtry=1,replace=TRUE, strata= train9[,1],
                       sampsize = c(round(length(train9[,1][which(train9[,1]==1)])*0.75,0),
                                    round(length(train9[,1][which(train9[,1]==1)])*0.75,0)),nodesize=2,
                       importance = TRUE, norm.votes = TRUE, keep.forest = TRUE)
  # Importance
  Gini9[[a]]<-importance(miskas)[,4]
  Acc9[[a]]<-(importance(miskas,scale = F)[,1]+importance(miskas,scale = F)[,2])/2
  # AUROC calculation
  pred1<-predict(miskas,fold9,type = "prob")
  predi1 <- prediction(pred1[,2],fold9[,1])
  perf1 <- performance(predi1, measure = "auc")
  AUC9[a]<-as.numeric(perf1@y.values)
  
  print((((a-1)*10+9)/300)*100)
  
  ### 10  
  miskas<-randomForest(x=train10[,-1],y=train10[,1],ntree=500,
                       mtry=1,replace=TRUE, strata= train10[,1],
                       sampsize = c(round(length(train10[,1][which(train10[,1]==1)])*0.75,0),
                                    round(length(train10[,1][which(train10[,1]==1)])*0.75,0)),nodesize=2,
                       importance = TRUE, norm.votes = TRUE, keep.forest = TRUE)
  # Importance
  Gini10[[a]]<-importance(miskas)[,4]
  Acc10[[a]]<-(importance(miskas,scale = F)[,1]+importance(miskas,scale = F)[,2])/2
  # AUROC calculation
  pred1<-predict(miskas,fold10,type = "prob")
  predi1 <- prediction(pred1[,2],fold10[,1])
  perf1 <- performance(predi1, measure = "auc")
  AUC10[a]<-as.numeric(perf1@y.values)
  
  print((((a-1)*10+10)/300)*100)
}
### proccesing, inspecting and obtaining relevant results
{
  Gini_neo_full<-rbind.data.frame(t(do.call(cbind.data.frame,Gini1)),t(do.call(cbind.data.frame,Gini2)),
                                  t(do.call(cbind.data.frame,Gini3)),t(do.call(cbind.data.frame,Gini4)),
                                  t(do.call(cbind.data.frame,Gini5)),t(do.call(cbind.data.frame,Gini6)),
                                  t(do.call(cbind.data.frame,Gini7)),t(do.call(cbind.data.frame,Gini8)),
                                  t(do.call(cbind.data.frame,Gini9)),t(do.call(cbind.data.frame,Gini10)))
  Acc_neo_full<-rbind.data.frame(t(do.call(cbind.data.frame,Acc1)),t(do.call(cbind.data.frame,Acc2)),
                                 t(do.call(cbind.data.frame,Acc3)),t(do.call(cbind.data.frame,Acc4)),
                                 t(do.call(cbind.data.frame,Acc5)),t(do.call(cbind.data.frame,Acc6)),
                                 t(do.call(cbind.data.frame,Acc7)),t(do.call(cbind.data.frame,Acc8)),
                                 t(do.call(cbind.data.frame,Acc9)),t(do.call(cbind.data.frame,Acc10)))
  ### plots 
  boxplot(Gini_neo_full,las=2, main="GINI index")
  boxplot(Acc_neo_full,las=2, main="Accuracy reduction after permutation")
  ### collecting relevant data
  Acc_neo_full_sim<-Acc_neo_full[,c(c(11,10,9,8,3,1,13,2,20,18,19,17,4,6,5,12,14,15,16,21,7)+21)]
  Acc_neo_full_org<-Acc_neo_full[,1:21]
  Acc_neo_full_mean_org<-apply(Acc_neo_full_org,2,mean)
  Acc_neo_full_mean_sim<-apply(Acc_neo_full_sim,2,mean)
  Acc_neo_full_uplim_all<-max(apply(Acc_neo_full[,1:42],2,mean)+apply(Acc_neo_full[,1:42],2,sd))
  Acc_neo_full_downlim_all<-min(apply(Acc_neo_full[,1:42],2,mean)-apply(Acc_neo_full[,1:42],2,sd))
  Acc_neo_full_uplim_sim<-Acc_neo_full_mean_sim+apply(Acc_neo_full_sim,2,sd)
  Acc_neo_full_downlim_sim<-Acc_neo_full_mean_sim-apply(Acc_neo_full_sim,2,sd)
  Acc_neo_full_uplim_org<-Acc_neo_full_mean_org+apply(Acc_neo_full_org,2,sd)
  Acc_neo_full_downlim_org<-Acc_neo_full_mean_org-apply(Acc_neo_full_org,2,sd)
  ###plots
  plot(1:21,Acc_neo_full_mean_org, ylim =
         c(Acc_neo_full_downlim_all,Acc_neo_full_uplim_all),xaxt='n',xlab="",
       type = "b", main= "ACCURACY DECREASE \n AFTER x VARIABLE PERMUTATION")
  polygon(c(1:21,rev(1:21)), c(Acc_neo_full_downlim_org,
                               rev(Acc_neo_full_uplim_org)), col = "grey75", border = T)
  abline(h=0)
  points(1:21, Acc_neo_full_mean_org)
  lines(1:21, Acc_neo_full_mean_org, lwd = 2)
  axis(1, at=1:21, labels=names(Acc_neo_full_mean_org), lty=2, col=1, las=2) 
  lines(1:21,Acc_neo_full_uplim_sim,type = "b",col=2)
  lines(1:21,Acc_neo_full_downlim_sim,type = "b",col=2)
  lines(1:21,Acc_neo_full_mean_sim,type = "b",col=3)
  ### selecting significantly important variables according to permutation importance
  
  Acc_var_logic_full_neo<-numeric()
  for (i in 1:21) {
    Acc_var_logic_full_neo[i]<-ecdf(Acc_neo_full_sim[,i])(mean(Acc_neo_full[,i]))
  }
  Acc_var_ind_full_neo<-which(Acc_var_logic_full_neo>=0.875) # variable selection criteria
  names(Acc_var_logic_full_neo)<-names(Acc_neo_full_org)
  neo_full_Vars_Acc<-names(Acc_neo_full_org)[Acc_var_ind_full_neo]
  ###
  #"pc7"        "pc2"        "topo_vid_d" "vand_top_p" "dist_krant" "slait_r10k"
  # "tri_reg_r1" "reg_wetnes" "reg_ls_eu
  ###
  
  #the same operations repeated for GINI Importance estimates:
  Gini_neo_full_sim<-Gini_neo_full[,c(c(11,10,9,8,3,1,13,2,20,18,19,17,4,6,5,12,14,15,16,21,7)+21)]
  Gini_neo_full_org<-Gini_neo_full[,1:21]
  Gini_neo_full_mean_org<-apply(Gini_neo_full_org,2,mean)
  Gini_neo_full_mean_sim<-apply(Gini_neo_full_sim,2,mean)
  Gini_neo_full_uplim_all<-max(apply(Gini_neo_full[,1:42],2,mean)+apply(Gini_neo_full[,1:42],2,sd))
  Gini_neo_full_downlim_all<-min(apply(Gini_neo_full[,1:42],2,mean)-apply(Gini_neo_full[,1:42],2,sd))
  Gini_neo_full_uplim_sim<-Gini_neo_full_mean_sim+apply(Gini_neo_full_sim,2,sd)
  Gini_neo_full_downlim_sim<-Gini_neo_full_mean_sim-apply(Gini_neo_full_sim,2,sd)
  Gini_neo_full_uplim_org<-Gini_neo_full_mean_org+apply(Gini_neo_full_org,2,sd)
  Gini_neo_full_downlim_org<-Gini_neo_full_mean_org-apply(Gini_neo_full_org,2,sd)
  ###
  plot(1:21,Gini_neo_full_mean_org, ylim =
         c(Gini_neo_full_downlim_all,Gini_neo_full_uplim_all),xaxt='n',xlab="",
       type = "b", main= "GINI INDEX \n AFTER x VARIABLE PERMUTATION")
  polygon(c(1:21,rev(1:21)), c(Gini_neo_full_downlim_org,
                               rev(Gini_neo_full_uplim_org)), col = "grey75", border = T)
  abline(h=0)
  points(1:21, Gini_neo_full_mean_org)
  lines(1:21, Gini_neo_full_mean_org, lwd = 2)
  axis(1, at=1:21, labels=names(Gini_neo_full_mean_org), lty=2, col=1, las=2) 
  lines(1:21,Gini_neo_full_uplim_sim,type = "b",col=2)
  lines(1:21,Gini_neo_full_downlim_sim,type = "b",col=2)
  lines(1:21,Gini_neo_full_mean_sim,type = "b",col=3)
  ###
  Gini_var_logic_full_neo<-numeric()
  for (i in 1:21) {
    Gini_var_logic_full_neo[i]<-ecdf(Gini_neo_full_sim[,i])(mean(Gini_neo_full[,i]))
  }
  Gini_var_ind_full_neo<-which(Gini_var_logic_full_neo>=0.875) # variable selection criteria
  names(Gini_var_logic_full_neo)<-names(Gini_neo_full_org)
  neo_full_Vars_Gini<-names(Gini_neo_full_org)[Gini_var_ind_full_neo]
  ###
  # "pc2"        "topo_vid_d" "vand_top_p" "dist_krant" "slait_r10k" "tri_reg_r1"
  # "reg_wetnes" "reg_ls_eu"
  ###
}
#combining names of significant variables and deleting duplicates
{
  neo_full_Vars<-c(neo_full_Vars_Acc,neo_full_Vars_Gini)
  neo_full_Vars<-neo_full_Vars[!duplicated(neo_full_Vars)]
}
### Final result of first round of variable importance estimation
neo_full_Vars<-c("pc7","pc2","topo_vid_d","vand_top_p","dist_krant","slait_r10k",
                 "tri_reg_r1","reg_wetnes","reg_ls_eu")
###
# Second round of variable importance estimation
# training data
{
  train_neo_set1<-dat_neo[,c("POINT_X","POINT_Y", "neo",neo_full_Vars)]
  #
  rownames(train_neo_set1)<-as.character(seq(nrow(train_neo_set1)))
  #
}

##Simulating spatially autocorrelated data
{
  attach(train_neo_set1)
  train_neo_set1[,"neo"]<-factor(neo)
  data_sp<-train_neo_set1
  library(gstat)
  library(sp)
  coordinates(data_sp)=~POINT_X+POINT_Y
  xy<-data.frame(POINT_X,POINT_Y)
  names(xy) <- c('x','y')
  #1
  a<-variogram(slait_r10k~POINT_X+POINT_Y, data_sp)
  b<-fit.variogram(a, vgm("Sph"))
  g.dummy <- gstat(formula=z~x+y, locations=~x+y, dummy=T, beta=1, model=vgm(psill=b[2,2], range=b[2,3], model='Exp'), nmax=20)
  sim_slait_r10k <- predict(g.dummy, newdata=xy, nsim=30)[,c(3:32)]
  #2
  a<-variogram(tri_reg_r1~POINT_X+POINT_Y, data_sp)
  b<-fit.variogram(a, vgm("Sph"))
  g.dummy <- gstat(formula=z~x+y, locations=~x+y, dummy=T, beta=1, model=vgm(psill=b[2,2], range=b[2,3], model='Exp'), nmax=20)
  sim_tri_reg_r1 <- predict(g.dummy, newdata=xy, nsim=30)[,c(3:32)]
  #3
  a<-variogram(dist_krant~POINT_X+POINT_Y, data_sp)
  b<-fit.variogram(a, vgm("Sph"))
  g.dummy <- gstat(formula=z~x+y, locations=~x+y, dummy=T, beta=1, model=vgm(psill=b[2,2], range=b[2,3], model='Exp'), nmax=20)
  sim_dist_krant <- predict(g.dummy, newdata=xy, nsim=30)[,c(3:32)]
  #4
  a<-variogram(prop_van_g~POINT_X+POINT_Y, data_sp)
  b<-fit.variogram(a, vgm("Sph"))
  g.dummy <- gstat(formula=z~x+y, locations=~x+y, dummy=T, beta=1, model=vgm(psill=b[2,2], range=b[2,3], model='Exp'), nmax=20)
  sim_prop_van_g <- predict(g.dummy, newdata=xy, nsim=30)[,c(3:32)]
  #5
  a<-variogram(reg_avht~POINT_X+POINT_Y, data_sp)
  b<-fit.variogram(a, vgm("Sph"))
  g.dummy <- gstat(formula=z~x+y, locations=~x+y, dummy=T, beta=1, model=vgm(psill=b[2,2], range=b[2,3], model='Exp'), nmax=20)
  sim_reg_avht <- predict(g.dummy, newdata=xy, nsim=30)[,c(3:32)]
  #6
  a<-variogram(vand_top_p~POINT_X+POINT_Y, data_sp)
  b<-fit.variogram(a, vgm("Sph"))
  g.dummy <- gstat(formula=z~x+y, locations=~x+y, dummy=T, beta=1, model=vgm(psill=b[2,2], range=b[2,3], model='Exp'), nmax=20)
  sim_vand_top_p <- predict(g.dummy, newdata=xy, nsim=30)[,c(3:32)]
  #7
  a<-variogram(topo_vid_d~POINT_X+POINT_Y, data_sp)
  b<-fit.variogram(a, vgm("Sph"))
  g.dummy <- gstat(formula=z~x+y, locations=~x+y, dummy=T, beta=1, model=vgm(psill=b[2,2], range=b[2,3], model='Exp'), nmax=20)
  sim_topo_vid_d <- predict(g.dummy, newdata=xy, nsim=30)[,c(3:32)]
  #8
  a<-variogram(litho_div~POINT_X+POINT_Y, data_sp)
  b<-fit.variogram(a, vgm("Sph"))
  g.dummy <- gstat(formula=z~x+y, locations=~x+y, dummy=T, beta=1, model=vgm(psill=b[2,2], range=b[2,3], model='Exp'), nmax=20)
  sim_litho_div <- predict(g.dummy, newdata=xy, nsim=30)[,c(3:32)]
  #9
  a<-variogram(pc1~POINT_X+POINT_Y, data_sp)
  b<-fit.variogram(a, vgm("Sph"))
  g.dummy <- gstat(formula=z~x+y, locations=~x+y, dummy=T, beta=1, model=vgm(psill=b[2,2], range=b[2,3], model='Exp'), nmax=20)
  sim_pc1 <- predict(g.dummy, newdata=xy, nsim=30)[,c(3:32)]
  #10
  a<-variogram(pc2~POINT_X+POINT_Y, data_sp)
  b<-fit.variogram(a, vgm("Sph"))
  g.dummy <- gstat(formula=z~x+y, locations=~x+y, dummy=T, beta=1, model=vgm(psill=b[2,2], range=b[2,3], model='Exp'), nmax=20)
  sim_pc2 <- predict(g.dummy, newdata=xy, nsim=30)[,c(3:32)]
  #11
  a<-variogram(pc3~POINT_X+POINT_Y, data_sp)
  b<-fit.variogram(a, vgm("Sph"))
  g.dummy <- gstat(formula=z~x+y, locations=~x+y, dummy=T, beta=1, model=vgm(psill=b[2,2], range=b[2,3], model='Exp'), nmax=20)
  sim_pc3 <- predict(g.dummy, newdata=xy, nsim=30)[,c(3:32)]
  #12
  a<-variogram(pc4~POINT_X+POINT_Y, data_sp)
  b<-fit.variogram(a, vgm("Sph"))
  g.dummy <- gstat(formula=z~x+y, locations=~x+y, dummy=T, beta=1, model=vgm(psill=b[2,2], range=b[2,3], model='Exp'), nmax=20)
  sim_pc4 <- predict(g.dummy, newdata=xy, nsim=30)[,c(3:32)]
  #13
  a<-variogram(pc5~POINT_X+POINT_Y, data_sp)
  b<-fit.variogram(a, vgm("Sph"))
  g.dummy <- gstat(formula=z~x+y, locations=~x+y, dummy=T, beta=1, model=vgm(psill=b[2,2], range=b[2,3], model='Exp'), nmax=20)
  sim_pc5 <- predict(g.dummy, newdata=xy, nsim=30)[,c(3:32)]
  #14
  a<-variogram(pc6~POINT_X+POINT_Y, data_sp)
  b<-fit.variogram(a, vgm("Sph"))
  g.dummy <- gstat(formula=z~x+y, locations=~x+y, dummy=T, beta=1, model=vgm(psill=b[2,2], range=b[2,3], model='Exp'), nmax=20)
  sim_pc6 <- predict(g.dummy, newdata=xy, nsim=30)[,c(3:32)]
  #15
  a<-variogram(pc7~POINT_X+POINT_Y, data_sp)
  b<-fit.variogram(a, vgm("Sph"))
  g.dummy <- gstat(formula=z~x+y, locations=~x+y, dummy=T, beta=1, model=vgm(psill=b[2,2], range=b[2,3], model='Exp'), nmax=20)
  sim_pc7 <- predict(g.dummy, newdata=xy, nsim=30)[,c(3:32)]
  #16
  a<-variogram(r_wdi_geo_~POINT_X+POINT_Y, data_sp)
  b<-fit.variogram(a, vgm("Sph"))
  g.dummy <- gstat(formula=z~x+y, locations=~x+y, dummy=T, beta=1, model=vgm(psill=b[2,2], range=b[2,3], model='Exp'), nmax=20)
  sim_r_wdi_geo_ <- predict(g.dummy, newdata=xy, nsim=30)[,c(3:32)]
  #17
  a<-variogram(s_litho~POINT_X+POINT_Y, data_sp)
  b<-fit.variogram(a, vgm("Sph"))
  g.dummy <- gstat(formula=z~x+y, locations=~x+y, dummy=T, beta=1, model=vgm(psill=b[2,2], range=b[2,3], model='Exp'), nmax=20)
  sim_s_litho <- predict(g.dummy, newdata=xy, nsim=30)[,c(3:32)]
  for (i in seq(30)){
    sim_s_litho[,i]<-17/(max(sim_s_litho[,i])-min(sim_s_litho[,i]))*(sim_s_litho[,i]-min(sim_s_litho[,i]))+1
  }
  sim_s_litho<-round(sim_s_litho,0)
  #18
  a<-variogram(reg_topex_~POINT_X+POINT_Y, data_sp)
  b<-fit.variogram(a, vgm("Sph"))
  g.dummy <- gstat(formula=z~x+y, locations=~x+y, dummy=T, beta=1, model=vgm(psill=b[2,2], range=b[2,3], model='Exp'), nmax=20)
  sim_reg_topex_<- predict(g.dummy, newdata=xy, nsim=30)[,c(3:32)]
  #19
  a<-variogram(reg_ls_eu~POINT_X+POINT_Y, data_sp)
  b<-fit.variogram(a, vgm("Sph"))
  g.dummy <- gstat(formula=z~x+y, locations=~x+y, dummy=T, beta=1, model=vgm(psill=b[2,2], range=b[2,3], model='Exp'), nmax=20)
  sim_reg_ls_eu <- predict(g.dummy, newdata=xy, nsim=30)[,c(3:32)]
  #20
  a<-variogram(reg_wetnes~POINT_X+POINT_Y, data_sp)
  b<-fit.variogram(a, vgm("Sph"))
  g.dummy <- gstat(formula=z~x+y, locations=~x+y, dummy=T, beta=1, model=vgm(psill=b[2,2], range=b[2,3], model='Exp'), nmax=20)
  sim_reg_wetnes <- predict(g.dummy, newdata=xy, nsim=30)[,c(3:32)]
  #
  kasyklos<-xy[round(runif(30,1,nrow(dat_neo)),0),]
  kasyklos1<-xy[round(runif(30,1,nrow(dat_neo)),0),]
  kasyklos2<-xy[round(runif(30,1,nrow(dat_neo)),0),]
  
  di1<-list()
  di2<-list()
  di3<-list()
  
  a<-rep(NA,nrow(dat_neo))
  sim_titnago_at<-data.frame(a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a)
  
  for (i in 1:30) {
    di1[[i]]<-sqrt((xy[,1]-kasyklos[i,1])^2+(xy[,2]-kasyklos[i,2])^2)
    di2[[i]]<-sqrt((xy[,1]-kasyklos1[i,1])^2+(xy[,2]-kasyklos1[i,2])^2)
    di3[[i]]<-sqrt((xy[,1]-kasyklos2[i,1])^2+(xy[,2]-kasyklos2[i,2])^2)
    for (b in 1:nrow(dat_neo)) {
      sim_titnago_at[,i][b]<-min(data.frame(di1[[i]],di2[[i]],di3[[i]])[b,])
    }
    print(i)
  }
  #Visual inspection of distance to flint simulations
  c = cut(sim_titnago_at[,4], breaks=64)
  cols = rainbow(64)[as.numeric(c)]
  plot(POINT_X,POINT_Y,col=cols, main="Atstumo iki titnago 'kasyklø' simuliacija",xlab="X",ylab="Y")
  #
  detach(train_neo_set1)
  sim_dat<-data.frame(sim_titnago_at=sim_titnago_at,sim_reg_wetnes=sim_reg_wetnes,
                      sim_reg_ls_eu=sim_reg_ls_eu,sim_reg_topex_=sim_reg_topex_,
                      sim_s_litho=sim_s_litho,sim_r_wdi_geo_=sim_r_wdi_geo_,
                      sim_pc7=sim_pc7,sim_pc6=sim_pc6,sim_pc5=sim_pc5,sim_pc4=sim_pc4,
                      sim_pc3=sim_pc3,sim_pc2=sim_pc2,sim_pc1=sim_pc1,sim_litho_div=
                        sim_litho_div,sim_topo_vid_d=sim_topo_vid_d,sim_vand_top_p=
                        sim_vand_top_p,sim_reg_avht=sim_reg_avht,sim_prop_van_g=
                        sim_prop_van_g,sim_dist_krant=sim_dist_krant,sim_tri_reg_r1=
                        sim_tri_reg_r1,sim_slait_r10k=sim_slait_r10k)
  sim_names<-c("sim_titnago_at", "sim_reg_wetnes", "sim_reg_ls_eu", "sim_reg_topex_", 
               "sim_s_litho", "sim_r_wdi_geo_", "sim_pc7", "sim_pc6", "sim_pc5", 
               "sim_pc4", "sim_pc3", "sim_pc2", "sim_pc1", "sim_litho_div", 
               "sim_topo_vid_d", "sim_vand_top_p", "sim_reg_avht", "sim_prop_van_g", 
               "sim_dist_krant", "sim_tri_reg_r1", "sim_slait_r10k")
}
#sampling simulations that correspond to selected variables
{
  norm_names<-sub(sim_names, pattern = "sim_", replacement = "")
  
  sim_neo_set1_names_ind<-numeric()
  for (i in 1:length(norm_names)) {
    sim_neo_set1_names_ind[i]<-length(which(grepl(norm_names[i],neo_full_Vars)))
  }
  
  sim_neo_set1_names_final_index<-which(sim_neo_set1_names_ind>0)
}


data<-train_neo_set1[,c("neo",neo_full_Vars)] #leaving only Neolithic sites and variables

lcp<-length(names(neo_full_Vars))*2
repeats<-30

#creating empty data structures
train<-vector(mode="list", repeats)

Gini1<-vector(mode="list", repeats)
Gini2<-vector(mode="list", repeats)
Gini3<-vector(mode="list", repeats)
Gini4<-vector(mode="list", repeats)
Gini5<-vector(mode="list", repeats)
Gini6<-vector(mode="list", repeats)
Gini7<-vector(mode="list", repeats)
Gini8<-vector(mode="list", repeats)
Gini9<-vector(mode="list", repeats)
Gini10<-vector(mode="list", repeats)

Acc1<-vector(mode="list", repeats)
Acc2<-vector(mode="list", repeats)
Acc3<-vector(mode="list", repeats)
Acc4<-vector(mode="list", repeats)
Acc5<-vector(mode="list", repeats)
Acc6<-vector(mode="list", repeats)
Acc7<-vector(mode="list", repeats)
Acc8<-vector(mode="list", repeats)
Acc9<-vector(mode="list", repeats)
Acc10<-vector(mode="list", repeats)

AUC1<-vector(mode="numeric", repeats)
AUC2<-vector(mode="numeric", repeats)
AUC3<-vector(mode="numeric", repeats)
AUC4<-vector(mode="numeric", repeats)
AUC5<-vector(mode="numeric", repeats)
AUC6<-vector(mode="numeric", repeats)
AUC7<-vector(mode="numeric", repeats)
AUC8<-vector(mode="numeric", repeats)
AUC9<-vector(mode="numeric", repeats)
AUC10<-vector(mode="numeric", repeats)
# repeated stratified 10-fold cross-validation & random forest
for (a in seq(repeats)){
  Acc1[[a]] <- NaN*seq(lcp)
  Acc2[[a]] <- NaN*seq(lcp)
  Acc3[[a]] <- NaN*seq(lcp)
  Acc4[[a]] <- NaN*seq(lcp)
  Acc5[[a]] <- NaN*seq(lcp)
  Acc6[[a]] <- NaN*seq(lcp)
  Acc7[[a]] <- NaN*seq(lcp)
  Acc8[[a]] <- NaN*seq(lcp)
  Acc9[[a]] <- NaN*seq(lcp)
  Acc10[[a]] <- NaN*seq(lcp)
  
  Gini1[[a]] <- NaN*seq(lcp)
  Gini2[[a]] <- NaN*seq(lcp)
  Gini3[[a]] <- NaN*seq(lcp)
  Gini4[[a]] <- NaN*seq(lcp)
  Gini5[[a]] <- NaN*seq(lcp)
  Gini6[[a]] <- NaN*seq(lcp)
  Gini7[[a]] <- NaN*seq(lcp)
  Gini8[[a]] <- NaN*seq(lcp)
  Gini9[[a]] <- NaN*seq(lcp)
  Gini10[[a]] <- NaN*seq(lcp)
  
  
  train[[a]]<-data.frame(data,sim_dat[,((sim_neo_set1_names_final_index-1)*30+a)])
  
  folds <- createFolds(data$neo, k = 10, list = FALSE,returnTrain = FALSE)
  #test samples
  fold1<-train[[a]][which(folds==1),]
  fold2<-train[[a]][which(folds==2),]
  fold3<-train[[a]][which(folds==3),]
  fold4<-train[[a]][which(folds==4),]
  fold5<-train[[a]][which(folds==5),]
  fold6<-train[[a]][which(folds==6),]
  fold7<-train[[a]][which(folds==7),]
  fold8<-train[[a]][which(folds==8),]
  fold9<-train[[a]][which(folds==9),]
  fold10<-train[[a]][which(folds==10),]
  #Samples used to build models
  train1<-train[[a]][-which(folds==1),]
  train2<-train[[a]][-which(folds==2),]
  train3<-train[[a]][-which(folds==3),]
  train4<-train[[a]][-which(folds==4),]
  train5<-train[[a]][-which(folds==5),]
  train6<-train[[a]][-which(folds==6),]
  train7<-train[[a]][-which(folds==7),]
  train8<-train[[a]][-which(folds==8),]
  train9<-train[[a]][-which(folds==9),]
  train10<-train[[a]][-which(folds==10),]
  
  
  
  ### 1 
  miskas<-randomForest(x=train1[,-1],y=train1[,1],ntree=500,
                       mtry=1,replace=TRUE, strata= train1[,1],
                       sampsize = c(round(length(train1[,1][which(train1[,1]==1)])*0.75,0),
                                    round(length(train1[,1][which(train1[,1]==1)])*0.75,0)),nodesize=2,
                       importance = TRUE, norm.votes = TRUE, keep.forest = TRUE)
  # Importance
  Gini1[[a]]<-importance(miskas)[,4]
  Acc1[[a]]<-(importance(miskas,scale = F)[,1]+importance(miskas,scale = F)[,2])/2
  # AUROC calculation
  pred1<-predict(miskas,fold1,type = "prob")
  predi1 <- prediction(pred1[,2],fold1[,1])
  perf1 <- performance(predi1, measure = "auc")
  AUC1[a]<-as.numeric(perf1@y.values)
  
  print((((a-1)*10+1)/300)*100)
  
  ### 2 
  miskas<-randomForest(x=train2[,-1],y=train2[,1],ntree=500,
                       mtry=1,replace=TRUE, strata= train2[,1],
                       sampsize = c(round(length(train2[,1][which(train2[,1]==1)])*0.75,0),
                                    round(length(train2[,1][which(train2[,1]==1)])*0.75,0)),nodesize=2,
                       importance = TRUE, norm.votes = TRUE, keep.forest = TRUE)
  # Importance
  Gini2[[a]]<-importance(miskas)[,4]
  Acc2[[a]]<-(importance(miskas,scale = F)[,1]+importance(miskas,scale = F)[,2])/2
  # AUROC calculation
  pred1<-predict(miskas,fold2,type = "prob")
  predi1 <- prediction(pred1[,2],fold2[,1])
  perf1 <- performance(predi1, measure = "auc")
  AUC2[a]<-as.numeric(perf1@y.values)
  
  print((((a-1)*10+2)/300)*100)
  
  ### 3  
  miskas<-randomForest(x=train3[,-1],y=train3[,1],ntree=500,
                       mtry=1,replace=TRUE, strata= train3[,1],
                       sampsize = c(round(length(train3[,1][which(train3[,1]==1)])*0.75,0),
                                    round(length(train3[,1][which(train3[,1]==1)])*0.75,0)),nodesize=2,
                       importance = TRUE, norm.votes = TRUE, keep.forest = TRUE)
  # Importance
  Gini3[[a]]<-importance(miskas)[,4]
  Acc3[[a]]<-(importance(miskas,scale = F)[,1]+importance(miskas,scale = F)[,2])/2
  # AUROC calculation
  pred1<-predict(miskas,fold3,type = "prob")
  predi1 <- prediction(pred1[,2],fold3[,1])
  perf1 <- performance(predi1, measure = "auc")
  AUC3[a]<-as.numeric(perf1@y.values)
  
  print((((a-1)*10+3)/300)*100)
  
  ### 4  
  miskas<-randomForest(x=train4[,-1],y=train4[,1],ntree=500,
                       mtry=1,replace=TRUE, strata= train4[,1],
                       sampsize = c(round(length(train4[,1][which(train4[,1]==1)])*0.75,0),
                                    round(length(train4[,1][which(train4[,1]==1)])*0.75,0)),nodesize=2,
                       importance = TRUE, norm.votes = TRUE, keep.forest = TRUE)
  # Importance
  Gini4[[a]]<-importance(miskas)[,4]
  Acc4[[a]]<-(importance(miskas,scale = F)[,1]+importance(miskas,scale = F)[,2])/2
  # AUROC calculation
  pred1<-predict(miskas,fold4,type = "prob")
  predi1 <- prediction(pred1[,2],fold4[,1])
  perf1 <- performance(predi1, measure = "auc")
  AUC4[a]<-as.numeric(perf1@y.values)
  
  print((((a-1)*10+4)/300)*100)
  
  ### 5  
  miskas<-randomForest(x=train5[,-1],y=train5[,1],ntree=500,
                       mtry=1,replace=TRUE, strata= train5[,1],
                       sampsize = c(round(length(train5[,1][which(train5[,1]==1)])*0.75,0),
                                    round(length(train5[,1][which(train5[,1]==1)])*0.75,0)),nodesize=2,
                       importance = TRUE, norm.votes = TRUE, keep.forest = TRUE)
  # Importance
  Gini5[[a]]<-importance(miskas)[,4]
  Acc5[[a]]<-(importance(miskas,scale = F)[,1]+importance(miskas,scale = F)[,2])/2
  # AUROC calculation
  pred1<-predict(miskas,fold5,type = "prob")
  predi1 <- prediction(pred1[,2],fold5[,1])
  perf1 <- performance(predi1, measure = "auc")
  AUC5[a]<-as.numeric(perf1@y.values)
  
  print((((a-1)*10+5)/300)*100)
  
  ### 6  
  miskas<-randomForest(x=train6[,-1],y=train6[,1],ntree=500,
                       mtry=1,replace=TRUE, strata= train6[,1],
                       sampsize = c(round(length(train6[,1][which(train6[,1]==1)])*0.75,0),
                                    round(length(train6[,1][which(train6[,1]==1)])*0.75,0)),nodesize=2,
                       importance = TRUE, norm.votes = TRUE, keep.forest = TRUE)
  # Importance
  Gini6[[a]]<-importance(miskas)[,4]
  Acc6[[a]]<-(importance(miskas,scale = F)[,1]+importance(miskas,scale = F)[,2])/2
  # AUROC calculation
  pred1<-predict(miskas,fold6,type = "prob")
  predi1 <- prediction(pred1[,2],fold6[,1])
  perf1 <- performance(predi1, measure = "auc")
  AUC6[a]<-as.numeric(perf1@y.values)
  
  print((((a-1)*10+6)/300)*100)
  
  ### 7  
  miskas<-randomForest(x=train7[,-1],y=train7[,1],ntree=500,
                       mtry=1,replace=TRUE, strata= train7[,1],
                       sampsize = c(round(length(train7[,1][which(train7[,1]==1)])*0.75,0),
                                    round(length(train7[,1][which(train7[,1]==1)])*0.75,0)),nodesize=2,
                       importance = TRUE, norm.votes = TRUE, keep.forest = TRUE)
  # Importance
  Gini7[[a]]<-importance(miskas)[,4]
  Acc7[[a]]<-(importance(miskas,scale = F)[,1]+importance(miskas,scale = F)[,2])/2
  # AUROC calculation
  pred1<-predict(miskas,fold7,type = "prob")
  predi1 <- prediction(pred1[,2],fold7[,1])
  perf1 <- performance(predi1, measure = "auc")
  AUC7[a]<-as.numeric(perf1@y.values)
  
  print((((a-1)*10+7)/300)*100)
  
  ### 8  
  miskas<-randomForest(x=train8[,-1],y=train8[,1],ntree=500,
                       mtry=1,replace=TRUE, strata= train8[,1],
                       sampsize = c(round(length(train8[,1][which(train8[,1]==1)])*0.75,0),
                                    round(length(train8[,1][which(train8[,1]==1)])*0.75,0)),nodesize=2,
                       importance = TRUE, norm.votes = TRUE, keep.forest = TRUE)
  # Importance
  Gini8[[a]]<-importance(miskas)[,4]
  Acc8[[a]]<-(importance(miskas,scale = F)[,1]+importance(miskas,scale = F)[,2])/2
  # AUROC calculation
  pred1<-predict(miskas,fold8,type = "prob")
  predi1 <- prediction(pred1[,2],fold8[,1])
  perf1 <- performance(predi1, measure = "auc")
  AUC8[a]<-as.numeric(perf1@y.values)
  
  print((((a-1)*10+8)/300)*100)
  
  ### 9  
  miskas<-randomForest(x=train9[,-1],y=train9[,1],ntree=500,
                       mtry=1,replace=TRUE, strata= train9[,1],
                       sampsize = c(round(length(train9[,1][which(train9[,1]==1)])*0.75,0),
                                    round(length(train9[,1][which(train9[,1]==1)])*0.75,0)),nodesize=2,
                       importance = TRUE, norm.votes = TRUE, keep.forest = TRUE)
  # Importance
  Gini9[[a]]<-importance(miskas)[,4]
  Acc9[[a]]<-(importance(miskas,scale = F)[,1]+importance(miskas,scale = F)[,2])/2
  # AUROC calculation
  pred1<-predict(miskas,fold9,type = "prob")
  predi1 <- prediction(pred1[,2],fold9[,1])
  perf1 <- performance(predi1, measure = "auc")
  AUC9[a]<-as.numeric(perf1@y.values)
  
  print((((a-1)*10+9)/300)*100)
  
  ### 10  
  miskas<-randomForest(x=train10[,-1],y=train10[,1],ntree=500,
                       mtry=1,replace=TRUE, strata= train10[,1],
                       sampsize = c(round(length(train10[,1][which(train10[,1]==1)])*0.75,0),
                                    round(length(train10[,1][which(train10[,1]==1)])*0.75,0)),nodesize=2,
                       importance = TRUE, norm.votes = TRUE, keep.forest = TRUE)
  # Importance
  Gini10[[a]]<-importance(miskas)[,4]
  Acc10[[a]]<-(importance(miskas,scale = F)[,1]+importance(miskas,scale = F)[,2])/2
  # AUROC calculation
  pred1<-predict(miskas,fold10,type = "prob")
  predi1 <- prediction(pred1[,2],fold10[,1])
  perf1 <- performance(predi1, measure = "auc")
  AUC10[a]<-as.numeric(perf1@y.values)
  
  print((((a-1)*10+10)/300)*100)
}
### proccesing, inspecting and obtaining relevant results
{
  Gini_neo_set1<-rbind.data.frame(t(do.call(cbind.data.frame,Gini1)),t(do.call(cbind.data.frame,Gini2)),
                                  t(do.call(cbind.data.frame,Gini3)),t(do.call(cbind.data.frame,Gini4)),
                                  t(do.call(cbind.data.frame,Gini5)),t(do.call(cbind.data.frame,Gini6)),
                                  t(do.call(cbind.data.frame,Gini7)),t(do.call(cbind.data.frame,Gini8)),
                                  t(do.call(cbind.data.frame,Gini9)),t(do.call(cbind.data.frame,Gini10)))
  Acc_neo_set1<-rbind.data.frame(t(do.call(cbind.data.frame,Acc1)),t(do.call(cbind.data.frame,Acc2)),
                                 t(do.call(cbind.data.frame,Acc3)),t(do.call(cbind.data.frame,Acc4)),
                                 t(do.call(cbind.data.frame,Acc5)),t(do.call(cbind.data.frame,Acc6)),
                                 t(do.call(cbind.data.frame,Acc7)),t(do.call(cbind.data.frame,Acc8)),
                                 t(do.call(cbind.data.frame,Acc9)),t(do.call(cbind.data.frame,Acc10)))
  ### plots 
  par(cex.axis=0.6)
  boxplot(Gini_neo_set1,las=2, main="GINI index")
  boxplot(Acc_neo_set1,las=2, main="Accuracy reduction after permutation")
  ###
  Acc_neo_set1_sim<-Acc_neo_set1[,names(Acc_neo_set1)[(length(neo_full_Vars)+1):(length(neo_full_Vars)*2)]]
  Acc_neo_set1_org<-Acc_neo_set1[,norm_names[sim_neo_set1_names_final_index]]
  Acc_neo_set1_mean_org<-apply(Acc_neo_set1_org,2,mean)
  Acc_neo_set1_mean_sim<-apply(Acc_neo_set1_sim,2,mean)
  Acc_neo_set1_uplim_all<-max(apply(Acc_neo_set1[,1:(length(neo_full_Vars)*2)],2,mean)+apply(Acc_neo_set1[,1:(length(neo_full_Vars)*2)],2,sd))
  Acc_neo_set1_downlim_all<-min(apply(Acc_neo_set1[,1:(length(neo_full_Vars)*2)],2,mean)-apply(Acc_neo_set1[,1:(length(neo_full_Vars)*2)],2,sd))
  Acc_neo_set1_uplim_sim<-Acc_neo_set1_mean_sim+apply(Acc_neo_set1_sim,2,sd)
  Acc_neo_set1_downlim_sim<-Acc_neo_set1_mean_sim-apply(Acc_neo_set1_sim,2,sd)
  Acc_neo_set1_uplim_org<-Acc_neo_set1_mean_org+apply(Acc_neo_set1_org,2,sd)
  Acc_neo_set1_downlim_org<-Acc_neo_set1_mean_org-apply(Acc_neo_set1_org,2,sd)
  ###
  plot(1:length(neo_full_Vars),Acc_neo_set1_mean_org, ylim =
         c(Acc_neo_set1_downlim_all,Acc_neo_set1_uplim_all),xaxt='n',xlab="",
       type = "b", main= "ACCURACY DECREASE \n AFTER x VARIABLE PERMUTATION")
  polygon(c(1:length(neo_full_Vars),rev(1:length(neo_full_Vars))), c(Acc_neo_set1_downlim_org,
                                                                     rev(Acc_neo_set1_uplim_org)), col = "grey75", border = T)
  abline(h=0)
  points(1:length(neo_full_Vars), Acc_neo_set1_mean_org)
  lines(1:length(neo_full_Vars), Acc_neo_set1_mean_org, lwd = 2)
  axis(1, at=1:length(neo_full_Vars), labels=names(Acc_neo_set1_mean_org), lty=2, col=1, las=2) 
  lines(1:length(neo_full_Vars),Acc_neo_set1_uplim_sim,type = "b",col=2)
  lines(1:length(neo_full_Vars),Acc_neo_set1_downlim_sim,type = "b",col=2)
  lines(1:length(neo_full_Vars),Acc_neo_set1_mean_sim,type = "b",col=3)
  ###
  Acc_var_logic_set1_neo<-numeric()
  for (i in 1:length(neo_full_Vars)) {
    Acc_var_logic_set1_neo[i]<-ecdf(Acc_neo_set1_sim[,i])(mean(Acc_neo_set1_org[,i]))
  }
  names(Acc_var_logic_set1_neo)<-names(Acc_neo_set1_org)
  Acc_var_ind_set1_neo<-which(Acc_var_logic_set1_neo>=0.95) # variable selection criteria
  #
  neo_set1_Vars_Acc<-names(Acc_neo_set1_org)[Acc_var_ind_set1_neo]
  ###
  #"reg_wetnes" "reg_ls_eu"  "pc2"        "topo_vid_d" "vand_top_p" "dist_krant"
  # "tri_reg_r1" "slait_r10k"
  ###
  ##the same operations repeated for GINI Importance estimates:
  Gini_neo_set1_sim<-Gini_neo_set1[,names(Gini_neo_set1)[(length(neo_full_Vars)+1):(length(neo_full_Vars)*2)]]
  Gini_neo_set1_org<-Gini_neo_set1[,norm_names[sim_neo_set1_names_final_index]]
  Gini_neo_set1_mean_org<-apply(Gini_neo_set1_org,2,mean)
  Gini_neo_set1_mean_sim<-apply(Gini_neo_set1_sim,2,mean)
  Gini_neo_set1_uplim_all<-max(apply(Gini_neo_set1[,1:(length(neo_full_Vars)*2)],2,mean)+apply(Gini_neo_set1[,1:(length(neo_full_Vars)*2)],2,sd))
  Gini_neo_set1_downlim_all<-min(apply(Gini_neo_set1[,1:(length(neo_full_Vars)*2)],2,mean)-apply(Gini_neo_set1[,1:(length(neo_full_Vars)*2)],2,sd))
  Gini_neo_set1_uplim_sim<-Gini_neo_set1_mean_sim+apply(Gini_neo_set1_sim,2,sd)
  Gini_neo_set1_downlim_sim<-Gini_neo_set1_mean_sim-apply(Gini_neo_set1_sim,2,sd)
  Gini_neo_set1_uplim_org<-Gini_neo_set1_mean_org+apply(Gini_neo_set1_org,2,sd)
  Gini_neo_set1_downlim_org<-Gini_neo_set1_mean_org-apply(Gini_neo_set1_org,2,sd)
  ###
  plot(1:length(neo_full_Vars),Gini_neo_set1_mean_org, ylim =
         c(Gini_neo_set1_downlim_all,Gini_neo_set1_uplim_all),xaxt='n',xlab="",
       type = "b", main= "GINI INDEX \n AFTER x VARIABLE PERMUTATION")
  polygon(c(1:length(neo_full_Vars),rev(1:length(neo_full_Vars))), c(Gini_neo_set1_downlim_org,
                                                                     rev(Gini_neo_set1_uplim_org)), col = "grey75", border = T)
  abline(h=0)
  points(1:length(neo_full_Vars), Gini_neo_set1_mean_org)
  lines(1:length(neo_full_Vars), Gini_neo_set1_mean_org, lwd = 2)
  axis(1, at=1:length(neo_full_Vars), labels=names(Gini_neo_set1_mean_org), lty=2, col=1, las=2) 
  lines(1:length(neo_full_Vars),Gini_neo_set1_uplim_sim,type = "b",col=2)
  lines(1:length(neo_full_Vars),Gini_neo_set1_downlim_sim,type = "b",col=2)
  lines(1:length(neo_full_Vars),Gini_neo_set1_mean_sim,type = "b",col=3)
  ###
  Gini_var_logic_set1_neo<-numeric()
  for (i in 1:length(neo_full_Vars)) {
    Gini_var_logic_set1_neo[i]<-ecdf(Gini_neo_set1_sim[,i])(mean(Gini_neo_set1_org[,i]))
  }
  names(Gini_var_logic_set1_neo)<-names(Gini_neo_set1_org)
  Gini_var_ind_set1_neo<-which(Gini_var_logic_set1_neo>=0.95) #variable selection criteria
  #
  neo_set1_Vars_Gini<-names(Gini_neo_set1_org)[Gini_var_ind_set1_neo]
  ###
  #reg_wetnes" "reg_ls_eu"  "pc2"        "topo_vid_d" "vand_top_p" "dist_krant"
  # "tri_reg_r1" "slait_r10k"
  ###
}
#combining names of significant variables and deleting duplicates
neo_set1_Vars<-c(neo_set1_Vars_Acc,neo_set1_Vars_Gini)
neo_set1_Vars<-neo_set1_Vars[!duplicated(neo_set1_Vars)]
###
for_cor_neo_set1<-data[,norm_names[sim_neo_set1_names_final_index]][,neo_set1_Vars]
correl_neo<-cor(for_cor_neo_set1,method = "spearman")
abs(correl_neo)>=0.79 # removing highly correlated variables
### leaving more important correlates. If GINI and Permutation importances disagree, then Permutation Importance used for final decision
Acc_neo_set1_mean_org-Acc_neo_set1_mean_sim
Gini_neo_set1_mean_org-Gini_neo_set1_mean_sim
Gini_var_logic_set1_neo
Acc_var_logic_set1_neo
### selecting final variables for predictive modelling with GAM
neo_set1_Vars<-c("pc2","vand_top_p","dist_krant",
                 "tri_reg_r1")

###GAM

# creating training data
train_neo_set2<-dat_neo[,c("POINT_X","POINT_Y", "neo",neo_set1_Vars)]
attach(train_neo_set2)
train_neo_set2[,3]<-factor(neo)
rownames(train_neo_set2)<-as.character(seq(nrow(train_neo_set2)))
train<-train_neo_set2[,-1:-2]
detach(train_neo_set2)
#loading libraries
library(scoring)
library(mgcv)
#
repeats<-30
#creating empty data structures
AUC1<-vector(mode="numeric", repeats)
AUC2<-vector(mode="numeric", repeats)
AUC3<-vector(mode="numeric", repeats)
AUC4<-vector(mode="numeric", repeats)
AUC5<-vector(mode="numeric", repeats)
AUC6<-vector(mode="numeric", repeats)
AUC7<-vector(mode="numeric", repeats)
AUC8<-vector(mode="numeric", repeats)
AUC9<-vector(mode="numeric", repeats)
AUC10<-vector(mode="numeric", repeats)

Score1_1<-vector(mode="numeric", repeats)
Score2_1<-vector(mode="numeric", repeats)
Score3_1<-vector(mode="numeric", repeats)
Score4_1<-vector(mode="numeric", repeats)
Score5_1<-vector(mode="numeric", repeats)
Score6_1<-vector(mode="numeric", repeats)
Score7_1<-vector(mode="numeric", repeats)
Score8_1<-vector(mode="numeric", repeats)
Score9_1<-vector(mode="numeric", repeats)
Score10_1<-vector(mode="numeric", repeats)

Score1_0<-vector(mode="numeric", repeats)
Score2_0<-vector(mode="numeric", repeats)
Score3_0<-vector(mode="numeric", repeats)
Score4_0<-vector(mode="numeric", repeats)
Score5_0<-vector(mode="numeric", repeats)
Score6_0<-vector(mode="numeric", repeats)
Score7_0<-vector(mode="numeric", repeats)
Score8_0<-vector(mode="numeric", repeats)
Score9_0<-vector(mode="numeric", repeats)
Score10_0<-vector(mode="numeric", repeats)

Total_Score1<-vector(mode="numeric", repeats)
Total_Score2<-vector(mode="numeric", repeats)
Total_Score3<-vector(mode="numeric", repeats)
Total_Score4<-vector(mode="numeric", repeats)
Total_Score5<-vector(mode="numeric", repeats)
Total_Score6<-vector(mode="numeric", repeats)
Total_Score7<-vector(mode="numeric", repeats)
Total_Score8<-vector(mode="numeric", repeats)
Total_Score9<-vector(mode="numeric", repeats)
Total_Score10<-vector(mode="numeric", repeats)

c<-mean(as.numeric(train$neo)-1)
# repeated stratified 10-fold cross-validation & GAM
for (a in seq(repeats)) {
  
  folds <- createFolds(factor(train$neo), k = 10, list = FALSE,returnTrain = FALSE)
  
  fold1<-train[which(folds==1),]
  fold2<-train[which(folds==2),]
  fold3<-train[which(folds==3),]
  fold4<-train[which(folds==4),]
  fold5<-train[which(folds==5),]
  fold6<-train[which(folds==6),]
  fold7<-train[which(folds==7),]
  fold8<-train[which(folds==8),]
  fold9<-train[which(folds==9),]
  fold10<-train[which(folds==10),]
  
  train1<-train[-which(folds==1),]
  train2<-train[-which(folds==2),]
  train3<-train[-which(folds==3),]
  train4<-train[-which(folds==4),]
  train5<-train[-which(folds==5),]
  train6<-train[-which(folds==6),]
  train7<-train[-which(folds==7),]
  train8<-train[-which(folds==8),]
  train9<-train[-which(folds==9),]
  train10<-train[-which(folds==10),]
  
  gam_model<-gam(neo ~ s(dist_krant) + s(pc2) + s(vand_top_p) +
                   s(tri_reg_r1), family = binomial, select=TRUE, method="REML", data = train1)
  #AUROC
  pred<-predict(gam_model,fold1,type = "response")
  pred1 <- prediction(as.vector(pred),fold1[,1])
  peri <- performance(pred1, measure = "auc")
  AUC1[a]<-as.numeric(peri@y.values)
  #Beta family Proper scoring rule with alpha = 1 and beta=1/c - 1
  Score1_1[a]<-mean(calcscore(as.numeric(fold1[,1])-1~pred,fam="beta", param= c(1,-1+1/c),bounds = c(0,1))[which(fold1[,1]==1)])
  Score1_0[a]<-mean(calcscore(as.numeric(fold1[,1])-1~pred,fam="beta", param= c(1,-1+1/c),bounds = c(0,1))[which(fold1[,1]==0)])
  Total_Score1[a]<-mean(calcscore(as.numeric(fold1[,1])-1~pred,fam="beta", param= c(1,-1+1/c),bounds = c(0,1)))
  ###
  print((((a-1)*10+1)/300)*100)
  ###
  gam_model<-gam(neo ~ s(dist_krant) + s(pc2) + s(vand_top_p) +
                   s(tri_reg_r1), family = binomial, select=TRUE, method="REML", data = train2)
  #AUROC
  pred<-predict(gam_model,fold2,type = "response")
  pred1 <- prediction(as.vector(pred),fold2[,1])
  peri <- performance(pred1, measure = "auc")
  AUC2[a]<-as.numeric(peri@y.values)
  #Beta family Proper scoring rule with alpha = 1 and beta=1/c - 1
  Score2_1[a]<-mean(calcscore(as.numeric(fold2[,1])-1~pred,fam="beta", param= c(1,-1+1/c),bounds = c(0,1))[which(fold2[,1]==1)])
  Score2_0[a]<-mean(calcscore(as.numeric(fold2[,1])-1~pred,fam="beta", param= c(1,-1+1/c),bounds = c(0,1))[which(fold2[,1]==0)])
  Total_Score2[a]<-mean(calcscore(as.numeric(fold2[,1])-1~pred,fam="beta", param= c(1,-1+1/c),bounds = c(0,1)))
  ###
  print((((a-1)*10+2)/300)*100)
  ###
  gam_model<-gam(neo ~ s(dist_krant) + s(pc2) + s(vand_top_p) +
                   s(tri_reg_r1), family = binomial, select=TRUE, method="REML", data = train3)
  #AUROC
  pred<-predict(gam_model,fold3,type = "response")
  pred1 <- prediction(as.vector(pred),fold3[,1])
  peri <- performance(pred1, measure = "auc")
  AUC3[a]<-as.numeric(peri@y.values)
  #Beta family Proper scoring rule with alpha = 1 and beta=1/c - 1
  Score3_1[a]<-mean(calcscore(as.numeric(fold3[,1])-1~pred,fam="beta", param= c(1,-1+1/c),bounds = c(0,1))[which(fold3[,1]==1)])
  Score3_0[a]<-mean(calcscore(as.numeric(fold3[,1])-1~pred,fam="beta", param= c(1,-1+1/c),bounds = c(0,1))[which(fold3[,1]==0)])
  Total_Score3[a]<-mean(calcscore(as.numeric(fold3[,1])-1~pred,fam="beta", param= c(1,-1+1/c),bounds = c(0,1)))
  ###
  print((((a-1)*10+3)/300)*100)
  ###
  gam_model<-gam(neo ~ s(dist_krant) + s(pc2) + s(vand_top_p) +
                   s(tri_reg_r1), family = binomial, select=TRUE, method="REML", data = train4)
  #AUROC
  pred<-predict(gam_model,fold4,type = "response")
  pred1 <- prediction(as.vector(pred),fold4[,1])
  peri <- performance(pred1, measure = "auc")
  AUC4[a]<-as.numeric(peri@y.values)
  #Beta family Proper scoring rule with alpha = 1 and beta=1/c - 1
  Score4_1[a]<-mean(calcscore(as.numeric(fold4[,1])-1~pred,fam="beta", param= c(1,-1+1/c),bounds = c(0,1))[which(fold4[,1]==1)])
  Score4_0[a]<-mean(calcscore(as.numeric(fold4[,1])-1~pred,fam="beta", param= c(1,-1+1/c),bounds = c(0,1))[which(fold4[,1]==0)])
  Total_Score4[a]<-mean(calcscore(as.numeric(fold4[,1])-1~pred,fam="beta", param= c(1,-1+1/c),bounds = c(0,1)))
  ###
  print((((a-1)*10+4)/300)*100)
  ###
  gam_model<-gam(neo ~ s(dist_krant) + s(pc2) + s(vand_top_p) +
                   s(tri_reg_r1), family = binomial, select=TRUE, method="REML", data = train5)
  #AUROC
  pred<-predict(gam_model,fold5,type = "response")
  pred1 <- prediction(as.vector(pred),fold5[,1])
  peri <- performance(pred1, measure = "auc")
  AUC5[a]<-as.numeric(peri@y.values)
  #Beta family Proper scoring rule with alpha = 1 and beta=1/c - 1
  Score5_1[a]<-mean(calcscore(as.numeric(fold5[,1])-1~pred,fam="beta", param= c(1,-1+1/c),bounds = c(0,1))[which(fold5[,1]==1)])
  Score5_0[a]<-mean(calcscore(as.numeric(fold5[,1])-1~pred,fam="beta", param= c(1,-1+1/c),bounds = c(0,1))[which(fold5[,1]==0)])
  Total_Score5[a]<-mean(calcscore(as.numeric(fold5[,1])-1~pred,fam="beta", param= c(1,-1+1/c),bounds = c(0,1)))
  ###
  print((((a-1)*10+5)/300)*100)
  ###
  gam_model<-gam(neo ~ s(dist_krant) + s(pc2) + s(vand_top_p) +
                   s(tri_reg_r1), family = binomial, select=TRUE, method="REML", data = train6)
  #AUROC
  pred<-predict(gam_model,fold6,type = "response")
  pred1 <- prediction(as.vector(pred),fold6[,1])
  peri <- performance(pred1, measure = "auc")
  AUC6[a]<-as.numeric(peri@y.values)
  #Beta family Proper scoring rule with alpha = 1 and beta=1/c - 1
  Score6_1[a]<-mean(calcscore(as.numeric(fold6[,1])-1~pred,fam="beta", param= c(1,-1+1/c),bounds = c(0,1))[which(fold6[,1]==1)])
  Score6_0[a]<-mean(calcscore(as.numeric(fold6[,1])-1~pred,fam="beta", param= c(1,-1+1/c),bounds = c(0,1))[which(fold6[,1]==0)])
  Total_Score6[a]<-mean(calcscore(as.numeric(fold6[,1])-1~pred,fam="beta", param= c(1,-1+1/c),bounds = c(0,1)))
  ###
  print((((a-1)*10+6)/300)*100)
  ###
  gam_model<-gam(neo ~ s(dist_krant) + s(pc2) + s(vand_top_p) +
                   s(tri_reg_r1), family = binomial, select=TRUE, method="REML", data = train7)
  #AUROC
  pred<-predict(gam_model,fold7,type = "response")
  pred1 <- prediction(as.vector(pred),fold7[,1])
  peri <- performance(pred1, measure = "auc")
  AUC7[a]<-as.numeric(peri@y.values)
  #Beta family Proper scoring rule with alpha = 1 and beta=1/c - 1
  Score7_1[a]<-mean(calcscore(as.numeric(fold7[,1])-1~pred,fam="beta", param= c(1,-1+1/c),bounds = c(0,1))[which(fold7[,1]==1)])
  Score7_0[a]<-mean(calcscore(as.numeric(fold7[,1])-1~pred,fam="beta", param= c(1,-1+1/c),bounds = c(0,1))[which(fold7[,1]==0)])
  Total_Score7[a]<-mean(calcscore(as.numeric(fold7[,1])-1~pred,fam="beta", param= c(1,-1+1/c),bounds = c(0,1)))
  ###
  print((((a-1)*10+7)/300)*100)
  ###
  gam_model<-gam(neo ~ s(dist_krant) + s(pc2) + s(vand_top_p) +
                   s(tri_reg_r1),family = binomial, select=TRUE, method="REML", data = train8)
  #AUROC
  pred<-predict(gam_model,fold8,type = "response")
  pred1 <- prediction(as.vector(pred),fold8[,1])
  peri <- performance(pred1, measure = "auc")
  AUC8[a]<-as.numeric(peri@y.values)
  #Beta family Proper scoring rule with alpha = 1 and beta=1/c - 1
  Score8_1[a]<-mean(calcscore(as.numeric(fold8[,1])-1~pred,fam="beta", param= c(1,-1+1/c),bounds = c(0,1))[which(fold8[,1]==1)])
  Score8_0[a]<-mean(calcscore(as.numeric(fold8[,1])-1~pred,fam="beta", param= c(1,-1+1/c),bounds = c(0,1))[which(fold8[,1]==0)])
  Total_Score8[a]<-mean(calcscore(as.numeric(fold8[,1])-1~pred,fam="beta", param= c(1,-1+1/c),bounds = c(0,1)))
  ###
  print((((a-1)*10+8)/300)*100)
  ###
  gam_model<-gam(neo ~ s(dist_krant) + s(pc2) + s(vand_top_p) +
                   s(tri_reg_r1), family = binomial, select=TRUE, method="REML", data = train9)
  #AUROC
  pred<-predict(gam_model,fold9,type = "response")
  pred1 <- prediction(as.vector(pred),fold9[,1])
  peri <- performance(pred1, measure = "auc")
  AUC9[a]<-as.numeric(peri@y.values)
  #Beta family Proper scoring rule with alpha = 1 and beta=1/c - 1
  Score9_1[a]<-mean(calcscore(as.numeric(fold9[,1])-1~pred,fam="beta", param= c(1,-1+1/c),bounds = c(0,1))[which(fold9[,1]==1)])
  Score9_0[a]<-mean(calcscore(as.numeric(fold9[,1])-1~pred,fam="beta", param= c(1,-1+1/c),bounds = c(0,1))[which(fold9[,1]==0)])
  Total_Score9[a]<-mean(calcscore(as.numeric(fold9[,1])-1~pred,fam="beta", param= c(1,-1+1/c),bounds = c(0,1)))
  ###
  print((((a-1)*10+9)/300)*100)
  ###
  gam_model<-gam(neo ~ s(dist_krant) + s(pc2) + s(vand_top_p) +
                   s(tri_reg_r1), family = binomial, select=TRUE, method="REML", data = train10)
  #AUROC
  pred<-predict(gam_model,fold10,type = "response")
  pred1 <- prediction(as.vector(pred),fold10[,1])
  peri <- performance(pred1, measure = "auc")
  AUC10[a]<-as.numeric(peri@y.values)
  #Beta family Proper scoring rule with alpha = 1 and beta=1/c - 1
  Score10_1[a]<-mean(calcscore(as.numeric(fold10[,1])-1~pred,fam="beta", param= c(1,-1+1/c),bounds = c(0,1))[which(fold10[,1]==1)])
  Score10_0[a]<-mean(calcscore(as.numeric(fold10[,1])-1~pred,fam="beta", param= c(1,-1+1/c),bounds = c(0,1))[which(fold10[,1]==0)])
  Total_Score10[a]<-mean(calcscore(as.numeric(fold10[,1])-1~pred,fam="beta", param= c(1,-1+1/c),bounds = c(0,1)))
  ###
  print((((a-1)*10+10)/300)*100)
}
### processing, inspecting and obtaining relevant results
{
  AUC_GAM_neo_org<-c(AUC1,AUC2,AUC3,AUC4,AUC5,AUC6,AUC8,AUC9,AUC10)
  Score_1_GAM_neo_org<-c(Score1_1,Score2_1,Score3_1,Score4_1,Score5_1,Score6_1,
                         Score7_1,Score8_1,Score9_1,Score10_1)
  Score_0_GAM_neo_org<-c(Score1_0,Score2_0,Score3_0,Score4_0,Score5_0,Score6_0,
                         Score7_0,Score8_0,Score9_0,Score10_0)
  Total_score_GAM_neo_org<-c(Total_Score1,Total_Score2,Total_Score3,
                             Total_Score4,Total_Score5,Total_Score6,
                             Total_Score7,Total_Score8,Total_Score9,
                             Total_Score10)
  mean(AUC_GAM_neo_org) #  0.8599881
  boxplot(AUC_GAM_neo_org)
  
  mean(Score_1_GAM_neo_org) # 0.1850581
  boxplot(Score_1_GAM_neo_org)
  
  mean(Score_0_GAM_neo_org) # 0.005661996
  boxplot(Score_0_GAM_neo_org)
  
  mean(Total_score_GAM_neo_org) # 0.01118553
  boxplot(Total_score_GAM_neo_org)
}
#Repeating GAM building proccess with simulated variables
#simulating selected variables variables
{
  attach(train_neo_set2)
  data_sp<-train_neo_set2
  library(gstat)
  library(sp)
  coordinates(data_sp)=~POINT_X+POINT_Y
  xy<-data.frame(POINT_X,POINT_Y)
  names(xy) <- c('x','y')
  #1
  a<-variogram(tri_reg_r1~POINT_X+POINT_Y, data_sp)
  b<-fit.variogram(a, vgm("Sph"))
  g.dummy <- gstat(formula=z~x+y, locations=~x+y, dummy=T, beta=1, model=vgm(psill=b[2,2], range=b[2,3], model='Exp'), nmax=20)
  sim_tri_reg_r1 <- predict(g.dummy, newdata=xy, nsim=30)[,c(3:32)]
  #2
  a<-variogram(vand_top_p~POINT_X+POINT_Y, data_sp)
  b<-fit.variogram(a, vgm("Sph"))
  g.dummy <- gstat(formula=z~x+y, locations=~x+y, dummy=T, beta=1, model=vgm(psill=b[2,2], range=b[2,3], model='Exp'), nmax=20)
  sim_vand_top_p <- predict(g.dummy, newdata=xy, nsim=30)[,c(3:32)]
  #3
  a<-variogram(pc2~POINT_X+POINT_Y, data_sp)
  b<-fit.variogram(a, vgm("Sph"))
  g.dummy <- gstat(formula=z~x+y, locations=~x+y, dummy=T, beta=1, model=vgm(psill=b[2,2], range=b[2,3], model='Exp'), nmax=20)
  sim_pc2 <- predict(g.dummy, newdata=xy, nsim=30)[,c(3:32)]
  #4
  a<-variogram(dist_krant~POINT_X+POINT_Y, data_sp)
  b<-fit.variogram(a, vgm("Sph"))
  g.dummy <- gstat(formula=z~x+y, locations=~x+y, dummy=T, beta=1, model=vgm(psill=b[2,2], range=b[2,3], model='Exp'), nmax=20)
  sim_dist_krant <- predict(g.dummy, newdata=xy, nsim=30)[,c(3:32)]
  #
  detach(train_neo_set2)
  #
}

repeats<-30
#creating empty data structures
AUC1<-vector(mode="numeric", repeats)
AUC2<-vector(mode="numeric", repeats)
AUC3<-vector(mode="numeric", repeats)
AUC4<-vector(mode="numeric", repeats)
AUC5<-vector(mode="numeric", repeats)
AUC6<-vector(mode="numeric", repeats)
AUC7<-vector(mode="numeric", repeats)
AUC8<-vector(mode="numeric", repeats)
AUC9<-vector(mode="numeric", repeats)
AUC10<-vector(mode="numeric", repeats)

Score1_1<-vector(mode="numeric", repeats)
Score2_1<-vector(mode="numeric", repeats)
Score3_1<-vector(mode="numeric", repeats)
Score4_1<-vector(mode="numeric", repeats)
Score5_1<-vector(mode="numeric", repeats)
Score6_1<-vector(mode="numeric", repeats)
Score7_1<-vector(mode="numeric", repeats)
Score8_1<-vector(mode="numeric", repeats)
Score9_1<-vector(mode="numeric", repeats)
Score10_1<-vector(mode="numeric", repeats)

Score1_0<-vector(mode="numeric", repeats)
Score2_0<-vector(mode="numeric", repeats)
Score3_0<-vector(mode="numeric", repeats)
Score4_0<-vector(mode="numeric", repeats)
Score5_0<-vector(mode="numeric", repeats)
Score6_0<-vector(mode="numeric", repeats)
Score7_0<-vector(mode="numeric", repeats)
Score8_0<-vector(mode="numeric", repeats)
Score9_0<-vector(mode="numeric", repeats)
Score10_0<-vector(mode="numeric", repeats)

Total_Score1<-vector(mode="numeric", repeats)
Total_Score2<-vector(mode="numeric", repeats)
Total_Score3<-vector(mode="numeric", repeats)
Total_Score4<-vector(mode="numeric", repeats)
Total_Score5<-vector(mode="numeric", repeats)
Total_Score6<-vector(mode="numeric", repeats)
Total_Score7<-vector(mode="numeric", repeats)
Total_Score8<-vector(mode="numeric", repeats)
Total_Score9<-vector(mode="numeric", repeats)
Total_Score10<-vector(mode="numeric", repeats)

c<-mean(as.numeric(train$neo)-1)
sim_train<-vector(mode="list", repeats)
#repeated stratified 10-fold cross-validation & GAM
for (a in seq(repeats)) {
  #creating training data
  sim_train[[a]]<-data.frame(neo=train$neo,dist_krant=sim_dist_krant[,a],pc2=sim_pc2[,a],
                             vand_top_p=sim_vand_top_p[,a],
                             tri_reg_r1=sim_tri_reg_r1[,a])
  folds <- createFolds(factor(sim_train[[a]]$neo), k = 10, list = FALSE,returnTrain = FALSE)
  
  fold1<-sim_train[[a]][which(folds==1),]
  fold2<-sim_train[[a]][which(folds==2),]
  fold3<-sim_train[[a]][which(folds==3),]
  fold4<-sim_train[[a]][which(folds==4),]
  fold5<-sim_train[[a]][which(folds==5),]
  fold6<-sim_train[[a]][which(folds==6),]
  fold7<-sim_train[[a]][which(folds==7),]
  fold8<-sim_train[[a]][which(folds==8),]
  fold9<-sim_train[[a]][which(folds==9),]
  fold10<-sim_train[[a]][which(folds==10),]
  
  train1<-sim_train[[a]][-which(folds==1),]
  train2<-sim_train[[a]][-which(folds==2),]
  train3<-sim_train[[a]][-which(folds==3),]
  train4<-sim_train[[a]][-which(folds==4),]
  train5<-sim_train[[a]][-which(folds==5),]
  train6<-sim_train[[a]][-which(folds==6),]
  train7<-sim_train[[a]][-which(folds==7),]
  train8<-sim_train[[a]][-which(folds==8),]
  train9<-sim_train[[a]][-which(folds==9),]
  train10<-sim_train[[a]][-which(folds==10),]
  
  gam_model<-gam(neo ~ s(dist_krant) + s(pc2) + s(vand_top_p) +
                   s(tri_reg_r1), family = binomial, select=TRUE, method="REML", data = train1)
  #AUROC
  pred<-predict(gam_model,fold1,type = "response")
  pred1 <- prediction(as.vector(pred),fold1[,1])
  peri <- performance(pred1, measure = "auc")
  AUC1[a]<-as.numeric(peri@y.values)
  #Beta family Proper scoring rule with alpha = 1 and beta=1/c - 1
  Score1_1[a]<-mean(calcscore(as.numeric(fold1[,1])-1~pred,fam="beta", param= c(1,-1+1/c),bounds = c(0,1))[which(fold1[,1]==1)])
  Score1_0[a]<-mean(calcscore(as.numeric(fold1[,1])-1~pred,fam="beta", param= c(1,-1+1/c),bounds = c(0,1))[which(fold1[,1]==0)])
  Total_Score1[a]<-mean(calcscore(as.numeric(fold1[,1])-1~pred,fam="beta", param= c(1,-1+1/c),bounds = c(0,1)))
  ###
  print((((a-1)*10+1)/300)*100)
  ###
  gam_model<-gam(neo ~ s(dist_krant) + s(pc2) + s(vand_top_p) +
                   s(tri_reg_r1), family = binomial, select=TRUE, method="REML", data = train2)
  #AUROC
  pred<-predict(gam_model,fold2,type = "response")
  pred1 <- prediction(as.vector(pred),fold2[,1])
  peri <- performance(pred1, measure = "auc")
  AUC2[a]<-as.numeric(peri@y.values)
  #Beta family Proper scoring rule with alpha = 1 and beta=1/c - 1
  Score2_1[a]<-mean(calcscore(as.numeric(fold2[,1])-1~pred,fam="beta", param= c(1,-1+1/c),bounds = c(0,1))[which(fold2[,1]==1)])
  Score2_0[a]<-mean(calcscore(as.numeric(fold2[,1])-1~pred,fam="beta", param= c(1,-1+1/c),bounds = c(0,1))[which(fold2[,1]==0)])
  Total_Score2[a]<-mean(calcscore(as.numeric(fold2[,1])-1~pred,fam="beta", param= c(1,-1+1/c),bounds = c(0,1)))
  ###
  print((((a-1)*10+2)/300)*100)
  ###
  gam_model<-gam(neo ~ s(dist_krant) + s(pc2) + s(vand_top_p) +
                   s(tri_reg_r1), family = binomial, select=TRUE, method="REML", data = train3)
  #AUROC
  pred<-predict(gam_model,fold3,type = "response")
  pred1 <- prediction(as.vector(pred),fold3[,1])
  peri <- performance(pred1, measure = "auc")
  AUC3[a]<-as.numeric(peri@y.values)
  #Beta family Proper scoring rule with alpha = 1 and beta=1/c - 1
  Score3_1[a]<-mean(calcscore(as.numeric(fold3[,1])-1~pred,fam="beta", param= c(1,-1+1/c),bounds = c(0,1))[which(fold3[,1]==1)])
  Score3_0[a]<-mean(calcscore(as.numeric(fold3[,1])-1~pred,fam="beta", param= c(1,-1+1/c),bounds = c(0,1))[which(fold3[,1]==0)])
  Total_Score3[a]<-mean(calcscore(as.numeric(fold3[,1])-1~pred,fam="beta", param= c(1,-1+1/c),bounds = c(0,1)))
  ###
  print((((a-1)*10+3)/300)*100)
  ###
  gam_model<-gam(neo ~ s(dist_krant) + s(pc2) + s(vand_top_p) +
                   s(tri_reg_r1), family = binomial, select=TRUE, method="REML", data = train4)
  #AUROC
  pred<-predict(gam_model,fold4,type = "response")
  pred1 <- prediction(as.vector(pred),fold4[,1])
  peri <- performance(pred1, measure = "auc")
  AUC4[a]<-as.numeric(peri@y.values)
  #Beta family Proper scoring rule with alpha = 1 and beta=1/c - 1
  Score4_1[a]<-mean(calcscore(as.numeric(fold4[,1])-1~pred,fam="beta", param= c(1,-1+1/c),bounds = c(0,1))[which(fold4[,1]==1)])
  Score4_0[a]<-mean(calcscore(as.numeric(fold4[,1])-1~pred,fam="beta", param= c(1,-1+1/c),bounds = c(0,1))[which(fold4[,1]==0)])
  Total_Score4[a]<-mean(calcscore(as.numeric(fold4[,1])-1~pred,fam="beta", param= c(1,-1+1/c),bounds = c(0,1)))
  ###
  print((((a-1)*10+4)/300)*100)
  ###
  gam_model<-gam(neo ~ s(dist_krant) + s(pc2) + s(vand_top_p) +
                   s(tri_reg_r1), family = binomial, select=TRUE, method="REML", data = train5)
  #AUROC
  pred<-predict(gam_model,fold5,type = "response")
  pred1 <- prediction(as.vector(pred),fold5[,1])
  peri <- performance(pred1, measure = "auc")
  AUC5[a]<-as.numeric(peri@y.values)
  #Beta family Proper scoring rule with alpha = 1 and beta=1/c - 1
  Score5_1[a]<-mean(calcscore(as.numeric(fold5[,1])-1~pred,fam="beta", param= c(1,-1+1/c),bounds = c(0,1))[which(fold5[,1]==1)])
  Score5_0[a]<-mean(calcscore(as.numeric(fold5[,1])-1~pred,fam="beta", param= c(1,-1+1/c),bounds = c(0,1))[which(fold5[,1]==0)])
  Total_Score5[a]<-mean(calcscore(as.numeric(fold5[,1])-1~pred,fam="beta", param= c(1,-1+1/c),bounds = c(0,1)))
  ###
  print((((a-1)*10+5)/300)*100)
  ###
  gam_model<-gam(neo ~ s(dist_krant) + s(pc2) + s(vand_top_p) +
                   s(tri_reg_r1), family = binomial, select=TRUE, method="REML", data = train6)
  #AUROC
  pred<-predict(gam_model,fold6,type = "response")
  pred1 <- prediction(as.vector(pred),fold6[,1])
  peri <- performance(pred1, measure = "auc")
  AUC6[a]<-as.numeric(peri@y.values)
  #Beta family Proper scoring rule with alpha = 1 and beta=1/c - 1
  Score6_1[a]<-mean(calcscore(as.numeric(fold6[,1])-1~pred,fam="beta", param= c(1,-1+1/c),bounds = c(0,1))[which(fold6[,1]==1)])
  Score6_0[a]<-mean(calcscore(as.numeric(fold6[,1])-1~pred,fam="beta", param= c(1,-1+1/c),bounds = c(0,1))[which(fold6[,1]==0)])
  Total_Score6[a]<-mean(calcscore(as.numeric(fold6[,1])-1~pred,fam="beta", param= c(1,-1+1/c),bounds = c(0,1)))
  ###
  print((((a-1)*10+6)/300)*100)
  ###
  gam_model<-gam(neo ~ s(dist_krant) + s(pc2) + s(vand_top_p) +
                   s(tri_reg_r1), family = binomial, select=TRUE, method="REML", data = train7)
  #AUROC
  pred<-predict(gam_model,fold7,type = "response")
  pred1 <- prediction(as.vector(pred),fold7[,1])
  peri <- performance(pred1, measure = "auc")
  AUC7[a]<-as.numeric(peri@y.values)
  #Beta family Proper scoring rule with alpha = 1 and beta=1/c - 1
  Score7_1[a]<-mean(calcscore(as.numeric(fold7[,1])-1~pred,fam="beta", param= c(1,-1+1/c),bounds = c(0,1))[which(fold7[,1]==1)])
  Score7_0[a]<-mean(calcscore(as.numeric(fold7[,1])-1~pred,fam="beta", param= c(1,-1+1/c),bounds = c(0,1))[which(fold7[,1]==0)])
  Total_Score7[a]<-mean(calcscore(as.numeric(fold7[,1])-1~pred,fam="beta", param= c(1,-1+1/c),bounds = c(0,1)))
  ###
  print((((a-1)*10+7)/300)*100)
  ###
  gam_model<-gam(neo ~ s(dist_krant) + s(pc2) + s(vand_top_p) +
                   s(tri_reg_r1),family = binomial, select=TRUE, method="REML", data = train8)
  #AUROC
  pred<-predict(gam_model,fold8,type = "response")
  pred1 <- prediction(as.vector(pred),fold8[,1])
  peri <- performance(pred1, measure = "auc")
  AUC8[a]<-as.numeric(peri@y.values)
  #Beta family Proper scoring rule with alpha = 1 and beta=1/c - 1
  Score8_1[a]<-mean(calcscore(as.numeric(fold8[,1])-1~pred,fam="beta", param= c(1,-1+1/c),bounds = c(0,1))[which(fold8[,1]==1)])
  Score8_0[a]<-mean(calcscore(as.numeric(fold8[,1])-1~pred,fam="beta", param= c(1,-1+1/c),bounds = c(0,1))[which(fold8[,1]==0)])
  Total_Score8[a]<-mean(calcscore(as.numeric(fold8[,1])-1~pred,fam="beta", param= c(1,-1+1/c),bounds = c(0,1)))
  ###
  print((((a-1)*10+8)/300)*100)
  ###
  gam_model<-gam(neo ~ s(dist_krant) + s(pc2) + s(vand_top_p) +
                   s(tri_reg_r1), family = binomial, select=TRUE, method="REML", data = train9)
  #AUROC
  pred<-predict(gam_model,fold9,type = "response")
  pred1 <- prediction(as.vector(pred),fold9[,1])
  peri <- performance(pred1, measure = "auc")
  AUC9[a]<-as.numeric(peri@y.values)
  #Beta family Proper scoring rule with alpha = 1 and beta=1/c - 1
  Score9_1[a]<-mean(calcscore(as.numeric(fold9[,1])-1~pred,fam="beta", param= c(1,-1+1/c),bounds = c(0,1))[which(fold9[,1]==1)])
  Score9_0[a]<-mean(calcscore(as.numeric(fold9[,1])-1~pred,fam="beta", param= c(1,-1+1/c),bounds = c(0,1))[which(fold9[,1]==0)])
  Total_Score9[a]<-mean(calcscore(as.numeric(fold9[,1])-1~pred,fam="beta", param= c(1,-1+1/c),bounds = c(0,1)))
  ###
  print((((a-1)*10+9)/300)*100)
  ###
  gam_model<-gam(neo ~ s(dist_krant) + s(pc2) + s(vand_top_p) +
                   s(tri_reg_r1), family = binomial, select=TRUE, method="REML", data = train10)
  #AUROC
  pred<-predict(gam_model,fold10,type = "response")
  pred1 <- prediction(as.vector(pred),fold10[,1])
  peri <- performance(pred1, measure = "auc")
  AUC10[a]<-as.numeric(peri@y.values)
  #Beta family Proper scoring rule with alpha = 1 and beta=1/c - 1
  Score10_1[a]<-mean(calcscore(as.numeric(fold10[,1])-1~pred,fam="beta", param= c(1,-1+1/c),bounds = c(0,1))[which(fold10[,1]==1)])
  Score10_0[a]<-mean(calcscore(as.numeric(fold10[,1])-1~pred,fam="beta", param= c(1,-1+1/c),bounds = c(0,1))[which(fold10[,1]==0)])
  Total_Score10[a]<-mean(calcscore(as.numeric(fold10[,1])-1~pred,fam="beta", param= c(1,-1+1/c),bounds = c(0,1)))
  ###
  print((((a-1)*10+10)/300)*100)
}
### proccessing, inspecting and obtaining relevant results
{
  AUC_GAM_neo_sim<-c(AUC1,AUC2,AUC3,AUC4,AUC5,AUC6,AUC8,AUC9,AUC10)
  Score_1_GAM_neo_sim<-c(Score1_1,Score2_1,Score3_1,Score4_1,Score5_1,Score6_1,
                         Score7_1,Score8_1,Score9_1,Score10_1)
  Score_0_GAM_neo_sim<-c(Score1_0,Score2_0,Score3_0,Score4_0,Score5_0,Score6_0,
                         Score7_0,Score8_0,Score9_0,Score10_0)
  Total_score_GAM_neo_sim<-c(Total_Score1,Total_Score2,Total_Score3,
                             Total_Score4,Total_Score5,Total_Score6,
                             Total_Score7,Total_Score8,Total_Score9,
                             Total_Score10)
}
#comparing results of null models with original data
{
  c(mean(AUC_GAM_neo_sim),mean(AUC_GAM_neo_org)) # 0.7557886 vs 0.8599881
  boxplot(AUC_GAM_neo_sim,AUC_GAM_neo_org)
  
  c(mean(Score_1_GAM_neo_sim),mean(Score_1_GAM_neo_org)) # 0.2714611 vs 0.1850581
  boxplot(Score_1_GAM_neo_sim,Score_1_GAM_neo_org)
  
  c(mean(Score_0_GAM_neo_sim),mean(Score_0_GAM_neo_org)) # 0.006654400 vs. 0.005661996
  boxplot(Score_0_GAM_neo_sim,Score_0_GAM_neo_org)
  
  c(mean(Total_score_GAM_neo_sim),mean(Total_score_GAM_neo_org)) # 0.01480764 vs. 0.01118553
  boxplot(Total_score_GAM_neo_sim,Total_score_GAM_neo_org)
}
