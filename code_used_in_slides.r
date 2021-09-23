## official ones from CRAN
if (!require(Morpho))
    install.packages("Morpho")
if (!require(car))
    install.packages("car")
if (!require(shapes))
    install.packages("shapes")

## TRANSFORMS
## Load packages and data
require(Rvcg);require(rgl);require(Morpho) 
data(dummyhead);data(humface)

## Visualize it
wire3d(dummyhead.mesh,col=2)
spheres3d(dummyhead.lm,col=2,radius = 2)
wire3d(humface,col=3)
spheres3d(humface.lm,col=5,radius = 2)

## Load packages and data
require(Rvcg);require(rgl);require(Morpho) 
data(dummyhead);data(humface)

## Visualize it
wire3d(dummyhead.mesh,col=2)
spheres3d(dummyhead.lm,col=2,radius = 2)
wire3d(humface,col=3)
spheres3d(humface.lm,col=5,radius = 2)


## OPA
require(shapes);require(Morpho)

## rotate second landmark config to first 
OPA <- rotonto(gorf.dat[,,1],gorf.dat[,,2]) 
deformGrid2d(OPA$X, OPA$Y,pch=19,cex1 = 2,cex2=2,
             ngrid=15, main="OPA")

## procSym 
require(Morpho); require(shapes)
gor.dat <- bindArr(gorf.dat,gorm.dat,along=3)
gorproc <- procSym(gor.dat)
gorproc


## align2procSym
procM <- procSym(gorm.dat)
## rotate newdata to existing GPA
newalign <-align2procSym(procM,newdata=gorf.dat)

## plot aligned data
deformGrid2d(procM$mshape,newalign[,,1],pch=19,
             cex1=2,cex2=2)
for (i in 1:dim(newalign)[3])
	points(newalign[,,i],cex=2,pch=21,
           col="orange")
points(procM$mshape,col="red",pch=19,cex=2)

## CVA
# create factors
groups <-as.factor(c(rep("female",30),rep("male",29)))
# perform CVA and test Mahalanobis distance
# between groups with permutation test by 10000 rounds)            
cvagor <- CVA(gorproc$orpdata,groups,rounds=10000,cv=T)
cvagor
classify(cvagor,cv=F)

## visualize a shape change from score -5 to 5:
CV1Vis <- restoreShapes(c(-5,5),cvagor$CVvis[,1],
                        cvagor$Grandm)

deformGrid2d(CV1Vis[,,1],CV1Vis[,,2],
             main="Shape Change of CV1",pch=19,
             cex1 = 2,cex2=2,ngrid=15)


## Compute CVA from first 10 PCs
cvagor1<-CVA(gorproc$PCscores[,1:10],groups,
             plot=FALSE,cv=T,rounds=10000)

## Compute shape changes
cv2PC <- cvagor1$CVvis[,1]+cvagor1$Grandm
cv2PC5 <- (c(-5,5))*rbind(cv2PC,cv2PC)
CV1aVis <- restoreShapes(cv2PC5,gorproc$PCs[,1:10],
                         gorproc$mshape)

## Visualize
deformGrid2d(CV1aVis[,,2],CV1aVis[,,1],
             main="Shape Change of CV1 from PC-Scores",
             pch=19,cex1 = 2, cex2=2,ngrid=15)

deformGrid2d(CV1aVis[,,2],CV1aVis[,,1],main="Shape Change of CV1 from PC-Scores",pch=19,cex1 = 2,
             cex2=2,ngrid=15)

## groupPCA
data(boneData)
proc <- procSym(boneLM)
pop_sex <- name2factor(boneLM, which=3:4)

gpca <- groupPCA(proc$orpdata, groups=pop_sex, 
                 rounds=10000, mc.cores=2)
gpca
grandmean <-gpca$Grandmean
## calculate landmarks from 1st and 2nd between-group PC 
## (-2 and +2 standard deviations)
                   
gpcavisPC1_2sd <- restoreShapes(c(-2,2)*sd(gpca$Scores[,1]),
                                 gpca$groupPCs[,1], grandmean) 
gpcavisPC2_2sd <- restoreShapes(c(-2,2)*sd(gpca$Scores[,2]), 
                                  gpca$groupPCs[,2], grandmean) 
deformGrid3d(gpcavisPC1_2sd[,,1], gpcavisPC1_2sd[,,2], ngrid = 0)
require(rgl)
## visualize grandmean mesh
 
grandm.mesh <- tps3d(skull_0144_ch_fe.mesh, 
                     boneLM[,,1],grandmean,threads=1)
wire3d(grandm.mesh, col="white")
spheres3d(grandmean, radius=0.005)

car::spm(gpca$Scores,groups=pop_sex,smooth=F)

## Typicalities
tp <- typprobClass(cvagor$CVscores,groups=groups)
tp
classify(tp, cv=F)
myoutlier <- which(tp$groupaffinCV == "none")
tp$probs[myoutlier,]
tp$probsCV[myoutlier,]

boxplot(cvagor$CVscores~groups)
points(2,cvagor$CVscores[myoutlier],pch=19,col="red",
       cex=2)

deformGrid2d(gorproc$rotated[,,myoutlier],gorproc$mshape,
             cex1=2,cex2=2,main="Outlier vs Mean",pch=19)

## Ordination on new data
data(boneData)

pop_sex <- name2factor(boneLM,which=3:4)
set.seed(42)
train <- sample(80,size=60)
trainSamp <- boneLM[,,train]
testSamp <- boneLM[,,-train]
procTrain <- procSym(trainSamp)
CVAtrain <- CVA(procTrain$PCscores[,1:10],
                groups=pop_sex[train], rounds=10000,cv=T)
CVAtrain

car::spm(CVAtrain$CVscores,groups=pop_sex[train],smooth=F)

## Align to GPA
test2GPA <- align2procSym(procTrain,testSamp)

## Get PC-scores of new data
test2GPA_PCs <- getPCscores(test2GPA,procTrain$PCs,procTrain$mshape)
## Feed the scores to the model
testCVAclass <- classify(CVAtrain,newdata=test2GPA_PCs[,1:10])
testCVAclass


## Permutation tests

## groupPCA
gpca$groupdists
gpca$probs

cvagor1$Dist

warpmovie2d(CV1aVis[,,2],CV1aVis[,,1],n=15,
            folder="gormovie",pch=20,cex=4,
            links=c(1,5,4:2,8:6,1),lwd=3)

## PLSR

## download the human data
if (!file.exists("PLSDemo.RData"))
download.file(url="https://bit.ly/3tSkYo3",
              destfile="PLSDemo.RData")

## load it into the workspace
load("PLSDemo.RData")

Todo_trunk_reslid_GPA <- procSym(Todo_trunk_reslid)

PLS_sapiens <- pls2B(Todo_trunk_reslid_GPA$orpdata[Thorax_indices,,1:30],
                     Todo_trunk_reslid_GPA$orpdata[Pelvis_indices,,1:30]
                      ,rounds=999,mc.cores=0,same.config = F)
PLS_sapiens

## plot 2nd latent dimension
Trunk_sapiens_sex <- name2factor(Todo_trunk_reslid[,,1:30],
                                 which=3)

## obtain latent variables and plot them
commShape_sapiens <- getPLSCommonShape(PLS_sapiens) 
plot(commShape_sapiens$XscoresScaled[,2],
    commShape_sapiens$YscoresScaled[,2],
    col=c("skyblue","dodgerblue3")
     [as.numeric(Trunk_sapiens_sex)],
    pch=c(16,17)[as.numeric(Trunk_sapiens_sex)],
    cex=2,    xlab="Block 1 (Thorax)",
    ylab="Block 2 (Pelvis)",
    xlim=c(-0.08,0.08),ylim=c(-0.10,0.10),
	main="2nd latent dimension")

## add legend
legend("topleft",inset=.01,legend=c("Females","Males"), 
      horiz=FALSE, cex=1.2, bty="n",pt.cex=1.7,         
      col=c("skyblue","dodgerblue3"),
      border=col_leg_sex,pch=c(16,17))

## Get shapes associated with +-2sd of the second latent dimension
pred_sapiens <- plsCoVarCommonShape(PLS_sapiens,2,2)
## Warp a reference specimen mesh to the shapes
ref_sapiens <- rbind(Trunk_data[Thorax_indices,,41], 
                     Trunk_data[Pelvis_indices,,41]) 
Torso_pred1 <- tps3d(Torso_Homo,ref_sapiens,
                     pred_sapiens[,,1])
Torso_pred2 <- tps3d(Torso_Homo,ref_sapiens,
                     pred_sapiens[,,2]) 

## warpmovie3d(Torso_pred1,Torso_pred2,n=15,
##             folder="ThoraxPelvisMovie",
##             xland=pred_sapiens[,,1],
##             yland=pred_sapiens[,,2])


ref <- Todo_trunk_reslid_GPA$orpdata[,,1]
pred <- predictPLSfromData(PLS_sapiens, y=ref[Pelvis_indices,],ncomp = 10) 
deformGrid3d(ref[Thorax_indices,],pred)
