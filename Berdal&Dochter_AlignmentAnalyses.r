library(lme4); library(MCMCglmm); library(dplyr); library(gtools);library(RMark);library(ggplot2)
source("C:/Users/Dochtermann/Dropbox/Working/Projects/Berdal&Dochter(AdaptiveAlignment)/1_Data&Analysis/evo.fun.r")
#### Data importation and reshaping ####
behavdata <- read.csv("C:/Users/Dochtermann/Dropbox/Working/Projects/Berdal&Dochter(AdaptiveAlignment)/1_Data&Analysis/FullData2017_2.csv", header=TRUE)

#remove columns not being used
sub_behavdata <- behavdata[,-c(1,15:17,22:57)]
#remove individuals with unknown/conflicting IDs
sub_behavdata <- sub_behavdata[sub_behavdata$Animal_ID!="?",]

#reshape data
behav_wide <- reshape(sub_behavdata, idvar = c("Date","Order","Animal_ID","Shelter","Sex",
                                               "Dev_stage","Recap._no","Trials","Tagger","Observer","Mass"),
                      timevar = "Trial", direction = "wide")

RMarkData=read.table(file="C:/Users/Dochtermann/Dropbox/Working/Projects/Berdal&Dochter(AdaptiveAlignment)/1. Data & Analysis/RMarkString.txt",
                     header=TRUE,colClasses=c("character","numeric","numeric","numeric",
                                              "numeric","numeric","numeric","numeric","numeric"))

#load("C:/Users/Dochtermann/Dropbox/Working/Projects/Berdal&Dochter(AdaptiveAlignment)/1_Data&Analysis/Alignment_revision.RData")

#### Calculate variances for use in priors ####
#MCMCglmm can get stuck at zero so initialize versus expected proportion based on Bell et al. meta-analysis
#compare MCMC estiamtes to lmer estimates to make sure this doesn't bias estiamtes
among_prior <- diag(c(var(behav_wide$Distance_moved.OF, na.rm=T)*0.37,var(behav_wide$TimeNear.AGG, na.rm=T)*0.37,var(behav_wide$Mean_dist_AP.AP, na.rm=T)*0.37))
within_prior <- diag(c(var(behav_wide$Distance_moved.OF, na.rm=T)*.63,var(behav_wide$TimeNear.AGG, na.rm=T)*.63,var(behav_wide$Mean_dist_AP.AP, na.rm=T)*.63))

#### MCMCglmm model #####
NITT=12000;THIN=10;BURN=2000
Mult=500

prior_3<-list(G=list(G1=list(V=among_prior,nu=.002,
                             alpha.mu = rep(0,3), alpha.V = diag(3)*1000)),
              R=list(V=within_prior,nu=.002))

MCMC.pema<-(MCMCglmm(
  cbind((sqrt(Distance_moved.OF)),(TimeNear.AGG),(Mean_dist_AP.AP))~
    (trait-1)*Shelter+(trait-1)*Sex+(trait-1)*Dev_stage+(trait-1)*scale(Mass,scale=FALSE),
  random=~us(trait):Animal_ID,rcov=~us(trait):units,
  family=rep("gaussian",3), nitt=NITT*Mult,thin=THIN*Mult,burnin=BURN*Mult,
  verbose=FALSE, prior=prior_3,pr=TRUE, data=behav_wide))

summary(MCMC.pema)

#plot(MCMC.pema$VCV)

# Calculate Repeatbilities ####

#Activity
Vi.act= MCMC.pema$VCV[,1]
Vw.act= MCMC.pema$VCV[,10]
tau.act=Vi.act/(Vi.act+Vw.act)


#Aggression
Vi.agg= MCMC.pema$VCV[,5]
Vw.agg= MCMC.pema$VCV[,14]
tau.agg=Vi.agg/(Vi.agg+Vw.agg)

#Anti-pred
Vi.AP= MCMC.pema$VCV[,9]
Vw.AP= MCMC.pema$VCV[,18]
tau.AP=Vi.AP/(Vi.AP+Vw.AP)

median(tau.act)
HPDinterval(tau.act)

median(tau.agg)
HPDinterval(tau.agg)

median(tau.AP)
HPDinterval(tau.AP)
# same for 2 other traits

# Compare mcmcglmm variances to those estimated by lme4 ####
OF_Dat <- behav_wide[which(behav_wide$Distance_moved.OF!="NA"),]
of_mm <- lmer(Distance_moved.OF~Shelter+Sex+Dev_stage+
                scale(Mass,scale=FALSE)+(1|Animal_ID),
              data=OF_Dat)
of_IV <- VarCorr(of_mm)$Animal_ID[1]
of_WV <- sigma(of_mm)^2
of_rep <- of_IV/(of_IV+of_WV)

AGG_Dat <- behav_wide[which(behav_wide$TimeNear.AGG!="NA"),]
agg_mm <- lmer(TimeNear.AGG~Shelter+Sex+Dev_stage+
                 scale(Mass,scale=FALSE)+(1|Animal_ID),
               data=AGG_Dat)
agg_IV <- VarCorr(agg_mm)$Animal_ID[1]
agg_WV <- sigma(agg_mm)^2
agg_rep <- agg_IV/(agg_IV+agg_WV)


AP_Dat <- behav_wide[which(behav_wide$Mean_dist_AP.AP!="NA"),]
ap_mm <- lmer(Mean_dist_AP.AP~Shelter+Sex+Dev_stage+
                scale(Mass,scale=FALSE)+(1|Animal_ID),
              data=AP_Dat)
ap_IV <- VarCorr(ap_mm)$Animal_ID[1]
ap_WV <- sigma(ap_mm)^2
ap_rep <- ap_IV/(ap_IV+ap_WV)

#### Correlations and Covariances ####
# Correlation matrices #
#Among-individual
post.among.cor.mat<-matrix(posterior.mode(posterior.cor(MCMC.pema$VCV[,1:9])),3)
post.cor.among <- posterior.cor(MCMC.pema$VCV[,1:9])
post.med.among <- matrix(apply(post.cor.among,2,median),3)
post.med.among
HPDinterval(posterior.cor(MCMC.pema$VCV[,1:9])) # CI for among-individual correlations
#pMCMC -- proportion of posterior on one side of zero
sum(post.cor.among[,2]<0)
sum(post.cor.among[,3]<0)
sum(post.cor.among[,6]>0)

#Within-individual
post.within.cor.mat<-matrix(posterior.mode(posterior.cor(MCMC.pema$VCV[,10:18])),3)
post.cor.within <- posterior.cor(MCMC.pema$VCV[,10:18])
post.med.within <- matrix(apply(post.cor.within,2,median),3)
post.med.within
HPDinterval(posterior.cor(MCMC.pema$VCV[,10:18])) # CI for within-individual correlations
#pMCMC -- proportion of posterior on one side of zero
sum(post.cor.within[,2]>0)
sum(post.cor.within[,3]<0)
sum(post.cor.within[,6]<0)

# Covariance matrices #

post.among.cov.mat<-matrix(posterior.mode((MCMC.pema$VCV[,1:9])),3)
post.within.cov.mat<-matrix(posterior.mode((MCMC.pema$VCV[,10:18])),3)

#posterior distribution of covariances diff from zero
colSums(MCMC.pema$VCV[,c(1:9)]<0)/1000
colSums(MCMC.pema$VCV[,c(10:18)]<0)/1000

#### Evolvability ####
among_autonomy <- as.mcmc(apply(MCMC.pema$VCV[,1:9],1,evo.fun))
within_autonomy <- as.mcmc(apply(MCMC.pema$VCV[,10:18],1,evo.fun))
median(among_autonomy)
HPDinterval(among_autonomy)
median(within_autonomy)
HPDinterval(within_autonomy)

#summarize repeatability estimates  for reporting and plotting ####
df_plot=data.frame(Trait=factor(c("Exp","Agg","AP")),
                   lmer_est=c(of_rep,
                              agg_rep,
                              ap_rep),
                   mcmc_median=c(median(tau.act),
                                 median(tau.agg),
                                 median(tau.AP)),
                   low.ci=c(HPDinterval(tau.act)[1],
                            HPDinterval(tau.agg)[1],
                            HPDinterval(tau.AP)[1]),
                   up.ci=c(HPDinterval(tau.act)[2],
                           HPDinterval(tau.agg)[2],
                           HPDinterval(tau.AP)[2]))

df_plot

#Repeatability plot ####
#add lme4 estimate for comparison to multiresponse MCMCglmm estimates

par(mar=c(5, 6, 4, 2),pty='s')
plot(NA,ylim=c(0,1),xlim=c(0.5,3.5),xaxt="n",
     xlab="Behavior",ylab="Repeatability",
     cex.lab=2,cex.axis=1.25)
arrows(x0=c(0.95,1.95,2.95),x1=c(0.95,1.95,2.95),
       y0=df_plot$low.ci,y1=df_plot$up.ci,
       angle=90,code=3,length=0.1)
points(x=c(0.95,1.95,2.95),y=df_plot$mcmc_median,
       pch=21,cex=2.5,
       bg=c("#55c667ff","#238a8dff","#404788FF"))
points(x=c(1.05,2.05,3.05),y=df_plot$lmer_est,
       pch=23,cex=1.5,
       bg=c("#55c667ff","#238a8dff","#404788FF"))
axis(1,at=c(1,2,3),labels=df_plot$Trait,cex.axis=1.25)

#Calculate P matrix after controlling for fixed effects ####
#just I + E
P_mat <- MCMC.pema$VCV[,1:9]+MCMC.pema$VCV[,10:18]

#### Calculate selection gradients via RMark ####

# Set covariates for survival (phi) and recapture probability (p)
Phi_F <- list(formula=~Sex1+Sex2+Dev+Mass+BLUPS.Act+BLUPS.Agg+BLUPS.AP)
p_F <- list(formula=~Sex1+Sex2+Dev+Mass+BLUPS.Act+BLUPS.Agg+BLUPS.AP)

Phi_dot <- list(formula=~Sex1+Sex2+Dev+Mass)
p_dot <- list(formula=~Sex1+Sex2+Dev+Mass)

subMARK <- RMarkData[,-c(7:9)]

Betastore <- matrix(NA,1000,3)
AICcStore <- matrix(NA,1000,2)

for(i in 1:1000){
  
  # extract one slice of the individual BLUBS for activity, aggression, and AP
  BLUP_slice <- data.frame(MCMC.pema$Sol[i,-c(1:18)])  
  # Create column for animal ID
  BLUP_slice$animal=factor(rownames(BLUP_slice))
  # Add headers
  names(BLUP_slice)=c("BLUPS","animal")
  
  # Add column with traits
  BLUP_slice$Trait=factor(c(rep("Act",72),
                            rep("Agg",143-71),
                            rep("AP",216-144)),
                          levels = c("Act","Agg","AP"))
  
  # Remove text from animal ID
  x <- gregexpr("[0-9]+", BLUP_slice$animal)
  BLUP_slice$animal <-  unlist(regmatches(BLUP_slice$animal, x))
  
  # Rescape dataset so each individual has only one row
  slice_wide <- reshape(BLUP_slice, idvar = "animal",timevar = "Trait", direction = "wide")
  # Make animal ID numeric and in increasing order
  slice_wide <- slice_wide[order(as.numeric(slice_wide$animal)),]  
  slice_wide[,2:4] <- apply(slice_wide[,2:4],2,scale) #interested in spatial orientation so scale to compare with correlation matrices
  
  # Merge with RMark dataset
  newMARK <- cbind(subMARK,slice_wide)
  
  # Create a list with the data and its attributes
  dp=process.data(newMARK,model="CJS", begin.time = 1, time.intervals = c(15,2,3,2,2,3,2,2,
                                                                          3,3,1,3,2,41,7,3,4,
                                                                          3,4,10,3,3,7,4,5))
  
  # Create design data
  ddl_F=make.design.data(dp)
  
  # Run MARK model
  DotModel <-  mark(dp,ddl,model.parameters=list(Phi=Phi_dot,p=p_dot))
  FullModel <- mark(dp,ddl,model.parameters=list(Phi=Phi_F,p=p_F))
  
  Betastore[i,1] <- coef(FullModel)[6,1]
  Betastore[i,2] <- coef(FullModel)[7,1]
  Betastore[i,3] <- coef(FullModel)[8,1]
  
  AICcStore[i,1] <- FullModel$results$AICc
  AICcStore[i,2] <- DotModel$results$AICc
  
}
Betastore <- as.data.frame(Betastore)
names(Betastore) <- c("Act","AGG","AP")
dAICc <- AICcStore[,2]-AICcStore[,1]
hist(dAICc)
median(dAICc)

#### Among x Within matrix alignment ####
#angle calculation functions
vec.out<-function(x,traits){y<-matrix(x,traits)
eigen(y)$vectors[,1]}

vec.cor<-function(z1=z1,z2=z2){
  abs(sum(z1*z2) / ( sqrt(sum(z1 * z1)) * sqrt(sum(z2 * z2)) ) )}

vec.angle<-function(x){acos(x)*(180/pi)}

#Calculate angle for paired estimates
among.cor.mat<-posterior.cor(MCMC.pema$VCV[,1:9])
within.cor.mat<-posterior.cor(MCMC.pema$VCV[,10:18])
p.cor.mat <- posterior.cor(P_mat[,])

among.eig<-t(apply(among.cor.mat,1,vec.out,traits=3))
within.eig<-t(apply(within.cor.mat,1,vec.out,traits=3))
p.eig<-t(apply(p.cor.mat,1,vec.out,traits=3))

vec.cor.out<-matrix(NA,nrow=1000,ncol=4)
for(i in 1:1000){
  z1<-among.eig[i,]
  z2<-within.eig[i,]
  z3 <- p.eig[i,]
  z4 <- Betastore[i,]
  vec.cor.out[i,1]<-vec.cor(z1,z2)
  vec.cor.out[i,2]<-vec.cor(z1,z4)
  vec.cor.out[i,3]<-vec.cor(z4,z2)
  vec.cor.out[i,4]<-vec.cor(z4,z3)
}

#Calculate vector correlation
apply(vec.cor.out,2,median)
HPDinterval(as.mcmc(vec.cor.out))

#convert to angles
angle.out<-vec.angle(vec.cor.out)
(ANGLES <- apply(angle.out,2,median))
HPDinterval(as.mcmc(angle.out))

#z-based significance test of significant misalignment
null.corr <- 0.975

z_AW <- (atan(median(angle.out[,1]))-atan(null.corr))/sqrt(2/length(unique(newMARK$ID)))
2*pnorm(-abs(z_AW))

z_AB <- (atan(median(angle.out[,2]))-atan(null.corr))/sqrt(2/length(unique(newMARK$ID)))
2*pnorm(-abs(z_AB))

z_WB <- (atan(median(angle.out[,3]))-atan(null.corr))/sqrt(2/length(unique(newMARK$ID)))
2*pnorm(-abs(z_WB))

z_PB <- (atan(median(angle.out[,4]))-atan(null.corr))/sqrt(2/length(unique(newMARK$ID)))
2*pnorm(-abs(z_PB))
#(all significant)

#Across the posterior:
z_AW_post <- (atan((angle.out[,1]))-atan(null.corr))/sqrt(2/length(unique(newMARK$ID)))
min(z_AW_post)
sum(z_AW_post<0)/1000

z_AB_post <- (atan((angle.out[,2]))-atan(null.corr))/sqrt(2/length(unique(newMARK$ID)))
min(z_AB_post)
sum(z_AB_post<0)/1000

z_WB_post <- (atan((angle.out[,3]))-atan(null.corr))/sqrt(2/length(unique(newMARK$ID)))
min(z_WB_post)
sum(z_WB_post<0)/1000

z_PB_post <- (atan((angle.out[,4]))-atan(null.corr))/sqrt(2/length(unique(newMARK$ID)))
min(z_PB_post)
sum(z_PB_post<0)/1000
#all significant across posterior

#save.image("C:/Users/Dochtermann/Dropbox/Working/Projects/Berdal&Dochter(AdaptiveAlignment)/1_Data&Analysis/Alignment_revision.RData")


#### Plot correlogram ####
library(ellipse);library(RColorBrewer);library(paletteer)

new_cor <- post.med.among
new_cor2 <- post.med.within
diag(new_cor) <- 1
diag(new_cor2) <- 1
new_mat <- matrix(NA,3,3)
new_mat[upper.tri(new_mat)] <- new_cor[upper.tri(new_cor)]
new_mat[lower.tri(new_mat)] <- new_cor2[lower.tri(new_cor2)]
row.names(new_mat)<-c("Exp","Agg","AP")
colnames(new_mat)<-c("","","")#c("Exp","Agg","AP") #c("Exploration(Exp)","Aggression (Agg)","Anti-predator response (AP)")
#diag(new_mat)<--1
paletteer_c("ggthemes::Red-Blue Diverging", 100)
my_colors <- c("white","#79A9CFFF","#ED6156FF",
"#A90C38FF","white","#ED6156FF",
"#A90C38FF","#ECB9ADFF","white")

par(pty='m',mar=c(6.5, 4, 1, 2))
plotcorr(new_mat,col=my_colors, mar=c(1,1,1,1),cex.lab=2,xpd=TRUE)
lines(x=c(1,3),y=c(3,1),lwd=2)
text(x=1,y=2,cex=1.75,labels=substitute(paste(bold('0.20'))))
text(x=1,y=1,cex=1.75,labels=substitute(paste(bold('-0.16'))))
text(x=2,y=1,cex=1.75,labels=substitute(paste(bold('-0.13'))))
text(x=2,y=3,cex=1.75,labels=substitute(paste(bold('-0.46'))),col="white")
text(x=3,y=3,cex=1.75,labels=substitute(paste(bold('-0.44'))),col="white")
text(x=3,y=2,cex=1.75,labels=substitute(paste(bold('-0.01'))))
text(x=2,y=.075,cex=2.5,labels="Behavior")

text(x=1,y=0.35,cex=2,labels="Exp")
text(x=2,y=0.35,cex=2,labels="Agg")
text(x=3,y=0.35,cex=2,labels="AP")
lines(x=c(0.25,3.67),y=c(.5,.5))
lines(x=c(0.25,0.25),y=c(.5,3.55))
lines(x=c(3.67,3.67),y=c(.5,3.55))
lines(x=c(0.25,3.67),y=c(3.55,3.55))

