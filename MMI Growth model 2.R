###############################################################################
# 20/02/2016
# Author: Alastair V Harry - alastair.harry@gmail.com
# Deptartment of Fisheries Western Australia & 
# Centre For Sustainable Tropical Fisheries & Aquaculture
# James Cook University, Townsville
#
# From Harry et al (2013) Age, growth, and reproductive biology of the spot-tail
# shark, Carcharhinus sorrah, and the Australian blacktip shark, C. tilstoni, 
# from the Great Barrier Reef World Heritage Area, north-eastern Australia. 
# Marine and Freshwater Research
#
#
# Work in progress
################################################################################

# Clear console and remove list objects
cat("\014");rm(list=ls())

# Specify libraries required for analysis
library(nlstools)
library(RCurl)
library(ggplot2)
library(dplyr)
library(propagate)
source("Models.R")

# Load example data from github
# If sss.verifypeer = TRUE doesn't work, set to false, but read this:
# http://ademar.name/blog/2006/04/curl-ssl-certificate-problem-v.html
data <- getURL("https://raw.githubusercontent.com/alharry/spot-tail/master/SSS.csv",ssl.verifypeer = FALSE)%>%
  textConnection()%>%
  read.csv(sep=",")%>%
  tbl_df()

# Skip to here and load your data if you don't want to use the 
# example dataset

# Specify sex: "m", "f" or "both"
sex="f"

# Prepare data for non-linear model fitting and plotting. Choose starting values
# and define extremes of plotting areas

# Specify an appropriate length at birth, L0, for use in 2 parameter models
Mini<-L0<-525

# Specify dimensions for plotting (used at the end of this script)
Max.Overall<-ceiling(max(na.omit(data$STL)))
Max.Plot<-ceiling(max(na.omit(data$AgeAgree)))

# Specify subset of data to be used (i.e. males or females or both) and keep
# an extra data frame with all the data (data1)
data1<-data
data<-if(sex=="both"){data}else{filter(data,Sex==sex)}
Max<-max(na.omit(data$STL)) 
Max.Age<-max(na.omit(data$AgeAgree))

# Define starting values for asymptotic growth models using Ford-Walford method
# for Linf and K. To estimate L0, a quadratic is fit to the mean length at
# age data, and the intercept extracted. This section I got from Derek Ogle's
# FishR. 
mean.age<-with(data,tapply(STL,round(AgeAgree),mean,na.rm=T)) # Mean STL at age
Lt1<-mean.age[2:length(mean.age)] # Length at time t+1
Lt<-mean.age[1:length(mean.age)-1] # Length at time t
pars<-lm(Lt1~Lt)$coef # Regression of Lt1 against Lt
start.k<- abs(-log(pars[2])) # Get starting value for k, Linf and L0
start.linf<-abs(pars[1]/(1-pars[2])) 
start.L0<-lm(mean.age~poly(as.numeric(names(mean.age)), 2, raw=TRUE))$coef[1]

# Remove all non-essential data from the data frame
data<-na.omit(dplyr::select(data,STL,AgeAgree))
names(data)<-c("STL","Age")
xy<-data
names(xy)<-c("y","x")

# Create a set of candidate models for MMI comparison. 
# Three commonly used non-linear growth models are listed below.
# The functions have each been solved (by someone else) in such a way that
# they explicitly include the y-intercept (L0, length at time zero) as a
# fitted model parameter, and differ from the more common form that explicitly
# includes the x-intercept (t0, hypothetical age at zero length). The L0 
# parameterisation makes more biological sense for sharks since they are born 
# live and fully formed. 
# In instances where you don't have a good data for very young individuals,
# it is probably best to use the 2 parameter versions, otherwise the size
# at birth will be overestimated by the model. 
# Constraining the model to have only 2 parameters instead of 3 causes a bias
# in the other parameters, so if you have good data for all sizes its best to
# use the 3 parameter version. 
# IMO it doesn't make much sense to include both 2 and 3 parameter versions of the
# same model in your set of candidate models.

# G1 - 3 parameter VB
mo<- VB3
P<-c(b1=start.linf[[1]],b2=start.L0[[1]],b3=start.k[[1]])
Test<-fit(mo,P,xy)


#G1 - 3 parameter VB
g1<-STL~abc2+(abc1-abc2)*(1-exp(-abc3*Age))
g1.Start<-list(abc1=start.linf,abc2=start.L0,abc3=start.k)
g1.name<-"VB3"
# G2 - 2 parameter VB
g2<- STL~ L0+(ab1-L0)*(1-exp(-ab2*Age))
g2.Start<-list(ab1=start.linf,ab2=start.k)
g2.name<-"VB2"
# G3 - 3 parameter Gompertz
g3<- STL~ abc2*exp(log(abc1/abc2)*(1-exp(-abc3*Age)))
g3.Start<-list(abc1=start.linf,abc2=start.L0,abc3=start.k)
g3.name<-"GOM3"
# G4 - 2 parameter Gompertz
g4<- STL~L0*exp(log(ab1/L0)*(1-exp(-ab2*Age)))
g4.Start<-list(ab1=start.linf,ab2=start.k)
g4.name<-"GOM2"
# G5 - 3 parameter Logistic model
g5<- STL~ (abc1*abc2*exp(abc3*Age))/(abc1+abc2*((exp(abc3*Age))-1))
g5.Start<-list(abc1=start.linf, abc2=start.L0, abc3=start.k)
g5.name<-"LOGI3"
# G6 - 2 parameter Logistic model
g6<- STL~ (ab1*L0*exp(ab2*Age))/(ab1+L0*((exp(ab2*Age))-1))
g6.Start<-list(ab1=start.linf,ab2=start.k)
g6.name<-"LOGI2"


# Select the list of candidate models to compare.
Models<-list(g1,g2,g3,g4,g5,g6)
Start<-list(g1.Start,g2.Start,g3.Start,g4.Start,g5.Start,g6.Start)
Mod.names<-c(g1.name,g2.name,g3.name,g4.name,g5.name,g6.name)
Mods<-length(Models)
Results<-list()

# Fit growth models, outputs are stored as a list object called Results. Models
# are fit using nonlinear least squares 
Results[[1]]<-nls(Models[[1]],data=data,start=Start[[1]],algorithm="port")
Results[[2]]<-nls(Models[[2]],data=data,start=Start[[2]],algorithm="port")
Results[[3]]<-nls(Models[[3]],data=data,start=Start[[3]],algorithm="port")
Results[[4]]<-nls(Models[[4]],data=data,start=Start[[4]],algorithm="port")
Results[[5]]<-nls(Models[[5]],data=data,start=Start[[5]],algorithm="port")
Results[[6]]<-nls(Models[[6]],data=data,start=Start[[6]],algorithm="port")

# Put best fit paramters into a matrix
par.matrix<-data.frame(matrix(NA,Mods,3))
names(par.matrix)<-c('abc1','abc2','abc3')
for(i in 1:Mods){
  npars=length(coef(Results[[i]]))
  for(j in 1:npars)
    par.matrix[i,j]<-coef(Results[[i]])[j]
}

# Look at diagnostic plots
for(i in 1:Mods){plot(nlsResiduals(Results[[i]]))}


# Run the multi-model inference comparison. The following section calculates
# AIC values, AIC differences, AIC weights and residual standard errors
AIC.vals<-rep(NA,Mods)
Residual.Standard.Error<-rep(NA,Mods)
for(s in 1:Mods){AIC.vals[s]<-AIC(Results[[s]])}
AIC.dif<-AIC.vals-min(AIC.vals)
AIC.weight= round((exp(-0.5*AIC.dif))/(sum(exp(-0.5*AIC.dif))),4)*100
for(t in 1:Mods){Residual.Standard.Error[t]<-summary(Results[[t]])$sigma}

# Store output of MMI analysis. This is the important bit. Beyond this section
# its all plotting and confidence intervals, which may have a tendency not 
# to work so well depending on the sparsity of your data. It's somewhat
# 'optional' beyond here anyway, since many people don't present
# confidence intervals with their growth models anyway. 
MMI.Analysis<-data.frame(AIC.vals,AIC.dif,AIC.weight,Residual.Standard.Error,
                         par.matrix, row.names=Mod.names)
names(MMI.Analysis)[1:4]<-c("AIC", "AIC dif", "w","RSE")

# Calculate bootstrap confidence intervals. First step is to bootstrap the data. 
# The nlsBoot function takes the data and randomly re-samples it a certain
# number of times. 'Iter' is the number of interations bootstraps to be run. 
# I usually do 10000 in the final analysis but this will take a while.
# The actual confidence intervals can be extracted from the 'Bootstrapping'
# list e.g. to extract CI for LOGI3 type Bootstrapping[[3]]$bootCI
Iter=1000
Bootstrapping<-list()
for(i in 1:length(Models)){Bootstrapping[[i]]<-nlsBoot(Results[[i]],Iter)}

# Create a matrix with 95% confidence intervals
Conf.Ints<-data.frame(matrix(NA,length(Models),6),row.names=Mod.names)
names(Conf.Ints)[1:length(Conf.Ints[1,])]<-c("2.5%","97.5%")
for(i in 1:length(Models)){
  npars=nrow(Bootstrapping[[i]]$bootCI)
  for(j in 1:(npars*2)){Conf.Ints[i,j]<-c(t(Bootstrapping[[i]]$bootCI[,-1]))[j]}}

# Store all necessary outputs from analysis together
MMI.Analysis<-cbind(MMI.Analysis,Conf.Ints)

# Round to four significant figures
MMI.Analysis<-signif(MMI.Analysis,4)

# Specify values of age that you want to predict length over
x.deter<-seq(0,Max.Age,0.1)

# Number of iterations is now the number of successful bootstraps (i.e.
# original Iter minus those bootstraps that didn't converge)

Iter<-nrow(Bootstrapping[[which(AIC.dif==0)]]$coefboot)

# This bit calculates the actual confidence intervals and prediction intervals
# for the 'best fit' model
for(k in which(AIC.dif==0)){
  plotmatrix<-matrix(NA,ncol=length(x.deter),nrow=Iter)
  
  Npars<-length(Start[[which(AIC.dif==0)]])
  
  for(l in (1:Iter)){
    abc1<-as.numeric(Bootstrapping[[k]]$coef[l,1])
    abc2<-if(Npars==3){as.numeric(Bootstrapping[[k]]$coef[l,2])}else{L0}
    abc3<-if(Npars==3){as.numeric(Bootstrapping[[k]]$coef[l,3])}
    else{as.numeric(Bootstrapping[[k]]$coef[l,2])}
    # Quick add on to do confidence intervals in 2 parameter models
    ab1<-abc1
    ab2<-abc3
    Age=x.deter
    
    plotmatrix[l,]<-formula(formula(Results[[k]])[[3]])
  }
  
  RSE<-Bootstrapping[[k]]$rse
  boundsmatrix<-matrix(NA,nrow=4,ncol=length(x.deter))
  
  for(m in 1:length(x.deter)){
    boundsmatrix[1,m]<-quantile(plotmatrix[,m],0.025)
    boundsmatrix[2,m]<-quantile(plotmatrix[,m],0.975)
    boundsmatrix[3,m]<-quantile((plotmatrix-RSE)[,m],0.025)
    boundsmatrix[4,m]<-quantile((plotmatrix+RSE)[,m],0.975)
  }
}

# Save outputs
y.deter<-predict(Results[[which(AIC.dif==0)]],list(Age=x.deter),
                 type="response")
Output<-data.frame(Age=x.deter,TL=y.deter,clower=boundsmatrix[1,],cupper=boundsmatrix[2,],
                   plower=boundsmatrix[3,],pupper=boundsmatrix[4,])
Output<-Test$confInts
# Plot data
p<-ggplot(Output,aes(x=x,y=Mean))+geom_line()+
  geom_line(aes(y=cUpper),linetype="dashed")+geom_line(aes(y=cLower),linetype="dashed")+
  geom_line(aes(y=pUpper),linetype="dotted")+geom_line(aes(y=pLower),linetype="dotted")+
  geom_point(data=Test$data,aes(x=x,y=y))+xlab("Age (Years)")+ylab("TL (cm)")+
  theme(panel.background=element_blank())+theme_classic()+xlim(0,Max.Plot)+ylim(Mini*0.9,Max.Overall*1.1)

# Save analysis
#write.csv(MMI.Analysis, "MMI.Analysis.csv")
Results<-list(MMI.Analysis,Output)

