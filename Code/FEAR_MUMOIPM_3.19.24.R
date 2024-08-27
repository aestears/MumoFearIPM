## Shapefile reading and combining for FEAR
## Add data files

install.packages("devtools")
library(tidyverse)
library(devtools)
install_github('aestears/PlantTracker')
library(plantTracker)
library(tidyverse)
library(devtools)
library(MuMIn)
library(sjPlot)
library(car)
library(data.table)
library(lme4)
library(regclass)
library(optimx)
library(lmerTest)
install_version("IPMpack", version = "2.1", repos = "http://cran.us.r-project.org")
library(IPMpack)
library(zoo)
library(ggeffects)
library("RColorBrewer")
library(cAIC4)
library(DescTools)
library(piecewiseSEM)

#The section below uses plant tracker to extract data for FEAR and MUMO with a 20cm DD buffer which can
#Include intra or interspecific density depending on what is wanted.  Adapted from Dave Atkins' code.

## Create new home directory
setwd("/Users/rmm528/Desktop/shapefile_work-main")
dir_home = getwd()

## (3) Generate a list of all of your files
path_list = list.files(path = paste0(dir_home,"/","shapefiles_2002_2022"), recursive = T, pattern = ".shp")
path_list = path_list[!grepl("xml", path_list)] ## drops duplicates that are unneeded

shps = lapply(paste0(dir_home,"/shapefiles_2002_2022/",path_list), sf::st_read, type = 3) ## read in shps from list
shps2 = lapply(shps, sf::st_cast, to="GEOMETRY") ## change the class to GEOMETRY
shp = sf::st_as_sf(data.table::rbindlist(shps2)) ## Bind together

unique(shp$z_Year) ## check for no 2008
unique(shp$species) ## see species list
shp = shp[shp$species != "All Other Species",] ## remove the all other category

shp = sf::st_cast(shp, to="MULTIPOLYGON") ## back to multipolygon

## run following line before coercing valid geometries to get the reason why each geometry is invalid
sf::st_is_valid(shp, reason = T)

shp_valid = sf::st_make_valid(shp) ## removes issues with any invalid geometries.

shp_valid$z_Year = as.integer(shp_valid$z_Year)
shp_valid$Plot = as.character(shp_valid$Plot)

## Loop runs each quadrat individually. Potential memory leak freezes computer if run all together.

PlotList_genet = list() ## empty list for trackspp output
i=1
for(i in 1:length(unique(shp_valid$Plot))){
  dat_i = shp_valid[shp_valid$Plot == unique(shp_valid$Plot)[i],]
  inv_i = list(c(2002:2007,2009:2022))
  names(inv_i) = unique(dat_i$Plot)
  PlotList_genet[[i]] = trackSpp(dat_i, 
                                 inv = inv_i,
                                 dorm = 1,
                                 buff = 0.05,
                                 clonal = T,
                                 buffGenet = 0.01,
                                 species = "species",
                                 site = "Site",
                                 quad = "Plot",
                                 year = "z_Year",
                                 geometry = "geometry")
  print(c("### Plot Number", i, "out of", length(unique(shp_valid$Plot))))
  
}

PlotList_genet = lapply(PlotList_genet,sf::st_cast, to="MULTIPOLYGON") ## change to multipolygon
pt_out_genet_plot = sf::st_as_sf(data.table::rbindlist(PlotList_genet, use.names = T)) ## bind together

saveRDS(pt_out_genet_plot, file = "pt_out_genet_plot_6.2.23.RDS") ## save genet

####intraspecific competition 
pt_out_genet<-readRDS(file.choose())
pt_out_genet = pt_out_genet_plot
#pt_out_ramet = readRDS("pt_out_ramet.RDS")


sf::st_crs(pt_out_genet) = NA ## removes coordinate reference system from files. IMPORTANT
## 5,10,15,20 cm buffers
i=1
for(i in 1:length(unique(pt_out_genet$Plot))){
  dat_i = pt_out_genet[pt_out_genet$Plot == unique(pt_out_genet$Plot)[i],]
  print(c("### Quad Number",i,"out of",length(unique(pt_out_genet$Plot))))
  GN_genettemp <<- getNeighbors(dat = dat_i , buff = 0.2 , method = "area", compType = "oneSpp", output = "summed",
                                trackID = "trackID",
                                species = "species",
                                quad = "Plot",
                                year= "z_Year",
                                site = "Site",
                                geometry = "geometry")
  if(i == 1){
    GN_genet_20cm_intra <<- GN_genettemp
    print("first")
  } else {
    GN_genet_20cm_intra <<- rbind(GN_genet_20cm_intra, GN_genettemp)
    print("working")
  }
}

saveRDS(GN_genet_20cm_intra, file = "GN_genet_20cm_intra.8.21.RDS") ## save genet

data_out = GN_genet_20cm

data_out$neighbors_area_20cm = GN_genet_20cm$neighbors_area
data_out$nBuff_area_20cm = GN_genet_20cm$nBuff_area

saveRDS(data_out, file = "gn_data.RDS") ## save final plant tracker data output


#The below section extracts desired environmental variables from data supplied by 
#Jeff Jenness.  These data are extracted and interpolated from PRISM.
################################
#Create environmental variables

#Get appropriate environmental data
HWB.clim<-read.csv("HWB.clim.csv")#Climate data provided by Jeff J.
HWB.clim<-setDT(HWB.clim)
HWB.clim.02.22<-subset(HWB.clim, z_Year>2001)#pull appropriate years.  This dataset ends in 2022.

#Create t-1 lagged climate data
HWB.t.one.clim<-subset(HWB.clim, z_Year>2000)
HWB.t.one.clim<-subset(HWB.t.one.clim, z_Year<2022)
names(HWB.t.one.clim) <- c('Key', 'Site', 'Plot', 'Month', 'z_Year','tminusPrecip', 
                           'tminusMaxT', 'tminusMeanT', 'tminusMinT', 'tminusMaxVPD', 'tminusMinVPD')
HWB.t.one.clim$z_Year<-HWB.t.one.clim[,5]+1#adds one to t-1 year to properly match years in climate data

HWB.full.clim<- merge(HWB.clim.02.22,HWB.t.one.clim,by=c("Site","Plot", "Month", "z_Year"))
HWB.full.clim<-setDT(HWB.full.clim)
HWB.annual.clim<-HWB.full.clim[ , .(Precip=sum(Precip),
                                    MaxT=max(MaxT), 
                                    MeanT=mean(MeanT), 
                                    MinT=min(MinT), 
                                    MaxVPD=max(MaxVPD), 
                                    MinVPD=min(MinVPD), 
                                    tminusPrecip=sum(tminusPrecip),
                                    tminusMaxT=max(tminusMaxT), 
                                    tminusMeanT=mean(tminusMeanT),
                                    tminusMinT=min(tminusMinT), 
                                    tminusMaxVPD=max(tminusMaxVPD), 
                                    tminusMinVPD=min(tminusMinVPD)), 
                                by = .(Site,Plot, z_Year)]
HWB.annual.clim$MeanVPD<- rowMeans(HWB.annual.clim[,c('MinVPD', 'MaxVPD')], na.rm=TRUE)
HWB.annual.clim$tminusMeanVPD<-rowMeans(HWB.annual.clim[,c('tminusMinVPD', 'tminusMaxVPD')], na.rm=TRUE)
HWB.annual.clim<-as.data.frame(HWB.annual.clim)

clim.labels<-HWB.annual.clim[,c(1:3)]
clim.only.scale<-HWB.annual.clim[,c("Precip","MaxT","MeanT","MinT","MaxVPD","MinVPD",
                                    "tminusPrecip", "tminusMaxT","tminusMeanT",
                                    "tminusMinT","tminusMaxVPD","tminusMinVPD", "MeanVPD", "tminusMeanVPD")]
clim.only.scale<-scale(clim.only.scale, center = FALSE, scale = apply(clim.only.scale, 2, sd, na.rm = T))
clim.only.scale<-as.data.frame(clim.only.scale)

HWB.annual.scale.<-cbind(clim.labels, clim.only.scale)#Scale, but do not center, climate data.



#The below section uses outputs from Plant Tracker to created vital rates regressions and 
#construct an IPM for both FEAR and MUMO
################################
#Working with FEAR and MUMO
#data_out.all.inter<-read_rds("GN_genet_20cm_inter.8.21.rds")# min year 2002, #max year is 2022.  the neighbors_area column is interspecific competition within the 20cm buffer.

data_out.all.intra<-read_rds("GN_genet_20cm_intra.8.21.rds")#the neighbors_area column is intraspecific competition within the 20cm buffer

#extract FEAR 
FEAR.data.intra<-data_out.all.intra[ which(data_out.all.intra$species=='Festuca arizonica'),]
FEAR.data.intra$geometry <- NULL

FEAR.IPM.data.intra<-merge(FEAR.data.intra, HWB.annual.scale, by=c("Site", "Plot", "z_Year"))#Merges plantTracker output with climate data.  Drops incomplete cases (e.g., plots that weren't measured in every year from 2002-2022)
FEAR.IPM.data.intra$Plot<-as.factor(FEAR.IPM.data.intra$Plot)
FEAR.IPM.data.intra$z_Year<-as.factor(FEAR.IPM.data.intra$z_Year)
FEAR.IPM.data.intra$basalArea_genet<-log(FEAR.IPM.data.intra$basalArea_genet)#log transform size data
FEAR.IPM.data.intra$neighbors_area.1<-NULL#drop this, this is non-intraspecific cover.  #confuse.
FEAR.IPM.data.intra$neighbors_area<-scale(FEAR.IPM.data.intra$neighbors_area, center = FALSE, scale = T )#scale density dependence to match other ind. var.

count.out<-count(FEAR.IPM.data.intra,Site,Plot,trackID)#Gives count of unique individuals.  count is row total
count.plot<-count(FEAR.IPM.data.intra, Site, Plot)#Give count of plots, n=row total


#Below section runs vital rates regressions and collects parameters for IPM
#########FEAR Parameter List for IPM

fear.params=data.frame (
  surv.int=NA,
  surv.slope=NA,
  surv.dd_beta=NA,
  surv.precip_beta=NA,
  surv.temp_beta=NA,
  surv.VPD_beta=NA,
  surv.tminusPrecip_beta=NA,
  surv.tminusMeanT_beta=NA,
  surv.tminusMeanVPD_beta=NA,
  #####
  growth.int=NA,
  growth.slope=NA,
  growth.sd=NA,
  growth.dd_beta=NA,
  growth.precip_beta=NA,
  growth.temp_beta=NA,
  growth.VPD_beta=NA,
  growth.tminusPrecip_beta=NA,
  growth.tminusMeanT_beta=NA,
  growth.tminusMeanVPD_beta=NA,
  
  #####
  flwr.int=NA,
  flwr.slope=NA,
  #####
  seed.int=NA,
  seed.slope=NA,
  #####
  recruit.size.mean=NA,
  recruit.size.sd=NA,
  #####
  establishment.prob=NA
)


######FEAR Survival Model
######I had run inter vs. intra in past iterations, can add that back if we think it is useful.

surv.fear.intra<- glmer(survives_tplus1~basalArea_genet  + 
                          neighbors_area + Precip+ MeanT +MeanVPD+
                          tminusPrecip+tminusMeanT+tminusMeanVPD+
                          (1|Plot)+(1|z_Year),
                        family=binomial,
                        data=FEAR.IPM.data.intra, control = glmerControl(
                         optimizer ='optimx', optCtrl=list(method='nlminb')))

vif(surv.fear.intra)# check variance inflation

summary(surv.fear.intra)#results sum
rsquared(surv.fear.intra)#extract marginal and conditional R2

#store betas in params list
fear.params$surv.int=fixef(surv.fear.intra)[1]
fear.params$surv.slope=fixef(surv.fear.intra)[2]
fear.params$surv.dd_beta=fixef(surv.fear.intra)[3]
fear.params$surv.precip_beta=fixef(surv.fear.intra)[4]
fear.params$surv.temp_beta=fixef(surv.fear.intra)[5]
fear.params$surv.VPD_beta=fixef(surv.fear.intra)[6]
fear.params$surv.tminusPrecip_beta=fixef(surv.fear.intra)[7]
fear.params$surv.tminusMeanT_beta=fixef(surv.fear.intra)[8]
fear.params$surv.tminusMeanVPD_beta=fixef(surv.fear.intra)[9]

#Plot
#set plotting theme
set_theme(base= theme_bw())#My prefered theme

#Plot against size at T
surv.fear.predict<-ggpredict(surv.fear.intra, terms="basalArea_genet[all]")

surv.fear.size<-ggplot(surv.fear.predict, aes(x, predicted))+
  geom_line(color="turquoise3")+
  geom_ribbon(aes(ymin=conf.low, ymax=conf.high, fill="band") ,alpha=0.3)+
  scale_color_manual("",values="turquoise3")+
  scale_fill_manual("", values="turquoise3")+
  xlab(bquote("log(Size)"[t]))+
  ylab("Probability of Survival")+
  ylim(0,1)+
  theme_bw()+
  theme(legend.position="none",
        axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank())
surv.fear.size

#Plot against significant environmental variables
surv.fear.predict<-ggpredict(surv.fear.intra, terms="tminusPrecip[all]")

surv.fear.precip<-ggplot(surv.fear.predict, aes(x, predicted))+
  geom_line(color="turquoise3")+
  geom_ribbon(aes(ymin=conf.low, ymax=conf.high, fill="band") ,alpha=0.3)+
                scale_color_manual("",values="turquoise3")+
                scale_fill_manual("", values="turquoise3")+
  xlab(bquote("Scaled Precipitation"[t-1]))+
  ylab("Probability of Survival")+
  ylim(0,1)+
  theme_bw()+
  theme(legend.position="none",
    axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank())
surv.fear.precip

######FEAR Growth Model
FEAR.IPM.data.intra<-FEAR.IPM.data.intra[-which(is.na(FEAR.IPM.data.intra$size_tplus1)),]#removes rows where the individual dies at t+1.  Obviously no growth, and they break the regression
FEAR.IPM.data.intra$size_tplus1<-log(FEAR.IPM.data.intra$size_tplus1)#put size on correct scale

growth.fear.intra<- lmer(size_tplus1~basalArea_genet  +
                     neighbors_area +
                            Precip+
                            MeanT+
                            MeanVPD+
                     tminusPrecip+tminusMeanT+tminusMeanVPD+
                            (1|Plot)+(1|z_Year),
                       data=FEAR.IPM.data.intra)

vif(growth.fear.intra)#Check VIF

summary(growth.fear.intra)
rsquared(growth.fear.intra)#extract marginal and conditional R2

#add to params list
fear.params$growth.int=fixef(growth.fear.intra)[1]
fear.params$growth.slope=fixef(growth.fear.intra)[2]
fear.params$growth.dd_beta=fixef(growth.fear.intra)[3]
fear.params$growth.precip_beta=fixef(growth.fear.intra)[4]
fear.params$growth.temp_beta=fixef(growth.fear.intra)[5]
fear.params$growth.VPD_beta=fixef(growth.fear.intra)[6]
fear.params$growth.tminusPrecip_beta=fixef(growth.fear.intra)[7]
fear.params$growth.tminusMeanT_beta=fixef(growth.fear.intra)[8]
fear.params$growth.tminusMeanVPD_beta=fixef(growth.fear.intra)[9]
fear.params$growth.sd=sd(resid(growth.fear.intra))

#######Plot
growth.fear.size<-ggpredict(growth.fear.intra, terms="basalArea_genet[all]")

growth.size<-ggplot(growth.fear.size, aes(x, predicted))+
  geom_line(color="turquoise3")+
  geom_ribbon(aes(ymin=conf.low, ymax=conf.high, fill="band") ,alpha=0.3)+
  scale_color_manual("",values="turquoise3")+
  scale_fill_manual("", values="turquoise3")+
  xlab(bquote("log(size)"[t]))+
  ylab(bquote("log(Size)"[t+1]))+
  ylim(-13,0)+
  theme_bw()+
  theme(legend.position="none",
        axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank())
growth.size


growth.fear.precip<-ggpredict(growth.fear.intra, terms="Precip[all]")
growth.precip<-ggplot(growth.fear.precip, aes(x, predicted))+
  geom_line(color="turquoise3")+
  geom_ribbon(aes(ymin=conf.low, ymax=conf.high, fill="band") ,alpha=0.3)+
  scale_color_manual("",values="turquoise3")+
  scale_fill_manual("", values="turquoise3")+
  xlab("Scaled Precipitation")+
  ylab(bquote("log(Size)"[t+1]))+
  ylim(-8,0)+
  theme_bw()+
  theme(legend.position="none",
        axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank())
growth.precip

##########tables of FEAR survival and growth model outputs
#These can be copied and pasted directly into the manuscript
#and additional data added (e.g. R^2)
##########################################################
tab.surv.fear<-tbl_regression(surv.fear.intra)
tab.surv.fear<-bold_p(tab.surv.fear, t=0.05, q=F)
tab.surv.fear<-add_glance_table(include = c(nobs, r.squared))
tab.grow.fear<-tbl_regression(growth.fear.intra)
tab.grow.fear<-bold_p(tab.grow.fear, t=0.05, q=F)
tab.fear<-tbl_merge(tbls=list(tab.surv.fear, tab.grow.fear), 
                    tab_spanner=c("**Survival**", "**Growth**"))
tab.fear

###################
### FEAR Fecundity
#Fully size-based, as we don't have multiple years of fec. data to relate to climate

fec.FEAR = read.csv(file.choose(), header=T)#exract fecundity data
fec.FEAR$basalArea_genet = log((fec.FEAR$basal.area.cm2/10000))#put on same size scale
fec.FEAR = fec.FEAR[!is.na(fec.FEAR$basalArea_genet),]#remove NAs
glume.FEAR = read.csv(file.choose(), header=T)#get glumes
glume.FEAR$basalArea_genet = log((glume.FEAR$BA.cm2/10000))#same size scale
glume.FEAR$basalArea_genet<-as.numeric(glume.FEAR$basalArea_genet)
fec.FEAR$flwr = ifelse(fec.FEAR$numb.stalks>0, 1, 0)

#Flower by size regression
flower.FEAR.linear = glm(flwr ~ basalArea_genet, data=fec.FEAR, family = "binomial")
summary(flower.FEAR.linear)
fear.params$flwr.int=flower.FEAR.linear$coefficients[1]
fear.params$flwr.slope=flower.FEAR.linear$coefficients[2]

#Function to convert flowering to probability, for use in reproduction kernel
p_bz <- function(basalArea_genet, params)
{
  linear.p <- params$flwr.int + params$flwr.slope * basalArea_genet      # linear predictor
  p <- 1/(1+exp(-linear.p))                                # logistic transformation to probability
  return(p)
}

#seeds by size regression
seeds.FEAR = glm(sum.glumes ~ basalArea_genet, data=glume.FEAR, family = "poisson")
summary(seeds.FEAR)

fear.params$seed.int=seeds.FEAR$coefficients[1]
fear.params$seed.slope=seeds.FEAR$coefficients[2]

newdat=data.frame(FEAR.IPM.data.intra$basalArea_genet)
names(newdat)[1]="basalArea_genet"

FEAR.IPM.data.intra$Seeds<-round(predict(seeds.FEAR, newdata=newdat, type="response"))

fear.params$recruit.size.mean=mean(FEAR.IPM.data.intra$basalArea_genet[FEAR.IPM.data.intra$recruit==1], na.rm =TRUE)
fear.params$recruit.size.sd=sd(FEAR.IPM.data.intra$basalArea_genet[FEAR.IPM.data.intra$recruit==1], na.rm =TRUE)

fear.params$establishment.prob=sum(FEAR.IPM.data.intra$recruit, na.rm = TRUE)/sum(FEAR.IPM.data.intra$Seeds,na.rm=TRUE)


#########Calculating IPMS
#used only statistically significant environmental variables in construction
#########################

s.x.fear=function(x,fear.params) {
  u=exp(fear.params$surv.int + 
          fear.params$surv.slope*x + 
          fear.params$surv.dd_beta*mean(FEAR.IPM.data.intra$neighbors_area)+#Set to average density of neighbors
          fear.params$surv.tminusPrecip_beta*mean(FEAR.IPM.data.intra$tminusPrecip))#set to average precipitation across measure period
  return(u/(1+u))
}

# 2. growth function 
g.yx.fear=function(xp,x,fear.params) {     	
  dnorm(xp,mean=fear.params$growth.int+
          fear.params$growth.slope*x+
          fear.params$growth.dd_beta*mean(FEAR.IPM.data.intra$neighbors_area)+ 
          fear.params$growth.precip_beta*mean(FEAR.IPM.data.intra$Precip),
         sd=fear.params$growth.sd)
}


#3. recruitment function.  This is Dave. A's development
p_bz <- function(basalArea_genet, fear.params)
{
  linear.p <- fear.params$flwr.int + fear.params$flwr.slope * basalArea_genet      # linear predictor
  p <- 1/(1+exp(-linear.p))                                # logistic transformation to probability
  return(p)
}

## Recruitment function (N.B - from birth in spring to first summer), logistic regression
pr_z <- function(fear.params) {
  linear.p <- fear.params$recruit.size.mean                             # linear predictor
  p <- 1/(1+exp(-linear.p))                                 # logistic transformation to probability
  return(p)
}

## Seed production function

b_z <- function(basalArea_genet, fear.params)
{
  N <- exp(fear.params$seed.int + fear.params$seed.slope * basalArea_genet)    # seed production of a size z plant
  return(N)
}

## Recruit size pdf

c_0z1 <- function(basalArea_genet, fear.params)
{
  mu <- fear.params$recruit.size.mean
  sig <- fear.params$recruit.size.sd
  p.deRecr <- dnorm(basalArea_genet, mean = mu, sd = sig)              # pdf of a size z1 recruit
  return(p.deRecr)
}


f.yx.fear=function (size_tplus1, basalArea_genet, fear.params) {
  
  return(p_bz(basalArea_genet, fear.params) * b_z(basalArea_genet, fear.params) * fear.params$establishment.prob * c_0z1(size_tplus1, fear.params))
  
}

##########Plot fecundity based on size using f.yx.fear
fec.out.FEAR<-as.data.frame(f.yx.fear(FEAR.IPM.data.intra$size_tplus1, FEAR.IPM.data.intra$basalArea_genet, fear.params))
fec.out.FEAR$basalArea_genet<-FEAR.IPM.data.intra$basalArea_genet
colnames(fec.out.FEAR)[1] <- "fec.out.FEAR"
fec.lm<-glm(fec.out.FEAR~basalArea_genet, data=fec.out.FEAR)

Fec.fear.size<-ggpredict(fec.lm, terms="basalArea_genet[all]")

FEAR.fecundity<-ggplot(Fec.fear.size, aes(x, predicted))+
  geom_line(color="turquoise3")+
  geom_ribbon(aes(ymin=conf.low, ymax=conf.high, fill="band") ,alpha=0.3)+
  scale_color_manual("",values="turquoise3")+
  scale_fill_manual("", values="turquoise3")+
  xlab(bquote("log(Size)"[t]))+
  ylab("Fecundity")+
  #ylim(-13,0)+
  theme_bw()+
  theme(legend.position="none",
        axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank())
FEAR.fecundity

#Run IPM to generate kernels
###################

min.size=0.9*(min(FEAR.IPM.data.intra$basalArea_genet)); max.size=1.1*(max(FEAR.IPM.data.intra$size_tplus1, na.rm=T)) # integration limits of -.9 to -1.1 of smallest and largest ever observed
n=200 # number of cells in the discretized kernel
b=min.size+c(0:n)*(max.size-min.size)/n #boundary points (the edges of the cells defining the kernel
y=0.5*(b[1:n]+b[2:(n+1)]) # mesh points (midpoints of the cells)
h=y[2]-y[1] # width of the cells

G=h*outer(y,y,g.yx.fear, fear.params) # growth kernel
S=s.x.fear(y,fear.params) # survival 
P=G
for(i in 1:n) P[,i]=G[,i]*S[i]  # growth/survival kernel
Fert=h*outer(y,y,f.yx.fear,fear.params)   # reproduction kernel
image(y,y,t(Fert),main='F. arizonica Fecundity kernel',col=hcl.colors(100,"Mako", rev=T))  # plot fertility kernel
contour(y,y,t(Fert), add = TRUE, drawlabels = TRUE)
K=P  #P kernel
image(y,y,t(K),main='F. arizonica Growth and Survival kernel',col=hcl.colors(100,"Mako", rev=T)) #Plot Growth and Survival kernel
contour(y,y,t(K), add = TRUE, drawlabels = TRUE)
image(y,y,t(K),main='F. arizonica Full kernel ^1; Lambda = 1.3', col=hcl.colors(100,"Mako", rev=T)) #Full kernel, not used because can't see F

# 1. get lamda,v,w  
lam=(eigen(K)$values[1]) 
w.eigen=(eigen(K)$vectors[,1])
stable.dist=w.eigen/sum(w.eigen) 
v.eigen=(eigen(t(K))$vectors[,1])
repro.val=v.eigen/v.eigen[1] 
lam#0.86

#Sensitivity and elasticity
v.dot.w=sum(stable.dist*repro.val)*h
sens=outer(repro.val,stable.dist)/v.dot.w
elas=matrix(as.vector(sens)*as.vector(K)/tam,nrow=n)
e<- apply(elas, 2, as.numeric)
s<- apply(sens, 2, as.numeric)

# 3. plot sensitivity and elasticity results
par(mfrow=c(1,2)) 
#image(y,y,t(K)^.3, xlab="Logsize (t)",ylab="Logsize (t+1)", main="F. arizonica full kernel")
#contour(y,y,t(K), add = TRUE, drawlabels = TRUE)
#plot(y,stable.dist,xlab="Logsize",type="l",main="Stable size distribution")
#plot(y,repro.val,xlab="Logize",type="l",main="Reproductive values") 
image(y,y,t(e),xlab=bquote("ln(size)"[t]),ylab=bquote("ln(size)"[t+1]),main="Elasticity",col=hcl.colors(100,"Mako", rev=T))  
contour(y,y,t(e), add = TRUE, drawlabels = TRUE)
image(y,y,t(s),xlab=bquote("ln(size)"[t]),ylab=bquote("ln(size)"[t+1]), main="Sensitivity", col=hcl.colors(100,"Mako", rev=T))
contour(y,y,t(s), add = TRUE, drawlabels = TRUE)

#########Analyze change in Lambda in response to -3 to +3 SD of temp, precip, and VPD

clim_sd<- 1:3 #Create range of SDs to perturb over

lamda.out.1<-list()#precip+
lamda.out.2<-list()#precip-
lamda.out.3<-list()#temp+
lamda.out.4<-list()#temp-
lamda.out.5<-list()#VPD+
lamda.out.6<-list()#VPD-

for (j in clim_sd){
  ########Perturb climate variables individually.  Change sign for negative and positive SD.
  #########This is currently done manually
  s.x.fear=function(x,fear.params) {
    u=exp(fear.params$surv.int + 
            fear.params$surv.slope*x + 
            fear.params$surv.dd_beta*mean(FEAR.IPM.data.intra$neighbors_area)+
            fear.params$surv.precip_beta*mean(HWB.annual.scale$Precip)+
            #fear.params$surv.precip_beta*(mean(HWB.annual.scale$Precip)-(j*sd(HWB.annual.scale$Precip))+
            fear.params$surv.temp_beta*mean(HWB.annual.scale$MeanT)+
          #fear.params$surv.temp_beta*(mean(HWB.annual.scale$MeanT)-(j*sd(HWB.annual.scale$MeanT))+                                
            #fear.params$surv.VPD_beta*mean(HWB.annual.scale$MeanVPD)+
          fear.params$surv.VPD_beta*(mean(HWB.annual.scale$MeanVPD)-(j*sd(HWB.annual.scale$MeanVPD))))
    #fear.params$surv.VPD_beta*mean(FEAR.IPM.data$MeanVPD))
    return(u/(1+u))
  }
  
  # 2. growth function #covariates are GrowthPrecip, MayVPD, logsize
  g.yx.fear=function(xp,x,fear.params) {     	
    dnorm(xp,mean=fear.params$growth.int+
            fear.params$growth.slope*x+
            fear.params$growth.dd_beta*mean(FEAR.IPM.data.intra$neighbors_area)+ 
           fear.params$growth.precip_beta*mean(HWB.annual.scale$Precip)+
            #fear.params$growth.precip_beta*(mean(HWB.annual.scale$Precip)-(j*sd(HWB.annual.scale$Precip))+
            fear.params$growth.temp_beta*mean(HWB.annual.scale$MeanT)+ 
         #fear.params$growth.temp_beta*(mean(HWB.annual.scale$MeanT)-(j*sd(HWB.annual.scale$MeanT)))+
            #fear.params$growth.VPD_beta*mean(HWB.annual.scale$MeanVPD),
          fear.params$growth.VPD_beta*(mean(HWB.annual.scale$MeanVPD)-(j*sd(HWB.annual.scale$MeanVPD))),
          sd=fear.params$growth.sd)
  }
  
  
  #3. recruitment function
  p_bz <- function(basalArea_genet, fear.params)
  {
    linear.p <- fear.params$flwr.int + fear.params$flwr.slope * basalArea_genet      # linear predictor
    p <- 1/(1+exp(-linear.p))                                # logistic transformation to probability
    return(p)
  }
  
  ## Recruitment function (N.B - from birth in spring to first summer), logistic regression
  pr_z <- function(fear.params) {
    linear.p <- fear.params$recruit.size.mean                             # linear predictor
    p <- 1/(1+exp(-linear.p))                                 # logistic transformation to probability
    return(p)
  }
  
  ## Seed production function
  
  b_z <- function(basalArea_genet, fear.params)
  {
    N <- exp(fear.params$seed.int + fear.params$seed.slope * basalArea_genet)    # seed production of a size z plant
    return(N)
  }
  
  ## Recruit size pdf
  
  c_0z1 <- function(basalArea_genet, fear.params)
  {
    mu <- fear.params$recruit.size.mean
    sig <- fear.params$recruit.size.sd
    p.deRecr <- dnorm(basalArea_genet, mean = mu, sd = sig)              # pdf of a size z1 recruit
    return(p.deRecr)
  }
  
  
  f.yx.fear=function (size_tplus1, basalArea_genet, fear.params) {
    
    return(p_bz(basalArea_genet, fear.params) * b_z(basalArea_genet, fear.params) * fear.params$establishment.prob * c_0z1(size_tplus1, fear.params))
    
  }
  
  ###################
  
  min.size=0.9*(min(FEAR.IPM.data.intra$basalArea_genet)); max.size=1.1*(max(FEAR.IPM.data.intra$size_tplus1, na.rm=T)) # integration limits
  n=100 # number of cells in the discretized kernel
  b=min.size+c(0:n)*(max.size-min.size)/n #boundary points (the edges of the cells defining the kernel
  y=0.5*(b[1:n]+b[2:(n+1)]) # mesh points (midpoints of the cells)
  h=y[2]-y[1] # width of the cells
  
  G=h*outer(y,y,g.yx.fear, fear.params) # growth kernel
  S=s.x.fear(y,fear.params) # survival 
  P=G
  for(i in 1:n) P[,i]=G[,i]*S[i]  # growth/survival kernel
  Fert=h*outer(y,y,f.yx.fear,fear.params)   # reproduction kernel
 # image(y,y,t(Fert),main='F. arizonica F kernel') # plot it
  
  K=P  #P kernel
 # image(y,y,t(K),main='F. arizonica P kernel')
  
 # image(y,y,t(K),main='F. arizonica Full kernel ^1; Lambda = 1.3') 
  
  # 1. get lamda,v,w  
  lam=(eigen(K)$values[1]) 
  w.eigen=(eigen(K)$vectors[,1])
  stable.dist=w.eigen/sum(w.eigen) 
  v.eigen=(eigen(t(K))$vectors[,1])
  repro.val=v.eigen/v.eigen[1] 
  
  lamda.out.6[j]<-lam# Have to change this label for each value of lambda.  
  
}





lambda.out<-c(lamda.out.1, lamda.out.2,lamda.out.3,lamda.out.4,lamda.out.5,lamda.out.6) #Bind lambdas for plotting
lambda.out<-unlist(as.numeric(lambda.out))
lambda.out.FEAR<-as.data.frame(lambda.out)
lambda.out.FEAR$StDev<-c(1,2,3,-1,-2,-3,1,2,3,-1,-2,-3,1,2,3,-1,-2,-3)#least elegant way to create new column of sd values
lambda.out.FEAR$Clim<-c("Precip", "Precip", "Precip", "Precip", "Precip", "Precip", "Temp", "Temp", "Temp", "Temp", "Temp", "Temp", "VPD", "VPD", "VPD", "VPD", "VPD", "VPD")

lam.sd.plot<-ggplot(data=lambda.out.FEAR, aes(x=StDev, y=lambda.out, group=Clim, color=Clim, shape=Clim)) + #Plot
  geom_line(size=1)+
  scale_color_manual(values=c("skyblue3", "darkolivegreen4", "darkgoldenrod"))+
  geom_point(size=3)+
  geom_hline(yintercept=c(0.86,1), color=c("grey48", "red"), linetype=c("dashed","dashed")) + 
  ylim(0.4,1.1)+
  labs(y="Lambda", x="Standard Deviations from Mean")+
  ggtitle("F. arizonica")+
  theme_classic()+
  theme(legend.title = element_blank())
lam.sd.plot



###################
###############MUMO code below
################## Not annotated because identical to FEAR analysis above.
MUMO.data.intra<-data_out.all.intra[ which(data_out.all.intra$species=='Muhlenbergia montana'),]      
MUMO.data$geometry <- NULL

MUMO.IPM.data.intra<-merge(MUMO.data.intra, HWB.annual.scale, by=c("Site", "Plot", "z_Year"))
MUMO.IPM.data.intra$Plot<-as.factor(MUMO.IPM.data.intra$Plot)
MUMO.IPM.data.intra$z_Year<-as.factor(MUMO.IPM.data.intra$z_Year)
MUMO.IPM.data.intra$basalArea_genet<-log(MUMO.IPM.data.intra$basalArea_genet)

count.out.mumo<-count(MUMO.IPM.data.intra,Site,Plot,trackID)#Gives count of unique individuals.  count is row total
count.plot.mumo<-count(MUMO.IPM.data.intra, Site, Plot)#Give count of plots, n=row total
total.plot.count<-merge(count.plot, count.plot.mumo, by=c("Site", "Plot"), all=T)

#########MUMO Parameter List
MUMO.params=data.frame (
  surv.int=NA,
  surv.slope=NA,
  surv.dd_beta=NA,
  surv.precip_beta=NA,
  surv.temp_beta=NA,
  surv.VPD_beta=NA,
  surv.tminusPrecip_beta=NA,
  surv.tminusMeanT_beta=NA,
  surv.tminusMeanVPD_beta=NA,
  #####
  growth.int=NA,
  growth.slope=NA,
  growth.sd=NA,
  growth.dd_beta=NA,
  growth.precip_beta=NA,
  growth.temp_beta=NA,
  growth.VPD_beta=NA,
  growth.tminusPrecip_beta=NA,
  growth.tminusMeanT_beta=NA,
  growth.tminusMeanVPD_beta=NA,
  
  #####
  flwr.int=NA,
  flwr.slope=NA,
  #####
  seed.int=NA,
  seed.slope=NA,
  #####
  recruit.size.mean=NA,
  recruit.size.sd=NA,
  #####
  establishment.prob=NA
)



######MUMO Survival Model
#MUMO.IPM.data.intra$basalArea_genet <- scale(MUMO.IPM.data.intra$basalArea_genet)

surv.MUMO.intra<- glmer(survives_tplus1~basalArea_genet  + 
                          neighbors_area + Precip+ MeanT +MeanVPD+
                          tminusPrecip+tminusMeanT+tminusMeanVPD+
                          (1|Plot)+(1|z_Year),
                        family=binomial,
                        data=MUMO.IPM.data.intra, control = glmerControl(
                          optimizer ='optimx', optCtrl=list(method='nlminb')))

vif(surv.MUMO.intra)

summary(surv.MUMO.intra)#Mean VPD, tminusPrecip
rsquared(surv.MUMO.intra)#extract marginal and conditional R2

MUMO.params$surv.int=fixef(surv.MUMO.intra)[1]
MUMO.params$surv.slope=fixef(surv.MUMO.intra)[2]
MUMO.params$surv.dd_beta=fixef(surv.MUMO.intra)[3]
MUMO.params$surv.precip_beta=fixef(surv.MUMO.intra)[4]
MUMO.params$surv.temp_beta=fixef(surv.MUMO.intra)[5]
MUMO.params$surv.VPD_beta=fixef(surv.MUMO.intra)[6]
MUMO.params$surv.tminusPrecip_beta=fixef(surv.MUMO.intra)[7]
MUMO.params$surv.tminusMeanT_beta=fixef(surv.MUMO.intra)[8]
MUMO.params$surv.tminusMeanVPD_beta=fixef(surv.MUMO.intra)[9]

#set plotting theme
set_theme(base= theme_bw())

#First plot against size
surv.MUMO.predict<-ggpredict(surv.MUMO.intra, terms="basalArea_genet[all]")

surv.MUMO.size<-ggplot(surv.MUMO.predict, aes(x, predicted))+
  geom_line(color="coral1")+
  geom_ribbon(aes(ymin=conf.low, ymax=conf.high, fill="band") ,alpha=0.3)+
  scale_color_manual("",values="coral1")+
  scale_fill_manual("", values="coral1")+
  xlab(bquote("log(Size)"[t]))+
  ylab("Probability of Survival")+
  ylim(0,1)+
  theme_bw()+
  theme(legend.position="none",
        axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank())
surv.MUMO.size

#then against climate variables

surv.MUMO.predict<-ggpredict(surv.MUMO.intra, terms="tminusPrecip[all]")

surv.MUMO.precip<-ggplot(surv.MUMO.predict, aes(x, predicted))+
  geom_line(color="coral1")+
  geom_ribbon(aes(ymin=conf.low, ymax=conf.high, fill="band") ,alpha=0.3)+
  scale_color_manual("",values="coral1")+
  scale_fill_manual("", values="coral1")+
  xlab(bquote("Scaled Precipitation"[t-1]))+
  ylab("Probability of Survival")+
  ylim(0,1)+
  theme_bw()+
  theme(legend.position = "none",
         axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank())
surv.MUMO.precip

surv.MUMO.predict<-ggpredict(surv.MUMO.intra, terms="MeanVPD[all]")

surv.MUMO.VPD<-ggplot(surv.MUMO.predict, aes(x, predicted))+
  geom_line(color="coral1")+
  geom_ribbon(aes(ymin=conf.low, ymax=conf.high, fill="band") ,alpha=0.3)+
  scale_color_manual("",values="coral1")+
  scale_fill_manual("", values="coral1")+
  xlab("Scaled Vapor Pressure Deficit")+
  ylab("Probability of Survival")+
  ylim(0,1)+
  theme_bw()+
  theme(legend.position="none",
        axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank())
surv.MUMO.VPD

######MUMO Growth Model
MUMO.IPM.data.intra<-MUMO.IPM.data.intra[-which(is.na(MUMO.IPM.data.intra$size_tplus1)),]
MUMO.IPM.data.intra$size_tplus1<-log(MUMO.IPM.data.intra$size_tplus1)

growth.MUMO.intra<- lmer(size_tplus1~basalArea_genet  +
                           neighbors_area +
                           Precip+
                           MeanT+
                           MeanVPD+
                           tminusPrecip+tminusMeanT+tminusMeanVPD+
                           (1|Plot)+(1|z_Year),
                         data=MUMO.IPM.data.intra)

vif(growth.MUMO.intra)#MaxT and Max VPD above 8, drop MaxT
summary(growth.MUMO.intra)#Precip, tminusMeanT
rsquared(growth.MUMO.intra)#extract marginal and conditional R2


MUMO.params$growth.int=fixef(growth.MUMO.intra)[1]
MUMO.params$growth.slope=fixef(growth.MUMO.intra)[2]
MUMO.params$growth.dd_beta=fixef(growth.MUMO.intra)[3]
MUMO.params$growth.precip_beta=fixef(growth.MUMO.intra)[4]
MUMO.params$growth.temp_beta=fixef(growth.MUMO.intra)[5]
MUMO.params$growth.VPD_beta=fixef(growth.MUMO.intra)[6]
MUMO.params$growth.tminusPrecip_beta=fixef(growth.MUMO.intra)[7]
MUMO.params$growth.tminusMeanT_beta=fixef(growth.MUMO.intra)[8]
MUMO.params$growth.tminusMeanVPD_beta=fixef(growth.MUMO.intra)[9]
MUMO.params$growth.sd=sd(resid(growth.MUMO.intra))

###First against size

growth.MUMO.size<-ggpredict(growth.MUMO.intra, terms="basalArea_genet[all]")

growth.size<-ggplot(growth.MUMO.size, aes(x, predicted))+
  geom_line(color="coral1")+
  geom_ribbon(aes(ymin=conf.low, ymax=conf.high, fill="band") ,alpha=0.3)+
  scale_color_manual("",values="coral1")+
  scale_fill_manual("", values="coral1")+
  xlab(bquote("log(size)"[t]))+
  ylab(bquote("log(Size)"[t+1]))+
  ylim(-19,0)+
  theme_bw()+
  theme(legend.position="none",
        axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank())
growth.size

growth.MUMO.precip<-ggpredict(growth.MUMO.intra, terms="Precip[all]")
growth.precip<-ggplot(growth.MUMO.precip, aes(x, predicted))+
  geom_line(color="coral1")+
  geom_ribbon(aes(ymin=conf.low, ymax=conf.high, fill="band") ,alpha=0.3)+
  scale_color_manual("",values="coral1")+
  scale_fill_manual("", values="coral1")+
  xlab("Scaled Precipitation")+
  ylab(bquote("log(Size)"[t+1]))+
  ylim(-8,0)+
  theme_bw()+
  theme(legend.position="none",
        axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank())
growth.precip

growth.MUMO.MeanT<-ggpredict(growth.MUMO.intra, terms="tminusMeanT[all]")
growth.tminusMeanT<-ggplot(growth.MUMO.precip, aes(x, predicted))+
  geom_line(color="coral1")+
  geom_ribbon(aes(ymin=conf.low, ymax=conf.high, fill="band") ,alpha=0.3)+
  scale_color_manual("",values="coral1")+
  scale_fill_manual("", values="coral1")+
  xlab(bquote("Scaled Temperature"[t-1]))+
  ylab(bquote("log(Size)"[t+1]))+
  ylim(-8,0)+
  theme_bw()+
  theme(legend.position="none",
        axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank())
growth.tminusMeanT

############MUMO Tables
tab.surv.mumo<-tbl_regression(surv.MUMO.intra)
tab.surv.mumo<-bold_p(tab.surv.mumo, t=0.05, q=F)
tab.grow.mumo<-tbl_regression(growth.MUMO.intra)
tab.grow.mumo<-bold_p(tab.grow.mumo, t=0.05, q=F)

tab.mumo<-tbl_merge(tbls=list(tab.surv.mumo, tab.grow.mumo), 
                    tab_spanner=c("**Survival**", "**Growth**"))



tab.mumo

########Publication table
table.1<-tbl_stack(tbls=list(tab.fear, tab.mumo), group_header=c("F. arizonica", "M. montana"))
table.1


### MUMO Fecundity
fec.MUMO = read.csv(file.choose(), header=T)
fec.MUMO$basalArea_genet = log((fec.MUMO$basal.area.cm2/10000))
#fec.MUMO$basalArea_genet<-scale(fec.MUMO$basalArea_genet)
fec.MUMO = fec.MUMO[!is.na(fec.MUMO$basalArea_genet),]
glume.MUMO = read.csv(file.choose(), header=T)
glume.MUMO$basalArea_genet = log((glume.MUMO$BA.cm2/10000))
#glume.MUMO$basalArea_genet<-scale(glume.MUMO$basalArea_genet)
glume.MUMO$basalArea_genet<-as.numeric(glume.MUMO$basalArea_genet)
fec.MUMO$flwr = ifelse(fec.MUMO$numb.stalks>0, 1, 0)

flower.MUMO.linear = glm(flwr ~ basalArea_genet, data=fec.MUMO, family = "binomial")
summary(flower.MUMO.linear)
MUMO.params$flwr.int=flower.MUMO.linear$coefficients[1]
MUMO.params$flwr.slope=flower.MUMO.linear$coefficients[2]

#Function to convert flowering to probability, for use in reproduction kernel
p_bz <- function(basalArea_genet, params)
{
  linear.p <- params$flwr.int + params$flwr.slope * basalArea_genet      # linear predictor
  p <- 1/(1+exp(-linear.p))                                # logistic transformation to probability
  return(p)
}

seeds.MUMO = glm(sum.glumes ~ basalArea_genet, data=glume.MUMO, family = "poisson")
summary(seeds.MUMO)

MUMO.params$seed.int=seeds.MUMO$coefficients[1]
MUMO.params$seed.slope=seeds.MUMO$coefficients[2]

newdat=data.frame(MUMO.IPM.data.intra$basalArea_genet)
names(newdat)[1]="basalArea_genet"

MUMO.IPM.data.intra$Seeds<-round(predict(seeds.MUMO, newdata=newdat, type="response"))

MUMO.params$recruit.size.mean=mean(MUMO.IPM.data.intra$basalArea_genet[MUMO.IPM.data.intra$recruit==1], na.rm =TRUE)
MUMO.params$recruit.size.sd=sd(MUMO.IPM.data.intra$basalArea_genet[MUMO.IPM.data.intra$recruit==1], na.rm =TRUE)

MUMO.params$establishment.prob=sum(MUMO.IPM.data.intra$recruit, na.rm = TRUE)/sum(MUMO.IPM.data.intra$Seeds,na.rm=TRUE)

#########Calculating IPMS
s.x.MUMO=function(x,MUMO.params) {
  u=exp(MUMO.params$surv.int + 
          MUMO.params$surv.slope*x + 
          MUMO.params$surv.dd_beta*mean(MUMO.IPM.data.intra$neighbors_area)+
          MUMO.params$surv.tminusPrecip_beta*mean(MUMO.IPM.data.intra$tminusPrecip)+
          MUMO.params$surv.VPD_beta*mean(MUMO.IPM.data.intra$MeanVPD))
  return(u/(1+u))
}

# 2. growth function 
g.yx.MUMO=function(xp,x,MUMO.params) {     	
  dnorm(xp,mean=MUMO.params$growth.int+
          MUMO.params$growth.slope*x+
          MUMO.params$growth.dd_beta*mean(MUMO.IPM.data.intra$neighbors_area)+ 
          MUMO.params$growth.precip_beta*mean(MUMO.IPM.data.intra$Precip)+
          MUMO.params$growth.precip_beta*mean(MUMO.IPM.data.intra$tminusMeanT),
        sd=MUMO.params$growth.sd)
}


#3. recruitment function
p_bz <- function(basalArea_genet, MUMO.params)
{
  linear.p <- MUMO.params$flwr.int + MUMO.params$flwr.slope * basalArea_genet      # linear predictor
  p <- 1/(1+exp(-linear.p))                                # logistic transformation to probability
  return(p)
}

## Recruitment function (N.B - from birth in spring to first summer), logistic regression
pr_z <- function(MUMO.params) {
  linear.p <- MUMO.params$recruit.size.mean                             # linear predictor
  p <- 1/(1+exp(-linear.p))                                 # logistic transformation to probability
  return(p)
}

## Seed production function

b_z <- function(basalArea_genet, MUMO.params)
{
  N <- exp(MUMO.params$seed.int + MUMO.params$seed.slope * basalArea_genet)    # seed production of a size z plant
  return(N)
}

## Recruit size pdf

c_0z1 <- function(basalArea_genet, MUMO.params)
{
  mu <- MUMO.params$recruit.size.mean
  sig <- MUMO.params$recruit.size.sd
  p.deRecr <- dnorm(basalArea_genet, mean = mu, sd = sig)              # pdf of a size z1 recruit
  return(p.deRecr)
}


f.yx.MUMO=function (size_tplus1, basalArea_genet, MUMO.params) {
  
  return(p_bz(basalArea_genet, MUMO.params) * b_z(basalArea_genet, MUMO.params) * MUMO.params$establishment.prob * c_0z1(size_tplus1, MUMO.params))
  
}

###Plot fecundity probability 

fec.out.MUMO<-as.data.frame(f.yx.MUMO(MUMO.IPM.data.intra$size_tplus1, MUMO.IPM.data.intra$basalArea_genet, MUMO.params))
fec.out.MUMO$basalArea_genet<-MUMO.IPM.data.intra$basalArea_genet
colnames(fec.out.MUMO)[1] <- "fec.out.MUMO"
fec.lm.MUMO<-glm(fec.out.MUMO~basalArea_genet, data=fec.out.MUMO)

fec.MUMO.size<-ggpredict(fec.lm.MUMO, terms="basalArea_genet[all]")

MUMO.fecundity<-ggplot(fec.MUMO.size, aes(x, predicted))+
  geom_line(color="coral1")+
  geom_ribbon(aes(ymin=conf.low, ymax=conf.high, fill="band") ,alpha=0.3)+
  scale_color_manual("",values="coral1")+
  scale_fill_manual("", values="coral1")+
  xlab(bquote("log(Size)"[t]))+
  ylab("Fecundity")+
  #ylim(-13,0)+
  theme_bw()+
  theme(legend.position="none",
        axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank())
MUMO.fecundity

###################




###################

min.size=0.9*(min(MUMO.IPM.data.intra$basalArea_genet)); max.size=1.1*(max(MUMO.IPM.data.intra$size_tplus1, na.rm=T)) # integration limits
n=200 # number of cells in the discretized kernel
b=min.size+c(0:n)*(max.size-min.size)/n #boundary points (the edges of the cells defining the kernel
y=0.5*(b[1:n]+b[2:(n+1)]) # mesh points (midpoints of the cells)
h=y[2]-y[1] # width of the cells

G=h*outer(y,y,g.yx.MUMO, MUMO.params) # growth kernel
S=s.x.MUMO(y,MUMO.params) # survival 
P=G
for(i in 1:n) P[,i]=G[,i]*S[i]  # growth/survival kernel
Fert=h*outer(y,y,f.yx.MUMO,MUMO.params)   # reproduction kernel
image(y,y,t(Fert),main='M. montana Fecundity kernel',col=hcl.colors(100,"Red-Purple", rev=T))  # plot it
contour(y,y,t(Fert), add = TRUE, drawlabels = TRUE)
K=P  #P kernel
image(y,y,t(K),main='M. montana Growth and Survival kernel',col=hcl.colors(100,"Red-Purple", rev=T)) 
contour(y,y,t(K), add = TRUE, drawlabels = TRUE)
image(y,y,t(K),main='M. montana Full kernel ^1; Lambda = 1.3', col=hcl.colors(100,"Red-Purple", rev=T)) 

# 1. get lamda,v,w  
lam=(eigen(K)$values[1]) 
w.eigen=(eigen(K)$vectors[,1])
stable.dist=w.eigen/sum(w.eigen) 
v.eigen=(eigen(t(K))$vectors[,1])
repro.val=v.eigen/v.eigen[1] 
lam#0.76

#Sensitivity and elasticity
v.dot.w=sum(stable.dist*repro.val)*h
sens=outer(repro.val,stable.dist)/v.dot.w
elas=matrix(as.vector(sens)*as.vector(K)/lam,nrow=n)
e<- apply(elas, 2, as.numeric)
s<- apply(sens, 2, as.numeric)

# 3. plot results
par(mfrow=c(1,2)) 
#image(y,y,t(K)^.3, xlab="Logsize (t)",ylab="Logsize (t+1)", main="M. montana full kernel")
#contour(y,y,t(K), add = TRUE, drawlabels = TRUE)
#plot(y,stable.dist,xlab="Logsize",type="l",main="Stable size distribution")
#plot(y,repro.val,xlab="Logize",type="l",main="Reproductive values") 
image(y,y,t(e),xlab=bquote("ln(size)"[t]),ylab=bquote("ln(size)"[t+1]),main="Elasticity",col=hcl.colors(100,"Red-Purple", rev=T))  
contour(y,y,t(e), add = TRUE, drawlabels = TRUE)
image(y,y,t(s),xlab=bquote("ln(size)"[t]),ylab=bquote("ln(size)"[t+1]), main="Sensitivity", col=hcl.colors(100,"Red-Purple", rev=T))
contour(y,y,t(s), add = TRUE, drawlabels = TRUE)

#########Analyze change in Lambda in response to -3 to +3 SD of temp, precip, and VPD

clim_sd<- 1:3

lamda.out.MUMO1<-list()#precip+
lamda.out.MUMO2<-list()#precip-
lamda.out.MUMO3<-list()#temp+
lamda.out.MUMO4<-list()#temp-
lamda.out.MUMO5<-list()#VPD+
lamda.out.MUMO6<-list()#VPD-

for (j in clim_sd){
  ########Perturb climate variables.  Change sign for negative and postive SD
  s.x.MUMO=function(x,MUMO.params) {
    u=exp(MUMO.params$surv.int + 
            MUMO.params$surv.slope*x + 
            MUMO.params$surv.dd_beta*mean(MUMO.IPM.data.intra$neighbors_area)+
            MUMO.params$surv.precip_beta*mean(HWB.annual.scale$Precip)+
            #MUMO.params$surv.precip_beta*(mean(HWB.annual.scale$Precip)-(j*sd(HWB.annual.scale$Precip)))+
            MUMO.params$surv.temp_beta*mean(HWB.annual.scale$MeanT)+
            #MUMO.params$surv.temp_beta*(mean(HWB.annual.scale$MeanT)-(j*sd(HWB.annual.scale$MeanT))+                                
            MUMO.params$surv.VPD_beta*(mean(HWB.annual.scale$MeanVPD)-(j*sd(HWB.annual.scale$MeanVPD))))
            #MUMO.params$surv.VPD_beta*mean(HWB.annual.scale$MeanVPD)))
    return(u/(1+u))
  }
  
  # 2. growth function #covariates are GrowthPrecip, MayVPD, logsize
  g.yx.MUMO=function(xp,x,MUMO.params) {     	
    dnorm(xp,mean=MUMO.params$growth.int+
            MUMO.params$growth.slope*x+
            MUMO.params$growth.dd_beta*mean(MUMO.IPM.data.intra$neighbors_area)+ 
            MUMO.params$growth.precip_beta*mean(HWB.annual.scale$Precip)+
            #MUMO.params$growth.precip_beta*(mean(HWB.annual.scale$Precip)-(j*sd(HWB.annual.scale$Precip))+
            MUMO.params$growth.temp_beta*mean(HWB.annual.scale$MeanT)+ 
            #MUMO.params$growth.temp_beta*(mean(HWB.annual.scale$MeanT)-(j*sd(HWB.annual.scale$MeanT))+
            #MUMO.params$growth.VPD_beta*mean(HWB.annual.scale$MeanVPD)),
            MUMO.params$growth.VPD_beta*(mean(HWB.annual.scale$MeanVPD)-(j*sd(HWB.annual.scale$MeanVPD))),
          sd=MUMO.params$growth.sd)
  }
  
  
  #3. recruitment function
  p_bz <- function(basalArea_genet, MUMO.params)
  {
    linear.p <- MUMO.params$flwr.int + MUMO.params$flwr.slope * basalArea_genet      # linear predictor
    p <- 1/(1+exp(-linear.p))                                # logistic transformation to probability
    return(p)
  }
  
  ## Recruitment function (N.B - from birth in spring to first summer), logistic regression
  pr_z <- function(MUMO.params) {
    linear.p <- MUMO.params$recruit.size.mean                             # linear predictor
    p <- 1/(1+exp(-linear.p))                                 # logistic transformation to probability
    return(p)
  }
  
  ## Seed production function
  
  b_z <- function(basalArea_genet, MUMO.params)
  {
    N <- exp(MUMO.params$seed.int + MUMO.params$seed.slope * basalArea_genet)    # seed production of a size z plant
    return(N)
  }
  
  ## Recruit size pdf
  
  c_0z1 <- function(basalArea_genet, MUMO.params)
  {
    mu <- MUMO.params$recruit.size.mean
    sig <- MUMO.params$recruit.size.sd
    p.deRecr <- dnorm(basalArea_genet, mean = mu, sd = sig)              # pdf of a size z1 recruit
    return(p.deRecr)
  }
  
  
  f.yx.MUMO=function (size_tplus1, basalArea_genet, MUMO.params) {
    
    return(p_bz(basalArea_genet, MUMO.params) * b_z(basalArea_genet, MUMO.params) * MUMO.params$establishment.prob * c_0z1(size_tplus1, MUMO.params))
    
  }
  
  ###################
  
  min.size=0.9*(min(MUMO.IPM.data.intra$basalArea_genet)); max.size=1.1*(max(MUMO.IPM.data.intra$size_tplus1, na.rm=T)) # integration limits
  n=100 # number of cells in the discretized kernel
  b=min.size+c(0:n)*(max.size-min.size)/n #boundary points (the edges of the cells defining the kernel
  y=0.5*(b[1:n]+b[2:(n+1)]) # mesh points (midpoints of the cells)
  h=y[2]-y[1] # width of the cells
  
  G=h*outer(y,y,g.yx.MUMO, MUMO.params) # growth kernel
  S=s.x.MUMO(y,MUMO.params) # survival 
  P=G
  for(i in 1:n) P[,i]=G[,i]*S[i]  # growth/survival kernel
  Fert=h*outer(y,y,f.yx.MUMO,MUMO.params)   # reproduction kernel
  # image(y,y,t(Fert),main='F. arizonica F kernel') # plot it
  
  K=P  #P kernel
  # image(y,y,t(K),main='F. arizonica P kernel')
  
  # image(y,y,t(K),main='F. arizonica Full kernel ^1; Lambda = 1.3') 
  
  # 1. get lamda,v,w  
  lam=(eigen(K)$values[1]) 
  w.eigen=(eigen(K)$vectors[,1])
  stable.dist=w.eigen/sum(w.eigen) 
  v.eigen=(eigen(t(K))$vectors[,1])
  repro.val=v.eigen/v.eigen[1] 
  
  lamda.out.MUMO6[j]<-lam
  
}





lambda.out.MUMO<-c(lamda.out.MUMO1, lamda.out.MUMO2,lamda.out.MUMO3,lamda.out.MUMO4,lamda.out.MUMO5,lamda.out.MUMO6)
lambda.out.MUMO<-unlist(as.numeric(lambda.out.MUMO))
lambda.out.MUMO<-as.data.frame(lambda.out.MUMO)
lambda.out.MUMO$StDev<-c(1,2,3,-1,-2,-3,1,2,3,-1,-2,-3,1,2,3,-1,-2,-3)#least elegant way to create new column of sd values
lambda.out.MUMO$Clim<-c("Precip", "Precip", "Precip", "Precip", "Precip", "Precip", "Temp", "Temp", "Temp", "Temp", "Temp", "Temp", "VPD", "VPD", "VPD", "VPD", "VPD", "VPD")

lam.sd.plot.MUMO<-ggplot(data=lambda.out.MUMO, aes(x=StDev, y=lambda.out.MUMO, group=Clim, color=Clim, shape=Clim)) +
  geom_line(size=1)+
  scale_color_manual(values=c("skyblue3", "darkolivegreen4", "darkgoldenrod"))+
  geom_point(size=3)+
  geom_hline(yintercept=c(0.76,1), color=c("grey48", "red"), linetype=c("dashed","dashed")) +
  ylim(0.4,1.1)+
  labs(y="Lambda", x="Standard Deviations from Mean")+
  ggtitle("M. montana")+
  theme_classic()+
  theme(legend.title = element_blank())
lam.sd.plot.MUMO
