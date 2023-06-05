############Appendix 3. R Code for analysis (To be put up on GitHub)############
################################################################################
################### SPECIES AND FUNCTIONAL INDICES #############################
################################################################################
remotes::install_github("paleolimbot/rbbt")


Taxa <- read.csv("~/Arctic Ph.D/Chapter 1. Functional analysis HMSC/Data/GitHub_NEG/Data for R/Abun_species_vegan.csv", 
                 header=TRUE, sep=",", row = 1)
Taxa <- Taxa[order(rownames(Taxa)),]
Traits <- read.csv("~/Arctic Ph.D/Chapter 1. Functional analysis HMSC/Data/GitHub_NEG/Data for R/Abun_traits_vegan.csv", 
                   header=TRUE, sep=",", row = 1)
Env <- read.csv("~/Arctic Ph.D/Chapter 1. Functional analysis HMSC/Data/GitHub_NEG/Data for R/Env.csv",  header=TRUE, sep=",", row = 1)
Env <- Env[order(rownames(Env)),]
# load libraries 
library(vegan)   #install.packages("devtools")
library(devtools)  #devtools::install_github("gavinsimpson/ggvegan")
library(ggvegan)
library(FD)
library(adiv)
library(RColorBrewer)

#Taxa indices----------------------------------------------------------------
indices <- data.frame(Abundance = rowSums(Taxa))
indices$Sp.Richness <- specnumber(Taxa)
indices$shannon <- diversity(Taxa, index = "shannon")
indices$simpson <- diversity(Taxa, index = "simpson")
indices$evenness <- diversity(Taxa, index = "shannon")/log(specnumber(Taxa)) #Evenness
indices$Habitat <- as.factor(Env$Habitat)

head(indices)

write.csv(Sp.index, "Sp_index.csv") # save species indices

#Vulnerability------------------------------------------------------------------
Vulnerability <- data.frame(Mean = colSums(Uni$V, na.rm = TRUE), 
                            SD = apply(Uni$V, 2, sd, na.rm=TRUE))

Vulnerability
write.csv(Vulnerability, "Vulnerability_abundance.csv") # save sp. vulnerability

#FUNCTIONAL INDICES-------------------------------------------------------------
# Load and check datasets
Traits <- traits 
Species <- Abundance.hmsc
#The funtional indices 
ex3<-dbFD(Traits, Species, corr="sqrt", m="max")  # m = max means to try with 
# the highest number of PCA axis. Low axes mean loss of information. 
ex3 

data<-as.data.frame(ex3$nbsp) # save the data
data$sing.sp<-(ex3$sing.sp)
data$FRic<-(ex3$FRic)
data$FEve<-(ex3$FEve)
data$FDiv<-(ex3$FDiv)
data$FDis<-(ex3$FDis)
data$RaoQ<-(ex3$RaoQ)
data$Station <- rownames(data) # write in the stations 
data$Habitat <- as.factor(Env$Habitat) # Write in the habitat type. 
data
write.csv(data, 'Functional_Indices.csv') # save functional indices
CWM <- ex3$CWM
write.csv(ex3$CWM,'CWM_traits.csv') # save CWM trait values. 

#Functional redundancy----------------------------------------------------------
#data frame with species as columns and stations are rows 
Species
#data frame with traits as columns and species as rows 
Traits

dist.MS <- function (Species, diag = FALSE, upper = FALSE, tol = 1e-07) 
{
  
  df <- data.frame(Species)
  if (!inherits(df, "data.frame")) 
    stop("df is not a data.frame")
  nlig <- nrow(df)
  d <- matrix(0, nlig, nlig)
  d.names <- row.names(df)
  fun1 <- function(x) {
    sum(abs(df[x[1], ] - df[x[2], ]))/sum(apply(cbind.data.frame(df[x[1], ], df[x[2], ]), 1, max))
  }
  df <- as.matrix(df)
  index <- cbind(col(d)[col(d) < row(d)], row(d)[col(d) < row(d)])
  
  d <- unlist(apply(index, 1, fun1))
  
  attr(d, "Size") <- nlig
  attr(d, "Labels") <- d.names
  attr(d, "Diag") <- diag
  attr(d, "Upper") <- upper
  attr(d, "method") <- "Marczewski–Steinhaus"
  attr(d, "call") <- match.call()
  class(d) <- "dist"
  return(d)
}

Dis <- dist.MS(Traits)
Uni <- uniqueness(Species, Dis, abundance = TRUE)

fac <- factor(Env$Habitat)
sapply(Uni$red, function(x) tapply(x, fac, mean))
sapply(Uni$red, function(x) tapply(x, fac, sd))

Uni$red$Station <- Env$Habitat
write.csv(Uni$red, "Functional_Redundancy.csv") # save functional redundancy 

#End----------------------------------------------------------------------------

#############################################################################
########################## STATISTICAL TESTING ##############################
#############################################################################
########################### ADONIS / PERMANOVA #############################
#############################################################################

#Species------------------------------------------------------------------------
## ADONIS: tests for sig differences in communities among the habitat groups 
## First species: NB* bray method used here
test.species.1 <- adonis(Species~Habitat,
                         data = Env, permutations = 999, method="bray")
test.species.1 # Yes, there is a sig difference between groups
# assumption: normally distributed. 
distances.Species<- vegdist(Species, method = "bray") 
test.species.2 <- anova(betadisper(distances.Species, Env$Habitat)) 
test.species.2 # F value is much smaller. P = 39% of the time we will get 
# greater than 1. Conclude: the dispersion/variance of the data is equal 

## Second Traits : NB* Gower method used here. 
test.traits.1 <- adonis(CWM~Habitat,
                        data = Env, permutations = 999, method="gower")
test.traits.1
# assumption: normally distributed. 
distances.CWM<- vegdist(CWM, method = "gower")
test.traits.2 <- anova(betadisper(distances.CWM, Env$Habitat))
# are there difference in the dispersion parameters that we've just created 
# (distance). Data dispersion is nested within the anova
test.traits.2 # P = 36% # Conclude: the dispersion/variance of the data is equal 

## SIMPER ANALYSIS 
# Species #
simper.species <- simper(Species, group = Env$Habitat,permutations =999)
simper.species
(summary(simper.species, ordered=TRUE))

# Traits #
simper.traits <- simper(CWM, group = Env$Habitat, permutations =999)
simper.traits
(summary(simper.traits, ordered=TRUE))

#End----------------------------------------------------------------------------

#############################################################################
#################################### HMSC ###################################
#############################################################################

# Libraries---------------------------------------------------------------------
library(Hmsc);library(parallel);library(corrplot)

######################### Import Data ##########################################
# Read in reduced data set(s) (see methods in main script)
# Environmental----------------------------------------------------------------- 
Env <- read.csv(file="Env.csv", sep = ";")
Env <- Env[,-c(1)]
rownames(Env) <- Env$Station
Env<-Env[order(Env$Station),]
Env$Habitat <- as.factor(Env$Habitat)
Env$logDepth <- log(Env$Depth)
Env

# Taxa abundance / Change for biomass ------------------------------------------------------------ 
Abundance.hmsc <- read.csv(file="Abundance_95_taxa.csv", sep = ";")
rownames(Abundance.hmsc) <- Abundance.hmsc[,1]
Abundance.hmsc <- Abundance.hmsc[,-c(1)]
Abundance.hmsc$Station <- rownames(Abundance.hmsc)
Abundance.hmsc<-Abundance.hmsc[order(Abundance.hmsc$Station),]
Abundance.hmsc <- Abundance.hmsc[,1:46] #This may vary according to dataset.
str(Abundance.hmsc)
Abundance.hmsc

# Species traits / Change for biomass ----------------------------------------------------------------
#traits <- read.csv(file="Abundance_95_traits.csv", sep=";")
traits <- read.csv(file="Abundance_95_traits.csv", sep = ";")
rownames(traits) <- traits$X
traits <- traits[,-c(1)]
str(traits)
traits


# Phylogenetic tree-------------------------------------------------------------
my.tree <- ape::read.tree(file="my_tree.phy")
plot(my.tree)

########################## SETTING UP THE MODEL ################################
Y <- Abundance.hmsc             # Taxa abundance (Y matrix)
logY <- log(Y+1)			# Log Y because of spread of data 
Ypa = 1*(Y > 0)             	# presence/absense model if wanted. 
Ycon <- logY                    	# If you want to do conditional model (part.1)
Ycon[Ycon==0] = NA  		# If you want to do conditional model (part.2)
Ycon                 			# Check data
my.tree                         		# Taxonomic Tree
traits                          		# Traits of taxa
TrData <- (traits)              	# TrData (trait data used in model)
Env                             		# Environmental data
XData <- (Env)                  	# XData (Env used in model)
n <- length(rownames(Env))     # number of stations
ns <- length(rownames(traits))  # number of taxa
rownames(Env)
rownames(logY)

studyDesign <- data.frame(sample=as.factor(XData$Station), 
                          habitat = as.factor(XData$Habitat))
#Random effect levels: 
rL1 <- HmscRandomLevel(units = levels(studyDesign$sample))
rL2 <- HmscRandomLevel(units = levels(studyDesign$habitat))
#Fixed effect levels: 
XFormula <- ~  logDepth + B.Temp + B.PSU + Oxygen + Habitat 
#Traits used in analysis
TrFormula<-~S1+ S2+ S3+ S4+ S5+ BF1+ BF2+ BF3+ BF4+ BF5+SK1+ SK2 +SK3+SK4+SK5+
  R1+R2+ R3+ R4+ LD1+ LD2+ LD3+ MV1 +MV2+ MV3+ MV4 +MO1 + MO2+MO3+MO4 +
  FH1 + FH2 + FH3 + FH4 +  Z1+ Z2+ Z3+ Z4

# NOTE: model is using Gaussian family /a normal distribution
model <- Hmsc(Y=logY, XData = XData, XFormula = XFormula, YScale=TRUE,
              TrData = TrData, TrFormula = TrFormula,
              phyloTree = my.tree,
              studyDesign = studyDesign,
              ranLevels = list(sample = rL1, habitat = rL2), 
              distr = "normal") 


# set MCMC parameters
nChains = 4
samples = 250
verbose = 100

for (thin in c(1,50,100,1000))
{  transient = 50*thin
model = sampleMcmc(model, thin = thin, samples = samples, transient = transient,
                   nChains = nChains, nParallel = nChains,
                   verbose = verbose, alignPost = TRUE, initPar = "fixed effects")
filename=file.path(paste0("Reduced_model_chains_",as.character(nChains),
                          "_samples_",as.character(samples),
                          "_thin_",as.character(thin),".Rdata"))
save(model,file=filename)

}

#End of model setting-----------------------------------------------------------

########################## EXPLORING THE HMSC MODEL #######################
#Model exploration--------------------------------------------------------------
model
model$studyDesign
model$ranLevelsUsed
model$TrFormula
model$XFormula
model$samples
model$thin
model$transient

# Alpha = the posterior (intercept/space)
# Beta = species niches / the slope, the relationship between species 
# and covariates
# Gamma =  the influence of traits on species niches / relationship between 
# species niches and traits
# omega = residual species associations 
# rho = strength of phylogenetic signal

mpost = convertToCodaObject(model)
postBeta = getPostEstimate(model, parName = "Beta") 
plotBeta(model, post = postBeta, param = "Sign", supportLevel = 0.90, 
         cex = c(0.3, 0.8, 0.8),
         spNamesNumbers = c(F,F), plotTree=TRUE, split = 0.4)
postGamma = getPostEstimate(model, parName = "Gamma") 
plotGamma(model, post = postGamma, param = "Sign", supportLevel = 0.90,
          cex = c(0.7, 0.3, 0.8))
par(mfrow=c(1,1))
OmegaCor <- computeAssociations(model)
supportLevel <- 0.95
toPlot <- ((OmegaCor[[1]]$support>supportLevel)
           + (OmegaCor[[1]]$support<(1-supportLevel))>0)*OmegaCor[[1]]$mean
corrplot(toPlot, method = "color",
         col=colorRampPalette(c("blue","white","red"))(200),
         tl.cex=.6, tl.col="black",
         title=paste("random effect level:", model$rLNames[1]), mar=c(0,0,1,0))

corrplot(toPlot,  type='lower', method = "shade", order = 'FPC', diag=FALSE, 
         tl.cex=.6, tl.col="black",
         title=paste("random effect level:", model$rLNames[1]), mar=c(0,0,1,0))

# Species residual associations 

############################## DIAGNOSTICS ##################################
# Alpha-------------------------------------------------------------------------
#plot(mpost$Alpha[[1]])
#hist(effectiveSize(mpost$Alpha[[1]]))
#summary(mpost$Alpha [[1]])[2]$quantiles [1,] #strong spatial signal if bounded
#away from zero
# Beta--------------------------------------------------------------------------
summary(mpost$Beta)
plot(mpost$Beta)
ess.beta = effectiveSize(mpost$Beta)
psrf.beta = gelman.diag(mpost$Beta, multivariate = FALSE, 
                        autoburnin = FALSE)$psrf
hist(ess.beta)
hist(psrf.beta, xlab="Potential scale reduction factor (beta)")
mean(gelman.diag(mpost$Beta, multivariate = FALSE, autoburnin = FALSE)$psrf)
sd(gelman.diag(mpost$Beta, multivariate = FALSE, autoburnin = FALSE)$psrf)
# Gamma-------------------------------------------------------------------------
ess.gamma = effectiveSize(mpost$Gamma)
psrf.gamma = gelman.diag(mpost$Gamma, multivariate = FALSE, 
                         autoburnin = FALSE)$psrf
hist(ess.gamma)
hist(psrf.gamma, xlab="Potential scale reduction factor (gamma)")
mean(gelman.diag(mpost$Gamma, multivariate = FALSE, autoburnin = FALSE)$psrf)
sd(gelman.diag(mpost$Gamma, multivariate = FALSE, autoburnin = FALSE)$psrf)
# Omega-------------------------------------------------------------------------
sppairs = matrix(sample(x = 1:ns^2, size = 100))
tmp = mpost$Omega[[1]]
for (chain in 1:length(tmp)){
  tmp[[chain]] = tmp[[chain]][,sppairs]
}
ess.omega = effectiveSize(tmp)
psrf.omega = gelman.diag(tmp, multivariate=FALSE, autoburnin = FALSE)$psrf
mean(psrf.omega)
sd(psrf.omega)
hist(ess.omega)
hist(psrf.omega)
# Rho---------------------------------------------------------------------------
summary(mpost$Rho) # posterior distribution of rho reveals strong evidence for 
#phylogenetic signal p = 1
effectiveSize(mpost$Rho)
gelman.diag(mpost$Rho, multivariate = FALSE, autoburnin = FALSE)$psrf


######################### VARIANCE PARTITIONING ############################
## To explain how much variation in species niches and occurrences 
head(model$X)
par(mfrow=c(1,1))
VP <- computeVariancePartitioning(model)
plotVariancePartitioning(model, VP = VP)
(VP$R2T$Beta) #The proportions of variance in species niches explained by 
# traits included in the model. A major part of the variation among species 
#niches related to their responses to habitat quality is explained by the traits. 

(VP$R2T$Y) #  the proportion of variance that the included traits explain 
# out of species abundances. The traits explain 11.8% of variance in species 
# abundances

########################## MODEL FIT #####################################

preds <- computePredictedValues(model)
MF <- evaluateModelFit(hM=model, predY=preds)
MF # RMSE and R2 for normal distribution 
# Hist of species-specific explanatory r2 variables.
hist(MF$R2, xlim = c(0,1), main=paste0("Mean = ", round(mean(MF$R2),2)))
hist(MF$RMSE) # Root mean squared errpr
#mean and standard deviation of R2
mean(MF$R2)
sd(MF$R2)

#Note: cross validation method cannot work here due to the sample size of data. 

########################### PLOTS ###########################################
# library
library(ggplot2);library(viridis);library(tidyr)

#ALL VARIANCE PARTITIONING -----------------------------------------------------
#Explained R2 variance partitioning
VP = computeVariancePartitioning(model)
plotVariancePartitioning(model, VP=VP)
v.part <- as.data.frame(VP$vals)
dim(v.part)
v.part <- as.data.frame(t(v.part))
v.part$Taxa <- rownames(v.part)
v.part$R2 <- MF$R2
v.part$RMSE <- MF$RMSE
(v.part)

#Make explained partition.
v.part$R2
Sp.var.part <- data.frame((v.part$R2/rowSums(v.part[1:6]))*v.part[1:6])
rowSums((v.part$R2/rowSums(v.part[1:6]))*v.part[1:6])
rowSums(Sp.var.part)
Sp.var.part$R2 <- v.part$R2
#Looks ok. now order to most explained: 
dim(Sp.var.part)
Sp.var.part <- Sp.var.part[order(Sp.var.part$R2),]
#Sp.var.part <- Sp.var.part[order(Sp.var.part$B.PSU),]
Sp.var.part$Taxa <- 30:1 #30 or 46 (Biomass or abundance data set respectively)
Sp.var.part$Taxa_Names <- rownames(Sp.var.part)
Sp.var.part

#Abundance 
mean(rowSums(v.part[1:4])*100) # 43% expl. by environmental covariates 
mean(rowSums(v.part[5:6])*100) # 20% expl. by random eff. in sample and habitat
mean(v.part$R2)*100 
sd(v.part$R2)*100 # 63% mean R2 - the model explains 63%. 
#Biomass
Sp.var.part_bio <- Sp.var.part
mean(rowSums(Sp.var.part_bio[1:4])*100) # 65.3% expl. by environmental covariates 
mean(rowSums(Sp.var.part_bio[5:6])*100) # 12,2% expl. by random eff.in sample & habitat
mean(Sp.var.part_bio$R2)*100            # 77.5% mean R2 - the model explains 77,5%. 
sd(Sp.var.part_bio$R2)*100              # 22,6% sd R2  

Sp.var.part2<- gather(Sp.var.part, Env, value,B.PSU,B.Temp,logDepth,Habitat, 
                      `Random..sample`,`Random..habitat`, factor_key=TRUE)

col <- c(#"dimgrey", 
  "darkslategray", "darkslategray4", "darkslategray3","darkslategray2", 
  "darkorange",  "darkorange2")

# Stacked : Abundance dataset --------------------------------------------------
Sp.var.part2
Sp.var.part2_plot <- ggplot(Sp.var.part2, aes(fill=forcats::fct_infreq(Env), 
                                              y=value, x=Taxa)) + 
  geom_bar(stat="identity", position = position_stack(reverse=TRUE)) + 
  guides(fill = guide_legend(reverse=TRUE, nrow=6,byrow=TRUE)) + 
  scale_fill_manual(values = col, name = " ", labels = 
                      c("Salinity (mean = 28.6)","Temperature (mean = 12.5)",
                        #"Oxygen (mean = 22.9)",
                        "logDepth (mean = 9)", "Habitat",
                        "Random: Habitat (mean = 8.9)",
                        "Random: Station (mean = 18.3)")) +
  xlab ("Number of Taxa") + 
  ylab("Taxa niche explained variance partitioned")+  
  ggtitle("   a) Abundance") +
  theme(legend.position=c(0.8,0.93), legend.title = element_blank(),
        legend.text = element_text(size = 12), 
        legend.background = element_rect(fill="white"),
        panel.grid.major = element_line(colour = "gray70", linetype = 2),
        panel.grid.minor = element_line(colour = "gray85", linetype = 2),
        panel.background = element_rect(fill = "white"),
        #legend.box = "horizontal", 
        axis.text=element_text(size = 16),
        axis.title=element_text(size=18),#face="bold"),
        plot.title = element_text(face= 2, size=20), 
        legend.text.align = 1)  

Sp.var.part2_plot

# Look at means and SD's for reference------------------------------------------
# Random Ham =  8.9 -- 6.3%
# Ran sample = 18.3 -- 21
# oxy = 22,9 -- 5,7
# Salinity = 28.6 -- 10,7
# Temperature = 12,5 -- 7,9
# Log Depth = 9 -- 8,3

(sd(Sp.var.part$Random..habitat/Sp.var.part$R2))
(sd(Sp.var.part$Random..sample/Sp.var.part$R2))
(sd(Sp.var.part$Oxygen/Sp.var.part$R2))
(sd(Sp.var.part$B.PSU/Sp.var.part$R2))
(mean(Sp.var.part$B.Temp/Sp.var.part$R2))
(sd(Sp.var.part$logDepth/Sp.var.part$R2))




# Stacked: Biomass dataset------------------------------------------------------
#Sp.var.part2_bio <- Sp.var.part2
Sp.var.part2_bio
Sp.var.part2_bio_plot <- 
  
  ggplot(Sp.var.part2_bio, aes(fill=forcats::fct_infreq(Env), y=value, x=Taxa))+ 
  geom_bar(stat="identity", position = position_stack(reverse=TRUE)) + 
  guides(fill = guide_legend(reverse=TRUE, nrow=6,byrow=TRUE)) + 
  scale_fill_manual(values = col, name = " ", labels = 
                      c("Salinity (mean = 30.2)","Temperature (mean = 19.3)",
                        "Oxygen (mean = 21)", "logDepth (mean = 14.6)",
                        "Random: Habitat (mean = 6,5)",
                        "Random: Station (mean = 8,5)")) +
  xlab ("Number of Taxa") + 
  ylab(" ")+  
  ggtitle("   b) Biomass") +
  theme(legend.position=c(0.8,0.93), legend.title = element_blank(),
        legend.text = element_text(size = 12), 
        legend.background = element_rect(fill="white"),
        panel.grid.major = element_line(colour = "gray70", linetype = 2),
        panel.grid.minor = element_line(colour = "gray85", linetype = 2),
        panel.background = element_rect(fill = "white"),
        #legend.box = "horizontal", 
        axis.text=element_text(size = 16),
        axis.title=element_text(size=18),#face="bold"),
        plot.title = element_text(face= 2, size=20), 
        legend.text.align = 1)  

Sp.var.part2_bio_plot


#means for reference---------------------------------------------------
#Random Ham = 7.3
#Ran sample = 9
#oxy = 19.4
#Salinity = 30.4
# Temperature = 19.1
# Log Depth = 14.9

#plot altogether---------------------------------------------------------------
library(patchwork)

Var_plots <- Sp.var.part2_plot + Sp.var.part2_bio_plot

ggsave("Var_plots.png",
       plot = Var_plots,
       width = 17,
       height = 7,
       units = "in",
       dpi = 400,
       limitsize = TRUE)

Var_plots_long <- Sp.var.part2_plot / Sp.var.part2_bio_plot

ggsave("Var_plots_long.png",
       plot = Var_plots_long,
       width = 8,
       height = 14,
       units = "in",
       dpi = 400,
       limitsize = TRUE)


#End----------------------------------------------------------------------------
