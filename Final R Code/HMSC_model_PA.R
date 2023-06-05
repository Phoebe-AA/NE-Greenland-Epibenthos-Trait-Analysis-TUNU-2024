##############################################################################
############### MOVING ON TO FITTING THE MODEL S2 ############################
##############################################################################
library(Hmsc); library(corrplot); library("vioplot"); library("colorspace"); library(parallel)
#getwd()
load(file = "allData_A.R")
localDir = "."
#data.directory = file.path(localDir, "data")
model.directory = file.path(localDir, "models")

set.seed (1)
range(colMeans(Y>0))
hist(colMeans(Y>0), main = "prevalence")
hist(log(Y[Y>0]), main = "log abundance conditional on presence")
#Data is skewed towards low abundance values. 

#Many species are rare (absent in many samples) - need zero inflated. 
# Therefore we need a hurdle model. 

hist(rowSums(Y>0))
#Species richness across samples. 

names(X)
#Explore the correlations between the variables. 
plot(X)

XFormula <- ~ Depth + B.Temp + B.PSU + Oxygen + Turbidity  + Flourescence #  + Habitat

TrFormula <- ~ Small + Small_med + Medium_S + Med_large + Large +
Globulose + Vermiform + Dorso_com + Later_com +  Upright +
  Calcareous + Siliceous + Chitinous + Cuticle + None_SK + 
  Asexual + Sexual_ext + Sexual_int + Sexual_brood + Pel_Plankto +
  Pel_Leci + Benthic_Dir + Sessile + Burrower + Crawler +
  Swimmer+ none_MV + Low + medium_MV + High +
  Deposit + Filter_Sus + Op_Sc + Pred + Arctic +
  Boreal + Cosmo


studyDesign <- data.frame(plot = as.factor(S$Habitat))
#xy <- data.frame(S[match(unique(S$Habitat), S$Habitat), 2:3]) # takes the first value 
#rownames(xy) <- unique(S$Habitat) 
#colnames(xy) <-c("Lat", "Long")
#head(xy)

#set the random effect. 
#rL.plot <- HmscRandomLevel(sData = xy, longlat = TRUE)
#rL.plot <- HmscRandomLevel()
studyDesign = data.frame(sample = as.factor(1:nrow(S)), habitat = as.factor(S$Habitat)) 
rL = HmscRandomLevel(units = studyDesign$sample) 
rL2 = HmscRandomLevel(units = studyDesign$habitat)

Ypa <- 1*(Y>0)
Yabu <- Y
Yabu[Y==0] =NA
Yabu = log(Yabu)

m1 <- Hmsc(Y=Ypa, XData = X, XFormula = XFormula,
           TrData = Tr, TrFormula = TrFormula,
           phyloTree = P,
           studyDesign = studyDesign,
           ranLevels = {list("sample" = rL, "habitat" = rL2)}, 
           distr = "probit")


models <- list(m1)
modelnames <- c("presence_absence")

save(models, modelnames, file = file.path(model.directory, "unfitted_models"))

load(file = "models/unfitted_models")
models
modelnames

# set MCMC parameters

samples_list = 250
thin_list = 1000
nChains = 4
for(Lst in 1:length(samples_list)){
  thin = thin_list[Lst]
  samples = samples_list[Lst]
  
  print(paste0("thin = ", as.character(thin),"; samples = ", as.character(samples)))
  nm = length(models)
  
  for(model in 1:nm) {
    #model = 1
    print(paste0("model = ",modelnames[model]))
    m = models[[model]]
    m = sampleMcmc(m, samples = samples, thin = thin, 
                   adaptNf=rep(ceiling(0.4*samples*thin),m$nr),
                   transient = ceiling(0.5*samples*thin), 
                   nChains = nChains, nParallel = nChains)
    models[[model]] = m
  }
  filename = paste("models/models_thin_", as.character(thin),
                   "_samples_", as.character(samples),
                   "_chains_",as.character(nChains),
                   ".Rdata",sep = "")
  save(models,modelnames,file=filename)
}


