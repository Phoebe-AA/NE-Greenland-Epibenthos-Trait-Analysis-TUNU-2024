###############################################################################
################ NON Metric multidimensional scaling ##########################
###############################################################################
rm(list=ls())
localDir = "."
indices.directory = file.path(localDir, "NMDS")

# libraries 
library(pairwiseAdonis)
library(vegan)
library(RColorBrewer)

# load taxa data 
load("allData_A.R")
#load("allData_B.R")

# re-assign names
Phy <- P
Species <- Y
Traits <- Tr 
Env <- X
Env[order(Env$Depth),]
Env[order(Env$B.PSU),]
Env[order(Env$B.Temp),]
Env[order(Env$Flourescence),]
Env[order(Env$Oxygen),]
Env[order(Env$Turbidity),]
max(Env$Flourescence) - min(Env$Flourescence) # 0.09287

CWM <- read.csv(file = file.path('indices//CWM_A.csv'),
                header = T, row.names = 1)
#CWM <- read.csv(file = file.path('indices//CWM_B.csv'))

#?metaMDS
# RUN NMDS On species 
nmds_sp <- metaMDS(Species, trymax=999, distance = "bray", 
                   autotransform = TRUE) # bray curtis is the normal plot 
plot(nmds_sp, type ="p", display = "sites")

## Start from previous best solution
nmds_sp <- metaMDS(Species, previous.best = nmds_sp)
plot(nmds_sp, type ="t", display = "sites")

# Save scores and make sure
nmds.sp.scores<- as.data.frame(scores(nmds_sp))
nmds.sp.scores


# RUN NMDS on CWM trait values 
nmds_traits <- metaMDS(CWM, trymax=999, distance = "gower", 
                       autotransform = F) # 

plot(nmds_traits, type ="t", display = "sites")
## Start from previous best solution
nmds_traits <- metaMDS(CWM, previous.best = nmds_traits)
plot(nmds_traits, type ="t", display = "sites")

nmds.traits.scores <- as.data.frame(scores(nmds_traits))
nmds.traits.scores


########################### PLOTTING ###########################################
png("NMDS//NMDS plots_A.png", width=10, height=6, units="in", res = 300, bg="white")
par(mfrow=c(1,2))
col_vec2 <- c("azure3", "cadetblue3", "aquamarine3", "antiquewhite3")
# ELLIPSES FOR SPECIES 
plot(nmds.sp.scores, col=brewer.pal(n = 4, name = "BrBG")[Env$Habitat], 
     ylim=c(-1, 1), xlim=c(-1.5,1.5), 
     type='n',
     cex=1.2, main = "a) Abundance: Taxa Compositon Distance", xlab="NMDS1", ylab="NMDS2")
ordiellipse(nmds.sp.scores, groups=S$Habitat, 
            kind="ehull", col=brewer.pal(n = 4, name = "BrBG"), lwd=2)
ordiellipse(nmds.sp.scores, groups=S$Habitat, 
            draw="polygon", col=brewer.pal(n = 4, name = "BrBG"), lwd=2)
ordispider(nmds.sp.scores, groups=S$Habitat, 
           col=brewer.pal(n = 4, name = "BrBG"), label=TRUE)
points(nmds.sp.scores, pch=21, col="black", bg=brewer.pal(n = 4, name = "BrBG")[S$Habitat], cex=0.7)

# ELLIPSES FOR TRAITS 
plot(nmds.traits.scores, col=brewer.pal(n = 3, name = "BrBG")[S$Habitat], 
     ylim=c(-0.6, 0.6), xlim=c(-1.2,1.2), cex=1.2, 
     main = "b) Abundance: Trait CWM distance", xlab="NMDS1", ylab="NMDS2", type='n')
ordiellipse(nmds_traits, groups=S$Habitat, 
            kind="ehull", col=brewer.pal(n = 4, name = "BrBG"), lwd=2)
ordiellipse(nmds_traits, groups=S$Habitat, 
            draw="polygon", col=brewer.pal(n = 4, name = "BrBG"), lwd=2)
ordispider(nmds_traits, groups=S$Habitat, 
           col=brewer.pal(n = 4, name = "BrBG"), label=TRUE)
points(nmds_traits, pch=21, col="black", bg=brewer.pal(n = 4, name = "BrBG")[S$Habitat], cex=0.7)
dev.off()

###################################################################
################### STATISTICAL TESTING ###########################
###################################################################
################### ADONIS / PERMANOVA ############################
###################################################################

# Distance matrices 
distances.Species<- vegdist(Species, method = "bray")
distances.CWM<- vegdist(CWM, method = "gower")

plot(distances.Species, #col = Env$Habitat, 
     pch = 19)
plot(distances.CWM, #col = Env$Habitat, 
     pch = 19)
plot(distances.Species, #col = Env$Habitat, 
     pch = 19, ylim=c(0,1), ylab = "Dissimiliarities distances")
points(distances.CWM, #col = Env$Habitat, 
     pch = 19, col = "grey")

#Test for normality 
normal.test1 <- anova(betadisper(distances.Species, Env$Habitat)) 
normal.test1 # F value is much smaller. P = 39% of the time we will get greater than 1
# Conclude: the dispersion/variance of the data is equal 
normal.test2 <- anova(betadisper(distances.CWM, Env$Habitat))
normal.test2 # P = 19%
# Conclude: the dispersion/variance of the data is equal 


## ADONIS: tests for sig differences in communities among the habitat groups 
test.species.1 <- adonis(Species~Habitat,
                         data = Env, permutations = 999, method="bray")
test.species.1 

test.traits.1 <- adonis(CWM~Habitat,
                        data = Env, permutations = 999, method="gower")
test.traits.1
# Yes, there is a sig difference between groups for both species and traits. 

# Pairwise comparison 
A_sp_pairwise <- pairwise.adonis(Species, # change to any dataframe you want. 
                Env$Habitat, p.adjust.m='bonferroni')
A_sp_pairwise <- as.data.frame(A_sp_pairwise)
A_sp_pairwise
A_tr_pairwise <- pairwise.adonis(CWM, Env$Habitat,sim.method = "gower", p.adjust.m='bonferroni') # change to any dataframe you want. 
A_tr_pairwise <- as.data.frame(A_tr_pairwise)
A_tr_pairwise
A_pairwise <- rbind(A_sp_pairwise, A_tr_pairwise)
A_pairwise
write.csv(A_pairwise, "A_pairwise.csv")

# Repeat for biomass. 
# END ----------------------------------------------------------------------