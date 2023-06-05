################################################################################
############################## INDICES #########################################
################################################################################
rm(list=ls())
localDir = "."
indices.directory = file.path(localDir, "indices")

# load taxa data 
load("allData_A.R")
load("allData_B.R")

# assign new names
Phy <- P
Species <- Y
Traits <- Tr
Env <- X
S
str(Y)

# load libraries 
library(vegan)
library(devtools)
library(FD)
library(RColorBrewer)
library(adiv)


#Species indices----------------------------------------------------------------

Sp.index <- data.frame(Station = S$Station,
                       Habitat = S$Habitat,
                       Sp.rich = specnumber(Y), Total.ind = rowSums(Y),
                       Shannon = diversity(Y, index = "shannon"),
                      Simpson = diversity(Y, index = "simpson"))
Sp.index

#Functional indices-------------------------------------------------------------

ex3<-dbFD(Tr, Y, corr="sqrt", m=4) 
# m = min means low number of PCA axis. Loses information. But crashes otherwise.
#20 out of 40 axes were removed. quality = 0,89
ex3  

data<-data.frame(FRic = ex3$FRic, 
                 FEve = ex3$FEve, FDiv = ex3$FDiv, FDis =ex3$FDis,
                 RaoQ =ex3$RaoQ, Station = S$Station)

Sp.index_A <- merge(Sp.index, data, by = "Station")
Sp.index_B <- merge(Sp.index, data, by = "Station")

CWM <- ex3$CWM
write.csv(ex3$CWM,file = file.path(indices.directory, 'CWM_A.csv'))
write.csv(ex3$CWM,file = file.path(indices.directory, 'CWM_B.csv'))



#Functional redundancy----------------------------------------------------------

#make distance matrix
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

#Then run uniqueness function from adiv
Uni <- uniqueness(Species, Dis, abundance = TRUE)
Uni$red$Habitat <- S$Habitat
Uni$red

#And then Vulnerability
Vulnerability <- data.frame(Mean = colSums(Uni$V, na.rm = TRUE), 
                            SD = apply(Uni$V, 2, sd, na.rm=TRUE))
Vulnerability2 <- data.frame(Mean = colMeans(Uni$V, na.rm = T),
                             SD = apply(Uni$V, 2, sd, na.rm=TRUE))
boxplot(Uni$V)
Vulnerability

Sp.index_A$Vul1 <- Vulnerability$Mean
Sp.index_A$Vul2 <- Vulnerability2$Mean
Sp.index_A$Red <- Uni$red$R

Sp.index_B$Vul1 <- Vulnerability$Mean
Sp.index_B$Vul2 <- Vulnerability2$Mean
Sp.index_B$Red <- Uni$red$R

write.csv(Sp.index_A, file = file.path(indices.directory, "spIndices_A.csv"))
write.csv(Sp.index_B, file = file.path(indices.directory, "spIndices_B.csv"))

############ READ THE DATA IN FOR BOTH ABUNDANCE AND BIOMASS ####################
Sp.index_A <- read.csv(file ="indices/spIndices_A.csv", row.names = 1)
Sp.index_B <- read.csv(file ="indices/spIndices_B.csv", row.names = 1)

Sp.index_A$Data <- 'Abundance' 
Sp.index_B$Data <- "Biomass"
All_indices <- rbind(Sp.index_A, Sp.index_B)
All_indices$Habitat <- as.factor(All_indices$Habitat)

### PLOTTING THE DATA-----------------------------------------------------------
#Boxplots-----------------------------------------------------------------------
col_vec2 <- c("azure3", "cadetblue3", "aquamarine3", "antiquewhite3")

Fjord_indices<-subset(All_indices,Habitat=="Fjord")
Shelf_indices<-subset(All_indices,Habitat=="Shelf")
Shelfbreak_indices<-subset(All_indices,Habitat=="Shelfbreak")
Slope_indices<-subset(All_indices,Habitat=="Slope")

means_fjord <- apply(Fjord_indices[,c(3,5,6,7,8,9,11,12, 13, 14)], 2, mean)
means_shelf <- apply(Shelf_indices[,c(3,5,6,7,8,9,11,12, 13, 14)], 2, mean)
means_shelfbreak <- apply(Shelfbreak_indices[,c(3,5,6,7,8,9,11,12, 13, 14)], 2, mean)
means_slope <- apply(Slope_indices[,c(3,5,6,7,8,9,11,12, 13, 14)], 2, mean)

means <- cbind(Fjord = means_fjord, Shelf = means_shelf, Shelfbreak = means_shelfbreak, Slope= means_slope)
means <- as.data.frame(t(means))
means$Habitat <- rownames(means)
means


################ STATISTICS ON INDICES #################################
library("mgcv")
library("fields")
library("gamair")
library('car')
library("dplyr")
library("agricolae")
library("ggplot2")
library("ggthemes")
# Levene test: P has to be over 0.05.
#leveneTest(Sp.rich ~ Habitat, center = mean, data = All_indices)
leveneTest(log(Sp.rich) ~ Habitat, center = mean, data = All_indices)#ok.
#leveneTest(sqrt(Sp.rich) ~ Habitat, center = mean, data = All_indices)
leveneTest(FRic ~ Habitat, center = mean, data = All_indices)
leveneTest(Simpson ~ Habitat, center = mean, data = All_indices)
leveneTest(FDiv ~ Habitat, center = mean, data = All_indices)
#leveneTest(Vul2 ~ Habitat, center = mean, data = All_indices)
leveneTest(Vul2^2 ~ Habitat, center = mean, data = All_indices)
leveneTest(Red ~ Habitat, center = mean, data = All_indices)

anova1 <- aov(log(All_indices$Sp.rich) ~ All_indices$Habitat)
anova2 <- aov(All_indices$FRic ~ All_indices$Habitat)
anova3 <- aov(All_indices$Simpson ~ All_indices$Habitat)
anova4 <- aov(All_indices$FDiv ~ All_indices$Habitat)
anova5 <- aov(All_indices$Vul2^2 ~ All_indices$Habitat)
anova6 <- aov(All_indices$Red ~ All_indices$Habitat)

par(mfrow=c(2,2))
plot(anova1) # sp rich --------- is a good model 
plot(anova2) # F rich  --------- is a good model 
plot(anova3) # simpson --------- is a good model 
plot(anova4) # f diversity 
plot(anova5) # vul
plot(anova6) # red

summary(anova1)
summary(anova2) 
summary(anova3) # not sig.
summary(anova4)
summary(anova5)
summary(anova6)

hsd1 = HSD.test(aov(log(Sp.rich)~Habitat,data=All_indices), "Habitat", group=T)
hsd2 = HSD.test(aov(sqrt(FRic)~Habitat,data=All_indices), "Habitat", group=T)
hsd3 = HSD.test(aov(Simpson~Habitat,data=All_indices), "Habitat", group=T)
hsd4 = HSD.test(aov(FDiv~Habitat,data=All_indices), "Habitat", group=T)
hsd5 = HSD.test(aov(log(-All_indices$Vul2+1)~Habitat,data=All_indices), "Habitat", group=T)
hsd6 = HSD.test(aov(Red~Habitat,data=All_indices), "Habitat", group=T)

hsd1
hsd2
hsd3
hsd4
hsd5
hsd6

# plot
png("Tukeys test results 2.png", width = 6, height = 10, units = 'in', res = 300)
par(mfrow=c(3,2))
plot(TukeyHSD(aov(log(Sp.rich)~Habitat,data=All_indices), "Habitat", group=T))
plot(TukeyHSD(aov(sqrt(FRic)~Habitat,data=All_indices), "Habitat", group=T))
plot(TukeyHSD(aov(Simpson~Habitat,data=All_indices), "Habitat", group=T))
plot(TukeyHSD(aov(FDiv~Habitat,data=All_indices), "Habitat", group=T))
plot(TukeyHSD(aov(log(-All_indices$Vul2+1)~Habitat,data=All_indices), "Habitat", group=T))
plot(TukeyHSD(aov(Red~Habitat,data=All_indices), "Habitat", group=T))
dev.off()

# allocate sig. group letters
groups1 <- hsd1$groups[order(row.names(hsd1$groups)),]
groups2 <- hsd2$groups[order(row.names(hsd2$groups)),]
groups3 <- hsd3$groups[order(row.names(hsd3$groups)),]
groups4 <- hsd4$groups[order(row.names(hsd4$groups)),]
groups5 <- hsd5$groups[order(row.names(hsd5$groups)),]
groups6 <- hsd6$groups[order(row.names(hsd6$groups)),]

lmF_SPric <- (lm(All_indices$FRic ~ log(All_indices$Sp.rich)))
plot(lmF_SPric)
summary(lmF_SPric)
lmF_SPric$fitted.values


GLM_Fric_all <- glm(formula = All_indices$FRic ~ (log(All_indices$Sp.rich)), 
    family = (gaussian(link= log)))
summary(GLM_Fric_all)
100*(GLM_Fric_all$null.deviance-GLM_Fric_all$deviance)/GLM_Fric_all$null.deviance
#72% is explained by our model. 
chisqr <- GLM_Fric_all$null.deviance - GLM_Fric_all$deviance
GLM_Fric_all$df.null - GLM_Fric_all$df.residual #1
#using online chi sqr p value calculator. P=0.0000. 
0.65021^2

p <- predict(GLM_Fric_all, type="response", se.fit = TRUE)
p <- data.frame(p.fit = p$fit, p.se = p$se.fit)
p <- p[order(p$p.fit),]
xpred <- data.frame(Sp.ric = (All_indices$Sp.rich))
xpred <- xpred[order(xpred$Sp.ric),]



# PLOTS --------------------------------------
png("indices//All_indices.png", width=9, height= 8, units ='in', pointsize=10, 
    bg="white", res= 300)
par(mfrow=c(3,3), mar=c(4,5,2,1))
#Species & Funct richness-------------------------------------------------------
summarized = All_indices %>% group_by(Habitat) %>% summarize(Sp.rich=max(Sp.rich))
boxplot(All_indices$Sp.rich ~ All_indices$Habitat, col=brewer.pal(n = 4, name = "BrBG"),
        xlab= " ", ylab="Taxa Richness", ylim=c(0,4+(max(summarized$Sp.rich))),
        cex.axis = 1.2, cex.lab = 1.3)
title("a)", adj=0.01, line = 1, cex.main = 1.5)
points(means$Sp.rich, col = "red", pch = 18, cex=1)
over <- 3+summarized$Sp.rich
text(y=over, x=1:4, labels = groups1[,2], font = 3, cex = 1.2)
text(0.5,0.07, labels= "F = 29.0, p<0.001", col="red", font = 3, adj=0, cex = 1.2)

# Funtional Richness ------------------
summarized = All_indices %>% group_by(Habitat) %>% summarize(FRic=max(FRic))
boxplot(All_indices$FRic~All_indices$Habitat, col=brewer.pal(n = 4, name = "BrBG"),
        xlab= " ", ylab="Functional Richness", ylim = c(0,1200),
        cex.axis = 1.2, cex.lab = 1.3)
title("b)", adj=0.01, line = 1, cex.main = 1.5)
points(means$FRic, col = "red", pch = 18, cex=1)
over <- 75+summarized$FRic
text(y=over, x=1:4, labels = groups2[,2], font = 3, cex = 1.2)
text(0.5,1, labels= "F = 17.6, p<0.001", col="red", font = 3, adj=0, cex = 1.2)

#FRIC OVER SP RICH ----------------------------------------------
plot(All_indices$FRic ~ (All_indices$Sp.rich), 
     col="black", pch=19,cex.axis = 1.2, cex.lab = 1.3, cex=1.5,
     xlab=" ", ylab="Functional Richness", xlim=c(0,60), ylim=c(0,800))
grid()
points(jitter(All_indices$FRic[19:27]) ~jitter(All_indices$Sp.rich[19:27]), cex=1.5,
       col="black", bg ="darkgrey", pch=21, cex.axis = 1.2, cex.lab = 1.3)
title("c)", adj=0.01, line = 1, cex.main = 1.5)
lines(p$p.fit ~ xpred, col="black", lwd=2)
lines((p$p.fit)+(p$p.se) ~ xpred, col="grey20", lwd=1.5, lty=2)
lines(p$p.fit-p$p.se ~ xpred, col="grey20", lwd=1.5, lty=2)
text(1,750, labels= "T = 8.7, p<0.001", col="red", font = 3, adj=0, cex=1.2)

#Diversity----------------------------------------------------------------------
summarized = All_indices %>% group_by(Habitat) %>% summarize(Simpson=max(Simpson))
boxplot(All_indices$Simpson~ All_indices$Habitat, col=brewer.pal(n = 4, name = "BrBG"),
        xlab= " ", ylab="Taxa Diversity",ylim=c(0,0.1+(max(summarized$Simpson))),
        cex.axis = 1.2, cex.lab = 1.3)
title("d)", adj=0.01, line = 1, cex.main = 1.5)
points(means$Simpson, col = "red", pch = 18, cex=1)
#over <- 0.1+summarized$Simpson
#text(y=over, x=1:4, labels = groups3[,2], font = 3, cex = 0.7)
text(0.5,0.01, labels= "F = 2.0, p>0.05", col="grey40", font = 3, adj=0, cex = 1.2)


summarized = All_indices %>% group_by(Habitat) %>% summarize(FDiv=max(FDiv))
boxplot(All_indices$FDiv~All_indices$Habitat, col=brewer.pal(n = 4, name = "BrBG"),
        xlab= " ", ylab="Functional Diversity", ylim=c(0,0.1+(max(summarized$FDiv))),
        cex.axis = 1.2, cex.lab = 1.3)
title("e)", adj=0.01, line = 1, cex.main = 1.5)
points(means$FDiv, col = "red", pch = 18, cex=1)
over <- 0.05+summarized$FDiv
text(y=over, x=1:4, labels = groups4[,2], font = 3, cex = 1.2)
text(0.5,0.005, labels= "F = 3.7, p<0.05", col="red", font = 3, adj=0, cex = 1.2)


plot(All_indices$FDiv[1:18] ~ All_indices$Sp.rich[1:18], ylim=c(0.3,1.0),
     col="black", pch=19,cex.axis = 1.2, cex.lab = 1.3, cex=1.5,
     xlab=" ", ylab="Functional Diversity", xlim=c(0,60))
grid()
points(All_indices$FDiv[19:27] ~All_indices$Sp.rich[19:27],cex=1.5,
       col="black", bg ="darkgrey", pch=21, cex.axis = 1.2, cex.lab = 1.3)
title("f)", adj=0.01, line = 1, cex.main = 1.5)

#Vulnerability and Redundancy---------------------------------------------------
summarized = All_indices %>% group_by(Habitat) %>% summarize(Vul2=max(Vul2))
boxplot(All_indices$Vul2~All_indices$Habitat, col=brewer.pal(n = 4, name = "BrBG"), 
        ylab="Taxa Vulnerability", xlab = "Habitat type", ylim=c(0,1), #4+(max(summarized$Vul2))),
        cex.axis = 1.2, cex.lab = 1.3)
title("g)", adj=0.01, line = 1, cex.main = 1.5)
points(means$Vul2, col = "red", pch = 18, cex=1)
over <- 0.05+summarized$Vul2
text(y=over, x=1:4, labels = groups5[,2], font = 3, cex = 1.2)
text(0.5,0.005, labels= "F = 29.6, p<0.001", col="red", font = 3, adj=0, cex = 1.2)

summarized = All_indices %>% group_by(Habitat) %>% summarize(Red=max(Red))
boxplot(All_indices$Red ~ All_indices$Habitat, col=brewer.pal(n = 4, name = "BrBG"),
        xlab= "Habitat Type", ylab="Functional Redundancy", ylim=c(0,0.05+(max(summarized$Red))),
        cex.axis = 1.2, cex.lab = 1.3)
title("h)", adj=0.01, line = 1, cex.main = 1.5)
points(means$Red, col = "red", pch = 18, cex=1)
over <- 0.025+summarized$Red
text(y=over, x=1:4, labels = groups6[,2], font = 3, cex = 1.2)
text(0.5,0.008, labels= "F = 4.1, p<0.05", col="red", font = 3, adj=0, cex = 1.2)


plot(All_indices$Red[1:18] ~ All_indices$Sp.rich[1:18], 
     col="black", pch=19,cex.axis = 1.2, cex.lab = 1.3, cex=1.5,
     xlab="Number of Taxa", ylab="Functional Redundancy", xlim=c(0,60))
grid()
points(All_indices$Red[19:27] ~ All_indices$Sp.rich[19:27],cex=1.5,
       col="black", bg ="darkgrey", pch=21, cex.axis = 1.2, cex.lab = 1.3)
title("i)", adj=0.01, line = 1, cex.main = 1.5)

dev.off()


