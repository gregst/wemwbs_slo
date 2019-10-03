###################################
#
# R source code for data analysis in
# "VALIDATION OF THE WARWICK-
# EDINBURGH MENTAL WELL-BEING SCALE 
# (WEMWBS) AMONG UNIVERSITY STUDENT 
# POPULATION IN SLOVENIA"
#
###################################

# Packages
library("foreign")
library(psych)
for (n in c('ggplot2', 'corrplot', 'ggpubr', 'SparseM', 'foreign', 'utils', 'relimp', 'ggplot2', 'ggdendro', 'psych', 'Hmisc', 'ltm', 'mirt', 'eRm', 'mokken', 'lavaan','semTools','semPlot', 'qgraph','sem','CTT','MBESS','cluster')) {if(!require(n,character.only=TRUE)){install.packages(n)}}
library(n,character.only=TRUE)

# Reading the data
mydata= read.spss("db_wemwbs_slo.sav", use.value.labels=FALSE, to.data.frame=TRUE)
names(mydata)
dim(mydata)

#####################  STEP 1 - descriptive statistics ##################### 
# Sum of all WEMWBS items
attach(mydata)
WEMWBS_total <- WEMWBSa + WEMWBSb + WEMWBSc + WEMWBSd + WEMWBSe + WEMWBSf + WEMWBSg + WEMWBSh + WEMWBSi + WEMWBSj + WEMWBSk + WEMWBSl + WEMWBSm + WEMWBSn
mydata$WEMWBS_total <- WEMWBS_total
detach(mydata)

WEMWBS_table <- mydata[,31:44]
lowerCor(WEMWBS_table, method = "spearman")

# Psychosomatic items
mydata[,23:30] <- 6 - mydata[,23:30]
attach(mydata)
mydata$psychosom_total <- Headache + Stomach_pain + Back_pain + Miserable + Irritable_or_badmood + Nervous + Sleeping_troubles + Dizzy 
mydata$somatic <- Headache + Stomach_pain + Back_pain + Dizzy 
mydata$psychological <- Miserable + Irritable_or_badmood + Nervous + Sleeping_troubles 
detach(mydata)

# Correlations between WEMWBS, psychological and somatic symptoms and self-related health
lowerCor(mydata[,c("WEMWBS_total", "psychological", "somatic", "General_health")], method = "spearman")
p.mat <- cor.mtest(mydata[,c("WEMWBS_total", "psychological", "somatic", "General_health")])[["p"]]
p.mat

##################### STEP 2 - IRT analyses (the mokken, ltm, and mirt packages) ##################### 
moscales.for.lowerbounds <- function( x, lowerbounds=seq(from=0.05,to=0.60,by=0.05) )
{
  ret.value <- NULL;
  for( lowerbound in lowerbounds )
  {
    tmp <- aisp( x,  lowerbound=lowerbound );
    if( is.null(ret.value) )
    {
      ret.value <- data.frame( "Item"=rownames(tmp), "Scales."=tmp[,1] );
    }
    else
    {
      ret.value <- cbind( ret.value, "Scales."=tmp[,1] );
    }
    names(ret.value)[ncol(ret.value)] <- paste("c=",sprintf("%.2f",lowerbound),sep="");
  }
  rownames(ret.value) <- NULL;
  ret.value;
}

# Compute scalability coefficients
WEMWBS_table <- WEMWBS_table[complete.cases(WEMWBS_table),]
coefH(WEMWBS_table)$H

# examine aisp for increasing c levels (run the function you defined above and give it a name)
motable.WEMWBS_table <- moscales.for.lowerbounds( WEMWBS_table )

# see the results
motable.WEMWBS_table

# save it as a data frame
WEMWBS_table2 <- as.data.frame(motable.WEMWBS_table)

##################### STEP 3 - Parametric IRT #####################
# Rating Scale model (equivalent of Rasch for ordinal items)

WEMWBS_table2 <- WEMWBS_table[,-c(6,14)]
fit1.WEMWBS_table2 <- PCM(WEMWBS_table) #, constrained = FALSE, Hessian=TRUE

# separation reliability (proportion of item variance not due to error - similar to C-alpha)
ppr1 <- person.parameter(fit1.WEMWBS_table2)

# item fit (between 0.6 and 1.4 acc to Wright BD, Linacre JM. Reasonable mean-square fit values. Rasch Meas Trans. 1994;8(2):370.)
itemfit.fit1.WEMWBS_table2 <- itemfit(ppr1)

# check min and max infit and outfit
min(itemfit.fit1.WEMWBS_table2$i.infitMSQ)
max(itemfit.fit1.WEMWBS_table2$i.infitMSQ)

min(itemfit.fit1.WEMWBS_table2$i.outfitMSQ)
max(itemfit.fit1.WEMWBS_table2$i.outfitMSQ)

##################### STEP 4 - Confirmatory factor analysis ##################### 
# Exploratory factor analysis (EFA)

# factor analysis via parallel analysis
p1 <- fa.parallel(WEMWBS_table,cor="poly")
png(filename="figure3.png", type="cairo", height = 6, width = 6, units = 'in', res=300)
plot(p1)
dev.off()

# very simple structure analysis
vss(WEMWBS_table, 2)
# default FA - 1 factor, min residual & principal axis
fa(WEMWBS_table,  nfactors=1, fm="minres", n.iter=10)
fa(WEMWBS_table,  nfactors=1, fm="pa")
# plot the fa solution
plot(fa(WEMWBS_table,  nfactors=1, fm="pa"))

# hierarchical cluster analysis using ICLUST (groups items)
iclust(WEMWBS_table, title="WEMWBS_table using Pearson correlations")
summary(iclust(WEMWBS_table))
iclust.diagram(iclust(WEMWBS_table, title="WEMWBS_table using Pearson correlations"))
# hierarchical factor solution to find omega coefficient
omega(WEMWBS_table, nfactors=1, sl=FALSE)

# omega with polychoric matrix
WEMWBS_table.poly <- polychoric(WEMWBS_table)
omega(WEMWBS_table.poly$rho, nfactors=1,  sl=FALSE)

# CFA
# specify the model
CFA.BES <- "WEMWBS_total =~ WEMWBSa + WEMWBSb + WEMWBSc + WEMWBSd + WEMWBSe + WEMWBSf + WEMWBSg + WEMWBSh + WEMWBSi + WEMWBSj + WEMWBSk + WEMWBSl + WEMWBSm + WEMWBSn"

# fit the model
fitCFA.BES <- lavaan::cfa(CFA.BES, data=WEMWBS_table)
# model summary
summary(fitCFA.BES, standardized=TRUE, fit.measures = TRUE)
# coefficients only
coef(fitCFA.BES)
# CFA diagram from psych package
lavaan.diagram(fitCFA.BES, errors=TRUE)

##################### STEP 5 - CCT ##################### 
# CTT for a single scale

# Alpha by bootstrapping
ci.reliability(data=WEMWBS_table, type="alpha", conf.level = 0.95, interval.type="perc", B=100)

# Guttman lambda 6 (G6) and Beta values
splitHalf(WEMWBS_table) 

# Omega
ci.reliability(data=WEMWBS_table, type="omega", conf.level = 0.95, interval.type="perc", B=100)

##################### STEP 6 ##################### 
# check everything about your scores

# check descriptives
summary(WEMWBS_table)

# Histograms
png(filename="figure1.png", type="cairo", height = 6, width = 6, units = 'in', res=300)
hist(WEMWBS_total, breaks=40 , border=F , col=rgb(0.1,0.8,0.3,0.5) , xlab="distribution of WEMWBS_total" , main="")
dev.off()

# Plot png image
png(filename="figure2.png", type="cairo", height = 8, width = 8, units = 'in', res=300)
cor.plot(lowerCor(WEMWBS_table, method = "spearman"), numbers=TRUE, main="Correlations between WEMWBS items", 
         cex=0.5, cex.axis=0.7, xlas = 2)
dev.off()