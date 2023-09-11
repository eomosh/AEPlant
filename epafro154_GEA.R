setwd("D:/LANDSCAPE_GENOMICS/EggplantCWR/EggplantAfro")
getwd()
library(vegan)    # Used to run PCA & RDA
library(lfmm)     # Used to run LFMM
library(qvalue)   # Used to post-process LFMM output
library(vcfR)
library(LEA)
library(adegenet)


epafro_vcf = read.vcfR("D:/LANDSCAPE_GENOMICS/EggplantCWR/EggplantAfro/epafro_cur80_LD50_numer_snps.vcf", 
                       verbose = FALSE)
head(epafro_vcf)

epafro_genind = vcfR2genind(epafro_vcf)
epafro_gen_df = genind2df(epafro_genind)
head(epafro_gen_df)
dim(epafro_gen_df)

epafro_gen_df[,c(1:15146)] = sapply(epafro_gen_df[,c(1:15146)], as.numeric)
epafro_gen_df[is.na(epafro_gen_df)] = 9
epafro_gen_df[epafro_gen_df == 00 ] = 0
epafro_gen_df[epafro_gen_df == 11] = 1


epafro_gen_df = as.matrix(epafro_gen_df)
str(epafro_gen_df)
class(epafro_gen_df)
write.csv(epafro_gen_df, 'epafro_gen_df.csv', row.names = TRUE)
write.table(epafro_gen_df, "epafro_gen_df.txt", sep = "\t", col.names = TRUE, row.names = TRUE)


#Environmental data matrix
epafro_env_info = read.csv("D:/LANDSCAPE_GENOMICS/EggplantCWR/EggplantAfro/epafro_envqced.csv")
epafro_env_info
head(epafro_env_info)
str(epafro_env_info)
epafro_env_info$country = as.factor(epafro_env_info$country)
epafro_env_info$species = as.factor(epafro_env_info$species)
epafro_env_info1 = epafro_env_info[,-c(1:7)]
head(epafro_env_info1)


coord = epafro_env_info[,c("lon", "lat")]
str(coord)
coord_spd = SpatialPointsDataFrame(coords = coord, data = epafro_env_info, 
                                   proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84"))
coord_spd

#afro_data_comb = cbind(epafro_env_info, epafro_gen_df)
sp  = c("aethiopicum", "anguivi", "anomalum", "campylacanthum", "cerasiferum", "coagulans", "dasyphyllum",
        "incanum", "macrocarpon")
pdata = epafro_env_info[,c(4, 7:17)]
head(pdata)
pdata = subset(pdata, species %in% sp)

library(ggplot2)
p = ggplot(data = pdata, aes(x= PWaQ_18, y= MTWeQ_8, fill = species)) +
  geom_boxplot()
p


#Genomic matrix
# Reading genotype matrix (molecular matrix: independent variables):
Y <- epafro_gen_df
head(epafro_gen_df)
gen.pca = rda(Y, scale=TRUE)
screeplot(gen.pca, bstick = TRUE, type = "barplot")

Y1 = read.csv("Y_epafrosnps.csv")
row.names(Y1) = Y1$taxon
Y1 = Y1[, -1]
dim(Y1)
colnames(Y1) = colnames(Y)
head(Y1)

#Q matrix
str_matrix = read.csv("D:/LANDSCAPE_GENOMICS/EggplantCWR/EggplantAfro/K_structure.csv")
head(str_matrix)
strmat = str_matrix
str_matrix = as.data.frame(str_matrix)

#str_k3 = read.csv('epafrok3_matrix.csv')
#####Moran's eigenvector maps (MEMs)######
#MEMs are used to model the spatial auto correlations in RDA
library(adespatial)
library(geosphere)
library(spdep)

## load geographic coordinates
geo = epafro_env_info[,1:7]
head(geo)
## calculate mean geographic coordinates of samples from each deme
grp.coord <- cbind(tapply(geo$lon, INDEX = geo$country, mean), tapply(geo$lat, INDEX = geo$country, mean))

## 1. generate MEMs
## 1.1 make neighborhood network
nb <- graph2nb( gabrielneigh(grp.coord, nnmult = 100), sym=TRUE) # Gabriel graph: neighbor definition

## 1.2 produce spatial weights matrix
listW   <- nb2listw(nb,style="W")   

## 1.3 Calculate the distances between neighbours -> 'longlat = TRUE' for longitude-latitude decimal degrees data
disttri <- nbdists(nb,grp.coord, longlat = TRUE) 

## 1.4 Use inverse distance as weights
fdist   <- lapply(disttri,function(x) x^(-1))               

## 1.5 Revised spatial weights matrix
listW   <- nb2listw(nb,glist=fdist,style="W")               

## 1.6 Eigen analysis of W 
MEM <- scores.listw(listW, MEM.autocor = "all")         
rownames(MEM) <- rownames(grp.coord)

## 1.7 Match MEMs of demes to samples
mem.full <- as.data.frame(apply(MEM, 2, function(x){x[match(geo$country,rownames(MEM))]}))
rownames(mem.full) <- geo$taxon

## 2. Forward selection
fwsel.mem <- forward.sel(Y = as.matrix(Y), X = as.matrix(mem.full)) # significant if p-value < 0.05
head(mem.full)
# extract the significant MEMs (p < 0.05)
sel.mem <- mem.full[,colnames(mem.full) %in% fwsel.mem$variables]
write.csv(sel.mem, "selected_mem.csv", row.names = TRUE)
getwd()

#all factors (env_preds + q_matrix + genomic matrix)
all.factors = cbind(epafro_env_info,str_matrix,sel.mem)
head(all.factors)
str(all.factors)
write.csv(all.factors, "all.factors_epafro.CSV", row.names=FALSE)
all.factors= sapply(all.factors[,c(8:28)], as.numeric)

all.factors = data.frame(all.factors)
class(all.factors)


#####RDA Full model####
##Environment+structure+space

rda.full = rda(Y1 ~ PDrM_14 + PWaQ_18 + PCoQ_19 + MTWeQ_8 + srad + phh2o + 
                 nitrogen + clay + silt + ocd + K1 + K2 + K3+ K4 + K5 + K6 + K7 + K8 + MEM1 + MEM2 + MEM5, 
               data = all.factors, scale = T)
rda.full

#adjusted R2 
RsquareAdj(rda.full)

#permutation test
rda.full.anova <- anova(rda.full, permutations = 5000)
rda.full.anova

## env+structure
rda.env_qmat = rda(Y ~ as.matrix(epafro_env_info1) + as.matrix(str_matrix), scale = T)
rda.env_qmat

#adjusted R2 
RsquareAdj(rda.env_qmat)

#permutation test
rda.env_qmat.anova <- anova(rda.env_qmat, permutations = 5000)
rda.env_qmat.anova

## environment+space
rda.env_sp = rda(Y ~ as.matrix(epafro_env_info1) + as.matrix(sel.mem), scale = T)
rda.env_sp

#adjusted R2 
RsquareAdj(rda.env_sp)

#permutation test
rda.env_sp.anova <- anova(rda.env_sp, permutations = 5000)
rda.env_sp.anova


####PARTITIONING SNP VARIATION####
#Simple RDA

###Contribution of population structure
#Ancestry coefficients (K=8) are used to model population structure in the model.

#Fit an RDA model

rda.snp.qmat <- rda(as.matrix(Y) ~ as.matrix(str_matrix), scale = TRUE)
rda.snp.qmat

#adjusted R2 
RsquareAdj(rda.snp.qmat)

#permutation test

rda.snp.qmat.anova <- anova(rda.snp.qmat, permutations = 5000)
rda.snp.qmat.anova


###Contribution of spatial autocorrelation

rda.snp.sp <- rda(as.matrix(Y) ~ as.matrix(sel.mem), scale = TRUE)
rda.snp.sp

#adjusted R2 
RsquareAdj(rda.snp.sp)

#Permutation test

rda.snp.sp.anova <- anova(rda.snp.sp, permutations = 5000)

rda.snp.sp.anova


####Partitioning SNPs1####
##Partial RDAs 
rda.env <- rda(Y ~ as.matrix(epafro_env_info1) + Condition(as.matrix(str_matrix) + as.matrix(sel.mem)), all.factors, 
               scale = T)
rda.env
rda.qmat <- rda(Y ~ as.matrix(str_matrix) + Condition(as.matrix(epafro_env_info1) + as.matrix(sel.mem)),all.factors, 
                scale = T)
rda.qmat
rda.sel.mem <- rda(Y ~ as.matrix(sel.mem) + Condition(as.matrix(str_matrix) + as.matrix(epafro_env_info1)), all.factors, 
                   scale = T)
rda.sel.mem

#adjusted R2 
RsquareAdj(rda.env)
RsquareAdj(rda.qmat)
RsquareAdj(rda.sel.mem)

#Permutation test

rda.env.anova <- anova(rda.env, permutations = 5000)
rda.env.anova
rda.qmat.anova <- anova(rda.qmat, permutations = 5000)
rda.qmat.anova
rda.sel.mem.anova <- anova(rda.sel.mem, permutations = 5000)
rda.sel.mem.anova


####Partitioning SNPs2####
##Partial RDA controlling population structure and spatial autocorrelation
rda.snp.env.qmat <- rda(Y ~ as.matrix(epafro_env_info1) + as.matrix(str_matrix) + Condition(as.matrix(sel.mem)), 
                        scale = T)
rda.snp.env.qmat

#adjusted R2 
RsquareAdj(rda.snp.env.qmat)

#Permutation test

rda.snp.env.qmat.anova <- anova(rda.snp.env.qmat, permutations = 5000)
rda.snp.env.qmat.anova

####Contribution of all environmental variables####

#simple #Individual effect of each environmental variable#

s_env.anova <- lapply(1:ncol(epafro_env_info1),function(x){anova(rda(Y ~ epafro_env_info1[,x]), permutations = 5000)})
s_env.anova
names(s_env.anova) <- colnames(epafro_env_info1)

s_tol.var <- sum(s_env.anova[[1]]$Variance)
s_envstat <- cbind(sapply(s_env.anova, function(x){x$F[1]}),
                 sapply(s_env.anova, function(x){x$Variance[1]}), 
                 sapply(s_env.anova, function(x){x$Variance[1]/s_tol.var})*100,
                 sapply(s_env.anova, function(x){x$`Pr(>F)`[1]}))
rownames(s_envstat) <- colnames(epafro_env_info1)
colnames(s_envstat) <- c("F","Var","PerceVar", "pvalue")
rda.snp.env.anova.ind <- s_envstat[order(s_envstat[,1],decreasing = T),]
rda.snp.env.anova.ind


#conditioned on population structure

env.anova <- lapply(1:ncol(epafro_env_info1),function(x){anova(rda(Y ~ epafro_env_info1[,x] + 
             Condition(as.matrix(str_matrix))), permutations = 5000)} )

env.anova
names(env.anova) <- colnames(epafro_env_info1)

tol.var <- sum(env.anova[[1]]$Variance)
envstat <- cbind(sapply(env.anova, function(x){x$F[1]}),
                 sapply(env.anova, function(x){x$Variance[1]}), 
                 sapply(env.anova, function(x){x$Variance[1]/tol.var})*100,
                 sapply(env.anova, function(x){x$`Pr(>F)`[1]}))
rownames(envstat) <- colnames(epafro_env_info1)
colnames(envstat) <- c("F","Var","PerceVar", "pvalue")
rda.snp.env.qmat.anova.ind <- envstat[order(envstat[,1],decreasing = T),]
rda.snp.env.qmat.anova.ind

#Marginal effect of each environmental variable
rda.env.strmat = rda(Y ~ as.matrix(epafro_env_info1) + Condition(as.matrix(str_matrix)), scale=FALSE)


rda.env.strmat.anova.bymar <- anova.cca(rda.env.strmat, model = "direct",
                                          by = "margin", parallel = 4, permutations = how(nperm = 5000))


##Partial RDA controlling spatial autocorrelation
#Total effect of all environmental variables

rda.snp.env.sp <- rda(as.matrix(Y) ~ as.matrix(epafro_env_info1) + Condition(as.matrix(sel.mem)))
rda.snp.env.sp

#adjusted R2 
RsquareAdj(rda.snp.env.sp)

#Permutation test

rda.snp.env.sp.anova <- anova(rda.snp.env.sp, permutations = 5000)
rda.snp.env.sp.anova


####Partitioning population structure####

#As the genetic clusters correspond to their eco-geographic regions, it is also interesting to understand 
#the relative contribution of environmental gradients and spatial autocorrelation.

#To do so, we take ancestry coefficients as response variables and treat environmental variables and MEMs 
#as explanatory variables in RDA models.

#Simple RDA

#Contribution of environment + space

rda.qmat.env_sp <- rda(as.matrix(scale(str_matrix)) ~ as.matrix(epafro_env_info1) + as.matrix(sel.mem), scale = TRUE )
rda.qmat.env_sp

#adjusted R2 
RsquareAdj(rda.qmat.env_sp)

#Permutation test
rda.qmat.env_sp.anova <- anova(rda.qmat.env_sp, permutations = 5000)
rda.qmat.env_sp.anova


#Contribution of environmental variables

rda.qmatenv <- rda(as.matrix(scale(str_matrix)) ~ as.matrix(epafro_env_info1) )
rda.qmatenv

#adjusted R2 
RsquareAdj(rda.qmatenv)

#Permutation test
rda.qmatenv.anova <- anova(rda.qmatenv, permutations = 5000)
rda.qmatenv.anova

#Contribution of spatial autocorrelation

rda.qmatsp <- rda(as.matrix(scale(str_matrix)) ~ as.matrix(sel.mem))
rda.qmatsp

#adjusted R2 
RsquareAdj(rda.qmatsp)

#Permutation test

rda.qmatsp.anova <- anova(rda.qmatsp, permutations = 5000)
rda.qmatsp.anova


###Partial RDA###
#Contribution of environmental variables
#Fit a model controlling spatial autocorrelation

rda.qmatenv.sp <- rda(as.matrix(scale(str_matrix)) ~ as.matrix(epafro_env_info1) + Condition(as.matrix(sel.mem)) )
rda.qmatenv.sp

#adjusted R2 
RsquareAdj(rda.qmatenv.sp)

#Permutation test

rda.qmatenv.sp.anova <- anova.cca(rda.qmatenv.sp, permutations = 5000)
rda.qmatenv.sp.anova

#Contribution of spatial autocorrelation
#Fit a model controlling environmental variables

rda.qmatsp.env <- rda(as.matrix(scale(str_matrix)) ~ as.matrix(sel.mem) + Condition(as.matrix(epafro_env_info1)))
rda.qmatsp.env

#adjusted R2 
RsquareAdj(rda.qmatsp.env)


#Permutation test
rda.qmatsp.env.anova <- anova(rda.qmatsp.env, permutations = 5000)
rda.qmatsp.env.anova





####GEA####

## simple RDA###
library(vegan)

rda.snp.env <-rda(Y ~ PDrM_14 + PWaQ_18 + PCoQ_19 + MTWeQ_8 + srad + 
                     phh2o + nitrogen + clay + silt + ocd, data=epafro_env_info1, scale= TRUE)
rda.snp.env

#adjusted R2 
RsquareAdj(rda.snp.env) #R = 0.0873 (9% variation explained): R2 = 0.0231 (3% variation explained)

summary(eigenvals(rda.snp.env, model = "constrained"))

#statistical test
#To conduct statistical tests, we use the R function rdadapt (see github.com/Capblancq/RDA-genome-scan). 
#This function returns both p-values and q-values.

library(robust)
library(qvalue)

rdadapt<-function(rda,K)
{
  loadings<-rda$CCA$v[,1:as.numeric(K)]
  resscale <- apply(loadings, 2, scale)
  resmaha <- covRob(resscale, distance = TRUE, na.action= na.omit, estim="pairwiseGK")$dist
  lambda <- median(resmaha)/qchisq(0.5,df=K)
  reschi2test <- pchisq(resmaha/lambda,K,lower.tail=FALSE)
  qval <- qvalue(reschi2test)
  q.values_rdadapt<-qval$qvalues
  return(data.frame(p.values=reschi2test, q.values=q.values_rdadapt))
}

k=8
srda.env.pq = rdadapt(rda.snp.env, K=k)
View(srda.env.pq)
sum(srda.env.pq$q.values <= 0.05) 
write.csv(srda.env.pq, "pvalues_srda_env.csv", row.names = FALSE)


#Permutation test. The permutation will take some time to finish.

rda.snp.env.anova <- anova(rda.snp.env, permutations = 5000)
rda.snp.env.anova

vif.cca(rda.snp.env) #all env factor VIFs are below 10 and most below 5. this indicates no problem with multicollinearity

#Plot samples and the environmental predictors
#levels(epafro_env_info$species) = c("aethiopicum","macrocarpon", "dasyphyllum", "anomalum", "anguivi" , 
                                         # "cerasiferum", "incanum", "coagulans", "aculeastrum", "aculeatissimum",
                                         # "arundo", "campylacanthum", "sp", "dasyanthum" , "mauense" , "nigriviolaceum",
                                         # "phoxocarpum", "setaceum")
#eco = epafro_env_info$species
#bg = c('#FF9933', '#FFFF33','#33FF33','#33FF99','#33FFFF','#3399FF','#3333FF','#9933FF','#FF33FF',
                #'#FF3399','#A0A0A0','#660000','#663300', '#666600', '#336600','#000066','#330066','#660066')
                
levels(epafro_env_info$country) = c('GHA', "KEN" ,"NGA", "SDN" ,"TZA", "UGA")
eco = epafro_env_info$country
bg = c('#FF33FF','#33a02c','#1f78b4','#FF8000', '#ffff33','#e31a1c' )


# axes 1 & 2
plot(rda.snp.env, type="n", scaling=3, xlim=c(-5,5), ylim=c(-10,10), xlab= 'RDA1 (27.6%)', ylab= 'RDA2 (15.0%)', main= "Wild African eggplant species sRDA, axes 1 and 2")
points(rda.snp.env, display="species", pch=20, cex=1, col="gray32", scaling=3)           # the SNPs
points(rda.snp.env, display="sites", pch=21, cex=1, col="gray32", scaling=3, bg=bg[eco]) # the species
text(rda.snp.env, scaling=2, display="bp", col="#0868ac", cex=1)                           # the predictors
legend("bottomright", legend=levels(eco), bty="n", col="gray32", pch=21, cex=1, pt.bg=bg)

# axes 1 & 3
plot(rda.snp.env, type="n", scaling=3, choices=c(1,3))
points(rda.snp.env, display="species", pch=20, cex=0.7, col="gray32", scaling=3, choices=c(1,3))
points(rda.snp.env, display="sites", pch=21, cex=1.3, col="gray32", scaling=3, bg=bg, choices=c(1,3))
text(rda.snp.env, scaling=3, display="bp", col="#0868ac", cex=1, choices=c(1,3))
legend("topleft", legend=levels(eco), bty="n", col="gray32", pch=21, cex=1, pt.bg=bg)

#identify the candidates
load.rda <- summary(rda.snp.env)$species[,1:3] # load the first 3 most informative axes
hist(load.rda[,1], main="Loadings on RDA1")
hist(load.rda[,2], main="Loadings on RDA2")
hist(load.rda[,3], main="Loadings on RDA3") 


#the function for outliers here, where x is the vector of loadings and z is the number of standard deviations to use:
outliers <- function(x,z){
  lims <- mean(x) + c(-1, 1) * z * sd(x) ## f.nd loadings +/- z SD from mean loading     
  x[x < lims[1] | x > lims[2]]           # locus names in these tails
}

#Apply the function to the first three constrained axes
cand1 <- outliers(load.rda[,1], 3) ## 51
cand2 <- outliers(load.rda[,2], 3) ## 43
cand3 <- outliers(load.rda[,3], 3) ## 42

epafro.rda.cand <- c(names(cand1), names(cand2), names(cand3)) ## j.st the names of the candidates
epafro.rda.cand

length(epafro.rda.cand[duplicated(epafro.rda.cand)]) ## 2 duplicate detections (detected on multiple RDA axes)

ncand <- length(cand1) + length(cand2) + length(cand3)
ncand
#remove duplicates
epafro.rda.cand <- epafro.rda.cand[!duplicated(epafro.rda.cand)] ## 102 unique candidates 


#Make one data frame with the axis, SNP name, loading, & correlation with each predictor:

cand1 <- cbind.data.frame(rep(1,times=length(cand1)), names(cand1), unname(cand1))
cand2 <- cbind.data.frame(rep(2,times=length(cand2)), names(cand2), unname(cand2))
cand3 <- cbind.data.frame(rep(3,times=length(cand3)), names(cand3), unname(cand3))

colnames(cand1) <- colnames(cand2) <- colnames(cand3) <- c("axis","snp","loading")

sRDAcand <- rbind(cand1, cand2, cand3)
sRDAcand$snp <- as.character(sRDAcand$snp)

#Add in the correlations of each candidate SNP with the eight environmental predictors:

foo <- matrix(nrow=(ncand), ncol=10)  # 10 columns for 10 predictors
colnames(foo) <- c('PDrM_14' , 'PWaQ_18' , 'PCoQ_19' , 'MTWeQ_8' , 'srad' , 
                     'phh2o' , 'nitrogen' , 'clay' , 'silt' , 'ocd')

for (i in 1:length(sRDAcand$snp)) {
  nam <- sRDAcand[i,2]
  snp.gen <- Y[,nam]
  foo[i,] <- apply(epafro_env_info1,2,function(x) cor(x,snp.gen))
}

sRDAcand <- cbind.data.frame(sRDAcand,foo)  
head(sRDAcand)

#Duplicate SNPs detection
length(sRDAcand$snp[duplicated(sRDAcand$snp)])  # 2 duplicate detections

foo <- cbind(sRDAcand$axis, duplicated(sRDAcand$snp)) 
table(foo[foo[,1]==1,2]) # no duplicates on axis 1
table(foo[foo[,1]==2,2]) # no duplicates on axis 2
table(foo[foo[,1]==3,2]) # 6 duplicates on axis 3

sRDAcand <- sRDAcand[!duplicated(sRDAcand$snp),] # remove duplicate detections
sRDAcand


#Which of the predictors each candidate SNP is most strongly correlated with:

for (i in 1:length(sRDAcand$snp)) {
  bar <- sRDAcand[i,]
  sRDAcand[i,14] <- names(which.max(abs(bar[4:13]))) # gives the variable
  sRDAcand[i,15] <- max(abs(bar[4:12]))              # gives the correlation
}

colnames(sRDAcand)[14] <- "predictor"
colnames(sRDAcand)[15] <- "correlation"

table(sRDAcand$predictor) 
#write.csv(sRDAcand, 'sRDAcandVsVars.csv', row.names = FALSE)

#visualize candidate SNPs in the ordination space

sel <- sRDAcand$snp
env <- sRDAcand$predictor
env[env=="PDrM_14"] <- '#1f78b4'
env[env=='PWaQ_18'] <- '#FF33FF'
env[env=="PCoQ_19"] <- '#a6cee3'
env[env=="MTWeQ_8"] <- '#6a3d9a'
env[env=="srad"] <- '#e31a1c'
env[env=="phh2o"] <- '#33a02c'
env[env=="nitrogen"] <- '#ffff33'
env[env=="clay"] <- '#fb9a99'
env[env=="silt"] <- '#b2df8a'
env[env=="ocd"] <- '#FF8000'
                  
# color by predictor:
col.pred <- rownames(rda.snp.env$CCA$v) # pull the SNP names
                
for (i in 1:length(sel)) {           # color code candidate SNPs
foo <- match(sel[i],col.pred)
col.pred[foo] <- env[i]
}
                
col.pred[grep(':',col.pred)] <- '#f1eef6' # non-candidate SNPs
empty <- col.pred
empty[grep("#f1eef6",empty)] <- rgb(0,1,0, alpha=0) # transparent
empty.outline <- ifelse(empty=="#00FF0000","#00FF0000","gray32")
bg <- c('#1f78b4','#FF33FF','#a6cee3','#6a3d9a','#e31a1c','#33a02c','#ffff33','#fb9a99','#b2df8a', '#FF8000')
                  
#Plot the SNPs
# axes 1 & 2
plot(rda.snp.env, type="n", scaling=3, xlim=c(-1,1), ylim=c(-1,1), xlab= 'RDA1 (27.6%)', ylab= 'RDA2 (15.0%)', main= "Wild African eggplant species sRDA, axes 1 and 2")
points(rda.snp.env, display="species", pch=19, cex=1, col="gray", bg=col.pred, scaling=3)
points(rda.snp.env, display="species", pch=21, cex=1, col=empty.outline, bg=empty, scaling=3)
text(rda.snp.env, scaling=3, display="bp", col="#0868ac", cex=1)
legend("bottomright", legend=c('PDrM_14','PWaQ_18', 'PCoQ_19','MTWeQ_8','srad','phh2o','nitrogen', 'clay',
                                                 'silt', 'ocd'), bty="n", col="gray", pch=21, cex=1, pt.bg=bg)
                
# axes 1 & 3
plot(rda.snp.env, type="n", scaling=3, xlim=c(-1,1), ylim=c(-1,1), xlab= 'RDA1 (23.6%)', ylab= 'RDA3 (12.3%)', choices=c(1,3), main= "Wild eggplant RDA, axes 1 and 3")
points(rda.snp.env, display="species", pch=21, cex=1, col="gray32", bg=col.pred, scaling=3, choices=c(1,3))
points(rda.snp.env, display="species", pch=21, cex=1, col=empty.outline, bg=empty, scaling=3, choices=c(1,3))
text(rda.snp.env, scaling=3, display="bp", col="#0868ac", cex=1, choices=c(1,3))
legend("bottomright", legend=c('PDrM_14', 'PWaQ_18', 'PCoQ_19','MTWeQ_8','srad','phh2o','nitrogen', 'clay',
                                                 'silt', 'ocd'), bty="n", col="gray32", pch=21, cex=1, pt.bg=bg) 
##
summary(rda.snp.env)  #simple RDA summary


## partialRDA###
library(vegan)
str(all.factors)

prda.env <- rda(Y1 ~  PDrM_14 + PWaQ_18 + PCoQ_19 + MTWeQ_8 + srad  +  phh2o + nitrogen  + clay + silt + ocd 
                + Condition (K1 + K2 + K3 + K4 + K5 + K6 + K7 + K8 + MEM1 + MEM2 + MEM5), data= all.factors,
                        scale = T)
prda.env

#adjusted R2 
RsquareAdj(prda.env)#R = 0.0752 (8 variation explained): adjR2 = 0.0126 (1.3% variation explained)

#Permutation test

prda.anova <- anova(prda.env, permutations = 3000) #0.04832 *
prda.anova

summary(eigenvals(prda.env, model = "constrained")) #RDA 51% cummulative proportion (1=28.29, 2=11.81, 3= 10.9)

#statistical test
#To conduct statistical tests, we use the R function rdadapt (see github.com/Capblancq/RDA-genome-scan). 
#This function returns both p-values and q-values.

library(robust)
library(qvalue)

rdadapt<-function(rda,K)
{
  loadings<-rda$CCA$v[,1:as.numeric(K)]
  resscale <- apply(loadings, 2, scale)
  resmaha <- covRob(resscale, distance = TRUE, na.action= na.omit, estim="pairwiseGK")$dist
  lambda <- median(resmaha)/qchisq(0.5,df=K)
  reschi2test <- pchisq(resmaha/lambda,K,lower.tail=FALSE)
  qval <- qvalue(reschi2test)
  q.values_rdadapt<-qval$qvalues
  return(data.frame(p.values=reschi2test, q.values=q.values_rdadapt))
}

k=8
rda.env.pq = rdadapt(prda.env, K=k)
View(rda.env.pq)
sum(rda.env.pq$q.values <= 0.05) 
write.csv(rda.env.pq, "pvalues_prda_env.csv", row.names = FALSE)

#Plot samples and the environmental predictors
#levels(epafro_env_info$species) = c("aethiopicum","macrocarpon", "dasyphyllum", "anomalum", "anguivi" , 
                                    #"cerasiferum", "incanum", "coagulans", "aculeastrum", "aculeatissimum",
                                    #"arundo", "campylacanthum", "sp", "dasyanthum" , "mauense" , "nigriviolaceum",
                                    #"phoxocarpum", "setaceum"   )
#eco = epafro_env_info$species
#bg = c('#FF9933', '#FFFF33','#33FF33','#33FF99','#33FFFF','#3399FF','#3333FF','#9933FF','#FF33FF',
                #'#FF3399','#A0A0A0','#660000','#663300', '#666600', '#336600','#000066','#330066','#660066')

str(epafro_env_info)
levels(epafro_env_info$country) = c('GHA', "KEN" ,"NGA", "SDN" ,"TZA", "UGA")
eco = epafro_env_info$country
bg = c('#FF33FF','#33a02c','#1f78b4','#FF8000', '#ffff33','#e31a1c' )


# axes 1 & 2
plot(prda.env, type="n", scaling=3, xlab= 'RDA1 (28.29%)', ylab= 'RDA2 (11.8%)', main= "African eggplant CWR pRDA, axes 1 and 2")
points(prda.env, display="species", pch=20, cex=0.7, col="gray32", scaling=2)           # the SNPs
points(prda.env, display="sites", pch=21, cex=1.3, col="gray32", scaling=3, bg= bg[eco]) # the sample
text(prda.env, scaling=3, display="bp", col="#0868ac", cex=1)                           # the predictors
legend("topleft", legend=levels(eco), bty="n", col="gray32", pch=21, cex=1, pt.bg=bg)

# axes 1 & 3
plot(prda.env, type="n", scaling=3, choices=c(1,3))
points(prda.env, display="species", pch=20, cex=0.7, col="gray32", scaling=3, choices=c(1,3))
points(prda.env, display="sites", pch=21, cex=1.3, col="gray32", scaling=3, bg=bg[eco], choices=c(1,3))
text(prda.env, scaling=3, display="bp", col="#0868ac", cex=1, choices=c(1,3))
legend("topleft", legend=levels(eco), bty="n", col="gray32", pch=21, cex=1, pt.bg=bg)

#identify the candidates
load.rda <- summary(prda.env)$species[,1:3] # load the first 3 most informative axes
hist(load.rda[,1], main="Loadings on RDA1")
hist(load.rda[,2], main="Loadings on RDA2")
hist(load.rda[,3], main="Loadings on RDA3") 


#the function for outliers here, where x is the vector of loadings and z is the number of standard deviations to use:
outliers <- function(x,z){
  lims <- mean(x) + c(-1, 1) * z * sd(x) ## f.nd loadings +/- z SD from mean loading     
  x[x < lims[1] | x > lims[2]]           # locus names in these tails
}

#Apply the function to the first three constrained axes
cand1 <- outliers(load.rda[,1], 3) ## 37
cand2 <- outliers(load.rda[,2], 3) ## 35
cand3 <- outliers(load.rda[,3], 3) ## 61

prda.env.cand <- c(names(cand1), names(cand2), names(cand3)) ## just the names of the candidates
prda.env.cand

length(prda.env.cand[duplicated(prda.env.cand)]) ## 1 duplicate detections (detected on multiple RDA axes)

ncand <- length(cand1) + length(cand2) + length(cand3)
ncand

#remove duplicates
prda.env.cand <- prda.env.cand[!duplicated(prda.env.cand)] ## 132 unique candidates 
prda.env.cand

#Make one data frame with the axis, SNP name, loading, & correlation with each predictor:

cand1 <- cbind.data.frame(rep(1,times=length(cand1)), names(cand1), unname(cand1))
cand2 <- cbind.data.frame(rep(2,times=length(cand2)), names(cand2), unname(cand2))
cand3 <- cbind.data.frame(rep(3,times=length(cand3)), names(cand3), unname(cand3))

colnames(cand1) <- colnames(cand2) <- colnames(cand3) <- c("axis","snp","loading")

cand <- rbind(cand1, cand2, cand3)
cand$snp <- as.character(cand$snp)

#Add in the correlations of each candidate SNP with the eight environmental predictors:

foo <- matrix(nrow=(ncand), ncol=10)  # 10 columns for 10 predictors
colnames(foo) <- c('PDrM_14' , 'PWaQ_18' , 'PCoQ_19' , 'MTWeQ_8' , 'srad' , 
                   'phh2o' , 'nitrogen' , 'clay' , 'silt' , 'ocd')

for (i in 1:length(cand$snp)) {
  nam <- cand[i,2]
  snp.gen <- Y[,nam]
  foo[i,] <- apply(epafro_env_info1,2,function(x) cor(x,snp.gen))
}

cand <- cbind.data.frame(cand,foo)  
head(cand)

#Duplicate SNPs detection
length(cand$snp[duplicated(cand$snp)])  # 1 duplicate detections

foo <- cbind(cand$axis, duplicated(cand$snp)) 
table(foo[foo[,1]==1,2]) # no duplicates on axis 1
table(foo[foo[,1]==2,2]) # no duplicates on axis 2
table(foo[foo[,1]==3,2]) # 1 duplicates on axis 3

cand <- cand[!duplicated(cand$snp),] # remove duplicate detections
cand
View(cand)
#write.csv( cand, "pRDA_Candidates.CSV", row.names = FALSE)

#Which of the predictors each candidate SNP is most strongly correlated with:


for (i in 1:length(cand$snp)) {
  bar <- cand[i,]
  cand[i,14] <- names(which.max(abs(bar[4:13]))) # gives the variable
  cand[i,15] <- max(abs(bar[4:12]))              # gives the correlation
}

colnames(cand)[14] <- "predictor"
colnames(cand)[15] <- "correlation"

table(cand$predictor) 
#write.csv(cand, "pRDAcandidatesVsvars.csv", row.names = FALSE)
#visualize candidate SNPs in the ordination space

sel <- cand$snp
env <- cand$predictor
env[env=="PDrM_14"] <- 'navyblue'
env[env=='PWaQ_18'] <- '#FF33FF'
env[env=="PCoQ_19"] <- 'dodgerblue'
env[env=="MTWeQ_8"] <- '#6a3d9a'
env[env=="srad"] <- '#e31a1c'
env[env=="phh2o"] <- '#33a02c'
env[env=="nitrogen"] <- 'yellow'
env[env=="clay"] <- 'hotpink'
env[env=="silt"] <- '#b2df8a'
env[env=="ocd"] <- '#FF8000'

          
# color by predictor:
col.pred <- rownames(prda.env$CCA$v) # pull the SNP names
            
for (i in 1:length(sel)) {           # color code candidate SNPs
foo <- match(sel[i],col.pred)
col.pred[foo] <- env[i]
}
                  
col.pred[grep(':',col.pred)] <- 'gray88' # non-candidate SNPs
empty <- col.pred
empty[grep("gray88",empty)] <- rgb(0,1,0, alpha=0) # transparent
empty.outline <- ifelse(empty=="#00FF0000","#00FF0000","gray32")
bg <- c('navyblue','#FF33FF','dodgerblue','#6a3d9a','#e31a1c','#33a02c','yellow','hotpink',
        '#b2df8a', '#FF8000')
                   
#Plot the SNPs
# axes 1 & 2
plot(prda.env, type="n", scaling=3, xlim=c(-1,1), ylim=c(-1,1), xlab= 'RDA1 (28.29%)', ylab= 'RDA2 (11.8%)', main= "Wild African eggplant species pRDA, axes 1 and 2")
points(prda.env, display="species", pch=19, cex=1, col="gray", bg=col.pred, scaling=3)
points(prda.env, display="species", pch=21, cex=1, col=empty.outline, bg=empty, scaling=3)
text(prda.env, scaling=3, display="bp", col="#0868ac", cex=1)
legend("bottomright", legend=c('PDrM_14','PWaQ_18', 'PCoQ_19','MTWeQ_8','srad','phh2o','nitrogen', 'clay',
          'silt', 'ocd'), bty="n", col="gray", pch=21, cex=1, pt.bg=bg)
                    
# axes 1 & 3
plot(prda.env, type="n", scaling=3, xlim=c(-1,1), ylim=c(-1,1), xlab= 'RDA1 (23.6%)', ylab= 'RDA3 (10.9%)', choices=c(1,3), main= "Wild eggplant RDA, axes 1 and 3")
points(prda.env, display="species", pch=21, cex=1, col="gray32", bg=col.pred, scaling=3, choices=c(1,3))
points(prda.env, display="species", pch=21, cex=1, col=empty.outline, bg=empty, scaling=3, choices=c(1,3))
text(prda.env, scaling=3, display="bp", col="#0868ac", cex=1, choices=c(1,3))
legend("bottomright", legend=c('PDrM_14', 'PWaQ_18', 'PCoQ_19','MTWeQ_8','srad','phh2o','nitrogen', 'clay',
         'silt', 'ocd'), bty="n", col="gray32", pch=21, cex=1, pt.bg=bg) 
                  
##
summary(prda.env)



##########################

pRDA_analysis_results = cand %>% dplyr::select(predictor, snp, axis, loading)
head(pRDA_analysis_results)
write.csv(pRDA_analysis_results, "pRDA_analysis_results.csv", row.names = FALSE)

rda_outliers_color <- function(rda_results = NULL) {
  snp_colors <- rep(NA, nrow(rda_results))
  colors <- colorRampPalette(c('#1f78b4','#FF33FF','#a6cee3','#6a3d9a','#e31a1c','#33a02c','yellow','#fb9a99',
                               '#b2df8a', '#FF8000'))
  unique_env <- unique(rda_results$predictor)
  my_col <- colors(length(unique_env))
  for (i in 1:length(unique_env)) {
    snp_colors[which(rda_results$predictor == unique_env[i])] <- my_col[i]
  }
  return(snp_colors)
}

cand_col <- rda_outliers_color(rda_results = pRDA_analysis_results)
cand_col

# rda_neutral_color()
rda_biplot <- function(rda_analysis = NULL, rda_analysis_results, 
                       candidate_color = NULL, neutral_color = "gray70", legend = FALSE) {
  snp_col <- rownames(rda_analysis$CCA$v)
  names(snp_col) <- snp_col
  snp_col[1:length(snp_col)] <- alpha(neutral_color, 0.4)
  
  plot(prda.env, type="n", scaling=3, xlim=c(-1,1), ylim=c(-1,1), xlab= 'RDA1 (28.29%)', 
       ylab= 'RDA2 (11.8%)', main= "African eggplant CWR RDA, axes 1 and 2")
  points(prda.env, display="species", pch=16, cex=1, col=snp_col, scaling=3)
  
  snp_col[1:length(snp_col)] <- rgb(0,1,0, alpha=0)
  for (i in 1:length(pRDA_analysis_results$snp)) {
    snp_col[match(pRDA_analysis_results$snp[i], names(snp_col))] <- candidate_color[i]
  }
  
  points(prda.env, display="species", pch=16, cex=1, col=snp_col, scaling=3)
  text(prda.env, scaling=3, display="bp", col="#0868ac", cex=1)
  tmp <- unique(cbind.data.frame(pRDA_analysis_results$predictor, candidate_color))
  colnames(tmp) <- c("ENV", "COL")
  if (legend) {
    legend("topleft", legend=tmp$ENV, bty="n", pch=16, cex=0.8, col=tmp$COL, bg="white")
  }
  
}

rda_biplot(rda_analysis = prda.env, 
           rda_analysis_results = pRDA_analysis_results, 
           candidate_color = cand_col, legend = TRUE)

##LFMM####
                  
library(lfmm)
fit.lfmm <- lfmm_ridge(Y = Y, X = epafro_env_info1, K = 8)
              
pv <- lfmm_test(Y = Y, X = epafro_env_info1, lfmm = fit.lfmm,calibrate = "gif")
                  
pvalues <- pv$calibrated.pvalue 
write.csv(pvalues, "lfmm_pvalues.csv", row.names = FALSE)                  
                  
# check if the distribution of p-value agree with the FDR assumption
hist(as.vector(pvalues), xlab = "Adjusted P-values", main = "Distribution of Adjusted P-values of LFMM (K = 8)")
                  
                  
# calculate q-value for each GEA test respectively
qvalues <-apply(pvalues,2, function(x){qvalue(x)$qvalues})
rownames(qvalues) <- rownames(pvalues)
colnames(qvalues) <- colnames(pvalues)
write.csv(qvalues, "lfmm_qvalues.csv", row.names = FALSE)                  
#Number of significant SNPs
#For all environmental variables
                  
sum(apply(qvalues, 1, function(x){any(x < 0.05)})) #There are 51 significant SNPs in total for all environmental variables.
                  
#For each environmental variable
                  
#We further look at the number of significant SNPs for each environmental variable.
                  
qvalues1 = qvalues
qvalues1 = as.data.frame(qvalues1)
length(which(qvalues1$PDrM_14 < 0.05)) #34
length(which(qvalues1$PWaQ_19 < 0.05)) #0
length(which(qvalues1$PCoQ_19 < 0.05)) #1
length(which(qvalues1$MTWeQ_8 < 0.05))#1
length(which(qvalues1$srad < 0.05)) #1
length(which(qvalues1$phh2o < 0.05))    #0
length(which(qvalues1$nitrogen < 0.05))#12
length(which(qvalues1$clay < 0.05))     #2
length(which(qvalues1$silt < 0.05))    #0
length(which(qvalues1$ocd < 0.05))    #0
                  
                  
PDrM_snps = which(qvalues1$PDrM_14 < 0.05)
PCoQ_snps = which(qvalues1$PCoQ_19 < 0.05)
MTWeQ_snps = which(qvalues1$MTWeQ_8 < 0.05)
srad_snps =which(qvalues1$srad < 0.05)
nitrogen_snps = which(qvalues1$nitrogen < 0.05)
clay_snps = which(qvalues1$clay < 0.05)

                  
(lfmm.snps1 <- row.names(qvalues1)[which(qvalues1$PDrM_14 < 0.05)])
(lfmm.snps2 <- row.names(qvalues1)[which(qvalues1$PCoQ_19 < 0.05)])
(lfmm.snps2 <- row.names(qvalues1)[which(qvalues1$MTWeQ_8 < 0.05)])
(lfmm.snps3 <- row.names(qvalues1)[which(qvalues1$srad < 0.05)])
(lfmm.snps3 <- row.names(qvalues1)[which(qvalues1$nitrogen < 0.05)])
(lfmm.snps4 <- row.names(qvalues1)[which(qvalues1$clay < 0.05)])
                  
PDrM_14snps = which(qvalues1$PDrM_14 < 0.05)
PDrM_14snps
PDrMsnps_colnames = subset(Y, TRUE, c(PDrM_14snps))
PDrM_outliers = colnames(PDrMsnps_colnames)
PDrM_outliers

PCoQ_snps = which(qvalues1$PCoQ_19 < 0.05)
PCoQ_snps
PCoQsnps_colnames = subset(Y, TRUE, c( PCoQ_snps))
PCoQ_outliers = colnames(PCoQsnps_colnames)
PCoQ_outliers

MTWeQ_snps = which(qvalues1$MTWeQ_8 < 0.05)
MTWeQ_snps
MTWeQsnps_colnames = subset(Y, TRUE, c( MTWeQ_snps))
MTWeQ_outliers = colnames(MTWeQsnps_colnames)
MTWeQ_outliers

srad_snps = which(qvalues1$srad < 0.05)
srad_snps
srad_colnames = subset(Y, TRUE, c(srad_snps))
srad_outliers = colnames(srad_colnames)
srad_outliers 

nitrogen_snps = which(qvalues1$nitrogen < 0.05)
nitrogen_snps
nitrogen_colnames = subset(Y, TRUE, c(nitrogen_snps))
nitrogen_outliers = colnames(nitrogen_colnames)
nitrogen_outliers 

clay_snps = which(qvalues1$clay < 0.05)
clay_snps
clay_colnames = subset(Y, TRUE, c(clay_snps))
clay_outliers = colnames(clay_colnames)
clay_outliers  



lfmm_outliers = c(PDrM_outliers,PCoQ_outliers, MTWeQ_outliers, srad_outliers, nitrogen_outliers, clay_outliers  )

lfmmlist = tibble::lst(PDrM_outliers,PCoQ_outliers, MTWeQ_outliers, srad_outliers, nitrogen_outliers, clay_outliers )
lfmmoutlier_df = data.frame(lapply(lfmmlist,`length<-`, max(lengths(lfmmlist))))
write.csv(lfmmoutlier_df, "lfmmoutlier_df.csv", row.names = FALSE)

####PCAdapt####
library(pcadapt)
genotype = Y

pca.genotype = read.pcadapt(t(genotype))
class(genotype)                
K = 12
x <- pcadapt(pca.genotype #, K = K
             )
plot(x, option = "screeplot") # 2 groups seems to be the correct value
                  
plot(x, option = "scores", pop = epafro_env_info$species) # how populations are shared among the 18 groups
                  
K <- 2
x <- pcadapt(pca.genotype #, K = K
             )
                  
summary(x) # numerical quantities obtained after performing a PCA
plot(x, option = "screeplot")
#plot(x, option = "manhattan")
#plot(x, option = "qqplot", threshold = 0.05)
plot(x, option = "stat.distribution") # Distribution of Mahalanobis distances.
                  
#padj1 = x$pvalues                  
padj <- p.adjust(x$pvalues,method="bonferroni")
write.csv(padj, "pcadapt_pvalues.csv", row.names = FALSE)
#write.csv(padj1, "pcadapt_pvalues without correction.csv", row.names = FALSE)
alpha <- 0.05
outliers_pcadapt <- which(padj < alpha)
length(outliers_pcadapt) #114 SNPs out of 15146
print(outliers_pcadapt)
                 
loc1 <- genotype[, outliers_pcadapt[1]]
#loc1 = loc1[(loc1 == 1) | (loc1 ==0)]
PDrM = epafro_env_info$PDrM_14
lociPDrM = data.frame(loc1, PDrM)                 
lociPDrM = lociPDrM[(lociPDrM$loc1 == 1) | (lociPDrM$loc1 == 0),]  
PDrM = lociPDrM$PDrM

mod <- glm(cbind(lociPDrM$loc1, 2 - lociPDrM$loc1) ~ PDrM, family = binomial) 
summary(mod) # This locus is significantly correlated to temperature
                  
outlier_names = subset(Y1, TRUE, c(outliers_pcadapt)) 
                  
outliers_pcadapt = colnames (outlier_names)
                  
#### all outliers####
lfmm_outliers   # 51   
epafro.rda.cand #130
prda.env.cand   #132
outliers_pcadapt #114
           

intersect(lfmm_outliers, outliers_pcadapt)#0
intersect(lfmm_outliers, epafro.rda.cand)#3
intersect(lfmm_outliers, prda.env.cand)#10
intersect(prda.env.cand, epafro.rda.cand)#14
intersect(prda.env.cand, outliers_pcadapt)#1
intersect(outliers_pcadapt, epafro.rda.cand)#3
Reduce(intersect, list(lfmm_outliers, epafro.rda.cand, prda.env.cand ))


all_outliers = c(lfmm_outliers, prda.env.cand, epafro.rda.cand, outliers_pcadapt) #427 outliers
all_outliers
write.csv(all_outliers, "all_outliers.csv")

gea_outliers = c(lfmm_outliers, prda.env.cand, epafro.rda.cand) #313 outliers
write.csv(gea_outliers,"gea_outliers.csv", row.names = FALSE)
          
#### Venn diagram ####

library("VennDiagram")
# all outliers
lfmm_outliers   # 51   
epafro.rda.cand #130
prda.env.cand   #132
outliers_pcadapt #114

intersect(lfmm_outliers, epafro.rda.cand)#3
intersect(prda.env.cand, epafro.rda.cand)#14
intersect(lfmm_outliers, prda.env.cand)#10
intersect(lfmm_outliers, outliers_pcadapt)#0
intersect(epafro.rda.cand,outliers_pcadapt )#0
intersect(prda.env.cand, outliers_pcadapt)#0

Reduce(intersect, list(lfmm_outliers, epafro.rda.cand, prda.env.cand ))#0
Reduce(intersect, list(lfmm_outliers, epafro.rda.cand, outliers_pcadapt ))#0
Reduce(intersect, list(epafro.rda.cand, prda.env.cand,outliers_pcadapt ))#0
Reduce(intersect, list(lfmm_outliers, prda.env.cand, outliers_pcadapt ))#0
Reduce(intersect, list(lfmm_outliers, epafro.rda.cand, prda.env.cand,outliers_pcadapt ))#0

all_outliers = c(lfmm_outliers, prda.env.cand, epafro.rda.cand, outliers_pcadapt) #427 outliers
all_outliers

write.csv(all_outliers, "all_outliers.csv")
# move to new plotting page
grid.newpage()
# create Venn diagram with three sets

dev.off()                 
draw.quad.venn(area1=51, area2=130, area3=132, area4=114,
n12=2, n23=14, n13=9, n14=0 ,n24=3, n34=0, n123=0, 
n124=0  , n234=0 , n134=0 , n1234=0 ,
category=c("LFMM","sRDA","pRDA","PCAdapt"),
col="gray24",fill=c("#1f78b4",'#b2df8a','#FF8000','#33a02c'), lty = 'dashed',
cex = 2,
cat.cex = 2)


####Manhattan plots####

library(qqman)
pvals = read.table("epafro_GEA_pvalues.txt", sep = '\t', header = TRUE)
pvals = pvals[-(1:223),] 
#pvals = pvals[224:15146,]
head(pvals)

#l
lfmm_outliers   # 51 

a = pvals[pvals$PDrM_14 < 0.05,]
b = pvals[pvals$PWaQ_18 < 0.05,]
c = pvals[pvals$PCoQ_19 < 0.05,]
d = pvals[pvals$MTWeQ_8 < 0.05,]
e = pvals[pvals$srad < 0.05,]
f = pvals[pvals$phh2o < 0.05,]
g = pvals[pvals$nitrogen < 0.05,]
h = pvals[pvals$clay < 0.05,]
i = pvals[pvals$silt < 0.05,]
j = pvals[pvals$ocd < 0.05,]

head(a)
NROW(a)
int_a = a$ID #34
int_b = b$ID #0
int_c = c$ID #1
int_d = d$ID #1
int_e = e$ID #1
int_f = f$ID #0
int_g = g$ID #12
int_h = h$ID #2
int_i = i$ID #0
int_j = j$ID #0
 
int_lfmm = c (int_a, int_b,int_c, int_d, int_e, int_f, int_g, int_h, int_i, int_j)

lfmm_man_PDrM = manhattan(pvals, chr = "CHROM", bp = "POS", p= 'PDrM_14', snp = "ID",col = c("grey20", "grey60"), 
                     logp = TRUE, #ylab = "lfmm_PDrM", 
                     xlab = 'CHROM', highlight = int_a, suggestiveline = FALSE,
                     genomewideline = FALSE, annotatePval = FALSE)

lfmm_man_PCoQ = manhattan(pvals, chr = "CHROM", bp = "POS", p= 'PCoQ_19', snp = "ID",col = c("grey20", "grey60"), 
                      logp = TRUE, #ylab = "lfmm_PDrM", 
                      xlab = 'CHROM', highlight = int_c, suggestiveline = FALSE,
                      genomewideline = FALSE, annotatePval = FALSE)

lfmm_man_MTWeQ = manhattan(pvals, chr = "CHROM", bp = "POS", p= 'MTWeQ_8', snp = "ID",col = c("grey20", "grey60"), 
                      logp = TRUE, #ylab = "lfmm_PDrM", 
                      xlab = 'CHROM', highlight = int_d, suggestiveline = FALSE,
                      genomewideline = FALSE, annotatePval = FALSE)

lfmm_man_srad = manhattan(pvals, chr = "CHROM", bp = "POS", p= 'srad', snp = "ID",col = c("grey20", "grey60"), 
                      logp = TRUE, #ylab = "lfmm_PDrM", 
                      xlab = 'CHROM', highlight = int_e, suggestiveline = FALSE,
                      genomewideline = FALSE, annotatePval = FALSE)

lfmm_man_N = manhattan(pvals, chr = "CHROM", bp = "POS", p= 'nitrogen', snp = "ID",col = c("grey20", "grey60"), 
                      logp = TRUE, #ylab = "lfmm_PDrM", 
                      xlab = 'CHROM', highlight = int_g, suggestiveline = FALSE,
                      genomewideline = FALSE, annotatePval = FALSE)

lfmm_man_clay = manhattan(pvals, chr = "CHROM", bp = "POS", p= 'clay', snp = "ID",col = c("grey20", "grey60"), 
                      logp = TRUE, #ylab = "lfmm_PDrM", 
                      xlab = 'CHROM', highlight = int_h, suggestiveline = FALSE,
                      genomewideline = FALSE, annotatePval = FALSE)

#srda manhattan
epafro.rda.cand #130

srda_snps = pvals[pvals$srda_p < 0.0027,]
head(srda_snps)
NROW(srda_snps)
int_srda = srda_snps$ID

srda_man = manhattan(pvals, chr = "CHROM", bp = "POS", p= 'srda_p', snp = "ID",col = c("grey20", "grey60"), 
                     logp = TRUE, #ylab = "srda_p", 
                     xlab = 'CHROM',
                     highlight = int_srda, suggestiveline = FALSE,
                     genomewideline = FALSE, annotatePval = FALSE)


#prda manhattan
prda.env.cand   #132

prda_snps = pvals[pvals$prda_p < 0.0027,]
head(prda_snps)
NROW(prda_snps)
int_prda = prda_snps$ID

prda_man = manhattan(pvals, chr = "CHROM", bp = "POS", p= 'prda_p', snp = "ID",col = c("grey20", "grey60"), 
                        logp = TRUE, #ylab = "prda_p", 
                     xlab = 'CHROM',highlight = int_prda, suggestiveline = FALSE,
                        genomewideline = FALSE, annotatePval = FALSE)


#pcadapt manhattan
outliers_pcadapt #114
pcadapt_snps = pvals[pvals$pcadapt_P < 0.05,]
head(pcadapt_snps)
NROW(pcadapt_snps)
int_pcadapt = pcadapt_snps$ID
pcadapt_man = manhattan(pvals, chr = "CHROM", bp = "POS", p= 'pcadapt_p', snp = "ID",col = c("grey20", "grey60"), 
                     logp = TRUE, #ylab = "pcadapt_P", 
                     xlab = 'CHROM',highlight = int_pcadapt, suggestiveline = FALSE,
                     genomewideline = FALSE, annotatePval = FALSE)


####Variation proportion plots####
# library
library(ggplot2)

var_data = read.csv("variance_datatable.csv")
var_data

# Stacked + percent
var.plot = ggplot(var_data, aes(fill=Factor, y=Value, x=Variation)) + 
  geom_bar(position="fill", stat="identity") + 
  scale_fill_brewer(palette = 'Set1')
var.plot 

ggplot(data=var_data, aes(x=Variation, y=Value, fill=Factor)) +
  geom_bar(stat="identity", position=position_dodge())+
  geom_text(aes(label=Value), vjust=1.6, color="white",
            position = position_dodge(0.9), size=3.5)+
  scale_fill_brewer(palette="Paired")+
  ggtitle("Plot of explainable SNP variation and Structure") +
  xlab("Variation") + ylab("Explained variation (%)")+
  theme_minimal()


###outlier genes associated to abiotic stress####
#import genes associated with abiotic stress
abiogenes = read.csv("OutlierGenes_abiotic stress.csv")
abiogenes
dim(abiogenes)

#filter the genomic data dataframe with the names of the genes associated with abiotic stress

epafro_abgen = epafro_gen_df[, which((names(epafro_gen_df) %in% abiogenes$Name) == TRUE) ]
head(epafro_abgen)
dim(epafro_abgen)
epafro_abgen$taxon = rownames(epafro_abgen)# create a column for the taxon
epafro_abgen = epafro_abgen[,c(65, 1:64)]
mkrlist = colnames(epafro_abgen[,-1])



pop.info = read.csv('epafro_envqced.csv') #import the population information
pop.info = pop.info[, c(1:7)]

epafro_abgeninfo = merge(pop.info, epafro_abgen, by= "taxon")
epafro_abgeninfo
dim(epafro_abgeninfo)
head(epafro_abgeninfo)

write.csv(epafro_abgeninfo, "Species_outliergenes.csv", row.names = FALSE)


####outlier genes and functions for every method####
#outlier names for every method
lfmm_outliers   # 51   
epafro.rda.cand #130
prda.env.cand   #132
outliers_pcadapt #114

#all genes
outliergenes = read.csv("Alloutlier_genes.csv")

#all genes and outlier method scores
lfmm_genes= subset(outliergenes, Name %in% c(lfmm_outliers)) #40
sRDA_genes = subset(outliergenes, Name %in% c(epafro.rda.cand)) #109
pRDA_genes = subset(outliergenes, Name %in% c(prda.env.cand)) #124
pcadapt_genes = subset(outliergenes, Name %in% c(outliers_pcadapt)) #91


# genes clearly associated with abiotic stress
abiogenes = read.csv("OutlierGenes_abiotic stress.csv")

#gene clearly associated with abiotic stress

lfmm_abgenes= subset(abiogenes, Name %in% c(lfmm_outliers)) #7
sRDA_abgenes = subset(abiogenes, Name %in% c(epafro.rda.cand)) #21
pRDA_abgenes = subset(abiogenes, Name %in% c(prda.env.cand)) #22
pcadapt_abgenes = subset(abiogenes, Name %in% c(outliers_pcadapt)) #20

#simple grouped plot method and gene frequencies 
method__df = as.matrix(data.frame(LFMM= c(40,7), sRDA= c(109,21), pRDA = c(142,22), PCAdapt = c(91,20)))
row.names(method__df) = c("All genes", "Abiotic stress genes")
method__df

barplot(method__df,
        col = c("black", "grey"),
        main = "# genes for every method outlier method",
        xlab = "Method",
        ylab = "Frequency",
        beside = TRUE)
legend("topright",
       legend = c("All genes", "Abiotic stress genes"),
       fill = c("black", "grey"))


#grouped ggplot method and gene frequencies 
method = c("LFMM", "LFMM", "sRDA", "sRDA", "pRDA", "pRDA", "PCAdapt", "PCAdapt")
group = c("All genes", "Abiotic stress genes","All genes", "Abiotic stress genes","All genes", "Abiotic stress genes",
          "All genes", "Abiotic stress genes")

#Frequency = c(40,7,109,21, 142,22,91,20)
proportion = c(10.9, 10.9, 29.9, 32.81, 38.8, 34.4, 24.9, 31.3)

df = data.frame(method, group, proportion)
df %>%
  ggplot(aes(method,proportion, fill=group))+
  geom_col(position="dodge") +
  labs(title="Proportion of genes detected by detection method",
       x="Method", y= "Proportion (%)") +
  geom_text(aes(label = round(proportion, 1)), 
            position = position_dodge(0.9),
            color="black",vjust = 0.5,hjust = 1)


####gene frequency plots for species####
setwd("D:/LANDSCAPE_GENOMICS/EggplantCWR/EggplantAfro")
all_mat = read.csv("allspecies_abiotic_genes.csv")
head(all_mat)
#c("aculeastrum", "aculetissimum", "aethiopiccum", "anguivi", "anomalum", "arundo", "campylacanthum", "cerasiferum", "coagulans", "dasyanthum", "dasyphyllum", "incanum", "macrocarpon", "mauense", "nigriviolaceum", "phoxocarpon", "setaceum", "sp")

#"X1065940.181.." gene
all_sub = all_mat[,c("species", "lon", "lat", "country","X1065940.181..")]
head(all_sub)
all_sub2 = all_sub[all_sub$species %in% c("aculeatissimum", "aethiopicum", "anguivi", "anomalum", "campylacanthum", 
                                        "cerasiferum", "coagulans", "dasyanthum", "dasyphyllum", "incanum", 
                                        "macrocarpon", "mauense", "setaceum"),]
head(all_sub2)
coord = all_sub2[,c(2:3)]
sample_sites = all_sub2[,4]
library(sp)
library(tidyr)
library(ggplot2)
coordinates(all_sub) = c("lon", "lat")
plot(all_sub)

mj = subset(all_sub2, X1065940.181.. == "MJ") 
mn = subset(all_sub2, X1065940.181.. == "MN")


a = table(mj$species)
b = table(mn$species)
a
b
ac = table(mj$country)
ac
bc = table(mn$country)
bc

spec = c("aculeatissimum", "aculeatissimum","aethiopicum","aethiopicum", "anguivi","anguivi", "anomalum", "anomalum",
         "campylacanthum","campylacanthum", "cerasiferum", "cerasiferum", "coagulans", "coagulans","dasyanthum",
         "dasyanthum", "dasyphyllum","dasyphyllum", "incanum", "incanum",  "macrocarpon", "macrocarpon", 
         "mauense","mauense", "setaceum","setaceum")
allele = c("mja", "mna","mja", "mna","mja", "mna","mja", "mna","mja", "mna","mja", "mna","mja", "mna",
          "mja", "mna","mja", "mna","mja", "mna","mja", "mna","mja", "mna","mja", "mna")

#Allele_freq = c(2,0,0,4,0,13,0,22,11,1,13,11,1,3,0,1,0,12,17,1,0,15,0,2,1,1) #RH38:1501545:284:+
Allele_freq = c(0,0, 0,4, 1,13, 1,21, 0,13, 1,22, 0,5, 2,0, 11,1, 0,18, 15,0, 0,2, 0,2)


allele_df = data.frame(spec,allele, Allele_freq)
allele_df


allele_df %>%
  ggplot(aes(spec,Allele_freq, fill=allele))+
  geom_col(position="dodge") +
  labs(title="Species Vs Allele frequency (protein FIP1:1065940:181:-)",
       x="Species", y= "Allele Frequency") +
  coord_flip()

###plot stress genes for populations
##major alleles
stg_df = all_mat[,c("species", "lon", "lat", "country","X93733.156..",'X224364.27..', 'X575767.196..','X695411.6..',
                    'X711581.46..','X1065940.15..','X1206322.139..', 'X1249948.100..','X1276278.21..','X1501545.284..'
                    ,'X1510703.8..', 'X1065940.181..')]
head(stg_df)

mja = subset(stg_df, X1065940.181.. == "MJ") 
mna = subset(stg_df, X1065940.181.. == "MN")


ac = table(mja$country)
ac
bc = table(mna$country)
bc

cntry = c("GHA","GHA", "NGA","NGA","UGA","UGA","SDN","SDN", "KEN","KEN","TZA","TZA")
allele = c("mja", "mna","mja", "mna","mja", "mna","mja", "mna","mja", "mna","mja", "mna")

al_freq = c(2,3, 26,32, 2,18, 0,43, 1,10, 0,2)


allele_cntry_df = data.frame(cntry,allele, al_freq)
allele_cntry_df

colr = c ("Blue", "Black")
allele_cntry_df %>%
  ggplot(aes(x = cntry, y=al_freq, fill = allele))+
  geom_col(position="dodge") +
  labs(title="Population Vs Allele frequency (protein FIP1:1065940:181:-)",
       x="Population", y= "Allele Frequency") +
  coord_flip()


##allele 1501545:284:+ plot (RH38 protein) #### 1065940:181:-(FIP1 protein)
str_matrix = read.csv("D:/LANDSCAPE_GENOMICS/EggplantCWR/EggplantAfro/K_structure.csv")

al = all_mat$X1065940.181..
al_df = all_mat[,c(1:7,76)]
dim(all_mat)
  
al_df =   cbind (al, str_matrix, coord, sample_sites)
al_df
al.admix = al_df[,c(2:9)]
al.coord = al_df[,c(5,6)]

library("rnaturalearth")
library("rnaturalearthdata")

world <- ne_countries(scale = "medium", returnclass = "sf")
class(world)
world$admin

ggplot(data = world) +
  geom_sf() +
  geom_point(data = all_mat[all_mat$X1065940.181..== "MJ", ], aes(x = lon, y = lat), size = 2, 
             shape = 21, fill = "red") +
  geom_point(data = all_mat[all_mat$X1065940.181..== "MN", ], aes(x = lon, y = lat), size = 2, 
             shape = 21, fill = "blue") +
  coord_sf(xlim = c(-10, 50), ylim = c(-10, 20), expand = FALSE) 
  theme(legend.position ="bottomleft")

####alleles on koppen climate####
# required packages 
library(raster); 
library(rasterVis); 
library(rworldxtra); 
data(countriesHigh)
library(latticeExtra)

# Read raster files
period='1986-2010'
r <- raster(paste('KG_', period, '.grd', sep=''))

# Color palette for climate classification
climate.colors=c("#960000", "#FF0000", "#FF6E6E", "#FFCCCC", "#CC8D14", "#CCAA54", "#FFCC00", "#FFFF64", "#007800", "#005000", "#003200", "#96FF00", "#00D700", "#00AA00", "#BEBE00", "#8C8C00", "#5A5A00", "#550055", "#820082", "#C800C8", "#FF6EFF", "#646464", "#8C8C8C", "#BEBEBE", "#E6E6E6", "#6E28B4", "#B464FA", "#C89BFA", "#C8C8FF", "#6496FF", "#64FFFF", "#F5FFFF")

# Legend must correspond to all climate classes, insert placeholders
r0 <- r[1:32]; r[1:32] <- seq(1,32,1)

# Converts raster field to categorical data
r <- ratify(r); rat <- levels(r)[[1]]

# Legend is always drawn in alphabetic order
rat$climate <- c('Af', 'Am', 'As', 'Aw', 'BSh', 'BSk', 'BWh', 'BWk', 'Cfa', 'Cfb','Cfc', 'Csa', 'Csb', 'Csc', 'Cwa','Cwb', 'Cwc', 'Dfa', 'Dfb', 'Dfc','Dfd', 'Dsa', 'Dsb', 'Dsc', 'Dsd','Dwa', 'Dwb', 'Dwc', 'Dwd', 'EF','ET', 'Ocean')

# Remove the placeholders
r[1:32] <- r0; levels(r) <- rat

# Select region (Australia)
# x1=80; x2=180; y1=-50; y2=20; xat=5; yat=5	
# Select region (Europe)
# x1=-20; x2=80; y1=30; y2=75; xat=5; yat=5		
# Select region (US)
# x1=-130; x2=-60; y1=20; y2=60; xat=5; yat=5
# Select region (Global)
# x1=-180; x2=180; y1=-90; y2=90; xat=20; yat=10

#select region
x1=-10; x2=50; y1=-10; y2=20
x = list(x = c(-10, 50), y = c(-10,20))
africa = extent(x)

####Select region 
#africa = readOGR("Africa.shp")
r_africa <- crop(r, africa)
plot(r_africa)
#loci data frame
al_df #see allele plots above

#allele dataframes containing the coordinates
mj.allele.df = all_mat[all_mat$X1065940.181..== "MJ", ]
mn.allele.df = all_mat[all_mat$X1065940.181..== "MN", ]

# Visualization		
if(.Platform$OS.type=="windows") {quartz<-function(...) windows(...)}
quartz(width=13, height=10)

allele_kop = rasterVis::levelplot(r_africa, col.regions=climate.colors, xlab="", ylab="", 
                                  scales=list(x=list(limits=c(xmin(r_africa), xmax(r_africa)), at=seq(xmin(r_africa), xmax(r_africa))), 
                                              y=list(limits=c(ymin(r_africa), ymax(r_africa)), at=seq(ymin(r_africa), ymax(r_africa))))) + 
latticeExtra::layer(panel.points(mn.allele.df[c("lon", "lat")], cex = 1, pch = 20, col = "blue")) +
latticeExtra::layer(panel.points(mj.allele.df[c("lon", "lat")], cex = 1, pch = 20, col = "black")) 

plot(allele_kop)

#####Adaptive Index####

study_area = readOGR("Africa.shp")
extent = c(xlim = c(-10, 50), ylim = c(-10, 20))
study_area  = raster::crop(study_area, extent)
plot(study_area)

#
Y1 = read.csv("outlier_snps.csv")
View(Y1)
row.names(Y1) = Y1$Accessions
Y1 = Y1[, -1]

outliers_df = read.csv("all_outliers.csv") #all methods outliers
gea_outliers

outlr_gea = read.csv("epafro_outlier_gea.csv") #only gea outliers
View(outlr_gea)

Y2 = outlr_gea
row.names(Y2) = Y2$Taxa
Y2 = Y2[,-1]
View(Y2)

#
head(epafroenv_info_qced) 
X = epafroenv_info_qced[,-c(1:7)]
head(X)

# Standardization of the environmental variables
X_scaled <- scale(X, center=TRUE, scale=TRUE) 
View(X_scaled)

# Now, we can create a table gathering all the predictors (environment + population structure + space)
predictors_rda <- data.frame(X_scaled, str_matrix, sel.mem)
View(predictors_rda)
# 2. Once the data set of candidate SNPs is created, the second step is to run the
# "adaptively enriched" RDA, this time without conditioning for population structure:

RDA_outliers1 <- rda(Y1 ~ PDrM_14 + PWaQ_18 + PCoQ_19 +  MTWeQ_8 + srad + phh2o + nitrogen
                    + clay + silt + ocd, predictors_rda)

RDA_outliers2 <- rda(Y2 ~ PDrM_14 + PWaQ_18 + PCoQ_19 +  MTWeQ_8 + srad + phh2o + nitrogen
                    + clay + silt + ocd, predictors_rda)
RDA_outliers1
RDA_outliers2

summary(RDA_outliers1)
summary(RDA_outliers2)
# Let's have a look at the RDA biplot ("adaptively enriched RDA space") to have an idea 
# about the association between markers and environment:

plot(RDA_outliers1, type="n", scaling=3, xlim=c(-1,1), ylim=c(-1,1), xlab = "RDA1 (69.2%)", ylab = "RDA2 (19.4%)", main = "Adaptively enriched RDA space")
points(RDA_outliers1, display="species", pch=16, cex=1, col="orange", scaling=3)
text(RDA_outliers1, scaling=3, display="bp", col="#0868ac", cex=1)

plot(RDA_outliers2, type="n", scaling=3, xlim=c(-1,1), ylim=c(-1,1), xlab = "RDA1 (17.2%)", ylab = "RDA2 (4.8%)", main = "Adaptively enriched RDA space")
points(RDA_outliers2, display="species", pch=16, cex=1, col="orange", scaling=3)
text(RDA_outliers2, scaling=3, display="bp", col="#0868ac", cex=1)

# Now, we need to identify the significant RDA axes to be used in the calculation 
# of the adaptive index:
signif.axis1 <- anova.cca(RDA_outliers1, by = "axis", parallel = getOption("mc.cores"))
signif.axis2 <- anova.cca(RDA_outliers2, by = "axis", parallel = getOption("mc.cores"))
signif.axis1
signif.axis2

# 3. Now, we can calculate the adaptive index across the landscape: a prediction of 
# adaptive similarity across the landscape based on the relative importance of the focal 
# drivers of divergent selection in shaping adaptive variation. 

# Raster stack where to store the environmental variables used:
#climate
library(raster)
library(ggplot2)
sel_bio = list.files(path = "D:/LANDSCAPE_GENOMICS/EggplantCWR/EggplantAfro/Sel_bioclim", pattern = ".tif",
                     all.files = TRUE, full.names = FALSE)
sel_bio <- stack(paste0(path = "D:/LANDSCAPE_GENOMICS/EggplantCWR/EggplantAfro/Sel_bioclim", 
                              "/", sel_bio))

srad_files= list.files(path = "D:/LANDSCAPE_GENOMICS/EggplantCWR/EggplantAfro/srad", pattern = ".tif",
                        all.files = TRUE, full.names = FALSE)
srad_stack = stack(paste0(path = "D:/LANDSCAPE_GENOMICS/EggplantCWR/EggplantAfro/srad", 
                          "/", srad_files))

srad = mean (srad_stack)

#soil
#phh2o_15_30cm
ph = raster("D:/LANDSCAPE_GENOMICS/Gap analysis/soil data/phh2o_15-30cm_SoilGrids.tif")
#nitrogen_15_30cm
nitro = raster("D:/LANDSCAPE_GENOMICS/Gap analysis/soil data/nitrogen_15-30cm_SoilGrids.tif")
#clay
clay = raster("D:/LANDSCAPE_GENOMICS/Gap analysis/soil data/clay_15-30cm_SoilGrids.tif")
#silt
silt = raster("D:/LANDSCAPE_GENOMICS/Gap analysis/soil data/silt_15-30cm_SoilGrids.tif")
#ocd
ocd = raster("D:/LANDSCAPE_GENOMICS/Gap analysis/soil data/ocd_15-30cm_SoilGrids.tif")

#stack soil files

soil_files = c(ph, nitro, clay, silt, ocd)
soil_rasters =  raster::stack(soil_files)
soil_rasters

#stack soil and bio rasters to make a predictors raster

soil_rasters2 = raster::resample(soil_rasters, sel_bio) # change soil rasters to same resolution as bio
res(soil_rasters2)
res(sel_bio)

sel_bio2 = raster::crop(sel_bio, extent(soil_rasters2)) # make wclm and soil raster to same extent
srad2 = raster::crop(srad, extent(soil_rasters2))

#combine all the env data
stack_current <- raster::stack(sel_bio2, srad2, soil_rasters2)
names(stack_current) = colnames(X)
# As we are going to multiply the scores of the environmental variables by  
# their standardized values, we need to retrieve the scaling coefficients of the 
# standardized environmental matrix we used in RDA previously (Practical 4):
scale_X <- attr(X_scaled, 'scaled:scale')
center_X <- attr(X_scaled, 'scaled:center')

# adaptive index calculation based on the 'adaptive_index()' function (credits: 
# Capblancq & Forester 2021:
adaptive_landscape_current <- adaptive_index(
  RDA = RDA_outliers, # adaptively enriched RDA
  K = 2, # n.er of significant RDA axes 
  env_pres = stack_current, # the  
  range = study_area, 
  method = "loadings", # here we specify to use the loadings of the env. variables 
  # on the RDA axes in the calculation of the adaptive index 
  scale_env = scale_X, # scaling coefficients 
  center_env = center_X)


RDA1_Adapt_index = adaptive_landscape_current$RDA1
RDA2_Adapt_index = adaptive_landscape_current$RDA2
adap_ind_rda1 = raster::extract(RDA1_Adapt_index, coord)
adap_ind_rda2 = raster::extract(RDA2_Adapt_index, coord)
adap_ind = cbind(adap_ind_rda1, adap_ind_rda2)


specie_adap_indexdf = cbind(epafro_env_info, adap_ind)
write.csv(specie_adap_indexdf, "epafro_adaptive index.csv", row.names = FALSE)

colpal=colorRampPalette(c("blue", "yellow"))
p1 <- ggplot() + 
  geom_raster(aes(x=raster_rda1[,1],y=raster_rda1[,2],fill=cut(raster_rda1[,3], breaks=c(-6,-3,0,3,6,9,12,15)))) + 
  scale_fill_manual(values=colpal(8)[-1]) +
  geom_point(aes(x= cities@coords[,1], y=cities@coords[,2]), size=1) +
  geom_text(aes(x= cities@coords[,1]+0.05, y=cities@coords[,2]+0.05, label=cities$CITY_NAME), size=3.5) +
  xlab("Longitude") + ylab("Latitude") +
  guides(fill=guide_legend(title="RDA scores")) +
  theme_bw() +
  theme(panel.grid.major = element_blank())

# Let's plot the adaptive landscape based on loadings from RDA1
plot(adaptive_landscape_current$RDA1, main = "Adaptive landscape (RDA1)")

# ... and RDA2
plot(adaptive_landscape_current$RDA2, main = "Adaptive landscape (RDA1)")

# As an exercise, let's finally try to outline evolutionary significant units (ESUs)
# together with adaptive groups among the studied populations for conservation 
# purposes (Funk et al. 2012):
adaptive_index_RDA1 <- extract(adaptive_landscape_current$RDA1, coordinates)
adaptive_index_RDA2 <- extract(adaptive_landscape_current$RDA2, coordinates)

jpeg("practical6 - adaptive index RDA1.jpeg", width = 7, height = 7, units = 'in', res = 800)
plot(study_area, ylim = c(-1.8, 4.3), xlab = "Longitude", ylab = "Latitude", 
     cex.lab = 1.5, col="lightgray", main = "Adaptive index based on RDA1", 
     cex.main = 1.5, font.main = 1)
plot(grid, add = TRUE)

# here we add individuals with symbol based on 'find.cluster()' analysis and colors 
# based on the adaptive index (RDA1)
plot(coordinates, add = TRUE, cex = 2, 
     pch = coordinates@data$symbol, 
     bg = dapc_color(4)[as.numeric(cut(adaptive_index_RDA1, breaks = 4))]) 

# add a legend 
legend("topleft", pch = c(21, 24), legend = c("ESU 1", "ESU 2"), pt.bg = "black", 
       pt.cex=1.3, cex = 1.3)

axis(1); axis(2)
addscalebar(plotepsg = "4326", plotunit = "Km", 
            pos = "bottomright", padin = c(0.3, 0.3), style = "bar",
            label.cex = 1.7, widthhint = 0.2, htin = 0.02,
            label.col = "black", bar.cols = "black", linecol = "black")
addnortharrow(pos = "topleft", scale = 1.05, padin = c(0.5, 0.9))
box()
dev.off()

# The same but with the adaptive index based on the RDA2 loadings:
jpeg("practical6 - adaptive index RDA2.jpeg", width = 7, height = 7, units = 'in', res = 800)
plot(study_area, ylim = c(-1.8, 4.3), xlab = "Longitude", ylab = "Latitude", 
     cex.lab = 1.5, col="lightgray", main = "Adaptive index based on RDA2", 
     cex.main = 1.5, font.main = 1)
plot(grid, add = TRUE)

# here we add individuals with symbol based on 'find.cluster()' analysis and colors 
# based on the adaptive index (RDA2)
plot(coordinates, add = TRUE, cex = 2, 
     pch = coordinates@data$symbol, 
     bg = dapc_color(4)[as.numeric(cut(adaptive_index_RDA2, breaks = 4))]) 
legend("topleft", pch = c(21, 24), legend = c("ESU 1", "ESU 2"), pt.bg = "black", 
       pt.cex=1.3, cex = 1.3)
axis(1); axis(2)
addscalebar(plotepsg = "4326", plotunit = "Km", 
            pos = "bottomright", padin = c(0.3, 0.3), style = "bar",
            label.cex = 1.7, widthhint = 0.2, htin = 0.02,
            label.col = "black", bar.cols = "black", linecol = "black")
addnortharrow(pos = "topleft", scale = 1.05, padin = c(0.5, 0.9))
box()
dev.off()
