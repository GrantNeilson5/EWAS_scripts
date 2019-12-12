
#Setting up parallel processors
library(doParallel)
library(lme4)
cl<-makeCluster(12)
registerDoParallel(cl)
clusterEvalQ(cl, library(lme4))
library(dplyr, warn.conflicts = FALSE)
library(tidyr, warn.conflicts = FALSE)
setwd("/gpfs/mrc0/projects/Research_Project-MRC193462")
load("BDR_combined_1221_QCd.rdat")

rownames(pheno) <- pheno$Basename
pheno <- pheno[order(rownames(pheno)),]
betas <- betas[,order(colnames(betas))]
print(identical(rownames(pheno), colnames(betas)))


# splitting basenmae into chip and postion as Excel often mucks up Chip
pheno$Basename2<-pheno$Basename
pheno<-separate(data = pheno, col = Basename2, into = c("Chip", "Position"), sep="_")


pheno$Chip <- as.factor(pheno$Chip)
pheno$Brain_ID <-as.factor(pheno$Brain_ID)
pheno$Basename <-as.character(pheno$Basename)

pheno <-pheno[-grep("202093120076_R05C01|202093120076_R08C01|203734300117_R02C01", rownames(pheno)),]
betas <- betas[,-grep("202093120076_R05C01|202093120076_R08C01|203734300117_R02C01", colnames(betas))]


testCpG<-function(row, pheno){
  
  modellmer<-lmer(betas[i,] ~ Age + Gender + Plate + prop + (1|BR), data = pheno, REML = FALSE, subset = which(is.na(pheno$Age) == FALSE))
  Beta<-fixef(modellmer)["Age"]
  SE<-summary(modellmer)$coefficients["Age", 2]
  model.null<-lmer(betas[i,] ~ Gender + Plate + prop + (1|BR), data= pheno, REML = FALSE, subset = which(is.na(pheno$Age) == FALSE))
  P<-anova(modellmer,model.null)["modellmer","Pr(>Chisq)"]
  return(c(Beta,SE,P))
}


res<-foreach(i= 1:nrow(betas), .combine=rbind) %dopar%{
  testCpG(betas[i,], pheno)	
}

rownames(res)<-rownames(betas)
colnames(res)<-c("Age_Beta", "Age_SE", "Age_P")

# As the beta values here represent proportion of DNA methylation (i.e. they lie between 0 and 1), the regression coefficients represent the change in proportion. 
# Typically we report our findings on the % scale therefore we will multiple the regression coefficients and SE by 100. 
# This will need to be done for all variables for which you saved the results.


res[,"Age_Beta"]<-res[,"Age_Beta"]*100
res[,"Age_SE"]<-res[,"Age_SE"]*100

### Annotate the output

# At this stage we will be interested in adding information regarding where each probe is located, and what genes or regulatory features it overlaps with. 
# There is a lot of annotation available, the code below only takes a subset of columns which I think may be most relevant. 


epicManifest<-read.csv("/gpfs/mrc0/projects/Research_Project-MRC193462/EWAS/epicManifest_hg38.csv", stringsAsFactors = F, header=T, sep=" ")


epicManifest<-epicManifest[match(rownames(res), epicManifest$Name),]
res<-cbind(res, epicManifest[,c("CHR", "MAPINFO", "UCSC_RefGene_Name", "UCSC_RefGene_Group", "Relation_to_UCSC_CpG_Island", "Regulatory_Feature_Name", "Regulatory_Feature_Group", "DNase_Hypersensitivity_NAME", "OpenChromatin_NAME", "TFBS_Evidence_Count", "Methyl450_Loci", "SNP_ID", "SNP_DISTANCE", "SNP_MinorAlleleFrequency")]) ## this can be edited to include additional columns or exclude as required

write.csv(res, file="/gpfs/mrc0/projects/Research_Project-MRC193462/EWAS/EWAS_Age_BDR.csv")
