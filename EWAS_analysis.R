#Setting up parallel processors
library(doParallel)
library(lme4)
cl<-makeCluster(8)
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
pheno$Brain_ID <-as.character(pheno$Brain_ID)
pheno$Basename <-as.character(pheno$Basename)

pheno <-pheno[-grep("202093120076_R05C01|202093120076_R08C01|203734300117_R02C01", rownames(pheno)),]
betas <- betas[,-grep("202093120076_R05C01|202093120076_R08C01|203734300117_R02C01", colnames(betas))]

# 
# ##Mixed effect model for sample type ####
# 
# res<-matrix(data = NA, ncol = 3, nrow = nrow(betas))
# rownames(res)<-rownames(betas)
# colnames(res)<-c("SampleType_Beta", "SampleType_SE", "SampleType_P")
# for(i in 1:nrow(betas)){
# 
#   modellmer<-lmer(betas[i,] ~ pheno$SampleType + pheno$Sex + pheno$Age + as.factor(pheno$Chip) + (1|pheno$ID), REML = FALSE, subset = which(is.na(pheno$SampleType) == FALSE))
#   res[i,1]<-fixef(modellmer)["pheno$SampleTypef"]
#   res[i,2]<-summary(modellmer)$coefficients["pheno$SampleTypef",2]
#   model.null<-lmer(betas[i,] ~  pheno$Sex + pheno$Age + as.factor(pheno$Chip) + (1|pheno$ID), REML = FALSE, subset = which(is.na(pheno$SampleType) == FALSE))
#   res[i,3]<-anova(modellmer,model.null)["modellmer","Pr(>Chisq)"]
# }
# write.csv(res, file ="EWAS_Brain_Regions.csv")
# 
# # As the beta values here represent proportion of DNA methylation (i.e. they lie between 0 and 1), the regression coefficients represent the change in proportion. 
# # Typically we report our findings on the % scale therefore we will multiple the regression coefficients and SE by 100. 
# # This will need to be done for all variables for which you saved the results.
# 
# 
# res[,"SampleType_Beta"]<-res[,"SampleType_Beta"]*100
# res[,"SampleType_SE"]<-res[,"SampleType_SE"]*100
# 
# ### Annotate the output
# 
# # At this stage we will be interested in adding information regarding where each probe is located, and what genes or regulatory features it overlaps with. 
# # There is a lot of annotation available, the code below only takes a subset of columns which I think may be most relevant. 
# 
# 
# epicManifest<-read.csv("/mnt/data1/IndiaWorkshop/epicManifest.csv", stringsAsFactors = F, header=T)
# 
# 
# epicManifest<-epicManifest[match(rownames(res), epicManifest$Name),]
# res<-cbind(res, epicManifest[,c("CHR", "MAPINFO", "UCSC_RefGene_Name", "UCSC_RefGene_Group", "Relation_to_UCSC_CpG_Island", "Regulatory_Feature_Name", "Regulatory_Feature_Group", "DNase_Hypersensitivity_NAME", "OpenChromatin_NAME", "TFBS_Evidence_Count", "Methyl450_Loci", "SNP_ID", "SNP_DISTANCE", "SNP_MinorAlleleFrequency")]) ## this can be edited to include additional columns or exclude as required
# 
# 
# ## QQ plot 
# pdf("/mnt/data1/NIMHANS/Down_stream_analysis/NIMHAMSg/EWAS/Plots/QQplot_of_EWAS_SampleType.pdf")
# qq(res$SampleType_P)
# dev.off()
# 
# 
# ##Manhattan Plot
# res<-res[which(res$CHR != "Y"),] ## to also exclude the X chromosome repeat and edit this line
# 
# res$CHR<-as.character(res$CHR)
# res$CHR[which(res$CHR == "X")]<-23 ## to also recode Y chromosome repeat and edit this line
# res$CHR<-as.numeric(res$CHR)
# res<-res[which(res$CHR != ""),]
# 
# bonfP<-0.05/nrow(res)
# 
# pdf("/mnt/data1/NIMHANS/Down_stream_analysis/NIMHAMSg/EWAS/Plots/Manhattan_plot_Mixed_Effect_Model_SampleType.pdf")
# manhattan(res, p = "Status_P", bp = "MAPINFO", chr = "CHR", genomewide = -log10(bonfP), suggestiveline = -log10(5e-5), logp=T, col=c("blue","yellow"), ylim=c(0,90))
# dev.off()
# 
# ### Mixed effect model of Sex #####
# 
# res1<-matrix(data = NA, ncol = 3, nrow = nrow(betas))
# rownames(res1)<-rownames(betas)
# colnames(res1)<-c("Sex_Beta", "Sex_SE", "Sex_P")
#   
# for(i in 1:nrow(betas)){
#   modellmer1<-lmer(betas[i,] ~ pheno$SampleType + pheno$Sex + pheno$Age + as.factor(pheno$Chip) + (1|pheno$ID), REML = FALSE, subset = which(is.na(pheno$Sex) == FALSE))
#   res1[i,1]<-fixef(modellmer1)["pheno$SexM"]
#   res1[i,2]<-summary(modellmer1)$coefficients["pheno$SexM",2]
#   model.null1<-lmer(betas[i,] ~  pheno$SampleType + pheno$Age + as.factor(pheno$Chip) + (1|pheno$ID), REML = FALSE, subset = which(is.na(pheno$Sex) == FALSE))
#   res1[i,3]<-anova(modellmer1,model.null1)["modellmer1","Pr(>Chisq)"]
# }
# 
# write.csv(res1, file="/mnt/data1/NIMHANS/Down_stream_analysis/NIMHAMSg/EWAS/EWAS_Sex_analysis.csv")
# 
# # As the beta values here represent proportion of DNA methylation (i.e. they lie between 0 and 1), the regression coefficients represent the change in proportion. 
# # Typically we report our findings on the % scale therefore we will multiple the regression coefficients and SE by 100. 
# # This will need to be done for all variables for which you saved the results.
# 
# 
# res[,"Sex_Beta"]<-res[,"Sex_Beta"]*100
# res[,"Sex_SE"]<-res[,"Sex_SE"]*100
# 
# ### Annotate the output
# 
# # At this stage we will be interested in adding information regarding where each probe is located, and what genes or regulatory features it overlaps with. 
# # There is a lot of annotation available, the code below only takes a subset of columns which I think may be most relevant. 
# 
# 
# #epicManifest<-read.csv("/mnt/data1/IndiaWorkshop/epicManifest.csv", stringsAsFactors = F, header=T)
# 
# 
# epicManifest<-epicManifest[match(rownames(res1), epicManifest$Name),]
# res1<-cbind(res1, epicManifest[,c("CHR", "MAPINFO", "UCSC_RefGene_Name", "UCSC_RefGene_Group", "Relation_to_UCSC_CpG_Island", "Regulatory_Feature_Name", "Regulatory_Feature_Group", "DNase_Hypersensitivity_NAME", "OpenChromatin_NAME", "TFBS_Evidence_Count", "Methyl450_Loci", "SNP_ID", "SNP_DISTANCE", "SNP_MinorAlleleFrequency")]) ## this can be edited to include additional columns or exclude as required
# 
# 
# ## QQ plot 
# pdf("/mnt/data1/NIMHANS/Down_stream_analysis/NIMHAMSg/EWAS/Plots/QQplot_of_EWAS_Sex.pdf")
# qq(res1$Sex_P)
# dev.off()
# 
# 
# ##Manhattan Plot
# res1<-res1[which(res1$CHR != "Y"),] ## to also exclude the X chromosome repeat and edit this line
# res1<-res1[which(res1$CHR != "X"),] ## to also exclude the X chromosome repeat and edit this line
# 
# 
# # res1$CHR<-as.character(res1$CHR)
# # res1$CHR[which(res1$CHR == "X")]<-23 ## to also recode Y chromosome repeat and edit this line
# res1$CHR<-as.numeric(res1$CHR)
# res1<-res1[which(res1$CHR != ""),]
# 
# 
# pdf("/mnt/data1/NIMHANS/Down_stream_analysis/NIMHAMSg/EWAS/Plots/Manhattan_plot_of_Sex.pdf")
# manhattan(res1, p = "Status_P", bp = "MAPINFO", chr = "CHR", genomewide = -log10(bonfP), suggestiveline = -log10(5e-5), logp=T, col=c("blue","yellow"), ylim=c(0,50))
# dev.off()

### Mixed effect model of Age #####

testCpG<-function(row, pheno){

res<-matrix(data = NA, ncol = 3, nrow = nrow(betas))
rownames(res)<-rownames(betas)
colnames(res)<-c("Age_Beta", "Age_SE", "Age_P")

for(i in 1:nrow(betas)){
  modellmer<-lmer(betas[i,] ~ pheno$Age + pheno$Gender + pheno$Plate + pheno$prop + (1|pheno$BR), REML = FALSE, subset = which(is.na(pheno$Age) == FALSE))
  res[i,1]<-fixef(modellmer)["pheno$Age"]
  res[i,2]<-summary(modellmer)$coefficients["pheno$Age", 2]
  model.null<-lmer(betas[i,] ~  pheno$Gender + pheno$Plate + pheno$prop + (1|pheno$BR), REML = FALSE, subset = which(is.na(pheno$Age) == FALSE))
  res[i,3]<-anova(modellmer,model.null)["modellmer","Pr(>Chisq)"]
}
}

res<-foreach(i=1:nrow(betas), .combine=rbind) %dopar%{
  testCpG(betas[i,], pheno)	
}



testCpG<-function(row, pheno){
  
    modellmer<-lmer(betas[1,] ~ Age + Gender + Plate + prop + (1|BR), data = pheno, REML = FALSE, subset = which(is.na(pheno$Age) == FALSE))
    Beta<-fixef(modellmer)["Age"]
    SE<-summary(modellmer)$coefficients["Age", 2]
    model.null<-lmer(betas[i,] ~ Gender + Plate + prop + (1|BR), data= pheno, REML = FALSE, subset = which(is.na(pheno$Age) == FALSE))
    P<-anova(modellmer,model.null)["modellmer","Pr(>Chisq)"]
    return(c(Beta,SE,P))
}


res<-foreach(i=1:10, .combine=rbind) %dopar%{
  testCpG(betas[i,], pheno)	
}

rownames(res)<-rownames(betas)
colnames(res)<-c("Age_Beta", "Age_SE", "Age_P")

res[,"Age_Beta"]<-res[,"Age_Beta"]*100
res[,"Age_SE"]<-res[,"Age_SE"]*100

### Annotate the output

# At this stage we will be interested in adding information regarding where each probe is located, and what genes or regulatory features it overlaps with. 
# There is a lot of annotation available, the code below only takes a subset of columns which I think may be most relevant. 


epicManifest<-read.csv("/gpfs/mrc0/projects/Research_Project-MRC193462/EWAS/epicManifest_hg38.csv", stringsAsFactors = F, header=T, sep=" ")


epicManifest<-epicManifest[match(rownames(res), epicManifest$Name),]
res<-cbind(res, epicManifest[,c("CHR", "MAPINFO", "UCSC_RefGene_Name", "UCSC_RefGene_Group", "Relation_to_UCSC_CpG_Island", "Regulatory_Feature_Name", "Regulatory_Feature_Group")]) ## this can be edited to include additional columns or exclude as required


write.csv(res, file="/gpfs/mrc0/projects/Research_Project-MRC193462/EWAS/EWAS_Age_analysis.csv")


# res<-matrix(data = NA, ncol = 3, nrow = nrow(betas))
# rownames(res)<-rownames(betas)
# colnames(res)<-c("Age_Beta", "Age_SE", "Age_P")
# 
# for(i in 1:nrow(betas)){
#   modellmer<-lmer(betas[i,] ~ pheno$Age + (1|pheno$Brain_ID) + pheno$Gender + pheno$Plate + pheno$Chip + pheno$BR + pheno$prop, REML = FALSE, subset = which(is.na(pheno$Age) == FALSE))
#   res[i,1]<-fixef(modellmer)["pheno$Age"]
#   res[i,2]<-summary(modellmer)$coefficients["pheno$Age", 2]
#   model.null<-lmer(betas[i,] ~  (1|pheno$Brain_ID) + pheno$Gender + pheno$Plate + pheno$Chip + pheno$BR + pheno$prop, REML = FALSE, subset = which(is.na(pheno$Age) == FALSE))
#   res[i,3]<-anova(modellmer,model.null)["modellmer","Pr(>Chisq)"]
# }
# 
# 
# # As the beta values here represent proportion of DNA methylation (i.e. they lie between 0 and 1), the regression coefficients represent the change in proportion. 
# # Typically we report our findings on the % scale therefore we will multiple the regression coefficients and SE by 100. 
# # This will need to be done for all variables for which you saved the results.
# 
# 
# res[,"Age_Beta"]<-res[,"Age_Beta"]*100
# res[,"Age_SE"]<-res[,"Age_SE"]*100

### Annotate the output

# At this stage we will be interested in adding information regarding where each probe is located, and what genes or regulatory features it overlaps with. 
# There is a lot of annotation available, the code below only takes a subset of columns which I think may be most relevant. 




# epicManifest<-epicManifest[match(rownames(res), epicManifest$Name),]
# res<-cbind(res, epicManifest[,c("CHR", "MAPINFO", "UCSC_RefGene_Name", "UCSC_RefGene_Group", "Relation_to_UCSC_CpG_Island", "Regulatory_Feature_Name", "Regulatory_Feature_Group", "DNase_Hypersensitivity_NAME", "OpenChromatin_NAME", "TFBS_Evidence_Count", "Methyl450_Loci", "SNP_ID", "SNP_DISTANCE", "SNP_MinorAlleleFrequency")]) ## this can be edited to include additional columns or exclude as required
# 
# write.csv(res, file="/gpfs/mrc0/projects/Research_Project-MRC193462/EWAS/EWAS_Age_individual.randomeffect.csv")

# ## QQ plot 
# pdf("/mnt/data1/NIMHANS/Down_stream_analysis/NIMHAMSg/EWAS/Plots/QQplot_of_EWAS_Age.pdf")
# qq(res1$Age_P)
# dev.off()
# 
# 
# ##Manhattan Plot
# res2<-res2[which(res2$CHR != "Y"),] ## to also exclude the X chromosome repeat and edit this line
# 
# res2$CHR<-as.character(res2$CHR)
# res2$CHR[which(res2$CHR == "X")]<-23 ## to also recode Y chromosome repeat and edit this line
# res2$CHR<-as.numeric(res2$CHR)
# res2<-res2[which(res2$CHR != ""),]
# 
# bonfP<-0.05/nrow(res2)
# 
# pdf("/mnt/data1/NIMHANS/Down_stream_analysis/NIMHAMSg/EWAS/Plots/Manhattan_Plots_of_Age.pdf")
# manhattan(res2, p = "Status_P", bp = "MAPINFO", chr = "CHR", genomewide = -log10(bonfP), suggestiveline = -log10(5e-5), logp=T, col=c("blue","yellow"), ylim=c(0,30))
# dev.off()

