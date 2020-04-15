library(qqman)
library(gplots)
library(missMethyl)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(limma)

setwd("/gpfs/mrc0/projects/Research_Project-MRC193462/EWAS/")
load("BDR_combined_1221_QCd.rdat")
res <- read.csv("EWAS_Age_BDR.csv", row.names = 1)

## QQ plot 
pdf("/gpfs/mrc0/projects/Research_Project-MRC193462/EWAS/Plots/QQplot_of_EWAS_Age.pdf")
lamda <- qchisq(1-median(res$phenotype_P),1)/qchisq(0.5,1)
qq(res$phenotype_P, main = "Phenotype Predicted Tissue")
mtext(paste("Lamda = ", signif(lamda,3)), side = 3, adj = 1)
dev.off()


##Manhattan Plot
res<-res[which(res$CHR != "Y"),] ## to also exclude the X chromosome repeat and edit this line
res<-res[which(res$CHR != "X"),] ## to also exclude the X chromosome repeat and edit this line


#res1$CHR<-as.character(res$CHR)
#res$CHR[which(res$CHR == "X")]<-23 ## to also recode Y chromosome repeat and edit this line
#res$CHR[which(res$CHR == "Y")]<-23 ## to also recode Y chromosome repeat and edit this line
res$CHR<-as.numeric(as.character(res$CHR))
res<-res[which(res$CHR != ""),]
res<-res[which(res$MAPINFO != ""),]

res$SNP <- rownames(res)
bonfP<-0.05/nrow(res)

pdf("/gpfs/mrc0/projects/Research_Project-MRC193462/EWAS/Plots/Manhattan_plot_of_Age.pdf")
manhattan(res, p = "Age_P", bp = "MAPINFO", chr = "CHR", genomewide = -log10(bonfP), suggestiveline = -log10(5e-5), logp=T, col=c("blue","yellow"), ylim=c(0,10))
dev.off()

##Heatmaps #######
##heatmap of top 500 significant probes

res<-res[order(res$Age_P),]
res_sig <- res[1:500,]
betas_r <- betas[rownames(res_sig),]

cell_rows <- as.data.frame(pheno$Gender, rownames(pheno))
colnames(cell_rows)<- "Gender"
cell_age <- as.data.frame(pheno$Age, rownames(pheno))
colnames(cell_age)<- "Age"
cell_BR <- as.data.frame(pheno$BR, rownames(pheno))
colnames(cell_BR)<- "BR"
cell_Plate <- as.data.frame(pheno$Plate, rownames(pheno))
colnames(cell_Plate)<- "Plate"
cell_Braak <- as.data.frame(pheno$BraakTangle_numeric, rownames(pheno))
colnames(cell_Braak)<- "Braak"
cell_Braak <- as.data.frame(pheno$BraakTangle_numeric, rownames(pheno))
colnames(cell_Braak)<- "Braak"
cell_Institue <- as.data.frame(pheno$Institute, rownames(pheno))
colnames(cell_Institue)<- "Institue"
cell_rows <- cbind(cell_rows, cell_age, cell_BR, cell_Plate,cell_Braak,cell_Institue)

pdf("/gpfs/mrc0/projects/Research_Project-MRC193462/EWAS/Plots/Heatmap_top_500_probes.pdf")
pheatmap(betas_r, 
         annotation_col = cell_rows,
         cluster_rows = F,
         show_rownames = F,
         show_colnames = F,
         main="Top significant probes between Ages")
dev.off()


####Go Analysis ##### 

ann <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
index<-which(res$Age_P<=0.00001)# this can be set depending on the p-value disribution of your dataset we normally set it to <p=1x10-5
background<-rownames(res)
rownames(res)[index]->test

gst <- gometh(sig.cpg=test, all.cpg=background, collection="GO", plot.bias = FALSE, prior.prob = TRUE)# this is using the GO term
## Warning in alias2SymbolTable(flat$symbol): Multiple symbols ignored for one
## or more aliases
gst<-gst[order(gst$P.DE),]
head(gst)


gst2 <- gst[1:5,] #select top 5 terms
gst2 <- gst2[order(gst2$P.DE, decreasing = T),] #order it so that the plot goes from most to least sig
gst2$P.DE <- signif(gst2$P.DE,3) #shorten the p.vals to three significant figures
gst2$TERM <- factor(gst2$TERM, levels = gst2$TERM) #lock it as factors so it does not go in alphabetical order
pdf("/gpfs/mrc0/projects/Research_Project-MRC193462/EWAS/Plots/GO_Age.pdf")
par(mar =c(15,15,15,30))
ggplot(gst2, aes(x =TERM, y = N)) +
  geom_bar(fill = 'royalblue', stat = 'identity') +
  coord_flip() +
  labs(x = "", y = "", title = 'Top 5 significant GO terms for Age') +
  geom_text(data=gst2,aes(x=TERM,y=N,label=P.DE),position = position_stack(vjust = 0.5)) + #adds the p.val text to GO terms
  theme_classic()
dev.off()
Collapse



