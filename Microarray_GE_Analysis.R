setwd("~/FGT/ICA")

# Affymetrix Microarray Analysis Basic (Skeleton) Workflow 
##Load the required libraries & load the files for the workflow
library(limma)
library(affy)
library(annotate)
library(mouse4302.db)# load chip-specific annotation
library(simpleaffy)
library(ggrepel)
library(scatterplot3d)
library(SummarizedExperiment)
library(org.Hs.eg.db)


## Load the main data- commented code is just for information
##system("tar -xvf GSE23556_RAW.tar")
#####system("cp  /shared_files/targets.txt .")

###### 1. LOADING CELL FILES #######################################################################################
# Load the target file into a table object
adf <-read.table("targets.txt",header=TRUE,as.is=TRUE)
# Load the expression values of all CEL files in the R working directory
raw_data <- ReadAffy()
# rename columns 
sampleNames(raw_data) <- adf$Name

# View a summary of the example data
raw_data

###### 2. QUALITY CONTROL #######################################################################################

### 2.1 Access Control Data (Added code)
control_idx <-grep("AFFX", rownames(raw_data))
control_idx

## 2.2 Build Quality Control Plots
colours <- c(rep("yellow",3),rep("red",3))

#QC Plot (Added code)
x.mas <- call.exprs(raw_data,"mas5")
qcs <-qc(raw_data,x.mas)
png("qc_plot.png")
plot(qcs, xlab="Fold Change", ylab="Array")
dev.off()

#Quality control histogram
png("quality_control_hist.png")
hist(raw_data, main="Quality control of raw data")
dev.off()

# Boxplot with different colour for the two sample group 
png("quality_control_box.png")
boxplot(raw_data, col = adf$Colour, las=2, main="Quality control of raw data", ylab="Log Intensity", xlab="Samples")
dev.off()

#Added code
#Build presence/ abscence call table
calls <-mas5calls(raw_data)
calls <- exprs(calls)
absent <- colSums(calls == 'A')/ nrow(calls)
present <-colSums(calls=='P')/ nrow(calls)
quality_control <- data.frame(absent, present)
write.table(quality_control,"quality_control.txt")

###### 3. DATA NORMALISATION #######################################################################################

## Normalise the data using RMA
normalised_data <- rma(raw_data)
normalised_data
# To obtain a matrix of the expression values, use exprs() 
exp_values <- exprs(normalised_data)


## Control Check - Plot Normalised Data
# Boxplot to observe the results of normalisation
# Notice differences with the boxplot from the raw data
png("data_normalisation.png")
boxplot(exp_values, col=colours,las=2, main="Quality Control After Normalisation", ylab="Log Intensity", xlab="Samples")
dev.off()


# MA plot of all samples 
png("MA_normalised.png")
mva.pairs(exp_values, cex=0.9, main = "MA Plot for Normalised Data")
dev.off()
# The same plot for the non-normalised raw data
png("MA_non_normalised.png")
mva.pairs(pm(raw_data), cex=0.9, main = "MA Plot for Non-Normalised Data")
dev.off()


#process the expression data
#save the data for use in the shiny app
variances <- apply(exp_values, 1, var)
#Rank 20 most expressed genes by highest variance
expression <- exprsvals[order(variances, decreasing=TRUE)[1:20],]

gene_names <- tmp$Symbol[match(rownames(expression), tmp$ID)]
rownames(expression)[! is.na(gene_names)] <- paste(gene_names[
  ! is.na(gene_names)], rownames(expression)[! is.na(gene_names)
  ], sep = ' / ')
save(expression,file="expression.Rdata")


###### 4. HIERARCHICAL CLUSTERING OF NORMALISED DATA####################################################################

# Performs hierarchical clustering with average linkage based on
# Pearsons Correlation Coefficient
hc<-hclust(as.dist(1-cor(exp_values, method="pearson")), method="average")

png("Hierarchical_plot.png")
plot(hc)
dev.off()

###### 5. PRINCIPAL COMPONENTS ANALYSIS ####################################################################

## Perform PCA of normalised data
pca <- prcomp(t(exp_values), scale=T)

# Plot the PCA results
png("PCA_plot.png")
legends <- c("SRF_control","SRF_KO")
percentVar <- round(100 * attr(pca, "percentVar"))
s3d<-scatterplot3d(pca$x[,1:3], pch=19, color=colours, main="Principal Components Analysis", xlab(paste0("PC1: ", percentVar[1], "% variance")))
legend("bottom", legend = adf$Group[1,4],
       col =  colours, pch = 16, 
       inset = -0.1, xpd = TRUE, horiz = TRUE)
s3d.coords <- s3d$xyz.convert(pca$x[,1:3])
text(s3d.coords$x, s3d.coords$y, labels = colnames(exp_values),pos = 3,offset = 0.5)
scores = as.data.frame(pca$x) 
summary(pca)
dev.off()
###### 6. FOLD FILTERING ####################################################################

## Perform fold filtering - Obtaining Expression Values#################################################################

#using expression matrix obtained in step 3 a matrix of expression values

#RMA outputs log2 data while MAS5 outputs linear data
#To convert from log
exprsvals10 <-2^exp_values
#check conversion
exp_values[1:10,]
#converted
exprsvals10[1:10,]

#More fold filtering
#check order of sample names
mysamples <- sampleNames(normalised_data)
#display the list
mysamples
#it is useful to obtain a vector of ProbeIDs here
probesets <- probeNames(raw_data)
#display the first 10 ProbeSets
probesets[1:10]

#Build final fold table
#Calculate the means
#Note mean of the log is not the same as the log of the mean!!
Srf_control.mean <- apply(exprsvals10[,adf$Name[1:3]],1,mean)
Srf_KO.mean <- apply(exprsvals10[,adf$Name[4:6]],1,mean)
#calculate some fold changes
Fold_Change <-Srf_control.mean /Srf_KO.mean
#build a summary table to hold all the data
fold_data= cbind(Srf_control.mean,Srf_KO.mean,Fold_Change)
#check the column names
colnames(fold_data)
#write the table of means as an output
write.table(fold_data,file="fold_data.txt", quote=F,
            sep="\t",col.names=NA)

###### 7. STATISTICAL ANALYSIS OF DIFFERENTIALLY EXPRESSED GENES####################################################################

#Check original sample order
sampleNames(normalised_data)
#Rename the samples
sampleNames(normalised_data) <-c("Srf_WT1","Srf_WT2","Srf_WT3","Srf_KO1","Srf_KO2","Srf_KO3")
#Check the samples have renamed
sampleNames(normalised_data)


##Building annotation for differential gene identification
#establish annotation for MOE430v2
#which annotation do we need
#modified from #http://gettinggeneticsdone.blogspot.co.uk/2012/01/annotating-limma-#results-with-gene.html

normalised_data@annotation

#packages in the annotation package
ls("package:mouse4302.db")

#build an annotation table
ID <- featureNames(normalised_data)
Symbol <- getSYMBOL(ID, "mouse4302.db")
Name <- as.character(lookUp(ID, "mouse4302.db", "GENENAME"))
tmp <- data.frame(ID=ID, Symbol=Symbol, Name=Name, stringsAsFactors=F)
tmp[tmp=="NA"] <- NA #fix padding with NA characters 
#assign as feature data of the current normalised_data
fData(normalised_data) <- tmp

## Statistical analysis using Limma

#Build the design matrix which indicates which RNA sample has been applied to each array
design <- model.matrix(~-1+factor(c(1,1,1,2,2,2)))
colnames(design) <- c("Srf_WT","Srf_KO")
#Check it makes sense
sampleNames(normalised_data)
#output the design matrix
design

#Contrast matrix instructs Limma which comparisons to make between the RNA samples
contrastmatrix <- makeContrasts(Srf_WT-Srf_KO,levels=design)
contrastmatrix

#Based on design matrix: linear model is fitted for each gene on the array 
#To model systematic part of data, so it can be distinguished from random variation 
#and make the contrasts
fit <- lmFit(normalised_data, design)

#Based on contrast matrix
#Obtain coefficients and standard errors for contrasts of the coefficients of the original model
fit2 <- contrasts.fit(fit, contrastmatrix)

#this last part essentially moderates the t-statistic using 
#the borrowed variance approach described in class
#To asses differential expression and compute moderated t-statistics, moderated Fstatistic, and log-odds of differential expression 
fit2 <- eBayes(fit2)

#make a venn diagram
#Multiple Testing Genewise Across Contrasts
#Does not adjust for multiple testing across genes
#nested F-test approach
png("venn_diagram.png")
clas <- classifyTestsF(fit2)
vennDiagram(clas)
dev.off()


#Top differentially expressed genes according to Limma (based on adjusted p-value)
topTable(fit2,coef=1,adjust="fdr")
limma_genes <-topTable(fit2,coef=1, adjust="fdr", number=nrow(normalised_data))
write.table(limma_genes,"limma_genes.txt")

#Top 10 differentially expressed genes according to Limma (based on adjusted p-value)
limma_genes10 <- topTable(fit2,coef=1, adjust="fdr", n=10)
write.table(limma_genes10,"limma_genes10.txt", sep = ",")


#Top differentially expressed genes based on Fold change
#order[-x] as want descending order
# abs(x) as need to ignore -ve and +ve sign for fold change, only interested in absolute value
fold_change_limma10 <- limma_genes[order(-abs(limma_genes$logFC)),]
fold_change_limma10 <- head(fold_change_limma10,10)
write.table(fold_change_limma10,"fold_change_top10.txt", sep = ",")
fold_change_limma10




###### 8. FUNCTIONAL ENRICHMENT ANALYSIS ####################################################################

Mm.H <- readRDS("/shared_files/MSigDB/Mm.h.all.v7.1.entrez.rds") 

#Check that you have the required objects
ls()

#Show the full contents of the annotation package
ls("package:mouse4302.db")

#Show the annotation keys in this database
keytypes(mouse4302.db) 
sampleNames(normalised_data)
## Process annotation for functional enrichment

#Here we select from the annotation a number of keys with the primary key being PROBEID
res <- select(mouse4302.db, keys = rownames(normalised_data), columns = c("ENTREZID", "ENSEMBL","SYMBOL"), keytype="PROBEID")
#View the top of the table
head(res)
#find the index of each row of the expression set in the #annotation object res
idx <- match(rownames(normalised_data), res$PROBEID)
#Use the index to set the phenotypic data in the ExpressionSet
fData(normalised_data) <- res[idx, ]
head(fData(normalised_data), 10)
#Find all rows that dont have an EntrezID and remove then
normalised_data_corrected <- normalised_data[is.na(fData(normalised_data)$ENTREZID)==0,]


## Functional Enrichment Analysis

#convert to indexes
H.indices <- ids2indices(Mm.H,fData(normalised_data_corrected)$ENTREZID)

#if you want to run mroast
enrichment_roast <-mroast(normalised_data_corrected,index=H.indices,design=design,contrast=contrastmatrix[,1],adjust.method = "BH")
#if you want to run camera
enrichment_camera <-camera(normalised_data_corrected,index=H.indices,design=design,contrast=contrastmatrix[,1],adjust.method = "BH")
#if you want to run romer
enrichment_romer <-romer(normalised_data_corrected,index=H.indices,design=design,contrast=contrastmatrix[,1],adjust.method = "BH")
#View the results

write.table(enrichment_camera,"enrichment_camera.txt",sep="\t")
#You can then examine the results in enrichment.txt.  It is a text file.  It can be downloaded to view in a spreadsheet such as Excel.

#results of functional enrichment analysis to be used for shiny app
save(enrichment_camera,file="enrichment_camera.Rdata")
save(enrichment_romer,file="enrichment_romer.Rdata")
save(enrichment_roast,file="enrichment_roast.Rdata")

## Session Information
sessionInfo()

###### 9. SHINY APP ####################################################################
library(shiny)
runApp("~/FGT/ICA/ICA_app/app.R")
###### 9. VOLCANO PLOT ####################################################################
#Code adapted from https://biocorecrg.github.io/CRG_RIntroduction/volcano-plots.html

# The significantly differentially expressed genes are the ones found in the upper-left and upper-right corners.
# Add a column to the data frame to specify if they are UP- or DOWN- regulated (log2FoldChange respectively positive or negative)

# add a column of NAs
limma_genes$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
limma_genes$diffexpressed[limma_genes$logFC > 0.6 & limma_genes$P.Value < 0.05] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
limma_genes$diffexpressed[limma_genes$logFC < -0.6 & limma_genes$P.Value < 0.05] <- "DOWN"

# Re-plot but this time color the points with "diffexpressed"
p <- ggplot(data=limma_genes, aes(x=logFC, y=-log10(P.Value), col=diffexpressed)) + geom_point() + theme_minimal()

# Add vertical lines for log2FoldChange thresholds, and one horizontal line for the p-value threshold 
p2 <- p + geom_vline(xintercept=c(-0.6, 0.6), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red")

## Change point color 

# 1. by default, it is assigned to the categories in an alphabetical order):
p3 <- p2 + scale_color_manual(values=c("blue", "black", "red"))

# 2. to automate a bit: ceate a named vector: the values are the colors to be used, the names are the categories they will be assigned to:
mycolors <- c("blue", "red", "black")
names(mycolors) <- c("DOWN", "UP", "NO")
p3 <- p2 + scale_colour_manual(values = mycolors)

# Now write down the name of genes beside the points...
# Create a new column "delabel" to de, that will contain the name of genes differentially expressed (NA in case they are not)
limma_genes$delabel <- NA
limma_genes$delabel[limma_genes$diffexpressed != "NO"] <- limma_genes$NGenes[limma_genes$diffexpressed != "NO"]

ggplot(data=limma_genes, aes(x=logFC, y=-log10(P.Value), col=diffexpressed, label=delabel)) + 
  geom_point() + 
  theme_minimal() +
  geom_text()

# Finally, we can organize the labels nicely using the "ggrepel" package and the geom_text_repel() function

# plot adding up all layers we have seen so far
ggplot(data=limma_genes, aes(x=logFC, y=-log10(P.Value), col=diffexpressed, label=delabel)) +
  geom_point() + 
  theme_minimal() +
  geom_text_repel() +
  scale_color_manual(values=c("blue", "black", "red")) +
  geom_vline(xintercept=c(-0.6, 0.6), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red")




