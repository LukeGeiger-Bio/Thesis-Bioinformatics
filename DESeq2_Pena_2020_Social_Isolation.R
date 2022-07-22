#See this website for further resources on importing stuff to/using R on a computing cluster/Unix:
#https://uwaterloo.ca/statistics-and-actuarial-science/research/resources/r-tutorial-unix-environment


if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.14") #silly program wants earlier version

BiocManager::install("DESeq2") #force u to install 
install.packages("data.frame")

#Library loading
library(tidyverse)
library(magrittr)
library(dplyr)
library(DESeq2)
library(RColorBrewer)
library(pheatmap)
library(ggplot2)

#Making metadata file to concatenate sample runs

cut_metadata <- All_Metadata[,c("Run","Sample.Name","acute.chronic","housing","region","sex")] 
cut_metadata #All_Metadata is a delimited file with all identifiers, metadata

##################.tsv file from Andrew w/combined all runs together

colnames(PenaData)[1] <- "Gene"
rownames(PenaData) <- PenaData$Gene
rownames(PenaData) <- sub("gene-", "", rownames(PenaData)) #Removes excess information so names in metadata match count matrix
CountMatrix <- PenaData[,2:639]
CountMatrix

rownames(cut_metadata) <-cut_metadata$Run #Make metadata,count matrix names into row and column names

cut_metadata_ordered <- cut_metadata[order(cut_metadata$Run),]#Consistent order or later steps fail
CountMatrix_ordered <-CountMatrix %>% select(order(colnames(.)))

flipped_CountMatrix <- t(CountMatrix_ordered)#To combine, so that rownames are same in metadata + matrix
#Don't have pseudoreplicates in your count matrices; combine all runs by individual! 

#######################do below in gen-comp2
merge_matrix <- merge(cut_metadata_ordered, flipped_CountMatrix, by ='row.names', all = FALSE)#Merges metadata and counts
merge_matrix
merge_grouped <-merge_matrix %>% group_by(Sample.Name)#Groups counts for runs by identifier
select_grouped <-merge_grouped[c(3:40206)]#Pulls everything but run names bc we've grouped by animal, double check column numbers
select_grouped
sum_matrix <-select_grouped %>% group_by(Sample.Name) %>% summarise_at(vars(c(6:40206)), sum)#A bit redundant to group again tbh
sum_matrix#This variation of summarise (summarise_at(vars,)) works on GC2 as of July 2022
#######################do this in gen-comp2^^


coldata <-cut_metadata_ordered#Making our standard variables for DESeq input
cts <-cts_sum_m
colnames(cts) <- cts[1,]
cts <- cts[2:40206,]#Cleaning redundant names

coldata <-coldata[,2:6]
coldata <- unique(coldata)
rownames(coldata) <-coldata[,1]



coldata$acute.chronic <-factor(coldata$acute.chronic) #Make metadata into factor levels for factor comparison later
coldata$housing <-factor(coldata$housing)
coldata$region <-factor(coldata$region)
coldata$sex <-factor(coldata$sex)

all(rownames(coldata) %in% colnames(cts))
rownames(coldata) %in% colnames(cts) #Checks for naming convention identicality 
colnames(cts) %in% rownames(coldata)

cts <- cts %>% select(order(colnames(.)))
coldata <- coldata %>% arrange(rownames(.)) #Useful way to fix out-of-order names

#Output from gen-comp2 will return a matrix with factor levels. I haven't found a way to undo this leveling in base R. 
#Best workaround is to download to loacal machine, then re-upload into the environment as a new delimited file. 
#Your output matrix is still unable to be used for DESeq2 ('negative' values == factor leveling) if you see spaces in matrix

write.table(cts, "C:/Users/trott/OneDrive/Desktop/cts_sum_m.txt")

cts <-cts_sum_m
coldata

#create dds object
dds <- DESeqDataSetFromMatrix(countData = cts,
                               colData = coldata,
                               design = ~sex + acute.chronic + region + housing)#Final factor is basis of comparison, 
                                                                                #other covariates come before
                                                                                # recall that X:Y is an interaction variable
dds
#count filtering in a row
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]


dds <- DESeq(dds)#Run DESeq2 with contrast of interest
res <- results(dds, contrast=c("housing","SI","G")) #LFC = log(SI/G)
resOrdered <- res[order(res$pvalue),]#Order by whatever you want
resOrdered

design(dds)
summary(res)
resultsNames(dds)#Did my comparisons run properly? 
cts
coldata


#######Graphing and exporting the results! 

hist(v_res$pvalue)#Check for issues with algorithm

###Sample heatmap, row = gene bins, color = normalized scaled expression:
ntd <- normTransform(dds)
select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:20]
df <- as.data.frame(colData(dds)[,c("acute.chronic","housing","region","sex")])
pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=TRUE, annotation_col=df)#Options for the dendrogram and sample clustering by similarity





vsd <- vst(dds, blind=FALSE)#Data variance stabilizing transformation
sampleDists <- dist(t(assay(vsd)))

###Sample-sample similarity heatmap:
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$condition, vsd$type, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

###Principal component analysis:
plotPCA(vsd, intgroup=c("housing", "region"))
plotPCA(vsd, intgroup=c("housing", "sex"))

###Volcano plots:
vp <-data.frame(res@rownames,res$pvalue,res$log2FoldChange)
vp
vp$diffexpressed <- "NO"
vp$diffexpressed[res$log2FoldChange > 0.6 & res$pvalue < 0.1] <- "UP"
vp$diffexpressed[res$log2FoldChange < -0.6 & res$pvalue < 0.1] <- "DOWN"
vp$vplabel <- NA
vp$vplabel[vp$diffexpressed != "NO"] <- res@rownames[vp$diffexpressed != "NO"]

ggplot(data=vp, aes(x=res$log2FoldChange, y=-log10(res$pvalue), col=diffexpressed, label=vplabel)) + 
  geom_point() + 
  theme_minimal() +
  geom_vline(xintercept=c(-0.6, 0.6), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red")



###Obtaining DEG list
up <-vp[vp$diffexpressed =='UP',]
down <-vp[vp$diffexpressed =='DOWN',]

write.csv(up, "C:/Users/trott/OneDrive/Desktop/up.csv")
write.csv(down, "C:/Users/trott/OneDrive/Desktop/down.csv")
write.csv(vp, "C:/Users/trott/OneDrive/Desktop/vp.csv")







#As above, so below:






##################
#amygdala analysis
amy_cts <-cts[,1:57]
amy_coldata <-coldata[1:57,]

a_dds <- DESeqDataSetFromMatrix(countData = amy_cts,
                              colData = amy_coldata,
                              design = ~sex + acute.chronic + housing)

#count filtering in a row
a_keep <- rowSums(counts(a_dds)) >= 10
a_dds <- a_dds[a_keep,]


a_dds <- DESeq(a_dds)#Run DESeq2 with contrast of interest
a_res <- results(a_dds, contrast=c("housing","SI","G"))
a_resOrdered <- a_res[order(a_res$pvalue),]
a_res_p <- a_res[a_res@listData$pvalue<.1,]
a_resOrdered
a_res_p
write.csv(a_res_p, "C:/Users/trott/OneDrive/Desktop/a_res_p.csv")


design(a_dds)
summary(a_res)
resultsNames(a_dds)

vp_a <-data.frame(a_res@rownames,a_res$pvalue,a_res$log2FoldChange)
vp_a
vp_a$diffexpressed <- "NO"
vp_a$diffexpressed[a_res$log2FoldChange > 0.6 & a_res$pvalue < 0.1] <- "UP"
vp_a$diffexpressed[a_res$log2FoldChange < -0.6 & a_res$pvalue < 0.1] <- "DOWN"
vp_a$vplabel <- NA
vp_a$vplabel[vp_a$diffexpressed != "NO"] <- a_res@rownames[vp_a$diffexpressed != "NO"]

library(ggplot2)
ggplot(data=vp_a, aes(x=a_res$log2FoldChange, y=-log10(a_res$pvalue), col=diffexpressed, label=vplabel)) + 
  geom_point() + 
  theme_minimal() +
  geom_vline(xintercept=c(-0.6, 0.6), col="red") +
  geom_hline(yintercept=-log10(0.1), col="red")

#################
#nucleus accumbens analysis
nac_cts <-cts[,58:112]
nac_coldata <-coldata[58:112,]

n_dds <- DESeqDataSetFromMatrix(countData = nac_cts,
                                colData = nac_coldata,
                                design = ~sex + acute.chronic + housing)

#count filtering in a row
n_keep <- rowSums(counts(n_dds)) >= 10
n_dds <- n_dds[n_keep,]


n_dds <- DESeq(n_dds)#Run DESeq2 with contrast of interest
n_res <- results(n_dds, contrast=c("housing","SI","G"))
n_resOrdered <- n_res[order(n_res$padj),]
n_res_p <- n_res[na.omit(n_res@listData$pvalue)<.1,]

n_resOrdered
n_res_p
write.csv(n_res_p, "C:/Users/trott/OneDrive/Desktop/n_res_p.csv")

design(n_dds)
summary(n_res)
resultsNames(n_dds)
?results

vp_n <-data.frame(n_res@rownames,n_res$pvalue,n_res$log2FoldChange)
vp_n
vp_n$diffexpressed <- "NO"
vp_n$diffexpressed[n_res$log2FoldChange > 0.6 & n_res$pvalue < 0.1] <- "UP"
vp_n$diffexpressed[n_res$log2FoldChange < -0.6 & n_res$pvalue < 0.1] <- "DOWN"
vp_n$vplabel <- NA
vp_n$vplabel[vp_n$diffexpressed != "NO"] <- n_res@rownames[vp_n$diffexpressed != "NO"]

library(ggplot2)
ggplot(data=vp_n, aes(x=n_res$log2FoldChange, y=-log10(n_res$pvalue), col=diffexpressed, label=vplabel)) + 
  geom_point() + 
  theme_minimal() +
  geom_vline(xintercept=c(-0.6, 0.6), col="red") +
  geom_hline(yintercept=-log10(0.1), col="red")

##################
#prefrontal cortex analysis
pfc_cts <-cts[,113:166]
pfc_coldata <-coldata[113:166,]

p_dds <- DESeqDataSetFromMatrix(countData = pfc_cts,
                                colData = pfc_coldata,
                                design = ~sex + acute.chronic + housing)

#count filtering in a row
p_keep <- rowSums(counts(p_dds)) >= 10
p_dds <- p_dds[p_keep,]


p_dds <- DESeq(p_dds)#Run DESeq2 with contrast of interest
p_res <- results(p_dds, contrast=c("housing","SI","G"))
p_resOrdered <- p_res[order(p_res$padj),]
p_res_p <- p_res[na.omit(p_res@listData$pvalue)<.1,]

p_resOrdered
p_res_p
write.csv(p_res_p, "C:/Users/trott/OneDrive/Desktop/p_res_p.csv")

design(p_dds)
summary(p_res)
resultsNames(p_dds)

vp_p <-data.frame(p_res@rownames,p_res$pvalue,p_res$log2FoldChange)
vp_p
vp_p$diffexpressed <- "NO"
vp_p$diffexpressed[p_res$log2FoldChange > 0.6 & p_res$pvalue < 0.1] <- "UP"
vp_p$diffexpressed[p_res$log2FoldChange < -0.6 & p_res$pvalue < 0.1] <- "DOWN"
vp_p$vplabel <- NA
vp_p$vplabel[vp_p$diffexpressed != "NO"] <- p_res@rownames[vp_p$diffexpressed != "NO"]

library(ggplot2)
ggplot(data=vp_p, aes(x=p_res$log2FoldChange, y=-log10(p_res$pvalue), col=diffexpressed, label=vplabel)) + 
  geom_point() + 
  theme_minimal() +
  geom_vline(xintercept=c(-0.6, 0.6), col="red") +
  geom_hline(yintercept=-log10(0.1), col="red")

############
#ventral tegmental area analysis
vta_cts <-cts[,166:223]
vta_coldata <-coldata[166:223,]

v_dds <- DESeqDataSetFromMatrix(countData = vta_cts,
                                colData = vta_coldata,
                                design = ~sex + acute.chronic + housing)

#count filtering in a row
v_keep <- rowSums(counts(v_dds)) >= 10
v_dds <- v_dds[v_keep,]


v_dds <- DESeq(v_dds)#Run DESeq2 with contrast of interest
v_res <- results(v_dds, contrast=c("housing","SI","G"))
v_resOrdered <- v_res[order(v_res$pvalue),]
v_res_p <- v_res[na.omit(v_res@listData$pvalue)<.1,]

v_resOrdered
v_res_p
write.csv(v_res_p, "C:/Users/trott/OneDrive/Desktop/v_res_p.csv")

design(v_dds)
summary(v_res)
resultsNames(v_dds)

vp_v <-data.frame(v_res@rownames,v_res$pvalue,v_res$log2FoldChange)
vp_v
vp_v$diffexpressed <- "NO"
vp_v$diffexpressed[v_res$log2FoldChange > 0.6 & v_res$pvalue < 0.1] <- "UP"
vp_v$diffexpressed[v_res$log2FoldChange < -0.6 & v_res$pvalue < 0.1] <- "DOWN"
vp_v$vplabel <- NA
vp_v$vplabel[vp_v$diffexpressed != "NO"] <- v_res@rownames[vp_v$diffexpressed != "NO"]

library(ggplot2)
ggplot(data=vp_v, aes(x=v_res$log2FoldChange, y=-log10(v_res$pvalue), col=diffexpressed, label=vplabel)) + 
  geom_point() + 
  theme_minimal() +
  geom_vline(xintercept=c(-0.6, 0.6), col="red") +
  geom_hline(yintercept=-log10(0.1), col="red")

############


