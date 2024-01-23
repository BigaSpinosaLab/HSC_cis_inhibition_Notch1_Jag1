################################################################################
##        PAPER : Thambyrajah et al. 2024 Nat Comm
##                "Cis inhibition of NOTCH1 through JAGGED1 sustains embryonic Hematopoietic stem cell fate"
##  Description : This script is for conducting a Differential Expression Analysis
##       Author : María Maqueda (mmaqueda@researchmar.net) at BigSpinLab
################################################################################

# NOTES: 
# a) Required input data for this analysis is an expression matrix (RAW counts).
#    This matrix can be directly downloaded from GEO repository (GSE230790) or
#    created/reproduced using original FASTQ files + shared shell scripts. 
# b) Comparisons are considering CompE.washed condition as the reference
# c) By executing this script, DEA results in Suppl Table T4 are obtained

################################################################################
## 1. Load packages
################################################################################

require("openxlsx")
require("DESeq2")
require("EnhancedVolcano")
require("vsn")
require(dplyr)
require("ggpubr")

################################################################################
## 2. Import relevant data and prepare it for analysis
################################################################################

# Raw counts matrix from GEO repo > GSE230790_RNAseq_Raw_counts.txt.gz
counts <- read.delim(file="/Users/mmaqueda/Downloads/GSE230790_RNAseq_Raw_counts.txt", 
                     header=TRUE, sep="\t")

# Create dataframe with samples metadata information (Embryo and Condition)
# colnames(counts)
# [1] "I_C1" "I_C2" "I_J1" "I_J3" "I_W1" "I_W2" "I_W3"

targets <- data.frame("Sample.id" = colnames(counts),
                      "Embryo" = c("e1", "e2", "e1", "e3", "e1", "e2", "e3"),
                      "Condition" = c(rep("CompE",2), 
                                      rep("Jag1",2), 
                                      rep("CompE.Washed",3)))

# Gene annotation (dim 55487x14): Created from corresponding GTF file from Ensembl
# Shown an entry example
gene.annot <- readRDS(file = "/Volumes/projectscomput/cancer/db_files/Genomes/Ensembl/mouse/mm10/release-102/Mus_musculus.GRCm38.102.ONLY_GENES.RDS")

# chrom  source feature   start     end score strand frame                                                                                               attribute
# 1     ensembl    gene 3102016 3102125     .      +     . gene_id ENSMUSG00000064842; gene_version 1; gene_name Gm26206; gene_source ensembl; gene_biotype snRNA;

# gene_id gene_type gene_name    Entrez gene_biotype
# ENSMUSG00000064842      <NA>   Gm26206 115487594        snRNA

################################################################################
## 3. Reduce count matrix: gene filtering. Genes must show a minimal expression in 
## at least two replicates. 
################################################################################

# Genes expressed in at least two replicates from the same condition
# By expressed we apply: >=10 raw counts in that replicate

conditions <- unique(targets$Condition)

keep.genes <- lapply(conditions, function(case){
   data <- counts[,targets$Condition %in% case]

   # Option A: two or more replicates per condition
   keep <- apply(data[,], 1, function(row) {
     a <- sum(row >= 10)
     return(ifelse(a>=2, TRUE, FALSE))
   })

   return(keep)
 })

kk <- Reduce(`|`, keep.genes)  # OR among all conditions
counts.red <- counts[kk,] # 10311 genes kept

# Reduce gene annotation to the one in reduced expression matrix according to the reduced
# count matrix. Remove discarded genes: they are directly in the same order

gene.annot.red <- gene.annot[which(gene.annot$gene_id %in% rownames(counts.red)),]

################################################################################
## 4. Create object raw and VST normalized (for plotting if required)
################################################################################

data.dds <- DESeq2::DESeqDataSetFromMatrix(countData = counts.red, 
                                           colData = targets, 
                                           design =  ~ Embryo + Condition, 
                                           tidy=FALSE)

#dim(data.dds) 
#colData(data.dds)
#head(assay(data.dds,"counts"))

# Transformation: for visualization and clustering purposes

#  Not blind to the experimental design: otherwise dispersion would be overestimated
vstData.local.dds <- vst(data.dds, blind=FALSE, fitType='local') 
#msdvstlocal <- meanSdPlot(assay(vstData.local.dds))
#msdvstlocal$gg + ggtitle("VST transformation - Local fT") # Plot it!

vstData.dds <- vstData.local.dds

################################################################################
## 5. DEA with DESeq2 and store results 
################################################################################

# Let's relevel the group factor (reference) to get the comparison of interest

colData(data.dds)$Condition <- relevel(colData(data.dds)$Condition, ref = "CompE.Washed")

# Apply DESeq2
dea.dds <- DESeq(data.dds)
resultsNames(dea.dds)  #Identify the name of the comparison of interest
# [1] "Intercept"                       "Embryo_e2_vs_e1"                 "Embryo_e3_vs_e1"                
# [4] "Condition_CompE_vs_CompE.Washed" "Condition_JAG1_vs_CompE.Washed" 

# Define comparisons of interest
comparisons <- c("Condition_CompE_vs_CompE.Washed", "Condition_JAG1_vs_CompE.Washed")

# IMPORTANT: Execute complete_dea() function in the next section prior to execute
# following instruction
all.comparisons <- lapply(comparisons, function(comp){
  dea <- complete.dea(DESeqSet = dea.dds,
                      comparison = comp,
                      fc_threshold  = 5, 
                      gene_annotations = gene.annot.red,
                      suffix.fig.files = "",
                      lim.ma.abs = 10,
                      norm.data= vstData.dds,
                      res.path = "../results/")
  return(dea)
})

names(all.comparisons) <- comparisons

# Generate a reduce list of DEA results by discarding those genes that are only expressed in one replicate

all.comparisons.reduced <- lapply(all.comparisons, function(res){
  
  # Select the integer columns (raw counts)
  res.raw <- res %>% dplyr::select(where(is.integer))  
  
  # Identify the different conditions under test
  conditions <- targets$Condition[which(targets$Sample.id %in% colnames(res.raw))]

  # Identify expressed genes per condition
  keep.genes <- lapply(unique(conditions), function(case){
    res.raw.cond <- res.raw[,conditions %in% case]
    
    # Option A: two or more replicates per condition
    keep <- apply(res.raw.cond[,], 1, function(row) {
      a <- sum(row >= 10)  # Sum up number of positions in that row with > 10 raw counts 
      return(ifelse(a>=2, TRUE, FALSE))
    })
    return(keep)
  })
  
  # OR among all conditions
  kk <- Reduce(`|`, keep.genes)  

 # Return logical vector for all genes to keep
  return(res[kk,])
  })

# Save results if needed
hs <- createStyle(textDecoration = "BOLD", fontColour = "#FFFFFF", fontSize = 12, fontName = "Arial Narrow", fgFill = "#4F80BD")
names(all.comparisons.reduced) <- sapply(names(all.comparisons.reduced), function(n) sub(pattern = "Condition_", replacement = "",x = n))

# (This is Suppl Table T4)
write.xlsx(all.comparisons.reduced,
           file = "../results/DEA_all_comparisons_REDUCED_CompEWashed_REF.xlsx", 
           headerStyle = hs)

# Function to compute DEA and lfcShrink and add annotations ---------------

################################################################################
## ANNEX. Function to compute DEA, lfcShrinkage, gather all information and add annotations
##    into the same data frame and generate a Volcano plot + excel file with results
################################################################################

# NOTE: Assumed to have ENSEMBL ids in the gene annotations. If not, function to
# be reviewed

complete.dea <- function(DESeqSet, #Obtain after computing DESeq function
                         comparison, 
                         gene_annotations, # File with the corresponding gene annotations 
                         suffix.fig.files,
                         fc_threshold, # Shrunken fold change threshold for Volcano Plot
                         lim.ma.abs= 10,  #Limit to consider in MA plot
                         res.path,  #Path to store MA plot and excel file with the results from DEA
                         norm.data) # VST normalized data 
{
  # 1. Obtain results from DEA
  results_dea <-results(DESeqSet,
                        name =  comparison,
                        alpha = 0.05)
  
  # 2. Print summary of results on the console
  print(comparison)
  summary(results_dea)
  
  # 3. LogFC shrinkage
  res.shr <- lfcShrink(dds = DESeqSet, 
                       coef= comparison,
                       res = results_dea, 
                       type = "apeglm")
  
  # 4. Store an MA-plot to see the effect
  pdf(file = file.path(res.path,paste("MAplot_",comparison,suffix.fig.files,".pdf",sep="")), width=16, height=10)
  par(mfrow=c(1,2))    
  DESeq2::plotMA(results_dea,ylim=c(-lim.ma.abs,lim.ma.abs),main=paste("DEA ",comparison,suffix.fig.files, sep=""))
  DESeq2::plotMA(res.shr, ylim=c(-lim.ma.abs,lim.ma.abs), main=paste("DEA ", comparison,suffix.fig.files," - Shrinkage performed", sep=""))
  dev.off()
  
  # 5. Collect results in the same data frame: DEA resuls and logFC shrinkage
  DEA_complete <- cbind(as.data.frame(results_dea), 
                        "Shrunken_lFC" = as.data.frame(res.shr)$log2FoldChange,
                        "Shrunken_lFCSE" = as.data.frame(res.shr)$lfcSE)
  
  # 6. Gene annotation to be included from the same gtf file used during alignment
  DEA_complete <- cbind(gene_annotations[,c("gene_id", "gene_name","Entrez")],DEA_complete)

  # 7. Volcano plot
  pdf(file = file.path(res.path,paste("Volcano_",comparison,suffix.fig.files,".pdf",sep="")), width=16, height=10)
  p <- EnhancedVolcano(toptable = DEA_complete,
                       lab= DEA_complete$gene_name,
                       x ='Shrunken_lFC',
                       y= 'padj',
                       xlab = "Shrunken fold change",
                       ylab = "-Log10(adj pval)",
                       title = comparison,
                       subtitle= "Differential Expression - AGM Notch1-Jag1 Hematopoietic",
                       pCutoff = 0.05,
                       max.overlaps = 100,
                       labSize=3,
                       drawConnectors = TRUE,
                       FCcutoff = fc_threshold,
                       caption = paste("FC cutoff: ", fc_threshold,"; adj p-val cutoff: 0.05",sep=""),
                       legendLabels=c("Not sign.", "Shrunken FC", "adj pval","adj pval & shrunken FC"),
                       legendPosition = "right",
                       legendIconSize = 5.0,
                       legendLabSize = 10)
  print(p)
  dev.off()
  
  
  # 8. Add raw and normalized counts for those samples
  samples <- unlist(strsplit(x = comparison, split = "_"))[c(2,4)]
  DEA_complete <- cbind(DEA_complete,
                        counts(DESeqSet)[,which(colData(DESeqSet)$Condition %in% samples)],  # Raw
                        assay(norm.data)[,which(colData(DESeqSet)$Condition %in% samples)])  # Norm
  
  # 7. Sort results by p-adj
  DEA_complete <- DEA_complete[order(DEA_complete$padj),]
  
  # 8. Store results in an excel file
  hs <- createStyle(textDecoration = "BOLD", fontColour = "#FFFFFF", fontSize = 12, fontName = "Arial Narrow", fgFill = "#4F80BD")
  write.xlsx(DEA_complete,file = file.path(res.path,paste("DEA_All_genes_",comparison,suffix.fig.files,".xlsx",sep="")), headerStyle = hs)
  
  # 8. Return complete dataframe
  return(DEA_complete)
}


