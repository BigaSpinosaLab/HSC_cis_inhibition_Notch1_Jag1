################################################################################
##        PAPER : Thambyrajah et al. 2024 Nat Comm
##                "Cis inhibition of NOTCH1 through JAGGED1 sustains embryonic Hematopoietic stem cell fate"
##  Description : This script is for plotting  a heatmap of selected genes
##       Author : María Maqueda (mmaqueda@researchmar.net) at BigSpinLab
################################################################################

# NOTES:
# For generating this plot, it is required to have the normalized expression matrix
# To be downloaded from GEO repositori (Accession ID GSE230790). Additionally,
# Suppl Table 4 is also required (with the results from DEA)

################################################################################
## 1. Load packages
################################################################################

require(ggplot2)
require(ComplexHeatmap)
require(openxlsx)

################################################################################
## 2. Import data: normalized expression matrix and gene annotation
################################################################################

mat <- read.delim(file="/Users/mmaqueda/Downloads/GSE230790_RNAseq_Norm_Counts_Filtered.txt", 
                  header=TRUE, sep="\t")


# Import DEA results to know if that particular gene is DEG from Jag1 vs Washed-out comparison
path <- "../downloads/Suppl_Table_T4.xlsx"

dea.final <- read_excel(path = path, 
                        sheet=excel_sheets(path)[2], # The second sheet is the comparison of interest
                        skip=3) # Skip the first 3 lines (table legend)

################################################################################
## 3. Define genes of interest and subset them from the complete matrix
################################################################################

# Genes of interest
path <-  "Genes_selection_Fig6F.xlsx"
GoI <- lapply(excel_sheets(path), read_excel, path = path)
names(GoI) <- excel_sheets(path)

#names(GoI)
# [1] "JAG1_vs_CompE.Washed_notch"     "JAG1_vs_CompE.Washed_MHC-Histo" "JAG1_VS_CompE_Immune Response" 

 # Rename cases
names(GoI) <- c("Notch Signaling", "MHC-Class I", "Immune Response")


GoI.df <- do.call(rbind.data.frame, GoI)
GoI.df$TOPIC <- c(rep(names(GoI)[1], unlist(lapply(GoI, function(g) nrow(g)))[1]),
                  rep(names(GoI)[2], unlist(lapply(GoI, function(g) nrow(g)))[2]),
                  rep(names(GoI)[3], unlist(lapply(GoI, function(g) nrow(g)))[3]))

# Add corresponding padj from DEA
GoI.df.extended <- merge(GoI.df, dea.final, by="gene_id", )

# Subset them from the matrix
mat.interest <- mat[which(rownames(mat) %in% GoI.df.extended$gene_id),] 

# Sort them in the same order
mat.interest <- mat.interest[match(GoI.df.extended$gene_id, rownames(mat.interest)),]

# Now Change row names to gene names (instead of Ensembl)
rownames(mat.interest) <- GoI.df.extended$gene_name.x

# Add an asterisk in the gene names if they show an adj pval between 0.05 and 0.1
# The rest are below 0.05 (adj pval)
rownames(mat.interest)[which(GoI.df.extended$padj>0.05)] <- paste(rownames(mat.interest)[which(GoI.df.extended$padj>0.05)], "*", sep="")

################################################################################
## 6. Conduct relative values between CompoundE <-> WashOut and JAG1 <-> WashOut
## to consider the embryo variability
################################################################################

mat.interest.relative <- as.data.frame(mat.interest)

mat.interest.relative$C1_W1 <- (mat.interest.relative$I_C1 - mat.interest.relative$I_W1)/mat.interest.relative$I_W1
mat.interest.relative$C2_W2 <- (mat.interest.relative$I_C2 - mat.interest.relative$I_W2)/mat.interest.relative$I_W2
  
mat.interest.relative$J1_W1 <- (mat.interest.relative$I_J1 - mat.interest.relative$I_W1)/mat.interest.relative$I_W1
mat.interest.relative$J3_W3 <- (mat.interest.relative$I_J3 - mat.interest.relative$I_W3)/mat.interest.relative$I_W3

mat.interest.relative <- mat.interest.relative[,c("C1_W1","C2_W2","J1_W1","J3_W3")]
colnames(mat.interest.relative) <- c("CompE_1", "CompE_2", "Fc-JAG1_1", "Fc-JAG1_3") 

# Matrix should be genes x samples
mat.interest.rel.scaled <- t(scale(t(mat.interest.relative)))

# Value range in the matrix
quantile(mat.interest.rel.scaled, c(0.1, 0.95))
# 10%       95% 
#   -0.9907647  1.4359453

# Top annotation: condition they belong

condition <- c("CompE", "CompE", "Fc-JAG1", "Fc-JAG1")
col_fun = circlize::colorRamp2(c(-1.1, 0, 1.5), c("skyblue3", "white", "#67001F"))

htmp <- Heatmap(mat.interest.rel.scaled, 
               name = "Scaled Exprs",  
               column_split = factor(condition, levels=c("CompE","Fc-JAG1")),
               cluster_columns = FALSE,
               cluster_column_slices = TRUE,
               column_title_gp = gpar(fontsize = 7, fontface="bold"),
               column_gap = unit(0.3, "mm"),
               row_split = factor(GoI.df.extended$TOPIC),
               row_gap = unit(0.9, "mm"),
               cluster_rows = TRUE,
               show_row_dend = TRUE,
               col = col_fun,
               top_annotation = HeatmapAnnotation(Condition = factor(condition, levels=c("CompE","Fc-JAG1")),
                                                  col = list(Condition = c("CompE" = "skyblue1", 
                                                                           "Fc-JAG1" = "seagreen2")),
                                                  simple_anno_size=unit(0.3,"cm"),
                                                  show_annotation_name = c(Condition=FALSE),
                                                  show_legend = FALSE),
               show_column_names = TRUE,
               column_names_rot = 45, 
               column_names_gp = gpar(fontsize=6),
               show_row_names = TRUE,
               row_names_gp = gpar(fontsize = 5),
               use_raster = TRUE,
               heatmap_legend_param = list(legend_direction = "horizontal"),
               raster_quality = 4)

pdf("../results/Fig6F_Heatmap_Selected_Genes_RT_Relative2WashOut.pdf", width=2, height=6)
draw(htmp, heatmap_legend_side="bottom")
dev.off()

