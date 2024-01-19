################################################################################
##        PAPER : Thambyrajah et al. 2024 
##                "Cis inhibition of NOTCH1 through JAGGED1 sustains embryonic Hematopoietic stem cell fate"
##  Description : Code for FIGURE 3 panels
##       Author : María Maqueda (mmaqueda@researchmar.net) at BigSpinLab
################################################################################

# Following panels from Figure 3 can be obtained from this script:

# A) UMAPs coloured by expression level of selected genes
# B) Venn Diagram with DEGs in E10.5 and E11.5 for GFI1+HE population Hey1/2+ vs Hey1/2- cells
# C) and D) Volcano plots from DEA for previous scenarios
# E) Barplot KEGG pathways from overrepresentation analysis
# F) Violin plots and boxplots of expression levels of selected genes 
# G) UMAP subset with HSC-HE and HSC cluster showing NOTCH1 and JAG1 MFI associated levels

################################################################################
## 1. Load required packages
################################################################################

require(ggplot2)
require(dplyr)
require(readxl)
require(ggvenn)
require(EnhancedVolcano)
require(reshape)
require(ggpubr)

################################################################################
## 2a. Import relevant data and prepare it for plotting
################################################################################

# Required data is deposited in GEO repository under the accession ID GSE230792 
# Please, download the following associated Suppl Files for executing this script
# and import them to R:

# GSE230792_scRNAseq_Cell_annotation.txt.gz
# GSE230792_scRNAseq_Norm_counts_filtered.txt.gz

norm_counts_GEO <- read.table(file = "GEO_upload/scRNAseq_Norm_counts_filtered.txt", 
                              header=TRUE, check.names=FALSE) # 30362 genes x 775 cells
cell_metadata_GEO <- read.table(file = "GEO_upload/GSE230792_scRNAseq_Cell_annotation_final.txt", 
                                header=TRUE, check.names=FALSE) # 860 cells x 22 obs

# NOTE: Cell metatada file includes information from additional cells that were 
#       discarded during QC. Those discarded cells show NAs in columns such as
#       Louvain cluster

# Merge data into one dataframe where cells are in rows and rest of information
# (metadata and expression levels) in columns to ease plotting

  # Discard cells not passing QC
cell_metadata_pass <- cell_metadata_GEO[cell_metadata_GEO$Library_name %in% colnames(norm_counts_GEO),]
norm_counts_GEO_transposed <- as.data.frame(t(norm_counts_GEO))
norm_counts_GEO_transposed$Cell_name <- rownames(norm_counts_GEO_transposed)

  # Merge data. This merged object is the required object for all plotting
sc.data.all <- merge(x = cell_metadata_pass, y = norm_counts_GEO_transposed, 
                     by.x = "Library_name", by.y= "Cell_name")

################################################################################
## 2b. Import Suppl Table T3 (refer to paper) to be used in these figures
################################################################################

# Import Suppl Table T3 with identified results from DEA and KEGG overreperesention analysis
path <- "../PAPER_review_NatComm/Final_Review_Jan24/REVISIONS_JAN24/Suppl_Table_T3_DEGs_Hey12_vs_rest_Gfi1_HE_REVIEWED_MM_Jan24_v2.xlsx"
suppT3 <- lapply(excel_sheets(path), read_excel, path = path)

# DEA results and select DEGs
dea <- suppT3[c(2,4)]  # DEA results
names(dea) <- excel_sheets(path)[c(2,4)]

# Reduce the list of genes to those with abs(log2FC) >2 and adj pval <0.05
dea.final.degs <- lapply(dea, function(df){
  #Remove NA entries (not tested)
  df <- df[!is.na(df[,"Expressed_p"]),]
  
  # Select DEGs
  df <- df %>%
    dplyr::filter(abs(Expressed_l) >2, adj.pval <0.05)
  
  # Return selection
  return(df)
})


# KEGG overrepresentation results
kegg <- suppT3[c(3,5)]  # KEGG overrepresentation results
names(kegg) <- excel_sheets(path)[c(3,5)]

# Skip first lines since they are header
kegg <- lapply(kegg, function(data) {
  colnames(data) <- data[5,]
  data <- data[-c(1:5),]
  return(data)
})

# Combine results from both stages
kegg$Overrepresent_KEGG_DEGs_E10.5$Stage <- rep("E10.5", nrow(kegg$Overrepresent_KEGG_DEGs_E10.5))
kegg$Overrepresent_KEGG_DEGs_E11.5$Stage <- rep("E11.5", nrow(kegg$Overrepresent_KEGG_DEGs_E11.5))

kegg.combined <- rbind(kegg$Overrepresent_KEGG_DEGs_E10.5, 
                       kegg$Overrepresent_KEGG_DEGs_E11.5)

################################################################################
## 3. Relabel cluster names
################################################################################

# Cluster relabeling for HSC and HE
# For the moment, we will keep the rest with just the numbers
sc.data.all <- sc.data.all %>%
  mutate(CLUSTER = case_when(
    Louvain_cluster == 4  ~   "HSC-HE",
    Louvain_cluster == 2  ~   "HSC",
    Louvain_cluster == 0  ~   "0",
    Louvain_cluster == 1  ~   "1",
    Louvain_cluster == 3  ~   "3",
    Louvain_cluster == 5  ~   "5",
    Louvain_cluster == 6  ~   "6",
    Louvain_cluster == 7  ~   "7",
    Louvain_cluster == 8  ~   "8",
    Louvain_cluster == 9  ~   "9",
    Louvain_cluster == 10  ~   "10",
    # TRUE ~ "Rest",  # in case we would like to join some in a 'Rest' cluster
  ))

################################################################################
## 4. FIGURE 3A
################################################################################

gene.sel<- c("Hes1", "Hey1", "Hey2", "Gata2")

for(gene in gene.sel)
{
  # Generate the individual plots 
  p <- ggplot(sc.data.all,aes_string(x="UMAP_C1",y="UMAP_C2", color=gene))+
    geom_point(color="black", size=0.5)+
    geom_point(aes_string(color=gene),size=0.3)+
    scale_color_gradient(low = "ivory", high = "red") +
    theme_void() +  # Empty theme (no axis at all)
    theme(legend.position = "bottom",  
          legend.title = element_blank(),
          legend.key.height = unit(0.2,"cm"),
          legend.text = element_text(size=4, face = "bold")) +
    annotate("text", x = -8, y = 4, label = gene, size=5)
  
  # Generate the pdf and store it
  pdf(paste0("tmp/Fig3A_UMAP_complete_", gene,"_expression_Red_scale.pdf"), width=2.5, height=2.5)
  print(p)
  dev.off()
}

################################################################################
## 4. FIGURE 3B
################################################################################

tovenn <- lapply(dea.final.degs, function(l) l$Expressed_n)
names(tovenn) <- c("E10.5","E11.5")

pdf(file="tmp/Fig3B_VennDiagram_DEGs_E10_E11_Hey1Hey2Exprs_vs_NotExprs.pdf",width=1.5,height=1.5)
ggvenn(
  tovenn, 
  fill_color = c('darkorange', 'darkorange4'),
  set_name_size = 3,
  text_size = 1.5,
  stroke_size = 0.25,fill_alpha = 0.8,
) + ggtitle("DEGs (abs logFC>2 & adj pval<0.05)")
dev.off()

################################################################################
## 5. FIGURE 3C and FIGURE 3D
################################################################################

# Specify genes to highlight (labelled) in the plots
highlight <- list("E10.5" =c("Mecom", "Hey1", "Hey2", "Smad6", "Cdkn1c", "Cd44",
                             "Smad6", "Jag1", "Ctnnb1", "Procr", "Gata3", "Notch1", 
                             "Fgd5", "Sox7", "Hoxa9", "Akt3","Prom1", "Fos", "Egfl7",
                             "Cyp26b1"),
                  "E11.5" =c("Hey1", "Hey2", "Mecom", "Jag1", "Ltbp4","Mllt3", "Smad7",
                             "Smad6", "Notch1", "Procr", "Fbln5", "Mcm2", "Prdm16",
                             "Myc", "Zeb2", "Foxc2", "Gata3", "Dnmt3a", "Prom1", 
                             "Hoxa9","Bcl2"))

# Prepare data and generate plots

for(i in 1:length(dea))
{
  # Not to plot the ones that were not be tested based on their cells percentage
  # presence
  toplot <- as.data.frame(dea[[i]][!is.na(dea[[i]][,3]),]) 
  
  # Create a new column indicating that gene is an outlier in terms of log2FC (abs>10)
  # or adj pval (<1e-12). This is just for the plot
  toplot$Outlier <-ifelse((abs(toplot$Expressed_l) > 10) | (toplot$adj.pval < 1e-12), "Yes", "No")
  
  # Saturate outlier values to 10 maximum (abs) or 1e-12 adj pval
  toplot$Expressed_l[toplot$Expressed_l>10] <- 10
  toplot$Expressed_l[toplot$Expressed_l< (-10)] <- -10
  toplot$adj.pval[toplot$adj.pval< 1e-12] <- 1e-12
  
  # Create custom key-value colors for logFC >1 or 2 and significant genes
  keyvals <- ifelse(
    (toplot$Expressed_l > 2) & (toplot$adj.pval < 0.05), 
    ifelse(names(highlight)[i] == "E10.5",'darkorange', 'darkorange4'),
    'grey')
  names(keyvals)[keyvals != 'grey'] = 'up'
  names(keyvals)[keyvals == 'grey'] <- 'rest'
  
  # Create custom key-value shape for those outlier genes
  keyvals.shape <- ifelse(toplot$Outlier == "Yes",17 , 19)
  
  names(keyvals.shape)[keyvals.shape == 19] <- 'NOoutlier'
  names(keyvals.shape)[keyvals.shape == 17] <- 'outlier'
  
  
  # Build the volcano
  p <- EnhancedVolcano(toplot,
                       lab = toplot[,1],
                       x = colnames(toplot)[2],  # Log2 fold change
                       y = colnames(toplot)[6],  # adj pval (FDR)
                       pCutoff = 0.05,
                       FCcutoff = 2,
                       xlim = c(-10, 10), #c(-6, 6) ZOOMed version; c(-32, 32) COMPLETE
                       ylim = c(0, 12),  #c(0, 15) This is a zoomed version; c(0, 25) COMPLETE
                       pointSize = 0.35, 
                       labSize = 2.5,
                       axisLabSize = 12,
                       captionLabSize = 5,
                       titleLabSize = 7,
                       borderWidth = 0.4,
                       #labFace="bold",
                       selectLab = highlight[[i]],
                       boxedLabels = FALSE,  # With TRUE is not possible to include so many labels
                       title = names(dea)[i],
                       subtitle = "",
                       caption = "FC cutoff: 2; adj p-val cutoff: 0.05",
                       legendPosition = "none",
                       shapeCustom = keyvals.shape,
                       colCustom = keyvals,
                       colAlpha = 1,
                       drawConnectors = TRUE,
                       max.overlaps = 30, # By default 15
                       widthConnectors = 0.2,
                       lengthConnectors = unit(0.005, "npc"),
                       typeConnectors = "open",
                       cutoffLineType = 'blank',
                       hline = c(0.05),
                       hlineCol = c('black'),
                       hlineType = c('longdash'),
                       hlineWidth = c(0.25),
                       vline = c(2),
                       vlineCol = c('black'),
                       vlineType = c('longdash'),
                       vlineWidth = c(0.25),
                       gridlines.major = FALSE,
                       gridlines.minor = FALSE)
  
  
  p <- p + ggplot2::theme(
    axis.title.y=element_text(size=7, colour="black"),
    axis.title.x=element_text(size=7, colour="black"),
    axis.text.x = element_text(size=7, face="bold", colour = "black"),
    axis.text.y = element_text(size=7, face="bold", colour = "black"),
    axis.ticks.y.left =  element_line(linewidth=0.25),
    axis.ticks.x.bottom =  element_line(linewidth=0.25))
  
  pdf(file=paste("tmp/Fig3CD_VOLCANO_",names(dea)[i],"_log2FC_2.pdf", sep=""),width=3,height=3.5)
  print(p)
  dev.off()
}


################################################################################
## 5. FIGURE 3E
################################################################################

# Remove from the pathway Description 'Mus musculus (house mouse)'
kegg.combined$Description <- sapply(kegg.combined$Description, function(des) unlist(strsplit(des, split=" - Mus musculus (house mouse)", fixed=TRUE)))

# Define pathways to show
paths <- c("PI3K-Akt signaling pathway", "ECM-receptor interaction", 
           "Notch signaling pathway", "TGF-beta signaling pathway",
           "Fluid shear stress and atherosclerosis", "DNA replication",
           "Wnt signaling pathway", "TNF signaling pathway", "Cell cycle",
           "Signaling pathways regulating pluripotency of stem cells",
           "Sphingolipid signaling pathway")

highlight <- kegg.combined[which(kegg.combined$Description %in% paths),]

pdf("tmp/Fig3E_Barplot_Selected_Overrepresented_KEGG_E10_E11_Hey1Hey2_DEA.pdf", width=2.5, height=3)
ggplot(data=highlight,aes(x=Description, y=-log10(as.numeric(p.adjust)), fill=Stage)) +
  geom_bar(stat="identity", position=position_dodge(preserve="single"),   # position=position_dodge(preserve="single")
           width = 0.8) +
  coord_flip() +
  theme_bw() +
  theme(legend.position = c(0.8,0.5),  
        legend.title = element_blank(),
        legend.key.size = unit(0.3,"cm"),
        legend.text = element_text(size=6, face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border=element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(),
        axis.text.y = element_text(size=6),
        axis.title.x =  element_text(size=7),
        axis.title.y =  element_blank()) +
  scale_fill_manual(values = c("E10.5" = "darkorange", 
                               "E11.5" = "darkorange4")) +
  scale_x_discrete(labels = function(x) 
    stringr::str_wrap(x, width = 23))+
  geom_hline(yintercept=-log10(0.05), linetype="dashed", color = "black",size=0.4)

dev.off()


################################################################################
## 6. Prepare data for Figure 3F and Figure 3G
################################################################################

# Create subset for the clusters of interest
sc.subsets <- list("HE" = sc.data.all %>% filter(Louvain_cluster==4), #HE  76 cells
                   "HSC" = sc.data.all %>% filter(Louvain_cluster==2))  #HSC  109 cells

# Genes to plot (Main Notch1 targets included)
main.targets <- c("Hes1","Gata2", "Hey1", "Hey2")


# Additional information: cell double positive for FACS JAG1 and NOTCH1
# This is required for Fig 3G

sc.subsets <- lapply(sc.subsets, function(subset){
  subset$Double_Pos_JAG1_NOTCH1 <- ifelse(subset$FACS_NOTCH1>500,
                                          ifelse(subset$FACS_JAG1 > 11000, "POSITIVE", "NEGATIVE"), 
                                          "NEGATIVE")
  subset$Double_Pos_JAG1_NOTCH1 <- factor(subset$Double_Pos_JAG1_NOTCH1, ordered=T)
  return(subset)
})

#table(sc.subsets$HE$Double_Pos_JAG1_NOTCH1)
# NEGATIVE POSITIVE 
# 54       21 

#table(sc.subsets$HSC$Double_Pos_JAG1_NOTCH1)
# NEGATIVE POSITIVE 
# 88       12


# Prepare data frame for plotting
genes <- main.targets

# Select expression matrix and prepare df
sc.subsets.genes <- lapply(seq(1, length(sc.subsets)), function(df){
  genes <- c(genes, "Double_Pos_JAG1_NOTCH1","UMAP_C1", "UMAP_C2")
  df.red <- sc.subsets[[df]][ , which(colnames(sc.subsets[[df]]) %in% genes)]
  df.red$Cell <- sc.subsets[[df]]$Row.names
  df.red$Cluster <- rep(names(sc.subsets)[df], nrow(sc.subsets[[df]]))
  df.red <- melt(df.red, id.vars=c("Cluster", "Double_Pos_JAG1_NOTCH1","UMAP_C1","UMAP_C2"),variable_name="Gene")
  return(df.red)
})
names(sc.subsets.genes) <- names(sc.subsets)

# Combine both clusters
complete.data <- do.call(rbind, sc.subsets.genes)

################################################################################
## 7. Figure 3F. Remember to execute previous section
################################################################################

for(i in 1:length(genes))
{
  toplot <- complete.data %>%
    filter(Gene==genes[i])
  
  set.seed(123) # For the jitter
  
  # Common characteristics for the plots
  
  p <- ggplot(toplot, aes(Cluster, value, fill=Cluster)) +
    geom_violin(alpha=0.5, lwd=0.2)  +
    geom_boxplot(width=0.05,lwd=0.1,outlier.size = 0.1) +
    geom_jitter(width=0.25,
                size=0.2,
                alpha=0.8,
                shape=21,
                color="black", stroke=0.1) +
    stat_compare_means(method = "wilcox.test",  label.x = 0.6, label.y = max(toplot$value) +0.5, label="p.format", size=1.5) +
    theme_classic() +
    theme(legend.position="none",
          plot.title = element_text(size=6),
          axis.title.y=element_text(size=4, colour="black"),
          axis.text.x = element_text(size=4, face="bold", colour = "black"),
          axis.text.y = element_text(size=4, face="bold", colour = "black"),
          axis.ticks.y.left =  element_line(linewidth=0.2),
          axis.ticks.x.bottom =  element_blank(),
          axis.line.x.bottom=element_line(linewidth = 0.2),
          axis.line.y.left = element_line(linewidth = 0.2)) +
    scale_fill_manual(values = c("HE" = "green4", 
                                 "HSC" = "darkorchid4")) +
    ggtitle(genes[i]) +
    ylab("Norm Expression") +
    xlab("")
  
  if(genes[i] != "Hey2")  # Add violin plot to boxplot
    p <- p + geom_violin(alpha=0.5, lwd=0.2)
  
  pdf(file=paste("tmp/Fig3F_",genes[i],"_expression_HE_HSC_Violins_Boxplots.pdf", sep=""), width=1.5, height=1.5)
  print(p)
  dev.off()
}

################################################################################
## 8. FIGURE 3G. Remember to execute section 6)
################################################################################

complete.data <- complete.data %>%
  mutate(CL_FACS = case_when(
    Cluster == "HE" & Double_Pos_JAG1_NOTCH1 == "POSITIVE" ~   "HE_pos",
    Cluster == "HSC" & Double_Pos_JAG1_NOTCH1 == "POSITIVE" ~   "HSC_pos",
    Cluster == "HE" & Double_Pos_JAG1_NOTCH1 == "NEGATIVE" ~   "HE_neg",
    Cluster == "HSC" & Double_Pos_JAG1_NOTCH1 == "NEGATIVE" ~   "HSC_neg",
    TRUE ~ "Unknown", 
  ))

pdf("tmp/Fig3G_UMAP_HSC_HE_Double_positives_JAG1_NOTCH1.pdf", width=3, height=2.5)
ggplot(complete.data,aes(x=UMAP_C1,y=UMAP_C2)) +
  geom_point(color= "black", size=0.5) +
  geom_point(aes(color=CL_FACS),size=0.3)+
  theme_void() +  # Empty theme (no axis at all)
  theme(legend.position = "right",  
        legend.key.size = unit(0.2,"cm"),
        legend.title = element_blank(),
        legend.text = element_text(size=3, face = "bold")) +
  scale_color_manual(values = c("HE_pos" = "green4", 
                                "HSC_pos" = "darkorchid4",
                                "HE_neg" = "snow2",
                                "HSC_neg" = "snow2",
                                "Unknown" = "white")) +
  guides(colour =guide_legend(ncol=1, byrow = FALSE))  +
  annotate("text", x = -3, y = -3, label = "HSC", size=2, color="darkorchid4")  +
  annotate("text", x = -0.5, y = 4, label = "HE", size=2, color="green4")
dev.off()

