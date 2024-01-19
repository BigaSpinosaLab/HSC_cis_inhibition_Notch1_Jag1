################################################################################
##        PAPER : Thambyrajah et al. 2024 Nat Comm
##                "Cis inhibition of NOTCH1 through JAGGED1 sustains embryonic Hematopoietic stem cell fate"
##  Description : Code for FIGURE 2 panels
##       Author : María Maqueda (mmaqueda@researchmar.net) at BigSpinLab
################################################################################

# Following panels from Figure 2 can be obtained from this script:

# B) UMAP coloured by population
# C) UMAP coloured by identified cluster (11 cluster)
# D) Heatmap with selected HSC-HE (Cluster 4) and HSC (Cluster 2) specific genes across all clusters
# E) Force Directed Graph layouts coloured by identified cluster or Jag1 expression 
# F) UMAP coloured by expression level of selected genes
# G) Dotplots showing FACS MFI for NOTCH1, DLL4 and JAG1 for HSC-HE cells or 
#    HSC CD45+ / CD45- cells.

################################################################################
## 1. Load required packages
################################################################################

require(ggplot2)
require(dplyr)
require(ComplexHeatmap)

################################################################################
## 2. Import relevant data and prepare it for plotting
################################################################################

# Required data is deposited in GEO repository under the accession ID GSE230792 
# Please, download the following associated Suppl Files for executing this script
# and import them to R:

# GSE230792_scRNAseq_Cell_annotation.txt.gz
# GSE230792_scRNAseq_Norm_counts_filtered.txt.gz

norm_counts_GEO <- read.table(file = "GEO_upload/GSE230792_scRNAseq_Norm_counts_filtered.txt", 
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
## 4. FIGURE 2B and FIGURE 2C
################################################################################

pdf("tmp/Fig2B_UMAP_complete_Population.pdf", width=2.5, height=2.5)
ggplot(sc.data.all,aes(x=UMAP_C1,y=UMAP_C2, color=population)) +
  geom_point(color= "black", size=0.5) +
  geom_point(aes(color=population),size=0.3)+
  theme_void() +  # Empty theme (no axis at all)
  theme(legend.position = "bottom",  
        legend.spacing.x = unit(0, 'cm'),
        legend.title = element_blank(),
        legend.text = element_text(size=5, face = "bold")) +
  scale_color_manual(values = c("Gfi1_HE" = "green", 
                                "Gfi1_pos_IAHC" = "darkorchid3",
                                "Gfi1_neg_IAHC" = "steelblue1"))  
dev.off()


pdf("tmp/Fig2C_UMAP_complete_Cluster.pdf", width=3, height=2.5)
ggplot(sc.data.all,aes(x=UMAP_C1,y=UMAP_C2)) +
  geom_point(color= "black", size=0.5) +
  geom_point(aes(color=CLUSTER),size=0.3)+
  theme_void() +  # Empty theme (no axis at all)
  theme(legend.position = "right",  
        legend.key.size = unit(0.2,"cm"),
        legend.title = element_blank(),
        legend.text = element_text(size=3, face = "bold")) +
  scale_color_manual(values = c("HSC-HE" = "green4", #green3? 
                                "HSC" = "darkorchid4",
                                "0" = "cadetblue2",
                                "1" = "cadetblue4",
                                "3" = "coral1" ,
                                "5" = "cornsilk2" ,
                                "6" = "darksalmon",
                                "7" = "plum",
                                "8" = "yellow2",
                                "9" = "snow3",
                                "10" = "slategray2")) +
  guides(colour =guide_legend(ncol=1, byrow = FALSE))  +
  annotate("text", x = -3, y = -3, label = "HSC", size=2, color="darkorchid4")  +
  annotate("text", x = -0.5, y = 4, label = "HSC-HE", size=2, color="green4")
dev.off()


################################################################################
## 5. FIGURE 2D
################################################################################

# HSC-HE and HSC specific genes to be plotted
gene.sel <- c("Gfi1", "Procr", "Hlf", "Ptprc", "Runx1", "Neurl3", "Mecom", 
                   "Flt3", "Rfng", "Smad6", "Jag1", "Dll4", "Notch1", "Notch2", 
                   "Dnmt3a", "Rara", "Tfeb", "Mllt3")
gene.sel.cols <- unlist(sapply(gene.sel, function(gene) which(colnames(sc.data.all) %in% gene)))


# Expression matrix with the subset of selected genes
mat <- sc.data.all[,gene.sel.cols]
rownames(mat) <- sc.data.all$Library_name

# Usually this matrix is genes x cell (not the inverse)
mat <- t(mat)

# Scale the rows (genes)
mat <- t(scale(t(mat)))

# Get the cluster annotation
cluster_anno <- sc.data.all$CLUSTER

# Value range in the matrix
quantile(mat, c(0.1, 0.95))
# 10%        95% 
#-0.6976231  2.1518349 

# Color scale with blue to red: make the white color map to 0. the red map to highest and the blue map to the lowest
col_fun = circlize::colorRamp2(c(-1, 0, 3), c("skyblue3", "white", "#67001F"))

htmp <- Heatmap(mat, name = "Scaled Exprs",  
                column_split = factor(cluster_anno, levels=c("HSC","HSC-HE","0","1","3","5","6","7","8","9","10")),
                cluster_columns = FALSE,
                show_column_dend = FALSE,
                cluster_column_slices = TRUE,
                column_title_gp = gpar(fontsize = 10, fontface="bold"),
                column_gap = unit(0.5, "mm"),
                cluster_rows = TRUE,
                show_row_dend = FALSE,
                row_names_side = "left",
                col = col_fun,
                row_names_gp = gpar(fontsize = 9, fontface="bold"),
                top_annotation = HeatmapAnnotation(Louvain_cluster = factor(cluster_anno, levels=c("0","1","3","5","6","7","8","9","10", "HSC-HE","HSC")),
                                                   Timepoint = factor(sc.data.all$developmental_stage),
                                                   col = list(Louvain_cluster = c("HSC-HE" = "green4", #green3? 
                                                                                  "HSC" = "darkorchid4",
                                                                                  "0" = "cadetblue2",
                                                                                  "1" = "cadetblue4",
                                                                                  "3" = "coral1" ,
                                                                                  "5" = "cornsilk2" ,
                                                                                  "6" = "darksalmon",
                                                                                  "7" = "plum",
                                                                                  "8" = "yellow2",
                                                                                  "9" = "snow3",
                                                                                  "10" = "slategray2"),
                                                              Timepoint = c("E10.5" = "darkorange", 
                                                                            "E11.5" = "darkorange4")),
                                                   simple_anno_size=unit(0.3,"cm"),
                                                   show_annotation_name = c(Louvain_cluster=FALSE, Timepoint=FALSE)),
                show_column_names = FALSE,
                use_raster = TRUE,
                heatmap_legend_param = list(legend_direction = "horizontal"),
                raster_quality = 4)


pdf("tmp/Fig2D_Heatmap_sel_genes_all_clusters.pdf", width=8, height=3.3)
draw(htmp, heatmap_legend_side="bottom")
dev.off()


################################################################################
## 6. FIGURE 2E
################################################################################

pdf("tmp/Fig2E_FDG_Cluster_annotation.pdf", width=3, height=2.5)
ggplot(sc.data.all,aes(x=FDG_fa1,y=FDG_fa2)) +
  geom_point(color= "black", size=0.5) +
  geom_point(aes(color=CLUSTER),size=0.3)+
  theme_void() +  # Empty theme (no axis at all)
  theme(legend.position = "right",  
        legend.key.size = unit(0.2,"cm"),
        legend.title = element_blank(),
        legend.text = element_text(size=3, face = "bold")) +
  scale_color_manual(values = c("HSC-HE" = "green4", 
                                "HSC" = "darkorchid4",
                                "0" = "cadetblue2",
                                "1" = "cadetblue4",
                                "3" = "coral1" ,
                                "5" = "cornsilk2" ,
                                "6" = "darksalmon",
                                "7" = "plum",
                                "8" = "yellow2",
                                "9" = "snow3",
                                "10" = "slategray2")) +
  guides(colour =guide_legend(ncol=1, byrow = FALSE))  
dev.off()


pdf("tmp/Fig2E_FDG_Jag1_expression.pdf", width=3, height=2.5)
ggplot(sc.data.all,aes(x=FDG_fa1,y=FDG_fa2)) +
  geom_point(color="black", size=0.5)+
  geom_point(aes_string(color="Jag1"),size=0.3)+
  scale_color_gradient(low = "ivory", high = "red") +
  theme_void() +  # Empty theme (no axis at all)
  theme(legend.position = "bottom",  
        legend.title = element_blank(),
        legend.key.height = unit(0.2,"cm"),
        legend.text = element_text(size=4, face = "bold")) +
  annotate("text", x = -1000, y = 2000, label = "Jag1", size=5)
dev.off()


################################################################################
## 7. FIGURE 2F
################################################################################

gene.sel2 <- c("Mecom", "Procr", "Notch1", "Dll4", "Jag1","Mllt3")

for(gene in gene.sel2)
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
  pdf(paste0("tmp/Fig2F_UMAP_complete_", gene,"_expression_Red_scale.pdf"), width=2.5, height=2.5)
  print(p)
  dev.off()
}

################################################################################
## 8. FIGURE 2G
################################################################################

# Obtain the FACS CD45 threshold to be considered (+/-)
#=====
# Consider CD45+ and CD45- for HSC, both dev stages. Threshold for positive or 
# negative is the Q1 value in HSC cluster cells E11.5

# CD45 FACS values in HSC for E10.5 and E11.5
toexplore <- sc.data.all[which(sc.data.all$CLUSTER == "HSC"),]
by(toexplore[,"FACS_CD45"], toexplore[ ,"developmental_stage"], summary)

# toexplore[, "developmental_stage"]: E10.5
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#   236.9   380.1   701.3  2407.7  1949.0 19826.3       7 
# --------------------------------------------------------------------------------------------------------------------------------- 
# toexplore[, "developmental_stage"]: E11.5
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
# 904.1  9310.4 16427.3 16342.2 22912.7 36117.8       2 

# CD45 to be considered positive for values > Q1 in E11 (9310.4)

sc.data.all <- sc.data.all %>% mutate(CLUSTER.TP = paste(CLUSTER,developmental_stage,sep=":"))
sc.data.all$th_CD45<-ifelse(sc.data.all$FACS_CD45 < 
                              quantile(toexplore$FACS_CD45[which(toexplore$developmental_stage == "E11.5")],
                                                        na.rm=TRUE)[2],
                            c("CD45-"),c("CD45+"))


sc.data.all <- sc.data.all %>% 
  mutate(CLUSTER.TP.CD45 = case_when(
    CLUSTER.TP == "HSC-HE:E10.5"  ~   "HSC-HE:E10.5",
    CLUSTER.TP == "HSC-HE:E11.5"  ~   "HSC-HE:E11.5",
    CLUSTER.TP == "HSC:E10.5"  ~   paste(CLUSTER.TP,th_CD45,sep=":"),
    CLUSTER.TP == "HSC:E11.5"  ~   paste(CLUSTER.TP,th_CD45,sep=":"),
    TRUE ~ "Rest"  # The rest are not of interest indeed for the plotting
  ))


# Cases to plot and order
#=====
FACS <- paste0("FACS_", c("NOTCH1", "JAG1", "DLL4"))

order <- c("HSC-HE:E10.5","HSC-HE:E11.5",
           "HSC:E10.5:CD45-", "HSC:E11.5:CD45-",
           "HSC:E11.5:CD45+")

order.lab <- c("HSC-HE:E10.5","HSC-HE:E11.5",
               "HSC:E10.5\nCD45-", "HSC:E11.5\nCD45-",
               "HSC:E11.5\nCD45+")


# Generate plots
#=====

# NOTE: Dots may be located in different position in x-axis which do not change
# any conclusion. This is due to the jitter

for(gene in FACS)
{
  # Reduce the plotting to HE and HSC cells OR  LEAVE all cells
  toplot <- sc.data.all[which(sc.data.all$CLUSTER %in% c("HSC-HE", "HSC")),]
  toplot$CLUSTER.TP.CD45 <- factor(toplot$CLUSTER.TP.CD45, levels= order)  # "HSC:E10.5:CD45+" Drop it since it only has 3 cells
  toplot <- toplot[!is.na(toplot$CLUSTER.TP.CD45),]
  
  #ggplot(toplot, aes_string(x="CLUSTER.TP", y=gene, fill="CLUSTER.TP")) +
  p <-  ggplot(toplot, aes_string(x="CLUSTER.TP.CD45", y=gene, color="CLUSTER.TP.CD45")) +
    geom_jitter(width=0.3, alpha=0.8, size=0.7, shape=19,  stroke=0.2) +  #color="black",shape=21,
    stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
                 geom="errorbar", color="black", width=0.2, size=0.2) +
    stat_summary(fun.y=mean, geom="point", size=0.2, color="black") +
    theme_classic() +
    theme(legend.position="none",
          plot.title = element_text(size=9),
          axis.title.y=element_text(size=7, colour="black"),
          axis.text.x = element_text(size=3, face="bold", colour = "black"),
          axis.text.y = element_text(size=5, face="bold", colour = "black"),
          axis.ticks.y.left =  element_line(linewidth=0.2),
          axis.ticks.x.bottom =  element_blank(),
          axis.line.x.bottom=element_line(linewidth = 0.2),
          axis.line.y.left = element_line(linewidth = 0.2)) +
    scale_color_manual(values = c("HSC-HE:E10.5" = "green3",
                                  "HSC-HE:E11.5" = "green4",
                                  "HSC:E10.5:CD45-" = "darkorchid2",
                                  "HSC:E11.5:CD45-" = "darkorchid4",
                                  "HSC:E11.5:CD45+" = "darkorchid4")) +
    ggtitle(gene) +
    ylab("FACS index (MFI)") +
    xlab("") +
    scale_x_discrete(labels = order.lab) +
    coord_fixed(ratio=ifelse(gene=="FACS_DLL4", 0.004, ifelse(gene=="FACS_JAG1",0.00016,0.0022))) 
  
  # Generate the pdf and store it
  pdf(paste("tmp/Fig2G_DotPlots_FACS_HE_HSC_CD45_", gene,".pdf", sep=""), width=3.5, height=2.5)
  print(p)
  dev.off()
  
}
