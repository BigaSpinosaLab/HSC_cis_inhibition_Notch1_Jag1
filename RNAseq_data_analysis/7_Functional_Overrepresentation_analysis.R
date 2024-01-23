################################################################################
##        PAPER : Thambyrajah et al. 2024 Nat Comm
##                "Cis inhibition of NOTCH1 through JAGGED1 sustains embryonic Hematopoietic stem cell fate"
##  Description : This script is for conducting Overrepresentation analysis over
##                the exclusively up- and down- regulated genes in Fc-JAG1 
##                (vs Washed Out) and the common DEGs between Fc-JAG1 and anti-JAG1 
##       Author : María Maqueda (mmaqueda@researchmar.net) at BigSpinLab
################################################################################

# NOTES:
# Overrepresentation analysis conducted over KEGG Pathways and GO BP databases
# For conducting this analysis, it is required to have the results from Differential
# Expression Analysis (either execute 6_DEA.R script or use corresponding sheets
# from Suppl Table T4)

################################################################################
## 1. Load packages
################################################################################

require(clusterProfiler)
require(org.Mm.eg.db)
require(readxl)
require(openxlsx)

################################################################################
## 2. Import data: DEA results (excel file) with genes to be enriched
################################################################################

# Indicate path where DEA results are stored
# In this case, results enclosed in Suppl Table T4 are imported as an example
path <- "../downloads/Suppl_Table_T4.xlsx"

dea.final <- lapply(excel_sheets(path)[1:2], function(comp)# The first two enclose the DEA results
                    read_excel(path = path, sheet=comp, 
                               skip=3)) # Skip the first 3 lines (table legend)
names(dea.final) <- excel_sheets(path)[1:2]

# Reduce the list of genes to those with adj pval <0.05
dea.final.subset <- lapply(dea.final, function(dea){
  de <- dea[which(dea$padj<0.05),]

  # Return selection
  return(de)
})

################################################################################
## 3. Create a list with three elements, being each element the list of genes to be
## analyzed
################################################################################

# Exclusively DEGs genes in Fc-JAG1 or Common with CompE vs Washed out comparison
exclusive_JAG1 <- setdiff(dea.final.subset$`Fc-Jag1_vs_WashedOut`$Entrez, dea.final.subset$CompE_vs_WashedOut$Entrez)

dea.final.degs <- list("UP_exclusive_Fc-JAG1" = dea.final.subset$`Fc-Jag1_vs_WashedOut` %>%
                         filter(Entrez %in% exclusive_JAG1, Shrunken_lFC > 0), # Sign is the same for shrunken or not,
                       "DOWN_exclusive_Fc-JAG1" = dea.final.subset$`Fc-Jag1_vs_WashedOut` %>%
                         filter(Entrez %in% exclusive_JAG1, Shrunken_lFC < 0),
                       "Common_Fc-JAG1_CompE" = data.frame("Entrez" = intersect(dea.final.subset$`Fc-Jag1_vs_WashedOut`$Entrez, 
                                                                                dea.final.subset$CompE_vs_WashedOut$Entrez)))

################################################################################
## 4. GO BP overrepresentation analysis
################################################################################

ora.GOBP <- lapply(dea.final.degs, function(g){
  enrichGO(gene   = g$Entrez,
           keyType = "ENTREZID",
           OrgDb     = "org.Mm.eg.db",
           ont= "BP",
           pvalueCutoff = 0.25, # To explore a longer list
           qvalueCutoff = 0.25, # To explore a longer list
           pAdjustMethod = "BH",
           readable=TRUE)
})

# Simplify results by semantic similarity
ora.GOBP.reduced <- lapply(ora.GOBP, function(l){
  simplify(l, cutoff=0.7, by="p.adjust", select_fun=min)})


################################################################################
## 5. KEGG overrepresentation analysis
################################################################################

# KEGG requires Entrez id entries

ora.kegg <- lapply(dea.final.degs, function(g){

  y <- enrichKEGG(gene   = g$Entrez,
                  organism     = "mmu", 
                  pvalueCutoff = 1, #Include all KEGG pathways
                  qvalueCutoff = 1,
                  pAdjustMethod = "BH")

  y <- setReadable(y, OrgDb = org.Mm.eg.db, keyType="ENTREZID")
  return(y)
})


################################################################################
## FINAL. Store all results in an excel file
################################################################################

hs <- createStyle(textDecoration = "BOLD", fontColour = "#FFFFFF", fontSize = 12, fontName = "Arial Narrow", fgFill = "#4F80BD")

sapply(seq(1,length(dea.final.degs)), function(stage){
  
  results <- list("GOBP"  = as.data.frame(ora.GOBP[[stage]]),
                  "GOBP.simplified"  = as.data.frame(ora.GOBP.reduced[[stage]]),
                  "KEGG" = as.data.frame(ora.kegg[[stage]]))
  
  write.xlsx(x = results, 
             file= paste("../results/Functional_Overrepresentation_DEGs_",names(dea.final.degs)[stage],".xlsx", sep=""), 
             headerStyle=hs)
})

