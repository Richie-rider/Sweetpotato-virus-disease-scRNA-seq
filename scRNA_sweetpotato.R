library(SeuratObject)
library(Seurat)
library(patchwork)
library(dplyr)
library(cowplot)
library(ggplot2)
library(openxlsx)
library(harmony)

set.seed(99999)


dir_path = "raw_data/"
dir_list = (list.files(path = dir_path))

my_list = {}
for(i in 1:length(dir_list)){
  # loading each sample
  data_dir = paste(dir_path,"/",dir_list[i],sep="")
  tmp = Read10X(data.dir = data_dir)
  
  # create seurat ojectives
  tmp <- CreateSeuratObject(counts =tmp, project = dir_list[i], min.cells=3, min.features=200)
  tmp <- SCTransform(tmp, verbose = TRUE, do.scale = TRUE)
  
  nExp_poi<- round(0.05*nrow(tmp@meta.data))
  # doublet filter
  tmp <- doubletFinder_v3(tmp, PCs = 1:10, pN = 0.25, pK = 0.09,
                          nExp = nExp_poi,reuse.pANN = FALSE, sct = TRUE)
  

  
  # remove MT, CH genes
  
  chloroplast_genes <- c("psbA","rps16","psbK","psbD","psbC","psaA",
                         "psaB","rbcL","psbB")
  
  mitochondria_genes <- c("trnS-UGA","cox1","cob","ccmFn","ccmFc2",
                          "cox2","rrnL","mttB","cob","rrnL")
  
  protoplast_genes <-  protoplast_genes&Gene_ID
  
  mito_pattern <- paste(mitochondria_genes, collapse = "|")
  chloro_pattern <- paste(chloroplast_genes, collapse = "|")
  protoplast_pattern <-paste(protoplast_genes, collapse = "|")
  
  tmp <-tmp%>% SetIdent(value = colnames( tmp@meta.data[7] ) )%>% subset( ident = "Singlet" )
  tmp [["percent.mt"]] <- PercentageFeatureSet( tmp, pattern = mito_pattern)
  tmp [["percent.ch"]] <- PercentageFeatureSet( tmp, pattern = chloro_pattern)
  tmp [["percent.proto"]] <- PercentageFeatureSet( tmp, pattern = protoplast_pattern)
  

# cut-off the low-qunality cell
  tmp<-subset(tmp,
            subset = nFeature_RNA > 500 & nFeature_RNA < 50000 &
              percent.mt < 5 & percent.ch<10 & nCount_RNA>400 &
              nCount_RNA<50000)
 
  tmp_list={}
  tmp_list[[1]]$data = tmp
  tmp_list[[1]]$name = dir_list[i]
  my_list[[i]]=tmp_list
}


# read the list of cell-cycle genes
cell_cyl<-read.csv("cell_cycle_genes.csv")
# read list of features
feature_list<-read_tsv("1_lib/features.tsv")
# correct the gene names
cycle_geme<-inner_join(cell_cyl, feature_list, by=c("Locus" ="locus"))

# define the list of sgenes and g2m genes
s.genes <- cycle_geme%>%
  filter(Phase=="G1-S")
g2m.genes<- cycle_geme%>%
  filter(Phase=="G2-M")

# create merged 
# work out the effects of cell cycle genes

# check if cc genes are all present in the genes we detected

cc.genes <- c(s.genes$gene_name, g2m.genes$gene_name)
cc.genes[!cc.genes %in% rownames(merged)]

merged <- CellCycleScoring(merged, 
                           s.features = s.genes$gene_name, 
                           g2m.features = g2m.genes$gene_name,
                           set.ident = TRUE)
merged<- SCTransform(merged,
                     vars.to.regress =  c("S.Score", "G2M.Score"),do.scale = TRUE,verbose = TRUE)

merged <- RunPCA(merged, npcs = 30, verbose = TRUE)

merged<- RunUMAP(merged, reduction = "pca", 
                 dims = 1:30, verbose = T)

SCT_result <- IntegrateLayers(object = merged, method = HarmonyIntegration,
                                       orig.reduction = "pca", new.reduction = "harmony",
                                       verbose = FALSE)

SCT_result <- FindNeighbors(SCT_result, reduction = "pca", dims = 1:30)
SCT_result <- FindClusters(SCT_result, resolution = 0.8, cluster.name = "clusters")

SCT_result <- RunUMAP(SCT_result, reduction = "pca", dims = 1:30)

SCT_result <- JoinLayers(SCT_result)

SCT_result <- SetIdent(SCT_result,value = "orig.ident")
SCT_result <- SCT_result %>% RunHarmony("orig.ident", plot_convergence = T)

SCT_result <- SCT_result  %>% 
  RunUMAP(reduction = "harmony", 
          dims = 1:30, 
          verbose = F,
          seed.use = 42) %>% 
  FindNeighbors(reduction = "harmony", k.param = 10, dims = 1:30) %>% 
  FindClusters(resolution = 0.8) %>% 
  identity()


# Optimize the Color panel
library(RColorBrewer)
number_of_clusters <-length(unique(Idents(SCT_result))) 

if (number_of_clusters > 9) {
  color_ramp <- colorRampPalette(brewer.pal(9, "Dark2"))
  colors <- color_ramp(number_of_clusters)
} else {
  colors <- brewer.pal(number_of_clusters, "Dark2")
}

DimPlot(object = SCT_result, reduction = 'umap', raster = FALSE,cols = colors,  
        label = T,
        label.box = T,
        repel = TRUE, pt.size = 0.01,
        shuffle = FALSE)&NoAxes()


# Plot the fractions of cell numbers in each cluster
library(tidyverse)

count_table <- table(Sweetpotato_harmony@meta.data$seurat_clusters, Sweetpotato_harmony@meta.data$sample)
df <- as.data.frame(count_table)
df$count_table <- rowSums(df[,c('JHY', 'VF_JHY')])

# Calculate the fractions for each sample within each cluster
df <- df %>%
  mutate(
    Fraction_VF_JHY = VF_JHY / total_cells,
    Fraction_JHY = JHY / total_cells
  )

# Gather the data to long format for plotting with ggplot2
df_long <- df %>%
  select(Cluster, starts_with("Fraction")) %>%
  pivot_longer(cols = -Cluster, names_to = "Sample", values_to = "Fraction")

df_long$Sample <- sub("Fraction_", "", df_long$Sample)
df_long$Sample <- factor(df_long$Sample, levels = c("VF_JHY", "JHY"))


# Define custom colors for the samples
custom_colors <- c("VF_JHY" = "grey","JHY" = "#499195")

# Plot the fractions as a stacked bar plot with the specified adjustments

ggplot(df_long, aes(x = Cluster, y = Fraction, fill = Sample)) +
  geom_bar(stat = "identity", position="stack") +
  scale_x_continuous(breaks = seq(0, 30, by = 5)) + # Set breaks at intervals of 5
  scale_y_continuous(labels = scales::percent_format()) + # Format the y-axis labels as percentages
  scale_fill_manual(values = custom_colors) + # Assuming you have a vector of colors
  labs(x = "Cluster", y = "Relative fraction of numbers of cells", fill = "Sample") +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(color="black"), 
    axis.text.y = element_text(color="black"),
    axis.line = element_line(colour = "black"),
    axis.ticks = element_line(colour = "black"),
    legend.position = "top",
    legend.text = element_text(size = 10)
  )



# Find cluster markers
SCT_result.markers <- FindAllMarkers(SCT_result, 
                                     only.pos = TRUE, 
                                     min.pct = 0.25, min.diff.pct = 0.25,# at least 25% cells are expressing the marker in one of the groups
                                     logfc.threshold = 0.25) # at least 10% more cells expesssing than the lower group
  
  
# select the top 3 from the markers

library(scCustomize)
library(qs)
library(viridis)
library(magrittr)
library(ComplexHeatmap)

pal <- viridis(n = 10, option = "D")


top_mark <-SCT_result.markers %>%
  group_by(cluster) %>%
  slice_max(n = 3, order_by = avg_log2FC)

# plot the top 10 markers
DotPlot(SCT_result, features = top_makr$gene)

# load cell type markers, and see if the clustering make sense
celltype<-read.csv("1_lib/final_marker.csv")

topmark<-SCT_result.markers%>%
  filter(avg_log2FC>0, p_val_adj<0.05)

DotPlot_scCustom(seurat_object = SCT_result,
                 features = unique(topmark$gene), 
                 flip_axes = T,
                 colors_use = rev(pal))

# need to output a table with reduced genes

allmarkers <- FindAllMarkers(SCT_result, 
                                 only.pos = FALSE, 
                                 min.pct = 0.1, slot = "data",min.diff.pct = 0.2,# at least 10% cells are expressing the marker in one of the groups
                                 logfc.threshold = 0.25) # at least 10% more cells expesssing than the lower group

DEGs_markers<-inner_join(test, celltype, by=c("gene"="gene"))%>%
  filter(avg_log2FC>0)%>%
  mutate(x=(pct.1-pct.2))%>%
  filter(x>0.2)

write.xlsx(top_mark,"output_files/SCT_marker.xlsx")
write.xlsx(mark_celltype, "output_files/marker_celltype.xlsx")

# Add modulescore
SA_response <- c("g30836", "g434", "g13599", "g54917", "g3944", "g61465", "g1137", "g43130", "g11800", "g43655", 
                 "g12899", "g2374", "g35701", "g13601", "g29775", "g25786", "g10000", "g59146", "g53695", "g5258", 
                 "g17257", "g59919", "g56260", "g28379", "g32956", "g46268", "g59229", "g11797", "g34483", "g16706", 
                 "g38090", "g24961", "g1152", "g59480", "g53691", "g47381", "g53696", "g38071", "g1159", "g53937", 
                 "g35683", "g51195", "g13600", "g2961", "g58391", "g61528", "g9104", "g60792", "g50211", "g4174", 
                 "g52724", "g42573", "g52942", "g29421", "g1259", "g11", "g32955", "g54069")
JA_response <- c("g8149", "g34519", "g20273", "g434", "g25973", "g3944", "g19780", "g19461", "g1137", "g43130", 
                 "g11633", "g43655", "g21630", "g46473", "g2374", "g59006", "g29775", "g55166", "g24785", "g20274", 
                 "g46618", "g37258", "g29282", "g24750", "g46309", "g17257", "g28379", "g5362", "g13640", "g43795", 
                 "g13402", "g24961", "g1152", "g44301", "g51475", "g39779", "g53342", "g16326", "g46572", "g54830", 
                 "g11634", "g44497", "g1386", "g21940", "g47381", "g19781", "g31288", "g47479", "g1159", "g53937", 
                 "g34766", "g51195", "g29259", "g29235", "g2961", "g23744", "g36877", "g58391", "g2864", "g14548", 
                 "g40689", "g21147", "g41187", "g50211", "g39764", "g23746", "g31299", "g64290", "g60960", "g16466", 
                 "g56407", "g29421", "g26200", "g37945", "g1259", "g42979", "g11", "g62652", "g55180", "g42719", "g59199")

SCT_result <- AddModuleScore(object = SCT_result,features = list(SA_response),
                                     ctrl = 5,name = "SA_response")
SCT_result <- AddModuleScore(object = SCT_result,features = list(JA_response),
                                     ctrl = 5,name = "JA_response")
roc1 <- c("grey","#EF008C")

FeaturePlot_scCustom (seurat_object= SCT_result, 
                      features = "SA_response1", colors_use = roc1)& NoAxes()&labs(title = "Response to SA")

FeaturePlot_scCustom (seurat_object= SCT_result, 
                      features = "JA_response1", colors_use = roc1)& NoAxes()&labs(title = "Response to JA")


# L_EC: 11,12,13,24
# U_EC: 7,8,17,28
# PP: 21
# Procambium (P) 14
# CC: 18 20
# BSC: 15
# GC: 16,19,23,27
# SMC: 4,5,25
# PMC: 0,2,3,6,9,10,22
# Unkonwn: 1,26

cell_type <- c("PMC","Unknown","PMC","PMC","SMC",
               "SMC","PMC","U_EC","U_EC","PMC",
               "PMC","L_EC","L_EC","L_EC","VC",
               "VC","GC","U_EC","VC","GC","VC",
               "VC","PMC","GC","L_EC","SMC",
               "Unknown","GC","U_EC")

names(cell_type) <-levels(SCT_result)
SCT_result_cell_type <- RenameIdents(SCT_result, cell_type)


# Find the infected cells and DEGs in infected cells vs uninfected cells

Infected_object <- subset(SCT_result_cell_type, subset = SPLCV>0.00000001)
Uninfected_object <- subset(SCT_result_cell_type, subset = SPLCV<0.00000001)

###  Differential Expression in Each cell types of Infected cells

cell_types <- levels(Idents(Infected_object))

# Create an empty list to store the results
de_results <- list()

# Loop through each cell type and perform DE analysis
for(cell_type in cell_types) {
  # Subset the Seurat object for a specific cell type
  subset_infected <- Infected_object
  subset_uninfected <- Uninfected_object
  
  # Run differential expression
  de_results[[cell_type]] <- FindMarkers(subset_infected, ident.1 = cell_type, 
                                         subset.2 = subset_uninfected@meta.data$cell_orig_ident, 
                                         logfc.threshold = 0.25,
                                         min.pct = 0.25,
                                         min.diff.pct =0.1)
}

# Specify the directory 
output_directory <- "~output_files/DEGs_viral_infected_vs_uninfected"

# Loop through the list of DEGs and write each to a set of excel files
for (cell_type in names(de_results)) {
  # Create a file name based on the cell type
  output_file_path <- file.path(output_directory, paste0("DEGs_", cell_type, ".xlsx"))
  
  # Write to xlsx
  write.xlsx(de_results[[cell_type]], output_file_path, rowNames = TRUE)
}



# Pesedutime analysis

library(monocle)
library(dplyr)
my_cds<-readRDS("monocle_dataset/all")

my_cds <- estimateSizeFactors(my_cds)
my_cds <- detectGenes(my_cds, min_expr = 0.1) # i

cds_ch_all<- my_cds 
cds_ch_all <- detectGenes( cds_ch_all, min_expr = 0.1)

# introduce a framework
expressed_genes <- row.names( subset(fData (my_cds), num_cells_expressed >= 10))


DEGs<-differentialGeneTest(cds_ch_all, 
                              fullModelFormulaStr = '~original_clusters', 
                              reducedModelFormulaStr = "~replicate",
                              cores = 8)

top300<-row.names(DEGs_wt )[order(DEGs$qval)][1:300]

# for combined

cds_ch_all<-cds_ch_all[expressed_genes, ]

DEGs_ch_all <- differentialGeneTest(cds_ch_all, 
                                    fullModelFormulaStr = '~treament+original_clusters', 
                                    reducedModelFormulaStr = "~replicate",
                                    cores = 8)



order_genes<-DEGs_ch_all[!DEGs_ch_all%in%top300]

top2000_ch_all <- row.names( order_genes )[ order( order_genes$qval )][1:2000]


cds_ch_all <- setOrderingFilter(cds_ch_all, 
                                ordering_genes = top2000_ch_all) # set the order using top 3000


#finish the introduction

cds_ch_all <- reduceDimension(cds_ch_all, 
                              method = 'DDRTree') # demension reduction

cds_ch_all <- orderCells(cds_ch_all) 


write.xlsx(top50_wt,"order.xlsx")
write.xlsx(top3000_ch_all,"order_gene.xlsx")

my_pseudotime_de <- differentialGeneTest(cds_ch_all,
                                         fullModelFormulaStr = "~sm.ns(Pseudotime)",
                                         cores = 8)

write.xlsx(my_pseudotime_de,"sodu_genes_all.xlsx")



#Plot the trajectory curve
# Extract UMAP data
library(ggplot2)
library(dplyr)

umap_data <- Embeddings(SCT_result, reduction = "umap")
umap_df <- as.data.frame(umap_data)
colnames(umap_df) <- c("UMAP_1", "UMAP_2")

# Add barcodes as a new column
umap_df$Barcode <- rownames(umap_df)
umap_df$Barcode <- sub("_1$", "", umap_df$Barcode)


pseudotime_table <- read.xlsx("~/sodu_genes_all.xlsx")

# Merge pseudotime data with UMAP coordinates
combined_data <- merge(umap_df, pseudotime_table, by = "Barcode")

combined_data <- combined_data %>%arrange(Pseudotime) 

# Calculate a smooth trajectory
library(splines)
combined_data$smoothed_UMAP_1 <- predict(loess(UMAP_1 ~ Pseudotime, data = combined_data, span = 0.1))
combined_data$smoothed_UMAP_2 <- predict(loess(UMAP_2 ~ Pseudotime, data = combined_data, span = 0.1))

# Plotting the smoothed trajectory
ggplot(combined_data, aes(x = Pseudotime)) +
  geom_density(fill = "#499195", alpha = 0.5) +  # Adjust color and transparency
  labs(title = "Density of Cells by Pseudotime",
       x = "Pseudotime",
       y = "Density") +
  theme_minimal()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())


# scRNA-seq coexpression analysis
library(tidyverse)
library(cowplot)
library(patchwork)
library(WGCNA)
library(hdWGCNA)
library(igraph)
theme_set(theme_cowplot())
set.seed(12345)
enableWGCNAThreads(nThreads = 8)

Sweet_WGCNA <- SetupForWGCNA(
  SCT_result_cell_type,
  gene_select = "fraction", # the gene selection approach
  fraction = 0.05, # fraction of cells that a gene needs to be expressed in order to be included
  wgcna_name = "Sweet_WGCNA" # the name of the hdWGCNA experiment
)


# construct metacells  in each group
Sweet_WGCNA <- MetacellsByGroups(
  seurat_obj = Sweet_WGCNA,
  group.by = c("Cell_Type", "sample"), # specify the columns in seurat_obj@meta.data to group by
  reduction = 'harmony', # select the dimensionality reduction to perform KNN on
  k = 25, # nearest-neighbors parameter
  max_shared = 10, # maximum number of shared cells between two metacells
  ident.group = 'Cell_Type' # set the Idents of the metacell seurat object
)

# normalize metacell expression matrix:
Sweet_WGCNA<- NormalizeMetacells(Sweet_WGCNA)

Sweet_WGCNA <- SetDatExpr(
  Sweet_WGCNA,
  group_name = "PMC", # the name of the group of interest in the group.by column
  group.by='Cell_Type', # the metadata column containing the cell type info. This same column should have also been used in MetacellsByGroups
  assay = 'RNA', # using RNA assay
  slot = 'data' # using normalized data
)


Sweet_WGCNA <- TestSoftPowers(
  Sweet_WGCNA,
  networkType = 'unsigned' # you can also use "unsigned" or "signed hybrid"
)

# plot the results:
plot_list <- PlotSoftPowers(Sweet_WGCNA)

# assemble with patchwork
wrap_plots(plot_list, ncol=2)

power_table <- GetPowerTable(Sweet_WGCNA)
head(power_table)


# construct co-expression network:
Sweet_WGCNA <- ConstructNetwork(
  Sweet_WGCNA,
  tom_name = 'PMC' # name of the topoligical overlap matrix written to disk
)

PlotDendrogram(Sweet_WGCNA , main='PMC hdWGCNA Dendrogram')

# compute all MEs in the full single-cell dataset
Sweet_WGCNA <- ModuleEigengenes(
  Sweet_WGCNA,
  group.by.vars="sample"
)

# harmonized module eigengenes:
hMEs <- GetMEs(Sweet_WGCNA )
# module eigengenes:
MEs <- GetMEs(Sweet_WGCNA, harmonized=FALSE)


# compute eigengene-based connectivity (kME):
Sweet_WGCNA  <- ModuleConnectivity(
  Sweet_WGCNA ,
  group.by = 'Cell_Type', group_name = 'PMC'
)

Sweet_WGCNA  <- ResetModuleNames(
  Sweet_WGCNA,
  new_name = "PMC-M"
)

PlotKMEs(Sweet_WGCNA, ncol=5)


library(UCell)
Sweet_WGCNA <- ModuleExprScore(
  Sweet_WGCNA,
  n_genes = 25,
  method='UCell'
)


# make a featureplot of hMEs for each module
plot_list <- ModuleFeaturePlot(
  Sweet_WGCNA,
  features='hMEs', # plot the hMEs
  order=TRUE # order so the points with highest hMEs are on top
)

# stitch together with patchwork
wrap_plots(plot_list, ncol=6)

ModuleCorrelogram(Sweet_WGCNA)

# Individual module network plots
ModuleNetworkPlot(
  Sweet_WGCNA,
  outdir = 'ModuleNetworks'
)


ModuleNetworkPlot(
  Sweet_WGCNA, 
  outdir='ModuleNetworks2', # new folder name
  n_inner = 20, # number of genes in inner ring
  n_outer = 30, # number of genes in outer ring
  n_conns = Inf, # show all of the connections
  plot_size=c(10,10), # larger plotting area
  vertex.label.cex=1 # font size
)

# Applying UMAP to co-expression networks
Sweet_WGCNA<- RunModuleUMAP(
  Sweet_WGCNA,
  n_hubs = 5, # number of hub genes to include for the UMAP embedding
  n_neighbors=10, # neighbors parameter for UMAP
  min_dist=0.1 # min distance between points in UMAP space
)


# get the hub gene UMAP table from the seurat object
umap_df <- GetModuleUMAP(Sweet_WGCNA)

# plot with ggplot
ggplot(umap_df, aes(x=UMAP1, y=UMAP2)) +
  geom_point(
    color=umap_df$color, # color each point by WGCNA module
    size=umap_df$kME*2 # size of each point based on intramodular connectivity
  ) +
  umap_theme()


ModuleUMAPPlot(
  Sweet_WGCNA,
  edge.alpha=0.25,
  sample_edges=TRUE,
  edge_prop=0.1, # proportion of edges to sample (20% here)
  label_hubs=2 ,# how many hub genes to plot per module?
  keep_grey_edges=FALSE
)




