library(circlize)
library(CARD)
library(ComplexHeatmap)
library(dplyr)
library(FNN)
library(ggplot2)
library(ggpubr)
library(Matrix)
library(patchwork)
library(reshape2)
library(rstatix)
library(Seurat)
library(SingleCellExperiment)
library(spacexr)
library(SpatialExperiment)
library(spdep)
library(SPOTlight)
library(SummarizedExperiment)
library(tidyr)
library(viridis)
####load spatial data####
dir="/Users/mingkewu/Documents/R/benchmark/sampleresult/"
base_dir <- "~/Documents/R/benchmark/data"
sample_dirs <- list.dirs(base_dir, full.names = TRUE, recursive = FALSE)
spe_list <- lapply(sample_dirs, function(d) {
  obj <- Load10X_Spatial(data.dir = d)
  obj$sample <- basename(d)
  return(obj)
})
names(spe_list) <- basename(sample_dirs)
spe_merged <- merge(
  x = spe_list[[1]],
  y = spe_list[-1],
  add.cell.ids = names(spe_list)
)
table(spe_merged$sample)
base_dir <- "/Users/mingkewu/Documents/R/benchmark/data1"
files <- list.files(base_dir)
samples <- unique(gsub("-filtered_feature_bc_matrix\\.h5|-spatial", "", files))
samples <- samples[nchar(samples) > 0]
spe_list <- list()
for(sample in samples) {
  h5_file <- list.files(base_dir, 
                        pattern = paste0(sample, ".*filtered_feature_bc_matrix\\.h5$"),
                        full.names = TRUE)
  spatial_dir <- list.files(base_dir,
                            pattern = paste0(sample, ".*-spatial$"),
                            full.names = TRUE)
  if(length(h5_file) > 0 && length(spatial_dir) > 0) {
    counts <- Read10X_h5(h5_file)
    image_dir <- spatial_dir[1]
    image <- Read10X_Image(image_dir)
    spe_list[[sample]] <- CreateSeuratObject(counts = counts)
    spe_list[[sample]][[paste0("slice1")]] <- image
    spe_list[[sample]]$sample <- sample
  }
}
spe_merged <- merge(spe_list[[1]], y = spe_list[-1], 
                    add.cell.ids = names(spe_list))
for(i in 1:length(names(spe_list))) {
  names(spe_merged@images)[i] <- names(spe_list)[i]
}
spe_merged[["Spatial"]] <- spe_merged[["RNA"]]
DefaultAssay(spe_merged) <- "Spatial"
spe_merged[["RNA"]] <- NULL
data_path <- "/Users/mingkewu/Documents/R/benchmark/data2/"
folders <- list.dirs(data_path, recursive = FALSE, full.names = TRUE)
all_objs <- list()
for (folder in folders) {
  data <- Read10X(file.path(folder, "filtered_feature_bc_matrix"))
  if (is.list(data)) {
    mat <- data$`Gene Expression`
  } else {
    mat <- data
  }
  h5_file <- file.path(folder, "filtered_feature_bc_matrix.h5")
  if (file.exists(h5_file)) {
    file.remove(h5_file)
  }
  write10xCounts(h5_file, mat)
  spatial_file <- file.path(folder, "spatial", "tissue_positions.csv")
  file.remove(spatial_file)
  obj <- Load10X_Spatial(
    data.dir = folder,
    filename = "filtered_feature_bc_matrix.h5"
  )
  sample_name <- basename(folder)
  obj$sample <- sample_name
  all_objs[[sample_name]] <- obj
}
merged_obj <- merge(all_objs[[1]], all_objs[-1], add.cell.ids = names(all_objs))
for(i in 1:length(names(all_objs))) {
  names(merged_obj@images)[i] <- names(all_objs)[i]
}
combined_obj <- merge(spe_merged, y = merged_obj)
spe=combined_obj 
DefaultAssay(spe)="Spatial"
table(spe$sample)
saveRDS(spe,paste0(dir,"merge.rds"))

####seurat-based label transfer####
path <- "lung_all/"
counts <- readMM(file.path(path, "counts.mtx"))
features <- fread(file.path(path, "features.tsv"), header = FALSE)
rownames(counts) <- features$V1
barcodes <- fread(file.path(path, "barcodes.tsv"), header = FALSE)
colnames(counts) <- barcodes$V1
metadata <- fread(file.path(path, "metadata.csv"))
metadata <- as.data.frame(metadata)
colnames(metadata)[1:6] <- c("barcode", "donor", "celltype", "disease", "sex", "tissue")
rownames(metadata) <- metadata[,1]  
metadata <- metadata[match(colnames(counts), rownames(metadata)), ]
counts <- counts[!duplicated(rownames(counts)), ]
adata <- CreateSeuratObject(counts = counts, meta.data = metadata)
obj=adata
obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
VlnPlot(obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0, raster = FALSE)
obj <- subset(obj, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & 
                 nCount_RNA > 500 & nCount_RNA < 30000 & 
                 percent.mt < 10)
obj$celltype=as.character(obj$celltype)
obj=JoinLayers(obj)
table(obj$donor)
obj[["RNA"]] <- split(obj[["RNA"]], f = obj$donor)
obj <- NormalizeData(obj)
obj <- FindVariableFeatures(obj)
obj <- ScaleData(obj)
obj <- RunPCA(obj)
obj <- IntegrateLayers(
  object = obj, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "harmony",
  verbose = T
)
obj <- FindNeighbors(obj, dims = 1:50, reduction = "harmony")
obj <- FindClusters(obj, resolution =0.2)
obj <- RunUMAP(obj, reduction = "harmony", dims = 1:50, reduction.name = "umap",min.dist = 0.1) 
obj=JoinLayers(obj)
obj <- FindClusters(obj, resolution =0.5)
table(obj$seurat_clusters)
obj$seurat_clusters <- as.character(obj$seurat_clusters)
cluster2type <- c(
  "7"  = "T_NK",
  "2"  = "T_NK",
  "0"  = "T_NK",
  "16" = "T_NK",
  "14" = "T_NK",
  "11" = "T_NK",
  "4"  = "B",
  "1"  = "Myeloid",
  "5"  = "Myeloid",
  "12" = "Mast",
  "9"  = "Stroma",
  "6"  = "Endothelial",
  "10"= "Plasma",
  "8"="AT2-like",
  "13"="AT1-like",
  "3"="Epi-like"
)
obj$anno <- ifelse(
  obj$seurat_clusters %in% names(cluster2type),
  cluster2type[obj$seurat_clusters],
  obj$seurat_clusters
)
obj <- subset(obj, (seurat_clusters %in% c(0:14,16)))
new_levels <- c("AT1-like", "AT2-like", "Epi-like","B", "Endothelial", 
                "Mast", "Myeloid", 
                "Plasma", "Stroma", "T_NK")
obj$anno <- factor(obj$anno, levels = new_levels)
saveRDS(obj,"~/Documents/R/lung/lung_all2.RDS")
sce <- readRDS("~/Documents/R/lung/lung_all2.RDS")
spe <- readRDS(paste0(dir,"merge.rds"))
sce <- NormalizeData(sce)
sce <- FindVariableFeatures(sce)
sce <- ScaleData(sce)
sce <- RunPCA(sce, npcs = 50)
spe.list <- SplitObject(spe, split.by = "sample")
predictions.list <- list()
for (i in seq_along(spe.list)) {
  sample_name <- names(spe.list)[i]
  cat("Processing sample:", sample_name, "\n")
  query <- spe.list[[i]]
  query <- NormalizeData(query)
  query <- FindVariableFeatures(query)
  anchors <- FindTransferAnchors(reference = sce, query = query, reduction = "rpca", dims = 1:30)
  predictions <- TransferData(anchorset = anchors, refdata = sce$anno, dims = 1:30)
  query <- AddMetaData(query, metadata = predictions)
  query$seurat_trans <- query$predicted.id
  predictions.list[[sample_name]] <- query
}
spe <- merge(predictions.list[[1]], predictions.list[-1])
seurat_matrix_list <- list()
for (i in seq_along(predictions.list)) {
  pred <- predictions.list[[i]]@meta.data
  pred_mat <- as.matrix(pred[, grep("prediction.score", colnames(pred))])
  colnames(pred_mat) <- gsub("prediction.score.", "", colnames(pred_mat))
  seurat_matrix_list[[i]] <- pred_mat
}
seurat_matrix <- do.call(rbind, seurat_matrix_list)
seurat_order <- c("AT1.like", "AT2.like", "B", "Endothelial", "Epi.like", "Mast", "Myeloid", "Plasma", "Stroma", "T_NK")
seurat_matrix <- seurat_matrix[, seurat_order]
cell_sums <- rowSums(seurat_matrix)
table(spe$seurat_trans)
saveRDS(spe,paste0(dir,"merge1.rds"))
saveRDS(seurat_matrix,paste0(dir,"matrix1.rds"))
seurat_trans_annotation <- spe$seurat_trans
saveRDS(seurat_trans_annotation, paste0(dir,"/seurat_trans_annotation.rds"))

####card-based deconvolution####
spe <- readRDS(paste0(dir,"merge.rds"))
spe.list <- SplitObject(spe, split.by = "sample")
spatial_count.list <- lapply(spe.list, function(x) {
  LayerData(x, assay = "Spatial", layer = "counts")
})
spatial_loc.list <- lapply(spe.list, function(x) {
  coords <- GetTissueCoordinates(x, image = names(x@images)[1])
  coords <- coords[colnames(LayerData(x, assay = "Spatial", layer = "counts")), c("x", "y")]
  coords
})
sc_count <- LayerData(sce, assay = "RNA", layer = "counts")
sce$ident <- sce$anno
metadata <- sce@meta.data
sc_meta <- metadata %>% dplyr::select(orig.ident, ident)
sc_meta$cellID <- rownames(sc_meta)
CARD.list <- list()
for (i in seq_along(spatial_count.list)) {
  CARD_obj <- createCARDObject(
    sc_count = sc_count,
    sc_meta = sc_meta,
    spatial_count = spatial_count.list[[i]],
    spatial_location = spatial_loc.list[[i]],
    ct.varname = "ident",
    ct.select = unique(sc_meta$ident),
    sample.varname = "orig.ident",
    minCountGene = 500,
    minCountSpot = 20
  )
  CARD_obj <- CARD_deconvolution(CARD_obj)
  CARD.list[[i]] <- CARD_obj
}
deconv.list <- lapply(CARD.list, function(x) x@Proportion_CARD)
names(deconv.list) <- names(spe.list)
for (i in seq_along(spe.list)) {
  mat <- deconv.list[[i]]
  spe.list[[i]] <- AddMetaData(spe.list[[i]], metadata = mat)
}
spe <- merge(spe.list[[1]], y = spe.list[-1])
mat <- spe@meta.data[, c("AT1.like","AT2.like","B","Endothelial","Epi.like",
                         "Mast","Myeloid","Plasma","Stroma","T_NK")]
spe$card <- colnames(mat)[max.col(as.matrix(mat), ties.method = "first")]
table(spe$card)
card_matrix <- mat
cell_sums <- rowSums(card_matrix)
saveRDS(spe,paste0(dir,"merge2.rds"))
saveRDS(card_matrix,paste0(dir,"matrix2.rds"))
card_annotation <- spe$card
saveRDS(card_annotation, paste0(dir,"/card_annotation.rds"))

####rctd-based deconvolution####
spe <- readRDS(paste0(dir,"merge.rds"))
ref_counts <- LayerData(sce, assay = "RNA", layer = "counts")
cell_types <- factor(sce$anno)
reference <- Reference(
  counts = ref_counts,
  cell_types = cell_types
)
spe.list <- SplitObject(spe, split.by = "sample")
spatial_count.list <- lapply(spe.list, function(x) {
  LayerData(x, assay = "Spatial", layer = "counts")
})
spatial_loc.list <- lapply(spe.list, function(x) {
  coords <- GetTissueCoordinates(x, image = names(x@images)[1])
  coords <- coords[colnames(LayerData(x, assay = "Spatial", layer = "counts")), c("x", "y")]
  coords
})
spatialRNA.list <- lapply(1:length(spatial_count.list), function(i) {
  SpatialRNA(
    coords = spatial_loc.list[[i]],
    counts = spatial_count.list[[i]]
  )
})
rctd.list <- lapply(spatialRNA.list, function(spatialRNA_i) {
  myRCTD <- create.RCTD(
    spatialRNA = spatialRNA_i,
    reference = reference,
    max_cores = 4
  )
  myRCTD <- run.RCTD(myRCTD, doublet_mode = "full")
  return(myRCTD)
})
rctd_weights.list <- lapply(rctd.list, function(r) {
  r@results$weights
})
names(rctd.list) <- names(spe.list)
rctd_main.list <- lapply(rctd.list, function(r) {
  w <- r@results$weights
  max_type <- colnames(w)[apply(w, 1, which.max)]
  data.frame(
    cell = rownames(w),
    rctd_celltype = max_type
  )
})
names(rctd_main.list) <- names(rctd.list)
for (s in names(spe.list)) {
  spe.list[[s]]$rctd <- rctd_main.list[[s]][
    match(colnames(spe.list[[s]]), rctd_main.list[[s]]$cell),
    "rctd_celltype"
  ]
}
spe <- merge(spe.list[[1]], y = spe.list[-1])
table(spe$rctd)
rctd_weights_all <- do.call(rbind,rctd_weights.list)
rctd_matrix=as.matrix(rctd_weights_all)
cell_sums <- rowSums(rctd_matrix)
saveRDS(spe,paste0(dir,"merge3.rds"))
saveRDS(rcrd_matrix,paste0(dir,"matrix3.rds"))
rctd_annotation <- spe$rctd
saveRDS(rctd_annotation, paste0(dir,"/rctd_annotation.rds"))

####spotlight-based annotation####
spe <- readRDS(paste0(dir,"merge.rds"))
spe.list <- SplitObject(spe, split.by = "sample")
sce_sce <- as.SingleCellExperiment(sce)
sce_sce$celltype <- as.factor(sce$anno)
seu <- sce
seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000)
Idents(seu) <- seu$anno
markers <- FindAllMarkers(
  seu,
  only.pos = TRUE,
  logfc.threshold = 0.5,
  min.pct = 0.25,
  min.diff.pct = 0.2
)
mgs_df <- markers %>%
  group_by(cluster) %>%
  slice_max(n = 50, order_by = avg_log2FC) %>%
  ungroup() %>%
  dplyr::select(gene, cluster, avg_log2FC) %>%
  dplyr::rename(weight = avg_log2FC)
hvg <- VariableFeatures(seu)
spotlight_results <- vector("list", length(spe.list))
for (i in seq_along(spe.list)) {
  current_spe <- spe.list[[i]]
  counts_matrix <- LayerData(current_spe, assay = "Spatial", layer = "counts")
  coords <- GetTissueCoordinates(current_spe, image = names(current_spe@images)[1])
  coords <- coords[colnames(counts_matrix), c("x", "y")]
  common_genes <- intersect(rownames(sce_sce), rownames(counts_matrix))
  sce_subset <- sce_sce[common_genes, ]
  counts_subset <- counts_matrix[common_genes, ]
  mgs_subset <- mgs_df[mgs_df$gene %in% common_genes, ]
  hvg_subset <- intersect(hvg, common_genes)
  res <- SPOTlight(
    x = sce_subset,
    y = counts_subset,
    groups = sce_subset$celltype,
    mgs = mgs_subset,
    hvg = hvg_subset,
    weight_id = "weight",
    group_id = "cluster",
    gene_id = "gene",
    verbose = TRUE
  )
  deconv_matrix <- res$mat
  predicted_celltype <- colnames(deconv_matrix)[max.col(deconv_matrix, ties.method = "first")]
  spe.list[[i]]$spotlight <- predicted_celltype
  spotlight_results[[i]] <- list(
    deconv_matrix = res$mat,
    nmf_model = res$NMF,
    coords = coords,
    sample_id = i
  )
  rm(current_spe, counts_matrix, coords, sce_subset, 
     counts_subset, mgs_subset, hvg_subset, res, 
     deconv_matrix, predicted_celltype)
  gc()  
}
spe <- merge(spe.list[[1]], y = spe.list[-1])
table(spe$spotlight)
spotlight_scores_list <- lapply(spotlight_results, function(x) {
  x$deconv_matrix
})
spotlight_matrix <- do.call(rbind, spotlight_scores_list)
cell_sums <- rowSums(spotlight_matrix)
saveRDS(spe,paste0(dir,"/merge4.rds"))
saveRDS(spotlight_matrix,paste0(dir,"/matrix4.rds"))
spotlight_annotation <- spe$spotlight
saveRDS(spotlight_annotation,paste0(dir,"/spotlight_annotation.rds"))

####data integration####
spe=readRDS("~/Documents/R/benchmark/merge.rds")
seurat_matrix <- readRDS(paste0(dir, "/matrix1.rds"))
card_matrix <- readRDS(paste0(dir, "/matrix2.rds"))
rctd_matrix <- readRDS(paste0(dir, "/matrix3.rds"))
spotlight_matrix <- readRDS(paste0(dir, "/matrix4.rds"))
seurat <- readRDS(paste0(dir, "seurat_trans_annotation.rds"))
card <- readRDS(paste0(dir, "card_annotation.rds"))
rctd <- readRDS(paste0(dir, "rctd_annotation.rds"))
spotlight <- readRDS(paste0(dir, "spotlight_annotation.rds"))
spe$seurat=seurat
spe$card=card
spe$rctd=rctd
spe$spotlight=spotlight
all_matrix <- list(
  rctd_matrix = rctd_matrix,
  card_matrix = card_matrix,
  seurat_matrix = seurat_matrix,
  spotlight_matrix = spotlight_matrix
)
saveRDS(spe, paste0(dir, "/merge_anno.rds"))
saveRDS(all_matrix, paste0(dir, "/all_cell_type_matrix.rds"))

####visualization####
spe <- readRDS(paste0(dir,"merge_anno.rds"))
spe$seurat_trans=spe$seurat
spe$seurat_trans <- gsub("-", ".", spe$seurat_trans)
spe$card <- gsub("-", ".", spe$card)
spe$rctd <- gsub("-", ".", spe$rctd)
spe$spotlight <- gsub("-", ".", spe$spotlight)
keep_cells <- !(is.na(spe$seurat_trans) | 
                  is.na(spe$card) | 
                  is.na(spe$rctd) | 
                  is.na(spe$spotlight))
spe <- spe[, keep_cells]
custom_colors <- c(
  "AT1.like" = "#ce6d87",
  "AT2.like" = "#a08582", 
  "B" = "#a48441",
  "Endothelial" = "#869846",
  "Epi.like" = "#465439",
  "Mast" = "#55A372",
  "Myeloid" = "#6191c2",
  "Plasma" = "#404f69",
  "Stroma" = "#9b79c0",
  "T_NK" = "#6a366d"
)  
dir_result <- "~/Documents/R/benchmark/result/"
p=SpatialDimPlot(spe, group.by = "seurat_trans", images = "P1_LUAD",  pt.size.factor = 2, 
               image.alpha = 0.5,cols = custom_colors)
ggsave(file.path(dir_result , "spatialdim_seurat.pdf"), p, width = 10, height = 8,dpi = 600)
p=SpatialDimPlot(spe, group.by = "card", images = "P1_LUAD",  pt.size.factor = 2, 
               image.alpha = 0.5,cols = custom_colors)
ggsave(file.path(dir_result , "spatialdim_card.pdf"), p, width = 10, height = 8,dpi = 600)
p=SpatialDimPlot(spe, group.by = "rctd", images = "P1_LUAD",  pt.size.factor = 2, 
               image.alpha = 0.5,cols = custom_colors)
ggsave(file.path(dir_result , "spatialdim_rctd.pdf"), p, width = 10, height = 8,dpi = 600)
p=SpatialDimPlot(spe, group.by = "spotlight", images = "P1_LUAD",  pt.size.factor = 2, 
               image.alpha = 0.5,cols = custom_colors)
ggsave(file.path(dir_result , "spatialdim_spotlight.pdf"), p, width = 10, height = 8,dpi = 600)

####celltype composition####
group.by_vars <- c("seurat_trans", "card", "rctd", "spotlight")
all_celltypes <- names(custom_colors)
all_results <- list()
for(group in group.by_vars) {
  spe@meta.data[[group]] <- factor(spe@meta.data[[group]], levels = all_celltypes)
  result <- as.data.frame(table(spe$sample, spe@meta.data[[group]]))
  colnames(result) <- c("Sample", "CellType", "Count")
  result <- result[result$Count > 0, ]
  result <- result %>%
    group_by(Sample) %>%
    mutate(Percentage = Count / sum(Count) * 100)
  result$Method <- group  
  all_results[[group]] <- result
}
df_all <- bind_rows(all_results)
df_all$CellType <- factor(df_all$CellType, levels = all_celltypes)
bar_data <- df_all %>%
  group_by(Method, CellType) %>%
  summarise(
    MeanPct = mean(Percentage, na.rm = TRUE),
    SDPct = sd(Percentage, na.rm = TRUE),
    .groups = "drop"
  )
method_order <- c("seurat_trans", "card", "rctd", "spotlight")
bar_data$Method <- factor(bar_data$Method, levels = method_order)
df_all$Method <- factor(df_all$Method, levels = method_order)
stat_results <- df_all %>%
  group_by(CellType) %>%
  pairwise_t_test(
    Percentage ~ Method,
    p.adjust.method = "bonferroni"  
  ) %>%
  add_xy_position(x = "CellType")  
stat_results <- stat_results %>%
  group_by(CellType) %>%
  mutate(y.position = c(49,51,53,55, 57, 59)) %>%
  ungroup()
stat_results_sig <- stat_results %>%
  filter(p.adj < 0.05)
stat_results <- stat_results %>%
  mutate(
    significance = case_when(
      p.adj < 0.001 ~ "***",
      p.adj < 0.01 ~ "**",
      p.adj < 0.05 ~ "*",
      TRUE ~ "ns"
    )
  ) %>%
  filter(significance != "ns")  
p=ggplot() +
  geom_col(data = bar_data, 
           aes(x = CellType, y = MeanPct, fill = Method),
           position = position_dodge(width = 0.9), 
           alpha = 0.6) +
  geom_errorbar(data = bar_data,
                aes(x = CellType, ymin = MeanPct - SDPct, ymax = MeanPct + SDPct,
                    group = Method),
                position = position_dodge(width = 0.9),
                width = 0.2, size = 0.5, color = "black") +
  #geom_point(data = df_all,
   #          aes(x = CellType, y = Percentage, group = Method),
    #         position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.9),
     #        size = 1, alpha = 0.4, color = "black") +
  stat_pvalue_manual(
    stat_results,
    label = "significance",
    tip.length = 0,
    size = 4,
    bracket.size = 0.3,
    hide.ns = TRUE,
    vjust = 0.9
  ) +
  scale_fill_manual(values = c("seurat_trans" = "#F8766D", 
                               "card" = "#7CAE00",
                               "rctd" = "#00BFC4", 
                               "spotlight" = "#C77CFF"),
                    name = "Method") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.3))) +ylim(-7,60)+
  theme_minimal()+
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 12, face = "bold"),
    legend.position = "bottom",
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 9),
    panel.grid = element_blank(),  
    panel.border = element_rect(fill = NA, color = "black", linewidth = 0.5),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14)
  )
ggsave(file.path(dir_result , "cellcompositionall.pdf"), p, width = 10, height = 8,dpi = 600)

####cor heatmap####
all_matrix=readRDS(paste0(dir,"/all_cell_type_matrix.rds"))
dataframes <- lapply(all_matrix, as.data.frame)
common_genes <- Reduce(intersect, lapply(dataframes, rownames))
dataframes <- lapply(dataframes, function(df) df[common_genes, ])
dataframes <- lapply(dataframes, function(df) {
  colnames(df) <- gsub("-", ".", colnames(df))
  df
})
methods <- list(
  seurat = dataframes$seurat_matrix,
  rctd = dataframes$rctd_matrix,
  card = dataframes$card_matrix,
  spotlight = dataframes$spotlight_matrix
)
celltypes=c("AT1.like", "AT2.like", "B", "Endothelial", "Epi.like", 
            "Mast", "Myeloid", "Plasma", "Stroma", "T_NK")
res_list <- list()
for (ct in celltypes) {
  for (m1 in names(methods)) {
    for (m2 in names(methods)) {
      x <- methods[[m1]][, ct]
      y <- methods[[m2]][, ct]
      cor_val <- cor(x, y, use = "complete.obs")
      res_list[[length(res_list) + 1]] <- data.frame(
        celltype = ct,
        method1 = m1,
        method2 = m2,
        cor = cor_val
      )
    }
  }
}
res_df <- do.call(rbind, res_list)
order_methods <- c("seurat",  "card","rctd", "spotlight")
res_df$method1 <- factor(res_df$method1, levels = order_methods)
res_df$method2 <- factor(res_df$method2, levels = order_methods)
p=ggplot(res_df, aes(x = method1, y = method2, fill = cor)) +
  geom_tile() +
  facet_wrap(~celltype, nrow = 2) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  theme_minimal() +
  coord_equal()
ggsave(file.path(dir_result, "corheatmap.pdf"), p, width = 10, height = 8,dpi = 600)

####rmse heatmap####
all_matrix=readRDS(paste0(dir,"/all_cell_type_matrix.rds"))
dataframes <- lapply(all_matrix, as.data.frame)
common_cells <- Reduce(intersect, lapply(dataframes, rownames))
dataframes <- lapply(dataframes, function(df) df[common_cells, ])
dataframes <- lapply(dataframes, function(df) {
  colnames(df) <- gsub("-", ".", colnames(df))
  df
})
seurat_mat <- dataframes$seurat_matrix
rctd_mat <- dataframes$rctd_matrix
card_mat <- dataframes$card_matrix
spotlight_mat <- dataframes$spotlight_matrix
common_cells_final <- Reduce(intersect, list(
  rownames(seurat_mat),
  rownames(rctd_mat),
  rownames(card_mat),
  rownames(spotlight_mat)
))
seurat_mat <- seurat_mat[common_cells_final, ]
rctd_mat <- rctd_mat[common_cells_final, ]
card_mat <- card_mat[common_cells_final, ]
spotlight_mat <- spotlight_mat[common_cells_final, ]
celltypes <- colnames(seurat_mat)  
coords_list <- lapply(1:length(spe@images), function(i) {
  img <- names(spe@images)[i]
  coords <- GetTissueCoordinates(spe, image = img)
  cells <- Cells(spe)[spe$sample == img]
  coords <- coords[cells, c("x", "y"), drop = FALSE]
  coords$cell <- rownames(coords)
  coords$sample <- img
  coords
})
coords_df <- do.call(rbind, coords_list)
common_all <- Reduce(intersect, list(
  rownames(coords_df),
  rownames(seurat_mat),
  rownames(rctd_mat),
  rownames(card_mat),
  rownames(spotlight_mat)
))
coords_df <- coords_df[common_all, ]
seurat_mat <- seurat_mat[common_all, ]
rctd_mat <- rctd_mat[common_all, ]
card_mat <- card_mat[common_all, ]
spotlight_mat <- spotlight_mat[common_all, ]
df_all <- data.frame(
  x = coords_df$x,
  y = coords_df$y,
  sample = coords_df$sample,
  row.names = rownames(coords_df)
)
for (ct in celltypes) {
  df_all[[paste0("seurat_", ct)]] <- seurat_mat[, ct]
  df_all[[paste0("rctd_", ct)]] <- rctd_mat[, ct]
  df_all[[paste0("card_", ct)]] <- card_mat[, ct]
  df_all[[paste0("spotlight_", ct)]] <- spotlight_mat[, ct]
}
compute_disagreement <- function(ct) {
  mat <- cbind(
    df_all[[paste0("seurat_", ct)]],
    df_all[[paste0("card_", ct)]],
    df_all[[paste0("rctd_", ct)]],
    df_all[[paste0("spotlight_", ct)]]
  )
  apply(mat, 1, var, na.rm = TRUE)
}
rmse <- function(x, y) sqrt(mean((x - y)^2, na.rm = TRUE))
methods <- c("seurat", "rctd", "card", "spotlight")
order_methods <- c("seurat", "card","rctd", "spotlight")
rmse_list <- list()
for (ct in celltypes) {
  rmse_matrix <- matrix(NA, nrow = length(methods), ncol = length(methods),
                        dimnames = list(methods, methods))
  for (i in 1:length(methods)) {
    for (j in 1:length(methods)) {
      rmse_val <- rmse(df_all[[paste0(methods[i], "_", ct)]], 
                       df_all[[paste0(methods[j], "_", ct)]])
      rmse_matrix[i, j] <- rmse_val
    }
  }
  rmse_df_ct <- as.data.frame(as.table(rmse_matrix))
  colnames(rmse_df_ct) <- c("method1", "method2", "rmse")
  rmse_df_ct$celltype <- ct
  rmse_list[[ct]] <- rmse_df_ct
}
rmse_df <- do.call(rbind, rmse_list)
rmse_df$method1 <- factor(rmse_df$method1, levels = order_methods)
rmse_df$method2 <- factor(rmse_df$method2, levels = order_methods)
p=ggplot(rmse_df, aes(method1, method2, fill = rmse)) +
  geom_tile() +
  facet_wrap(~ celltype,nrow=2) +
  scale_fill_viridis_c() +
  coord_equal() +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(dir_result, "rmseheatmap.pdf"), p, width = 10, height = 8,dpi = 600)

####Moran's I####
compute_moran <- function(df, value_col) {
  coords <- as.matrix(df[, c("x", "y")])
  nb <- knn2nb(knearneigh(coords, k = 6))
  lw <- nb2listw(nb)
  moran.test(df[[value_col]], lw)
}
moran_list <- list()
for (s in unique(df_all$sample)) {
  sub <- df_all[df_all$sample == s, ]
  for (ct in celltypes) {
    mat <- cbind(
      sub[[paste0("seurat_", ct)]],
      sub[[paste0("card_", ct)]],
      sub[[paste0("rctd_", ct)]],
      sub[[paste0("spotlight_", ct)]]
    )
    sub$dis <- apply(mat, 1, var, na.rm = TRUE)
    test <- compute_moran(sub, "dis")
    moran_list[[length(moran_list)+1]] <- data.frame(
      sample = s,
      celltype = ct,
      moran_I = test$estimate[1],
      p_value = test$p.value
    )
  }
}
moran_df <- do.call(rbind, moran_list)
moran_wide <- moran_df %>%
  pivot_wider(id_cols = sample, names_from = celltype, values_from = moran_I)
moran_matrix <- as.matrix(moran_wide[, -1])
rownames(moran_matrix) <- moran_wide$sample
col_fun <- colorRamp2(seq(min(moran_matrix, na.rm = TRUE), 
                          max(moran_matrix, na.rm = TRUE), 
                          length = 100),
                      viridis(100))
p=Heatmap(moran_matrix,
        name = "Moran's I",
        col = col_fun,
        cluster_rows = TRUE,     
        cluster_columns = F,  
        show_row_names =F,
        show_column_names = TRUE,
        row_names_gp = gpar(fontsize = 8),
        column_names_gp = gpar(fontsize = 8, rot = 45),
        row_dend_width = unit(2, "cm"),
        column_dend_height = unit(2, "cm"))
pdf(file.path(dir_result, "moranheatmap.pdf"), width = 12, height = 12)
draw(p)
dev.off()
compute_pairwise_moran <- function(df, ct, method1, method2) {
  col1 <- paste0(method1, "_", ct)
  col2 <- paste0(method2, "_", ct)
  df$diff <- abs(df[[col1]] - df[[col2]])
  coords <- as.matrix(df[, c("x", "y")])
  valid_idx <- !is.na(df$diff)
  coords_valid <- coords[valid_idx, ]
  diff_valid <- df$diff[valid_idx]
  if(length(diff_valid) < 3 | sd(diff_valid, na.rm = TRUE) == 0) {
    return(data.frame(moran_I = NA, p_value = NA))
  }
  nb <- knn2nb(knearneigh(coords_valid, k = min(6, nrow(coords_valid)-1)))
  lw <- nb2listw(nb, style = "W", zero.policy = TRUE)
  test <- moran.test(diff_valid, lw, zero.policy = TRUE)
  return(data.frame(moran_I = test$estimate[1], p_value = test$p.value))
}
methods <- c("seurat", "card", "rctd", "spotlight")
method_pretty <- c("Seurat", "CARD", "RCTD", "SPOTlight")
pairwise_results <- list()
for (s in unique(df_all$sample)) {
  sub <- df_all[df_all$sample == s, ]
  for (ct in celltypes) {
    for (i in 1:(length(methods)-1)) {
      for (j in (i+1):length(methods)) {
        result <- compute_pairwise_moran(sub, ct, methods[i], methods[j])
        pairwise_results[[length(pairwise_results)+1]] <- data.frame(
          sample = s,
          celltype = ct,
          method1 = method_pretty[i],
          method2 = method_pretty[j],
          moran_I = result$moran_I,
          p_value = result$p_value
        )
      }
    }
  }
}
pairwise_df <- do.call(rbind, pairwise_results)
pairwise_df$method_pair <- paste0(pairwise_df$method1, " vs ", pairwise_df$method2)
pairwise_df$significant <- ifelse(!is.na(pairwise_df$p_value) & pairwise_df$p_value < 0.05, "p < 0.05", "p >= 0.05")
pairwise_df$method_pair <- factor(pairwise_df$method_pair, 
                                  levels = c("Seurat vs CARD", "Seurat vs RCTD", "Seurat vs SPOTlight",
                                             "CARD vs RCTD", "CARD vs SPOTlight", "RCTD vs SPOTlight"))

summary_df <- pairwise_df %>%
  group_by(celltype, method_pair) %>%
  summarise(mean_moran = mean(moran_I, na.rm = TRUE),
            sd_moran = sd(moran_I, na.rm = TRUE),
            .groups = "drop")
p=ggplot() +
  geom_bar(data = summary_df, aes(x = method_pair, y = mean_moran), 
           stat = "identity", fill = "steelblue", alpha = 0.5) +
  geom_errorbar(data = summary_df, aes(x = method_pair, 
                                       ymin = mean_moran - sd_moran, ymax = mean_moran + sd_moran), 
                width = 0.2, size = 0.5) +
  geom_point(data = pairwise_df, aes(x = method_pair, y = moran_I), 
             size = 1, alpha = 0.4, position = position_jitter(width = 0.2)) +
  facet_wrap(~ celltype, ncol = 5) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        panel.background = element_blank(),
        strip.background = element_blank())
ggsave(file.path(dir_result, "moranpair.pdf"), p, width = 10, height = 8,dpi = 600)

####lisa####
tumor_types <- c("AT1.like","AT2.like","Epi.like")
methods <- c("seurat_trans","card","rctd","spotlight")
spe.list <- SplitObject(spe, split.by = "sample")
spatial_loc.list <- lapply(spe.list, function(x) {
  coords <- GetTissueCoordinates(x, image = names(x@images)[1])
  coords <- coords[colnames(LayerData(x, assay = "Spatial", layer = "counts")), c("x","y")]
  coords
})
compute_lisa <- function(df, value_col) {
  coords <- as.matrix(df[,c("x","y")])
  nb <- knn2nb(knearneigh(coords, k = 4))
  listw <- nb2listw(nb, style = "W")
  lisa <- localmoran(df[[value_col]], listw)
  df$Ii <- lisa[,1]
  df$Pvalue <- lisa[,5]
  m <- mean(df[[value_col]], na.rm = TRUE)
  df$LISA_cluster <- "NS"
  df$LISA_cluster[df[[value_col]] > m & df$Ii > 0 & df$Pvalue < 0.05] <- "High-High"
  df$LISA_cluster[df[[value_col]] < m & df$Ii > 0 & df$Pvalue < 0.05] <- "Low-Low"
  df$LISA_cluster[df$Ii < 0 & df$Pvalue < 0.05] <- "Outlier"
  df
}
all_results <- list()
for (m in methods) {
  tmp_list <- list()
  for (s in names(spe.list)) {
    obj <- spe.list[[s]]
    coords <- spatial_loc.list[[s]]
    df <- data.frame(
      cell = colnames(obj),
      x = coords[,1],
      y = coords[,2],
      label = obj[[m]][,1],
      stringsAsFactors = FALSE
    )
    df$value <- ifelse(df$label %in% tumor_types, 1, 0)
    df$sample <- s
    df <- compute_lisa(df, "value")
    df$method <- m
    tmp_list[[s]] <- df
  }
  all_results[[m]] <- bind_rows(tmp_list)
}
lisa_results <- bind_rows(all_results)
summary_df <- lisa_results %>%
  group_by(sample, method) %>%
  summarise(
    HH_ratio = mean(LISA_cluster == "High-High"),
    LL_ratio = mean(LISA_cluster == "Low-Low"),
    Outlier_ratio = mean(LISA_cluster == "Outlier"),
    .groups = "drop"
  )
summary_df$method <- factor(summary_df$method, 
                            levels = c("seurat_trans", "card", "rctd", "spotlight"))
method_colors <-c("seurat_trans" = "#F8766D", 
                  "card" = "#7CAE00",
                  "rctd" = "#00BFC4", 
                  "spotlight" = "#C77CFF")
compare_list <- list(
  c("seurat_trans", "card"),
  c("seurat_trans", "rctd"),
  c("seurat_trans", "spotlight"),
  c("card", "rctd"),
  c("card", "spotlight"),
  c("rctd", "spotlight")
)
sig_comparisons <- list()
sig_y_positions <- c()
y_positions <- c(0.42, 0.45, 0.48, 0.51, 0.54, 0.57)
for(i in 1:length(compare_list)) {
  group1 <- summary_df$HH_ratio[summary_df$method == compare_list[[i]][1]]
  group2 <- summary_df$HH_ratio[summary_df$method == compare_list[[i]][2]]
  test <- wilcox.test(group1, group2)
  p_val <- test$p.value
  if(p_val < 0.05) {
    sig_comparisons[[length(sig_comparisons) + 1]] <- compare_list[[i]]
    sig_y_positions <- c(sig_y_positions, y_positions[i])
  }
}
p=ggplot(summary_df, aes(x = method, y = HH_ratio, fill = method)) +
  geom_violin(alpha = 0.6, width = 0.7, trim = FALSE, color = "black", linewidth = 0.8) +
  geom_point(position = position_jitter(width = 0, height = 0),
             size = 2.5, alpha = 0.7, color = "black", shape = 16) +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 3, color = "black") +
  stat_summary(fun.data = mean_cl_normal, geom = "errorbar", 
               width = 0.15, linewidth = 0.7, color = "black") +
  scale_fill_manual(values = method_colors) +
  theme_classic(base_size = 13) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 30, hjust = 1, size = 11, color = "black"),
        axis.text.y = element_text(size = 11, color = "black"),
        axis.title.y = element_text(size = 12, face = "bold")) +
  stat_compare_means(comparisons = sig_comparisons,  
                     method = "wilcox.test",
                     label = "p.signif",
                     size = 4,
                     bracket.size = 0.3,vjust = 0.8,
                     tip.length = 0.0,
                     label.y = sig_y_positions) +
  coord_cartesian(ylim = c(-0.1, max(y_positions) * 1)) 
ggsave(file.path(dir_result, "lisa_tumor.pdf"), p, width = 5, height = 4,dpi = 600)

immu_types <- c("B","T_NK","Mast","Myeloid","Plasma")
methods <- c("seurat_trans","card","rctd","spotlight")
spe.list <- SplitObject(spe, split.by = "sample")
spatial_loc.list <- lapply(spe.list, function(x) {
  coords <- GetTissueCoordinates(x, image = names(x@images)[1])
  coords <- coords[colnames(LayerData(x, assay = "Spatial", layer = "counts")), c("x","y")]
  coords
})
compute_lisa <- function(df, value_col) {
  coords <- as.matrix(df[,c("x","y")])
  nb <- knn2nb(knearneigh(coords, k = 4))
  listw <- nb2listw(nb, style = "W")
  lisa <- localmoran(df[[value_col]], listw)
  df$Ii <- lisa[,1]
  df$Pvalue <- lisa[,5]
  m <- mean(df[[value_col]], na.rm = TRUE)
  df$LISA_cluster <- "NS"
  df$LISA_cluster[df[[value_col]] > m & df$Ii > 0 & df$Pvalue < 0.05] <- "High-High"
  df$LISA_cluster[df[[value_col]] < m & df$Ii > 0 & df$Pvalue < 0.05] <- "Low-Low"
  df$LISA_cluster[df$Ii < 0 & df$Pvalue < 0.05] <- "Outlier"
  df
}
all_results <- list()
for (m in methods) {
  tmp_list <- list()
  for (s in names(spe.list)) {
    obj <- spe.list[[s]]
    coords <- spatial_loc.list[[s]]
    df <- data.frame(
      cell = colnames(obj),
      x = coords[,1],
      y = coords[,2],
      label = obj[[m]][,1],
      stringsAsFactors = FALSE
    )
    df$value <- ifelse(df$label %in% immu_types, 1, 0)
    df$sample <- s
    df <- compute_lisa(df, "value")
    df$method <- m
    tmp_list[[s]] <- df
  }
  all_results[[m]] <- bind_rows(tmp_list)
}
lisa_results <- bind_rows(all_results)
summary_df <- lisa_results %>%
  group_by(sample, method) %>%
  summarise(
    HH_ratio = mean(LISA_cluster == "High-High"),
    LL_ratio = mean(LISA_cluster == "Low-Low"),
    Outlier_ratio = mean(LISA_cluster == "Outlier"),
    .groups = "drop"
  )
summary_df$method <- factor(summary_df$method, 
                            levels = c("seurat_trans", "card", "rctd", "spotlight"))
method_colors <- c("seurat_trans" = "#F8766D", 
                   "card" = "#7CAE00",
                   "rctd" = "#00BFC4", 
                   "spotlight" = "#C77CFF")
compare_list <- list(
  c("seurat_trans", "card"),
  c("seurat_trans", "rctd"),
  c("seurat_trans", "spotlight"),
  c("card", "rctd"),
  c("card", "spotlight"),
  c("rctd", "spotlight")
)
sig_comparisons <- list()
sig_y_positions <- c()
y_positions <- c(0.32, 0.35, 0.38, 0.31, 0.34, 0.37)
for(i in 1:length(compare_list)) {
  group1 <- summary_df$HH_ratio[summary_df$method == compare_list[[i]][1]]
  group2 <- summary_df$HH_ratio[summary_df$method == compare_list[[i]][2]]
  test <- wilcox.test(group1, group2)
  p_val <- test$p.value
  if(p_val < 0.05) {
    sig_comparisons[[length(sig_comparisons) + 1]] <- compare_list[[i]]
    sig_y_positions <- c(sig_y_positions, y_positions[i])
  }
}
p=ggplot(summary_df, aes(x = method, y = HH_ratio, fill = method)) +
  geom_violin(alpha = 0.6, width = 0.7, trim = FALSE, color = "black", linewidth = 0.8) +
  geom_point(position = position_jitter(width = 0, height = 0),
             size = 2.5, alpha = 0.7, color = "black", shape = 16) +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 3, color = "black") +
  stat_summary(fun.data = mean_cl_normal, geom = "errorbar", 
               width = 0.15, linewidth = 0.7, color = "black") +
  scale_fill_manual(values = method_colors) +
  theme_classic(base_size = 13) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 30, hjust = 1, size = 11, color = "black"),
        axis.text.y = element_text(size = 11, color = "black"),
        axis.title.y = element_text(size = 12, face = "bold")) +
  stat_compare_means(comparisons = sig_comparisons,  
                     method = "wilcox.test",
                     label = "p.signif",
                     size = 4,
                     bracket.size = 0.3,vjust = 0.8,
                     tip.length = 0.0,
                     label.y = sig_y_positions) +
  coord_cartesian(ylim = c(-0.1, max(y_positions) * 1)) 
ggsave(file.path(dir_result, "lisa_immu.pdf"), p, width = 5, height = 4,dpi = 600)

####permutation####
methods <- c("seurat_trans","card","rctd","spotlight")
tumor_types <- c("AT1.like","AT2.like","Epi.like")
method_colors <- c("seurat_trans" = "#F8766D", 
                   "card" = "#7CAE00",
                   "rctd" = "#00BFC4", 
                   "spotlight" = "#C77CFF")
get_coords <- function(obj) {
  img_name <- names(obj@images)[1]
  coords <- GetTissueCoordinates(obj, image = img_name)
  coords <- coords[match(colnames(obj), rownames(coords)), c("x","y")]
  coords
}
get_group <- function(x) {
  if (is.na(x)) return("Other")
  if (x %in% tumor_types) return("Tumor")
  if (x == "B") return("B")
  if (x == "T_NK") return("T_NK")
  if (x == "Mast") return("Mast")
  if (x == "Myeloid") return("Myeloid")
  if (x == "Plasma") return("Plasma")
  if (x == "Stroma") return("Stroma")
  if (x == "Endothelial") return("Endothelial")
  return("Other")
}
compute_neighborhood_enrichment <- function(df, n_perm = 50) {
  coords <- as.matrix(df[,c("x","y")])
  nb <- knn2nb(knearneigh(coords, k = 4))
  edges <- do.call(rbind, lapply(1:length(nb), function(i) {
    if(length(nb[[i]]) == 0) return(NULL)
    data.frame(from = i, to = nb[[i]])
  }))
  edges$from_type <- df$group[edges$from]
  edges$to_type   <- df$group[edges$to]
  obs <- table(edges$from_type, edges$to_type)
  perm_array <- array(0, dim = c(nrow(obs), ncol(obs), n_perm))
  rownames(perm_array) <- rownames(obs)
  colnames(perm_array) <- colnames(obs)
  for (p in 1:n_perm) {
    perm_label <- sample(df$group)
    perm_from <- perm_label[edges$from]
    perm_to   <- perm_label[edges$to]
    perm_array[,,p] <- table(
      factor(perm_from, levels=rownames(obs)),
      factor(perm_to, levels=colnames(obs))
    )
  }
  perm_mean <- apply(perm_array, c(1,2), mean)
  perm_sd   <- apply(perm_array, c(1,2), sd)
  z <- (obs - perm_mean) / (perm_sd + 1e-6)
  z
}
all_results <- list()
for (m in methods) {
  sample_list <- list()
  for (s in names(spe.list)) {
    obj <- spe.list[[s]]
    coords <- get_coords(obj)
    df <- data.frame(
      x = coords[,1],
      y = coords[,2],
      label = obj[[m]][,1],
      stringsAsFactors = FALSE
    )
    df <- df[!is.na(df$label), ]
    df$group <- sapply(df$label, get_group)
    zmat <- compute_neighborhood_enrichment(df, n_perm = 50)
    sample_list[[s]] <- zmat
  }
  all_results[[m]] <- sample_list
}
extract_pairs <- function(zmat) {
  required <- c("Tumor","B","T_NK","Mast","Myeloid","Plasma","Stroma","Endothelial")
  if (!all(required %in% rownames(zmat))) return(NULL)
  data.frame(
    Tumor_B = zmat["Tumor","B"],
    Tumor_TNK = zmat["Tumor","T_NK"],
    Tumor_Mast = zmat["Tumor","Mast"],
    Tumor_Myeloid = zmat["Tumor","Myeloid"],
    Tumor_Plasma = zmat["Tumor","Plasma"],
    Tumor_Stroma = zmat["Tumor","Stroma"],
    Tumor_Endothelial = zmat["Tumor","Endothelial"]
  )
}
summary_list <- list()
for (m in methods) {
  tmp <- lapply(all_results[[m]], extract_pairs)
  tmp <- bind_rows(tmp)
  tmp$method <- m
  summary_list[[m]] <- tmp
}
summary_df <- bind_rows(summary_list)
pairs_to_plot <- c("Tumor_Stroma", "Tumor_Endothelial", "Tumor_Myeloid", 
                   "Tumor_TNK", "Tumor_B", "Tumor_Plasma", "Tumor_Mast")
plot_list <- list()
for (pair in pairs_to_plot) {
  plot_df <- summary_df[, c("method", pair)]
  colnames(plot_df) <- c("method", "Z_score")
  plot_df <- plot_df[!is.na(plot_df$Z_score), ]
  method_levels <- c("seurat_trans", "card", "rctd", "spotlight")
  plot_df$method <- factor(plot_df$method, levels = method_levels)
  p_values <- list()
  for (comp in list(
    c("seurat_trans", "card"),
    c("seurat_trans", "rctd"),
    c("seurat_trans", "spotlight"),
    c("card", "rctd"),
    c("card", "spotlight"),
    c("rctd", "spotlight")
  )) {
    group1 <- plot_df$Z_score[plot_df$method == comp[1]]
    group2 <- plot_df$Z_score[plot_df$method == comp[2]]
    if (length(group1) >= 2 & length(group2) >= 2) {
      p_val <- wilcox.test(group1, group2)$p.value
      if (p_val < 0.05) {
        p_values[[length(p_values) + 1]] <- list(comp = comp, p = p_val)
      }
    }
  }
  sig_comparisons <- lapply(p_values, function(x) x$comp)
  if (length(sig_comparisons) > 0) {
    y_max <- max(plot_df$Z_score, na.rm = TRUE)
    y_min <- min(plot_df$Z_score, na.rm = TRUE)
    label_y_positions <- seq(y_max + 1, y_max + 2 * length(sig_comparisons), length.out = length(sig_comparisons))
    p <- ggplot(plot_df, aes(x = method, y = Z_score, fill = method)) +
      geom_violin(alpha = 0.6, width = 0.7, trim = FALSE, color = "black", linewidth = 0.8) +
      geom_point(position = position_jitter(width = 0, height = 0),
                 size = 2.5, alpha = 0.7, color = "black", shape = 16) +
      stat_summary(fun = mean, geom = "point", shape = 18, size = 3, color = "black") +
      stat_summary(fun.data = mean_cl_normal, geom = "errorbar", 
                   width = 0.15, linewidth = 0.7, color = "black") +
      scale_fill_manual(values = method_colors) +
      theme_classic(base_size = 13) +
      theme(legend.position = "none",
            axis.text.x = element_text(angle = 30, hjust = 1, size = 11, color = "black"),
            axis.text.y = element_text(size = 11, color = "black"),
            axis.title.y = element_text(size = 12, face = "bold")) +
      stat_compare_means(comparisons = sig_comparisons,  
                         method = "wilcox.test",
                         label = "p.signif",
                         size = 4,
                         bracket.size = 0.3,
                         vjust = 0.8,
                         tip.length = 0.0,
                         label.y = label_y_positions) +
      coord_cartesian(ylim = c(y_min - 0.1, max(label_y_positions) + 10)) +
      labs(title = gsub("_", "-", pair), y = "Z-score", x = NULL)
  } else {
    y_max <- max(plot_df$Z_score, na.rm = TRUE)
    y_min <- min(plot_df$Z_score, na.rm = TRUE)
    p <- ggplot(plot_df, aes(x = method, y = Z_score, fill = method)) +
      geom_violin(alpha = 0.6, width = 0.7, trim = FALSE, color = "black", linewidth = 0.8) +
      geom_point(position = position_jitter(width = 0, height = 0),
                 size = 2.5, alpha = 0.7, color = "black", shape = 16) +
      stat_summary(fun = mean, geom = "point", shape = 18, size = 3, color = "black") +
      stat_summary(fun.data = mean_cl_normal, geom = "errorbar", 
                   width = 0.15, linewidth = 0.7, color = "black") +
      scale_fill_manual(values = method_colors) +
      theme_classic(base_size = 13) +
      theme(legend.position = "none",
            axis.text.x = element_text(angle = 30, hjust = 1, size = 11, color = "black"),
            axis.text.y = element_text(size = 11, color = "black"),
            axis.title.y = element_text(size = 12, face = "bold")) +
      coord_cartesian(ylim = c(y_min - 0.1, y_max + 10)) +
      labs(title = gsub("_", "-", pair), y = "Z-score", x = NULL)
  }
  plot_list[[pair]] <- p
}
p<- ggarrange(plotlist = plot_list, ncol = 7, nrow = 1)
ggsave(file.path(dir_result, "neighbor.pdf"), p, width = 25, height = 4,dpi = 600)

####cell type score####
Tcell_genes <- c("CD3D","CD3E","CD3G","TRAC","IL7R","CCR7","LTB","TCF7")
Bcell_genes <- c("MS4A1","CD79A","CD79B","CD19","CD74","HLA-DRA")
Myeloid_genes <- c("LYZ","LST1","S100A8","S100A9","CD14","FCGR3A","CST3")
Endo_genes <- c("PECAM1","VWF","CLDN5","KDR","FLT1","ENG")
Stromal_genes <- c("COL1A1","COL1A2","DCN","LUM","COL3A1","FAP","PDGFRA")
Mast_genes <- c("TPSAB1","TPSB2","CPA3","KIT","HDC")
Plasma_genes <- c("MZB1","XBP1","PRDM1","SDC1","JCHAIN")
gene_list <- list(
  Tcell = Tcell_genes,
  Bcell = Bcell_genes,
  Myeloid = Myeloid_genes,
  Endo = Endo_genes,
  Stromal = Stromal_genes,
  Mast = Mast_genes,
  Plasma=Plasma_genes 
)
methods <- c("seurat_trans","card","rctd","spotlight")
DefaultAssay(spe) <- "Spatial"
spe <- NormalizeData(spe)
for (m in methods) {
  obj <- spe
  for (ct in names(gene_list)) {
    obj <- AddModuleScore(
      obj,
      features = list(gene_list[[ct]]),
      name = paste0(ct, "_", m, "_score"),
      assay = "Spatial"
    )
  }
  spe <- obj
}
methods <- c("seurat_trans", "card", "rctd", "spotlight")
score_to_annotation_map <- list(
  Tcell = "T_NK",
  Bcell = "B",
  Plasma = "Plasma",
  Myeloid = "Myeloid",
  Endo = "Endothelial",
  Stromal = "Stroma",
  Mast = "Mast"
)
all_data <- list()
for (method in methods) {
  celltype_col <- method
  df_list <- list()
  for (score_prefix in names(score_to_annotation_map)) {
    annotation_name <- score_to_annotation_map[[score_prefix]]
    score_col <- paste0(score_prefix, "_", method, "_score1")
    cells_of_interest <- spe@meta.data[[celltype_col]] == annotation_name
    cells_of_interest[is.na(cells_of_interest)] <- FALSE
    if (sum(cells_of_interest, na.rm = TRUE) > 0) {
      temp_df <- data.frame(
        CellType = score_prefix,
        Score = spe@meta.data[[score_col]][cells_of_interest],
        Method = method
      )
      df_list[[score_prefix]] <- temp_df
    }
  }
  if (length(df_list) > 0) {
    all_data[[method]] <- do.call(rbind, df_list)
  }
}
combined_data <- do.call(rbind, all_data)
combined_data$Method <- factor(combined_data$Method, levels = c("seurat_trans", "card", "rctd", "spotlight"))
method_comparisons <- list(
  c("seurat_trans", "card"),
  c("seurat_trans", "rctd"),
  c("seurat_trans", "spotlight"),
  c("card", "rctd"),
  c("card", "spotlight"),
  c("rctd", "spotlight")
)
mast_comparisons <- list(
  c("seurat_trans", "spotlight")
)
cell_types <- c("Tcell", "Bcell", "Plasma", "Myeloid", "Endo", "Stromal","Mast")
plot_list <- list()
for (ct in cell_types) {
  subset_data <- combined_data[combined_data$CellType == ct, ]
  if (length(unique(subset_data$Method)) >= 2) {
    p <- ggplot(subset_data, aes(x = Method, y = Score, fill = Method)) +
      geom_violin(trim = FALSE, alpha = 0.7) +
      geom_boxplot(width = 0.1, fill = "white", alpha = 0.5) +
      stat_compare_means(comparisons = method_comparisons,
                         method = "wilcox.test",
                         label = "p.signif",
                         size = 3,vjust =0.3,tip.length = 0,
                         bracket.size = 0.3) +
      theme_bw() +
      labs(title = ct, y = "Module Score", x = "") +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "none")
  } else {
    p <- ggplot(subset_data, aes(x = Method, y = Score, fill = Method)) +
      geom_violin(trim = FALSE, alpha = 0.7) +
      geom_boxplot(width = 0.1, fill = "white", alpha = 0.5) +
      theme_bw() +
      labs(title = ct, y = "Module Score", x = "") +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "none")
  }
  plot_list[[ct]] <- p
}

combined_plot <- wrap_plots(plot_list, ncol = 4)
print(combined_plot)
p=combined_plot <- wrap_plots(plot_list, ncol = 4)
ggsave(file.path(dir_result, "tmecellscore.pdf"), p, width = 8, height = 8,dpi = 600)

####tumor score####
methods <- c("seurat_trans", "card", "rctd", "spotlight")
gene_list <- list(
  Lung_epi= c("EPCAM","KRT7","KRT8","KRT18","KRT19","CLDN3","CLDN4","CDH1","KRT5"),
  Epi_proliferative = c("MKI67","TOP2A","PCNA","UBE2C","CDK1","CCNB1","CENPF","HMGB2"),
  EMT = c("VIM","FN1","SNAI1","SNAI2","ZEB1","ZEB2","TWIST1","CDH2","ITGA5"),
  MYC_program = c("MYC","E2F1","E2F2","E2F3","NPM1","RPL5","RPL11"),
  PI3K_AKT = c("PIK3CA","AKT1","MTOR","PTEN","TSC1","TSC2"),
  MAPK_EGFR = c("EGFR","KRAS","BRAF","MAPK1","MAPK3","FOS","JUN"),
  UPR_ER_stress = c("XBP1","HSPA5","DDIT3","ATF4")
)
DefaultAssay(spe) <- "Spatial"
for (m in methods) {
  obj <- spe
  for (ct in names(gene_list)) {
    obj <- AddModuleScore(
      obj,
      features = list(gene_list[[ct]]),
      name = paste0(ct, "_", m, "_score"),
      assay = "Spatial"
    )
  }
  spe <- obj
}
cell_types_of_interest <- c("AT1.like", "AT2.like", "Epi.like")
pathways <- names(gene_list)
all_data <- list()
for (method in methods) {
  for (pathway in pathways) {
    score_col <- paste0(pathway, "_", method, "_score1")
    cells_of_interest <- spe@meta.data[[method]] %in% cell_types_of_interest
    cells_of_interest[is.na(cells_of_interest)] <- FALSE
    if (sum(cells_of_interest, na.rm = TRUE) > 0) {
      temp_df <- data.frame(
        Pathway = pathway,
        Score = spe@meta.data[[score_col]][cells_of_interest],
        Method = method
      )
      all_data[[paste0(method, "_", pathway)]] <- temp_df
    }
  }
}
plot_data <- do.call(rbind, all_data)
plot_data$Method <- factor(plot_data$Method, levels = c("seurat_trans", "card", "rctd", "spotlight"))
method_comparisons <- list(
  c("seurat_trans", "card"),
  c("seurat_trans", "rctd"),
  c("seurat_trans", "spotlight"),
  c("card", "rctd"),
  c("card", "spotlight"),
  c("rctd", "spotlight")
)
plot_list <- list()
for (pathway in pathways) {
  subset_data <- plot_data[plot_data$Pathway == pathway, ]
  if (nrow(subset_data) > 0 && length(unique(subset_data$Method)) >= 2) {
    p <- ggplot(subset_data, aes(x = Method, y = Score, fill = Method)) +
      geom_violin(trim = FALSE, alpha = 0.7) +
      geom_boxplot(width = 0.1, fill = "white", alpha = 0.5) +
      stat_compare_means(comparisons = method_comparisons,
                         method = "wilcox.test",
                         label = "p.signif",
                         size = 3,vjust = 0.3,tip.length = 0,
                         bracket.size = 0.3) +
      theme_bw() +
      labs(title = pathway, y = "Module Score", x = "") +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "none")
  } else if (nrow(subset_data) > 0) {
    p <- ggplot(subset_data, aes(x = Method, y = Score, fill = Method)) +
      geom_violin(trim = FALSE, alpha = 0.7) +
      geom_boxplot(width = 0.1, fill = "white", alpha = 0.5) +
      theme_bw() +
      labs(title = pathway, y = "Module Score", x = "") +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "none")
  } else {
    next
  }
  
  plot_list[[pathway]] <- p
}
combined_plot <- wrap_plots(plot_list, ncol = 4)
p=combined_plot
ggsave(file.path(dir_result, "tumorpathscore.pdf"), p, width = 8, height = 8,dpi = 600)

####tme score####
gene_list <- list(
  Immune_escape = c("CD274","PDCD1","CTLA4","LAG3","TIGIT","HAVCR2","IDO1","ARG1","HLA-A","HLA-B","HLA-C","B2M","TAP1","TAP2"),
  ECM_remodeling = c("COL1A1","COL1A2","COL3A1","COL5A1","FN1","SPARC","LOX","LOXL2","MMP2","MMP9","MMP14","ITGA5","ITGB1"),
  Immune_infiltration = c("PTPRC","CD3D","CD3E","CD4","CD8A","MS4A1","CD79A","NKG7","GNLY","GZMB","PRF1","LYZ","CD14","LST1","S100A8","S100A9"),
  Fibroblast_activation = c("COL1A1","COL1A2","COL3A1","ACTA2","TAGLN","PDGFRA","PDGFRB","FAP","THY1","CXCL12"),
  Antigen_presentation = c("HLA-A","HLA-B","HLA-C","B2M","HLA-DRA","HLA-DRB1","TAP1","TAP2"),
  Angiogenesis = c("VEGFA","KDR","FLT1","ANGPT1","ANGPT2","TEK","PECAM1","VWF"),
  Hypoxia = c("HIF1A","EPAS1","VEGFA","LDHA","SLC2A1")
)
methods <- c("seurat_trans", "card", "rctd", "spotlight")
DefaultAssay(spe) <- "Spatial"
for (m in methods) {
  obj <- spe
  for (ct in names(gene_list)) {
    obj <- AddModuleScore(
      obj,
      features = list(gene_list[[ct]]),
      name = paste0(ct, "_", m, "_score"),
      assay = "Spatial"
    )
  }
  spe <- obj
}
cell_types_of_interest <- c( "B", "Endothelial",  "Mast", "Myeloid", "Plasma", "Stroma", "T_NK")
pathways <- names(gene_list)
all_data <- list()
for (method in methods) {
  for (pathway in pathways) {
    score_col <- paste0(pathway, "_", method, "_score1")
    cells_of_interest <- spe@meta.data[[method]] %in% cell_types_of_interest
    cells_of_interest[is.na(cells_of_interest)] <- FALSE
    if (sum(cells_of_interest, na.rm = TRUE) > 0) {
      temp_df <- data.frame(
        Pathway = pathway,
        Score = spe@meta.data[[score_col]][cells_of_interest],
        Method = method
      )
      all_data[[paste0(method, "_", pathway)]] <- temp_df
    }
  }
}
plot_data <- do.call(rbind, all_data)
plot_data$Method <- factor(plot_data$Method, levels = c("seurat_trans", "card", "rctd", "spotlight"))
method_comparisons <- list(
  c("seurat_trans", "card"),
  c("seurat_trans", "rctd"),
  c("seurat_trans", "spotlight"),
  c("card", "rctd"),
  c("card", "spotlight"),
  c("rctd", "spotlight")
)
plot_list <- list()
for (pathway in pathways) {
  subset_data <- plot_data[plot_data$Pathway == pathway, ]
  if (nrow(subset_data) > 0 && length(unique(subset_data$Method)) >= 2) {
    p <- ggplot(subset_data, aes(x = Method, y = Score, fill = Method)) +
      geom_violin(trim = FALSE, alpha = 0.7) +
      geom_boxplot(width = 0.1, fill = "white", alpha = 0.5) +
      stat_compare_means(comparisons = method_comparisons,
                         method = "wilcox.test",
                         label = "p.signif",
                         size = 3,vjust = 0.3,tip.length = 0,
                         bracket.size = 0.3) +
      theme_bw() +
      labs(title = pathway, y = "Module Score", x = "") +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "none")
  } else if (nrow(subset_data) > 0) {
    p <- ggplot(subset_data, aes(x = Method, y = Score, fill = Method)) +
      geom_violin(trim = FALSE, alpha = 0.7) +
      geom_boxplot(width = 0.1, fill = "white", alpha = 0.5) +
      theme_bw() +
      labs(title = pathway, y = "Module Score", x = "") +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "none")
  } else {
    next
  }
  plot_list[[pathway]] <- p
}
combined_plot <- wrap_plots(plot_list, ncol = 4)
print(combined_plot)
p=combined_plot
ggsave(file.path(dir_result, "tmepathscore.pdf"), p, width = 8, height = 8,dpi = 600)

####gradient####
tumor_types <- c("AT1.like", "AT2.like", "Epi.like")
methods <- c("seurat_trans", "card", "rctd", "spotlight")
get_coords <- function(spe) {
  coords <- GetTissueCoordinates(spe, image = names(spe@images)[1])
  coords <- coords[colnames(spe), c("x", "y")]
  coords
}
get_cell_distance <- function(spe, label_col, coords, cell_types) {
  labels <- spe@meta.data[[label_col]]
  target_cells <- colnames(spe)[labels %in% cell_types]
  if(length(target_cells) < 3) return(NULL)
  target_coords <- coords[target_cells, , drop = FALSE]
  nn <- get.knnx(as.matrix(target_coords), as.matrix(coords), k = 1)
  dist <- nn$nn.dist[, 1]
  names(dist) <- rownames(coords)
  dist
}
run_continuous_DE_fast <- function(spe, label_col, target_types, sample_size = 5000) {
  coords <- get_coords(spe)
  dist <- get_cell_distance(spe, label_col, coords, target_types)
  if(is.null(dist)) return(NULL)
  expr <- GetAssayData(spe, layer = "counts")
  keep <- !is.na(dist[colnames(spe)])
  if(sum(keep) < 50) return(NULL)
  expr <- expr[, keep, drop = FALSE]
  dist <- dist[colnames(expr)]
  celltype <- spe[[label_col]][colnames(expr), 1]
  if(ncol(expr) > sample_size) {
    set.seed(123)
    samp_idx <- sample(ncol(expr), sample_size)
    expr <- expr[, samp_idx]
    dist <- dist[samp_idx]
    celltype <- celltype[samp_idx]
  }
  genes <- rownames(expr)
  gene_present <- rowSums(expr > 0) >= 10
  genes <- genes[gene_present]
  if(length(genes) == 0) return(NULL)
  expr_mat <- as.matrix(expr[genes, , drop = FALSE])
  expr_log <- log1p(expr_mat)
  dist_vec <- dist
  celltype_factor <- as.factor(celltype)
  run_lm_fast <- function(gene_expr) {
    df <- data.frame(y = gene_expr, dist = dist_vec, celltype = celltype_factor)
    fit <- lm(y ~ dist + celltype, data = df)
    coef_sum <- summary(fit)$coefficients
    if("dist" %in% rownames(coef_sum)) {
      return(c(coef_sum["dist", "Estimate"], coef_sum["dist", "Pr(>|t|)"]))
    } else {
      return(c(NA, NA))
    }
  }
  lm_res <- t(apply(expr_log, 1, run_lm_fast))
  colnames(lm_res) <- c("beta", "pval")
  res <- data.frame(gene = genes, beta = lm_res[, "beta"], pval = lm_res[, "pval"])
  res <- res[!is.na(res$beta), ]
  if(nrow(res) == 0) return(NULL)
  res$FDR <- p.adjust(res$pval, method = "BH")
  res
}
all_celltype_res <- list()
for (method in methods) {
  cat("\n========================================\n")
  cat("Processing method:", method, "\n")
  cat("========================================\n")
  celltype_res <- list()
  all_celltypes <- unique(unlist(lapply(spe.list, function(spe) {
    unique(spe@meta.data[[method]])
  })))
  all_celltypes <- all_celltypes[!all_celltypes %in% tumor_types]
  cat("Found", length(all_celltypes), "cell types to process\n")
  for (ct in all_celltypes) {
    cat("\n  [", which(all_celltypes == ct), "/", length(all_celltypes), "] Cell type:", ct, "\n")
    start_time <- Sys.time()
    res_list <- list()
    for(j in seq_along(spe.list)) {
      spe <- spe.list[[j]]
      sample_name <- names(spe.list)[j]
      if(ct %in% unique(spe@meta.data[[method]])) {
        cat("    Processing sample:", sample_name, "... ")
        res <- run_continuous_DE_fast(spe, method, ct, sample_size = 3000)
        if(!is.null(res) && nrow(res) > 0) {
          res_list[[sample_name]] <- res
          cat(nrow(res), "genes\n")
        } else {
          cat("no significant results\n")
        }
      }
      rm(spe)
      gc()
    }
    if(length(res_list) > 0) {
      combined <- bind_rows(res_list, .id = "sample")
      final <- combined %>%
        group_by(gene) %>%
        summarise(
          mean_beta = mean(beta, na.rm = TRUE),
          mean_p = mean(pval, na.rm = TRUE),
          freq = sum(FDR < 0.05, na.rm = TRUE),
          n_samples = n()
        ) %>%
        mutate(FDR = p.adjust(mean_p, method = "BH")) %>%
        arrange(FDR)
      celltype_res[[ct]] <- final
      cat("    Saved", nrow(final), "genes for", ct, "\n")
    } else {
      cat("    No results for", ct, "\n")
    }
    end_time <- Sys.time()
    cat("    Time elapsed:", round(difftime(end_time, start_time, units = "mins"), 2), "minutes\n")
    rm(res_list, combined, final)
    gc()
  }
  all_celltype_res[[method]] <- celltype_res
  saveRDS(all_celltype_res, paste0("results_", method, ".rds"))
  cat("\nSaved results for", method, "to results_", method, ".rds\n")
  rm(celltype_res, all_celltypes)
  gc()
}
methods <- c("seurat_trans","card","rctd","spotlight")
get_significant_genes <- function(res_list, method_name, celltype, fdr_cutoff = 0.05) {
  if(!method_name %in% names(res_list)) return(character(0))
  if(!celltype %in% names(res_list[[method_name]])) return(character(0))
  df <- res_list[[method_name]][[celltype]]
  if(is.null(df) || nrow(df) == 0) return(character(0))
  df$gene[df$FDR < fdr_cutoff & !is.na(df$FDR)]
}
all_celltypes <- c("Myeloid", "Endothelial", "Stroma", "T_NK", "Plasma", "Mast", "B")
sig_genes_2methods <- c()
for(ct in all_celltypes) {
  for(i in 1:(length(methods)-1)) {
    for(j in (i+1):length(methods)) {
      genes_i <- get_significant_genes(all_celltype_res, methods[i], ct)
      genes_j <- get_significant_genes(all_celltype_res, methods[j], ct)
      common <- intersect(genes_i, genes_j)
      if(length(common) > 0) {
        sig_genes_2methods <- union(sig_genes_2methods, common)
      }
    }
  }
}
beta_matrix <- matrix(NA, 
                      nrow = length(sig_genes_2methods),
                      ncol = length(methods) * length(all_celltypes))
neglogp_matrix <- matrix(NA,
                         nrow = length(sig_genes_2methods),
                         ncol = length(methods) * length(all_celltypes))
col_idx <- 1
col_names <- c()
for(ct in all_celltypes) {
  for(m in methods) {
    col_name <- paste0(ct, "_", m)
    col_names <- c(col_names, col_name)
    if(m %in% names(all_celltype_res) && 
       ct %in% names(all_celltype_res[[m]]) &&
       !is.null(all_celltype_res[[m]][[ct]])) {
      df <- all_celltype_res[[m]][[ct]]
      if(nrow(df) > 0) {
        for(i in seq_along(sig_genes_2methods)) {
          gene <- sig_genes_2methods[i]
          idx <- which(df$gene == gene)
          if(length(idx) > 0) {
            if("mean_beta" %in% colnames(df)) {
              beta_matrix[i, col_idx] <- df$mean_beta[idx[1]]
            }
            if("mean_p" %in% colnames(df)) {
              neglogp_matrix[i, col_idx] <- -log10(df$mean_p[idx[1]])
            }
          }
        }
      }
    }
    col_idx <- col_idx + 1
  }
}
rownames(beta_matrix) <- sig_genes_2methods
colnames(beta_matrix) <- col_names
rownames(neglogp_matrix) <- sig_genes_2methods
colnames(neglogp_matrix) <- col_names
keep_rows <- apply(!is.na(beta_matrix), 1, any)
beta_matrix <- beta_matrix[keep_rows, ]
neglogp_matrix <- neglogp_matrix[keep_rows, ]
beta_clean <- beta_matrix
beta_clean[is.na(beta_clean)] <- 0
row_order <- order(apply(beta_clean, 1, function(x) max(abs(x), na.rm = TRUE)), decreasing = TRUE)
beta_matrix <- beta_matrix[row_order, ]
neglogp_matrix <- neglogp_matrix[row_order, ]
col_split <- rep(all_celltypes, each = length(methods))
col_fun_beta <- colorRamp2(c(-0.0001, -0.00005, 0, 0.00005, 0.0001), 
                           c("darkblue", "blue", "white", "red", "darkred"))
col_fun_p <- colorRamp2(c(0, 5, 10, 15, 20, 30, 40),
                        c("pink", "lightcoral", "red", "firebrick", "darkred", "darkred", "darkred"))
pdf(file.path(dir_result, "heatmap_beta.pdf"), width = 20, height = 16)
Heatmap(beta_matrix,
        name = "beta",
        col = col_fun_beta,
        na_col = "gray90",
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        show_row_names = TRUE,
        show_column_names = TRUE,
        row_names_gp = gpar(fontsize = 8),
        column_names_gp = gpar(fontsize = 8, rot = 45),
        column_split = col_split,
        row_title = "Genes",
        column_title = "Beta Value",
        row_names_max_width = unit(10, "cm"),
        heatmap_legend_param = list(title = "beta", direction = "horizontal"))
dev.off()
pdf(file.path(dir_result, "heatmap_neglog10p.pdf"), width = 20, height = 16)
Heatmap(neglogp_matrix,
        name = "-log10(p)",
        col = col_fun_p,
        na_col = "gray90",
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        show_row_names = TRUE,
        show_column_names = TRUE,
        row_names_gp = gpar(fontsize = 8),
        column_names_gp = gpar(fontsize = 8, rot = 45),
        column_split = col_split,
        row_title = "Genes",
        column_title = "-log10(p)",
        row_names_max_width = unit(10, "cm"),
        heatmap_legend_param = list(title = "-log10(p)", direction = "horizontal"))
dev.off()
saveRDS(all_celltype_res, file.path(dir_result, "all_celltype_results.rds"))
