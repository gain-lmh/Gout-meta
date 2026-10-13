library(Seurat)
library(clustree)
library(cowplot)
library(data.table)
library(ggplot2)
library(patchwork)
library(stringr)
library(qs)
library(Matrix)
library(dplyr)

####批量读取Seurat对象
setwd("/data3/home/hu/Tongfeng/data")
dir()

samples <- c("A1","A2","A3","Con1","Con2","Con3","GR1","GR2","GR3","GS1","GS2","GS3")
data_dir <- "/data3/home/hu/Tongfeng/data"  # 你的整理后的文件夹所在目录

seu_list <- lapply(samples, function(s){
  # Read10X 读取10X矩阵
  data_path <- file.path(data_dir, s)
  counts <- Read10X(data.dir = data_path)
  
  # 创建 Seurat 对象
  seu <- CreateSeuratObject(counts = counts, project = s)
  
  # 添加样本信息
  seu$sample <- s
  
  return(seu)
})

names(seu_list) <- samples

seu_all <- merge(seu_list[[1]], y = seu_list[-1], add.cell.ids = samples)

seu_all$group <- case_when(
  grepl("Con", seu_all$sample) ~ "Control",
  grepl("GR", seu_all$sample)  ~ "Int gout remission",
  grepl("GS", seu_all$sample)  ~ "Adv gout remission",
  grepl("A", seu_all$sample)   ~ "Acute flare"
)

seu_all <- JoinLayers(seu_all)
table(seu_all$group)
table(seu_all$sample)

#####V5 单细胞标准分析流程
{
  ##Step 0、加载需要的R包
  library(Seurat)
  library(dplyr)
  library(patchwork)
  #library(mindr)
  library(Matrix)
  
  ##Step 1、数据准备
  sce.all <- seu_all
  
  #####查看一下数据状况
  as.data.frame(sce.all@assays$RNA$counts[1:10, 1:2])
  head(sce.all@meta.data, 10)
  summary(sce.all@meta.data$nCount_RNA)
  
  
  
  #Step2 ####质控QC
  {
    # Add number of genes per UMI for each cell to metadata
    sce.all$log10GenesPerUMI <- log10(sce.all$nFeature_RNA) / log10(sce.all$nCount_RNA)
    head(sce.all@meta.data)
    #### 线粒体比例
    # Compute percent mito ratio
    sce.all$mitoRatio <- PercentageFeatureSet(object = sce.all, pattern = "^mt-")
    sce.all$mitoRatio <- sce.all@meta.data$mitoRatio / 100
    summary(sce.all$mitoRatio)
    boxplot(sce.all$mitoRatio)
    
    ###质控前小提琴图
    VlnPlot(sce.all, 
            features = c("nFeature_RNA", "nCount_RNA", "log10GenesPerUMI", "mitoRatio"), 
            ncol = 4,
            pt.size = 0)
    
    ####质控过滤
    # Filter out low quality cells using selected thresholds - these will change with experiment
    filtered_seurat <- subset(x = sce.all, 
                              subset= (nCount_RNA >= 200) & 
                                (nFeature_RNA >= 150) & 
                                (nFeature_RNA <= 4000) & 
                                (log10GenesPerUMI > 0.80))
    paste0("质控过滤掉了", ncol(sce.all) - ncol(filtered_seurat), "个细胞")
    
    ####过滤基因
    
    ###质控后小提琴图
    
    # 经典框框
    col <- c("#866DA6",  #BAT
                   "#38A8D7",  #Colon
                   "#E8C1C7",  #Duodenum
                   "#CD3E39",  #gWAT
                   "#EEC684",  #Heart
                   "#BAD09E",  #Ileum
                   "#D77DA7",  #iWAT
                   "#CEC1D7",  #Jejunum
                   "#A7CCE5", #Kidney
                   "#B6773D",  #Liver
                   "#DDBE90",  #Lung
                   "#9BC35A",  #Muscle
                   "#E7AFC9",  #Stomach
                   "#67B1C9",  #
                   "#E1D572"   #
    )
    
    
    pdf("/data3/home/hu/Tongfeng/Result/Fig1/nFeature_RNA_vlnplot.pdf",height = 8,width = 12)
    VlnPlot(filtered_seurat, 
            features = c("nFeature_RNA", "nCount_RNA", "log10GenesPerUMI"), 
            ncol = 4,cols = col,
            pt.size = 0)
    dev.off()
  }
  
  #Step3  #标准化
  seurat_phase <- NormalizeData(filtered_seurat)
  seurat_phase  
  
  ##Step 4、评估细胞周期的影响
  # Load cell cycle markers
  s.genes <- cc.genes$s.genes
  g2m.genes <- cc.genes$g2m.genes
  
  library(biomaRt)
  # 将人类基因转换为小鼠基因
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl",host = "https://dec2021.archive.ensembl.org/")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl",host = "https://dec2021.archive.ensembl.org/")
  
  # 查询基因的同源性
  s.genes = getLDS(attributes = c("hgnc_symbol"), 
                   filters = "hgnc_symbol", 
                   values = s.genes, 
                   mart = human, 
                   attributesL = c("mgi_symbol"), 
                   martL = mouse, 
                   uniqueRows=T)
  
  g2m.genes = getLDS(attributes = c("hgnc_symbol"), 
                     filters = "hgnc_symbol", 
                     values = g2m.genes, 
                     mart = human, 
                     attributesL = c("mgi_symbol"), 
                     martL = mouse, 
                     uniqueRows=T)
  
  # Score cells for cell cycle
  seurat_phase <- CellCycleScoring(seurat_phase, 
                                   g2m.features = g2m.genes$MGI.symbol, 
                                   s.features = s.genes$MGI.symbol)
  
  # 现在的meta.data中多出了细胞周期评分“S.Score”和“G2M.Score”，以及推断的细胞所处的周期“Phase”
  head(seurat_phase@meta.data)
  
  # 查看一下细胞周期的分布情况
  table(seurat_phase$Phase)
  
  # ##Step 5、鉴定高变基因 Identify the most variable genes
  seurat_phase <- FindVariableFeatures(seurat_phase, 
                                       selection.method = "vst", # 默认值
                                       nfeatures = 2000) # 默认值
  
  # Scale the counts
  seurat_phase <- ScaleData(seurat_phase)
  Layers(seurat_phase)
  
  ###Step 6 降维 #Perform PCA。如果没有指定features，RunPCA默认使用FindVariableFeatures找到的高变基因作为PCA输入.
  seurat_phase <- RunPCA(seurat_phase)
  
  # Plot the PCA colored by cell cycle phase
  p1 <- DimPlot(seurat_phase,
                reduction = "pca",
                group.by= "Phase")
  p2 <- DimPlot(seurat_phase,
                reduction = "pca",
                group.by= "Phase",
                split.by = "Phase")
  library(cowplot)
  plot_grid(p1, p2, ncol = 2, labels = "AUTO")
  
  #####Step7 harmony 降维
  library(harmony)
  seuratObj <- RunHarmony(seurat_phase, "orig.ident")
  names(seuratObj@reductions)
  # Set seed
  set.seed(123456)
  # Run UMAP
  # 分群
  #umap降维聚类
  scRNA_harmony <- seuratObj
  
  #确定pc数量
  ElbowPlot(object = scRNA_harmony,ndims = 40)
  scRNA_harmony <- RunUMAP(scRNA_harmony, reduction = "harmony", dims = 1:30)
  scRNA_harmony <- FindNeighbors(scRNA_harmony, reduction = "harmony", dims = 1:20) %>% FindClusters(resolution = 0.5)
  ##作图
  #group_by_cluster
  ##作图##可视化
  ncluster <- length(unique(scRNA_harmony$seurat_clusters))
  library(RColorBrewer)
  mycol <- colorRampPalette(brewer.pal(12, "Set3"))(ncluster)
  
  plot1 = DimPlot(scRNA_harmony, reduction = "umap", label=T,pt.size = 0.5,cols = mycol,group.by = "seurat_clusters") 
  #group_by_sample
  plot2 = DimPlot(scRNA_harmony, reduction = "umap", group.by='orig.ident',pt.size = 0.5) 
  #combinate
  plotc <- plot1+plot2
  
  pdf("/data3/home/hu/Tongfeng/Result/Fig1/umap_cluster.pdf",width = 10,height = 6)
  plotc
  dev.off()
  
  pdf("/data3/home/hu/Tongfeng/Result/Fig1/umap_cluster_0.pdf",width = 7,height = 6)
  plot1
  dev.off()
  
  save.image("/data3/home/hu/Tongfeng/Result/Fig1_0331.RData")
  
  ######step 7
  #细胞注释我粗略的整理出来三种方法，可以供大家选择
  #寻找每一个cluster里与其它所有细胞相比之后的差异marker，并且每个亚群只展示FC最大的两位
  pbmc.markers <- FindAllMarkers(scRNA_harmony, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  
  library(tidyverse)
  top6 <- pbmc.markers %>%
    group_by(cluster) %>%
    top_n(n = 5, wt = avg_log2FC)
  
  write.csv(pbmc.markers,"/data3/home/hu/Tongfeng/Result/Fig1/pbmc.markers.csv")
  write.csv(top6,"top6_marker.csv")
  
  #####使用已知的细胞marker用来注释
  {
    #####single R注释
    library(SingleR)
    library(celldex)
    
    # human
    ref <- HumanPrimaryCellAtlasData() 
    
    sce1 <- scRNA_harmony ###赋值给其他变量，避免修改原变量。
    sce_for_SingleR <- GetAssayData(sce1, slot="data")
    sce_for_SingleR
    clusters <- sce1@meta.data$seurat_clusters
    pred.hesc <- SingleR(test = sce_for_SingleR, 
                         ref = ref, 
                         labels = ref$label.main,
                         method = "cluster", 
                         clusters = clusters, 
                         assay.type.test = "logcounts", 
                         assay.type.ref = "logcounts")
    
    table(pred.hesc$labels)  
    
    
    celltype = data.frame(ClusterID=rownames(pred.hesc), 
                          celltype=pred.hesc$labels, 
                          stringsAsFactors = F) 
    sce1@meta.data$singleR = celltype[match(clusters,celltype$ClusterID),'celltype']
    P9 <- DimPlot(sce1, reduction = "umap", group.by = "singleR",label = T)
    P9
    P10 <- DimPlot(sce1, reduction = "umap", group.by = "seurat_clusters",label = T)
    P10
    #combinate
    plotc <- P9+P10
  }
  
  
  # 创建cluster对应细胞类型
  celltype <- c(
    "T_cells","T_cells","T_cells","NK_cell","Monocyte",
    "Neutrophils","B_cell","B_cell","NK_cell","Monocyte",
    "T_cells","Dendritic_cells","Dendritic_cells","T_cells",
    "T_cells","B_cell","Monocyte","Plasma_cells","T_cells",
    "Plasma_cells","Monocyte","Monocyte","T_cells"
  )
  
  # 确保顺序对应cluster编号
  Idents(scRNA_harmony) <- scRNA_harmony$seurat_clusters
  names(celltype) <- levels(scRNA_harmony)
  
  # 赋值
  scRNA_harmony <- RenameIdents(scRNA_harmony, celltype)
  
  # 保存到meta.data（推荐）
  scRNA_harmony$celltype <- Idents(scRNA_harmony)
  
  scRNA_harmony <- readRDS("/data3/home/hu/Aging/pbmc/scRNA-fina.rds")
  Idents(scRNA_harmony) <- scRNA_harmony$cell_type2
  ####绘制marker基因 dotplot图 ####
  marker_list=list()
  
  marker_list[['Monocyte']]=c('CD14',"S100A8")
  marker_list[['NK_cell']]=c('NKG7', 'GNLY')
  marker_list[['Megakaryocytes']]= c('PPBP', 'ITGA2B')
  marker_list[['T_cell']]=c('CD3D', 'CD3E')
  marker_list[['B_cell']]= c('CD79B', 'MS4A1')
  marker_list[['cDC']]= c( 'FCER1A')
  marker_list[['Plasma_cell']]= c('MZB1', 'JCHAIN')
  marker_list[['pDC']]= c('LILRA4', 'CLEC4C')

  
  p <- DotPlot(object = scRNA_harmony, features=marker_list)
  p
  
  library(RColorBrewer)
  library(tidyverse)
  library(magrittr)
  library(ggh4x)
  pal1 = c("#409079","#52a5c1","#c65341","#425785","#d6873b","#92b8da","#9f2b39",
           "#b5aa82","#de9d3d","#347852","#ca8399","#296097","#564b84","#817f7e")
  #pal1 <- mycol[1:14] 
  
  
  v <- unlist(marker_list, use.names=FALSE)
  p_data = p$data %>% 
    # mutate(id=factor(id,levels=rev(levels(id))))%>% 
    mutate(features.plot=factor(features.plot,levels=rev(v)))
  
  
  dput(levels(p_data$feature.groups))
  
  labels <- scRNA_harmony$cell_type2
  col_use = data.frame(group =
                         unique(labels)) %>% 
    mutate(mycol = pal1[1:nrow(.)]) %>% 
    arrange(group)
  
  col_map_y = p_data %>% 
    distinct(features.plot,feature.groups) %>% 
    set_colnames(c('name','group')) %>% 
    left_join(col_use ,by='group')
  col_map_y
  
  p11 <- ggplot(p_data, aes(x = id, y = features.plot)) + 
    geom_point(aes(size = pct.exp, color = avg.exp.scaled)) + 
    scale_radius(breaks = c(25, 50, 75, 100), range = c(0,6)) + 
    theme_test() + 
    scale_color_gradient2(low = "steelblue", mid = "white", high = "firebrick")+
    theme(axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1,face = 1, size = 12, 
                                     color = "black"), 
          
          axis.text.y = element_text(size = 10, face = 1, color = "black"),
          legend.text = element_text(size = 10, face = 1), 
          legend.title = element_text(size =10, face = 1),
          legend.position = 'top')+
    labs(x=NULL,y=NULL,
         size = "Fraction of cells", color = "Mean Expression")+
    guides(y = "axis_nested")+
    theme(
      ggh4x.axis.nestline.y = element_line(size = 2,color = rev(col_use$mycol)),
      ggh4x.axis.nesttext.y = element_text(colour =rev( col_use$mycol)))+
    scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = 100))
  p11
  ggsave(filename = "/data3/home/hu/Tongfeng/new_Result/Fig1/marker_dotplot.pdf",plot = p11,height = 5.5,width = 8)
  
  
  #####Umap图美化
  P10 <- DimPlot(scRNA_harmony, reduction = "umap",label = T,group.by = "cell_type2")
  P10
  
  library(ggplot2)
  library(dplyr)
  library(tidydr)
  library(mascarade)
  #---- 01 获取数据，设置配色 ---------------------------------------
  plot_df <- cbind(scRNA_harmony@reductions$umap@cell.embeddings, scRNA_harmony@meta.data)
  
  marker_list <- list()
  marker_list[['Monocyte']]=c('CD14',"S100A8")
  marker_list[['NK_cell']]=c('NKG7', 'GNLY')
  marker_list[['Megakaryocytes']]= c('PPBP', 'ITGA2B')
  marker_list[['T_cell']]=c('CD3D', 'CD3E')
  marker_list[['B_cell']]= c('CD79B', 'MS4A1')
  marker_list[['cDC']]= c( 'FCER1A')
  marker_list[['Plasma_cell']]= c('MZB1', 'JCHAIN')
  marker_list[['pDC']]= c('LILRA4', 'CLEC4C')
  
  
  # 基因表达均值热图
  scRNA_harmony$celltype <- factor(scRNA_harmony$cell_type,levels = names(marker_list))
  DefaultAssay(scRNA_harmony) <- "RNA"  # 可根据实际情况切换为 "SCT"
  scRNA_v4 <- Convert_Assay(seurat_object = scRNA_harmony, convert_to = "V3")
  
  
  ht1 <- GroupHeatmap(scRNA_v4,
                      features = unlist(marker_list),
                      group.by = c("celltype"),
                      show_row_names = TRUE,
                      row_names_side = 'left',
                      #legend_title = 'z-score',
                      column_title = ''
                      
  )
  ht1$plot
  
  
  # 面板大小固定
  pdf("/data3/home/hu/Tongfeng/new_Result/Fig1/card_marker.pdf",height = 5,width = 7)
  panel_fix(ht1$plot, height = 4, width = 6)
  dev.off()
  
  
  color_mapping <-  c("T_cell" = "#B44D51","NK_cell" =  "#DA884F","Megakaryocytes" = "#8AB6C3","B_cell" = "#D0CC48",
                      "cDC" = "#989EBE","Plasma_cell" = "#B8A2C4","Monocyte" = "#386388","pDC" = "#564b84")
  #---- 02 计算质心距离 ---------------------------------------------
  maskTable <- generateMask(
    dims = plot_df[, c("umap_1", "umap_2")], 
    cluster = plot_df$cell_type2, 
    minDensity = 1.0, 
    smoothSigma = 0.05
  )
  #---- 03 作图 -----------------------------------------------------
  p1 <- ggplot(plot_df, aes(x = umap_1, y = umap_2)) +
    geom_point(aes(color = cell_type2), size = 0.05) +  # 使用 Celltype 映射颜色
    geom_path(data = maskTable, aes(group = group), 
              linewidth = 0.3, linetype = 2, color = "black") +
    coord_fixed(ratio = 1) +
    theme_classic() +
    scale_color_manual(values = color_mapping,guide = guide_legend(override.aes = list(size = 4)) ) +  # 手动指定颜色映射
    theme_dr(xlength = 0.22, ylength = 0.22, 
             arrow = grid::arrow(length = unit(0.15, "inches"), type = "closed")) +
    theme(
      panel.grid = element_blank(),
      legend.position = "right",  # 图例在右侧
      legend.title = element_text(face = "bold")  # 图例标题加粗
    )
  
  ggsave(filename = "/data3/home/hu/Tongfeng/new_Result/Fig1/scRNA_celltype.pdf",plot = p1,
         width = 7,height = 5)
  
  
  
  # 导入数据
  library(SCP)
  scRNA_harmony@meta.data$celltype <- scRNA_harmony@meta.data$cell_type2
  sce.all.int <- scRNA_harmony
  
  
  # 统一颜色映射
  cluster_levels <- unique(sce.all.int$celltype)
  cluster_levels <- cluster_levels[order(cluster_levels)]
  cluster_colors <- setNames(c("#946A85","#989EBE","#8AB6C3","#386388","#DA884F","#564b84","#B8A2C4","#B44D51"), cluster_levels)    
  
  color_mapping <-  c("T_cell" = "#B44D51","NK_cell" =  "#DA884F","Megakaryocytes" = "#8AB6C3","B_cell" = "#946A85",
                      "cDC" = "#989EBE","Plasma_cell" = "#B8A2C4","Monocyte" = "#386388","pDC" = "#564b84")
  

  
  
  sce.all.int$Ecotyper <- factor(sce.all.int$celltype,levels =cluster_levels )
  # 经典框框
  p1 <- CellDimPlot(sce.all.int, group.by = "celltype", pt.alpha = 1,
                    reduction = "UMAP", show_stat = F,label_insitu = T,pt.size = 0.005,
                    palcolor = cluster_colors,raster = TRUE,label = T );p1
  
  sce.all.int$group <- factor(sce.all.int$group,levels = c("Control","Disease","Remission"))
  p2 <- CellDimPlot(sce.all.int, group.by = "celltype", split.by = "group",
                    reduction = "UMAP", show_stat = F,label_insitu = T,pt.alpha = 1,
                    palcolor = cluster_colors,raster = TRUE,label = F,pt.size = 0.5,bg_color = "#FFFFFF");p2

  
  
  
  ggsave(filename = "/data3/home/hu/Tongfeng/new_Result//Fig1/scRNA_umap.pdf",plot = p1,
         width = 7,height = 5)
  
  ggsave(filename = "/data3/home/hu/Tongfeng/new_Result//Fig1/scRNA_umap_split.pdf",plot = p2,
         width = 20,height = 5)
  
  saveRDS(scRNA_harmony,"/data3/home/hu/Tongfeng/new_Result//Fig1/Fig1_scRNA_harmony.RDS")
  
  scRNA_harmony <- readRDS("/data3/home/hu/Tongfeng/Result/Fig1/Fig1_scRNA_harmony.RDS")
  
  
  
  #####一些feature plot的绘制
  library(scCustomize)

  gene <- c("CD3D","NKG7","S100A8","MS4A1")
  
  plots <- list()
  for (i in 1:length(gene)){
    
    plots[[i]] <- FeaturePlot_scCustom(seurat_object = scRNA_harmony, features = gene[i],raster = T,
                         colors_use = colorRampPalette(c("#EFEFEF", "#D1E5F0", "#67A9CF", "#2166AC", "#08306B"))(50))+NoAxes()+NoLegend()+NoAxes()

  }
  
  
  library(patchwork)
  p<-wrap_plots(plots, ncol = 2);p
  
  library(ComplexHeatmap)
  library(circlize)
  #给图例定义范围和颜色
  col_fun = colorRamp2(c(0,1,2,3,4), c("#EFEFEF", "#D1E5F0", "#67A9CF", "#2166AC", "#08306B"))
  lgd3 = Legend(at = c(0,2,4), col_fun = col_fun, 
                nrow = 1,#border = "black",direction = "horizontal",
                legend_height = unit(4, "cm")
  )
  
  draw(lgd3,x = unit(0.99, "npc"), y = unit(0.4, "npc"),
       just = c("right", "bottom"))
  
  #可以调整图例的方向，高度，还有加不加border。
  #用draw把图例画在合适的位置。
  p1 <- p+draw(lgd3,x = unit(0.99, "npc"), y = unit(0.4, "npc"),
               just = c("right", "bottom"))
  
  library(colorspace)
  hcl_palettes(plot = TRUE)
  
  ggsave(p,file="/data3/home/hu/Tongfeng/new_Result/Fig1/featureplot_subset4.pdf",width =6 ,height = 6)
  


  scRNA_harmony$group <- factor(scRNA_harmony$group,levels = c("Control","Disease",
                                                               "Remission"))
  
  # 细胞比例绘制
  pdf("/data3/home/hu/Tongfeng/new_Result/Fig1/cell_radio.pdf",width = 5,height = 5)
  CellStatPlot(scRNA_harmony,
               stat.by = "celltype",palcolor = cluster_colors,
               group.by = "group",plot_type = "trend",
               label = TRUE) %>% panel_fix(height = 16, width = 3)
  dev.off()
} 


#####绘制炎症和代谢在不同细胞类型总体趋势变化
{
  ###炎症基因集和代谢基因集
  # 加载代谢相关基因
  library(readxl)
  
  # 读取文件
  data <- read_excel("/data3/home/hu/Tongfeng/Result/Fig2/mmc2.xlsx",skip = 1)
  kegg_data <- data[grepl("KEGG", data$From), ]  
  KEGG_metabolism_gene <- kegg_data$`Gene symbol`
  
  library(msigdbr)
  library(dplyr)
  
  hallmark <- msigdbr(species = "Homo sapiens", category = "H")
  inflam_genes <- hallmark %>%
    filter(gs_name == "HALLMARK_INFLAMMATORY_RESPONSE") %>%
    pull(gene_symbol) %>%
    unique()
  
  fatty_acid_genes <- hallmark %>%
    filter(gs_name == "HALLMARK_FATTY_ACID_METABOLISM") %>%
    pull(gene_symbol) %>%
    unique()
  
  glycolysis_genes <- hallmark %>%
    filter(gs_name == "HALLMARK_GLYCOLYSIS") %>%
    pull(gene_symbol) %>%
    unique()
  
  gene_list <- list(inflam_genes,fatty_acid_genes,glycolysis_genes)
  names(gene_list) <- c("Inflammatory","Fatty acid metabolism","Glycolysis")
  
  ###对单细胞进行伪bulk处理
  library(qs)
  library(ggplot2)
  library(DESeq2)
  library(Seurat)
  library(tidyverse)
  library(dplyr)


  # 使用 AggregateExpression 进行加和处理
  # assays: 选择要聚合的 assay
  # return.seurat: TRUE 会返回一个简化版的 Seurat 对象，方便后续操作
  pb_seurat <- AggregateExpression(scRNA_harmony, 
                                   group.by = "orig.ident", 
                                   assays = "RNA", 
                                   slot = "counts", 
                                   return.seurat = TRUE)
  
  # 此时 pb_seurat 的每一个 "细胞" 实际上就是一个组的加和结果
  # 基因表达矩阵如下：
  pb_matrix <- pb_seurat@assays$RNA$counts
  # 提取聚合后的 count 矩阵
  data_pb <- as.matrix(pb_seurat@assays$RNA$counts)
  
  # 规范化行名（基因名）和列名（组名）
  data_pb_rounded <- round(data_pb)
  
  # 根据列名（组名）创建对应的条件向量
  # 假设你的组别是有顺序的，比如从健康到严重
  group_names <- colnames(data_pb_rounded)
  group_data <- scRNA_harmony@meta.data[,c("orig.ident","group")]
  group_data <- unique(group_data)
  rownames(group_data) <- group_data$orig.ident
  condition <- group_data[group_names,]$group
  condition <- factor(condition,levels = c("Control","Remission","Disease"))
  
  condition <- as.numeric(condition)
  # 创建 colData
  col_data <- data.frame(
    row.names = group_names,
    condition = condition
  )
  
  
  library(DESeq2)
  
  # 构建 DESeqDataSet
  dds <- DESeqDataSetFromMatrix(countData = data_pb_rounded, 
                                colData = col_data, 
                                design = ~ condition)
  
  # 运行归一化
  dds <- estimateSizeFactors(dds)
  resdata <- as.data.frame(counts(dds, normalized = TRUE))
  
  # Z-score 标准化
  scale_data <- as.data.frame(t(scale(t(resdata))))
  scale_data <- as.data.frame(t(scale_data))
  
  gene_name <- names(gene_list)
  
  for (i in 1:length(gene_list)) {
    
    immue_gene <- intersect(gene_list[[i]], colnames(scale_data))
    
    immue_data <- scale_data[, immue_gene]
    immue_data$time <- condition
    median_data <- aggregate(immue_data[,-ncol(immue_data)], list(immue_data$time), median)
    
    genes <- colnames(median_data)[-1]
    
    df_long <- tidyr::pivot_longer(median_data, cols = all_of(genes),
                                   names_to = "Gene", values_to = "Expression")
    df_long <- na.omit(df_long)
    
    # MK检验
    library(trend)
    mk_test <- mk.test(df_long$Expression, continuity = TRUE)
    P_value <- mk_test$p.value
    Z_score <- as.numeric(mk_test$statistic)
    
    # 绘图
    col <- c("#EEBF6D", "#d0b9d5", "#B8B7A3", "#B8B7A3", "#9E3150",  
             "#60b4e2", "#FC8C5A", "#D44C3C", "#B1EAEF", "#8EB69C", "#EEBF6D",
             "#44757A", "#c8e27f", "#fca1c2")
    
    library(ggplot2)
    library(ggthemes)
    
    p <- ggplot(df_long) +
      theme_par() +
      guides(colour = "none") +
      geom_smooth(aes(x = Group.1, y = Expression, color = Gene),
                  method = "loess", se = FALSE, size = 0.2) +
      scale_x_continuous(
        breaks = c(1:3),
        labels = c(1:3)
      ) +
      ylim(-1.2, 1.2) +
      labs(x = "Group", y = "Zscore") +
      ggtitle(gene_name[i]) +
      theme(plot.title = element_text(size = 12, hjust = 0.5)) +
      scale_color_manual(values = rep(col[i], length(genes)), labels = "None") +
      geom_smooth(aes(x = Group.1, y = Expression),
                  method = "loess", se = TRUE, size = 1.3,
                  fill = "white", color = col[i], alpha = 1) +
      annotate("text", x = 1, y = 1.05,
               label = paste0("Z = ", round(Z_score, 2),
                              "\tP = ", formatC(P_value, format = "e", digits = 2)),
               size = 4, color = "black", hjust = 0);p
    
    # 保存图像
    file_name <- paste0("/data3/home/hu/Tongfeng/new_Result//Fig1/Trend/",gene_name[i],"total.pdf")
    ggsave(file_name, plot = p, width = 7, height = 4)
    
  }
  
  
  #####按照细胞类型再进行划分查看
  celltype <- unique(scRNA_harmony$celltype)
  
  for (j in 1:length(celltype)) {
    temp_scRNA <- subset(scRNA_harmony,subset=celltype == celltype[j])
    
    pb_seurat <- AggregateExpression(temp_scRNA, 
                                     group.by = "orig.ident", 
                                     assays = "RNA", 
                                     slot = "counts", 
                                     return.seurat = TRUE)
    
    # 此时 pb_seurat 的每一个 "细胞" 实际上就是一个组的加和结果
    # 基因表达矩阵如下：
    pb_matrix <- pb_seurat@assays$RNA$counts
    # 提取聚合后的 count 矩阵
    data_pb <- as.matrix(pb_seurat@assays$RNA$counts)
    
    # 规范化行名（基因名）和列名（组名）
    data_pb_rounded <- round(data_pb)
    
    # 根据列名（组名）创建对应的条件向量
    # 假设你的组别是有顺序的，比如从健康到严重
    group_names <- colnames(data_pb_rounded)
    group_data <- scRNA_harmony@meta.data[,c("orig.ident","group")]
    group_data <- unique(group_data)
    rownames(group_data) <- group_data$orig.ident
    condition <- group_data[group_names,]$group
    condition <- factor(condition,levels = c("Control","Int gout remission",
                                             "Adv gout remission","Acute flare"))
    
    condition <- as.numeric(condition)
    # 创建 colData
    col_data <- data.frame(
      row.names = group_names,
      condition = condition
    )
    
    
    library(DESeq2)
    
    # 构建 DESeqDataSet
    dds <- DESeqDataSetFromMatrix(countData = data_pb_rounded, 
                                  colData = col_data, 
                                  design = ~ condition)
    
    # 运行归一化
    dds <- estimateSizeFactors(dds)
    resdata <- as.data.frame(counts(dds, normalized = TRUE))
    
    # Z-score 标准化
    scale_data <- as.data.frame(t(scale(t(resdata))))
    scale_data <- as.data.frame(t(scale_data))
    
    gene_name <- names(gene_list)
    
    for (i in 1:length(gene_list)) {
      
      immue_gene <- intersect(gene_list[[i]], colnames(scale_data))
      
      immue_data <- scale_data[, immue_gene]
      immue_data$time <- condition
      median_data <- aggregate(immue_data[,-ncol(immue_data)], list(immue_data$time), median)
      
      genes <- colnames(median_data)[-1]
      
      df_long <- tidyr::pivot_longer(median_data, cols = all_of(genes),
                                     names_to = "Gene", values_to = "Expression")
      df_long <- na.omit(df_long)
      
      # MK检验
      library(trend)
      mk_test <- mk.test(df_long$Expression, continuity = TRUE)
      P_value <- mk_test$p.value
      Z_score <- as.numeric(mk_test$statistic)
      
      # 绘图
      col <- c("#EEBF6D", "#d0b9d5", "#B8B7A3", "#B8B7A3", "#9E3150",  
               "#60b4e2", "#FC8C5A", "#D44C3C", "#B1EAEF", "#8EB69C", "#EEBF6D",
               "#44757A", "#c8e27f", "#fca1c2")
      
      library(ggplot2)
      library(ggthemes)
      
      p <- ggplot(df_long) +
        theme_par() +
        guides(colour = "none") +
        geom_smooth(aes(x = Group.1, y = Expression, color = Gene),
                    method = "loess", se = FALSE, size = 0.2) +
        scale_x_continuous(
          breaks = c(0:8),
          labels = c(0:8)
        ) +
        ylim(-1.2, 1.2) +
        labs(x = "Group", y = "Zscore") +
        ggtitle(gene_name[i]) +
        theme(plot.title = element_text(size = 12, hjust = 0.5)) +
        scale_color_manual(values = rep(col[i], length(genes)), labels = "None") +
        geom_smooth(aes(x = Group.1, y = Expression),
                    method = "loess", se = TRUE, size = 1.3,
                    fill = "white", color = col[i], alpha = 1) +
        annotate("text", x = 3, y = 1.05,
                 label = paste0("Z = ", round(Z_score, 2),
                                "\tP = ", formatC(P_value, format = "e", digits = 2)),
                 size = 4, color = "black", hjust = 0);p
      
      # 保存图像
      file_name <- paste0("/data3/home/hu/Tongfeng/new_Result//Fig1/Trend/",celltype[j],"_",gene_name[i],".pdf")
      ggsave(file_name, plot = p, width = 7, height = 4)
      
    }
    
  }
  
}


####差异分析
scRNA_harmony$celltype <-  factor(scRNA_harmony$celltype,levels =names(marker_list) )
Idents(scRNA_harmony) <- scRNA_harmony$celltype
deg <- FindAllMarkers(scRNA_harmony, min.pct = 0.25, logfc.threshold =0.25)
View(deg)
library(scRNAtoolVis)

#用title.col参数修改矩形的填充颜色
p1 <- jjVolcano(diffData = deg,         
          tile.col = cluster_colors)

library(ggplot2)

# ==============================================================================
# 1. 调用 jjVolcano 进行核心参数瘦身与精调
# ==============================================================================
p1 <- jjVolcano(
  diffData = deg,         
  tile.col = cluster_colors,
  
  # --- 【核心瘦身：减少大小】 ---
  pt.size = 0.5,             # 1. 缩小火山图散点大小（默认一般是1或1.5，改成0.5最精致）
  gene.size = 2.5,           # 2. 缩小基因文本标签的字号（从默认的3-4降到2.5，防止撑大图幅）
  tile.height = 0.08,        # 3. 压扁顶部分组彩色条的高度，给下方火山预留更多纵向空间
  
  # --- 【标签密度与排版控制】 ---
  topGeneN = 5,              # 4. 每组只标记最显著的前5个基因（默认可能过多，导致文字爆满）
  col.alpha = 0.75,          # 5. 赋予散点轻微透明度，叠在一起时更有晶格层次感
  
  # --- 【黑科技：防止文字重叠】 ---
  # 如果 jjVolcano 提示支持 repel 或 force 参数，可显式指定（若不支持可删掉这两行）：
  # max.overlaps = 20,       # 限制重叠标签的最大显色数
  # force = 1.5,             # 加大标签之间的排斥力
)

# ==============================================================================
# 2. 叠加 ggplot2 顶级科研硬朗主题修饰
# ==============================================================================
p1_fixed <- p1 + 
  theme_bw(base_family = "sans") + # 使用中性无衬线字体，极具现代感
  theme(
    # 核心绘图区外边框：加上极具质感的细黑框，拒绝发飘
    panel.border = element_rect(color = "#333333", fill = NA, size = 0.8),
    panel.grid.major = element_blank(), # 抠掉背景凌乱的网格线
    panel.grid.minor = element_blank(),
    
    # 坐标轴字体调小加粗，迎合整体瘦身风格
    axis.text.x = element_text(size = 9, color = "black", face = "bold"),
    axis.text.y = element_text(size = 9, color = "black"),
    axis.title = element_text(size = 10, face = "bold", color = "#222222"),
    
    # 图例精细化缩放
    legend.position = "right",
    legend.title = element_text(size = 8.5, face = "bold"),
    legend.text = element_text(size = 8)
  )

# ==============================================================================
# 3. 【最关键一步】通过锁定 ggsave 物理尺寸彻底“减少大小”
# ==============================================================================
# 直接在屏幕上查看
print(p1_fixed)

# 保存为 PDF 时，通过精准控制 width 和 height，将其压缩成紧凑的精致小图
ggsave(
  filename = "Premium_Mini_jjVolcano.pdf", 
  plot = p1_fixed, 
  width = 7.5,      # 【核心控制】横向宽度收窄，防止由于细胞群太多导致火山图拉得巨长
  height = 4,     # 黄金纵向高度
  dpi = 600
)

ggsave("/data3/home/hu/Tongfeng/new_Result/Fig1/volcano.pdf", plot = p1, width = 7, height = 4)

qsave(scRNA_harmony,"/data3/home/hu/Tongfeng/new_Result//Fig1/scRNA_harmony.qs")



#####进行一个Hallmarker富集分析
{
  library(msigdbr)
  
  # 1. 从 msigdbr 提取人类 (Homo sapiens) 的 Hallmark 基因集 (Category 为 H)
  m_df <- msigdbr(species = "Homo sapiens", category = "H")
  
  # 2. 将数据框格式转换为 AUCell 需要的 list 格式
  hallmark_list <- split(x = m_df$gene_symbol, f = m_df$gs_name)
  
  # 3. 检查一下
  length(hallmark_list) # 应该返回 50
  names(hallmark_list)[1:3] # 查看前三个通路的名称
  
  # 过滤掉不在你 scRNA 对象中的基因
  hallmark_list <- lapply(hallmark_list, function(x) intersect(x, rownames(scRNA)))
  
  # 过滤掉基因数太少（比如少于 5 个基因）的通路
  hallmark_list <- hallmark_list[lapply(hallmark_list, length) > 5]
  
  library(AUCell)
  # 构建排名（非常耗时，建议开启多核心）
  cells_rankings <- AUCell_buildRankings(scRNA[["RNA"]]$counts, nCores = 8) 
  
  # 计算 AUC
  cells_AUC <- AUCell_calcAUC(hallmark_list, cells_rankings)
  
  # 提取结果到 Seurat 对象的 meta.data
  auc_results <- as.data.frame(t(getAUC(cells_AUC)))
  scRNA <- AddMetaData(scRNA, metadata = auc_results)
  
  library(ComplexHeatmap)
  library(dplyr)
  library(circlize)
  
  # 1. 提取 meta.data 中的 Hallmark 得分 (假设列名以 HALLMARK_ 开头)
  pathway_names <- grep("^HALLMARK_", colnames(scRNA@meta.data), value = TRUE)
  df_score <- scRNA@meta.data[, c("group","cell_type2","orig.ident", pathway_names)]
  saveRDS(df_score,"/data3/home/hu/Tongfeng/new_Result/Fig2/df_score_hallmarker.RDS")
  
  df_score <- readRDS("/data3/home/hu/Tongfeng/new_Result/Fig2/df_score_hallmarker.RDS")
  # 2. 按组别计算平均分
  avg_score <- df_score[,c(-2,-3)] %>%
    dplyr::group_by(group) %>%
    dplyr::summarise(across(everything(), mean)) %>%
    as.data.frame()
  
  # 3. 整理矩阵：行名为通路，列名为组别
  rownames(avg_score) <- avg_score$group
  plot_mtx <- t(avg_score[, -1]) # 转置，让通行为行
  colnames(plot_mtx) <- avg_score$group
  
  # 4. 重新排序组别（确保按照 Control -> Disease -> Remission 显示）
  plot_mtx <- plot_mtx[, c("Control", "Disease", "Remission")]
  
  # 计算每行（通路）的标准差，选出变化最明显的 Top 20
  top_gap <- names(sort(apply(plot_mtx, 1, sd), decreasing = TRUE))[1:20]
  
  {
    # 定义筛选列表
    metabolic_related <- c(
      "HALLMARK_GLYCOLYSIS", "HALLMARK_OXIDATIVE_PHOSPHORYLATION", 
      "HALLMARK_FATTY_ACID_METABOLISM", "HALLMARK_CHOLESTEROL_HOMEOSTASIS",
      "HALLMARK_XENOBIOTIC_METABOLISM", "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY",
      "HALLMARK_MTORC1_SIGNALING", "HALLMARK_ADIPOGENESIS"
    )
    
    immune_inflammatory_related <- c(
      "HALLMARK_INFLAMMATORY_RESPONSE", "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
      "HALLMARK_IL6_JAK_STAT3_SIGNALING", "HALLMARK_COMPLEMENT",
      "HALLMARK_INTERFERON_GAMMA_RESPONSE", "HALLMARK_INTERFERON_ALPHA_RESPONSE",
      "HALLMARK_IL2_STAT5_SIGNALING", "HALLMARK_ALLOGRAFT_REJECTION",
      "HALLMARK_HYPOXIA", "HALLMARK_COAGULATION"
    )
    
    # 合并关注的通路
    target_pathways <- c(metabolic_related, immune_inflammatory_related)
    
    # 确保这些通路在你的数据中存在
    target_pathways <- intersect(target_pathways, rownames(plot_mtx))
    
  }
  plot_mtx_sub <- plot_mtx[target_pathways, ]
  
  # 数据标准化（Z-score），让颜色对比更明显
  plot_mtx_z <- t(scale(t(plot_mtx_sub)))
  
  # 定义颜色：蓝色（低评分）到白色到红色（高评分）
  col_fun = colorRamp2(c(-2, 0, 2), c("#313695", "white", "#a50026"))
  
  library(ComplexHeatmap)
  library(circlize)
  
  # 1. 优化颜色：使用更专业的双极色板 (RdYlBu 或 RdBu)
  col_fun = colorRamp2(c(-2, 0, 2), c("#4575b4", "white", "#d73027"))
  
  # 2. 准备侧栏分类（代谢 vs 炎症）
  # 假设 target_pathways 是你筛选出的通路名
  row_anno_df <- data.frame(
    Type = ifelse(grepl("METABOLISM|GLYCOLYSIS|PHOSPHORYLATION|LIPID|CHOLESTEROL", rownames(plot_mtx_z)), 
                  "Metabolism", "Immune/Inflam")
  )
  rownames(row_anno_df) <- rownames(plot_mtx_z)
  
  # 3. 定义侧栏颜色
  anno_colors = list(Type = c("Metabolism" = "#9C6678", "Immune/Inflam" = "#BD6C5C"))
  
  # 开始绘图
  library(ggplot2)
  library(reshape2)
  library(dplyr)
  
  # ==============================================================================
  # 1. 数据重构（针对纯 AUC 打分值进行清洗与排序）
  # ==============================================================================
  # 将您的 AUCell 矩阵转换为长表格格式
  melt_df <- melt(plot_mtx_z)
  colnames(melt_df) <- c("Pathway", "Group", "AUC_score")
  
  # 名字净化：去掉大写 HALLMARK_ 缀，并将下划线换为空格
  melt_df$Pathway <- gsub("HALLMARK_", "", melt_df$Pathway)
  melt_df$Pathway <- gsub("_", " ", melt_df$Pathway)
  
  # 补全分类标签
  melt_df$Type <- rep(row_anno_df$Type, times = length(unique(melt_df$Group)))
  melt_df$Type <- factor(melt_df$Type, levels = c("Immune/Inflam", "Metabolism"))
  
  # 锁定 X 轴分组顺序
  melt_df$Group <- factor(melt_df$Group, levels = c("Control", "Disease", "Remission"))
  
  # 【智能趋势排序】按平均 AUC 评分从低到高排序，使气泡在视觉上呈现完美的对角线倾斜，告别凌乱
  pathway_order <- melt_df %>%
    dplyr::group_by(Pathway) %>%
    dplyr::summarize(mean_auc = mean(AUC_score)) %>%
    dplyr::arrange(mean_auc) %>%
    pull(Pathway)
  melt_df$Pathway <- factor(melt_df$Pathway, levels = pathway_order)
  
  # ==============================================================================
  # 2. 轴引导线与核心图表配置
  # ==============================================================================
  p_auc_bubble <- ggplot(melt_df, aes(x = Group, y = Pathway)) +
    
    # 背景水平引导线，辅助视线水平对齐
    geom_hline(aes(yintercept = Pathway), color = "#EAEAEA", linetype = "solid", size = 0.4) +
    
    # 绘制核心气泡：颜色和大小双维度均绑定 AUC 评分
    # 颜色代表高低变化，大小代表该通路在该组中整体的平均激活丰度
    geom_point(aes(size = AUC_score, color = AUC_score), alpha = 0.95) +
    
    # 按分类进行物理分块
    facet_grid(Type ~ ., scales = "free_y", space = "free_y") + 
    
    # --- 顶级学术调色盘 ---
    # 颜色使用高对比度的反转红蓝（高表达深红，低表达深蓝）
    scale_color_distiller(palette = "RdBu", direction = -1, oob = scales::squish) +
    # 科学控制气泡圆点的大小范围
    scale_size_continuous(range = c(2.5, 7.5)) +
    
    # --- 硬朗科研边框主题 ---
    theme_bw(base_family = "sans") + 
    labs(
      title = "Metabolic & Immune Landscape (AUCell Activity)", 
      x = NULL, 
      y = NULL, 
      size = "Mean AUC Score",       # 【修正】图例名称改为平均AUC得分
      color = "Pathway Activity\n(AUC score)"
    ) +
    
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 12, color = "#222222", margin = margin(b = 12)),
      
      # 核心绘图区外边框：加上极具质感的细黑框
      panel.border = element_rect(color = "#333333", fill = NA, size = 0.8),
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(),
      panel.spacing = unit(5, "mm"), # 两个大分块之间的物理间隔，避免黏连
      
      # 坐标轴美化
      axis.text.x = element_text(angle = 30, hjust = 1, face = "bold", size = 10.5, color = "black"),
      axis.text.y = element_text(size = 9.5, color = "black", face = "plain"),
      axis.ticks = element_line(color = "#333333", size = 0.5),
      
      # 右侧条带文字设置（初始背景会被下一步的黑科技动态覆盖）
      strip.text.y = element_text(angle = -90, face = "bold", size = 10, color = "white"),
      strip.background = element_rect(fill = "grey20", color = "#333333", size = 0.8), 
      
      # 图例排版
      legend.position = "right",
      legend.title = element_text(size = 9, face = "bold"),
      legend.text = element_text(size = 8.5),
      legend.box.spacing = unit(6, "mm")
    )
  
  # ==============================================================================
  # 3. 【黑科技】强行将您的专属 anno_colors 注入右侧条带
  # ==============================================================================
  g <- ggplotGrob(p_auc_bubble)
  strips <- grep("strip-r|strip-t|strip_layout", names(g$grobs))
  if(length(strips) == 0) strips <- grep("strip", g$layout$name)
  
  # 映射您指定的：Immune/Inflam = #BD6C5C, Metabolism = #9C6678
  strip_colors <- c("#BD6C5C", "#9C6678") 
  j <- 1
  for (i in strips) {
    if (g$grobs[[i]]$layout$name == "strip") {
      g$grobs[[i]]$grobs[]$children[]$gp$fill <- strip_colors[j]
      j <- j + 1
    }
  }
  
  # ==============================================================================
  # 4. 打印并安全导出矢量 PDF (AI秒开，完美保留彩色条带)
  # ==============================================================================
  # 在 R 控制台或 RStudio 中预览
  grid::grid.draw(g)
  
  # 导出为 PDF 文件
  pdf("Premium_AUCell_ColorStrip_Plot.pdf", width = 5.8, height = 6.2)
  grid::grid.draw(g)
  dev.off()
  
  
  
  
  
}

####绘制按照组别的箱线图
{
  # 确保安装并加载核心可视化与统计包
  if (!requireNamespace("ggpubr", quietly = TRUE)) install.packages("ggpubr")
  library(ggplot2)
  library(ggpubr)
  library(patchwork)
  
  # ==============================================================================
  # 1. 定义需要分析的通路和组别顺序
  # ==============================================================================
  target_pathways <- c(
    "HALLMARK_FATTY_ACID_METABOLISM", 
    "HALLMARK_GLYCOLYSIS", 
    "HALLMARK_INFLAMMATORY_RESPONSE", 
    "HALLMARK_TNFA_SIGNALING_VIA_NFKB"
  )
  
  # 确保 group 列是因子，并严格按照实验逻辑排序
  df_score$group <- factor(df_score$group, levels = c("Control", "Disease", "Remission"))
  
  # 定义显著性两两比较的组合列表
  my_comparisons <- list(
    c("Control", "Disease"),
    c("Disease", "Remission"),
    c("Control", "Remission")
  )
  
  # 使用您之前指定的标志性分组配色（可根据您的喜好微调）
  group_colors <- c("Control" = "#77ABC5", "Disease" = "#BD6C5C", "Remission" = "#E2B18E")
  
  # ==============================================================================
  # 2. 自动化循环批量建图
  # ==============================================================================
  box_plots <- list()
  
  for (pathway in target_pathways) {
    
    # 净化图表大标题：去掉 HALLMARK_ 缀，使版面高级
    clean_title <- gsub("HALLMARK_", "", pathway)
    clean_title <- gsub("_", " ", clean_title)
    
    p <- ggboxplot(
      data = df_score, 
      x = "group", 
      y = pathway, 
      fill = "group",
      palette = group_colors,
      
      # --- 视觉润色参数 ---
      alpha = 0.85,           # 填充轻微透明度，极具呼吸感
      width = 0.5,            # 箱体宽度，不要太粗，0.5 最优雅
      outlier.shape = 16,     # 异常点显示为实心圆点
      outlier.size = 0.8,     # 异常点缩小，防止喧宾夺主
      outlier.color = "grey60",
      
      # --- 边框与极简背景 ---
      ggtheme = theme_classic()
    ) +
      # --- 自动化添加显著性标记 ---
      stat_compare_means(
        comparisons = my_comparisons, 
        method = "wilcox.test",        # 使用经典的非参数 Wilcoxon 检验
        label = "p.signif",            # 仅显示星号（*, **, ***, ns），界面更干净
        step.increase = 0.12,          # 自动调整多重比较显著性连线的高度间隔，防止重叠
        tip.length = 0.02
      ) +
      # 全局整体显著性（Kruskal-Wallis 检验，显示在图表顶部）
      stat_compare_means(label.y.npc = "top", label.x.npc = "left", size = 3.5) +
      
      # 细节美化
      labs(title = clean_title, x = NULL, y = "AUCell Score") +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 10, color = "#222222"),
        axis.text.x = element_text(face = "bold", size = 10, color = "black"),
        axis.text.y = element_text(size = 9, color = "black"),
        axis.title.y = element_text(size = 10, face = "bold"),
        legend.position = "none"       # 隐藏多余的下方图例，因为 X 轴已有标签
      )
    
    box_plots[[pathway]] <- p
  }
  
  # ==============================================================================
  # 3. 拼图与高质量矢量导出 (AI秒开)
  # ==============================================================================
  # 将 4 张图拼成 2x2 的精致矩阵图
  final_merged_box <- wrap_plots(box_plots, ncol = 2)
  
  # 在 Rstudio 界面中直接查看
  print(final_merged_box)
  
  # 导出为高清晰度矢量 PDF，文件仅几十 KB，AI 可以完美拆解所有线条、箱体和文字
  ggsave(
    filename = "Gout_PBMC_Pathway_Boxplots.pdf", 
    plot = final_merged_box, 
    width = 7.5,    # 黄金宽比例
    height = 7.0,   # 留足纵向显著性横线的伸展空间
    dpi = 600
  )
  
}
####绘制细胞类型的箱线图
{
  library(ggplot2)
  library(reshape2)
  library(dplyr)
  library(ggpubr)
  
  # ==============================================================================
  # 1. 数据重构：按【样本 + 细胞类型 + 组别】计算平均得分（Pseudo-bulk）
  # ==============================================================================
  target_pathways <- c(
    "HALLMARK_FATTY_ACID_METABOLISM", 
    "HALLMARK_GLYCOLYSIS", 
    "HALLMARK_INFLAMMATORY_RESPONSE", 
    "HALLMARK_TNFA_SIGNALING_VIA_NFKB"
  )
  
  # 计算伪 bulk 均值
  df_cell_bulk <- df_score %>%
    dplyr::group_by(orig.ident, cell_type2, group) %>%
    dplyr::summarize(across(all_of(target_pathways), mean, .names = "{col}"), .groups = 'drop') %>%
    dplyr::filter(!is.na(cell_type2)) # 剔除没有明确分类的细胞
  
  # 将宽表格转换为 ggplot2 喜欢的长表格格式（方便一键切分 4 个通路）
  melt_cell_df <- melt(
    df_cell_bulk, 
    id.vars = c("orig.ident", "cell_type2", "group"),
    measure.vars = target_pathways,
    variable.name = "Pathway",
    value.name = "AUC_Score"
  )
  
  # 净化通路名字（去掉 HALLMARK_ 并将下划线换为空格）
  melt_cell_df$Pathway <- gsub("HALLMARK_", "", melt_cell_df$Pathway)
  melt_cell_df$Pathway <- gsub("_", " ", melt_cell_df$Pathway)
  
  # ==============================================================================
  # 2. 绘制按细胞类型分组的硬朗科研风箱线图
  # ==============================================================================
  # 调用您之前定义好的专属细胞类型颜色映射向量（cluster_colors）
  # 确保它的名字和 melt_cell_df$cell_type2 里的细胞类型完全对齐
  
  p_cell_type_box <- ggplot(melt_cell_df, aes(x = cell_type2, y = AUC_Score, fill = cell_type2)) +
    
    # 1. 绘制基础箱线（收窄箱体，更精致）
    geom_boxplot(alpha = 0.85, width = 0.5, outlier.shape = NA, size = 0.5) +
    
    # 2. 叠加病人真实的伪 bulk 散点，展现样本间的异质性
    geom_jitter(color = "#333333", size = 1.0, alpha = 0.6, width = 0.15) +
    
    # 3. 【核心物理分块】一张图分成 4 个面板，分别代表 4 个通路，公用一套 Y 轴坐标系
    facet_wrap(~ Pathway, scales = "free_y", ncol = 2) +
    
    # 4. 自动添加全局组间显著性统计（Kruskal-Wallis 检验，看不同细胞间是否有差异）
    stat_compare_means(label.y.npc = "top", label.x.npc = "left", size = 3, color = "#444444") +
    
    # 5. 配色与经典学术边框
    scale_fill_manual(values = cluster_colors) + # 注入您的专属细胞颜色
    theme_bw(base_family = "sans") +
    labs(title = "Pathway Activity Across Cell Types", x = NULL, y = "Mean AUC Score") +
    
    theme(
      # 大标题居中加粗
      plot.title = element_text(hjust = 0.5, face = "bold", size = 13, color = "#222222", margin = margin(b = 12)),
      
      # 绘图区黑色硬朗外框
      panel.border = element_rect(color = "#444444", fill = NA, size = 0.8),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.spacing = unit(4, "mm"), # 加大 4 个分块面板之间的物理间距，极具呼吸感
      
      # 坐标轴美化：细胞类型标签倾斜 45 度，防止过长重叠
      axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 9.5, color = "black"),
      axis.text.y = element_text(size = 9, color = "black"),
      axis.title.y = element_text(size = 11, face = "bold"),
      axis.ticks = element_line(color = "#444444", size = 0.5),
      
      # 顶部 4 个通路条带栏（Strip）润色：高级淡灰底 + 纯黑粗体字
      strip.background = element_rect(fill = "grey93", color = "#444444", size = 0.8),
      strip.text = element_text(face = "bold", size = 10, color = "black"),
      
      # 隐藏右侧多余图例（因为 X 轴已经写明了细胞类型）
      legend.position = "none"
    )
  
  # ==============================================================================
  # 3. 完美适配 Adobe Illustrator (AI) 的安全导出
  # ==============================================================================
  # 在 RStudio 中预览
  print(p_cell_type_box)
  
  # 导出为高质量矢量 PDF
  ggsave(
    filename = "Pathway_Activity_By_CellType.pdf", 
    plot = p_cell_type_box, 
    width = 8.5,     # 2x2 排版，横向留足空间给多个细胞类型
    height = 7.5,    # 纵向高度
    dpi = 600
  )
  
}
