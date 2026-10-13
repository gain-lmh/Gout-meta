library(qs)
library(ggplot2)
library(DESeq2)
library(Seurat)
library(tidyverse)
library(dplyr)

scRNA_Monocyte2 <- qread("/data3/home/hu/Tongfeng/new_Result/Fig3/scRNA_Monocyte.qs")
{
  library(Seurat)
  library(sctransform)
  library(glmGamPoi)
  
  # ==================== 1. 提取目标细胞并“彻底洗干净” ====================
  
  # 1.1 提取你想重聚类的细胞（这里以分群为 1 和 2 的细胞为例，你可以换成你的目标细胞或群）
  # 比如：sce_sub <- subset(sce, idents = "Monocyte") 或者具体的群
  sce_sub <- scRNA_Monocyte2
  
  # 1.2 【核心绝招】强行把分样本的图层重新合并，并将默认 Assay 换回原始 RNA
  # 这一步能彻底清洗掉之前错乱的 SCT、PCA 和 UMAP 等旧信息，退回到最干净的初始计数矩阵
  DefaultAssay(sce_sub) <- "RNA"
  sce_sub <- JoinLayers(sce_sub)
  
  # 1.3 顺手删掉不需要的旧低维空间和旧 Assay，让对象变成一个纯净的“新对象”
  sce_sub[["SCT"]] <- NULL
  sce_sub@reductions <- list()
  
  # ==================== 2. 重新按 v5 标准分批次 ====================
  
  # 2.1 重新计算这些细胞的线粒体比例
  sce_sub[["percent.mt"]] <- PercentageFeatureSet(sce_sub, assay = "RNA", pattern = "^MT-")
  
  # 2.2 【关键步骤】将洗干净后的数据，重新按样本(orig.ident)切分成干净的 v5 图层
  sce_sub <- split(sce_sub, f = sce_sub$orig.ident)
  
  # ==================== 3. 重新运行高级标准化与 PCA ====================
  
  # 3.1 重新进行单独的 SCT 标准化
  sce_sub <- SCTransform(sce_sub, assay = "RNA", vars.to.regress = "percent.mt", verbose = FALSE)
  
  # 3.2 运行重聚类所需的全新 PCA
  sce_sub <- RunPCA(sce_sub, assay = "SCT", verbose = FALSE)
  
  # ==================== 4. 重新进行 CCA 整合去批次（此时绝对不会报错） ====================
  
  # 此时 SCT 的图层是基于最新提取的细胞全新生成的，没有任何历史残留 Bug，可以完美运行
  # ==================== 4. 重新进行 CCA 整合去批次（加入核心修正参数） ====================
  sce_sub <- IntegrateLayers(
    object = sce_sub, 
    method = CCAIntegration,          
    orig.reduction = "pca",           
    new.reduction = "integrated.cca", 
    assay = "SCT",                    
    normalization.method = "SCT",     # 【核心修正】必须死死加上这一行！告诉 v5 这是 SCT 数据
    verbose = FALSE
  )
  
  # ==================== 5. 下游全新的降维与分群 ====================
  
  # 基于重聚类整合后的新低维空间进行降维
  sce_sub <- RunUMAP(sce_sub, reduction = "integrated.cca", dims = 1:10, verbose = FALSE)
  sce_sub <- RunTSNE(sce_sub, reduction = "integrated.cca", reduction.name = "tsne", dims = 1:10, verbose = FALSE)
  
  # 重新构建近邻并寻找更细致的亚群（重聚类）
  sce_sub <- FindNeighbors(sce_sub, reduction = "integrated.cca", dims = 1:10, verbose = FALSE)
  sce_sub <- FindClusters(sce_sub, verbose = FALSE, resolution = 0.4) # 可以适当调小或调大分辨率查看细分群
  
  # ==================== 6. 重聚类后的差异表达分析 ====================
  
  # 重新合并图层以计算新亚群的 Marker 基因
  sce_sub <- JoinLayers(sce_sub) 
  sce_sub <- PrepSCTFindMarkers(sce_sub, assay = "SCT", verbose = TRUE)
  
  sub_markers <- FindAllMarkers(
    object = sce_sub, 
    test.use = "wilcox",
    only.pos = TRUE,
    logfc.threshold = 0.25
  )
  
  # 检查组别是否融合良好
  p1 <- DimPlot(sce_sub, reduction = "umap", group.by = "group")
  # 检查分簇情况
  p2 <- DimPlot(sce_sub, reduction = "umap", label = TRUE)
  
  library(patchwork)
  p1 + p2
  
  library(dplyr)
  
  # 1. 筛选每群前 10 个最显著的 Marker 基因
  # 依据：先按细胞群(cluster)分组，再按 p_val_adj 升序（越小越显著）或 avg_log2FC 降序（越大越特异）排序
  top10_markers <- all_markers %>%
    group_by(cluster) %>%
    slice_max(order_by = avg_log2FC, n = 10) # 如果你的列名是 avg_logFC，请把 2 去掉
  
  # 2. 在控制台查看结果
  print(top10_markers)
  
  # 3. 提取出这前 10 基因的纯名字列表（去重后用于后续画图）
  top10_genes <- unique(top10_markers$gene)
  
  # ==================== 1. 挑选最能代表四大分类的经典 Marker 基因 ====================
  # 基于你之前运行出来的 Top 基因以及髓系免疫学经典标志物进行精选：
  classic_markers <- c(
    # 1. 经典单核 (Cluster 0, 1): 基础强表达
    "CD14","LGALS2", "THBS1", "GPR183", "MMP9",
    # 2. 过渡单核 (Cluster 2, 3, 7): 基质调节与补体
    "PARP2", "CLU", "HTRA1", "C1QB",
    # 3. 炎症单核 (Cluster 4, 5, 6): 急性趋化、细胞因子与干扰素响应
    "CXCL8", "CXCL1", "CCL3", "CCL4", "RSAD2", "IFIT1", "ISG15",
    # 4. 非经典单核 (Cluster 8): 强抗原递呈特征
    "CD1C", "CLEC10A", "HLA-DQA1", "HLA-DPB1",
    "FCGR3A","HLA-DRA"
  )
  
  #####查看这些基因在不同簇中的表达
  library(Seurat)
  library(ggplot2)
  
  # 1. 确保当前激活的分群是 0-10 的原始分群
  Idents(sce_sub) <- "SCT_snn_res.0.4"
  
  # 2. 绘制经典 Marker 基因在 0-10 簇中的气泡图
  p_dot <- DotPlot(
    object = sce_sub, 
    features = classic_markers, 
    assay = "SCT"
  ) + 
    RotatedAxis() +              # 将横坐标基因名旋转 45 度，防止重叠
    coord_flip() +               # 强烈推荐：颠倒横纵坐标（纵轴为基因，横轴为簇），更符合文献阅读习惯
    theme(
      axis.text.x = element_text(size = 12, face = "bold"),
      axis.text.y = element_text(size = 10, face = "bold")
    )
  
  # 3. 展现图表
  print(p_dot)
  
  
  library(Seurat)
  
  # ==================== 1. 暴力过滤低质量细胞 ====================
  # 切换到你原本的 0-10 分群结果
  Idents(sce_sub) <- "SCT_snn_res.0.4"
  # 只提取高质量的单核细胞群体（彻底抛弃 1, 9, 10）
  sce_pure <- subset(sce_sub, idents = c("0", "2", "3", "4", "5", "6", "7", "8"))
  
  # ==================== 2. 核心：直接复用旧的 CCA 空间进行重聚类 ====================
  # 【关键点】我们不重新算 CCA，但必须重新运行以下三步，让画面重新铺满：
  
  # 2.1 重新基于现有的 integrated.cca 空间计算近邻图
  sce_pure <- FindNeighbors(sce_pure, reduction = "integrated.cca", dims = 1:12, verbose = FALSE)
  
  # 2.2 重新寻找新亚群（因为细胞变少了，它们会重新凝聚并获得从 0 开始的新编号）
  # 这里可以保持 resolution = 0.4 
  sce_pure <- FindClusters(sce_pure, verbose = FALSE, resolution = 0.4)
  
  # 2.3 重新运行 UMAP 投影（这样剩下的细胞会自动拉伸、填满整张图，没有任何空白黑洞）
  sce_pure <- RunUMAP(sce_pure, reduction = "integrated.cca", dims = 1:12, verbose = FALSE)
  
  # ==================== 3. 可视化查看全新结果 ====================
  
  # 此时你会得到一张全新的、没有任何空缺、细胞分布异常丝滑漂亮的 UMAP 图
  p_new <- DimPlot(sce_pure, reduction = "umap", label = TRUE, repel = TRUE) + 
    ggtitle("Pure Monocytes: Re-clustered on Existing CCA Space")
  
  print(p_new)
  
  # 1. 确保当前激活的分群是 0-10 的原始分群
  Idents(sce_pure) <- "SCT_snn_res.0.4"
  
  # 2. 绘制经典 Marker 基因在 0-10 簇中的气泡图
  p_dot <- DotPlot(
    object = sce_pure, 
    features = classic_markers, 
    assay = "SCT"
  ) + 
    RotatedAxis() +              # 将横坐标基因名旋转 45 度，防止重叠
    coord_flip() +               # 强烈推荐：颠倒横纵坐标（纵轴为基因，横轴为簇），更符合文献阅读习惯
    theme(
      axis.text.x = element_text(size = 12, face = "bold"),
      axis.text.y = element_text(size = 10, face = "bold")
    )
  
  # 3. 展现图表
  print(p_dot)

  library(Seurat)
  
  # 1. 确保激活的是最新重新聚类出来的分群结果
  # 如果你不放心，可以直接指定为默认 Idents
  # Idents(sce_pure) <- "seurat_clusters" 
  
  # 2. 针对这全新的 11 个高质量群，建立严格的经典分类映射
  new_pure_annotation <- c(
    "0"  = "Classical_Monocytes",
    "1"  = "Classical_Monocytes",
    "2"  = "Classical_Monocytes",
    "3"  = "Inflammatory_Monocytes", # CXCL8+
    "4"  = "Inflammatory_Monocytes", # CCL3/4+
    "5"  = "Inflammatory_Monocytes", # ISG+ 干扰素
    "6"  = "Non_classical_Monocytes",# CD16+ 高表达
    "7"  = "Non_classical_Monocytes",# HLA+ 抗原递呈
    "8"  = "Intermediate_Monocytes", # CLU+ 过渡
    "9"  = "Classical_Monocytes",
    "10" = "Classical_Monocytes"
  )
  
  # 3. 智能替换标签并写入元数据
  sce_pure <- RenameIdents(sce_pure, new_pure_annotation)
  sce_pure$cell_type_classic <- Idents(sce_pure)
  
  # 4. 打印最终完美单核细胞亚型的细胞计数
  print("最终高质量单核细胞经典分类计数：")
  print(table(sce_pure$cell_type_classic))
  
  # 5. 重新绘制没有任何空白、没有任何低质量孤岛、铺满整个画面的完美 UMAP 图
  DimPlot(sce_pure, reduction = "umap", label = TRUE, repel = TRUE, label.size = 4) + 
    ggtitle("Final Purified Monocyte Classical Lineages")
  
  
  #####基因marker表达###########
  sce_pure$main_celltype <- factor(sce_pure$cell_type_classic,
                                   levels = c("Classical_Monocytes","Intermediate_Monocytes",
                                              "Inflammatory_Monocytes","Non_classical_Monocytes"))
  
  cluster_colors <- c("#4DBBD5", "#989EBE", "#B44D51","#3C5488FF")   
  # 经典框框
  p1 <- CellDimPlot(sce_pure, group.by = "main_celltype", alpha.highlight = 0.5,
                    reduction = "UMAP", show_stat = F,label_insitu = T,pt.size = 0.01,
                    palcolor = cluster_colors,raster = FALSE,label = T );p1
  
  ggsave(filename = "/data3/home/hu/Tongfeng/new_Result/Fig4/New result/mono_seurat.pdf",plot = p1,height = 5,width = 7)
  
  classic_markers <- c(
    # 1. 经典单核 (Cluster 0, 1): 基础强表达
    "CD14","LGALS2", "THBS1",  
    # 2. 过渡单核 (Cluster 2, 3, 7): 基质调节与补体
    "PARP2", "CLU", 
    # 3. 炎症单核 (Cluster 4, 5, 6): 急性趋化、细胞因子与干扰素响应
    "CXCL8", "CXCL1", "CCL3", "CCL4", "RSAD2", "IFIT1", "ISG15",
    # 4. 非经典单核 (Cluster 8): 强抗原递呈特征
    "CD1C", "CLEC10A", "HLA-DQA1", "HLA-DPB1",
    "FCGR3A","HLA-DRA"
  )
  # 基因表达均值热图
  ht1 <- GroupHeatmap(sce_pure,
                      features = classic_markers,
                      group.by = c("main_celltype"),
                      show_row_names = TRUE,
                      row_names_side = 'left',
                      #legend_title = 'z-score',
                      column_title = ''
                      
  )
  ht1$plot
  
  # 面板大小固定
  pdf("/data3/home/hu/Tongfeng/new_Result/Fig4/New result/mono_marker.pdf",height = 5,width = 7)
  panel_fix(ht1$plot, height = 4, width = 6)
  dev.off()
  
  
  library(Seurat)
  library(scCustomize)
  library(patchwork)
  library(ggplot2)
  
  # 1. 精选四大核心经典 Marker 基因
  gene <- c("CD14", "CXCL1", "FCGR3A", "CLU")
  
  # 2. 统一定制顶级渐变配色（从静息浅灰，到过渡浅蓝，再到强表达深海蓝）
  my_blue_palette <- colorRampPalette(c("#EAEAEA", "#08306B"))(100)
  
  # 3. 循环高效建图
  plots <- list()
  for (i in 1:length(gene)){
    
    plots[[i]] <- FeaturePlot_scCustom(
      seurat_object = sce_pure, 
      features = gene[i],
      raster = TRUE,                 # 【体积控制核心 1】开启光栅化
      raster.dpi = c(300, 300),      # 【体积控制核心 2】锁定 300 DPI 分辨率，让 PDF 瞬间缩小 90%
      pt.size = 0.6,                 # 适度调大点的大小，弥补光栅化后的颗粒感
      colors_use = my_blue_palette
    ) + 
      NoAxes() +                     # 抹去刻度轴
      labs(title = gene[i]) +        # 将基因名作为大标题
      theme(
        plot.title = element_text(face = "bold.italic", size = 14, hjust = 0.5), # 基因名粗斜体居中
        legend.position = "right",   # 暂时让图例靠右，方便后续 patchwork 自动捕捉合并
        legend.direction = "vertical"
      )
  }
  
  # 4. 【核心美化与缩量】用 patchwork 强行合并图例并排版
  # guides = "collect" 会自动识别 4 个图例完全一样，从而把它们合并成唯一的一个！
  final_p <- wrap_plots(plots, ncol = 2) + 
    plot_layout(guides = "collect") & 
    theme(
      legend.position = "bottom",       # 将合并后的唯一图例优雅地锁在整个大图的正下方
      legend.box = "horizontal",
      legend.title = element_blank(),   # 隐去多余的 "Expression" 字样，保持严谨简洁
      legend.text = element_text(size = 9, face = "bold"),
      legend.key.width = unit(0.8, "cm"), # 控制横向渐变条的宽度
      legend.key.height = unit(0.2, "cm") # 控制横向渐变条的高度
    )
  
  # 5. 【高清保存】以标准小巧尺寸（6x6.5英寸）安全导出
  # 这样导出的 PDF 文件一般只有 300KB 左右，图例绝对不会错位，且可以直接放大无限倍
  ggsave(
    filename = "/data3/home/hu/Tongfeng/new_Result/Fig1/featureplot_subset4.pdf", 
    plot = final_p, 
    width = 6, 
    height = 6.5,  # 稍微拉高 0.5 英寸给底部的横向单图例留出完美空间
    dpi = 300
  )
  
  
  library(Seurat)
  library(SCP)
  # 细胞比例绘制
  pdf("/data3/home/hu/Tongfeng/new_Result/Fig4/New result/cell_radio.pdf",width = 5,height = 5)
  CellStatPlot(sce_pure,
               stat.by = "main_celltype",palcolor = cluster_colors,
               group.by = "group",plot_type = "trend",
               label = TRUE) %>% panel_fix(height = 4, width = 3)
  dev.off()
  
  
  
}


######绘制在不同单核细胞类型中的代谢情况
{
  wgcna_score <- readRDS("/data3/home/hu/Tongfeng/new_Result/Fig4/New result/WGCNA_module_score.RDS")
  sub_wgcna <- wgcna_score[colnames(sce_pure),]
  sce_pure@meta.data <- cbind(sce_pure@meta.data,sub_wgcna[,c(4:9)])
  qsave(sce_pure,"/data3/home/hu/Tongfeng/new_Result/Fig4/New result/sce_sub.qs")
  
  table(sce_pure$group,sce_pure$main_celltype)
  aggregate(sce_pure$WGCNA,list(sce_pure$group),mean)
  aggregate(sce_pure$WGCNA,list(sce_pure$main_celltype),mean)
  
  #####查看不同细胞类型代谢通路的表达
  meta.data <- sce_pure@meta.data
  
  library(tidyverse)
  library(ggpubr)
  
  # 1. 整理数据：选择需要的列并转换为长格式
  plot_df <- meta.data %>%
    dplyr::select(main_celltype, WGCNA, `Urea cycle`, 
                  `Arachidonic acid metabolism`, `Eicosanoid metabolism`, 
                  `Nitrogen metabolism`, `Fatty acid biosynthesis`) %>%
    # 转换为长格式以便绘图
    pivot_longer(cols = -main_celltype, names_to = "Pathway", values_to = "Score") %>%
    # 设置细胞类型顺序
    mutate(main_celltype = factor(main_celltype, 
                                  levels = c("Classical_Monocytes","Intermediate_Monocytes",
                                             "Inflammatory_Monocytes","Non_classical_Monocytes")))
  
  
  library(tidyverse)
  library(ggpubr)
  library(ggsci) # 引入学术配色包
  
  # 1. 数据预处理
  plot_df <- meta.data %>%
    dplyr::select(main_celltype, WGCNA, `Urea cycle`, 
                  `Arachidonic acid metabolism`, `Eicosanoid metabolism`, 
                  `Nitrogen metabolism`, `Fatty acid biosynthesis`) %>%
    pivot_longer(cols = -main_celltype, names_to = "Pathway", values_to = "Score") %>%
    mutate(
      main_celltype = factor(main_celltype, levels = c("Classical_Monocytes","Intermediate_Monocytes",
                                                       "Inflammatory_Monocytes","Non_classical_Monocytes")),
      # 将 WGCNA_score 排在第一位，突出核心结果
      Pathway = factor(Pathway, levels = c("WGCNA", "Arachidonic acid metabolism", 
                                           "Eicosanoid metabolism", "Fatty acid biosynthesis", 
                                           "Nitrogen metabolism", "Urea cycle"))
    )
  
  library(ggplot2)
  library(ggpubr)
  
  # 1. 锁定您的专属 4 色方案
  cluster_colors <- c("#4DBBD5", "#989EBE", "#B44D51", "#3C5488FF")
  
  # 2. 绘制顶级艺术质感学术图
  p_polished <- ggplot(plot_df, aes(x = main_celltype, y = Score, fill = main_celltype)) +
    
    # 【微调 1】柔化小提琴图：增加细黑边（linewidth=0.2），让半透明的小提琴轮廓更精致，不再模糊
    geom_violin(
      scale = "width", 
      alpha = 0.55,          # 稍微提高一点饱和度
      color = "#333333",     # 给予极细的深灰色边缘线，增加立体感
      linewidth = 0.2,
      trim = TRUE,
      position = position_dodge(0.8)
    ) +
    
    # 【微调 2】精致内嵌箱线：加粗中位线，箱体边框使用更高级的深黑色
    geom_boxplot(
      width = 0.14,          # 保持优雅的极窄设计
      color = "#111111",     # 稍深一点的黑框
      linewidth = 0.5,
      outlier.shape = NA, 
      position = position_dodge(0.8)
    ) +
    
    # 3. 分面逻辑（保持完美形状）
    facet_wrap(~ Pathway, scales = "free_y", nrow = 2) + 
    
    # 4. 映射配色
    scale_fill_manual(values = cluster_colors) + 
    
    # 5. 【核心修复 3】完美居中 P 值星号
    # 强制让统计星号挪到 4 个小提琴的正中央（x = 2.5），并且上移一些，看起来非常工整
    stat_compare_means(
      label = "p.signif", 
      label.x = 2.5,         # 【核心修改】1.5 改为 2.5，死死居中在 4 组正上方
      method = "kruskal.test",
      size = 4.5,            # 稍微放大一点点星号
      color = "black"
    ) +
    
    # 6. 【终极细节美化】打造轻奢感科研主题
    theme_bw(base_size = 13) + 
    theme(
      # 彻底抹去多余网格线
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(),
      panel.border = element_rect(linewidth = 0.8, color = "#222222"), # 更严谨的框线
      
      # 【高级感核心】分面标题栏（Strip）注入柔和的高级灰，不再是一片惨白，文字加粗
      strip.background = element_rect(fill = "#EAEAEA", color = "#222222", linewidth = 0.8),
      strip.text = element_text(face = "bold.italic", size = 11, color = "#111111"), # 通路名用粗斜体，更符合学术规范
      
      # 彻底隐藏多余的横坐标标签（因为顶部的横向 Legend 已经说得很清楚了）
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.x = element_blank(),
      
      # 纵坐标文字加粗变黑，可读性更强
      axis.text.y = element_text(color = "black", size = 10),
      axis.title.y = element_text(face = "bold", size = 12, color = "black", margin = margin(r = 12)),
      
      # 图例无框、无标题，横向居中排列
      legend.position = "top",
      legend.title = element_blank(),
      legend.text = element_text(size = 11, face = "bold", color = "#222222"),
      legend.background = element_blank(),
      legend.box.margin = margin(b = -5) # 稍微让图例往图表方向靠拢，减少上方留白
    ) +
    
    labs(
      y = "Pathway Activation Score (UCell)"
    ) +
    
    # 留出 18% 的上方空间，配合居中的星号，绝对不会让触须和星号撞车
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.18)))
  
  # 渲染查看
  print(p_polished)
  
  ggsave(plot = p_polished,filename = "/data3/home/hu/Tongfeng/new_Result/Fig4/New result/WGCNA_voin.pdf",
         width = 9,height = 5)
}


#####根据分数中位数划分高低分组
sce_pure$class <- ifelse(sce_pure$WGCNA > median(sce_pure$WGCNA),"High","Low")
library(ggplot2)
library(dplyr)
library(patchwork)

# 1. 提取元数据并计算每个细胞类型中 High/Low 的比例
ratio_df <- sce_pure@meta.data %>%
  group_by(cell_type_classic, class) %>%
  summarise(Count = n(), .groups = 'drop') %>%
  group_by(cell_type_classic) %>%
  mutate(Percentage = Count / sum(Count) * 100)

# 2. 绘制精美的百分比堆叠柱状图
p_ratio <- ggplot(ratio_df, aes(x = cell_type_classic, y = Percentage, fill = class)) +
  geom_bar(stat = "identity", width = 0.6, color = "black", linewidth = 0.4) +
  scale_fill_manual(values = c("High" = "#B44D51", "Low" = "#4DBBD5")) + # 红色代表 High，蓝色代表 Low
  theme_bw(base_size = 12) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, color = "black", face = "bold"),
    axis.text.y = element_text(color = "black"),
    legend.position = "top",
    legend.title = element_blank()
  ) +
  labs(y = "Cell Proportion (%)", x = "")

print(p_ratio)
ggsave(plot = p_ratio,filename = "/data3/home/hu/Tongfeng/new_Result/Fig4/New result/WGCNA_High_Low_p_ratio.pdf",
       width = 7,height = 5)


#####拟时序分析
{
  
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  #library(monocle3)
  library(Matrix)
  library(stringr)
  library(harmony)
  library(monocle)
  
  scRNA = qread("/data3/home/hu/Tongfeng/new_Result/Fig4/New result/sce_sub.qs")
  DefaultAssay(scRNA) <- "SCT"
  
  ##提取表型信息--细胞信息(建议载入细胞的聚类或者细胞类型鉴定信息、实验条件等信息)
  #expr_matrix <- as(as.matrix(scRNA@assays$RNA@layers$counts), 'sparseMatrix')
  expr_matrix <- LayerData(scRNA, assay = "SCT", layer = "counts")
  
  ##提取表型信息到p_data(phenotype_data)里面 
  p_data <- scRNA@meta.data[,c(1,9,13:21)]
  p_data$celltype <- scRNA$main_celltype  ##整合每个细胞的细胞鉴定信息到p_data里面。如果已经添加则不必重复添加
  ##提取基因信息 如生物类型、gc含量等
  f_data <- data.frame(gene_short_name = row.names(scRNA),row.names = row.names(scRNA))
  ##expr_matrix的行数与f_data的行数相同(gene number), expr_matrix的列数与p_data的行数相同(cell number)
  
  #构建CDS对象
  pd <- new('AnnotatedDataFrame', data = p_data) 
  fd <- new('AnnotatedDataFrame', data = f_data)
  #将p_data和f_data从data.frame转换AnnotatedDataFrame对象。
  cds <- newCellDataSet(expr_matrix,
                        phenoData = pd,
                        featureData = fd,
                        lowerDetectionLimit = 0.5,
                        expressionFamily = negbinomial.size())
  
  #估计size factor和离散度
  cds <- estimateSizeFactors(cds)
  cds <- estimateDispersions(cds)
  
  
  
  ##使用seurat选择的高变基因⚠️
  scRNA <- FindVariableFeatures(scRNA, selection.method = "vst", nfeatures = 1000)
  express_genes <- VariableFeatures(scRNA)
  #cds <- setOrderingFilter(cds, express_genes)
  #plot_ordering_genes(cds)
  
  #这一步输入的expressed_genes来自于步骤4。
  #⚠️⚠️后续分析使用的是该方法
  #也可输入seurat筛选出的高变基因：expressed_genes <- VariableFeatures(pbmc) 
  diff <- differentialGeneTest(cds[express_genes[1:500],],fullModelFormulaStr="~celltype",cores=1) 
  #~后面是表示对谁做差异分析的变量，理论上可以为p_data的任意列名
  head(diff)
  
  ##差异表达基因作为轨迹构建的基因,差异基因的选择标准是qval<0.01,decreasing=F表示按数值增加排序
  deg <- subset(diff, qval < 0.001) #选出2829个基因
  deg <- deg[order(deg$qval,decreasing=F),]
  head(deg)
  
  ##差异基因的结果文件保存
  #write.table(deg,file="1.1train.monocle.DEG.xls",col.names=T,row.names=F,sep="\t",quote=F)
  
  ## 轨迹构建基因可视化
  ordergene <- rownames(deg) 
  cds <- setOrderingFilter(cds, ordergene)  
  #这一步是很重要的，在我们得到想要的基因列表后，我们需要使用setOrderingFilter将它嵌入cds对象，后续的一系列操作都要依赖于这个list。
  #setOrderingFilter之后，这些基因被储存在cds@featureData@data[["use_for_ordering"]]，可以通过table(cds@featureData@data[["use_for_ordering"]])查看
  #pdf("train.ordergenes.pdf")
  plot_ordering_genes(cds)
  dev.off()
  #出的图黑色的点表示用来构建轨迹的差异基因，灰色表示背景基因。红色的线是根据第2步计算的基因表达大小和离散度分布的趋势(可以看到，找到的基因属于离散度比较高的基因)
  
  
  ## 降维
  cds <- reduceDimension(cds, 
                         max_components = 2, 
                         method = 'DDRTree', 
                         num_dim = 6,           # 尝试减小维度（如 5 或 6），有助于简化主干
                         ncenter = 60)        # 适当设置节点数，节点越少轨迹越简单
  
  # 重新排序
  cds <- orderCells(cds)
  cds <- orderCells(cds, root_state = 2)
  
  cds$class <- ifelse(cds$WGCNA > median(cds$WGCNA),"High","Low")
  
  setwd("/data3/home/hu/Tongfeng/new_Result/Fig4")
  # 以pseudotime值上色
  pdf("1.4train.monocle.pseudotime.pdf",width = 6.5,height = 4)
  plot_cell_trajectory(cds,color_by="Pseudotime", size=1,show_backbone=TRUE) +
    scale_color_gradientn(
      colors = c("#2E3192","#EBEABE","#DE9EBE")) + theme(legend.position = "right") 
  dev.off()
  
  # 以细胞类型上色
  pdf("1.5train.monocle.celltype.pdf",width = 7.5,height = 4)
  plot_cell_trajectory(cds,color_by="celltype", size=0.1,show_backbone=TRUE) +
    scale_color_manual(breaks = c("Classical_Monocytes","Intermediate_Monocytes",
                                  "Inflammatory_Monocytes","Non_classical_Monocytes"), 
                       values=c("#4DBBD5", "#989EBE", "#B44D51", "#3C5488FF"))+
    theme(legend.position = "right")
  dev.off()
  
  
  
  pdf("1.6train.monocle.stage.pdf",width = 7.5,height = 4)
  plot_cell_trajectory(cds,color_by="group", size=0.1,show_backbone=TRUE) +
    scale_color_manual(breaks = c("Control", "Remission","Disease"), values=c('#a6c896','#f2b9b2',"#FC8D62")) + theme(legend.position = "right")
  dev.off()
  
  
  pdf("1.8train.monocle.state.pdf",width = 7.5,height = 4)
  plot_cell_trajectory(cds,color_by="State", size=0.1,show_backbone=TRUE) +
    scale_color_manual(breaks = c("1", "2", "3","4","5"), values=c("#0078a3", "#FC8D62", "#66C2A5","#989DDE","#BDDD51")) + theme(legend.position = "right")
  dev.off()
  
  pdf("1.8train.monocle.state.pdf",width = 7.5,height = 4)
  plot_cell_trajectory(cds,color_by="class", size=0.1,show_backbone=TRUE) +
    scale_color_manual(breaks = c("High", "Low"), values=c("#f2b9b2", "#BDDD51")) + theme(legend.position = "right")
  dev.off()
  
  
  library(Monocle)
  library(ggplot2)
  
  # 1. 挑选您想展示的 4 个目标基因（请替换为您自己数据中真实存在的基因名）
  my_select_genes <- c("CD5L", "CD68", "SPP1", "VCAM1")
  
  # 2. 过滤掉数据集中不存在的基因，防止报错
  available_genes <- my_select_genes[my_select_genes %in% row.names(cds)]
  
  # 3. 打开 PDF 画布，由于是纵向单列排布，建议将高度(height)设长一些
  pdf("point1_genes_branched_pseudotime.pdf", width = 5, height = 7)
  
  # 4. 核心绘图函数
  library(Monocle)
  library(ggplot2)
  
  # 1. 挑选您想要展示的 4 个标志性基因
  target_genes <- c("CD5L", "CD68", "SPP1", "VCAM1")
  
  # 2. 从您的单细胞对象 (cds) 中筛选出这些基因
  # 注意：cds 必须是已经运行完 reduceDimension(..., method = 'DDRTree') 和 orderCells() 的对象
  cds_subset <- cds[available_genes, ]
  
  # 3. 核心绘图函数
  plot_genes_in_pseudotime(
    cds_subset, 
    color_by = "State",            # 按照分化状态（或你的分组/细胞类型）上色
    nrow = 4,                      # 强行排成 4 行（即单列纵向排布）
    ncol = 1,                      
    cell_size = 1                 # 调整散点的大小
   # trendline_type = "gam"         # 拟合曲线类型，可以是 "deg" 或 "gam"
  ) + 
    # 4. 叠加 ggplot2 样式精修，使其完美对齐原图
    scale_y_log10() +                # 对应原图的非线性 y 轴表达量
    theme_classic() +                # 变成干净的白色背景
    theme(
      legend.position = "none",       # 隐藏右侧多余的图例
      strip.background = element_blank(), # 去掉基因名字外面的灰色方框
      strip.text = element_text(size = 14, face = "plain") # 调整基因字体
    )
  
  
  dev.off()
  
  
  ####为每个分支绘制饼图 ----细胞类型
  {
    library(tidyverse)
    
    # ==================== 1. 【核心修复】使用 tidyverse 绝对精准计算比例 ====================
    cds_metadata <- cds@phenoData@data
    cds_metadata$point <- cds_metadata$State
    
    # 抛弃繁琐的两次 aggregate，用一行流全自动搞定计数与组内百分比
    cell_number <- cds_metadata %>%
      filter(!is.na(point) & !is.na(celltype)) %>% # 剔除空值增强健壮性
      group_by(point, celltype) %>%
      summarise(cell_number = n(), .groups = "drop") %>% # 自动统计每个 Point 下各细胞类型的数量
      group_by(point) %>%
      mutate(
        point_total = sum(cell_number),
        cell_radio = cell_number / point_total,          # 绝对精准的组内比例
        label_text = paste0(round(cell_radio * 100, 1), "%") # 生成优雅的百分比标签（如 35.4%）
      ) %>%
      ungroup() %>%
      rename(cell_anno = celltype, Point = point) # 保持你后续代码的变量名一致
    
    # ==================== 2. 锁定您的专属 4 色方案 ====================
    colour <- c("#4DBBD5", "#989EBE", "#B44D51", "#3C5488FF")
    
    # ==================== 3. 绘制顶级学术分面饼图 ====================
    p_pie <- ggplot(cell_number, aes(x = "", y = cell_radio, fill = cell_anno)) + 
      # 绘制饼图扇区：加入极细的白边或淡黑边，比粗黑边高级得多
      geom_col(color = "white", linewidth = 0.4, show.legend = TRUE) + 
      
      # 【核心美化】在扇区中央极其工整地叠加上百分比标签
      geom_text(
        aes(label = label_text), 
        position = position_stack(vjust = 0.5), # 完美的极坐标堆叠居中
        size = 3.8,                             # 稍微收紧字体
        fontface = "bold", 
        color = "white"                         # 在深色底色上用白色字，视觉对比度极佳
      ) + 
      
      # 强制使用你的专属 4 种单核细胞颜色
      scale_fill_manual(values = colour) +
      
      # 转换为极坐标饼图样式
      coord_polar("y", start = 0) + 
      
      # 动态分面：根据 Point 的数量自动决定排版（自动转行，防止硬塞在一行导致压扁）
      facet_wrap(~ Point, nrow = 1) + 
      
      # 极简背景
      theme_void(base_size = 14) + 
      theme(
        # 分面顶部条框美化：粗斜体增加文章质感
        strip.text = element_text(size = 15, face = "bold.italic", margin = margin(b = 10)),
        
        # 【核心美化】图例重置：置于右侧，垂直居中排列，不再悬空错位
        legend.position = "right",
        legend.title = element_blank(),
        legend.text = element_text(size = 11, face = "bold"),
        legend.background = element_blank(),
        
        # 调整四周留白，防止 coord_polar 裁剪边缘文字
        plot.margin = margin(t = 10, r = 10, b = 10, l = 10)
      )
    
    p_pie
    # ==================== 4. 安全、小巧、高清保存 ====================
    # 根据 Point 的数量动态调整宽度。如果是 4 个 Point，长宽比 14x4.5 是极美的黄金比例
    ggsave(
      plot = p_pie,
      filename = "/data3/home/hu/Tongfeng/new_Result/Fig4/New result/pie_plot_point.pdf",
      width = 14, 
      height = 4.5,
      dpi = 300
    )
    
}


#### 给每个细胞进行纤维化、炎症和衰老基因集打分
{
  ###纤维化基因集
  Fib_gene <- c(
    "AKR1B1", "AL035706.1", "ARL4C", "ARRDC2", "BTG2", "COL4A1", "COL4A2", 
    "CYTOR", "EHD4", "ERVW-1", "FTOP1", "GSN", "HTR2A", "IER5", "IL27RA", 
    "INMT", "LINC01725", "LPAL2", "NFKB2", "PKM", "S100A4", "SOX5", 
    "TPM4", "TRBC2", "VIM", "XYLB"
  )
  
  ###炎症相关基因集
  library(msigdbr)
  msig_h <- msigdbr(species = "Homo sapiens", category = "H")
  hall_inflam <- msig_h %>%
    filter(gs_name == "HALLMARK_INFLAMMATORY_RESPONSE")
  Inf_gene <- unique(hall_inflam$gene_symbol)
  
  
  ###衰老相关
  hall_inflam <- msig_h %>%
    filter(gs_name == "HALLMARK_TNFA_SIGNALING_VIA_NFKB")
  SIGNALING_VIA_NFKB_gene <- unique(hall_inflam$gene_symbol)
  
  DNA_damage_gene <- read.csv("/data3/home/hu/Tongfeng/new_Result/Fig4/DDR.csv") 
  
  SenMayo <- read.csv("/data3/home/hu/Tongfeng/new_Result/Fig4/SenMayo.csv")
  
  
  s.genes <- cc.genes$s.genes
  g2m.genes <- cc.genes$g2m.genes
  Cell_cycle_gene <- c(s.genes)
  
  gene_list <- list(Fib_gene,Inf_gene,Cell_cycle_gene,SIGNALING_VIA_NFKB_gene,
                    DNA_damage_gene$Gene.symbol)
  names(gene_list) <- c("Fib_gene","Inf_gene","Cell_cycle_gene","SIGNALING_VIA_NFKB_gene",
                        "DNA_damage_gene")
  
  ####### AUCell
  library(AUCell)
  library(clusterProfiler) 
  library(ggsci)
  ###纤维化基因集
  exprMat <- GetAssayData(scRNA, slot = "data")
  cells_rankings <- AUCell_buildRankings(exprMat)  
  
  ###纤维化基因集
  geneSets = gene_list
  
  cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.1)
  geneSet <- c("Fib_gene","Inf_gene","Cell_cycle_gene","SIGNALING_VIA_NFKB_gene",
               "DNA_damage_gene")
  
  auc_matrix <- getAUC(cells_AUC)[geneSet, ]
  
  metadata2 <- cds@phenoData@data
  
  metadata3 <- cbind(metadata2, t(auc_matrix))
  
  cor.test(metadata3$WGCNA,metadata3$Inf_gene)
  {
    library(ggplot2)
    library(ggpubr)
    
    # 高级学术美化版
    p1 <- ggplot(metadata3, aes(x = WGCNA, y = Inf_gene, color = main_celltype)) +
      # 1. 缩小点半径(size=0.5)，微调透明度，大幅缓解2万细胞的重叠死黑感
      geom_point(alpha = 0.15, size = 0.5, shape = 16) +
      # 2. 将拟合线统一设为深灰色/黑色，置信区间设为浅灰，让趋势更突出
      geom_smooth(method = "lm", se = TRUE, color = "#222222", size = 1, fill = "grey80", alpha = 0.4) +
      # 3. 分面展示
      facet_wrap(~main_celltype) +
      # 4. 美化相关系数标签，改用斜体，调整精确位置
      stat_cor(
        method = "pearson", 
        label.x.npc = 0.05, 
        label.y.npc = 0.95, 
        size = 4.5, 
        fontface = "italic",
        color = "black"
      ) +
      # 5. 严格同步你前两张图的红、蓝、绿高级学术配色
      scale_color_manual(values = c(
        "Classical Monocytes" = "#E64B35FF",     # 经典型：珊瑚红
        "Intermediate Monocytes" = "#4DBBD5FF",  # 中间型：天空蓝
        "Non-classical Monocytes" = "#00A087FF"  # 非经典：孔雀绿
      )) +
      # 6. 补全符合英文文献规范的完整坐标轴标签
      labs(
        x = "WGCNA Module Activation Score",
        y = "Inflammatory Gene Expression"
      ) +
      # 7. 极致定制主题
      theme_bw(base_size = 12) + 
      theme(
        # 移除右侧冗余图例（分面标题已有细胞名）
        legend.position = "none",
        # 剔除背景网格线，画面更干净
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        # 坐标轴文字加黑
        axis.text = element_text(color = "black", size = 10),
        axis.title = element_text(face = "bold", size = 12),
        # 分面外框条美化：改用极其高雅的淡灰色背景和细边框
        strip.background = element_rect(fill = "#F5F5F5", color = "grey70", size = 0.6),
        strip.text = element_text(face = "bold", size = 11, color = "black"),
        # 绘图区外框微调
        panel.border = element_rect(color = "grey70", fill = NA, size = 0.6)
      )
    
    ggsave(plot = p1,filename = "clltype_type_module_Inla.pdf",width = 12,height = 5)
  }
  
  
  
  plot_cell_trajectory(cds,color_by="WGCNA", size=1,show_backbone=TRUE) + scale_color_gsea()
  {
    library(ggplot2)
    library(viridis)
    
    # 1. 生成原始的拟时序图（此时由于官方 Bug，出来可能是全黑的）
    p_trajectory <- plot_cell_trajectory(
      cds = cds,
      color_by = "WGCNA",
      size = 0.5,               # 适度收窄细胞点，增加精致感
      show_backbone = TRUE
    )
    
    # 2. 【核心解封 Bug】强行重构图层映射，把细胞点的颜色死死绑定到 WGCNA 变量上
    # Monocle2 的第 1 个图层（[[1]]）就是细胞点。我们要对它的 mapping 强制重新赋予 aes(color = WGCNA)
    p_trajectory$layers[[1]]$mapping <- aes(by = NULL, color = WGCNA)
    
    # 3. 【视觉升级】微调主干骨架线（Backbone）
    # 第 2 个图层（[[2]]）是分化主干线，我们把它加粗，让分化树更清晰立体
    p_trajectory$layers[[2]]$aes_params$size <- 1.1
    p_trajectory$layers[[2]]$aes_params$colour <- "#111111" # 纯黑色骨架线
    
    # 4. 注入高级学术配色与轻奢主题（抹去多余网格）
    p_final <- p_trajectory +
      
      # 使用经典的 magma 熔岩色（低分暗紫，高分暴亮，文章故事性极强）
      # 如果你更喜欢原本的 GSEA 蓝红配色，可以换成：scale_color_gsea(name = "WGCNA")
      scale_color_viridis_c(
        option = "magma", 
        direction = 1,
        name = "WGCNA\nScore" # 图例标题
      ) +
      
      # 极致严谨的单细胞学术主题
      theme_classic(base_size = 14) +
      theme(
        panel.grid = element_blank(), # 彻底移除方格网格线
        axis.line = element_line(linewidth = 0.8, color = "#333333"),
        axis.text = element_text(color = "black", size = 11),
        axis.title = element_text(face = "bold", size = 13, color = "black"),
        
        # 图例精致修饰
        legend.position = "right",
        legend.title = element_text(size = 11, face = "bold"),
        legend.text = element_text(size = 10),
        legend.background = element_blank(),
        legend.key.height = unit(1.1, "cm"),
        legend.key.width = unit(0.35, "cm")
      ) +
      labs(
        x = "Component 1",
        y = "Component 2"
      )
    
    # 5. 【高清、极小体积】保存为 PDF 
    ggsave(
      plot = p_final,
      filename = "/data3/home/hu/Tongfeng/new_Result/Fig1/WGCNA_Trajectory_Perfect.pdf",
      width = 7,
      height = 5.5,
      dpi = 300
    )
    
    
  }
  
  score_data <- metadata3[,c("WGCNA","Fib_gene","Inf_gene","SIGNALING_VIA_NFKB_gene",
                             "DNA_damage_gene","Cell_cycle_gene")]
  ####计算相关性
  library(Hmisc)
  cor_res <- rcorr(as.matrix(score_data), type = "pearson")
  
  cor_mat <- cor_res$r
  p_mat <- cor_res$P
  
  cor_mat <- cor_mat["WGCNA", -1, drop = FALSE]
  p_mat <- p_mat["WGCNA", -1, drop = FALSE]
  
  ##生成显著性标记
  sig_label <- matrix("", nrow=nrow(p_mat), ncol=ncol(p_mat))
  sig_label[p_mat < 0.05] <- "*"
  sig_label[p_mat < 0.01] <- "**"
  sig_label[p_mat < 0.001] <- "***"
  
  #合并相关系数和显著性
  display_mat <- matrix(paste0(sprintf("%.2f", cor_mat), sig_label),
                        nrow=nrow(cor_mat))
  
  rownames(display_mat) <- rownames(cor_mat)
  colnames(display_mat) <- colnames(cor_mat)
  
  
  library(pheatmap)
  pdf("D:\\博二上\\NPJ二修\\xinxiutu\\AGS_cor.pdf",height = 4,width = 7)
  pheatmap(cor_mat,
           cluster_rows = FALSE,
           cluster_cols = FALSE,
           color = colorRampPalette(c("#3B4CC0","white","#B40426"))(100),
           display_numbers = display_mat,
           fontsize_number = 12,
           fontsize = 12,
           border_color = NA,
           cellwidth = 70,
           cellheight = 40)
  dev.off()
  
}


####绘制热图
{
  ## 拟时差异基因热图
  #这里是把排序基因（ordergene）提取出来做回归分析，来找它们是否跟拟时间有显著的关系
  #如果不设置，就会用所有基因来做它们与拟时间的相关性
  
  Time_diff <- differentialGeneTest(cds[ordergene,], cores = 1, 
                                    fullModelFormulaStr = "~sm.ns(Pseudotime)")
  Time_diff <- Time_diff[,c(5,2,3,4,1,6)] #把gene放前面，也可以不改
  #write.csv(Time_diff, "2.0Time_diff_all.csv", row.names = F)
  Time_genes <- Time_diff %>% pull(gene_short_name) %>% as.character()
  p=plot_pseudotime_heatmap(cds[Time_genes[1:100],], num_clusters=4, show_rownames=T, return_heatmap=T)
  ggsave("2.1Time_heatmapAll.pdf", p, width = 5, height = 10)
  
  # 提取每个cluster
  p$tree_row
  clusters <- cutree(p$tree_row, k = 4)
  clustering <- data.frame(clusters)
  clustering[,1] <- as.character(clustering[,1])
  colnames(clustering) <- "Gene_Clusters"
  table(clustering)
  write.csv(clustering, "2.2Time_clustering_all.csv", row.names = F)
  
  #前100
  Time_diff = Time_diff[order(Time_diff$qval),]
  #Time_genes <- top_n(Time_diff, n = 100, desc(qval)) %>% pull(gene_short_name) %>% as.character()
  Time_genes = Time_diff$gene_short_name[1:100]
  p = plot_pseudotime_heatmap(cds[Time_genes,], num_clusters=3, show_rownames=T, return_heatmap=T)
  ggsave("2.3Time_heatmapTop100.pdf", p, width = 5, height = 10)

  
  ## 分支BEAM分析
  # 3. 重新运行你的 BEAM 命令
  cds <- qread("/data3/home/hu/Tongfeng/new_Result/Fig4/New result/cds.qs")
  saveRDS(ordergene,"/data3/home/hu/Tongfeng/new_Result/Fig4/New result/ordergene.rds")
  BEAM_res <- BEAM(cds[ordergene,], branch_point = 1, cores = 10, progenitor_method = 'duplicate')
  BEAM_res <- BEAM_res[order(BEAM_res$qval),]
  BEAM_res <- BEAM_res[,c("gene_short_name", "pval", "qval")]
  
  
  plot_genes_branched_heatmap(cds[row.names(subset(BEAM_res, qval < 1e-20)),],
                              branch_point = 1,
                              num_clusters = 3,
                              cores = 1,
                              use_gene_short_name = T,
                              show_rownames = T)
  
  write.csv(BEAM_res, '2.5BEAM_res.csv', quote = F, row.names = F)
  
  #选前100个基因可视化
  #BEAM_genes <- top_n(BEAM_res, n = 100, desc(qval)) %>% pull(gene_short_name) %>% as.character()
  BEAM_res = BEAM_res[order(BEAM_res$qval),]
  BEAM_genes = BEAM_res$gene_short_name[1:100]
  
  
  {
    library(monocle)
    
    # 1. 动态修复构建分支CDS对象时的矩阵类判断漏洞（解决你眼前的报错）
    trace("newCellDataSet", where = asNamespace("monocle"), tracer = quote({
      class <- function(x) base::class(x)[1]
    }), at = 1, print = FALSE)
    
    # 2. 顺便修复后续 orderCells() 可能会触发的同类漏洞（防止下步报错）
    trace("project2MST", where = asNamespace("monocle"), tracer = quote({
      class <- function(x) base::class(x)[1]
    }), at = 1, print = FALSE)
    
    # 1. 提取绘图核心数据集
    heatmap_output <- plot_genes_branched_heatmap(
      cds[row.names(subset(BEAM_res, qval < 1e-20)), ],
      branch_point = 1,
      num_clusters = 3,
      cores = 1,
      use_gene_short_name = TRUE,
      show_rownames = TRUE,
      return_heatmap = TRUE  # ◄ 核心改动：只返回数据，不直接画图
    )
    
    # 2. 从返回的对象中解出矩阵和上下游注释框
    heatmap_matrix <- heatmap_output$heatmap_matrix
    annotation_row <- heatmap_output$annotation_row
    annotation_col <- heatmap_output$annotation_col
    
    # 3. 使用 pheatmap 官方标准库进行任意细节的高度美化
    library(pheatmap)
    pheatmap(
      heatmap_matrix,
      cluster_rows = TRUE,          # 对基因进行聚类
      cluster_cols = FALSE,         # 拟时序的列保持时间轴，不聚类
      annotation_col = annotation_col, # 顶部自动添加分支Fate1/Fate2标签
      annotation_row = annotation_row, # 左侧自动添加基因Cluster分组标签
      show_rownames = FALSE,         # 显示基因名
      show_colnames = FALSE,        # 隐藏密集的细胞名
      fontsize_row = 6,             # 调整基因名称字体大小避免重叠
      color = colorRampPalette(c("#313695", "white", "#a50026"))(100) # 改用高级学术红蓝配色
    )
    
  }
  
  
  # 提取每个cluster
  p$annotation_row
  clustering <- data.frame(annotation_row)
  clustering[,1] <- as.character(clustering[,1])
  colnames(clustering) <- "Gene_Clusters"
  table(clustering)
  write.csv(clustering, "2.6BEAM_clustering_top100.csv", row.names = T)
  
  
  
  ## 富集
  library(clusterProfiler)#用来做富集分析
  library(org.Hs.eg.db)
  library(ggplot2)
  ###### DEG转换为ENTREZID
  for (i in 1:3) {
    DEG.gene_symbol = as.character(rownames(clustering)[clustering$Gene_Clusters==i])
    DEG.entrez_id = mapIds(x = org.Hs.eg.db,
                           keys = DEG.gene_symbol,
                           keytype = "SYMBOL",
                           column = "ENTREZID")
    #有NA则去掉
    DEG.entrez_id = na.omit(DEG.entrez_id)
    
    ####### GO
    erich.go.all = enrichGO(gene = DEG.entrez_id,
                            OrgDb = org.Hs.eg.db,
                            keyType = "ENTREZID",
                            ont = "all",
                            #pAdjustMethod = "BH",
                            pvalueCutoff = 0.05,
                            qvalueCutoff = 0.05,
                            readable = TRUE)
    View(erich.go.all@result)
    
    path <- paste0("2.7.GO_BEAM_c",i,".csv")
    write.csv(erich.go.all@result,path)
  }
  
}


}
