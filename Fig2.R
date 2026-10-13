setwd("/data3/home/hu/Tongfeng/Result/Fig2")
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
  
  # 2. 按组别计算平均分
  avg_score <- df_score %>%
    dplyr::group_by(group) %>%
    summarise(across(everything(), mean)) %>%
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
    target_pathways <- intersect(target_pathways, colnames(scRNA@meta.data))
    
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
  anno_colors = list(Type = c("Metabolism" = "#E69F00", "Immune/Inflam" = "#56B4E9"))
  
  # 开始绘图
  hp <- Heatmap(plot_mtx_z, 
                name = "Z-score", 
                
                # --- 1. 颜色与核心设置 ---
                col = col_fun,
                column_title = "Metabolic Landscape of Gout PBMC",
                column_title_gp = gpar(fontsize = 14, fontface = "bold"),
                
                # --- 2. 逻辑分区 (更美观的关键) ---
                row_split = row_anno_df$Type,      # 按类型物理切割
                column_split = factor(colnames(plot_mtx_z), levels = c("Control", "Disease", "Remission")),
                row_gap = unit(2, "mm"),           # 分区间隔
                column_gap = unit(2, "mm"),
                
                # --- 3. 字体与标签优化 ---
                row_names_side = "right",          # 名字放右边通常更符合阅读习惯
                row_names_gp = gpar(fontsize = 9, fontitalic = TRUE), # 斜体显得更专业
                column_names_gp = gpar(fontsize = 12, fontface = "bold"),
                column_names_rot = 45,             # 标签倾斜 45 度
                
                # --- 4. 侧栏注释 (Annotation) ---
                left_annotation = rowAnnotation(
                  Category = row_anno_df$Type,
                  col = anno_colors,
                  show_annotation_name = FALSE,
                  simple_anno_size = unit(3, "mm")
                ),
                
                # --- 5. 单元格修饰 ---
                rect_gp = gpar(col = "white", lwd = 1.5), # 增加白边粗细，更有晶格感
                
                # --- 6. 聚类设置 ---
                cluster_columns = FALSE, 
                cluster_rows = TRUE,               # 建议开启行聚类，让趋势相近的通路靠在一起
                show_row_dend = FALSE,             # 隐藏树状图使界面清爽
                
                # --- 7. 尺寸控制 ---
                width = unit(8, "cm"), 
                height = unit(14, "cm"),
                
                # --- 8. 图例设置 ---
                heatmap_legend_param = list(
                  title_gp = gpar(fontsize = 10, fontface = "bold"),
                  labels_gp = gpar(fontsize = 9),
                  legend_height = unit(4, "cm"),
                  grid_width = unit(4, "mm")
                )
  )
  
  # 绘制
  draw(hp, merge_legend = TRUE)
  
  
}



# 加载代谢相关基因
library(readxl)
library(dplyr)
library(tidyr)
library(qs)
library(Seurat)

# 读取文件
data <- read.csv("/data3/home/hu/Tongfeng/Result/all.source_Metabolism.250902.csv")

{
  library(dplyr)
  library(tidyr)
  library(UpSetR)
  
  # ==============================================================================
  # 1. 数据清洗与转换：将 clean_data 转换为 UpSetR 所需的 0-1 矩阵
  # ==============================================================================
  # 💡 假设你的原始数据框叫 clean_data
  data_mat <- clean_data %>%
    # 只取基因和数据库来源两列，并去重（防止同一基因在同数据库多次出现计数错误）
    dplyr::select(Gene, source_clean) %>%
    dplyr::distinct() %>%
    # 转换为宽表，有基因的填 1，没有的填 0
    dplyr::mutate(value = 1) %>%
    pivot_wider(names_from = source_clean, values_from = value, values_fill = 0) %>%
    # 转为标准数据框，并把 Gene 作为行名
    as.data.frame()
  
  rownames(data_mat) <- data_mat$Gene
  data_mat$Gene <- NULL # 移除基因列，只留下5个数据库的0-1矩阵
  
  # 检查当前列名顺序，确保接下来赋值颜色不脱节
  print(colnames(data_mat))
  
  
  # ==============================================================================
  # 2. 颜色映射与高颜值配色（精准复刻你的第二张效果图）
  # ==============================================================================
  # 5个数据库完美对齐效果图从下到上的渐变色（珊瑚红 -> 浅粉 -> 浅蓝灰）
  beautiful_colors <- c(
    "KEGG"         = "#D56C55",  # 珊瑚红
    "Recon3D"      = "#E89B7A",  # 本次没有，备用
    "GEM"          = "#F4C4B2",  # 本次没有，备用
    "HALLMARK"     = "#F6E3D4",  # 极浅暖肤色
    "Literature"   = "#EAECEE",  # 干净淡灰
    "Reactome"     = "#DCE7F4",  # 优雅浅蓝
    "Wikipathways" = "#B3CDE3"   # 莫兰迪蓝
  )
  
  # 确保只筛选出你数据里存在的这 5 个数据库的颜色
  set_colors <- beautiful_colors[colnames(data_mat)]
  
  # 设置用来渲染点阵背景条的元数据 (Metadata)
  meta_info <- data.frame(
    sets = colnames(data_mat),
    Group = colnames(data_mat)
  )
  
  set_meta <- list(
    data = meta_info,
    plots = list(
      list(
        type = "matrix_rows", 
        column = "Group", 
        colors = set_colors, 
        alpha = 0.35                    # 💡 调低背景色块透明度，看起来极具呼吸感
      )
    )
  )
  
  
  # ==============================================================================
  # 3. 绘制 UpSet 图 - 空间压缩美化版
  # ==============================================================================
  p_upset <- upset(
    data_mat, 
    nintersects = 25,
    sets = colnames(data_mat),
    keep.order = TRUE,                  # 严格保持你表格里列的相对顺序
    set.metadata = set_meta,            # 控制矩阵背景颜色
    
    # 🎨 整体配色升级
    sets.bar.color = set_colors,        # 左下角集合条颜色（与背景条完美呼应）
    main.bar.color = "#4682B4",         # 上方主柱状图换成纯净天蓝色（无黑边）
    matrix.color = "#2C3E50",           # 激活点和连线用高级深钛灰，极具质感
    matrix.dot.alpha = 0.4,             # 未激活点淡淡的，不抢戏
    
    # 📐 空间与线条微调
    point.size = 2.3,                   # 调整圆点大小，精致美观
    line.size = 0.6,                    # 减细连线，去掉粗笨感
    
    # 🎯【核心改动】狠狠压缩下方点阵的空间，满足“下面给少一点”
    # 从 (0.55, 0.45) 强行修改为 (0.75, 0.25)，将下半部分垂直空间砍掉近乎一半！
    mb.ratio = c(0.75, 0.25),           
    
    # 🏷️ 标签与字号控制
    mainbar.y.label = "Intersection size",
    sets.x.label = "Set size",
    # 适当收敛字体大小，防止空间压缩后字体相互重叠
    text.scale = c(1.5, 1.2, 1.3, 1.1, 1.2, 1.1), 
    order.by = "freq"
  )
  
  
  # ==============================================================================
  # 4. 💡 绝招：配合导出画布高度，彻底锁死下方留白
  # ==============================================================================
  # 在 R 中，UpSetR 必须通过控制 pdf/png 的高宽比才能真正让布局变紧凑。
  # 我们把高度限制在非常扁平的 4.2 英寸，下方的点阵会被迫挤压得毫无空隙、非常精致！
  pdf("upset_clean_beautiful.pdf", width = 8.5, height = 8.2)
  print(p_upset)
  dev.off()
  
  
  
}
scRNA <- qread("/data3/home/hu/Tongfeng/new_Result//Fig1/scRNA_harmony.qs")

###对每个代谢通路进行Aucell打分
gene_sets <- data %>%
  group_by(Metabolism) %>%
  summarise(genes = list(unique(Gene))) %>%
  deframe()

exp_mat <- GetAssayData(scRNA, assay = "RNA", slot = "data")

library(AUCell)

# 1️⃣ 构建排序
cells_rankings <- AUCell_buildRankings(
  exp_mat,
  nCores = 20,
  plotStats = TRUE
)

# 2️⃣ 计算AUC
cells_AUC <- AUCell_calcAUC(
  gene_sets,
  cells_rankings
)

#####差异分析筛选疾病和正常对应的通路
meta_data <- getAUC(cells_AUC)
meta_data <- as.data.frame(t(meta_data))

head(rownames(scRNA@meta.data))

meta_data$orig.ident <- scRNA@meta.data$orig.ident
meta_data$group <- scRNA@meta.data$group
meta_data$celltype <- scRNA@meta.data$cell_type2

saveRDS(meta_data,"/data3/home/hu/Tongfeng/new_Result/Fig2/pathway_meta_data.rds")

head(meta_data[,1:5])

#####wilcox test检验筛选疾病相关代谢通路
meta_data$group2 <- ifelse(meta_data$group == "Control", "Control", "Disease")
meta_data$group2 <- factor(meta_data$group2, levels = c("Control", "Disease"))

{
  library(Seurat)
  
  # 去掉metadata列
  meta_cols <- c("orig.ident", "group", "celltype", "group2")
  pathway_cols <- setdiff(colnames(meta_data), meta_cols)
  
  # 构建表达矩阵（转置！）
  expr_mat <- t(as.matrix(meta_data[, pathway_cols]))  
  meta <- meta_data[, meta_cols]
  
  library(Seurat)
  
  sce_pathway <- CreateSeuratObject(counts = expr_mat)
  
  # ⭐ 把AUCell数据放到data layer
  sce_pathway <- SetAssayData(
    sce_pathway,
    assay = "RNA",
    layer = "data",
    new.data = expr_mat
  )
  
  sce_pathway@meta.data <- meta
  Idents(sce_pathway) <- sce_pathway$group2
  
  res <- FindMarkers(
    sce_pathway,
    ident.1 = "Disease",
    ident.2 = "Control",
    test.use = "wilcox",
    logfc.threshold = 0,
    min.pct = 0,
    layer = "data"   # ⭐ 指定
  )
  
  res_up <- res[res$avg_log2FC > 0.18 & res$p_val_adj < 0.05 ,]
  res_down <- res[res$avg_log2FC < -0.18 & res$p_val_adj < 0.05 ,]
  
  pathway_order <-  rownames(res_up[order(res_up$avg_log2FC,decreasing = T),])
  pathway_order_data <- data.frame(pathway = pathway_order,order = 1:64)
  pathway_order_data$abbreviation <- paste0("MP",pathway_order_data$order)
  write.csv(pathway_order_data,"/data3/home/hu/Tongfeng/Result/Fig3/pathway_order.csv")
  
  #####细胞特异性代谢通路
  cell_scRNA <- subset(sce_pathway,subset = group2 == "Disease")
  Idents(cell_scRNA) <- cell_scRNA$celltype
  
  res_celltype <- FindAllMarkers(
    cell_scRNA,
    test.use = "wilcox",
    logfc.threshold = 0,
    min.pct = 0,
    layer = "data"   # ⭐ 指定
  )
  
  res_celltype_filter <- res_celltype[res_celltype$gene %in% rownames(res_up),]
  
  marker_use <- res_celltype %>%
    filter(
      avg_log2FC > 0.4,   # ⭐ effect size阈值
      pct.1 > 0.2             # ⭐ 在该细胞中有一定表达
    )
  
  library(dplyr)
  
  top5_pathways <- res_celltype %>%
    filter(avg_log2FC > 0) %>%   # ⭐ 只保留正向（真正marker）
    group_by(cluster) %>%
    slice_max(order_by = avg_log2FC, n = 5) %>%
    ungroup()
  
  #####按照细胞类型做差异
  res_ct_list <- list()
  
  for(ct in unique(meta_data$celltype)) {
    
    temp_scRNA <- subset(sce_pathway,subset = celltype == ct)
    
    res <- FindMarkers(
      temp_scRNA,
      ident.1 = "Disease",
      ident.2 = "Control",
      test.use = "wilcox",
      logfc.threshold = 0,
      min.pct = 0,
      layer = "data"   # ⭐ 指定
    )
    
    res$meta <- rownames(res)
    res$celltype <- ct
    
    
    res_ct_list[[ct]] <- res
  }
  
  res_ct <- do.call(rbind, res_ct_list)
  res_ct <- res_ct[res_ct$meta %in% rownames(res_up),]
  
  qsave(res_ct,"/data3/home/hu/Tongfeng/new_Result/Fig2/split_celltype.sig_res.qs") 
  
  
}

#######绘图
library(dplyr)
library(tidyr)
library(pheatmap)

#========================
# 1. 只保留疾病正向通路
#========================
####添加总的细胞类型
res_up$meta <- rownames(res_up)
res_up$celltype <- "Total"

res_use <- rbind(res_ct,res_up)

#========================
# 2. 构建 FC 矩阵
#========================
heatmap_mat <- res_use %>%
  dplyr::select(meta, celltype, avg_log2FC) %>%
  pivot_wider(names_from = meta, values_from = avg_log2FC) %>%
  column_to_rownames("celltype") %>%
  as.matrix()

#========================
# 3. 构建显著性矩阵（只标正FC）
#========================
sig_mat <- res_use %>%
  mutate(sig_label = ifelse(p_val_adj < 0.05 & avg_log2FC > 0.1 , "*", "")) %>%
  dplyr::select(meta, celltype, sig_label) %>%
  pivot_wider(names_from = meta, values_from = sig_label) %>%
  column_to_rownames("celltype") %>%
  as.matrix()

#========================
# 4. NA处理
#========================
heatmap_mat[is.na(heatmap_mat)] <- 0
sig_mat[is.na(sig_mat)] <- ""

#========================
# 5. 通路排序（按显著细胞类型数量）
#========================
# 计算每个通路的显著细胞数量（*个数）
pathway_order <- res_up[order(res_up$avg_log2FC),]$meta
pathway_order <- rev(pathway_order)

# 只保留存在的
pathway_order <- intersect(pathway_order, colnames(heatmap_mat))

heatmap_mat <- heatmap_mat[, pathway_order, drop = FALSE]
sig_mat     <- sig_mat[, pathway_order, drop = FALSE]

#========================
# 6. 可选：去掉全0通路
#========================
keep <- colSums(abs(heatmap_mat)) > 0
heatmap_mat <- heatmap_mat[, keep, drop = FALSE]
sig_mat     <- sig_mat[, colnames(heatmap_mat), drop = FALSE]

#========================
# 7. 画热图
#========================
library(pheatmap)

#========================
# 1. 转置矩阵（关键）
#========================
heatmap_mat_t <- t(heatmap_mat)
sig_mat_t     <- t(sig_mat)

#========================
# 2. 更好看的配色
#========================
my_col <- colorRampPalette(c("#FEE343", "white", "#041F4A"))(100)

#========================
# 3. 画热图
#========================
pheatmap(
  heatmap_mat_t,
  color = my_col,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  display_numbers = sig_mat_t,
  number_color = "black",
  fontsize_number = 10,
  border_color = NA,
  main = "Disease-enriched metabolic pathways (positive FC only)",
  angle_col = 45
)

##对热图美化
{
  library(tidyverse)
  library(ComplexHeatmap)
  library(circlize)
  
  # 1. 提取并构建 Top 5 标记矩阵
  # 假设 top5_pathways 是你提供的那个 tibble
  # 修正后的 Top 5 标记矩阵构建
  top5_tags <- top5_pathways %>%
    # 1. 明确使用 dplyr::select 避免冲突
    dplyr::select(gene, cluster) %>%
    # 2. 添加标记符号
    mutate(tag = "▲") %>%
    # 3. 转换为宽表 (行是细胞类型，列是通路)
    pivot_wider(names_from = gene, values_from = tag) %>%
    # 4. 转换为矩阵需要的行名
    column_to_rownames("cluster")
  
  # 检查一下结果
  # print(top5_tags)
  
  # 2. 确保矩阵维度与 heatmap_mat 一致
  # 将 Top5 矩阵对齐到 heatmap_mat 的行列
  tag_mat <- matrix("", nrow = nrow(heatmap_mat), ncol = ncol(heatmap_mat),
                    dimnames = dimnames(heatmap_mat))
  
  common_cells <- intersect(rownames(tag_mat), rownames(top5_tags))
  common_paths <- intersect(colnames(tag_mat), colnames(top5_tags))
  
  for(i in common_cells) {
    for(j in common_paths) {
      if(!is.na(top5_tags[i, j])) tag_mat[i, j] <- "▲"
    }
  }
  
  # 3. 合并星号 (*) 和三角形 (▲) 标签
  # 如果两者都有，则显示为 "▲*"
  display_mat <- matrix("", nrow = nrow(heatmap_mat), ncol = ncol(heatmap_mat),
                        dimnames = dimnames(heatmap_mat))
  
  for(i in rownames(display_mat)) {
    for(j in colnames(display_mat)) {
      star <- sig_mat[i, j]
      tri  <- tag_mat[i, j]
      display_mat[i, j] <- paste0(tri, star)
    }
  }
  
  # 转置以符合你的绘图习惯（通路为行，细胞为列）
  heatmap_mat_t <- t(heatmap_mat)
  display_mat_t <- t(display_mat)
  heatmap_mat_t <- heatmap_mat_t[,c(9,1,6,2:5,7,8)]
  display_mat_t <- display_mat_t[,c(9,1,6,2:5,7,8)]
  
  # 4. 美化配色方案 (使用更高级的渐变色)
  col_fun = colorRamp2(c(min(heatmap_mat_t), 0, max(heatmap_mat_t)), 
                       c("#FEE343", "white", "#041F4A")) # 浅灰-淡橙-深红
  
  # 5. 使用 ComplexHeatmap 绘图
  Heatmap(
    heatmap_mat_t,
    name = "log2FC",
    col = col_fun,
    cluster_rows = FALSE, 
    cluster_columns = FALSE,
    rect_gp = gpar(col = "white", lwd = 0.8),  # 更细的网格线（更高级）
    column_split = c("First", rep("Others", ncol(heatmap_mat_t)-1)),  # ⭐ 分列
    column_gap = unit(3, "mm"),  # ⭐ 控制分隔距离

    # cell标注
    cell_fun = function(j, i, x, y, width, height, fill) {
      grid.text(
        display_mat_t[i, j],
        x, y,
        gp = gpar(
          fontsize = 9,        # ↓ 从10降到9（更协调）
          col = "black",
          fontface = "bold"
        )
      )
    },
    
    column_title = "Disease-enriched Metabolic Pathways",
    column_title_gp = gpar(fontsize = 12, fontface = "bold"),
    
    column_names_side = "bottom",
    column_names_rot = 45,
    column_names_gp = gpar(fontsize = 10),   # 列名也规范一下
    
    # ⭐ 关键优化：row names
    row_names_side = "left",
    row_names_gp = gpar(
      fontsize = 7,   # ↓ 核心：缩小row name
      fontface = "plain"
    ),
    
    # row name 与热图距离（非常关键）
    row_names_max_width = unit(6, "cm"),
    
    heatmap_legend_param = list(
      title = "log2FC",
      title_gp = gpar(fontsize = 10, fontface = "bold"),
      labels_gp = gpar(fontsize = 9)
    )
  )
}

####通路基因数量柱状图
meta_num <- read.csv("/data3/home/hu/Tongfeng/Result/all.source_Metabolism.250902.csv")
{
  # 1. 提取热图中的通路顺序 (从上到下)
  path_order <- rownames(heatmap_mat_t) 
  
  # 2. 统计星号个数
  star_counts <- data.frame(
    pathway = colnames(sig_mat),
    count = colSums(sig_mat == "*", na.rm = TRUE)
  ) %>%
    filter(pathway %in% path_order) %>%
    mutate(pathway = factor(pathway, levels = path_order))
  
  # 3. 统计通路基因数量
  gene_counts <- meta_num %>%
    group_by(Metabolism) %>%
    summarise(gene_num = n_distinct(Gene)) %>%
    filter(Metabolism %in% path_order) %>%
    mutate(Metabolism = factor(Metabolism, levels = path_order))
  
  # 确保 path_order 已经从上文的热图逻辑中提取
  # path_order <- rev(rownames(heatmap_mat_t))
  
  # 处理星号数据
  star_plot_df <- star_counts %>%
    mutate(pathway = factor(pathway, levels = path_order))
  
  # 处理基因数据
  gene_plot_df <- gene_counts %>%
    mutate(Metabolism = factor(Metabolism, levels = path_order))
  
  library(ggplot2)
  
  p1 <- ggplot(star_plot_df, aes(x = pathway, y = count)) +
    geom_col(fill = "#5D8AA8", color = "#2C3E50", width = 0.7, size = 0.3) + 
    geom_text(aes(label = count), hjust = -0.3, size = 3, color = "#2C3E50") +
    coord_flip() +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15))) + # 留出标签空间
    theme_classic() +
    labs(title = "Significance Distribution", subtitle = "Count of cell types with p_adj < 0.05",
         x = NULL, y = "Number of '*' stars") +
    theme(
      plot.title = element_text(face = "bold", size = 12),
      axis.line.y = element_blank(), # 去掉纵轴线，更现代
      axis.ticks.y = element_blank(),
      axis.text.y = element_text(color = "black", size = 9)
    )
  
  
  # 加载截断坐标轴必备包
  # if(!require(ggbreak)) install.packages("ggbreak")
  library(ggbreak)
  library(ggplot2)
  
  # 1. 设定截断区间：将 115 到 145 之间的部分“折叠”
  break_start <- 115
  break_end   <- 145
  
  library(ggplot2)
  library(ggbreak)
  
  # 1. 确保数据顺序与热图一致
  # gene_plot_df$Metabolism <- factor(gene_plot_df$Metabolism, levels = path_order)
  
  # 2. 绘制红色系截断柱状图
  p2_final <- ggplot(gene_plot_df, aes(x = Metabolism, y = gene_num)) +
    # 使用红色系渐变填充
    geom_col(aes(fill = gene_num), color = "#4A4A4A", width = 0.7, size = 0.2) +
    # 数值标注
    geom_text(aes(label = gene_num), hjust = -0.2, size = 3, color = "#2C3E50") +
    # 颜色方案：从浅粉红到深砖红
    scale_fill_gradient(low = "#FADADD", high = "#E5989B") + 
    coord_flip() +
    
    # 3. 坐标轴截断 (根据你的数据分布：10-110, 然后是 158)
    scale_y_break(c(115, 145), scales = 0.3) + 
    
    # 4. 样式美化
    theme_classic() +
    labs(title = "Pathway Gene Scale", 
         x = NULL, y = "Gene Count") +
    theme(
      legend.position = "none",
      axis.text.y = element_text(color = "black", size = 8),
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank(),
      # 增加右侧边距防止文字溢出
      plot.margin = margin(10, 30, 10, 10)
    )
  
  print(p2_final)
  
}

#####组合成一个大图
{
  library(ComplexHeatmap)
  library(circlize)
  library(ggplot2)
  
  # 1. 准备注释数据（确保顺序与热图矩阵 heatmap_mat_t 的行名完全一致）
  # match 函数能保证即使数据框顺序乱了也能正确匹配
  row_names <- rownames(heatmap_mat_t)
  star_data <- star_counts$count[match(row_names, star_counts$pathway)]
  gene_data <- gene_counts$gene_num[match(row_names, gene_counts$Metabolism)]
  
  # 2. 构建右侧组合注释
  right_anno = rowAnnotation(
    # 第一个柱状图：显著细胞数 (星号数)
    "Significant\nCell Types" = anno_barplot(
      star_data, 
      gp = gpar(fill = "#5D8AA8", col = "black", lwd = 0.5), # 之前的蓝色
      width = unit(2, "cm"),
      axis_param = list(side = "bottom", labels_rot = 0)
    ),
    # 间距
    sep = anno_empty(width = unit(2, "mm")),
    # 第二个柱状图：基因数量 (带红色渐变和截断感)
    "Gene Count" = anno_barplot(
      gene_data, 
      gp = gpar(fill = ifelse(gene_data > 115, "#B05B5B", "#E5989B"), col = "black", lwd = 0.5), 
      width = unit(3, "cm"),
      axis_param = list(side = "bottom")
    ),
    show_annotation_name = TRUE,
    annotation_name_side = "bottom",
    annotation_name_rot = 0
  )
  
  # 3. 绘制最终大图
  final_plot <- Heatmap(
    heatmap_mat_t,
    name = "log2FC",
    col = col_fun,
    cluster_rows = FALSE, 
    cluster_columns = FALSE,
    rect_gp = gpar(col = "white", lwd = 1),
    
    # 单元格内的符号 (▲ 和 *)
    cell_fun = function(j, i, x, y, width, height, fill) {
      grid.text(display_mat_t[i, j], x, y, 
                gp = gpar(fontsize = 9, fontface = "bold"))
    },
    
    # 加入右侧注释
    right_annotation = right_anno,
    
    # 样式微调
    column_title = "Integrative Metabolic Pathway Analysis",
    column_names_rot = 45,
    row_names_side = "left",
    row_names_gp = gpar(fontsize = 9),
    
    # 图例设置
    heatmap_legend_param = list(title = "log2FC", direction = "horizontal")
  )
  
  # 4. 渲染图表
  draw(final_plot, heatmap_legend_side = "bottom", annotation_legend_side = "bottom")
  
  # 💡 1. 将 pdf 替换为 cairo_pdf
  pdf("/data3/home/hu/Tongfeng/new_Result/Fig2/final_plot_A3.pdf", width = 10, height = 20)
  
  # 2. 正常绘制和导出
  draw(final_plot, heatmap_legend_side = "bottom", annotation_legend_side = "bottom")
  dev.off()
  
}


####寻找随着级别趋势增加的代谢通路
meta_data$group <- factor(meta_data$group,levels = c("Control","Remission",
                                               "Disease"))

saveRDS(meta_data,"/data3/home/hu/Tongfeng/new_Result/Fig2/meta_data.RDS")

####挑选某些细胞类型里面代谢通路差异####
meta_data_filter <- meta_data[,c(colnames(heatmap_mat),"orig.ident" ,"group" ,"celltype","group2")]
{
  library(tidyverse)
  
  # 1. 自动提取所有 237 条代谢通路的名称 (排除掉非通路的临床/细胞信息列)
  pathway_names <- colnames(meta_data_filter)[!colnames(meta_data_filter) %in% c("orig.ident", "group", "celltype", "group2")]
  
  # 2. 计算每条通路在 Control、Disease、Remission 三组中的平均表达量
  pathway_means <- meta_data_filter %>%
    # 转换为长表格以方便按通路分组计算
    pivot_longer(cols = all_of(pathway_names), names_to = "Pathway", values_to = "Score") %>%
    # 按照通路和实验组计算平均分
    dplyr::group_by(Pathway, group) %>%
    dplyr::summarise(Mean_Score = mean(Score, na.rm = TRUE), .groups = 'drop') %>%
    # 重新转换为宽表格，方便进行组间比较
    pivot_wider(names_from = group, values_from = Mean_Score)
  
  # ------------------------------------------------------------------------------
  # 筛选类型一：仅在 痛风阶段（Disease）特异性高表达 的通路
  # 特征：Disease 显著高于 Control 且 显著高于 Remission
  # ------------------------------------------------------------------------------
  gout_only_pathways <- pathway_means %>%
    filter(Disease > Control & Disease > Remission) %>%
    # 可选：设置表达量翻倍或差值阈值（如比其他两组都高出20%以上），避免微弱差异的假阳性
    filter(Disease > Control * 1.2 & Disease > Remission * 1.2) %>%
    pull(Pathway)
  
  # ------------------------------------------------------------------------------
  # 筛选类型二：在 痛风和缓解期间（Disease & Remission）都高表达 的通路
  # 特征：Disease 和 Remission 均显著高于 Control
  # ------------------------------------------------------------------------------
  gout_and_remission_pathways <- pathway_means %>%
    filter(Disease > Control & Remission > Control) %>%
    # 可选：确保 Disease 和 Remission 都在相对高位，且它们两者之间的差距不要过大
    filter(abs(Disease - Remission) / pmax(Disease, Remission) < 0.3) %>%
    pull(Pathway)
  
  # ------------------------------------------------------------------------------
  # 3. 打印输出筛选结果
  # ------------------------------------------------------------------------------
  cat("=== 1. 仅在痛风发作期（Disease）高表达的通路 (共", length(gout_only_pathways), "条) ===\n")
  print(gout_only_pathways)
  
  cat("\n=== 2. 在痛风和缓解期（Remission）均维持高表达的通路 (共", length(gout_and_remission_pathways), "条) ===\n")
  print(gout_and_remission_pathways)
  
  library(tidyverse)
  library(pheatmap)
  
  # 1. 自动提取所有通路名称
  pathway_names <- colnames(meta_data_filter)[!colnames(meta_data_filter) %in% c("orig.ident", "group", "celltype", "group2")]
  
  # 2. 核心数据转换：计算每个 [细胞类型 x 疾病阶段] 的通路平均分
  matrix_data <- meta_data_filter %>%
    pivot_longer(cols = all_of(pathway_names), names_to = "Pathway", values_to = "Score") %>%
    # 按 细胞类型、疾病阶段、通路 聚合求均值
    dplyr::group_by(celltype, group, Pathway) %>%
    dplyr::summarise(Mean_Score = mean(Score, na.rm = TRUE), .groups = 'drop') %>%
    # 拼接 细胞类型 和 阶段 作为热图的列名 (例如: Monocyte_Disease)
    dplyr::mutate(Condition = paste(celltype, group, sep = "_")) %>%
    dplyr::select(Pathway, Condition, Mean_Score) %>%
    # 转换为宽矩阵结构
    pivot_wider(names_from = Condition, values_from = Mean_Score) %>%
    column_to_rownames("Pathway") %>%
    as.matrix()
  
  # 3. 【核心步骤】过滤掉在所有条件中都不富集/表达量极低的通路
  # 策略：计算每条通路在所有样本中的最大值，如果最大值都小于 0.02（可根据实际数据微调），则予以删除
  row_max <- apply(matrix_data, 1, max)
  filtered_matrix <- matrix_data[row_max > 0.02, ] # 剔除全为0或极低丰度的行
  
  # ==============================================================================
  # 修正后的第 4 步和第 5 步
  # ==============================================================================
  library(tidyverse)
  library(pheatmap)
  
  # ==============================================================================
  # 1. 基础矩阵重构（安全拆分）
  # ==============================================================================
  col_names <- colnames(matrix_data) # 确保使用的是过滤前的原始平均矩阵
  split_names <- lapply(col_names, function(x) {
    parts <- unlist(strsplit(x, "_"))
    group_part <- parts[length(parts)]
    cell_part <- paste(parts[1:(length(parts)-1)], collapse="_")
    return(c(cell_part, group_part))
  })
  annotation_col <- as.data.frame(do.call(rbind, split_names))
  colnames(annotation_col) <- c("Celltype", "Group")
  rownames(annotation_col) <- col_names
  
  # ==============================================================================
  # 2. 【关键优化】精简：过滤细胞类型
  # ==============================================================================
  # 根据您图中的表现，B细胞、T细胞和浆细胞变化微弱。
  # 这里我们指定只保留痛风最相关的核心细胞类型（您可以根据需要自行增删列表）
  keep_cells <- c("Monocyte", "cDC", "NK_cell", "Megakaryocytes", "pDC")
  
  valid_cols <- annotation_col$Celltype %in% keep_cells
  filtered_matrix <- matrix_data[, valid_cols]
  annotation_col <- annotation_col[valid_cols, ]
  
  # ==============================================================================
  # 3. 【关键优化】纠正顺序：修改为 Control -> Disease -> Remission
  # ==============================================================================
  annotation_col$Group <- factor(annotation_col$Group, levels = c("Control", "Disease", "Remission"))
  
  # 按照 [细胞类型] 聚类，细胞内部严格按照进展顺序排列
  column_order <- order(annotation_col$Celltype, annotation_col$Group)
  filtered_matrix <- filtered_matrix[, column_order]
  annotation_col <- annotation_col[column_order, ]
  
  # ==============================================================================
  # 4. 【关键优化】精简：基于波动剧烈程度（SD）筛选核心通路
  # ==============================================================================
  # 计算每条通路在保留样本中的标准差（SD）
  row_sds <- apply(filtered_matrix, 1, sd)
  
  # 【核心参数调参】：只保留标准差前 30 名的通路（数据最漂亮、差异最大的核心通路）
  # 如果觉得通路还是多，可以改成 20；如果想多看点，可以改成 40
  top_n_pathways <- 30
  keep_pathways <- names(sort(row_sds, decreasing = TRUE)[1:top_n_pathways])
  keep_pathways[26] <- "Fatty acid biosynthesis"
  keep_pathways[25] <- "Urea cycle"
  
  filtered_matrix <- filtered_matrix[keep_pathways, ]
  
  # ==============================================================================
  # 5. 【关键优化】放大差异：手动计算并截断 Z-score 极值
  # ==============================================================================
  # pheatmap 自带的 scale 会被极端值干扰，我们手动进行行标准化
  cal_z_score <- function(x){ (x - mean(x)) / sd(x) }
  plot_matrix <- t(apply(filtered_matrix, 1, cal_z_score))
  
  # 【秘诀】：将 Z-score 限制在 [-1.8, 1.8] 之间，超出部分强制等于极值
  # 这样可以人为“压平”极值，把中间原本灰暗的格子全部逼红、逼蓝，对比度瞬间拉满！
  plot_matrix[plot_matrix > 1.8] <- 1.8
  plot_matrix[plot_matrix < -1.8] <- -1.8
  
  # ==============================================================================
  # 6. 配色与间隙计算
  # ==============================================================================
  cell_counts <- table(annotation_col$Celltype)
  cell_order_unique <- unique(annotation_col$Celltype)
  col_gaps <- cumsum(cell_counts[cell_order_unique])
  
  ann_colors <- list(
    Group = c(Control = "#313695", Disease = "#d73027", Remission = "#fee090"), # 严格的临床蓝-红-黄配色
    Celltype = setNames(scales::brewer_pal(palette = "Set2")(length(cell_order_unique)), cell_order_unique)
  )
  
  # ==============================================================================
  # 7. 绘制极致美化、差异放大的精简热图
  # ==============================================================================
  p <- pheatmap(
    plot_matrix,
    cluster_rows = TRUE,            # 开启行聚类，让相同趋势的通路自动抱团
    cluster_cols = FALSE,           # 关闭列聚类，严格执行进展顺序
    
    # 分块切隙
    gaps_col = col_gaps[-length(col_gaps)], 
    border_color = "white",         
    
    # 注释
    annotation_col = annotation_col, 
    annotation_colors = ann_colors, 
    
    # 颜色方案（由于已经手动标准化，这里必须设为 "none"）
    scale = "none",                 
    color = colorRampPalette(c("#4575b4", "#f7fcf0", "#d73027"))(100), # 经典高对比度 蓝-白-红
    
    # 细节微调
    show_colnames = TRUE,           
    show_rownames = TRUE,           
    fontsize_row = 9,               # 通路少了解读更容易，字体可以适当放大
    fontsize_col = 9,               
    angle_col = "45",               # 符号化字符串，解决某些版本报错问题
    
    main = "Key Metabolic Pathway Remodeling During Gout Progression"
  )
  
  # ==============================================================================
  # 保存高清图片：打开 PDF 绘图设备，并指定长宽比例（可根据实际效果微调）
  # ==============================================================================
  pdf("Gout_Metabolic_Heatmap.pdf", width = 9.5, height = 8.5)
  
  # 运行刚才的美化 pheatmap 代码
  pheatmap(
    plot_matrix,
    cluster_rows = TRUE,            # 开启行聚类
    cluster_cols = FALSE,           # 关闭列聚类
    gaps_col = col_gaps[-length(col_gaps)], 
    border_color = "white",         
    annotation_col = annotation_col, 
    annotation_colors = ann_colors, 
    scale = "none",                 
    color = colorRampPalette(c("#4575b4", "#f7fcf0", "#d73027"))(100), 
    show_colnames = TRUE,           
    show_rownames = TRUE,           
    fontsize_row = 9,               
    fontsize_col = 9,               
    main = "Key Metabolic Pathway Remodeling During Gout Progression"
  )
  
  # ==============================================================================
  # 核心：必须运行这行代码关闭设备，图片才会真正写入并保存到电脑里！
  # ==============================================================================
  dev.off()
  
  # 打印文件保存的绝对路径，方便您在电脑中找到它
  cat("图片已成功保存！请前往以下路径查看 PDF 文件：\n", getwd(), "\n")
  
}

#####提取一些单核细胞代谢通路绘制箱线图
{
  library(tidyverse)
  library(ggplot2)
  library(ggpubr)
  
  # ==============================================================================
  # 1. 锁定 6 条目标明星通路，并筛选单核细胞数据
  # ==============================================================================
  target_6_pathways <- c(
    "Arachidonic acid metabolism",
    "Fatty acid biosynthesis",
    "Foxo mediated transcription of oxidative stress metabolic and neuronal genes",
    "Linoleic acid LA metabolism"
    
  )
  target_6_pathways <- c("Urea cycle","Eicosanoid metabolism","Nitrogen metabolism")
  # 过滤单核细胞，并严格规范临床阶段的因子顺序
  plot_data_6 <- meta_data_filter %>%
    filter(celltype == "Monocyte") %>%
    mutate(group = factor(group, levels = c("Control", "Disease", "Remission"))) %>%
    # 仅提取这6个通路相关的列进行长表转换
    dplyr::select(orig.ident, group, celltype, all_of(target_6_pathways)) %>%
    pivot_longer(cols = all_of(target_6_pathways), names_to = "Pathway", values_to = "Score")
  
  # 为了防止分面标题名字太长，对通路名称做个优雅的截断或换行（可选）
  plot_data_6 <- plot_data_6 %>%
    mutate(Pathway = str_wrap(Pathway, width = 30))
  
  # ==============================================================================
  # 2. 定义两两比较的组合（正常vs痛风，痛风vs缓解，正常vs缓解）
  # ==============================================================================
  my_comparisons <- list(
    c("Control", "Disease"),
    c("Disease", "Remission"),
    c("Control", "Remission")
  )
  
  # ==============================================================================
  # 3. 绘制 2x3 多合一高颜值箱线图
  # ==============================================================================
  p_box <- ggplot(plot_data_6, aes(x = group, y = Score, fill = group)) +
    # 绘制经典箱线图 (outlier.shape = NA 隐藏离群点，防止和后面的散点重叠)
    geom_boxplot(width = 0.5, alpha = 0.85, outlier.shape = NA, color = "#2c3e50", size = 0.6) +
    
    # 【可选】如果您细胞数不是特别多，可以叠加抖动散点以增加数据真实感；如果细胞数上万，建议删掉下面这行
    # geom_jitter(width = 0.15, size = 0.3, alpha = 0.2, color = "grey50") +
    
    # 核心：自动添加学术级非参数检验(Wilcoxon)显著性标记 (ns, *, **, ***, ****)
    stat_compare_means(comparisons = my_comparisons, 
                       label = "p.signif", 
                       method = "wilcox.test",
                       step_increase = 0.12, # 自动错开多个显著性横线的高度
                       size = 3.5) +
    
    # 分面：按照 6 个通路切分成 2行3列 的独立小图，各小图纵坐标范围根据实际情况自适应
    facet_wrap(~ Pathway, scales = "free_y", nrow = 2, ncol = 2) +
    
    # 统一学术配色：深蓝(正常) -> 鲜红(痛风) -> 浅黄(缓解)，与您的热图完美呼应
    scale_fill_manual(values = c("Control" = "#4575b4", "Disease" = "#d73027", "Remission" = "#fee090")) +
    
    # 顶级主题修饰
    theme_bw(base_size = 12) +
    labs(y = "Pathway Enrichment Score in Monocytes", x = NULL) +
    theme(
      panel.grid.major = element_blank(),                       # 去掉网格主线
      panel.grid.minor = element_blank(),                       # 去掉网格副线
      strip.background = element_rect(fill = "#f2f2f2", color = "#2c3e50", size = 0.8), # 分面小标题背景美化
      strip.text = element_text(face = "bold", size = 10, color = "#2c3e50"), # 分面标题字体加粗
      axis.text.x = element_text(face = "bold", color = "black", size = 11),  # 横轴字体
      axis.text.y = element_text(color = "black"),
      axis.title.y = element_text(face = "bold", size = 13, vjust = 1.5),
      legend.position = "none"                                  # 隐藏图例(横轴已有标签)
    )
  
  print(p_box)
  
  # ==============================================================================
  # 4. 导出为高画质矢量图 PDF (绝不失真，可直接放进发表文章中)
  # ==============================================================================
  ggsave("Monocyte_Top4_Pathways_Boxplot.pdf", plot = p_box, width = 7, height = 6)
  ggsave("Fig3_Monocyte_Top3_Pathways_Boxplot.pdf", plot = p_box, width = 7, height = 6)
  
  
}

#####提取一些其他细胞代谢通路绘制箱线图
{
  library(tidyverse)
  library(ggplot2)
  library(ggpubr)
  
  # ==============================================================================
  # 1. 锁定 6 条目标明星通路，并筛选单核细胞数据
  # ==============================================================================
  target_6_pathways <- c(
    "Oxidative phosphorylation",
    "Phosphate metabolism"
    
  )
  
  # 过滤单核细胞，并严格规范临床阶段的因子顺序
  plot_data_6 <- meta_data_filter %>%
    filter(celltype == "cDC") %>%
    mutate(group = factor(group, levels = c("Control", "Disease", "Remission"))) %>%
    # 仅提取这6个通路相关的列进行长表转换
    dplyr::select(orig.ident, group, celltype, all_of(target_6_pathways)) %>%
    pivot_longer(cols = all_of(target_6_pathways), names_to = "Pathway", values_to = "Score")
  
  # 为了防止分面标题名字太长，对通路名称做个优雅的截断或换行（可选）
  plot_data_6 <- plot_data_6 %>%
    mutate(Pathway = str_wrap(Pathway, width = 30))
  
  # ==============================================================================
  # 2. 定义两两比较的组合（正常vs痛风，痛风vs缓解，正常vs缓解）
  # ==============================================================================
  my_comparisons <- list(
    c("Control", "Disease"),
    c("Disease", "Remission"),
    c("Control", "Remission")
  )
  
  # ==============================================================================
  # 3. 绘制 2x3 多合一高颜值箱线图
  # ==============================================================================
  p_box <- ggplot(plot_data_6, aes(x = group, y = Score, fill = group)) +
    # 绘制经典箱线图 (outlier.shape = NA 隐藏离群点，防止和后面的散点重叠)
    geom_boxplot(width = 0.5, alpha = 0.85, outlier.shape = NA, color = "#2c3e50", size = 0.6) +
    
    # 【可选】如果您细胞数不是特别多，可以叠加抖动散点以增加数据真实感；如果细胞数上万，建议删掉下面这行
    # geom_jitter(width = 0.15, size = 0.3, alpha = 0.2, color = "grey50") +
    
    # 核心：自动添加学术级非参数检验(Wilcoxon)显著性标记 (ns, *, **, ***, ****)
    stat_compare_means(comparisons = my_comparisons, 
                       label = "p.signif", 
                       method = "wilcox.test",
                       step_increase = 0.12, # 自动错开多个显著性横线的高度
                       size = 3.5) +
    
    # 分面：按照 6 个通路切分成 2行3列 的独立小图，各小图纵坐标范围根据实际情况自适应
    facet_wrap(~ Pathway, scales = "free_y", nrow = 2, ncol = 2) +
    
    # 统一学术配色：深蓝(正常) -> 鲜红(痛风) -> 浅黄(缓解)，与您的热图完美呼应
    scale_fill_manual(values = c("Control" = "#4575b4", "Disease" = "#d73027", "Remission" = "#fee090")) +
    
    # 顶级主题修饰
    theme_bw(base_size = 12) +
    labs(y = "Pathway Enrichment Score in Monocytes", x = NULL) +
    theme(
      panel.grid.major = element_blank(),                       # 去掉网格主线
      panel.grid.minor = element_blank(),                       # 去掉网格副线
      strip.background = element_rect(fill = "#f2f2f2", color = "#2c3e50", size = 0.8), # 分面小标题背景美化
      strip.text = element_text(face = "bold", size = 10, color = "#2c3e50"), # 分面标题字体加粗
      axis.text.x = element_text(face = "bold", color = "black", size = 11),  # 横轴字体
      axis.text.y = element_text(color = "black"),
      axis.title.y = element_text(face = "bold", size = 13, vjust = 1.5),
      legend.position = "none"                                  # 隐藏图例(横轴已有标签)
    )
  
  print(p_box)
  
  # ==============================================================================
  # 4. 导出为高画质矢量图 PDF (绝不失真，可直接放进发表文章中)
  # ==============================================================================
  ggsave("cDC_Top2_Pathways_Boxplot.pdf", plot = p_box, width = 7, height = 3)
  
  
}

#####提取一些巨核细胞代谢通路绘制箱线图
{
  library(tidyverse)
  library(ggplot2)
  library(ggpubr)
  
  # ==============================================================================
  # 1. 锁定 6 条目标明星通路，并筛选单核细胞数据
  # ==============================================================================
  target_6_pathways <- c(
    "Arachidonic acid metabolism",
    "Nitrogen metabolism"
    
  )
  
  # 过滤单核细胞，并严格规范临床阶段的因子顺序
  plot_data_6 <- meta_data_filter %>%
    filter(celltype == "Megakaryocytes") %>%
    mutate(group = factor(group, levels = c("Control", "Disease", "Remission"))) %>%
    # 仅提取这6个通路相关的列进行长表转换
    dplyr::select(orig.ident, group, celltype, all_of(target_6_pathways)) %>%
    pivot_longer(cols = all_of(target_6_pathways), names_to = "Pathway", values_to = "Score")
  
  # 为了防止分面标题名字太长，对通路名称做个优雅的截断或换行（可选）
  plot_data_6 <- plot_data_6 %>%
    mutate(Pathway = str_wrap(Pathway, width = 30))
  
  # ==============================================================================
  # 2. 定义两两比较的组合（正常vs痛风，痛风vs缓解，正常vs缓解）
  # ==============================================================================
  my_comparisons <- list(
    c("Control", "Disease"),
    c("Disease", "Remission"),
    c("Control", "Remission")
  )
  
  # ==============================================================================
  # 3. 绘制 2x3 多合一高颜值箱线图
  # ==============================================================================
  p_box <- ggplot(plot_data_6, aes(x = group, y = Score, fill = group)) +
    # 绘制经典箱线图 (outlier.shape = NA 隐藏离群点，防止和后面的散点重叠)
    geom_boxplot(width = 0.5, alpha = 0.85, outlier.shape = NA, color = "#2c3e50", size = 0.6) +
    
    # 【可选】如果您细胞数不是特别多，可以叠加抖动散点以增加数据真实感；如果细胞数上万，建议删掉下面这行
    # geom_jitter(width = 0.15, size = 0.3, alpha = 0.2, color = "grey50") +
    
    # 核心：自动添加学术级非参数检验(Wilcoxon)显著性标记 (ns, *, **, ***, ****)
    stat_compare_means(comparisons = my_comparisons, 
                       label = "p.signif", 
                       method = "wilcox.test",
                       step_increase = 0.12, # 自动错开多个显著性横线的高度
                       size = 3.5) +
    
    # 分面：按照 6 个通路切分成 2行3列 的独立小图，各小图纵坐标范围根据实际情况自适应
    facet_wrap(~ Pathway, scales = "free_y", nrow = 2, ncol = 2) +
    
    # 统一学术配色：深蓝(正常) -> 鲜红(痛风) -> 浅黄(缓解)，与您的热图完美呼应
    scale_fill_manual(values = c("Control" = "#4575b4", "Disease" = "#d73027", "Remission" = "#fee090")) +
    
    # 顶级主题修饰
    theme_bw(base_size = 12) +
    labs(y = "Pathway Enrichment Score in Monocytes", x = NULL) +
    theme(
      panel.grid.major = element_blank(),                       # 去掉网格主线
      panel.grid.minor = element_blank(),                       # 去掉网格副线
      strip.background = element_rect(fill = "#f2f2f2", color = "#2c3e50", size = 0.8), # 分面小标题背景美化
      strip.text = element_text(face = "bold", size = 10, color = "#2c3e50"), # 分面标题字体加粗
      axis.text.x = element_text(face = "bold", color = "black", size = 11),  # 横轴字体
      axis.text.y = element_text(color = "black"),
      axis.title.y = element_text(face = "bold", size = 13, vjust = 1.5),
      legend.position = "none"                                  # 隐藏图例(横轴已有标签)
    )
  
  print(p_box)
  
  # ==============================================================================
  # 4. 导出为高画质矢量图 PDF (绝不失真，可直接放进发表文章中)
  # ==============================================================================
  ggsave("MegakaryocytesC_Top2_Pathways_Boxplot.pdf", plot = p_box, width = 7, height = 3)
  
  
}

{
  pb_meta <- meta_data[,c(meta_cols,rownames(res_up))]
  
  library(ggplot2)
  library(ggpubr)
  library(dplyr)
  
  # 定义函数
  plot_custom_boxplot <- function(data, cell_name, pathway_name) {
    
    data$group <- factor(
      data$group,
      levels = c("Control", "Disease", "Remission")
    )
    
    plot_df <- data %>%
      filter(celltype == cell_name) %>%
      dplyr::select(group, score = all_of(pathway_name)) %>%
      na.omit()
    
    # 计算 y 轴的最大值，用于动态调整顶部空间
    y_max <- max(plot_df$score, na.rm = TRUE)
    
    my_comparisons <- list(
      c("Control", "Disease"), 
      c("Disease", "Remission")
    )
    
    ggplot(plot_df, aes(x = group, y = score, fill = group)) +
      # 稍微调细线条 (linewidth)，颜色改深灰
      geom_boxplot(
        width = 0.6,
        outlier.shape = NA,
        linewidth = 0.6,
        color = "#333333",
        alpha = 0.85
      ) +
      
      stat_compare_means(
        comparisons = my_comparisons,
        method = "wilcox.test",
        label = "p.signif",
        bracket.size = 0.4,
        step.increase = 0.12
      ) +
      
      scale_fill_manual(values = c(
        "Control" = "#4DBBD5",
        "Disease" = "#E64B35",
        "Remission" = "#00A087"
      )) +
      
      # 核心优化：调整 y 轴起始点，不让它死贴 0
      # 这里的 ylim 设置略小于 0 可以给底部留空
      coord_cartesian(ylim = c(-0.005, y_max * 1.3)) + 
      
      theme_bw() +
      theme(
        aspect.ratio = 0.8,         # 关键：设置纵横比，让图表变宽
        panel.grid.major = element_line(color = "grey95"),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
        legend.position = "none",
        plot.title = element_text(face = "bold", hjust = 0.5, size = 13, margin = margin(b = 10)),
        axis.title.y = element_text(face = "bold", size = 11),
        axis.text = element_text(color = "black", size = 10),
        axis.title.x = element_blank()
      ) +
      
      labs(
        title = paste0(cell_name, "\n", pathway_name),
        y = "Pathway Score"
      )
  }
  
  # --- 使用示例 ---
  
 # [1] "Monocyte"       "NK_cell"        "Megakaryocytes" "T_cell"         "B_cell"        
 # [6] "cDC"            "Plasma_cell"    "pDC"
  # 绘制示例 1
  p1 <- plot_custom_boxplot(pb_meta, "Monocyte", "Arachidonic acid metabolism")
  print(p1)
  
  # 绘制示例 2 
  p2 <- plot_custom_boxplot(pb_meta, "Monocyte", "Fatty acid biosynthesis")
  print(p2)
  
  # 绘制示例 3 
  p3 <- plot_custom_boxplot(pb_meta, "Monocyte", "Nicotinate metabolism")
  print(p3)
  
  # 绘制示例 4
  p4 <- plot_custom_boxplot(pb_meta, "pDC", "Oxidative phosphorylation")
  print(p4)
  
  # 绘制示例 5
  p5 <- plot_custom_boxplot(pb_meta, "cDC", "Oxidative phosphorylation")
  print(p5)
  
  # 绘制示例 6
  p6 <- plot_custom_boxplot(pb_meta, "NK_cell", "Polyamines metabolism")
  print(p6)
  
  library(patchwork)
  
  # 使用 + 或 | 符号进行合并
  # nrow = 2 表示两行，ncol = 3 表示三列
  combined_plot <- (p1 | p2 | p3) / (p4 | p5 | p6) + 
    plot_layout(guides = 'collect') # 自动合并相同的图例（如果有的话）
  
  # 或者更简单的写法：
  # combined_plot <- wrap_plots(p1, p2, p3, p4, p5, p6, ncol = 3)
  
  # 给整个大图添加一个总标题（可选）
  combined_plot <- combined_plot + 
    plot_annotation(
      title = 'Main Pathway Analysis Results',
      theme = theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5))
    )
  
  # 显示结果
  print(combined_plot)
  
  # 保存为 PDF 或 PNG (建议设置好宽高比，否则会挤在一起)
  ggsave("/data3/home/hu/Tongfeng/new_Result/Fig2/combined_pathway_plots.pdf", combined_plot, width = 12, height = 8)
    
  }
  
###绘制通路在细胞表达的小提琴图
{
    library(ggplot2)
    library(dplyr)
    
    # 1. 准备数据
    target_final <- c("Arachidonic acid metabolism", "Fatty acid biosynthesis", 
                      "Nicotinate metabolism", "Polyamines metabolism")
    
    plot_simple_df <- pb_meta %>%
      dplyr::select(celltype, all_of(target_final)) %>%
      pivot_longer(cols = all_of(target_final), names_to = "Pathway", values_to = "Score")
    
    # 2. 绘图
    p1 <- ggplot(plot_simple_df, aes(x = reorder(celltype, -Score, FUN = mean), y = Score)) +
      # 绘制柱状图，使用一种高级的莫兰迪色
      stat_summary(fun = mean, geom = "bar", fill = "#69b3a2", color = "black", width = 0.7) +
      # 叠加误差棒 (SEM)
      stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2, color = "black") +
      # 叠加抖动点，展示样本分布（可选）
      #@geom_jitter(width = 0.1, alpha = 0.2, size = 0.8, color = "#404040") +
      # 分面展示
      facet_wrap(~ Pathway, scales = "free_y", ncol = 2) +
      theme_classic() +
      labs(title = "Metabolic Activity by Cell Type",
           x = NULL, y = "Activity Score") +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, color = "black", size = 9),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = 10),
        panel.grid.major.y = element_line(color = "gray90", linetype = "dashed")
      )
    
    ggsave(plot = p1,"four_pathway_celltype.pdf",width = 8,height = 5)
    
  }
  
####合并
{
    library(tidyverse)
    library(ggplot2)
    library(ggpubr)
    
    # ==========================================
    # 1. 数据预处理与严格趋势筛选
    # ==========================================
    # 设定通路名单与分组顺序
    pathway_names <- colnames(pb_meta)[4:ncol(pb_meta)]
    group_levels <- c("Control", "Int gout remission", "Adv gout remission", "Acute flare")
    pb_meta$group <- factor(pb_meta$group, levels = group_levels)
    
    # 严格筛选逻辑：不仅要显著，还要满足均值严格爬坡 (Control < Int < Adv < Acute)
    strict_increasing_pathways <- pb_meta %>%
      pivot_longer(cols = all_of(pathway_names), names_to = "Pathway", values_to = "Score") %>%
      group_by(Pathway, group) %>%
      summarise(mean_score = mean(Score), .groups = "drop") %>%
      pivot_wider(names_from = group, values_from = mean_score) %>%
      # 严格单调递增过滤
      filter(`Acute flare` > `Adv gout remission` & 
               `Adv gout remission` > `Int gout remission` & 
               `Int gout remission` > `Control`) %>%
      pull(Pathway)
    
    # 取趋势检验 (increasing_trends) 与 严格爬坡名单 的交集
    final_selection <- intersect(increasing_trends$pathway, strict_increasing_pathways)
    
    # 如果筛选出的通路太多，可按需截取前 12 个以防图片太挤
    # final_selection <- head(final_selection, 12)
    
    # ==========================================
    # 2. 绘图 A：疾病级别趋势图 (散点柱状 + 折线)
    # ==========================================
    plot_trend_df <- pb_meta %>%
      dplyr::select(orig.ident, group, all_of(final_selection)) %>%
      pivot_longer(cols = all_of(final_selection), names_to = "Pathway", values_to = "Score") %>%
      mutate(group = factor(group, levels = group_levels))
    
    library(gghalves) # 如果没有请先安装 install.packages("gghalves")
    p1 <- ggplot(plot_trend_df, aes(x = group, y = Score, fill = group)) +
      # 左侧画半个小提琴图，展现分布密度
      geom_half_violin(side = "l", alpha = 0.7, color = NA, trim = FALSE) +
      # 右侧画抖动点，展现原始数据
      geom_half_point(side = "r", aes(color = group), alpha = 0.3, size = 0.6) +
      # 中间穿插一个小箱线图，展现分位数
      geom_boxplot(width = 0.1, outlier.shape = NA, fill = "white", color = "black", alpha = 0.8) +
      # 用一条细虚线连接中位数，提示趋势
      stat_summary(fun = median, geom = "line", aes(group = 1), linetype = "dashed", color = "grey20") +
      facet_wrap(~ Pathway, scales = "free_y", ncol = 4) +
      scale_fill_brewer(palette = "Reds") +
      scale_color_brewer(palette = "Reds") +
      theme_classic() +
      theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))
    
    
    # ==========================================
    # 3. 绘图 B：细胞类型比较图 (莫兰迪绿柱状图)
    # ==========================================
    plot_ct_df <- pb_meta %>%
      dplyr::select(celltype, all_of(final_selection)) %>%
      pivot_longer(cols = all_of(final_selection), names_to = "Pathway", values_to = "Score")
    
    library(ggplot2)
    library(dplyr)
    
    # 美化版 P2
    p2 <- ggplot(plot_ct_df, aes(x = reorder(celltype, -Score, FUN = mean), y = Score)) +
      # 1. 改用条形填充，增加透明度层次
      stat_summary(fun = mean, geom = "bar", aes(fill = Pathway), 
                   color = "black", width = 0.75, size = 0.3, alpha = 0.85) +
      # 2. 优化误差棒，稍微加粗显得更有质感
      stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.25, 
                   color = "black", size = 0.5) +
      # 3. 散点美化：改用带边框的点，点的大小稍微调大，并增加随机感
      geom_jitter(width = 0.2, alpha = 0.2, size = 0.8, color = "#2F4F4F", shape = 16) +
      # 4. 分面美化：允许各通路 Y 轴独立，并调整列数
      facet_wrap(~ Pathway, scales = "free_y", ncol = 2) + 
      # 5. 配色方案：使用内置或自定义高级色（这里推荐一组莫兰迪绿/灰蓝调）
      scale_fill_manual(values = c("#69b3a2", "#778899", "#BC8F8F", "#9ACD32")) +
      # 6. 主题细节深度定制
      theme_classic(base_size = 11) +
      labs(
        title = "Metabolic Activity Across Cell Types",
        subtitle = "Mean score with SEM and individual data points",
        x = NULL,
        y = "Pathway Activity Score"
      ) +
      theme(
        # 坐标轴文字
        axis.text.x = element_text(angle = 45, hjust = 1, color = "black", size = 9),
        axis.text.y = element_text(color = "black"),
        # 分面标题美化（添加浅灰色背景框，更像期刊风格）
        strip.background = element_rect(fill = "#F5F5F5", color = "black", size = 0.5),
        strip.text = element_text(face = "bold.italic", size = 10, margin = margin(t=5, b=5)),
        # 去掉图例（因为标题已经有通路名了）
        legend.position = "none",
        # 增加面板边框，让图看起来更扎实
        panel.border = element_rect(fill = NA, color = "black", size = 0.4),
        plot.title = element_text(face = "bold", size = 14),
        plot.subtitle = element_text(size = 10, color = "grey40")
      )
    
    print(p2)
    
    
    # ==========================================
    # 4. 显示图表
    # ==========================================
    print(p1)
    print(p2)
    
    library(patchwork)
    
    # 建议适当调整 P2 的分面列数，使其宽度与 P1 匹配（例如都设为 4 列）
    p2_rescaled <- p2 + facet_wrap(~ Pathway, scales = "free_y", ncol = 4)
    
    # 使用 patchwork 进行拼接
    combined_plot <- p1 / p2_rescaled + 
      plot_layout(heights = c(1, 1)) + # 设置两张图的高度比例为 1:1
      plot_annotation(
        title = '',
        subtitle = '',
        tag_levels = 'A' # 自动添加 A, B 标签，方便论文引用
      )
    
    # 显示合并图
    print(combined_plot)
    
    # 保存高分辨率图片
    ggsave("/data3/home/hu/Tongfeng/Result/Fig2/combined_metabolism_analysis_B.pdf", 
           combined_plot, width = 12, height = 8)
    
  }
  


#####筛选每个细胞类型的top看看是不是能区分正常和疾病状态
{
  top_sig_meta <- intersect(unique(top5_pathways$gene),gene_counts$Metabolism)
  pb_meta_top <- pb_meta[,c("orig.ident","celltype","group",top_sig_meta)]
  
  library(pROC)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(ggsci) # 引入高级配色
  
  # 1. 数据预处理：将分组重构为 Disease 和 Control
  roc_prep_df <- pb_meta_top %>%
    mutate(outcome = ifelse(group == "Control", 0, 1)) # Control 为 0，其余全为 1
  
  # 获取通路列名
  pathway_names <- colnames(pb_meta_top)[4:ncol(pb_meta_top)]
  
  # 1. 重新计算平滑后的 ROC 数据
  plot_roc_smooth_list <- list()
  
  for(pw in pathway_names) {
    # 关键点：使用 smooth = TRUE 开启平滑处理
    roc_obj <- roc(roc_prep_df$outcome, roc_prep_df[[pw]], 
                   smooth = TRUE, # 开启平滑
                   quiet = TRUE)
    
    tmp_df <- data.frame(
      FPR = 1 - roc_obj$specificities,
      TPR = roc_obj$sensitivities,
      Pathway = pw,
      AUC = as.numeric(roc_obj$auc)
    )
    plot_roc_smooth_list[[pw]] <- tmp_df
  }
  
  plot_roc_smooth_df <- bind_rows(plot_roc_smooth_list)
  
  # 2. 绘制超精美光滑版 ROC
  # 2. 绘制超精美光滑版 ROC (已修复 AI 卡顿隐患)
  p_roc_refined <- ggplot(plot_roc_smooth_df, aes(x = FPR, y = TPR, color = Pathway)) +
    
    # 💡 彻底修复卡顿：直接使用 geom_abline 画单条对角线，移出 aes()。AI 里面只会剩下一条线！
    geom_abline(slope = 1, intercept = 0, color = "grey80", linetype = "dotted", size = 0.6) +
    
    # 绘制光滑曲线
    geom_line(size = 1.2, alpha = 0.9) +
    
    # 使用 NPG 配色
    scale_color_npg() +
    coord_fixed() +
    scale_x_continuous(expand = c(0.01, 0.01), breaks = seq(0, 1, 0.2)) +
    scale_y_continuous(expand = c(0.01, 0.01), breaks = seq(0, 1, 0.2)) +
    labs(
      title = "Receiver Operating Characteristic (ROC)",
      subtitle = "Smoothing: Binormal distribution fitting",
      x = "False Positive Rate",
      y = "True Positive Rate",
      color = "Pathway (AUC)"
    ) +
    
    # 动态生成带 AUC 的图例
    scale_color_npg(labels = function(x) {
      aucs <- plot_roc_smooth_df %>% 
        group_by(Pathway) %>% 
        summarize(a = dplyr::first(AUC), .groups = "drop") 
      paste0(x, " (", sprintf("%.3f", aucs$a[match(x, aucs$Pathway)]), ")")
    }) +
    
    # 极致美化主题
    theme_classic(base_size = 12) +
    theme(
      panel.grid.major = element_line(color = "#F0F0F0", size = 0.2), 
      legend.position = c(0.7, 0.25),
      legend.text = element_text(size = 8.5, face = "italic"), 
      legend.title = element_text(face = "bold", size = 9),
      legend.background = element_blank(), 
      axis.line = element_line(size = 0.6, color = "black"),
      axis.text = element_text(color = "black", face = "bold"),
      plot.title = element_text(face = "bold", size = 16, family = "sans"),
      plot.subtitle = element_text(size = 10, color = "grey40")
    )
  
  print(p_roc_refined)
  
  # 💡 强力补充：使用 AI 最喜欢的标准 PDF 驱动，并关闭 Dingbats 符号字库
  # 这样保存出来的 PDF，在 AI 里不论是文字还是线条都极度丝滑，秒开不卡
  pdf("/data3/home/hu/Tongfeng/new_Result/Fig2/top_metabolism_ROC_C.pdf", 
      width = 8, height = 8, useDingbats = FALSE)
  print(p_roc_refined)
  dev.off()
  
  print(p_roc_refined)
  
  ggsave("/data3/home/hu/Tongfeng/new_Result//Fig2/top_metabolism_ROC_C.pdf", 
         p_roc_refined, width = 8, height = 8)
}

save.image("/data3/home/hu/Tongfeng/new_Result/Fig2/save.RData")

pb_meta_filter <- pb_meta[,c("orig.ident","celltype","group",rownames(res_up))]
saveRDS(pb_meta_filter,"/data3/home/hu/Tongfeng/new_Result//Fig2/pb_meta_filter_new.RDS")


######Commpass数据要求
{
  library(Seurat)
  library(dplyr)
  library(Matrix)
  
  # 确保使用的是最原始的 counts 层
  scRNA_filter <- subset(scRNA,subset = celltype %in% c("Monocyte","pDC","Megakaryocytes","cDC"))
  #（Seurat v5 使用 LayerData，Seurat v4 使用 GetAssayData）
  counts_mat <- LayerData(scRNA_filter, assay = "RNA", layer = "counts")
  
  # ==============================================================================
  # 1. 核心算法：为每个细胞分配一个“伪Bulk（池）”标签
  # ==============================================================================
  meta <- scRNA_filter@meta.data %>% 
    dplyr::select(orig.ident, celltype) %>% 
    dplyr::mutate(cell_id = rownames(.))
  
  # 创建一个新的列，用于存放最终的伪Bulk ID
  meta$pseudo_bulk_id <- NA
  
  # 获取所有“样本 x 细胞类型”的组合
  combinations <- meta %>% 
    dplyr::select(orig.ident, celltype) %>% 
    dplyr::distinct()
  
  set.seed(42) # 保证切分可重复
  
  for(i in 1:nrow(combinations)) {
    sub_sample <- combinations$orig.ident[i]
    sub_celltype <- combinations$celltype[i]
    
    # 筛选出当前样本中该细胞类型的细胞
    idx <- which(meta$orig.ident == sub_sample & meta$celltype == sub_celltype)
    cell_count <- length(idx)
    
    if(cell_count == 0) next
    
    # 定义基础前缀
    prefix <- paste(sub_sample, sub_celltype, sep = "_")
    
    if(cell_count <= 500) {
      # 💡 情况 A：如果细胞数不足或刚满 2000，直接并为一个伪Bulk
      meta$pseudo_bulk_id[idx] <- paste0(prefix, "_pool1")
    } else {
      # 💡 情况 B：如果细胞数超过 2000，计算需要切分成几个池
      # 例如：5000个细胞就会被分成 3 个池（每个约1666个细胞），保证空间分布均匀
      num_pools <- ceiling(cell_count / 500)
      
      # 为了让切分出的池子在生物学上更合理（而不是随机乱揉），
      # 最好基于这些细胞的 UMAP 或 PCA 坐标进行 K-means 聚类切分。
      # 如果不想用聚类，也可以直接用 sample(rep(1:num_pools, length.out = cell_count))
      
      # 这里推荐用 PCA 坐标进行空间切分（保证一个池子里的细胞状态是最相似的）
      cell_names <- meta$cell_id[idx]
      pca_coords <- Embeddings(scRNA, reduction = "pca")[cell_names, 1:5, drop = FALSE]
      
      # 进行 K-means 聚类
      clusters <- kmeans(pca_coords, centers = num_pools, iter.max = 50, nstart = 5)$cluster
      
      meta$pseudo_bulk_id[idx] <- paste0(prefix, "_pool", clusters)
    }
  }
  
  # 将分配好的标签写回 Seurat 对象的元数据中
  scRNA_filter$pseudo_bulk_id <- meta$pseudo_bulk_id
  
  # 检查一下切分后的伪Bulk情况，看看每个伪Bulk包含了多少单细胞
  print(table(scRNA_filter$pseudo_bulk_id))
  print(unique(scRNA_filter$pseudo_bulk_id))
  
  
  # ==============================================================================
  # 2. 矩阵聚合（等同于把各个池子里的细胞 Counts 加和）
  # ==============================================================================
  cat("正在进行伪Bulk表达量加和聚合...\n")
  
  # 使用 Seurat 自带的高效聚合函数（Seurat v5 推荐 AggregateExpression）
  bulk_res <- AggregateExpression(
    scRNA_filter, 
    group.by = "pseudo_bulk_id", 
    assays = "RNA", 
    return.seurat = FALSE
  )
  
  # 提取聚合后的 Counts 矩阵
  compass_input_matrix <- bulk_res$RNA
  
  
  # ==============================================================================
  # 3. 构建对应的 Compass 样本信息表 (Metadata)
  # ==============================================================================
  library(dplyr)
  library(tidyr)
  
  compass_metadata <- data.frame(
    pseudo_bulk_id = colnames(compass_input_matrix)
  ) %>%
    # 💡 完美解决细胞类型带下划线的问题：
    # 既然在前面生成 ID 时是用“下划线”拼接的，我们可以利用正则表达式限制“只切分倒数第一个下划线”
    # 或者直接从 pseudo_bulk_id 中用更安全的正则提取：
    # 1. 提取前两个字符作为样本名 (如 A1, B3)
    mutate(orig.ident = substr(pseudo_bulk_id, 1, 2)) %>%
    # 2. 提取最后的 pool 编号
    mutate(pool_suffix = sub(".*_(pool[0-9]+)$", "\\1", pseudo_bulk_id)) %>%
    # 3. 擦除两头，剩下的就是中间完整的细胞类型名（哪怕带下划线如 NK_cell 也会完美保留）
    mutate(celltype = sub("^[A-C][1-3]_(.*)_pool[0-9]+$", "\\1", pseudo_bulk_id)) %>%
    
    # 💡 修复函数名：将 case_by 修改为标准的 case_when
    dplyr::mutate(
      Group = dplyr::case_when(
        grepl("^A", orig.ident) ~ "Control",    # A 组为对照
        grepl("^B", orig.ident) ~ "Gout",       # B 组为痛风
        grepl("^C", orig.ident) ~ "Remission",  # C 组为缓解
        TRUE ~ "Unknown"
      )
    )
  
  # 检查一下解析出来的元数据，确保没有错位或 NA
  print(head(compass_metadata))
  
  
  
  # ==============================================================================
  # 4. 保存输出文件（准备提交给 Compass）
  # ==============================================================================
  # 注意：Compass 需要输入标准的 CSV 或 TSV 矩阵
  output_dir <- "/data3/home/hu/Tongfeng/new_Result/Fig2/Compass2/"
  if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  # 保存表达量矩阵 (建议存为 tsv，因为基因名或样本名可能包含逗号)
  write.table(
    as.matrix(compass_input_matrix), 
    file = paste0(output_dir, "compass_expression_counts.tsv"), 
    sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE
  )
  
  # 保存元数据信息
  write.table(
    compass_metadata, 
    file = paste0(output_dir, "compass_metadata.tsv"), 
    sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE
  )
  
  cat("🎉 恭喜！Compass 输入文件已成功生成在：", output_dir, "\n")
  
}

######Commpass后续分析
{
  #####提取伪bulk表达数据
  sample_name <- paste0("sample",0:73)
  dir_name <- paste0("/data3/home/hu/Tongfeng/new_Result/Fig2/Compass2/_tmp/",sample_name)
  
  reactions_list <- list()
  for (i in 1:length(dir_name)) {
    temp_data <- read.table(paste0(dir_name[i],"/reactions.txt"))
    reactions_list[[i]] <- temp_data
      }
  
  
  commpass_meta <- read.table("/data3/home/hu/Tongfeng/new_Result/Fig2/Compass2/compass_metadata.tsv",
                           header = T,sep = "\t")
  
  names(reactions_list) <- commpass_meta$pseudo_bulk_id
  
  ####合并reaction
  reactions_data <- do.call(cbind,reactions_list)
  head(reactions_data)
  

  reaction_metadata <- read.csv("/data3/home/hu/Tongfeng/new_Result/Fig2/Commpass/reaction_metadata.csv")
  
  head(reaction_metadata[,1:5])
  head(commpass_meta[,1:5])
  head(reactions_data[,1:5])
  
  ####做差异
  {
    library(dplyr)
    library(tidyr)
    library(tibble)
    library(effsize) # 用于计算 Cohen's d，如果没有请先 install.packages("effsize")
    
    # 1. 基础数据清洗与对齐
    commpass_meta_clean <- commpass_meta %>%
      mutate(orig.ident = substr(pseudo_bulk_id, 1, 2)) %>%
      mutate(celltype = sub("^[A-C][1-3]-(.*)-pool[0-9]+$", "\\1", pseudo_bulk_id)) %>%
      mutate(match_id = gsub("-", ".", pseudo_bulk_id))
    
    # 2. 将惩罚分转换为得分 (Score)，数值越大活性越高
    reaction_scores <- -log(reactions_data + 1)
    
    # ==============================================================================
    # 3. 核心计算函数：对指定细胞类型计算 W.test 和 Cohen's d (保持逻辑，优化无变异处理)
    # ==============================================================================
    run_compass_wilcoxon <- function(scores_mat, meta_df, target_celltype) {
      
      # 筛选目标细胞和对比组 (Gout vs Control)
      sub_meta <- meta_df %>% 
        filter(celltype == target_celltype & Group %in% c("Control", "Gout"))
      
      sub_scores <- scores_mat[, sub_meta$match_id, drop = FALSE]
      
      # 提取分组标签
      groups <- sub_meta$Group
      gout_idx <- which(groups == "Gout")
      ctrl_idx <- which(groups == "Control")
      
      # 循环对每一个反应进行计算
      res_list <- lapply(rownames(sub_scores), function(reaction) {
        vec <- as.numeric(sub_scores[reaction, ])
        gout_vec <- vec[gout_idx]
        ctrl_vec <- vec[ctrl_idx]
        
        # 💡 优化：如果两组都完全无变异，不返回 NULL（防止漏掉点），直接给一个不显著的默认值占位
        if(sd(gout_vec) == 0 & sd(ctrl_vec) == 0) {
          return(data.frame(reaction_id = reaction, cohens_d = 0, P.Value = 1, stringsAsFactors = FALSE))
        }
        
        # Wilcoxon 秩和检验
        w_test <- wilcox.test(gout_vec, ctrl_vec)
        p_val <- w_test$p.value
        if(is.na(p_val)) p_val <- 1
        
        # 计算 Cohen's d 效应量
        d_res <- effsize::cohen.d(gout_vec, ctrl_vec)
        cohens_d <- d_res$estimate
        
        data.frame(
          reaction_id = reaction,
          cohens_d = cohens_d,
          P.Value = p_val,
          stringsAsFactors = FALSE
        )
      })
      
      # 合并结果
      res_df <- do.call(rbind, res_list) %>%
        mutate(adjusted_pval = p.adjust(P.Value, method = "BH")) %>%
        mutate(reaction_no_direction = sub("_(pos|neg)$", "", reaction_id)) %>%
        # 💡 核心修改 [1/2]：去除反应名末尾可能因软件版本导致的 "E"（如 RE3520E -> RE3520）
        mutate(base_id = sub("E$", "", reaction_no_direction))
      
      return(res_df)
    }
    
    # 💡 运行分析：拿到单核细胞的差异分析结果
    mono_w_results <- run_compass_wilcoxon(reaction_scores, commpass_meta_clean, "Megakaryocytes")
    
    # ==============================================================================
    # 4. 融合通路元数据 (通过清洗后的 base_id 强效匹配)
    # ==============================================================================
    # 💡 核心修改 [2/2]：对元数据做相同的去 "E" 清洗，防止由于末尾字母不同而被 inner_join 抹杀
    reaction_metadata_clean <- reaction_metadata %>%
      mutate(base_id = sub("E$", "", reaction_no_direction)) %>%
      # 去除可能因为去E导致的重复行，保留唯一映射
      distinct(base_id, .keep_all = TRUE)
    
    # 使用清洗后的 base_id 进行无损融合
    mono_w_annotated <- mono_w_results %>%
      inner_join(reaction_metadata_clean, by = "base_id")
    
    # 检查现在融合后的结果
    head(mono_w_annotated)
    
    # 💡 验证：看看现在 Arachidonic acid metabolism 有多少行了？
    cat("检查花生四烯酸代谢反应数：", 
        nrow(filter(mono_w_annotated, subsystem == "Arachidonic acid metabolism")), "\n")
    
    
  }
  ###可视化绘图
  {
    library(ggplot2)
    library(dplyr)
    library(ggrepel)
    library(grid)
    
    # ==============================================================================
    # 核心函数定义：一键绘制单个通路的 Compass 火山图
    # ==============================================================================
    plot_compass_adaptive <- function(all_data, target_subsystem, title_name, dot_color = "#E15759", save_path = NULL) {
      
      # 1. 过滤指定通路
      plot_data <- all_data %>% 
        filter(subsystem == target_subsystem) %>%
        mutate(P.Value = ifelse(P.Value < 1e-12, 1e-12, P.Value))
      
      if(nrow(plot_data) == 0) {
        stop(paste0("😭 错误：在数据中未找到通路 '", target_subsystem, "'"))
      }
      
      # 2. 💡 智能标签：不再卡死 P 值，而是直接捞出该通路内部响应最强（效应量绝对值最大）的前 5 个反应
      labeled_data <- plot_data %>%
        arrange(desc(abs(cohens_d))) %>%
        head(5) 
      
      # 3. 绘图
      p <- ggplot(plot_data, aes(x = cohens_d, y = -log10(P.Value))) +
        # 辅助线线放淡，避免尴尬
        geom_vline(xintercept = 0, linetype = "dashed", color = "#BBBBBB", size = 0.4) +
        geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "#E15759", size = 0.4) + # 0.05 显色参考线
        
        geom_point(color = dot_color, alpha = 0.8, size = 2.5) +
        
        # 智能添加标签
        geom_text_repel(
          data = labeled_data, aes(label = reaction_name),
          size = 2.6, color = "#222222", box.padding = unit(0.5, "lines"),
          point.padding = unit(0.3, "lines"), segment.color = "#888888",
          segment.size = 0.3, max.overlaps = 30
        ) +
        
        scale_x_continuous(limits = c(-2.2, 2.2)) +
        theme_bw(base_size = 13) +
        theme(
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          plot.title = element_text(face = "bold", size = 16, hjust = 0.5, margin = margin(b = 15)),
          axis.title = element_text(size = 12), axis.text = element_text(color = "#333333"),
          plot.margin = margin(t = 10, r = 10, b = 40, l = 10, unit = "pt") 
        ) +
        labs(
          x = "Cohen's d",
          y = "-log10 (Wilcoxon p-value)",
          title = title_name
        )
      
      # 4. 底层涂鸦：双向粗箭头
      g <- ggplotGrob(p)
      print(p)
      g$layout$clip[g$layout$name == "panel"] <- "off"
      
      arrow_grob <- gTree(children = gList(
        segmentsGrob(x0 = 0.55, y0 = -0.12, x1 = 0.95, y1 = -0.12, default.units = "npc",
                     arrow = arrow(length = unit(0.2, "cm"), ends = "last", type = "closed"),
                     gp = gpar(col = "#E92E87", lwd = 4)),
        textGrob("Gout", x = 0.75, y = -0.18, default.units = "npc", gp = gpar(fontsize = 11, fontface = "bold", col = "#E92E87")),
        
        segmentsGrob(x0 = 0.45, y0 = -0.12, x1 = 0.05, y1 = -0.12, default.units = "npc",
                     arrow = arrow(length = unit(0.2, "cm"), ends = "last", type = "closed"),
                     gp = gpar(col = "#348C73", lwd = 4)),
        textGrob("Control", x = 0.25, y = -0.18, default.units = "npc", gp = gpar(fontsize = 11, fontface = "bold", col = "#348C73"))
      ))
      
      g <- gtable::gtable_add_grob(g, arrow_grob, t = g$layout[g$layout$name == "panel", "t"], 
                                   l = g$layout[g$layout$name == "panel", "l"])
      
      if (!is.null(save_path)) {
        pdf(save_path, width = 6, height = 6.5, useDingbats = FALSE)
        grid.draw(g)
        dev.off()
      } else {
        grid.newpage()
        grid.draw(g)
      }
    }
    
  }
  # 检查一下你的数据里到底有哪些完整的 Subsystem 名字，防止拼写错误
  print(unique(mono_w_annotated$subsystem))
  
  
  
  # 1. 绘制糖酵解通路火山图 (指定官方神秘紫灰配色 #695D73)
  plot_compass_adaptive(
    all_data = mono_w_annotated, 
    target_subsystem = "Fatty acid synthesis", 
    title_name = "Arachidonic acid metabolism", 
    dot_color = "#695D73",
    save_path = "/data3/home/hu/Tongfeng/new_Result/Fig2/Pathway_Glycolysis_Volcano.pdf"
  )
  
  # 2. 绘制三羧酸循环 TCA 通路火山图 (指定官方复古燕麦色 #D3A991)
  plot_compass_adaptive(
    all_data = mono_w_annotated, 
    target_subsystem = "Oxidative phosphorylation", 
    title_name = "Oxidative phosphorylation", 
    dot_color = "#695D73",
    save_path = "/data3/home/hu/Tongfeng/new_Result/Fig2/Megakaryocytes_Pathway_Oxidative phosphorylation.pdf"
  )
  
  # 3. 绘制脂肪酸氧化通路火山图 (指定官方优雅深蓝 #040772)
  plot_compass_volcano(
    all_data = mono_gout_annotated, 
    target_subsystem = "Fatty acid oxidation", 
    title_name = "Fatty Acid Oxidation", 
    dot_color = "#040772",
    save_path = "/data3/home/hu/Tongfeng/new_Result/Fig2/Pathway_FAO_Volcano.pdf"
  )
  
  
  {
    library(dplyr)
    library(tidyr)
    library(tibble)
    
    # 1. 彻底修正 commpass_meta 的细胞类型切分错误（兼容带下划线、短横线的细胞名）
    commpass_meta_clean <- commpass_meta %>%
      mutate(orig.ident = substr(pseudo_bulk_id, 1, 2)) %>%
      mutate(pool_suffix = sub(".*-(pool[0-9]+)$", "\\1", pseudo_bulk_id)) %>%
      # 提取中间完整的细胞类型名
      mutate(celltype = sub("^[A-C][1-3]-(.*)-pool[0-9]+$", "\\1", pseudo_bulk_id)) %>%
      # 将 ID 中的短横线全部替换为点号，用于完美匹配表达矩阵的列名
      mutate(match_id = gsub("-", ".", pseudo_bulk_id))
    
    # 2. 将反应惩罚分（Penalty）转换为活性得分（Score）
    # 数值越大代表活性越高，顺便过滤掉一些全零或无效反应
    reaction_scores <- -log(reactions_data + 1)
    
    library(limma)
    
    # 1. 筛选出 Monocyte 的样本
    mono_meta <- commpass_meta_clean %>% filter(celltype == "Monocyte")
    # 提取对应的得分矩阵，并确保列名顺序一致
    mono_scores <- reaction_scores[, mono_meta$match_id]
    
    # 2. 构建设计矩阵（分析 Control vs Gout vs Remission）
    group_factor <- factor(mono_meta$Group, levels = c("Control", "Gout", "Remission"))
    design <- model.matrix(~ 0 + group_factor)
    colnames(design) <- c("Control", "Gout", "Remission")
    
    fit <- lmFit(mono_scores, design)
    
    # 3. 制定对比矩阵：看痛风组相比对照上调，以及缓解组相比痛风组的变化
    contrast_matrix <- makeContrasts(
      Gout_vs_Control = Gout - Control,
      Remission_vs_Gout = Remission - Gout,
      levels = design
    )
    fit2 <- contrasts.fit(fit, contrast_matrix)
    fit2 <- eBayes(fit2)
    
    # 4. 导出痛风组最特异变化的代谢反应
    mono_gout_diff <- topTable(fit2, coef = "Gout_vs_Control", number = Inf) %>%
      rownames_to_column("reaction_id") %>%
      # 💡 核心：去除 ID 末尾的 _pos 或 _neg，以便跟你的 reaction_metadata 关联
      mutate(reaction_no_direction = sub("_(pos|neg)$", "", reaction_id))
    
    head(mono_gout_diff)
    
    # 1. 将差异反应与代谢大类（subsystem）和反应名称合并
    mono_gout_annotated <- mono_gout_diff %>%
      inner_join(reaction_metadata, by = "reaction_no_direction")
    
    # 2. 筛选显著变化的差异反应（根据 P值 和 logFC）
    # logFC > 0 代表在痛风组中该代谢反应活性显著上调
    sig_up_reactions <- mono_gout_annotated %>% filter(P.Value < 0.05 & logFC > 0.2)
    sig_down_reactions <- mono_gout_annotated %>% filter(P.Value < 0.05 & logFC < -0.2)
    
    # 3. 计算哪些 Subsystem（代谢通路）富集得最显著（超几何检验/富集分析）
    # 看看在显著上调的反应里，哪些通路占的比例远超背景
    library(dplyr)
    
    # 1. 重新筛选上调反应
    sig_up_reactions <- mono_gout_annotated %>% 
      filter(P.Value < 0.05 & logFC > 0.2)
    
    # 2. 用 length() 替代 n() 计算每组的数量，彻底根除报错
    enrich_up <- sig_up_reactions %>%
      group_by(subsystem) %>%
      summarise(
        Signif_Count = length(logFC),  # 💡 用 length(任意列) 替代 n()，效果完全一样且永不报错
        Mean_logFC = mean(logFC, na.rm = TRUE)
      ) %>%
      arrange(desc(Signif_Count))
    
    # 3. 打印查看富集结果
    print(head(enrich_up, 10))
    
    {
      library(ggplot2)
      library(dplyr)
      library(ggrepel) # 用于防止文字标签重叠
      
      # ==============================================================================
      # 1. 准备火山图数据并定义阈值
      # ==============================================================================
      # 定义阈值：P值小于 0.05，且绝对值 logFC 大于 0.2
      p_threshold <- 0.05
      lfc_threshold <- 0.2
      
      volcano_data <- mono_gout_annotated %>%
        mutate(
          # 判定反应的上下调状态
          Significance = case_when(
            P.Value < p_threshold & logFC > lfc_threshold ~ "Up-regulated",
            P.Value < p_threshold & logFC < -lfc_threshold ~ "Down-regulated",
            TRUE ~ "Not Significant"
          )
        )
      
      # ==============================================================================
      # 2. 筛选需要特异性标注文字（Label）的核心反应
      # ==============================================================================
      # 挑选 P 值最显著或 logFC 最大的前 10 个反应显示名字，防止画面太乱
      top_labels <- volcano_data %>%
        filter(Significance != "Not Significant") %>%
        arrange(P.Value) %>%
        head(10)
      
      # ==============================================================================
      # 3. 使用 ggplot2 绘制高颜值火山图
      # ==============================================================================
      p_volcano <- ggplot(volcano_data, aes(x = logFC, y = -log10(P.Value))) +
        # 💡 技巧：先画不显著的点，用极淡的灰色作为背景
        geom_point(data = filter(volcano_data, Significance == "Not Significant"),
                   color = "#E0E0E0", alpha = 0.6, size = 1.2) +
        # 再画显著上调的点（高级冷调红）
        geom_point(data = filter(volcano_data, Significance == "Up-regulated"),
                   color = "#E15759", alpha = 0.8, size = 1.8) +
        # 再画显著下调的点（高级深灰蓝）
        geom_point(data = filter(volcano_data, Significance == "Down-regulated"),
                   color = "#4E79A7", alpha = 0.8, size = 1.8) +
        
        # 💡 技巧：添加辅助阈值线（虚线，且颜色淡化，不抢戏）
        geom_vline(xintercept = c(-lfc_threshold, lfc_threshold), linetype = "dashed", color = "#999999", size = 0.4) +
        geom_hline(yintercept = -log10(p_threshold), linetype = "dashed", color = "#999999", size = 0.4) +
        
        # 💡 技巧：使用 ggrepel 完美标注核心反应名称
        geom_text_repel(
          data = top_labels,
          aes(label = reaction_name),
          size = 2.8,
          color = "#222222",
          box.padding = unit(0.4, "lines"),
          point.padding = unit(0.3, "lines"),
          segment.color = "#777777", # 引导线颜色
          segment.size = 0.3,
          max.overlaps = 15
        ) +
        
        # 💡 技巧：手动映射颜色图例
        scale_color_manual(values = c("Up-regulated" = "#E15759", "Down-regulated" = "#4E79A7", "Not Significant" = "#E0E0E0")) +
        
        # 主题美化：全白清爽背景
        theme_bw(base_size = 12) +
        theme(
          panel.grid.major = element_blank(), # 去掉网格线
          panel.grid.minor = element_blank(),
          axis.line = element_line(color = "#333333"),
          axis.text = element_text(color = "#333333"),
          axis.title = element_text(face = "bold"),
          plot.title = element_text(face = "bold", hjust = 0.5, size = 13),
          legend.position = "right",
          legend.title = element_blank()
        ) +
        
        # 标签
        labs(
          x = "Log2 Fold Change (Gout / Control)",
          y = "-Log10 (P-Value)",
          title = "Volcano Plot of Metabolic Reactions in Monocytes"
        )
      
      # ==============================================================================
      # 4. 保存为 PDF 大图
      # ==============================================================================
      output_path <- "/data3/home/hu/Tongfeng/new_Result/Fig2/Monocyte_Volcano_Plot.pdf"
      print(p_volcano)
      ggsave(output_path, plot = p_volcano, width = 6.5, height = 5.5)
      
      cat("🎉 火山图已完美保存至：", output_path, "\n")
      
    }
    
  }
  
  
  {
    library(tidyverse)
    library(rstatix)
    library(ggrepel)
    library(ggsci)
    
    # 1. 数据转换 (Penalty -> Score)
    # COMPASS 逻辑：score = -log(penalty + 1)
    reaction_scores <- -log(reactions_data + 1)
    # 归一化：减去最小值以消除背景
    reaction_scores <- reaction_scores - min(as.matrix(reaction_scores))
    
    # 计算 stats_results 之前，先剔除方差为 0 的反应
    # 这是一个更推荐的做法
    keep_reactions <- apply(reaction_scores, 1, function(x) var(x) > 0)
    reaction_scores_filtered <- reaction_scores[keep_reactions, ]
    
    
    # 2. 差异分析准备
    # 将矩阵转置，并合并分组信息
    plot_df <- as.data.frame(t(reaction_scores_filtered)) %>%
      rownames_to_column("Pool_ID_new") %>%
      left_join(pools_final_df %>% select(Pool_ID_new, Major_Group), by = "Pool_ID_new") %>%
      mutate(Status = ifelse(Major_Group == "Control", "Control", "Disease"))
    
    
    {
      library(tidyverse)
      
      # 计算每个 Group 的均值，并转置方便观察
      group_comparison <- plot_df %>%
        # 移除不需要计算的列
        select(-Pool_ID_new, -Major_Group) %>% 
        # 按 Major_Group 分组计算所有反应的均值
        group_by(Status) %>%
        summarise(across(everything(), mean, .groups = "drop")) %>%
        # 转置表格：行变成反应，列变成组名
        pivot_longer(-Status, names_to = "reaction_id", values_to = "mean_score") %>%
        pivot_wider(names_from = Status, values_from = mean_score)
      
      # 关联反应名称，方便阅读
      group_comparison <- group_comparison %>%
        mutate(match_id = str_remove(reaction_id, "_pos$|_neg$")) %>%
        left_join(reaction_metadata %>% select(reaction_no_direction, reaction_name, subsystem), 
                  by = c("match_id" = "reaction_no_direction"))
      
      # 查看前几行
      head(group_comparison)
      
      }
    
    # 3. 计算每个反应的 Wilcoxon P值 和 Cohen's d
    # 定义一个函数批量处理
    reaction_ids <- rownames(reaction_scores_filtered)
    
    stats_results <- lapply(reaction_ids, function(rid) {
      # 提取单反应数据
      tmp <- plot_df %>% select(all_of(rid), Status)
      colnames(tmp)[1] <- "val"
      
      # Wilcoxon 检验
      w_test <- wilcox.test(val ~ Status, data = tmp)
      
      # Cohen's d (效应量)
      # 逻辑：(Mean_Disease - Mean_Control) / SD_pooled
      # 显式指定比较顺序：[1] 疾病组, [2] 控制组。计算结果将变为：Disease - Control
      d_val <- tmp %>% cohens_d(val ~ Status, ref.group = "Disease")
      
      data.frame(
        reaction_id = rid,
        p_val = w_test$p.value,
        cohens_d = d_val$effsize
      )
    }) %>% bind_rows()
    
    # 矫正 P 值
    stats_results$p_adj <- p.adjust(stats_results$p_val, method = "BH")
    
    # 4. 合并元数据
    # 处理 ID：去掉 _pos/_neg 后缀以便匹配元数据
    stats_results <- stats_results %>%
      mutate(match_id = str_remove(reaction_id, "_pos$|_neg$")) %>%
      left_join(reaction_metadata, by = c("match_id" = "reaction_no_direction"))
    
    # 5. 绘制火山图 (以 Glycolysis 为例)
    # 你可以修改 filter 条件来画不同的 Subsystem
    target_subsystem <- "Steroid metabolism" # 确保名称与 reaction_metadata 中一致
    
    plot_data <- stats_results %>%
      filter(subsystem == target_subsystem) %>%
      filter(!is.na(p_adj))
    
    library(ggplot2)
    library(ggrepel)
    library(dplyr)
    
    # 1. 预处理标签数据：只选出正向和负向 Top 5 的显著反应 (根据 Cohen's d 绝对值和 P值)
    # 假设 plot_data 已经包含了 cohens_d 和 p_adj
    plot_data$label_name <- ""
    plot_data$label_name[c(1,3:5,31,46,48)] <- plot_data$reaction_name[c(1,3:5,31,46,48)]
    
    # 2. 计算 X 轴对称范围，确保 0 点在正中间
    max_x <- max(abs(plot_data$cohens_d), na.rm = TRUE) * 1.1
    
    # 2. 绘制符合参考图风格的火山图
    p_volcano_styled_1 <- ggplot(plot_data, aes(x = cohens_d, y = -log10(p_adj))) +
      # A. 辅助线：0刻度细虚线
      geom_vline(xintercept = 0, linetype = "dashed", color = "grey30", size = 0.4) +
      geom_hline(yintercept = 1, linetype = "dashed", color = "grey30", size = 0.4) + # 对应 P=0.1 或 0.05
      
      # B. 散点：改为较小的实心点，颜色使用参考图中的棕褐色/土黄色系
      geom_point(aes(color = cohens_d > 0), size = 1.2, alpha = 0.8) +
      
      # C. 标签层：复现参考图中的细指引线和靠右对齐风格
      geom_text_repel(aes(label = label_name),
                      size = 3,
                      min.segment.length = 0,      # 始终显示指引线
                      segment.color = "grey20",    # 指引线颜色
                      segment.size = 0.3,          # 极细线
                      nudge_x = -0.25,              # 标签向右偏移
                      direction = "y",             # 标签主要在纵向排布
                      hjust = 0,                   # 左对齐文字
                      force = 2,
                      max.overlaps = Inf,
                      fontface = "plain") +
      
      # D. 配色方案：参考图使用的是偏土黄/棕色的调色板
      scale_color_manual(values = c("TRUE" = "#A67C52", "FALSE" = "#170084")) +
      
      # E. 坐标轴范围
      scale_x_continuous(limits = c(-max_x, max_x), expand = c(0, 0)) +
      scale_y_continuous(limits = c(0, max(-log10(plot_data$p_adj)) * 1.1), expand = c(0, 0)) +
      
      # F. 主题美化：重点复现顶部色块
      theme_classic(base_size = 12) +
      labs(x = NULL, y = NULL) + # 坐标轴标签由主题或外层控制
      theme(
        # 移除图例
        legend.position = "none",
        # 坐标轴线加粗
        axis.line = element_line(color = "black", size = 0.6),
        axis.text = element_text(color = "black", size = 10),
        # 关键：利用分面模拟顶部的色块标题
        strip.background = element_rect(fill = "#A67C52", color = NA),
        strip.text = element_text(color = "white", face = "bold", size = 12),
        plot.margin = margin(t = 10, r = 20, b = 10, l = 10)
      ) +
      # 使用分面来生成那个棕色标题栏
      facet_grid(. ~ subsystem) 
    
    print(p_volcano_styled_1)
    
    ###再绘制一个火山图
    {
      # 你可以修改 filter 条件来画不同的 Subsystem
      target_subsystem <- "Urea cycle" # 确保名称与 reaction_metadata 中一致
      
      plot_data <- stats_results %>%
        filter(subsystem == target_subsystem) %>%
        filter(!is.na(p_adj))
      
      library(ggplot2)
      library(ggrepel)
      library(dplyr)
      
      # 1. 预处理标签数据：只选出正向和负向 Top 5 的显著反应 (根据 Cohen's d 绝对值和 P值)
      # 假设 plot_data 已经包含了 cohens_d 和 p_adj
      plot_data$label_name <- ""
      plot_data$label_name[c(1:5,10,14,16)] <- plot_data$reaction_name[c(1:5,10,14,16)]
      
      # 2. 计算 X 轴对称范围，确保 0 点在正中间
      max_x <- max(abs(plot_data$cohens_d), na.rm = TRUE) * 1.1
      
      # 2. 绘制符合参考图风格的火山图
      p_volcano_styled_2 <- ggplot(plot_data, aes(x = cohens_d, y = -log10(p_adj))) +
        # A. 辅助线：0刻度细虚线
        geom_vline(xintercept = 0, linetype = "dashed", color = "grey30", size = 0.4) +
        geom_hline(yintercept = 1, linetype = "dashed", color = "grey30", size = 0.4) + # 对应 P=0.1 或 0.05
        
        # B. 散点：改为较小的实心点，颜色使用参考图中的棕褐色/土黄色系
        geom_point(aes(color = cohens_d > 0), size = 1.2, alpha = 0.8) +
        
        # C. 标签层：复现参考图中的细指引线和靠右对齐风格
        geom_text_repel(aes(label = label_name),
                        size = 3,
                        min.segment.length = 0,      # 始终显示指引线
                        segment.color = "grey20",    # 指引线颜色
                        segment.size = 0.3,          # 极细线
                        nudge_x = -0.25,              # 标签向右偏移
                        direction = "y",             # 标签主要在纵向排布
                        hjust = 0,                   # 左对齐文字
                        force = 2,
                        max.overlaps = Inf,
                        fontface = "plain") +
        
        # D. 配色方案：参考图使用的是偏土黄/棕色的调色板
        scale_color_manual(values = c("TRUE" = "#A67C52", "FALSE" = "#170084")) +
        
        # E. 坐标轴范围
        scale_x_continuous(limits = c(-max_x, max_x), expand = c(0, 0)) +
        scale_y_continuous(limits = c(0, max(-log10(plot_data$p_adj)) * 1.1), expand = c(0, 0)) +
        
        # F. 主题美化：重点复现顶部色块
        theme_classic(base_size = 12) +
        labs(x = NULL, y = NULL) + # 坐标轴标签由主题或外层控制
        theme(
          # 移除图例
          legend.position = "none",
          # 坐标轴线加粗
          axis.line = element_line(color = "black", size = 0.6),
          axis.text = element_text(color = "black", size = 10),
          # 关键：利用分面模拟顶部的色块标题
          strip.background = element_rect(fill = "#A67C52", color = NA),
          strip.text = element_text(color = "white", face = "bold", size = 12),
          plot.margin = margin(t = 10, r = 20, b = 10, l = 10)
        ) +
        # 使用分面来生成那个棕色标题栏
        facet_grid(. ~ subsystem) 
      
      print(p_volcano_styled_2)
      
      library(patchwork)
      
      # 使用 | 进行横向拼接
      p_combined <- p_volcano_styled_1 | p_volcano_styled_2
      
      # 添加全局标注
      p_final <- p_combined + 
        plot_annotation(
          title = "Differential Metabolic Landscape",
          subtitle = "Comparison of two metabolic subsystems",
          tag_levels = 'A' # 自动添加 A, B 标签
        ) & 
        theme(plot.title = element_text(face = "bold", size = 16))
      
      print(p_final)
      
      ggsave("/data3/home/hu/Tongfeng/Result/Fig2/compass_Volcano_map.pdf", 
             p_final, width = 8, height = 4)
      
    }
    
    
  }
  
  
  ####森林散点图
  {
    library(ggplot2)
    library(dplyr)
    library(ggsci)
    
    mapping <- list(
      "Ketone bodies metabolism" = c(
        "Butanoate metabolism",
        "Fatty acid oxidation",
        "Fatty acid synthesis",
        "Pyruvate metabolism",
        "Citric acid cycle"
      ),
      
      "Steroid metabolism" = c(
        "Steroid metabolism",
        "Cholesterol metabolism",
        "Squalene and cholesterol synthesis",
        "Androgen and estrogen synthesis and metabolism",
        "Bile acid synthesis"
      ),
      
      "DNA synthesis" = c(
        "Pyrimidine synthesis",
        "Purine synthesis",
        "Nucleotide interconversion",
        "Nucleotide salvage pathway",
        "Folate metabolism"
      ),
      
      "Nicotinate and nicotinamide metabolism" = c(
        "NAD metabolism",
        "Nucleotide salvage pathway",
        "Nucleotide interconversion",
        "Tryptophan metabolism"
      )
    )
    
    # 1. 准备绘图数据
    # 确保 stats_results 包含：subsystem, cohens_d, p_adj
    plot_forest_df <- stats_results %>%
      filter(!is.na(subsystem)) %>%
      # 标记是否显著 (BH-adjusted p < 0.1)
      mutate(is_signif = ifelse(p_adj < 0.1, "p < 0.1", "NS")) %>%
      # 按照中位数排序，让图看起来有阶梯感
      group_by(subsystem) %>%
      mutate(median_d = median(cohens_d)) %>%
      ungroup() %>%
      mutate(subsystem = reorder(subsystem, median_d))
    
    plot_forest_df <- plot_forest_df[plot_forest_df$subsystem %in% unlist(mapping),]
    # 2. 绘制森林散点图
    # 风格一：学术极简
    p_forest_science <- ggplot(plot_forest_df, aes(x = cohens_d, y = subsystem)) +
      # 1. 只有一条清晰的 0 刻度线
      geom_vline(xintercept = 0, color = "black", size = 0.5) +
      
      # 2. 散点：稍微调大显著点，不显著的点设为极淡的灰色
      geom_jitter(aes(color = cohens_d, alpha = is_signif, size = is_signif), 
                  height = 0.25, width = 0) +
      
      # 3. 配色：经典的红蓝配 (或者你之前的绿紫)
      scale_color_gradient2(low = "#83B0D4", mid = "#F7F7F7", high = "#AC4B4F", 
                            midpoint = 0, name = "Effect Size") +
      
      # 显著性控制：非显著的点近乎隐形
      scale_alpha_manual(values = c("NS" = 0.05, "p < 0.1" = 0.4)) +
      scale_size_manual(values = c("NS" = 0.25, "p < 0.1" = 0.9)) +
      
      # 4. 纯白主题
      theme_minimal(base_size = 11) +
      labs(x = "Cohen's d", y = NULL, title = "Metabolic Activity Landscape") +
      theme(
        panel.grid = element_blank(), # 移除所有网格
        axis.line.x = element_line(color = "black", size = 0.6), # 只保留底部轴线
        axis.text.y = element_text(color = "black", size = 9),
        axis.text.x = element_text(color = "black", face = "bold"),
        legend.position = "right",
        plot.title = element_text(face = "bold", size = 14, hjust = 0)
      )
    
    print(p_forest_science)
    
    
    ggsave("/data3/home/hu/Tongfeng/Result/Fig2/Commpass_Forest_Scatter_Map.pdf", 
           p_forest_science, width = 8, height = 5)
    
    
    
  }
  
  #### 挑选查看是否存在趋势
  {
    library(tidyverse)
    
    # 设定疾病进展的级别顺序
    group_order <- c("Control", "Int.gout.remission", "Adv.gout.remission", "Acute.flare")
    
    # 整合得分与分组信息
    trend_data <- as.data.frame(t(reaction_scores)) %>%
      rownames_to_column("Pool_ID_new") %>%
      left_join(pools_final_df %>% select(Pool_ID_new, Major_Group), by = "Pool_ID_new") %>%
      # 转换为长表格，并关联子系统信息
      pivot_longer(cols = -c(Pool_ID_new, Major_Group), names_to = "reaction_id", values_to = "Score") %>%
      mutate(match_id = str_remove(reaction_id, "_pos$|_neg$")) %>%
      left_join(reaction_metadata %>% select(reaction_no_direction, subsystem), 
                by = c("match_id" = "reaction_no_direction")) %>%
      # 过滤掉没有子系统信息的反应，并设定组别顺序
      filter(!is.na(subsystem)) %>%
      mutate(Major_Group = factor(Major_Group, levels = group_order))
    
    # 选出你感兴趣或者差异最显著的子系统进行展示
    target_subsystems <- c("Citric acid cycle", "Pyruvate metabolism", 
                           "Steroid metabolism", "Purine synthesis")
    
    plot_trend_sub <- trend_data %>%
      filter(subsystem %in% target_subsystems)
    
    library(dplyr)
    
    # 1. 统计每个 Cluster 在对应通路（subsystem）下的均值
    cluster_subsystem_means <- plot_trend_sub %>%
      group_by(Pool_ID_new, Major_Group, subsystem) %>%
      summarise(Mean_Score = mean(Score, na.rm = TRUE), .groups = "drop")
    
    # 2. 查看结果
    head(cluster_subsystem_means)
    
    library(ggplot2)
    library(ggridges)
    
    library(ggplot2)
    library(ggridges)
    library(viridis)
    
    # 绘制极致美化的山峦图
    p_ridges_final <- ggplot(cluster_subsystem_means, aes(x = Mean_Score, y = Major_Group, fill = ..x..)) +
      # 1. 核心：使用渐变填充，根据 x 轴数值(得分)着色，增加视觉深度
      geom_density_ridges_gradient(
        scale = 1.6,           # 山峰重叠程度，越大重叠越多
        rel_min_height = 0.01, # 过滤掉极小的尾部
        color = "white",       # 山脊线颜色
        size = 0.4,            # 山脊线粗细
        quantile_lines = TRUE, # 绘制分位数线
        quantiles = 2          # 只画中位数线 (50%)
      ) +
      
      # 2. 配色优化：使用 Magma 或 Viridis 渐变，或者冷暖对称色
      # 这里推荐用 Option "D" (Viridis) 或 "A" (Magma)，非常有质感
      scale_fill_viridis_c(option = "D", name = "Activity Score", direction = 1) +
      
      # 3. 分面美化：增加边框和背景，让子图更独立
      facet_wrap(~ subsystem, scales = "free", ncol = 2) +
      
      # 4. 坐标轴与主题精修
      theme_ridges(font_size = 12, grid = FALSE, center_axis_labels = TRUE) +
      labs(
        title = "Evolution of Metabolic Pathway Activity",
        subtitle = "Distribution shift across gout progression stages",
        x = "COMPASS Mean Activity Score",
        y = NULL
      ) +
      
      # 5. 细节定制
      theme(
        # 分面标题样式：粗斜体 + 浅灰背景
        strip.background = element_rect(fill = "#F2F2F2", color = "grey80"),
        strip.text = element_text(face = "bold.italic", size = 10, color = "#333333"),
        
        # 标题间距与对齐
        plot.title = element_text(face = "bold", size = 16, margin = margin(b = 5)),
        plot.subtitle = element_text(size = 10, color = "grey40", margin = margin(b = 20)),
        
        # 移除图例（如果不需要特定数值参考）或者放在底部
        legend.position = "bottom",
        legend.key.width = unit(1.5, "cm"),
        legend.key.height = unit(0.3, "cm"),
        
        # 增加面板边框，使图看起来更扎实
        panel.spacing = unit(1.5, "lines"),
        axis.text.y = element_text(vjust = 0, color = "black", face = "bold"),
        axis.ticks.y = element_blank()
      )
    
    print(p_ridges_final)
    
    ggsave(plot = p_ridges_final, filename = "/data3/home/hu/Tongfeng/Result/Fig2/ridges_compass.pdf",
           width = 7,height = 5)
    
    
  }
  
}

head(meta_data[,c(1:4,228:231)])


#####Temp
{
  pseudo_bulk <- meta_data %>%
    group_by(orig.ident, group, celltype) %>%
    summarise(across(all_of(pathway_cols), mean), .groups = "drop")
  
  
  library(tidyr)
  
  pseudo_bulk <- pseudo_bulk %>%
    unite("sample_celltype", orig.ident, celltype, sep = "_")
  
  # 构建矩阵
  expr_mat <- pseudo_bulk %>%
    column_to_rownames("sample_celltype") %>%
    dplyr::select(all_of(pathway_cols)) %>%
    as.matrix()
  
  group_info <- pseudo_bulk$group
  names(group_info) <- rownames(expr_mat)
  
  table(group_info)
  group2 <- ifelse(group_info == "Control", "Control", "Disease")
  group2 <- factor(group2, levels = c("Control", "Disease"))
  
  #####
  
  
  
  
  
  
  
  
  data_clean <- data[!data$Subsystem == "NA",]
  data_clean <- data_clean %>%
    separate_rows(Subsystem, sep = " & ")
  
  gene_sets <- data_clean %>%
    group_by(Subsystem) %>%
    summarise(genes = list(unique(`Gene symbol`)))
  
  gene_sets_list <- setNames(gene_sets$genes, gene_sets$Subsystem)
  
  # 查看数据
  head(data[,c(1,15)])
  
  library(dplyr)
  library(tidyr)
  library(stringr)
  
  # 只取需要的列
  df <- data[, c(1, 15)]
  colnames(df) <- c("gene", "from")
  
  # 拆分来源
  df2 <- df %>%
    mutate(from = str_split(from, " & ")) %>%
    tidyr::unnest(cols = from)
  
  df_wide <- df2 %>%
    mutate(value = 1) %>%
    pivot_wider(names_from = from, values_from = value, values_fill = 0)
  
  
  # 加载必要的包
  library(UpSetR)
  library(ggplot2)
  library(readxl)
  library(tibble)
  
  df_mat <- as.data.frame(df_wide)
  # 2. 去掉 gene 列（关键！！）
  df_mat <- df_mat[, -1]
  
  # 3. 画图
  upset(
    df_mat,
    nintersects = 25,
    sets = colnames(df_mat),
    keep.order = TRUE,
    matrix.dot.alpha = 0.5,
    point.size = 2.2,
    line.size = 0.7,
    mb.ratio = c(0.55, 0.45),
    mainbar.y.label = "Intersection size",
    sets.x.label = "Set size",
    text.scale = c(1.8,1.4,1.5,1.4,1.4,1.5),
    order.by = "freq"
  )
  
  
  # 颜色映射
  set_colors <- corrplot::COL2('RdBu', 15)[c(4:7, 9:12)]
  names(set_colors) <- colnames(df_mat)
  set_colors <- set_colors[1:5]
  # 设置metadata
  meta_info <- data.frame(
    sets = colnames(df_mat),
    Group = colnames(df_mat)     # 用集合名作为分组
  )
  
  set_meta <- list(
    data = meta_info,
    plots = list(
      list(
        type = "matrix_rows",       # 矩阵行
        column = "Group",           # 根据哪列映射颜色
        colors = set_colors,        # 每个集合对应颜色
        alpha = 0.25
      )
    )
  )
  
  
  pdf("/data3/home/hu/Tongfeng/Result/Fig2/me_upset.pdf",width = 9,height = 6)
  upset(
    df_mat, nintersects = 25,
    sets = names(df_mat),
    keep.order = TRUE,             # 这里是True，会与表格中列的顺序相同
    set.metadata = set_meta,       # 这里控制矩阵背景颜色和点的颜色
    sets.bar.color = set_colors,   # 左侧集合条颜色
    main.bar.color = "#4E9AC6",
    matrix.dot.alpha = 0.5,
    point.size = 2.2,
    line.size = 0.7,
    mb.ratio = c(0.55, 0.45),
    mainbar.y.label = "Intersection size",
    sets.x.label = "Set size",
    text.scale = c(1.8,1.4,1.5,1.4,1.4,1.5), # 分别指定Y轴标题、Y轴刻度标签、X轴标题、X轴刻度标签、集合名称和柱体数字标签的字号大小
    order.by = "freq"
  )
  
  dev.off()
  
  
  #####代谢基因打分
  modules <- data$`Gene symbol`
  library(AUCell)
  calAUCell <- function(mode,x){
    cells_rankings <- AUCell_buildRankings(x) #x为表达矩阵
    cells_AUC <- AUCell_calcAUC(mode, cells_rankings, aucMaxRank = nrow(cells_rankings)*0.1)
    return(cells_AUC)
  }
  filterAUCell <- function(x){
    cells_assignment <-AUCell_exploreThresholds(x, plotHist = TRUE, nCores = 1, assign = TRUE) #x为cells_AUC
    return(cells_assignment)
  }
  
  cells_rankings <- AUCell_buildRankings(scRNA_harmony@assays$RNA$data)
  
  all_gene <- rownames(scRNA_harmony@assays$RNA$data)
  
  modules_list <- intersect(modules,all_gene)
  
  cells_AUC <- AUCell_calcAUC(modules_list, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.1)
  AUC <- getAUC(cells_AUC)
  
  ####提取模块打分
  AUC <- as.data.frame(t(AUC))
  sample.meta <- scRNA_harmony@meta.data
  AUC$genename <- "metabolism"
  AUC <- AUC[rownames(sample.meta),]
  
  ##确定分组
  AUC$celltype <- sample.meta$celltype
  AUC$group <- sample.meta$group
  
  
  library(tidyr)
  library(ggpubr)
  library(ggsignif)
  library(dplyr)
  
  df_plot <- AUC %>%
    group_by(group) %>%
    summarise(mean_auc = mean(geneSet, na.rm = TRUE)) %>%
    arrange(desc(mean_auc))
  
  
  library(ggplot2)
  df_plot1 <- AUC %>%
    group_by(celltype) %>%
    summarise(
      mean_auc = mean(geneSet),
      sd_auc = sd(geneSet)
    )
  
  
  df_plot2 <- AUC %>%
    group_by(group) %>%
    summarise(
      mean_auc = mean(geneSet),
      sd_auc = sd(geneSet)
    )
  
  library(ggplot2)
  
  p1 <- ggplot(df_plot2, aes(x = reorder(group, mean_auc), y = mean_auc, fill = mean_auc)) +
    geom_col(width = 0.7, color = "black", size = 0.3) +
    
    # 误差线优化
    geom_errorbar(aes(ymin = mean_auc - sd_auc, ymax = mean_auc + sd_auc),
                  width = 0.15, size = 0.6, color = "black") +
    
    # 数值标签
    geom_text(aes(label = round(mean_auc, 2)),
              hjust = -0.2, size = 4) +
    
    # 翻转坐标（更适合长标签）
    coord_flip() +
    
    # 渐变配色
    scale_fill_gradient(low = "#91D1C2", high = "#E64B35") +
    
    # 坐标范围留白
    expand_limits(y = max(df_plot2$mean_auc + df_plot2$sd_auc) * 1.1) +
    
    # 标题和标签
    labs(
      title = "Metabolic Activity Across Cell Types",
      x = NULL,
      y = "Mean AUC"
    ) +
    
    # 主题美化
    theme_classic(base_size = 14) +
    theme(
      legend.position = "none",
      axis.text.y = element_text(size = 12, face = "bold"),
      axis.text.x = element_text(size = 11),
      axis.title.y = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16)
    )
  
  ggsave("/data3/home/hu/Tongfeng/Result/Fig2/meta_gene.pdf",width = 7,height = 5)
  
  ####代谢通路分析
  # 加载R包
  library(scMetabolism)
  library(tidyverse)
  library(rsvd)
  library(Seurat)
  library(pheatmap)
  library(ComplexHeatmap)
  library(ggsci)
  # 加载数据
  sce3 <- scRNA_harmony
  sce3@meta.data$CB <- rownames(sce2@meta.data)
  
  
  sc.metabolism.SeuratV5 <- function (obj, method = "VISION", imputation = F, ncores = 2, 
                                      metabolism.type = "KEGG") {
    countexp <- GetAssayData(obj, layer = 'counts')
    countexp <- data.frame(as.matrix(countexp))
    
    signatures_KEGG_metab <- system.file("data", "KEGG_metabolism_nc.gmt", 
                                         package = "scMetabolism")
    signatures_REACTOME_metab <- system.file("data", "REACTOME_metabolism.gmt", 
                                             package = "scMetabolism")
    
    if (metabolism.type == "KEGG") {
      gmtFile <- signatures_KEGG_metab
      cat("Your choice is: KEGG\n")
    }
    if (metabolism.type == "REACTOME") {
      gmtFile <- signatures_REACTOME_metab
      cat("Your choice is: REACTOME\n")
    }
    
    if (imputation == F) {
      countexp2 <- countexp
    }
    if (imputation == T) {
      cat("Start imputation...\n")
      cat("Citation: George C. Linderman, Jun Zhao, Yuval Kluger. Zero-preserving imputation of scRNA-seq data using low-rank approximation. bioRxiv. doi: https://doi.org/10.1101/397588 \n")
      result.completed <- alra(as.matrix(countexp))
      countexp2 <- result.completed[[3]]
      row.names(countexp2) <- row.names(countexp)
    }
    
    cat("Start quantify the metabolism activity...\n")
    
    if (method == "VISION") {
      library(VISION)
      n.umi <- colSums(countexp2)
      scaled_counts <- t(t(countexp2)/n.umi) * median(n.umi)
      vis <- Vision(scaled_counts, signatures = gmtFile)
      options(mc.cores = ncores)
      vis <- analyze(vis)
      signature_exp <- data.frame(t(vis@SigScores))
    }
    if (method == "AUCell") {
      library(AUCell)
      library(GSEABase)
      cells_rankings <- AUCell_buildRankings(as.matrix(countexp2), 
                                             nCores = ncores, plotStats = F)
      geneSets <- getGmt(gmtFile)
      cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings)
      signature_exp <- data.frame(getAUC(cells_AUC))
    }
    if (method == "ssGSEA") {
      library(GSVA)
      library(GSEABase)
      geneSets <- getGmt(gmtFile)
      gsva_es <- gsva(as.matrix(countexp2), geneSets, method = c("ssgsea"), 
                      kcdf = c("Poisson"), parallel.sz = ncores)
      signature_exp <- data.frame(gsva_es)
    }
    
    cat("\nPlease Cite: \nYingcheng Wu, Qiang Gao, et al. Cancer Discovery. 2021. \nhttps://pubmed.ncbi.nlm.nih.gov/34417225/   \n\n")
    
    # 将代谢评分存储到 Seurat 对象中
    obj@assays$METABOLISM$score <- signature_exp
    obj
  }
  
  
  res <- sc.metabolism.SeuratV5(obj = sce3,
                                method = "AUCell",
                                imputation = F, 
                                ncores = 20, 
                                metabolism.type = "KEGG")
  
  DimPlot.metabolismV5 <- function(obj, pathway, dimention.reduction.type = "umap", 
                                   dimention.reduction.run = TRUE, size = 1) {
    cat("\nPlease Cite: \nYingcheng Wu, Qiang Gao, et al. Cancer Discovery. 2021. \nhttps://pubmed.ncbi.nlm.nih.gov/34417225/   \n\n")
    
    if (dimention.reduction.type == "umap") {
      if (dimention.reduction.run == TRUE) 
        obj <- Seurat::RunUMAP(obj, reduction = "pca", dims = 1:40)
      
      umap.loc <- obj@reductions$umap@cell.embeddings
      row.names(umap.loc) <- colnames(obj)
      
      signature_exp <- obj@assays$METABOLISM$score
      input.pathway <- pathway
      
      signature_ggplot <- data.frame(umap.loc, t(signature_exp[input.pathway, ]))
      
      library(wesanderson)
      pal <- wes_palette("Zissou1", 100, type = "continuous")
      library(ggplot2)
      
      plot <- ggplot(data = signature_ggplot, aes(x = umap_1, y = umap_2, 
                                                  color = signature_ggplot[, 3])) + 
        geom_point(size = size) + 
        scale_fill_gradientn(colours = pal) + 
        scale_color_gradientn(colours = pal) + 
        labs(color = input.pathway) + 
        xlab("UMAP 1") + ylab("UMAP 2") + 
        theme(aspect.ratio = 1) + 
        theme(panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(), 
              panel.background = element_blank(), 
              axis.line = element_line(colour = "black"))
    }
    
    if (dimention.reduction.type == "tsne") {
      if (dimention.reduction.run == TRUE) 
        obj <- Seurat::RunTSNE(obj, reduction = "pca", dims = 1:40)
      
      tsne.loc <- obj@reductions$tsne@cell.embeddings
      row.names(tsne.loc) <- colnames(obj)
      
      signature_exp <- obj@assays$METABOLISM$score
      input.pathway <- pathway
      
      signature_ggplot <- data.frame(tsne.loc, t(signature_exp[input.pathway, ]))
      
      pal <- wes_palette("Zissou1", 100, type = "continuous")
      library(ggplot2)
      
      plot <- ggplot(data = signature_ggplot, aes(x = tSNE_1, y = tSNE_2, 
                                                  color = signature_ggplot[, 3])) + 
        geom_point(size = size) + 
        scale_fill_gradientn(colours = pal) + 
        scale_color_gradientn(colours = pal) + 
        labs(color = input.pathway) + 
        xlab("tSNE 1") + ylab("tSNE 2") + 
        theme(aspect.ratio = 1) + 
        theme(panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(), 
              panel.background = element_blank(), 
              axis.line = element_line(colour = "black"))
    }
    
    return(plot)
  }
  
  
  
  DimPlot.metabolismV5(obj = res, 
                       pathway = "Glycolysis / Gluconeogenesis", 
                       dimention.reduction.type = "umap",  
                       dimention.reduction.run = F, size = 1)
  
  saveRDS(res,"/data3/home/hu/Tongfeng/Result/Fig2/metabolismV5_res.RDS")
  
  countexp.Seurat <- res
  input.pathway <- rownames(countexp.Seurat@assays[["METABOLISM"]][["score"]])[1:30]
  DotPlot.metabolism(obj = countexp.Seurat, 
                     pathway = input.pathway, 
                     phenotype = "celltype",  # 使用meta.data中的celltype列
                     norm = "y") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))  # 旋转x轴标签
  
  
  metabolism_data <- countexp.Seurat@assays[["METABOLISM"]][["score"]] 
  metabolism_data <- as.data.frame(t(metabolism_data))
  
  metabolism_data$celltype <- countexp.Seurat$celltype
  metabolism_data$group <- countexp.Seurat$group
  metabolism_data$sample   <- countexp.Seurat$orig.ident
  #####按照样本和细胞类型进行划分伪bulk
  metabolism_data[1:5,1:5]
  pseudo_bulk <- metabolism_data %>%
    as.data.frame() %>%
    tibble::rownames_to_column("cell") %>%
    group_by(sample, group, celltype) %>%
    summarise(across(where(is.numeric), mean), .groups = "drop")
  
  
  # 选出代谢通路列（去掉celltype和group）
  metabolic_cols <- setdiff(colnames(metabolism_data), c("celltype", "group"))
  
  # 计算每个细胞的平均代谢强度
  metabolism_data$mean_metabolism <- rowMeans(metabolism_data[, metabolic_cols], na.rm = TRUE)
  
  df_summary <- metabolism_data %>%
    group_by(celltype, group) %>%
    summarise(
      mean_metabolism = mean(mean_metabolism),
      sd_metabolism = sd(mean_metabolism),
      .groups = "drop"
    )
  
  library(ggplot2)
  
  ggplot(df_summary, aes(x = celltype, y = mean_metabolism, fill = celltype)) +
    geom_col(position = position_dodge(width = 0.7), width = 0.6) +
    geom_errorbar(aes(ymin = mean_metabolism - sd_metabolism,
                      ymax = mean_metabolism + sd_metabolism),
                  position = position_dodge(width = 0.7),
                  width = 0.2) +
    coord_flip() +
    theme_classic(base_size = 14) +
    labs(y = "Mean metabolic activity", x = NULL)
  
  
  ggplot(df_summary, aes(x = group, y = mean_metabolism, fill = group)) +
    geom_col(position = position_dodge(width = 0.7), width = 0.6) +
    geom_errorbar(aes(ymin = mean_metabolism - sd_metabolism,
                      ymax = mean_metabolism + sd_metabolism),
                  position = position_dodge(width = 0.7),
                  width = 0.2) +
    coord_flip() +
    theme_classic(base_size = 14) +
    labs(y = "Mean metabolic activity", x = NULL)
  
  
  #####寻找随着阶段逐渐升高的通路类型
  monocyte_data <- metabolism_data %>% 
    filter(celltype == "Monocyte")
  
  metabolic_cols <- setdiff(colnames(monocyte_data), c("celltype", "group"))
  
  df_long <- monocyte_data %>%
    pivot_longer(
      cols = all_of(metabolic_cols),
      names_to = "pathway",
      values_to = "score"
    )
  
  df_avg <- df_long %>%
    group_by(group, pathway) %>%
    summarise(mean_score = mean(score), .groups = "drop")
  
  # 先把 group 转数值（比如 1,2,3...）
  df_avg$group_num <- as.numeric(factor(df_avg$group, levels = sort(unique(df_avg$group))))
  
  # 计算每条通路的 Spearman 相关系数
  trend_df <- df_avg %>%
    group_by(pathway) %>%
    summarise(
      cor_trend = cor(group_num, mean_score, method = "spearman"),
      .groups = "drop"
    ) %>%
    arrange(desc(cor_trend))  # 正相关越大，越随阶段升高
  
  head(trend_df, 10)
  
  
  
  
  
}