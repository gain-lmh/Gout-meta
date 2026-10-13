#数据过滤
rm(list = ls())
options(stringsAsFactors = F)

#####按照样本进行WGCNA分析
pb_meta_filter <- readRDS("/data3/home/hu/Tongfeng/new_Result//Fig2//pb_meta_filter_new.RDS")
pb_meta_filter <- as.data.frame(pb_meta_filter)
rownames(pb_meta_filter) <- paste(pb_meta_filter$orig.ident,pb_meta_filter$celltype,sep = "_")


data <- pb_meta_filter[,-c(1:3)]
dim(data)
#3325 147

#构建基因共表达网络
library(WGCNA)
#1.判断数据质量
gsg <- goodSamplesGenes(data, verbose = 3)
gsg$allOK  #这个数据结果是TRUE

#2.设置间隔
powers <- c(c(1:10), seq(from = 12, to=30, by=2)) 
dat <- as.data.frame(data) #注意倒置；
sft <- pickSoftThreshold(dat, powerVector = powers, verbose = 5) 

#3.画βvalue图,确定阈值
pdf("/data3/home/hu/Tongfeng/new_Result//Fig3/beta-value.pdf",width = 9,height = 4)
par(mfrow = c(1,2))
cex1 <- 0.95
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="signed R^2",type="n",main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red")
abline(h=0.95,col="red")
plot(sft$fitIndices[,1], sft$fitIndices[,5],xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()

options(stringsAsFactors = FALSE)
allowWGCNAThreads()
library(WGCNA)
cor <- WGCNA::cor 

#4.构建共表达网络
#设置一下power和maxBlockSize两个数值，power就是上面得到的R^2达到0.9的数值6；
net <- blockwiseModules(dat, 
                        power = 10, 
                        maxBlockSize = 2000,
                        TOMType ='unsigned', 
                        minModuleSize = 4, 
                        corType = "pearson",   # ← 加这一句
                        
                        ## 关键调整部分 ##
                        deepSplit = 3,          # 增加切割深度，设为最大值4
                        mergeCutHeight = 0.,   # 降低合并阈值，减少模块合并
                        ################
                        
                        reassignThreshold = 0, 
                        numericLabels = TRUE, 
                        pamRespectsDendro = FALSE,
                        saveTOMs = TRUE, 
                        saveTOMFileBase = "drought",
                        verbose = 3)

table(net$colors)


# 0    1   2   3   4
#2582 190 188 184 181

#5.绘制基因聚类和模块颜色组合
moduleLabels <- net$colors
moduleColors <- labels2colors(moduleLabels)
pdf("/data3/home/hu/Tongfeng/new_Result//Fig3/genes-modules.pdf",width = 8,height = 4)
plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

{
  # 假设 net$colors 中有数字标签 0, 1, 2, 3...
  # 我们保留 0 为灰色 (grey)，其余替换为 mycol 中的颜色
  num_modules <- length(table(net$colors)) - 1 # 减去 grey 模块
  my_module_colors <- rep("#989EBE", length(net$colors)) # 默认全是灰色
  
  # 循环映射：将模块 1 设为 mycol[1]，模块 2 设为 mycol[2]，依此类推
  mycol <- c("#989EBE","#386388","#876EC1","#8AB6C3","#DA884F","#E19F8C","#D0CC48","#668BAC")
  mycol_1 <- mycol[-1]
  for(i in 1:num_modules) {
    my_module_colors[net$colors == i] <- mycol[i]
  }
  
  # 此时 my_module_colors 就是使用了你自定义颜色的向量
  # 绘图
  pdf("/data3/home/hu/Tongfeng/new_Result//Fig3/genes-modules.pdf",width = 8,height = 4)
  plotDendroAndColors(
    net$dendrograms[[1]], 
    my_module_colors[net$blockGenes[[1]]], # 使用自定义颜色向量
    "Custom Modules",                      # 这里的标题可以改
    dendroLabels = FALSE, 
    hang = 0.03,
    addGuide = TRUE, 
    guideHang = 0.05,
    main = "Cluster Dendrogram with Custom Colors" # 添加一个主标题
  )
  dev.off()
  
}

###模块Count绘制柱状图
{
  library(ggplot2)
  # 1. 统计每个模块的数量
  module_counts <- as.data.frame(table(net$colors))
  colnames(module_counts) <- c("Module_ID", "Count")
  
  # 2. 将数字 ID 转换为 ME 标签 (例如 0 变成 ME0)
  module_counts$Module <- paste0("ME", module_counts$Module_ID)
  
  # 3. 排序：确保按照 ME0, ME1, ME2... 的顺序排列
  module_counts$Module <- factor(module_counts$Module, 
                                 levels = paste0("ME", 0:(nrow(module_counts)-1)))
  
  library(ggplot2)
  library(dplyr)
  
  # 1. 整理数据并进行排序
  module_counts_sorted <- module_counts %>%
    # 按照 Count 从大到小排序
    arrange(desc(Count)) %>%
    # 重新定义因子水平，使得绘图时从上到下按数量排列
    mutate(Module = factor(Module, levels = rev(Module)))
  
  # 2. 准备颜色映射 (需要根据排序后的 Module 顺序来匹配颜色)
  # 注意：这里需要确保颜色与 Module ID 对应准确
  # 建议先建立一个全量的颜色对应表
  all_colors <- c("ME0" = "#989EBE", setNames(mycol_1[1:(nrow(module_counts)-1)], paste0("ME", 1:(nrow(module_counts)-1))))
  current_colors <- all_colors[levels(module_counts_sorted$Module)]
  
  # 3. 绘制横向柱状图
  p1 <- ggplot(module_counts_sorted, aes(x = Module, y = Count, fill = Module)) +
    # coord_flip() 将纵向改为横向
    # width 设置为 0.9 或更大，可以使柱子之间的间隙变小
    geom_bar(stat = "identity", width = 0.85, color = "white", size = 0.2) +
    # 在柱子末端添加数值
    geom_text(aes(label = Count), hjust = -0.2, size = 3.5, fontface = "bold") +
    # 反转坐标轴
    coord_flip() +
    # 应用颜色
    scale_fill_manual(values = current_colors) +
    # 细节美化
    labs(title = "Pathway Counts per Metabolic Module",
         x = NULL, 
         y = "Number of Pathways") +
    theme_classic() +
    theme(
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      axis.text.y = element_text(face = "bold", size = 10),
      axis.line.y = element_blank(), # 去掉纵轴线让图更清爽
      axis.ticks.y = element_blank()
    ) +
    # 适当放宽 y 轴上限给标签留位
    scale_y_continuous(expand = expansion(mult = c(0, 0.15)))
  
  ggsave("/data3/home/hu/Tongfeng/new_Result//Fig3/Count_Module_order.pdf", 
         p1, width = 3, height = 4)  
  
}
#1.导入临床信息
data_clin <- pb_meta_filter
data_clin <- as.data.frame(data_clin)
rownames(data_clin) <- paste(data_clin$orig.ident,data_clin$celltype,sep = "_")
data_clin <- data_clin[,1:3]

# 1. 转换分类变量为 0/1 矩阵 (One-hot encoding)
library(magrittr)
design <- model.matrix(~0 + celltype, data = data_clin)
# 简化列名，去掉前缀 "group"
colnames(design) <- unique(data_clin$celltype)

# 2. 计算相关性和 P 值
MEs <- net$MEs
# 确保样本顺序一致
MEs <- MEs[rownames(design), ]
moduleTraitCor <- cor(MEs, design, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples = nrow(design))

library(WGCNA)

# 设置热图显示的文本（相关系数及P值）
textMatrix <- paste(round(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) <- dim(moduleTraitCor)

# 绘图
library(ComplexHeatmap)
library(circlize) # 用于精确控制颜色映射

# 1. 定义更高级的配色（红蓝经典，但色调更柔和）
col_fun = colorRamp2(c(-1, 0, 1), c("#2166AC", "white", "#B2182B"))

# 2. 绘制热图
Heatmap(moduleTraitCor, 
        name = "Correlation", 
        col = col_fun,
        
        # 单元格内的文字（显著性/相关系数）
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(textMatrix[i, j], x, y, gp = gpar(fontsize = 8))
        },
        
        # 样式调整
        cluster_rows = TRUE,         # 开启行聚类，让模块关系更清晰
        cluster_columns = FALSE,      # 通常性状不聚类，保持设计矩阵顺序
        row_names_side = "left",
        column_names_rot = 45,        # 标签倾斜 45 度，防止重叠
        
        # 边框与间隙
        rect_gp = gpar(col = "white", lwd = 1), # 单元格白边，增加呼吸感
        
        # 标题设置
        column_title = "Module-trait relationships",
        column_title_gp = gpar(fontsize = 14, fontface = "bold"),
        
        # 图例美化
        heatmap_legend_param = list(
          title = "Corr",
          at = c(-1, -0.5, 0, 0.5, 1),
          labels = c("-1", "-0.5", "0", "0.5", "1")
        ))



{
  # 1. 转换分类变量为 0/1 矩阵 (One-hot encoding)
  library(magrittr)
  design <- model.matrix(~0 + group, data = data_clin)
  # 简化列名，去掉前缀 "group"
  colnames(design) <- levels(data_clin$group)
  
  # 2. 计算相关性和 P 值
  MEs <- net$MEs
  # 确保样本顺序一致
  MEs <- MEs[rownames(design), ]
  moduleTraitCor <- cor(MEs, design, use = "p")
  moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples = nrow(design))
  
  library(WGCNA)
  
  # 设置热图显示的文本（相关系数及P值）
  textMatrix <- paste(round(moduleTraitCor, 2), "\n(",
                      signif(moduleTraitPvalue, 1), ")", sep = "")
  dim(textMatrix) <- dim(moduleTraitCor)
  
  # 绘图
  labeledHeatmap(Matrix = moduleTraitCor,
                 xLabels = colnames(design),
                 yLabels = names(MEs),
                 ySymbols = names(MEs),
                 colorLabels = FALSE,
                 colors = blueWhiteRed(50),
                 textMatrix = textMatrix,
                 setStdMargins = FALSE,
                 cex.text = 0.5,
                 zlim = c(-1,1),
                 main = "Module-trait relationships")
  
}

# 绘制模块中包含的通路
# 注意：power 必须和你 blockwiseModules 里的设置保持一致，即 10
softPower = 10
adjacency = adjacency(dat, power = softPower, type = "unsigned")
TOM = TOMsimilarity(adjacency, TOMType = "unsigned")

circle_list <- list()

for (i in 0:6) {
  
  mycol <- c("#989EBE","#386388","#876EC1","#8AB6C3","#DA884F","#E19F8C","#D0CC48","#668BAC")
  # 2. 确定模块颜色和通路名称
  # 在你的代码中，ME3 对应的数字标签是 3，我们找到它对应的颜色名
  module_label = i 
  module_color = labels2colors(module_label) # 得到的应该是绿色或其他颜色
  pathway_names = colnames(dat)
  
  # 3. 提取模块内的索引
  inModule = (moduleLabels == module_label)
  modTOM = TOM[inModule, inModule]
  modPathways = pathway_names[inModule]
  
  # 给矩阵命名，这样导出的网络图节点才有名字
  dimnames(modTOM) = list(modPathways, modPathways)
  
 library(igraph)
 library(ggraph)
 library(tidygraph)
  
  pathway_order <- read.csv("/data3/home/hu/Tongfeng/Result/Fig3/pathway_order.csv")
  # 1. 确保 modTOM 的行名和列名是通路名称
  # modPathways 是你之前提取的通路名称向量
  dimnames(modTOM) = list(modPathways, modPathways)
  rownames(modTOM) <- pathway_order$abbreviation[match(rownames(modTOM),pathway_order$pathway)]
  colnames(modTOM) <- pathway_order$abbreviation[match(colnames(modTOM),pathway_order$pathway)]
  
  # 2. 转换为 graph 对象并强制转换节点属性
  graph <- as_tbl_graph(modTOM, directed = FALSE) %>%
    activate(nodes) %>%
    mutate(name = rownames(modTOM) ) %>% # 强制赋予 name 属性
    mutate(degree = centrality_degree()) %>%
    activate(edges) 
  
  # 3. 绘图
  p1 <- ggraph(graph, layout = 'circle') + 
    geom_edge_link(aes(edge_alpha = weight), edge_colour = "#A1A1A1", show.legend = FALSE) + 
    geom_node_point(aes(size = degree), color = mycol[i+1], show.legend = FALSE) + 
    # 使用 vjust 和 hjust 让文字稍微偏离圆点，或者设置 repel = TRUE
    geom_node_text(aes(label = name), size = 3, color = "black", repel = TRUE) + 
    theme_graph() +
    theme(plot.margin = margin(20, 20, 20, 20))
  
  circle_list[[i+1]] <- p1
}
# install.packages("patchwork")
library(patchwork)
# 使用 wrap_plots 函数自动排版
combined_plot <- wrap_plots(circle_list, ncol = 4, nrow = 2) 

# 显示图像
combined_plot

# 保存高分辨率图片
circle_list[[3]]
ggsave("/data3/home/hu/Tongfeng/new_Result//Fig3/combined_circle_C_MP.pdf", 
       combined_plot, width = 10, height = 6)


#####计算每个代谢模块中这些通路的贡献
# 计算所有通路与所有模块的相关性
geneModuleMembership = as.data.frame(cor(dat, MEs, use = "p"))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples = nrow(dat)))

# 整理每个模块的贡献排名
library(dplyr)
contribution_list <- list()

for(i in 1:ncol(MEs)){
  mod_name <- colnames(MEs)[i]
  # 提取属于该模块的通路
  # 注意：net$colors 的标签可能是数字 0,1,2... 需要根据你的 MEs 顺序对应
  # 假设 MEs 的列名是 ME1, ME2... 对应 net$colors 的 1, 2...
  mod_idx <- as.numeric(gsub("ME", "", mod_name)) 
  pathways_in_mod <- colnames(dat)[net$colors == mod_idx]
  
  # 提取 MM 值并排序
  df <- data.frame(
    Pathway = pathways_in_mod,
    Contribution_Score = geneModuleMembership[pathways_in_mod, mod_name]
  ) %>% arrange(desc(Contribution_Score))
  
  contribution_list[[mod_name]] <- df
}

# 查看 ME3 模块中通路贡献排名
head(contribution_list[["ME3"]])
#####可视化
{
  library(ggplot2)
  library(ggrepel)
  library(dplyr)
  
  # 1. 将 list 转换为 dataframe
  plot_df <- bind_rows(contribution_list, .id = "Module")
  
  # 2. 标记出每个模块贡献度最高的通路 (Top 1)
  plot_df <- plot_df %>%
    group_by(Module) %>%
    mutate(is_top = ifelse(Contribution_Score == max(Contribution_Score), "yes", "no")) %>%
    ungroup()
  
  # 3. 设置模块顺序 (ME1, ME2...)
  plot_df$Module <- factor(plot_df$Module, levels = paste0("ME", 0:(length(unique(plot_df$Module))-1)))
  plot_df_abs <- plot_df
  plot_df_abs$Contribution_Score <- abs(plot_df_abs$Contribution_Score)
  
  library(ggplot2)
  library(ggrepel)
  library(dplyr)
  
  top_moudel <- plot_df_abs[plot_df_abs$is_top == "yes" ,]
  top_moudel$abbreviation <-  pathway_order$abbreviation[match(top_moudel$Pathway,pathway_order$pathway)]
  plot_df_abs$abbreviation <-  pathway_order$abbreviation[match(plot_df_abs$Pathway,pathway_order$pathway)]
  write.csv(top_moudel,"/data3/home/hu/Tongfeng/new_Result//Fig3/top_moudel_abbreviation.csv")
  write.csv(plot_df_abs,"/data3/home/hu/Tongfeng/new_Result//Fig3/plot_df_abs_abbreviation.csv")
  
  
  # 1. 绘图：箱线图 + 抖动散点
  p1 <- ggplot(plot_df_abs, aes(x = Module, y = Contribution_Score, color = Module)) +
    # 底层：绘制箱线图，去掉填充色，只留线条
    geom_boxplot(aes(fill = Module), alpha = 0.1, outlier.shape = NA, width = 0.5) +
    # 中层：绘制所有通路的抖动散点
    geom_jitter(width = 0.2, alpha = 0.8, size = 1) +
    # 顶层：突出显示 Top 1 的点（更大的实心圆或菱形）
    geom_point(data = filter(plot_df, is_top == "yes"), 
               size = 2, shape = 16) + 
    # 颜色应用
    scale_color_manual(values = mycol) +
    scale_fill_manual(values = mycol) +
    # 调整坐标轴范围：不强制从0开始，根据数据自动调整
    # 如果你想手动固定范围，可以使用 coord_cartesian(ylim = c(0.4, 1))
    
    # 细节修饰
    labs(title = "Contribution of Pathways in each Module",
         x = NULL, 
         y = "Module Membership (Correlation)") +
    theme_bw() +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10, face = "bold"),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(), # 关闭垂直网格线，更像参考图
      plot.title = element_text(hjust = 0.5, face = "bold", size = 12)
    )
  
  # 保存高分辨率图片
  ggsave("/data3/home/hu/Tongfeng/new_Result//Fig3/Contribution_Module.pdf", 
         p1, width = 6, height = 4)  
}

#####计算WGCNA模块与Mfuzz模块交集
{
  mfuzz_clusters <- readRDS("/data3/home/hu/Tongfeng/new_Result/Fig3/mfuzz_clusters.RDS")
  mfuzz_pathway <- mfuzz_clusters[mfuzz_clusters$Cluster == 6,]
  
  moduleLabels_data <- data.frame(pathway = names(moduleLabels),MP = moduleLabels)
  moduleLabels_data <- moduleLabels_data[order(moduleLabels_data$MP),]
  wgcna_pathway <- moduleLabels_data[moduleLabels_data$MP == 2,]
  
  inter_module <- intersect(wgcna_pathway$pathway,mfuzz_pathway$Pathway)
  inter_module <- inter_module[1:5]
}

#####计算模块在每个细胞得分
meta_data <- readRDS("/data3/home/hu/Tongfeng/new_Result/Fig2/pathway_meta_data.rds")

{
  # 找出属于该模块的通路名称
  pathways_in_mod <- inter_module
  
  # 提取这些通路在各样本中的表达值 (Score)
  sub_dat <- meta_data[, c("orig.ident","group","celltype",pathways_in_mod), drop = FALSE]
  
  # 提取这些通路的贡献值 (MM值)
  # 确保从之前算的 geneModuleMembership 中提取对应的 MM 列
  weights <- geneModuleMembership[pathways_in_mod, "ME2"]
  
  # 计算：(通路表达值 * 对应权重) 之和 / 模块通路数量
  # 使用矩阵乘法提高效率
  weighted_sum <- as.matrix(sub_dat[,-(1:3)]) %*% as.matrix(weights)
  weighted_module_scores <- weighted_sum / length(pathways_in_mod)
  
  sub_dat$WGCNA <- weighted_module_scores[,1]
}

####计算髓系细胞中WGCNA得分与免疫炎症的相关性
{
  Hallmarker_score <- readRDS("/data3/home/hu/Tongfeng/new_Result/Fig2/df_score_hallmarker.RDS")
  Monocyte_score <- Hallmarker_score[Hallmarker_score$cell_type2 == "Monocyte", ]
  sub_dat_Monocyte <- sub_dat[sub_dat$celltype == "Monocyte", ]
  
  library(ggplot2)
  library(ggpubr)
  library(dplyr)
  
  # 定义封装函数
  # score_matrix: 包含通路得分的矩阵 (如 Monocyte_score)
  # dat_matrix: 包含对比得分的矩阵 (如 sub_dat_Monocyte)
  # path_x: x轴的通路名
  # path_y: y轴的通路名
  plot_metacell_cor <- function(score_matrix, dat_matrix, path_x, path_y, cells_per_spot = 10) {
    
    # 1. 提取并匹配数据
    common_cells <- intersect(rownames(score_matrix), rownames(dat_matrix))
    
    df <- data.frame(
      hallmark = score_matrix[common_cells, path_x],
      metabolism = dat_matrix[common_cells, path_y]
    ) %>% na.omit()
    
    # 2. 随机聚合为 Metacell (Spot)
    set.seed(123)
    df_spot <- df %>%
      mutate(spot_id = sample(rep(1:ceiling(n()/cells_per_spot), 
                                  each = cells_per_spot, length.out = n()))) %>%
      group_by(spot_id) %>%
      summarise(
        hallmark = mean(hallmark),
        metabolism = mean(metabolism)
      )
    
    # 3. 绘图
    p <- ggplot(df_spot, aes(x = hallmark, y = metabolism)) +
      # 回归线
      geom_smooth(method = "lm", color = "#555555", fill = "#EEEEEE", 
                  se = TRUE, linewidth = 0.8, alpha = 0.5) +
      
      # 散点 (莫兰迪蓝)
      geom_point(size = 1.5, shape = 21, fill = "#8491B4FF", 
                 color = "#333333", stroke = 0.5, alpha = 0.8) +
      
      # 统计标注
      stat_cor(method = "spearman", 
               label.x.npc = "left", 
               label.y.npc = "top",
               size = 4.5, 
               fontface = "italic") +
      
      # 主题美化
      theme_bw(base_size = 14) +
      theme(
        aspect.ratio = 1,
        panel.grid.major = element_line(color = "grey95"),
        panel.grid.minor = element_blank(),
        axis.title = element_text(face = "bold", size = 12),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14)
      ) +
      
      # 动态标签
      labs(
        title = "Metacell Correlation",
        x = path_x,
        y = path_y
      )
    
    return(p)
  }
  
  # --- 使用方法演示 ---
  
  # 比如你想画当前的通路：
  p1 <- plot_metacell_cor(
    score_matrix = Monocyte_score, 
    dat_matrix = sub_dat_Monocyte, 
    path_x = "HALLMARK_INFLAMMATORY_RESPONSE", 
    path_y = "WGCNA"
  )
  print(p1)
  
  # 如果想换一个通路，直接改名字即可：
  p2 <- plot_metacell_cor(
    score_matrix = Monocyte_score, 
    dat_matrix = sub_dat_Monocyte, 
    path_x = "HALLMARK_TNFA_SIGNALING_VIA_NFKB", 
    path_y = "WGCNA"
  )
  print(p2)
  
  library(patchwork)
  
  # 1. 左右合并
  combined_p <- p1 + p2
  
  # 2. 如果想上下合并，可以使用 / 符号
  # combined_p <- p1 / p2
  
  # 3. 细节美化：收集相同图例，并添加子图标签 (A, B)
  combined_p <- (p1 + p2) + 
    plot_layout(guides = 'collect') + 
    plot_annotation(tag_levels = 'A') # 给图片左上角打上 A, B 标签
  
  print(combined_p)
  
  # 4. 保存为高分辨率图片
  ggsave("/data3/home/hu/Tongfeng/new_Result/Fig3/WGCNA_Correlation_Plots.pdf", combined_p, width = 10, height = 5)
  
  
}

####查看模块得分对于疾病状态打分
weighted_module_scores_df <- weighted_module_scores
library(dplyr)

# 按照样本和细胞类型进行合并，计算平均分
sub_dat_pb <- sub_dat %>%
  group_by(orig.ident, group, celltype) %>%
  summarise(
    # 计算所有数值列（通路得分）的平均值
    across(where(is.numeric), \(x) mean(x, na.rm = TRUE)),
    .groups = "drop" # 合并后解除分组状态
  )

# 查看合并后的结果
head(sub_dat_pb)
sub_dat_pb <- as.data.frame(sub_dat_pb)
rownames(sub_dat_pb) <- paste(sub_dat_pb$orig.ident,sub_dat_pb$celltype,sep = "_")
#sub_dat_pb <- sub_dat_pb[sub_dat_pb$group != "Remission",]

{
  library(pROC)
  library(ggplot2)
  
  # 1. 提取样本分组信息 (假设你的 data_clin 中有 group 列)
  # 确保顺序与 weighted_module_scores_df 一致
  sample_info <- data_clin[rownames(sub_dat_pb), ]
  sample_info <- sample_info[sample_info$celltype %in% c("Monocyte","cDC"),]
  
  # 2. 定义二分类标签：1 为 Disease (痛风相关), 0 为 Control (正常)
  # 请根据你实际的分组列名修改，这里假设列名是 'group'，正常组叫 'Normal'
  y_true <- ifelse(sample_info$group == "Control", 0, 1) 
  
  # 存储 ROC 结果
  roc_list <- list()
  auc_values <- c()
  
  weighted_module_scores_df_filter <- sub_dat_pb[rownames(sample_info),]
  #weighted_module_scores_df_filter <- weighted_module_scores_df
  # 循环计算每个 ME 的 ROC
  for(me in colnames(weighted_module_scores_df_filter)[4:9]) {
    res_roc <- roc(y_true, weighted_module_scores_df_filter[,me], quiet = TRUE)
    roc_list[[me]] <- res_roc
    auc_values[me] <- auc(res_roc)
  }
  
  # 打印 AUC 排名
  print(sort(auc_values, decreasing = TRUE))
  
  # 将颜色映射到模块名
  # ME0 灰色，其余用 mycol
  all_me_cols <- mycol
  names(all_me_cols) <- names(auc_values)
  
  library(pROC)
  library(ggplot2)
  library(dplyr)
  
  # 重新计算光滑后的 ROC，并生成包含 AUC 的图例名称
  smooth_roc_list <- list()
  legend_labels <- c()
  
  for(me in names(roc_list)) {
    # 开启光滑处理 (smooth = TRUE)
    s_roc <- smooth(roc_list[[me]], method = "binormal") 
    smooth_roc_list[[me]] <- s_roc
    
    # 创建带 AUC 值的图例标签，例如 "ME3 (AUC: 0.649)"
    legend_labels[me] <- paste0(me, " (AUC: ", round(auc(roc_list[[me]]), 3), ")")
  }
  
  # 确保颜色映射与新的图例标签对应
  names(all_me_cols) <- names(legend_labels)
  
  ggroc(smooth_roc_list, size = 1.2, legacy.axes = TRUE) +
    # 绘制对角基准线
    geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "grey50", size = 0.8) +
    # 应用自定义颜色并更新图例标签
    scale_color_manual(values = all_me_cols, labels = legend_labels) +
    # 坐标轴美化：强制比例 1:1，增加刻度密度
    coord_equal() +
    scale_x_continuous(expand = c(0.01, 0.01), breaks = seq(0, 1, 0.2)) +
    scale_y_continuous(expand = c(0.01, 0.01), breaks = seq(0, 1, 0.2)) +
    # 标题与标签
    labs(title = "ROC Analysis of Metabolic Modules",
         subtitle = "Weighted Scores Differentiating Gout vs. Control",
         x = "False Positive Rate (1 - Specificity)",
         y = "True Positive Rate (Sensitivity)",
         color = "Metabolic Modules") +
    # 主题深度美化
    theme_bw() + 
    theme(
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold", size = 16, hjust = 0),
      plot.subtitle = element_text(size = 10, color = "grey30"),
      legend.title = element_text(face = "bold", size = 10),
      legend.text = element_text(size = 9),
      legend.background = element_rect(fill = alpha("white", 0.5)),
      legend.position = c(0.75, 0.3), # 将图例移入图中空白处，更美观
      axis.title = element_text(face = "bold"),
      axis.text = element_text(color = "black")
    ) +
    # 突出显示 AUC 最高的模块
    annotate("label", x = 0.15, y = 0.85, 
             label = paste0("Leading Module: ", names(which.max(auc_values))),
             fill = all_me_cols[names(which.max(auc_values))],
             color = "white", fontface = "bold", size = 4)
  
}

######不同模块在级别中的打分情况
{
  library(ggplot2)
  library(tidyr)
  library(dplyr)
  library(ggpubr) # 用于添加统计显著性
  
  # 1. 合并评分与临床分组
  # 确保 rownames 匹配
  plot_data <- sub_dat_pb %>%
    mutate(SampleID = rownames(.)) %>%
    left_join(data_clin %>% mutate(SampleID = rownames(data_clin)), by = "SampleID")
  
  # 2. 转换为长格式 (针对 8 个模块)
  plot_long <- plot_data %>%
    pivot_longer(cols = starts_with("WGCNA"), 
                 names_to = "Module", 
                 values_to = "Score")
  
  data <- sub_dat_Monocyte[sub_dat_Monocyte$celltype == "Monocyte",]
  data$group <- factor(
    data$group,
    levels = c("Control", "Disease", "Remission")
  )
  
  plot_df <- data %>%
    dplyr::select(group, score = all_of("WGCNA")) %>%
    na.omit()
  
  # 计算 y 轴的最大值，用于动态调整顶部空间
  y_max <- max(plot_df$score, na.rm = TRUE)
  
  my_comparisons <- list(
    c("Control", "Disease"), 
    c("Disease", "Remission")
  )
  
  p1 <- ggplot(plot_df, aes(x = group, y = score, fill = group)) +
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
      title = paste0("WGCNA"),
      y = "Pathway Score"
    )
  
  ggsave(filename = "/data3/home/hu/Tongfeng/new_Result/Fig3/WGCNA_group.pdf",p1,
         width = 5,height = 3)
}

#####绘制模块在不同细胞类型总体趋势变化
{
  data <- read.csv("/data3/home/hu/Tongfeng/Result/all.source_Metabolism.250902.csv")
  module_list <- as.list(inter_module)
  # 创建一个空列表存储结果
  module_gene_list <- list()
  
  names(module_list) <- inter_module
  
  # 遍历每个模块
  for (mod in names(module_list)) {
    # 1. 获取该模块包含的所有通路名称
    pathways_in_mod <- module_list[[mod]]
    
    # 2. 在你的 data 表中筛选这些通路对应的基因
    # 确保 data$Metabolism 里的名称和 module_list 里的通路名完全一致
    genes <- data %>%
      filter(Metabolism %in% pathways_in_mod) %>%
      pull(Gene) %>%
      unique() # 去除在不同通路中重复出现的基因
    
    # 3. 存储结果
    module_gene_list[[mod]] <- genes
  }
  
  # 查看 ME1 模块包含的基因数量
  module_gene <- unique(unlist(module_gene_list))
  write.csv(module_gene,"/data3/home/hu/Tongfeng/new_Result/Fig3/module_gene.csv")
  
  
  #####绘制基因共表达网络
  gene_inter <- read.table("/data3/home/hu/Tongfeng/new_Result/Fig3/string_interactions_short .tsv",
                           header = T)
  
  gene_inter <- gene_inter[,c(1,2,13)]
  gene_inter <- gene_inter[gene_inter$node1 %in% module_gene & gene_inter$node2 %in% module_gene,]
  
  gene_inter_filter <- gene_inter[gene_inter$combined_score > 0.8,]
  
  {
    library(tidyverse)
    
    # 1. 将 list 转换为宽表 (One-hot encoding)
    # 这一步会为每个通路创建一列
    node_attr_wide <- stack(module_gene_list) %>%
      rename(gene = values, pathway = ind) %>%
      mutate(present = 1) %>%
      pivot_wider(names_from = pathway, 
                  values_from = present, 
                  values_fill = 0)
    
    # 2. 补全网络中存在但不在列表里的基因（设为 0）
    all_network_genes <- unique(c(gene_inter_filter$node1, gene_inter_filter$node2))
    node_attr_final <- data.frame(gene = all_network_genes) %>%
      left_join(node_attr_wide, by = "gene") %>%
      replace(is.na(.), 0)
    
    node_attr_final
    
    # 3. 导出给 Cytoscape
    write.csv(node_attr_final, "/data3/home/hu/Tongfeng/new_Result/Fig3/Cytoscape_Node_Charts.csv", 
              row.names = FALSE, quote = FALSE)
    
    write.csv(gene_inter_filter, "/data3/home/hu/Tongfeng/new_Result/Fig3/Cytoscape_Edge_Table.csv", 
              row.names = FALSE, quote = FALSE)
    
    library(tidyverse)
    
    # 假设你的 5 个通路列名如下（请根据实际情况微调列名顺序）
    pathway_cols <- c("Urea cycle", 
                      "Arachidonic acid metabolism", 
                      "Eicosanoid metabolism", 
                      "Nitrogen metabolism", 
                      "Fatty acid biosynthesis")
    
    # 指定对应的 ID 值
    pathway_ids <- c(1, 3, 5, 7, 9)
    
    node_attr_final <- node_attr_final %>%
      mutate(
        # 计算每个基因参与的通路总数
        pathway_sum = rowSums(select(., all_of(pathway_cols))),
        
        # 核心逻辑判断
        Group_ID = case_when(
          pathway_sum >= 2 ~ 11,                   # 两个及以上通路设为 11
          `Urea cycle` == 1 ~ 1,
          `Arachidonic acid metabolism` == 1 ~ 3,
          `Eicosanoid metabolism` == 1 ~ 5,
          `Nitrogen metabolism` == 1 ~ 7,
          `Fatty acid biosynthesis` == 1 ~ 9,
          TRUE ~ 0                                # 不在列表中的基因设为 0
        )
      ) %>%
      select(-pathway_sum) # 删除中间计算列
    
    # 重新导出（不带引号，确保导入顺畅）
    write.csv(node_attr_final, "Cytoscape_Node_Charts_with_ID.csv", row.names = FALSE, quote = FALSE)
    
    # 检查一下前几行结果
    head(node_attr_final)
    
    
  }
  
}
#####绘制模块基因在不同细胞类型总体趋势变化
  {
   
    ###对单细胞进行伪bulk处理
    library(qs)
    library(ggplot2)
    library(DESeq2)
    library(Seurat)
    library(tidyverse)
    library(dplyr)
    library(qs)
  
    scRNA_harmony <- readRDS("/data3/home/hu/Tongfeng/new_Result//Fig1/Fig1_scRNA_harmony.RDS")
    scRNA_Monocyte <- subset(scRNA_harmony,subset = cell_type2 == "Monocyte")
    qsave(scRNA_Monocyte,"/data3/home/hu/Tongfeng/new_Result/Fig3/scRNA_Monocyte.qs")
    
    library(Seurat)
    library(tidyverse)
    library(ComplexHeatmap)
    
    # 1. 提取 180 个基因的平均表达量
    # 确保你的组别列名是 group，级别顺序是 Control, Disease, Remission
    scRNA_Monocyte$group <- factor(scRNA_Monocyte$group, levels = c("Control", "Disease", "Remission"))
    
    # 计算每个组的平均表达值
    avg_exp <- AverageExpression(scRNA_Monocyte, 
                                 features = module_gene, # 替换为你那180个基因列表
                                 group.by = "group", 
                                 slot = "data")$RNA
    
    # 2. Z-score 标准化（为了展示趋势，必须标准化）
    exp_scaled <- t(scale(t(avg_exp)))
    exp_scaled <- na.omit(exp_scaled)
  
    # 将矩阵转为长表格
    plot_df <- as.data.frame(exp_scaled) %>%
      rownames_to_column("Gene") %>%
      pivot_longer(-Gene, names_to = "Group", values_to = "Expression")
    
    plot_df$Group <- factor(plot_df$Group, levels = c("Control", "Disease", "Remission"))
    
    # 绘图
    ggplot(plot_df, aes(x = Group, y = Expression, group = Gene)) +
      geom_line(alpha = 0.1, color = "grey") + # 绘制 180 条细线
      geom_smooth(aes(group = 1), method = "loess", color = "#E64B35", size = 1.5) + # 绘制整体拟合趋势
      theme_classic() +
      labs(title = "Overall Expression Trend of 180 Genes",
           y = "Relative Expression (Z-score)", x = "")
    
  
}
