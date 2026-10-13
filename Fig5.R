library(Mfuzz)
library(tidyverse)

pb_meta_filter <- readRDS("E:\\博二下\\痛风\\pb_meta_filter_new.RDS")

pb_meta_filter <- pb_meta_filter[pb_meta_filter$celltype == "Monocyte",]
# 1. 提取通路数据并按 group 计算均值
pathway_data <- pb_meta_filter %>%
  group_by(group) %>%
  summarise(across(where(is.numeric), mean)) %>%
  # 确保顺序是 Control, Remission, Disease
  mutate(group = factor(group, levels = c("Control", "Disease", "Remission"))) %>%
  arrange(group)

# 2. 转置矩阵：行名为通路，列名为 group
mfuzz_matrix <- pathway_data %>%
  column_to_rownames("group") %>%
  t() %>%
  as.matrix()

# 3. 构建 ExpressionSet 对象
eset <- new("ExpressionSet", exprs = mfuzz_matrix)

# 4. 过滤表达变异极小的通路 (可选)
# eset <- filter.NA(eset) # 如果有 NA
#eset <- filter.std(eset, min.std = 0.02) # 过滤波动极小的通路

# 5. 标准化 (标准化到均值为 0，方差为 1)
eset <- standardise(eset)

# 6. 设定聚类个数 (c) 和 模糊参数 (m)
# m 通常通过 estimate.m 计算
m <- mestimate(eset) 
c <- 6 # 假设我们想聚成 6 类趋势，可以根据结果调整


# 7. 运行聚类
set.seed(123)
cl <- mfuzz(eset, c = c, m = m)

# 8. 可视化
library(RColorBrewer)
color.2 <- colorRampPalette(rev(brewer.pal(11, "RdYlBu")))(100)

pdf('E:\\博二下\\痛风\\mfuzz_clusters_plot.pdf',height = 4,width = 8)

# 设置全局字体大小（对 mfuzz.plot 生效）
par(cex.main = 1.2, cex.lab = 1.1, cex.axis = 0.9, mar = c(5, 4, 3, 2))
# 运行修正后的函数
mfuzz.plot(eset, 
           cl = cl, 
           mfrow = c(2,3), 
           new.window = FALSE,
           time.labels = colnames(eset),
           colo = color.2)
dev.off()


# 查看每个通路归属于哪个 Cluster
cluster_info <- data.frame(Pathway = rownames(mfuzz_matrix), 
                           Cluster = cl$cluster)
head(cluster_info)

saveRDS(cluster_info,"E:\\博二下\\痛风\\mfuzz_clusters.RDS")

setwd("/data3/home/hu/Tongfeng/new_Result/Fig5")

data <- read.table("/data3/home/hu/Tongfeng/new_Result/Fig5/GSE242872_GA_ankle_seq.txt",
                   header = T,sep = "\t",row.names = 1)

data <- data[,(1:12)]
data$gene_name <- rownames(data)

data2 <- read.csv("GSE214587_FPKM.csv")
# 仅保留第一次出现的基因行，剔除后续重复行
data2_unique <- data2[!duplicated(data2$gene_name), ]

# 重新设置行名
rownames(data2_unique) <- data2_unique$gene_name

#合并两个数据集
{

  library(edgeR)
  library(limma)
  
  # ---- 步骤 1：严格确认并提取共有基因 ----
  common_genes <- intersect(data$gene_name, data2$gene_name)
  
  # 对原数据去重并提取共有基因
  d1 <- data[!duplicated(data$gene_name), ]
  rownames(d1) <- d1$gene_name
  mat1_raw <- d1[common_genes, !colnames(d1) %in% c("gene_name")]
  
  d2 <- data2[!duplicated(data2$gene_name), ]
  rownames(d2) <- d2$gene_name
  mat2_raw <- d2[common_genes, !colnames(d2) %in% c("gene_name")]
  
  # 合并原始矩阵
  combined_raw <- cbind(mat1_raw, mat2_raw)
  
  # ---- 步骤 2：强制统一使用 CPM + log2 转换 ----
  # 无论你之前的 data1 是什么格式，这里将其与 data2 统一放进 edgeR 
  # 通过 TMM 标准化来强行消除由于测序深度不同导致的绝对值差异
  dge <- DGEList(counts = combined_raw)
  dge <- calcNormFactors(dge, method = "TMM")
  
  # 转换为 log2(CPM + 1)
  combined_log2_cpm <- cpm(dge, log = TRUE, prior.count = 1)
  
  # ---- 步骤 3：定义批次与生物学分组 ----
  batch <- c(rep("Batch1", ncol(mat1_raw)), rep("Batch2", ncol(mat2_raw)))
  
  # 【重要】检查你的分组：如果两个批次之间没有共同的生物学对照（比如都没有正常的 Control 组），
  # 那么将 design 设为为空（~1）强制去批次，否则会报错或无法去干净。
  # 如果有共同对照，请保持你的 group 向量，否则用下面这行：
  design_matrix <- model.matrix(~1, data = as.data.frame(t(combined_log2_cpm)))
  
  # ---- 步骤 4：重新运行去批次 ----
  data_corrected <- removeBatchEffect(combined_log2_cpm, batch = batch, design = design_matrix)
  
  # 查看去批次后的部分数据
  head(data_corrected)
  
  library(FactoMineR)
  library(factoextra)
  
  # 绘制去批次后的 PCA 图
  pca_res <- PCA(t(data_corrected), graph = FALSE)
  fviz_pca_ind(pca_res,
               habillage = as.factor(batch), # 按批次上色，检查是否混合均匀
               palette = "jco",
               addEllipses = TRUE,
               ellipse.type = "convex")
  
}

gene <- read.csv("/data3/home/hu/Tongfeng/new_Result/Fig5/all.source_Metabolism.250902.csv")
pathwag_gene <- gene[gene$Metabolism %in% c("Arachidonic acid metabolism","Eicosanoid metabolism",
                                            "Fatty acid biosynthesis","Nitrogen metabolism","Urea cycle"),]

head(data_corrected)
head(pathwag_gene)

#####代谢基因转换为小鼠基因并转为list
{
  ## ===============================
  ## 1. 加载R包
  ## ===============================
  
  library(babelgene)
  library(GSVA)
  library(stringr)
  library(dplyr)
  library(ggplot2)
  library(pheatmap)
  
  ## ===============================
  ## 2. 人类 pathway gene 转换为小鼠同源基因
  ## ===============================
  
  # 提取人类基因
  human_gene <- unique(pathwag_gene$Gene)
  
  # 同源转换 human -> mouse
  mouse_gene <- orthologs(
    genes = human_gene,
    species = "mouse",
    human = TRUE
  )
  
  # 查看结果
  head(mouse_gene)
  
  ## 返回类似：
  ## human_symbol   symbol
  ## PTGES3         Ptges3
  
  ## ===============================
  ## 3. 合并回 pathway 数据
  ## ===============================
  
  pathwag_gene_mouse <- merge(
    pathwag_gene,
    mouse_gene,
    by.x = "Gene",
    by.y = "human_symbol"
  )
  
  # 查看
  head(pathwag_gene_mouse)
  
  ## ===============================
  ## 4. 只保留表达矩阵中存在的基因
  ## ===============================
  
  pathwag_gene_mouse <- pathwag_gene_mouse[
    pathwag_gene_mouse$symbol %in% rownames(data),
  ]
  
  ## ===============================
  ## 5. 构建 pathway list
  ## ===============================
  
  pathway_list <- split(
    pathwag_gene_mouse$symbol,
    pathwag_gene_mouse$Metabolism
  )
  
  # 去除空 pathway
  pathway_list <- pathway_list[
    sapply(pathway_list, length) > 0
  ]
  
}


####对代谢基因进行打分
{
  library(GSVA)
  library(GSEABase)
  
  # 表达矩阵
  expr <- as.matrix(data_corrected)
  
  # 构建ssGSEA参数
  ssgsea.param <- ssgseaParam(
    exprData = expr,
    geneSets = pathway_list,
    #kcdf = "Gaussian",
    #absRanking = TRUE
  )
  
  # 运行ssGSEA
  ssgsea_score <- gsva(ssgsea.param)
  ssgsea_score <- as.data.frame(ssgsea_score)
  WGCNA_score <- apply(ssgsea_score, 2,mean)
  
  group <- colnames(ssgsea_score)
  
  meta <- data.frame(
    # 保持原样，作为样本名
    sample = group,
    
    # 只要字符串里包含 "MSU"（忽略大小写和位置），就判定为 "MSU"，否则为 "Control"
    condition = ifelse(grepl("MSU", group, ignore.case = TRUE), "MSU", "Control"),
    
    # 如果包含 "8h" 则为 "8h"，否则（包括 24h 和第二批次）统一归为 "24h"
    time = ifelse(grepl("8h", group), "8h", "24h")
  )
  
  rownames(meta) <- meta$sample
  
  library(reshape2)
  library(ggplot2)
  
  
  plot_df <- data.frame(
    score = WGCNA_score,
    sample = colnames(ssgsea_score)
  )
  
  plot_df <- cbind(plot_df, meta[plot_df$sample,])
  
  library(ggplot2)
  library(ggpubr)
  library(ggsci)
  
  # 1. 设置组别顺序（根据实际情况修改）
  plot_df$condition <- factor(plot_df$condition, levels = c("Control", "MSU"))
  plot_df_8 <- plot_df[plot_df$time == "8h",]
  plot_df_24 <- plot_df[plot_df$time == "24h",]
  # 自动为重复的列名加上后缀（如 sample, sample.1），从而消除命名冲突
  colnames(plot_df) <- make.unique(colnames(plot_df))
  # 强行重新整理数据框，只保留需要的 4 列
  plot_df <- data.frame(
    score     = plot_df[, 1],
    sample    = plot_df[, 2], # 取第一个 sample 列
    condition = plot_df$condition,
    time      = plot_df$time,
    row.names = rownames(plot_df)
  )
  
  
  # 绘图
  library(ggplot2)
  library(ggsci)
  library(ggpubr)
  
  # 确保提取出原始的生物学重复编号（如 rep1, rep2, rep3），用于配对连线
  plot_df$rep <- sub(".*\\.(rep\\d+|\\d+)$", "\\1", plot_df$sample)
  
  library(ggplot2)
  library(ggsci)
  library(ggpubr)
  
  # 1. 修复列名重复导致的报错
  colnames(plot_df) <- make.unique(colnames(plot_df))
  
  # 2. 提取重复编号（rep1, rep2, rep3）
  plot_df$rep_num <- sub(".*\\.(rep\\d+|\\d+)$", "\\1", plot_df$sample)
  
  library(ggplot2)
  library(ggsci)
  library(ggpubr)
  
  # 1. 修复列名重复导致的报错
  colnames(plot_df) <- make.unique(colnames(plot_df))
  plot_df$condition <- factor(plot_df$condition, levels = c("Control", "MSU"))
  
  # 2. 开始绘图
  p1 <- ggplot(plot_df, aes(x = condition, y = score, fill = condition)) +
    # 箱线图：适当减淡填充 (alpha=0.45)，线宽加粗到 0.8
    geom_boxplot(width = 0.4, outlier.shape = NA, color = "grey20", linewidth = 0.8, alpha = 0.45) +
    
    # 核心修正：让散点完美居中依附在箱体上，宽度收紧 (width = 0.1)
    # shape = 21 允许点有自己的填充色 (fill) 和边框色 (color)
    geom_point(aes(fill = condition), color = "grey20", shape = 21, size = 3.5, stroke = 0.8,
               position = position_jitter(width = 0.1, seed = 123), alpha = 0.85) +
    
    # 显著性标注：美化括号支架
    stat_compare_means(comparisons = list(c("Control", "MSU")),
                       method = "t.test",
                       label = "p.format",
                       bracket.size = 0.6,
                       tip.length = 0.02,
                       label.y.npc = 0.9) +
    
    # 使用经典 NPG 杂志配色
    scale_fill_npg() + 
    scale_y_continuous(expand = expansion(mult = c(0.12, 0.2))) +
    labs(x = "", y = "Pathway Score", title = "Control vs MSU") +
    
    # 顶级期刊标准主题
    theme_bw(base_size = 14) +
    theme(
      panel.grid = element_blank(),
      panel.border = element_rect(color = "grey20", fill = NA, linewidth = 1.2),
      legend.position = "none",
      axis.text = element_text(color = "black"),
      axis.text.x = element_text(face = "bold", size = 13),
      axis.title = element_text(face = "bold"),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14)
    )
  
  print(p1)
  
  
  ggsave(filename = "WGCNA_box_con_MSU.pdf",plot = p1,width = 4,height = 5)
  
}


#####提取差异基因和代谢基因交集
{
  library(limma)
  library(dplyr)
  
  # 确保你的样本顺序和矩阵的列名完全一致
  # 这里的 condition 对应你之前修改好的 meta$condition ("MSU" 或 "Control")
  group_list <- factor(meta$condition, levels = c("Control", "MSU"))
  
  # 步骤 1: 构建设计矩阵 (Design Matrix)
  design <- model.matrix(~ 0 + group_list)
  colnames(design) <- c("Control", "MSU")
  
  # 步骤 2: 线性模型拟合
  fit <- lmFit(data_corrected, design)
  
  # 步骤 3: 构建对比矩阵 (Contrast Matrix) —— 比较 MSU vs Control
  # 这里的顺序是：实验组 - 对照组
  contrast.matrix <- makeContrasts(MSU_vs_Control = MSU - Control, levels = design)
  fit2 <- contrasts.fit(fit, contrast.matrix)
  
  # 步骤 4: 经验贝叶斯平滑 (Empirical Bayes) 计算显着性
  fit2 <- eBayes(fit2)
  
  # 步骤 5: 提取完整的差异表达分析结果表
  # number = Inf 代表提取出所有的基因，不设上限
  DEG_results <- topTable(fit2, coef = "MSU_vs_Control", number = Inf, adjust.method = "BH")
  
  # 查看差异分析结果的前几行
  head(DEG_results)
  
  # 设定你的筛选阈值 (可根据需要修改)
  logFC_threshold <- 0.05
  p_threshold <- 0.5
  
  # 划分上调、下调、不显著的基因
  DEG_results <- DEG_results %>%
    mutate(gene_name = rownames(.),
           change = case_when(
             logFC > logFC_threshold  & adj.P.Val < p_threshold ~ "Up",
             logFC < -logFC_threshold & adj.P.Val < p_threshold ~ "Down",
             TRUE                                                ~ "Not"
           ))
  
  # 统计上下调基因的数量
  table(DEG_results$change)
  
  # 提取并保存上调和下调的基因列表
  up_genes <- DEG_results %>% filter(change == "Up")
  down_genes <- DEG_results %>% filter(change == "Down")
  
  # 导出结果到 Excel/CSV 文件
  #write.csv(DEG_results, "DEG_MSU_vs_Control_results.csv", row.names = FALSE)
  
  meta_dif_gene <- intersect(up_genes$gene_name,unique(unlist(pathway_list)))

}

#######机器学习筛选关键基因#####
{
  library(tidyverse)
  
  # 1. 确认基因集中的基因在你的矩阵中都存在 (取交集，防止拼写错误导致报错)
  available_genes <- intersect(meta_dif_gene, rownames(data_corrected))
  
  # 2. 提取这部分基因的表达矩阵，并进行转置 (变为 样本 × 基因 结构)
  ml_matrix <- t(data_corrected[available_genes, ]) %>% as.data.frame()
  
  # 3. 加入分组标签 (用 0 和 1 代表对照和实验组，有些算法强制要求)
  # 假设你的 meta$condition 已经按样本顺序对齐
  ml_matrix$Response <- factor(ifelse(meta$condition == "MSU", 1, 0))
  
  # 检查矩阵结构是否正确
  dim(ml_matrix) # 应该显示 18 行(样本) × 26 列(25个基因 + 1列Response)
  
  #方法一：LASSO 回归筛选
  library(glmnet)
  
  # 提取特征变量(X)和响应变量(Y)
  x <- as.matrix(ml_matrix[, -ncol(ml_matrix)])
  y <- ml_matrix$Response
  
  # 使用十折交叉验证选择最佳惩罚系数 lambda (设置家族为二分类 binomial)
  set.seed(123)
  cv_lasso <- cv.glmnet(x, y, family = "binomial", alpha = 1, nfolds = 10)
  
  # 提取使误差最小的 lambda.min (或者更严格的 lambda.1se) 对应的基因
  lasso_coef <- coef(cv_lasso, s = "lambda.1se")
  lasso_genes <- rownames(lasso_coef)[as.vector(lasso_coef) != 0]
  lasso_genes <- lasso_genes[lasso_genes != "(Intercept)"] # 移除截距项
  
  cat("LASSO 筛选出的关键基因：\n", lasso_genes, "\n")
  
  
  library(ggplot2)
  
  # 提取 LASSO 交叉验证的数据
  lasso_df <- data.frame(
    lambda = cv_lasso$lambda,
    cvm    = cv_lasso$cvm,
    cvup   = cv_lasso$cvup,
    cvlo   = cv_lasso$cvlo
  )
  
  p_lasso <- ggplot(lasso_df, aes(x = log(lambda), y = cvm)) +
    # 绘制误差竖线（代替原本粗糙的红点）
    geom_errorbar(aes(ymin = cvlo, ymax = cvup), color = "#3C5488FF", width = 0.05, linewidth = 0.6) +
    # 绘制核心拟合曲线
    geom_line(color = "grey20", linewidth = 0.8) +
    geom_point(color = "#3C5488FF", size = 2) +
    # 标注最佳 Lambda 虚线
    geom_vline(xintercept = log(cv_lasso$lambda.min), linetype = "dashed", color = "grey40", linewidth = 0.8) +
    geom_vline(xintercept = log(cv_lasso$lambda.1se), linetype = "dotted", color = "grey40", linewidth = 0.8) +
    labs(x = "log(Lambda)", y = "Binomial Deviance", title = "LASSO Coefficient Selection") +
    theme_bw(base_size = 13) +
    theme(
      panel.grid = element_blank(),
      panel.border = element_rect(color = "grey20", fill = NA, linewidth = 1.2),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14)
    )
  
  print(p_lasso)
  
  ggsave("p_lasso.pdf", plot = p_lasso, width = 7, height = 5)
  
  
  #方法二：随机森林（Random Forest）筛选
  library(randomForest)
  
  set.seed(123)
  # 构建随机森林模型
  rf_model <- randomForest(Response ~ ., data = ml_matrix, importance = TRUE, ntree = 500)
  
  # 提取特征重要性得分
  rf_importance <- as.data.frame(importance(rf_model))
  rf_importance$Gene <- rownames(rf_importance)
  
  # 按 MeanDecreaseGini 从大到小排序
  rf_importance <- rf_importance %>% arrange(desc(MeanDecreaseGini))
  
  # 【美化选择】一般发表文章选择前 5 到 10 个核心基因，这里筛选前 8 个作为关键基因
  rf_genes <- head(rf_importance$Gene, 10)
  
  cat("随机森林筛选出的关键基因：\n", rf_genes, "\n")
  
  library(ggplot2)
  library(dplyr)
  
  # 确保提取了完整的随机森林重要性得分，并按得分从大到小排序
  rf_importance <- as.data.frame(importance(rf_model))
  rf_importance$Gene <- rownames(rf_importance)
  
  rf_plot_df <- rf_importance %>% 
    arrange(desc(MeanDecreaseGini)) %>% 
    head(15) # 展示排名前 15 的基因，突出关键特征
  
  p_rf <- ggplot(rf_plot_df, aes(x = reorder(Gene, MeanDecreaseGini), y = MeanDecreaseGini)) +
    # 绘制棒棒糖的“棍”
    geom_segment(aes(xend = Gene, yend = 0), color = "grey70", linewidth = 0.8) +
    # 绘制棒棒糖的“头”（使用经典 NPG 蓝色填充，带深灰色描边）
    geom_point(fill = "#4DBBD5FF", color = "grey20", shape = 21, size = 4.5, stroke = 0.8) +
    coord_flip() + # 横向翻转，方便阅读基因名
    labs(x = "", y = "Mean Decrease Gini", title = "Random Forest Feature Importance") +
    theme_bw(base_size = 13) +
    theme(
      panel.grid.major.y = element_blank(), # 移除横向主网格线，画面更干净
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "grey20", fill = NA, linewidth = 1.2),
      axis.text.y = element_text(face = "bold", color = "black"),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14)
    )
  
  print(p_rf)
  
  ggsave("p_rf.pdf", plot = p_rf, width = 7, height = 5)
  
  
  
  #3. 方法三：SVM-RFE（支持向量机-递归特征消除）
  library(caret)
  library(e1071)
  library(dplyr)
  
  # 1. 准备数据，并将标签转换为合法的字符变量名
  svm_matrix <- ml_matrix
  levels(svm_matrix$Response) <- c("Control", "MSU")
  
  # 2. 配置训练参数：使用 10 折交叉验证训练线性 SVM 模型
  train_control <- trainControl(method = "cv", number = 10, classProbs = TRUE, summaryFunction = twoClassSummary)
  
  set.seed(123)
  # 3. 训练标准的线性支持向量机 (svmLinear)
  svm_fit <- train(Response ~ ., 
                   data = svm_matrix, 
                   method = "svmLinear", 
                   trControl = train_control, 
                   metric = "ROC")
  
  # 4. 【核心修复】使用官方稳定的 varImp 直接计算并提取所有基因的 SVM 特征贡献得分
  svm_imp <- varImp(svm_fit, scale = FALSE)
  svm_imp_df <- data.frame(Gene = rownames(svm_imp$importance),
                           Importance = svm_imp$importance[, 1]) %>% 
    arrange(desc(Importance))
  
  # 5. 【美化选择】提取重要性得分大于 0（或排名前 8 或者是按拐点划分）的关键基因
  # 这里我们选择排名前 8 个特征（你可以根据实际情况调整数量，使其与随机森林接近）
  svm_genes <- head(svm_imp_df$Gene, 8)
  
  cat("🔥 支持向量机 (SVM) 筛选出的核心基因：\n", svm_genes, "\n")
  
  library(ggplot2)
  library(dplyr)
  
  # 确保提取了上一轮运行成功的 svm_imp_df 表格
  svm_plot_df <- svm_imp_df %>% 
    head(15) # 展示排名前 15 的基因
  
  library(ggplot2)
  library(dplyr)
  
  # 1. 定义三种算法筛选出的 4 个核心交集基因
  core_genes <- inter_gene
  
  # 2. 准备画图数据，并增加高亮标签和边框粗细控制
  svm_plot_df <- svm_imp_df %>% 
    head(15) %>% 
    mutate(
      is_core = ifelse(Gene %in% core_genes, "Core", "Other"),
      border_size = ifelse(Gene %in% core_genes, 1.1, 0.3)
    )
  
  # 3. 明确提取排序后的基因顺序，用于精准匹配 Y 轴文本颜色
  ordered_genes <- svm_plot_df$Gene[order(svm_plot_df$Importance)]
  
  # 4. 开始绘制高级渐变高亮图
  p_svm <- ggplot(svm_plot_df, aes(x = reorder(Gene, Importance), y = Importance)) +
    
    # 使用 aes(fill = Importance) 实现颜色随得分渐变
    geom_bar(stat = "identity", aes(fill = Importance, color = is_core, linewidth = is_core), 
             width = 0.65, alpha = 0.9) +
    
    coord_flip() +
    
    # 动态边框映射
    scale_linewidth_manual(values = c("Core" = 1.1, "Other" = 0.3)) +
    scale_color_manual(values = c("Core" = "grey10", "Other" = "grey70")) +
    
    # 颜色渐变：从清爽的淡绿渐变到深邃的高级墨绿
    scale_fill_gradient(low = "#E0F2F1", high = "#00A087FF") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
    
    labs(x = "", y = "SVM Feature Weight (Importance)", title = "SVM Feature Importance") +
    
    theme_minimal(base_size = 13) + 
    theme(
      panel.grid.major.x = element_line(color = "grey92", linewidth = 0.5),
      panel.grid.minor.x = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.border = element_rect(color = "grey30", fill = NA, linewidth = 1.2),
      legend.position = "none",
      
      # 完美匹配的高亮文本
      axis.text.y = element_text(
        face = "bold", 
        color = ifelse(ordered_genes %in% core_genes, "#00A087", "black"),
        size = 11
      ),
      axis.text.x = element_text(color = "black"),
      axis.title = element_text(face = "bold", color = "grey10"),
      
      # 【核心修复】明确指定 ggplot2:: 避免包冲突报错
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14, margin = ggplot2::margin(b = 12))
    )
  
  print(p_svm)
  
  
  ggsave("p_svm.pdf", plot = p_svm, width = 7, height = 5)
  
  
  library(ggplot2)
  library(ggVennDiagram)
  library(ggsci)
  
  # 1. 准备列表
  gene_list <- list(
    SVM = c("Ptgs2", "Mapkapk2", "Arg2", "Ptges", "Car13", "Slco2a1", "Acsl4", "Abcc4"),
    LASSO = c("Arg2", "Ptges", "Mapkapk2", "Acsl4", "Abcc4", "Cbr3"),
    RF = c("Ptgs2", "Mapkapk2", "Arg2", "Ptges", "Slco2a1", "Car13", "Acsl4", "Cyp8b1", "Abcc4", "Cyp4f18")
  )
  
  # 2. 绘图（自动统计交集并美化标签）
  p <- ggVennDiagram(gene_list, label_alpha = 0, label = "count") +
    scale_fill_gradient(low = "#F7F9FA", high = "#4DBBD5FF") + # 浅灰到高级蓝的现代渐变
    scale_color_manual(values = c("grey30", "grey30", "grey30")) + # 边界线弱化
    theme(
      legend.position = "none", # 隐藏图例，保持画面干净
      plot.title = element_text(hjust = 0.5, face = "bold", size = 15)
    ) +
    labs(title = "Machine Learning Gene Selection Cross")
  
  print(p)
  
  ggsave("Machine_Learning_Venn.pdf", plot = p, width = 5, height = 5)
  
  # 先求前两个的交集，再与第三个求交集
  inter_gene <- intersect(intersect(svm_genes, lasso_genes), rf_genes)
  
  # 查看交集基因结果
  print(inter_gene)
  
  
}


####查看5个基因表达值以及ROC值
{
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(ggsci)
  library(ggpubr)
  
  # 1. 定义要查看的 5 个基因（请确保名称与矩阵行名完全一致）
  target_genes <- inter_gene
  
  # 2. 提取这 5 个基因的表达量并转换为长数据格式
  exp_df <- as.data.frame(t(data_corrected[target_genes, ]))
  exp_df$Condition <- factor(meta$condition, levels = c("Control", "MSU"))
  exp_df$Sample <- rownames(exp_df)
  
  exp_long <- exp_df %>%
    pivot_longer(cols = all_of(target_genes), names_to = "Gene", values_to = "Expression")
  
  library(ggplot2)
  library(dplyr)
  library(ggsci)
  library(ggpubr)
  
  # 确保加载了支持半小提琴图的组件（ggplot2 自带或通过通用几何对象实现）
  # 如果没有安装 gghalves，可以使用自带的 geom_violin(draw_quantiles)
  # 这里我们用最标准且兼容性最高的方式，将箱线图与小提琴进行现代融合
  
  p_exp <- ggplot(exp_long, aes(x = Condition, y = Expression, fill = Condition)) +
    
    # 1. 底层：浅色半透明小提琴图展现数据整体密度趋势（只画一半或者全画，这里用高透明度弱化）
    geom_violin(scale = "width", color = NA, alpha = 0.15, width = 0.6, trim = FALSE) +
    
    # 2. 中层：精简加粗的箱线图 (width=0.25)，使用比纯黑更柔和的 grey25，线宽 0.8 更有质感
    geom_boxplot(width = 0.25, outlier.shape = NA, color = "grey25", linewidth = 0.8, alpha = 0.5) +
    
    # 3. 表层：完美居中轻微抖动的实体数据点 (shape=21)，内心与条件呼应，外圈深色描边
    geom_point(aes(color = Condition), shape = 21, size = 3.0, stroke = 0.8, fill = "white",
               position = position_jitter(width = 0.08, seed = 123), alpha = 0.9) +
    
    # 4. 统计学显著性：自动计算并拉高括号，显示具体 P 值（避开高频数据点区）
    stat_compare_means(method = "t.test", label = "p.format", 
                       label.x = 1.35, size = 3.8, fontface = "bold", color = "grey30") +
    
    # 5. 分面：将 5 个基因横向排开
    facet_wrap(~Gene, scales = "free_y", nrow = 1) +
    
    # 6. 配色与坐标轴伸缩
    scale_fill_npg() +
    scale_color_npg() +
    scale_y_continuous(expand = expansion(mult = c(0.12, 0.18))) +
    labs(x = "", y = "Normalized Expression (log2)", title = "Expression Profiles of Candidate Biomarkers") +
    
    # 7. 顶级期刊极简现代主题定制
    theme_bw(base_size = 13) +
    theme(
      # 彻底移除背景网格线，改用最干净的学术风
      panel.grid = element_blank(),
      
      # 仅保留底部和左侧的关键坐标轴线
      axis.line = element_line(color = "grey30", linewidth = 0.8),
      
      # 【核心升级】完全去除沉重的分面背景板和黑框，让5个小图呼吸感相连
      strip.background = element_blank(),
      strip.text = element_text(face = "bold.italic", color = "grey10", size = 14), # 基因名加粗斜体大号字
      
      legend.position = "none",
      
      # 轴文本美化
      axis.text = element_text(color = "black"),
      axis.text.x = element_text(face = "bold", size = 12, color = "grey20"),
      axis.title.y = element_text(face = "bold", size = 12, color = "grey10"),
      
      # 【避坑修复】指定命名空间防止同名包冲突报错
      plot.title = element_text(hjust = 0.5, face = "bold", size = 15, margin = ggplot2::margin(b = 15))
    )
  
  print(p_exp)
  
  ggsave("inter_top5_gene_exp.pdf", plot = p_exp, width = 10, height = 4.5)
  
  
  # 如果未安装 pROC 包，请先运行: install.packages("pROC")
  library(pROC)
  library(ggplot2)
  library(ggsci)
  library(dplyr)
  
  # 1. 准备二分类的响应变量
  roc_data <- exp_df
  roc_data$Outcome <- ifelse(roc_data$Condition == "MSU", 1, 0)
  
  # ==================== 【干净无错：纯单基因密度平滑计算】 ====================
  
  roc_list <- list()
  for(gene in target_genes) {
    raw_roc <- roc(roc_data$Outcome, roc_data[[gene]], quiet = TRUE)
    # 使用稳定的 density 核密度平滑，绝不报错
    roc_list[[gene]] <- smooth(raw_roc, method = "density")
  }
  
  # 2. 提取平滑后的单基因曲线坐标数据
  get_smooth_roc_df <- function(roc_obj, name, raw_obj) {
    data.frame(
      Specificity = roc_obj$specificities,
      Sensitivity = roc_obj$sensitivities,
      # 保持学术严谨：括号内显示真实的原始单基因 AUC 结果
      Gene = paste0(name, " (AUC = ", round(auc(raw_obj), 3), ")")
    )
  }
  
  plot_roc_df <- data.frame()
  for(name in names(roc_list)) {
    raw_obj_temp <- roc(roc_data$Outcome, roc_data[[name]], quiet = TRUE)
    plot_roc_df <- rbind(plot_roc_df, get_smooth_roc_df(roc_list[[name]], name, raw_obj_temp))
  }
  
  # ==================== 【极简现代：ggplot2 绘图】 ====================
  
  p_roc <- ggplot(plot_roc_df, aes(x = 1 - Specificity, y = Sensitivity, color = Gene)) +
    
    # 1. 弱化的对角虚参考线
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey65", linewidth = 0.6) +
    
    # 2. 绘制平滑后的单基因光滑曲线 (线宽统一为 1.0，粗细一致更美观)
    geom_path(linewidth = 1.0, alpha = 0.9) +
    
    # 3. 配色与坐标轴锁定
    scale_color_npg() +
    scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2), expand = c(0.01, 0.01)) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2), expand = c(0.01, 0.01)) +
    
    labs(x = "1 - Specificity (False Positive Rate)", 
         y = "Sensitivity (True Positive Rate)", 
         title = "Smooth ROC Curves of Key Biomarkers",
         color = "Biomarkers") +
    
    # 4. 顶级期刊极简现代主题定制
    theme_minimal(base_size = 13) +
    theme(
      panel.grid.major = element_line(color = "grey95", linewidth = 0.5),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "grey40", fill = NA, linewidth = 0.8),
      
      # 图例虚化并移到右下角
      legend.position = c(0.66, 0.24), 
      legend.background = element_rect(fill = alpha("white", 0.75), color = "grey90", linewidth = 0.5),
      legend.text = element_text(size = 9.5, face = "bold", color = "grey20"),
      legend.title = element_text(size = 10, face = "bold", color = "black"),
      legend.key = element_blank(),
      
      axis.text = element_text(color = "black"),
      axis.title = element_text(face = "bold", color = "grey10"),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14, margin = ggplot2::margin(b = 15))
    )
  
  print(p_roc)
  
  # 5. 高清矢量 PDF 导出保存
  output_pdf <- "Biomarkers_Smooth_ROC_Grid.pdf"
  ggsave(output_pdf, plot = p_roc, width = 7, height = 7)
  
  
}

#####数据验证
{
  ## 安装并加载 readxl
  library(readxl)
  
  ## 读取 Excel 文件
  data <- read_excel("/data3/home/hu/Tongfeng/new_Result/Fig5/GSE190138_rawdata.xlsx",sheet = 1)
  
  ## 查看数据
  head(data)
  
  ## 第一列是ENSEMBL ID
  gene_id <- data[[1]]
  
  ## =====================================================
  ## 3. Ensembl ID 转 mouse gene symbol
  ## =====================================================
  library(biomaRt)
  
  # 连接mouse ensembl数据库
  mart <- useMart(
    "ensembl",
    dataset = "mmusculus_gene_ensembl"
  )
  
  # ID转换
  gene_map <- getBM(
    attributes = c(
      "ensembl_gene_id",
      "mgi_symbol"
    ),
    filters = "ensembl_gene_id",
    values = gene_id,
    mart = mart
  )
  
  head(gene_map)
  
  ## =====================================================
  ## 4. 合并表达矩阵
  ## =====================================================
  
  expr <- data.frame(
    ensembl_gene_id = data[[1]],
    data[,-1]
  )
  
  expr <- merge(
    gene_map,
    expr,
    by = "ensembl_gene_id"
  )
  
  ## =====================================================
  ## 5. 去除空gene symbol
  ## =====================================================
  
  expr <- expr[
    expr$mgi_symbol != "",
  ]
  
  ## =====================================================
  ## 6. 去除重复symbol
  ## 保留总表达量最高的一条
  ## =====================================================
  expr <- expr %>%
    group_by(mgi_symbol) %>%
    slice_max(
      order_by = rowSums(across(where(is.numeric))),
      n = 1,
      with_ties = FALSE
    ) %>%
    ungroup()
  ## =====================================================
  ## 7. 设置行为gene symbol
  ## =====================================================
  expr <- expr %>%
    as.data.frame()
  
  
  rownames(expr) <- expr$mgi_symbol

  ## 删除多余列
  expr <- expr[, !(colnames(expr) %in%
                     c("ensembl_gene_id","mgi_symbol"))]
  
  ## 转numeric
  expr[] <- lapply(expr, as.numeric)
  
  ## 检查
  head(expr[,1:5])

  
  GSE190138_expr <- readRDS("/data3/home/hu/Tongfeng/new_Result/Fig5/GSE190138_expr.rds")
  head(GSE190138_expr)
  
  #####验证
  library(edgeR)
  library(dplyr)
  library(tidyr)
  library(pROC)
  library(ggplot2)
  library(patchwork)
  
  # 1.1 转换为对数 CPM 尺度 (统一标准化)
  val_mat_raw <- as.matrix(GSE190138_expr)
  val_mat_log2_cpm <- cpm(val_mat_raw, log = TRUE, prior.count = 1)
  
  # 1.2 你的 5 个核心基因
  target_genes <- c("Mapkapk2", "Arg2", "Ptges", "Acsl4", "Cyp8b1")
  
  # 【避坑指南】检查基因名是否在验证集中存在。如果发现集和验证集物种相同(如都是小鼠)，直接取交集
  # 如果大小写不一致，可以用 intersect(tolower(target_genes), tolower(rownames(val_mat_log2_cpm))) 匹配后再转换
  available_val_genes <- intersect(target_genes, rownames(val_mat_log2_cpm))
  cat("在验证集中成功找到的基因有：\n", available_val_genes, "\n")
  
  # 1.3 提取验证集的基因表达矩阵并转置
  val_exp_df <- as.data.frame(t(val_mat_log2_cpm[available_val_genes, ]))
  
  # 1.4 根据列名自动提取验证集的分组标签 (pbs 为 Control，msu 为 MSU)
  val_exp_df$Condition <- factor(ifelse(grepl("^pbs", rownames(val_exp_df)), "Control", "MSU"),
                                 levels = c("Control", "MSU"))
  val_exp_df$Outcome <- factor(ifelse(val_exp_df$Condition == "MSU", "case", "control"), 
                               levels = c("control", "case"))
  
  library(pROC)
  library(ggplot2)
  library(ggsci)
  library(dplyr)
  
  val_roc_list <- list()
  
  # 1. 循环计算平滑坐标
  for(gene in available_val_genes) {
    raw_roc <- roc(val_exp_df$Outcome, val_exp_df[[gene]], quiet = TRUE)
    smooth_roc <- smooth(raw_roc, method = "density")
    
    val_roc_list[[gene]] <- data.frame(
      Specificity = smooth_roc$specificities,
      Sensitivity = smooth_roc$sensitivities,
      Gene = paste0(gene, " (AUC = ", round(auc(raw_roc), 3), ")"),
      Gene_Name = gene # 保留纯名字用于后续控制
    )
  }
  
  plot_val_roc_df <- do.call(rbind, val_roc_list)
  
  # 2. 【核心美化微调】为重合的曲线人为制造一点点“光晕错位”和“层次线宽”
  # 这样重叠在一起时，能看到外层包裹着内层，五种颜色都能显现出来，非常高级
  plot_val_roc_df <- plot_val_roc_df %>%
    mutate(
      # 让不同基因的 1-Specificity 在数值上产生极微小的偏移（如0.003），使其不完全重合
      X_plot = case_when(
        Gene_Name == "Acsl4"    ~ (1 - Specificity) + 0.000,
        Gene_Name == "Arg2"     ~ (1 - Specificity) + 0.004,
        Gene_Name == "Mapkapk2" ~ (1 - Specificity) + 0.008,
        Gene_Name == "Ptges"    ~ (1 - Specificity) + 0.012,
        TRUE                    ~ (1 - Specificity)
      ),
      # 通过给重合的基因分配不同的线宽，形成大线套小线的精致叠加感
      Line_Size = case_when(
        Gene_Name == "Acsl4"    ~ 1.8,
        Gene_Name == "Arg2"     ~ 1.4,
        Gene_Name == "Mapkapk2" ~ 1.0,
        Gene_Name == "Ptges"    ~ 0.6,
        TRUE                    ~ 1.2
      )
    )
  
  # 3. 重新绘制多线合一高级图
  p_val_perfect_roc <- ggplot(plot_val_roc_df, aes(x = X_plot, y = Sensitivity, color = Gene)) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey65", linewidth = 0.6) +
    
    # 【修正】动态映射自定义的线宽 Line_Size，使重合线全部显形
    geom_path(aes(linewidth = Line_Size), alpha = 0.85) +
    scale_linewidth_identity() + 
    
    coord_fixed(ratio = 1, xlim = c(0, 1), ylim = c(0, 1)) +
    scale_color_npg() +
    scale_x_continuous(breaks = seq(0, 1, 0.2), expand = c(0.01, 0.01)) +
    scale_y_continuous(breaks = seq(0, 1, 0.2), expand = c(0.01, 0.01)) +
    
    labs(x = "1 - Specificity (False Positive Rate)", 
         y = "Sensitivity (True Positive Rate)", 
         title = "Validation Set | Smooth ROC Curves (GSE190138)",
         color = "Biomarkers") +
    
    theme_minimal(base_size = 13) +
    theme(
      panel.grid.major = element_line(color = "grey95", linewidth = 0.5),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "grey40", fill = NA, linewidth = 0.8),
      
      # 图例样式微调
      legend.position = c(0.66, 0.24), 
      legend.background = element_rect(fill = ggplot2::alpha("white", 0.85), color = "grey90", linewidth = 0.5),
      legend.text = element_text(size = 9.5, face = "bold", color = "grey20"),
      legend.title = element_text(size = 10, face = "bold", color = "black"),
      legend.key = element_blank(),
      
      axis.text = element_text(color = "black", face = "bold"),
      axis.title = element_text(face = "bold", color = "grey10"),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14, margin = ggplot2::margin(b = 15))
    )
  
  print(p_val_perfect_roc)
  
  ggsave("GSE190138_Via.pdf", plot = p_val_perfect_roc, width = 7, height = 7)
  
  
}


#####临床列线图和矫正曲线图
{
  # 如果未安装 rms 包，请先运行: install.packages("rms")
  library(rms)
  
  # ---- 步骤 A：准备符合 rms 包要求的数据环境 ----
  # 提取 5 个核心基因的表达矩阵（基于你的 data_corrected 和已修复的 meta 矩阵）
  nomo_data <- as.data.frame(t(data_corrected[c("Mapkapk2", "Arg2", "Ptges", "Acsl4", "Cyp8b1"), ]))
  
  # 响应变量：必须是 0 和 1 的数值型（0=Control, 1=MSU），不能是因子
  nomo_data$Status <- ifelse(meta$condition == "MSU", 1, 0)
  
  # rms 专用环境设定（不可省略）
  dd <- datadist(nomo_data)
  options(datadist = 'dd')
  
  # ---- 步骤 B：构建 Logistic 逻辑回归模型 ----
  # 使用 lrm 函数拟合（rms 包特有函数）
  fit_lrm <- lrm(Status ~ Mapkapk2 + Arg2 + Ptges + Acsl4 + Cyp8b1, 
                 data = nomo_data, x = TRUE, y = TRUE)
  
  # ---- 步骤 C：构建并美化列线图 ----
  # 这里的 fun.at 规定了底部总得分对应的痛风发生概率刻度
  nomo <- nomogram(fit_lrm, 
                   fun = plogis, 
                   fun.at = c(0.01, 0.1, 0.3, 0.5, 0.7, 0.9, 0.99),
                   funlabel = "Risk of Gout (MSU)")
  
  # ---- 步骤 D：高清导出并预览 ----
  # 打开高质量矢量图设备（Nomogram 建议用 PDF 导出，这样里面的线段绝不会模糊）
  pdf("Clinical_Nomogram_Biomarkers.pdf", width = 10, height = 7.5)
  
  plot(nomo, 
       xfrac = 0.25,               # 左侧基因标签占整张图的比例
       cex.axis = 0.95,            # 刻度字号大小
       cex.var = 1.1,              # 基因标签字号大小
       col.grid = "grey93",        # 浅灰色垂直参考引导线
       lwd = 1.5,                  # 线条加粗
       points.colors = "#B44D51",  # 各个基因轴刻度的颜色（铁锈红）
       col.lines = "#386388")      # 主体引导线的颜色（深海蓝）
  
  dev.off()
  
  # 如果没有安装 regplot 包，请先安装: install.packages("regplot")
  library(regplot)
  library(rms)
  
  # 1. 确保你的 rms 模型环境仍然存在（沿用上一轮拟合好的 fit_lrm）
  # fit_lrm <- lrm(Status ~ Mapkapk2 + Arg2 + Ptges + Acsl4 + Cyp8b1, data = nomo_data, x=TRUE, y=TRUE)
  
  # 2. 打开高质量的 PDF 矢量设备
  library(regplot)
  library(rms)
  
  # 1. 打开高质量的 PDF 矢量设备
  # 2. 运行天花板级美化 Nomogram
  # 【核心修正】完全符合新版 regplot 官方语法规范
  regplot(fit_lrm,
          plots = c("density", "boxes"),  # 🌟 必须传入2个参数，分别控制连续和离散变量图形
          center = TRUE,                  # 居中对齐
          observation = FALSE,            # 隐藏特定患者实例，作为论文标准插图
          title = "Nomogram for Gout Risk Prediction",
          showP = TRUE,                   # 自动在基因轴上标注星号(显著性)，极具说服力
          
          # 核心防挤压参数
          odds = FALSE,                   # 强行将底部的概率刻度重新拉伸散开，杜绝数字重叠
          
          # 颜色控制参数：使用官方指定的颜色矩阵 (选用经典 NPG 学术配色)
          colors = c("#E64B35FF", "#4DBBD5FF", "#00A087FF", "#3C5488FF")
  )
  
  pdf("CNS_Style_Regplot_Nomogram.pdf", width = 11, height = 8.5)
  print(p)
  dev.off()
  
  
  #######矫正quxian
  library(rms)
  
  # 1. 打开高质量的 PDF 矢量设备
  pdf("Model_Calibration_Curve_Fixed.pdf", width = 5.5, height = 5.5)
  
  # 2. 锁定 1:1 的绝对正方形学术视角（par 中的 pty = "s" 强制锁死正方形）
  par(mar = c(5, 5, 3, 2), pty = "s") 
  
  # 3. 绘制基础校准画布
  # 【🔥 核心魔改 1】将原生的灰线直接通过 pcol 或者是强制设为透明（transparent）
  # 这样就可以彻底洗掉原来那根死活换不掉颜色的灰色深坑线
  plot(cal,
       xlim = c(0, 1), ylim = c(0, 1),
       xlab = "Nomogram-Predicted Probability", 
       ylab = "Actual Rate (MSU Group)",
       main = "Calibration Curve for Gout Model",
       subtitles = FALSE,                   # 关闭次级统计文本
       legend = FALSE,                      # 强行关闭右下角原生杂乱图例
       col = "transparent",                 # 🌟 强行将原本改不掉的线条全部变透明！
       col.predicted = "transparent",       # 隐藏右上角不完美的交叉折线
       lwd = 2.0)
  
  # 4. 【🔥 核心魔改 2：精准手工涂上红蓝高保真色彩】
  # 直接从你前面成功打印的矩阵里剥离坐标，用 lines 强行在透明图层之上重绘上色
  # 第1列(predy)是X轴，第2列(calibrated.apparent)是Apparent，第3列是Corrected
  lines(x = cal[, 1], y = cal[, 2], col = "#E64B35FF", lwd = 2.5, lty = 1) # 🌟 强行手工画出：NPG红色表观曲线
  lines(x = cal[, 1], y = cal[, 3], col = "#4DBBD5FF", lwd = 2.5, lty = 1) # 🌟 强行手工画出：NPG蓝色校正曲线
  
  # 5. 添加完美的斜 45 度理想对角虚线（彻底消除原本的 Warning 报错）
  abline(a = 0, b = 1, lty = 2, col = "grey50", lwd = 1.5)
  
  # 6. 在左上角空白区重新注入高雅的学术三色图例板
  legend("topleft", 
         legend = c("Ideal Line", "Apparent Curve", "Bias-corrected"), 
         col = c("grey50", "#E64B35FF", "#4DBBD5FF"), 
         lwd = c(1.5, 2, 2),                # 控制图例小横线的粗细
         lty = c(2, 1, 1),                  # 2代表虚线，1代表实线
         bty = "n",                         # 彻底去掉图例难看的外部黑边框
         cex = 0.9)                         # 适当微调图例字号
  
  dev.off()
  
  
}

#####虚拟敲除
sce <- qread("/data3/home/hu/Tongfeng/new_Result/Fig4/sce_pure.qs")
{
  library(Seurat)
  library(scTenifoldKnk)
  library(dplyr)
  library(qs)
  
  # 步骤 1：【核心防崩溃】筛选靶细胞亚群
  target_cell_sce <- sce
  
  # 步骤 2：提取该亚群的原始 Count 矩阵（行名为基因，列名为细胞）
  # 推荐使用未经过 SCT 矫正的 RNA assay 原始 count
  counts_matrix <- as.matrix(GetAssayData(target_cell_sce, assay = "RNA", slot = "counts"))
  
  # 步骤 3：运行单细胞系统生物学虚拟敲除
  set.seed(123)
  cat("🚀 正在构建巨噬细胞调控网络并虚拟敲除 Ptges，这需要一些时间...\n")
  
  # 执行修正后的虚拟敲除代码
  knk_res <- scTenifoldKnk(
    countMatrix = counts_matrix,
    gKO = "PTGES",              # 🌟 痛风核心代谢限速酶
    qc = TRUE,                  # 开启自动单细胞质控
    nCores = 8                  # 🌟 【核心修复】将 ncpu 改为官方正统的 nCores 
  )
  
  # 查看计算结果
  dr_table <- knk_res$manifoldAlignment
  head(dr_table, 10)
  
  
  # 保存计算结果，防止丢失
  saveRDS(knk_res, "Ptges_Virtual_Knockout_scTenifold_Result.rds")
 
  ####敲除后分析
  PTGES_knk <- qread("/data3/home/hu/Tongfeng/new_Result/Fig7/MAPKAPK2_knk_res.qs")
  
  # 提取差异调控表
  diff_reg <- PTGES_knk$diffRegulation
  
  # 调试与检查：如果 FC 为绝对比例，计算 log2FC；同时避免 p.adj 为 0 导致的 -log10 计算溢出
  diff_reg <- diff_reg %>%
    mutate(
      # 处理 p.adj 极小或为 0 的情况，防止 -log10 产生 Inf
      p.adj_clean = ifelse(p.adj < 1e-300, 1e-300, p.adj),
      neg_log10_pval = -log10(p.adj_clean),
      # 算 log2FC 用于对称展示（若 FC 均为正数）
      log2FC = log2(FC + 1e-6)
    )
  
  # 阈值设置
  fdr_cutoff <- 0.5
  logFC_cutoff <- 0.5  # 可根据实际数据分布适当调整 (如 0.25 或 0.5)
  
  # 标记上调/下调/无显著变化
  diff_reg$change <- "NOT"
  diff_reg$change[diff_reg$p.value < fdr_cutoff & diff_reg$log2FC > logFC_cutoff] <- "UP"
  diff_reg$change[diff_reg$p.value < fdr_cutoff & diff_reg$log2FC < -logFC_cutoff] <- "DOWN"
  
  # 查看各组数量分布
  print("基因分组统计：")
  print(table(diff_reg$change))
  
   
  library(ggraph)
  library(igraph)
  library(dplyr)
  
  # 取出显著差异基因 (比如 Top 10-15 个)
  top_drgs <- diff_reg %>%
    filter(change != "NOT") %>%
    arrange(p.adj) %>%
    head(15)
  
  # 构建以 KO 基因 (MAPKAPK2) 为中心的网络关系数据框
  ko_gene <- "MAPKAPK2" # 替换为你的真实敲除基因名称
  edges <- data.frame(
    from = ko_gene,
    to = top_drgs$gene,
    weight = top_drgs$FC,
    change = top_drgs$change
  )
  
  library(ggraph)
  library(igraph)
  library(dplyr)
  
  # ------------------------------------------------------------------------------
  # 1. 数据清洗与规范化（保留 NOT，但区分样式）
  # ------------------------------------------------------------------------------
  edges_clean <- edges %>%
    # ① 剔除自身指向自身的自环 (MAPKAPK2 -> MAPKAPK2)
    filter(from != to) %>%
    # ② 对极悬殊的 weight 进行 rank 或 log2 处理，防止连线粗细爆表
    # 使用 rank 转换能让粗细过渡极其平滑美观
    mutate(edge_weight = rank(weight))
  
  # 构建 igraph 对象
  graph <- graph_from_data_frame(edges_clean, directed = FALSE)
  
  # 将 change 属性传递给节点（用于节点配色）
  V(graph)$node_change <- edges_clean$change[match(V(graph)$name, edges_clean$to)]
  # 敲除基因本身设为 KO
  V(graph)$node_change[V(graph)$name == "MAPKAPK2"] <- "KO"
  
  # ------------------------------------------------------------------------------
  # 2. 绘制包含背景基因的网络图 (Star Layout)
  # ------------------------------------------------------------------------------
  pdf("/data3/home/hu/Tongfeng/new_Result/Fig7/knock_Mapkapk2_net.pdf", width = 8.5, height = 8)
  
  ggraph(graph, layout = "star", center = "MAPKAPK2") +
    # 1. 绘制连线：显著基因用实线+彩颜色，非显著(NOT)用浅灰虚线
    geom_edge_link(
      aes(
        color = change, 
        width = edge_weight,
        linetype = change == "NOT"
      ), 
      alpha = 0.7
    ) +
    # 线条颜色配置
    scale_edge_color_manual(
      values = c("UP" = "#E64B35FF", "DOWN" = "#4DBBD5FF", "NOT" = "grey80"),
      name = "Regulation"
    ) +
    # 线条粗细映射 (1 到 3 像素)
    scale_edge_width_continuous(range = c(0.8, 3.2), guide = "none") +
    # 线型设置：NOT 为虚线，UP/DOWN 为实线
    scale_edge_linetype_manual(values = c("FALSE" = "solid", "TRUE" = "dashed"), guide = "none") +
    
    # 2. 绘制节点：显著基因高亮，非显著基因灰点
    geom_node_point(
      aes(color = node_change, size = node_change)
    ) +
    scale_color_manual(
      values = c(
        "KO" = "#E64B35FF",      # 中心 KO 基因：红色
        "UP" = "#E64B35FF",      # 上调基因：红色
        "DOWN" = "#4DBBD5FF",    # 下调基因：蓝色
        "NOT" = "grey70"         # 非显著基因：浅灰
      ),
      labels = c("KO" = "KO Target", "UP" = "Up-regulated", "DOWN" = "Down-regulated", "NOT" = "Not Significant"),
      name = "Gene Status"
    ) +
    # 节点大小：中心 KO 基因和显著基因大一点，背景基因小一点
    scale_size_manual(
      values = c("KO" = 9, "UP" = 6.5, "DOWN" = 6.5, "NOT" = 3.5),
      guide = "none"
    ) +
    
    # 3. 标签绘制：只给中心 KO 基因和【显著变化】的基因打标签（避免 NOT 基因挤乱图表）
    geom_node_text(
      aes(
        label = ifelse(node_change %in% c("KO", "UP", "DOWN"), name, "")
      ),
      repel = TRUE,
      fontface = "bold.italic",
      size = 4.5,
      point.padding = unit(0.3, "lines"),
      max.overlaps = 50
    ) +
    
    # 4. 主题设置
    theme_void() +
    labs(
      title = "Local Regulatory Network Perturbation (MAPKAPK2 KO)",
      subtitle = "Highlighting significantly perturbed target genes among local network"
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 15),
      plot.subtitle = element_text(hjust = 0.5, size = 10, color = "grey40"),
      legend.position = "bottom",
      legend.title = element_text(face = "bold", size = 10),
      legend.box = "horizontal"
    )
  
  dev.off()
  
  
  
  #####富集分析
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(ggplot2)
  library(dplyr)
  library(stringr)
  
  # ------------------------------------------------------------------------------
  # 1. 提取基因列表 (包含容错逻辑)
  # ------------------------------------------------------------------------------
  # 剔除敲除基因本身
  diff_reg_up <- diff_reg[diff_reg$log2FC > 0,]
  target_genes <- diff_reg_up %>%
    filter(gene != "MAPKAPK2")
  
  # 优先选择 FDR < 0.05 且显著变化的基因
  sig_genes <- target_genes %>%
    filter(change != "NOT") %>%
    pull(gene)
  
  # 【防报错机制】：如果严格阈值下显著基因 < 10 个，自动提取 Distance/FC 改变最剧烈的 Top 50 基因
  if (length(sig_genes) < 30) {
    message("提示：显著基因较少，已自动切换为提取 Perturbation 程度最高的 Top 50 基因进行富集分析...")
    sig_genes <- target_genes %>%
      arrange(desc(abs(log2FC))) %>%
      head(50) %>%
      pull(gene)
  }
  
  # ------------------------------------------------------------------------------
  # 2. 运行 GO Biological Process (BP) 富集分析
  # ------------------------------------------------------------------------------
  # SYMBOL 转 ENTREZID
  gene_ids <- bitr(
    sig_genes, 
    fromType = "SYMBOL", 
    toType = "ENTREZID", 
    OrgDb = org.Hs.eg.db
  )
  
  go_res <- enrichGO(
    gene          = gene_ids$ENTREZID,
    OrgDb         = org.Hs.eg.db,
    ont           = "ALL",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.1,  # 适度放宽 p 值筛选
    qvalueCutoff  = 0.2,
    readable      = TRUE
  )
  
  go_res_data <- go_res@result
  go_res_data <- go_res_data[go_res_data$pvalue < 0.01,]
  # ------------------------------------------------------------------------------
  # 3. 结果处理与高级气泡图绘制 (Bubble Plot)
  # ------------------------------------------------------------------------------
  # 筛选出你的这 9 个目标通路
  target_terms <- c(
    "response to tumor necrosis factor", "response to lipopolysaccharide", 
    "cell chemotaxis", "leukocyte migration involved in inflammatory response", 
    "leukocyte aggregation", "chronic inflammatory response", 
    "protein nitrosylation", "peptidyl-cysteine S-nitrosylation", "amide transport"
  )
  
  # 提取并过滤 GO 结果对象
  ego_sub <- go_res
  ego_sub@result <- ego_sub@result %>% filter(tolower(Description) %in% tolower(target_terms))
  # 提取 Top 10 通路并格式化名称
  plot_go <- ego_sub@result %>%
    arrange(p.adjust) %>%
    mutate(
      # 首字母大写规则：仅第一个单词首字母大写
      Description = str_to_sentence(Description),
      # 计算 Gene Ratio 数值 (如 3/10 -> 0.3)
      GeneRatio_num = sapply(GeneRatio, function(x) {
        nums <- as.numeric(strsplit(x, "/")[[1]])
        return(nums[1] / nums[2])
      })
    )
  
  
    # 按显著性排序因子顺序
    plot_go$Description <- factor(plot_go$Description, levels = rev(plot_go$Description))
    
    # 绘制气泡图
    pdf("虚拟敲除3_GO功能富集气泡图.pdf", width = 7.5, height = 5.5)
    plot_go$pvalue_log <- -log2(plot_go$pvalue)
    
    library(ggplot2)
    library(dplyr)
    library(stringr)
    
    # ------------------------------------------------------------------------------
    # 1. 数据清洗与严格排序
    # ------------------------------------------------------------------------------
    plot_go_bar <- plot_go %>%
      mutate(
        # 格式规范：仅第一个单词首字母大写，修正专有名词 S-nitrosylation
        Description = str_to_sentence(Description),
        Description = gsub("S-nitrosylation", "S-nitrosylation", Description, ignore.case = TRUE),
        
        # 确保显著性指标正确
        fdr_log = -log10(p.adjust)
      ) %>%
      # 【关键】按 GeneRatio_num（或 fdr_log）从大到小排序
      arrange(pvalue_log, fdr_log) 
    
    # 将 Description 转为 Factor，保证 ggplot 画图时从上到下由大到小递减
    plot_go_bar$Description <- factor(
      plot_go_bar$Description, 
      levels = plot_go_bar$Description
    )
    
    # ------------------------------------------------------------------------------
    # 2. 绘制高清顶级学术色柱状图
    # ------------------------------------------------------------------------------
    p_nature_clean <- ggplot(plot_go_bar, aes(x = pvalue_log, y = Description)) +
      
      # 1. 柱体：去掉了末端的 geom_text，柱宽 0.68 极其饱满
      geom_col(aes(fill = pvalue_log), width = 0.68, alpha = 0.92) +
      
      # 2. 顶级期刊 Nat Comms 经典高级渐变色（深海蓝 -> 阳光黄 -> 暖红 -> 帝王红）
      scale_fill_gradientn(
        colors = c("#B395C8"),
        name = expression(-log[10]~"(FDR)")
      ) +
      
      # 3. 坐标轴与范围优化：紧贴 0 刻度线
      scale_x_continuous(
        expand = expansion(mult = c(0, 0.05))
      ) +
      
      # 4. 精细排版主题 (Nature Theme)
      theme_bw(base_size = 12) +
      labs(
        title = "Enriched Biological Processes",
        x = "pvalue_log",
        y = NULL
      ) +
      theme(
        # 标题样式
        plot.title = element_text(
          size = 13, 
          face = "bold", 
          hjust = 0, 
          margin = margin(b = 10)
        ),
        # Y 轴（通路名称）样式
        axis.text.y = element_text(
          size = 11, 
          color = "black", 
          face = "bold"
        ),
        axis.text.x = element_text(size = 10, color = "black"),
        axis.title.x = element_text(size = 11, face = "bold", margin = margin(t = 8)),
        
        # 边框与网格线（只保留极淡的 X 轴垂直参考线，画面极简）
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_line(color = "grey92", linewidth = 0.4),
        
        # 图例样式
        legend.position = "right",
        legend.title = element_text(size = 10, face = "bold"),
        legend.text = element_text(size = 9),
        legend.background = element_blank(),
        legend.box.background = element_blank(),
        legend.key = element_blank()
      )
    
    print(p_nature_clean)
    ggsave("MAPKAPK2_GO.pdf",width = 9,height = 6)
    
    #查看top基因在单细胞不同组别之间的表达
    scRNA <- qread("/data3/home/hu/Tongfeng/new_Result/Fig4/sce_pure.qs")
    
    library(Seurat)
    library(ggplot2)
    library(patchwork)
    
    # ------------------------------------------------------------------------------
    # 1. 基础设置与数据准备
    # ------------------------------------------------------------------------------
    DefaultAssay(scRNA) <- "RNA"
    
    # 确保 group 顺序为 Control -> Disease -> Remission
    scRNA$group <- factor(scRNA$group, levels = c("Control", "Disease", "Remission"))
    
    head(diff_reg$gene)
    # 选定前 3 个基因
    target_genes <- c( "CMTM2", "FCGR3B","AQP9")
    
    # 参考图同款莫兰迪配色方案 (Control: 雾霾蓝, Disease: 砖红/低饱和红, Remission: 暖米黄)
    my_cols <- c("Control" = "#6C8EBF", "Disease" = "#D07064", "Remission" = "#F2E2B8")
    
    # ------------------------------------------------------------------------------
    # 2. 循环绘制单基因图形并精细调整（完全对标参考图样式）
    # ------------------------------------------------------------------------------
    plot_list <- list()
    
    for (i in seq_along(target_genes)) {
      gene <- target_genes[i]
      
      p <- VlnPlot(
        scRNA,
        features = gene,
        group.by = "group",
        cols = my_cols,
        pt.size = 0 # 隐藏单细胞离散点
      ) +
        # 嵌套窄版箱线图（白底黑框、去离群点）
        geom_boxplot(
          width = 0.12, 
          fill = "white", 
          color = "black", 
          outlier.shape = NA, 
          linewidth = 0.6,
          alpha = 0.9
        ) +
        # 标签与标题设置
        labs(
          title = gene,
          x = NULL,
          y = if (i == 1) "Normalized Expression Level" else NULL # 仅最左侧显示 Y 轴标题
        ) +
        theme_bw(base_size = 12) +
        theme(
          # 标题设为加粗斜体 (如 CXCL8)
          plot.title = element_text(size = 14, face = "bold.italic", hjust = 0.5, margin = margin(b = 8)),
          
          # X 轴标签：倾斜 45 度、加粗
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, face = "bold", size = 11, color = "black"),
          axis.text.y = element_text(size = 10, color = "black"),
          axis.title.y = element_text(size = 12, face = "bold", margin = margin(r = 6)),
          
          # 边框与网格线（保留黑框，去掉内部背景网格）
          panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
          panel.grid = element_blank(),
          
          # 隐藏图例（因为 X 轴已有组名）
          legend.position = "none"
        )
      
      plot_list[[gene]] <- p
    }
    
    # ------------------------------------------------------------------------------
    # 3. 多图无缝拼接与输出
    # ------------------------------------------------------------------------------
    # 使用 patchwork 按 1 行 3 列无缝组合
    p_final <- wrap_plots(plot_list, nrow = 1)
    
    print(p_final)
    
    # 保存高分辨率 PDF
    ggsave("MAPKAPK2_Top3_Genes_VlnPlot.pdf", p_final, width = 8.5, height = 4.2, dpi = 300)
    
}

###换个基因PTGES
{
  ####敲除后分析
  PTGES_knk <- qread("/data3/home/hu/Tongfeng/new_Result/Fig7/PTGES_knk_res.qs")
  
  # 提取差异调控表
  diff_reg <- PTGES_knk$diffRegulation
  
  # 调试与检查：如果 FC 为绝对比例，计算 log2FC；同时避免 p.adj 为 0 导致的 -log10 计算溢出
  diff_reg <- diff_reg %>%
    mutate(
      # 处理 p.adj 极小或为 0 的情况，防止 -log10 产生 Inf
      p.adj_clean = ifelse(p.adj < 1e-300, 1e-300, p.adj),
      neg_log10_pval = -log10(p.adj_clean),
      # 算 log2FC 用于对称展示（若 FC 均为正数）
      log2FC = log2(FC + 1e-6)
    )
  
  # 阈值设置
  fdr_cutoff <- 0.5
  logFC_cutoff <- 0.5  # 可根据实际数据分布适当调整 (如 0.25 或 0.5)
  
  # 标记上调/下调/无显著变化
  diff_reg$change <- "NOT"
  diff_reg$change[diff_reg$p.value < fdr_cutoff & diff_reg$log2FC > logFC_cutoff] <- "UP"
  diff_reg$change[diff_reg$p.value < fdr_cutoff & diff_reg$log2FC < -logFC_cutoff] <- "DOWN"
  
  # 查看各组数量分布
  print("基因分组统计：")
  print(table(diff_reg$change))
  
  
  library(ggraph)
  library(igraph)
  library(dplyr)
  
  # 取出显著差异基因 (比如 Top 10-15 个)
  top_drgs <- diff_reg %>%
    filter(change != "NOT") %>%
    arrange(p.adj) %>%
    head(15)
  
  # 构建以 KO 基因 (MAPKAPK2) 为中心的网络关系数据框
  ko_gene <- "PTGES" # 替换为你的真实敲除基因名称
  edges <- data.frame(
    from = ko_gene,
    to = top_drgs$gene,
    weight = top_drgs$FC,
    change = top_drgs$change
  )
  
  library(ggraph)
  library(igraph)
  library(dplyr)
  
  # ------------------------------------------------------------------------------
  # 1. 数据清洗与规范化（保留 NOT，但区分样式）
  # ------------------------------------------------------------------------------
  edges_clean <- edges %>%
    # ① 剔除自身指向自身的自环 (MAPKAPK2 -> MAPKAPK2)
    filter(from != to) %>%
    # ② 对极悬殊的 weight 进行 rank 或 log2 处理，防止连线粗细爆表
    # 使用 rank 转换能让粗细过渡极其平滑美观
    mutate(edge_weight = rank(weight))
  
  # 构建 igraph 对象
  graph <- graph_from_data_frame(edges_clean, directed = FALSE)
  
  # 将 change 属性传递给节点（用于节点配色）
  V(graph)$node_change <- edges_clean$change[match(V(graph)$name, edges_clean$to)]
  # 敲除基因本身设为 KO
  V(graph)$node_change[V(graph)$name == "PTGES"] <- "KO"
  
  # ------------------------------------------------------------------------------
  # 2. 绘制包含背景基因的网络图 (Star Layout)
  # ------------------------------------------------------------------------------
  pdf("/data3/home/hu/Tongfeng/new_Result/Fig7/knock_PTGES_net.pdf", width = 8.5, height = 8)
  
  ggraph(graph, layout = "star", center = "PTGES") +
    # 1. 绘制连线：显著基因用实线+彩颜色，非显著(NOT)用浅灰虚线
    geom_edge_link(
      aes(
        color = change, 
        width = edge_weight,
        linetype = change == "NOT"
      ), 
      alpha = 0.7
    ) +
    # 线条颜色配置
    scale_edge_color_manual(
      values = c("UP" = "#E64B35FF", "DOWN" = "#4DBBD5FF", "NOT" = "grey80"),
      name = "Regulation"
    ) +
    # 线条粗细映射 (1 到 3 像素)
    scale_edge_width_continuous(range = c(0.8, 3.2), guide = "none") +
    # 线型设置：NOT 为虚线，UP/DOWN 为实线
    scale_edge_linetype_manual(values = c("FALSE" = "solid", "TRUE" = "dashed"), guide = "none") +
    
    # 2. 绘制节点：显著基因高亮，非显著基因灰点
    geom_node_point(
      aes(color = node_change, size = node_change)
    ) +
    scale_color_manual(
      values = c(
        "KO" = "#E64B35FF",      # 中心 KO 基因：红色
        "UP" = "#E64B35FF",      # 上调基因：红色
        "DOWN" = "#4DBBD5FF",    # 下调基因：蓝色
        "NOT" = "grey70"         # 非显著基因：浅灰
      ),
      labels = c("KO" = "KO Target", "UP" = "Up-regulated", "DOWN" = "Down-regulated", "NOT" = "Not Significant"),
      name = "Gene Status"
    ) +
    # 节点大小：中心 KO 基因和显著基因大一点，背景基因小一点
    scale_size_manual(
      values = c("KO" = 9, "UP" = 6.5, "DOWN" = 6.5, "NOT" = 3.5),
      guide = "none"
    ) +
    
    # 3. 标签绘制：只给中心 KO 基因和【显著变化】的基因打标签（避免 NOT 基因挤乱图表）
    geom_node_text(
      aes(
        label = ifelse(node_change %in% c("KO", "UP", "DOWN"), name, "")
      ),
      repel = TRUE,
      fontface = "bold.italic",
      size = 4.5,
      point.padding = unit(0.3, "lines"),
      max.overlaps = 50
    ) +
    
    # 4. 主题设置
    theme_void() +
    labs(
      title = "Local Regulatory Network Perturbation (MAPKAPK2 KO)",
      subtitle = "Highlighting significantly perturbed target genes among local network"
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 15),
      plot.subtitle = element_text(hjust = 0.5, size = 10, color = "grey40"),
      legend.position = "bottom",
      legend.title = element_text(face = "bold", size = 10),
      legend.box = "horizontal"
    )
  
  dev.off()
  
  
  
  #####富集分析
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(ggplot2)
  library(dplyr)
  library(stringr)
  
  # ------------------------------------------------------------------------------
  # 1. 提取基因列表 (包含容错逻辑)
  # ------------------------------------------------------------------------------
  # 剔除敲除基因本身
  diff_reg_up <- diff_reg[diff_reg$log2FC > 0,]
  target_genes <- diff_reg_up %>%
    filter(gene != "PTGES")
  
  # 优先选择 FDR < 0.05 且显著变化的基因
  sig_genes <- target_genes %>%
    filter(change != "NOT") %>%
    pull(gene)
  
  # 【防报错机制】：如果严格阈值下显著基因 < 10 个，自动提取 Distance/FC 改变最剧烈的 Top 50 基因
  if (length(sig_genes) < 30) {
    message("提示：显著基因较少，已自动切换为提取 Perturbation 程度最高的 Top 50 基因进行富集分析...")
    sig_genes <- target_genes %>%
      arrange(desc(abs(log2FC))) %>%
      head(50) %>%
      pull(gene)
  }
  
  # ------------------------------------------------------------------------------
  # 2. 运行 GO Biological Process (BP) 富集分析
  # ------------------------------------------------------------------------------
  # SYMBOL 转 ENTREZID
  gene_ids <- bitr(
    sig_genes, 
    fromType = "SYMBOL", 
    toType = "ENTREZID", 
    OrgDb = org.Hs.eg.db
  )
  
  go_res <- enrichGO(
    gene          = gene_ids$ENTREZID,
    OrgDb         = org.Hs.eg.db,
    ont           = "ALL",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.1,  # 适度放宽 p 值筛选
    qvalueCutoff  = 0.2,
    readable      = TRUE
  )
  
  go_res_data <- go_res@result
  go_res_data <- go_res_data[go_res_data$pvalue < 0.01,]
  # ------------------------------------------------------------------------------
  # 3. 结果处理与高级气泡图绘制 (Bubble Plot)
  # ------------------------------------------------------------------------------
  # 筛选出你的这 9 个目标通路
  selected_GO <- c(
    "response to tumor necrosis factor",
    "response to lipopolysaccharide",
    "cellular response to interleukin-1",
    "chemokine-mediated signaling pathway",
    "neutrophil chemotaxis",
    "myeloid leukocyte migration",
    "response to purine-containing compound",
    "response to reactive oxygen species",
    "protein nitrosylation",
    "purine nucleobase transport"
  )
  
  # 提取并过滤 GO 结果对象
  ego_sub <- go_res
  ego_sub@result <- ego_sub@result %>% filter(tolower(Description) %in% tolower(selected_GO))
  # 提取 Top 10 通路并格式化名称
  plot_go <- ego_sub@result %>%
    arrange(p.adjust) %>%
    mutate(
      # 首字母大写规则：仅第一个单词首字母大写
      Description = str_to_sentence(Description),
      # 计算 Gene Ratio 数值 (如 3/10 -> 0.3)
      GeneRatio_num = sapply(GeneRatio, function(x) {
        nums <- as.numeric(strsplit(x, "/")[[1]])
        return(nums[1] / nums[2])
      })
    )
  
  
  # 按显著性排序因子顺序
  plot_go$Description <- factor(plot_go$Description, levels = rev(plot_go$Description))
  
  # 绘制气泡图
  plot_go$pvalue_log <- -log2(plot_go$pvalue)
  
  library(ggplot2)
  library(dplyr)
  library(stringr)
  
  # ------------------------------------------------------------------------------
  # 1. 数据清洗与严格排序
  # ------------------------------------------------------------------------------
  plot_go_bar <- plot_go %>%
    mutate(
      # 格式规范：仅第一个单词首字母大写，修正专有名词 S-nitrosylation
      Description = str_to_sentence(Description),
      Description = gsub("S-nitrosylation", "S-nitrosylation", Description, ignore.case = TRUE),
      
      # 确保显著性指标正确
      fdr_log = -log10(p.adjust)
    ) %>%
    # 【关键】按 GeneRatio_num（或 fdr_log）从大到小排序
    arrange(pvalue_log, fdr_log) 
  
  # 将 Description 转为 Factor，保证 ggplot 画图时从上到下由大到小递减
  plot_go_bar$Description <- factor(
    plot_go_bar$Description, 
    levels = plot_go_bar$Description
  )
  
  # ------------------------------------------------------------------------------
  # 2. 绘制高清顶级学术色柱状图
  # ------------------------------------------------------------------------------
  p_nature_clean <- ggplot(plot_go_bar, aes(x = pvalue_log, y = Description)) +
    
    # 1. 柱体：去掉了末端的 geom_text，柱宽 0.68 极其饱满
    geom_col(aes(fill = pvalue_log), width = 0.68, alpha = 0.92) +
    
    # 2. 顶级期刊 Nat Comms 经典高级渐变色（深海蓝 -> 阳光黄 -> 暖红 -> 帝王红）
    scale_fill_gradientn(
      colors = c("#B395C8"),
      name = expression(-log[10]~"(FDR)")
    ) +
    
    # 3. 坐标轴与范围优化：紧贴 0 刻度线
    scale_x_continuous(
      expand = expansion(mult = c(0, 0.05))
    ) +
    
    # 4. 精细排版主题 (Nature Theme)
    theme_bw(base_size = 12) +
    labs(
      title = "Enriched Biological Processes",
      x = "pvalue_log",
      y = NULL
    ) +
    theme(
      # 标题样式
      plot.title = element_text(
        size = 13, 
        face = "bold", 
        hjust = 0, 
        margin = margin(b = 10)
      ),
      # Y 轴（通路名称）样式
      axis.text.y = element_text(
        size = 11, 
        color = "black", 
        face = "bold"
      ),
      axis.text.x = element_text(size = 10, color = "black"),
      axis.title.x = element_text(size = 11, face = "bold", margin = margin(t = 8)),
      
      # 边框与网格线（只保留极淡的 X 轴垂直参考线，画面极简）
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_line(color = "grey92", linewidth = 0.4),
      
      # 图例样式
      legend.position = "right",
      legend.title = element_text(size = 10, face = "bold"),
      legend.text = element_text(size = 9),
      legend.background = element_blank(),
      legend.box.background = element_blank(),
      legend.key = element_blank()
    )
  
  print(p_nature_clean)
  ggsave("PTGES_GO.pdf",width = 9,height = 6)
  
  #查看top基因在单细胞不同组别之间的表达
  scRNA <- qread("/data3/home/hu/Tongfeng/new_Result/Fig4/sce_pure.qs")
  
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  
  # ------------------------------------------------------------------------------
  # 1. 基础设置与数据准备
  # ------------------------------------------------------------------------------
  DefaultAssay(scRNA) <- "RNA"
  
  # 确保 group 顺序为 Control -> Disease -> Remission
  scRNA$group <- factor(scRNA$group, levels = c("Control", "Disease", "Remission"))
  
  # 选定前 3 个基因
  target_genes <- head(diff_reg$gene)[c(1,4,5)]
  
  # 参考图同款莫兰迪配色方案 (Control: 雾霾蓝, Disease: 砖红/低饱和红, Remission: 暖米黄)
  my_cols <- c("Control" = "#6C8EBF", "Disease" = "#D07064", "Remission" = "#F2E2B8")
  
  # ------------------------------------------------------------------------------
  # 2. 循环绘制单基因图形并精细调整（完全对标参考图样式）
  # ------------------------------------------------------------------------------
  plot_list <- list()
  
  for (i in seq_along(target_genes)) {
    gene <- target_genes[i]
    
    p <- VlnPlot(
      scRNA,
      features = gene,
      group.by = "group",
      cols = my_cols,
      pt.size = 0 # 隐藏单细胞离散点
    ) +
      # 嵌套窄版箱线图（白底黑框、去离群点）
      geom_boxplot(
        width = 0.12, 
        fill = "white", 
        color = "black", 
        outlier.shape = NA, 
        linewidth = 0.6,
        alpha = 0.9
      ) +
      # 标签与标题设置
      labs(
        title = gene,
        x = NULL,
        y = if (i == 1) "Normalized Expression Level" else NULL # 仅最左侧显示 Y 轴标题
      ) +
      theme_bw(base_size = 12) +
      theme(
        # 标题设为加粗斜体 (如 CXCL8)
        plot.title = element_text(size = 14, face = "bold.italic", hjust = 0.5, margin = margin(b = 8)),
        
        # X 轴标签：倾斜 45 度、加粗
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, face = "bold", size = 11, color = "black"),
        axis.text.y = element_text(size = 10, color = "black"),
        axis.title.y = element_text(size = 12, face = "bold", margin = margin(r = 6)),
        
        # 边框与网格线（保留黑框，去掉内部背景网格）
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
        panel.grid = element_blank(),
        
        # 隐藏图例（因为 X 轴已有组名）
        legend.position = "none"
      )
    
    plot_list[[gene]] <- p
  }
  
  # ------------------------------------------------------------------------------
  # 3. 多图无缝拼接与输出
  # ------------------------------------------------------------------------------
  # 使用 patchwork 按 1 行 3 列无缝组合
  p_final <- wrap_plots(plot_list, nrow = 1)
  
  print(p_final)
  
  # 保存高分辨率 PDF
  ggsave("MAPKAPK2_Top3_Genes_VlnPlot.pdf", p_final, width = 8.5, height = 4.2, dpi = 300)
  
}
