library(readr)
library(tidyverse)
clinical=read_tsv("C:\\Users\\lenovo\\Desktop\\Survival analysis\\clinical.project-tcga-brca.2024-03-11\\clinical.tsv")

############################提取基本的临床信息
#提取样品ID
ID=clinical$case_submitter_id
#提取年龄
age=clinical$age_at_index
#提取性别
gender=clinical$gender
#提取生存时间
time=clinical$days_to_death
#提取生存状态
status=clinical$vital_status
#提取TMN分期
pathologicT=clinical$ajcc_pathologic_t
pathologicM=clinical$ajcc_pathologic_m
pathologicN=clinical$ajcc_pathologic_n
#提取stage分期
pathologicStage=clinical$ajcc_pathologic_stage

############################合并临床信息
#合并信息
TCGA_merge=cbind(ID,
                 age,
                 gender,
                 time,
                 status,
                 pathologicT,
                 pathologicM,
                 pathologicN,
                 pathologicStage)

############################删除缺失值并导出
#删除缺失值
TCGA_merge[which(TCGA_merge=="'--")]=NA
TCGA_clinical=na.omit(TCGA_merge)
TCGA_clinical=as.data.frame(TCGA_clinical)
#删除重复ID
duplicated(TCGA_clinical$ID)
TCGA_clinical<-TCGA_clinical[!duplicated(TCGA_clinical$ID),]
duplicated(TCGA_clinical)
#导出文件
rownames(TCGA_clinical)=TCGA_clinical$ID
TCGA_clinical=TCGA_clinical[,2:ncol(TCGA_clinical)]
write.csv(TCGA_clinical,file = "TCGA_clinical.csv",quote = F)

############################尝试生存分析
rm(list = ls())
#install.packages(c("survival", "survminer"))
library("survival")
library("survminer")
TCGA_clinical <- read.csv("C:\\Users\\lenovo\\Desktop\\Survival analysis\\TCGA_clinical.csv")
library(dplyr)
TCGA_clinical <- rename(TCGA_clinical, sex = gender)
# 并且已经提取了生存时间（time）、生存状态（status）和性别（sex）
TCGA_clinical$time <- as.numeric(TCGA_clinical$time) #将时间转化为数值类型
time <- TCGA_clinical$time
TCGA_clinical$status <- ifelse(TCGA_clinical$status == "Alive", FALSE , TRUE)  #将生存转化为TRUE, FALSE
status <- TCGA_clinical$status
sex <- TCGA_clinical$gender
# 构建Surv对象
surv_object <- Surv(time = TCGA_clinical$time, event = TCGA_clinical$status)
# 创建一个cox比例风险回归模型
library(survival)
cox_model <- coxph(Surv(time, status) ~ sex, data = TCGA_clinical)

# 或者进行生存曲线拟合
library(survminer)
ggsurvplot(survfit(Surv(time, status) ~ sex, data = TCGA_clinical))


############################################下载基因数据
if (!requireNamespace("BiocManager", quietly = TRUE))
install.packages("BiocManager")
library(BiocManager)

#BiocManager::install("TCGAbiolinks")
library(TCGAbiolinks)

#install.packages("survminer")
library(survminer)

#BiocManager::install("SummarizedExperiment")
library(SummarizedExperiment)

project="TCGA-BRCA"

#查看有哪些分类数据
TCGAbiolinks:::getProjectSummary("TCGA-BRCA")

#定义查询数据为基因表达数据
data.category = "Gene expression"

#指定了基因表达量化数据
data.type = "Gene expression quantification"

#工作流类型为 HTSeq 计数
workflow.type="HTSeq - Counts"

#测序平台为 Illumina HiSeq
platform = "Illumina HiSeq"

file.type = "results"

#列出了具体的样本 ID
listSamples <- c("TCGA-E9-A1NG-11A-52R-A14M-07","TCGA-BH-A1FC-11A-32R-A13Q-07",
                 "TCGA-A7-A13G-11A-51R-A13Q-07","TCGA-BH-A0DK-11A-13R-A089-07",
                 "TCGA-E9-A1RH-11A-34R-A169-07","TCGA-BH-A0AU-01A-11R-A12P-07",
                 "TCGA-C8-A1HJ-01A-11R-A13Q-07","TCGA-A7-A13D-01A-13R-A12P-07",
                 "TCGA-A2-A0CV-01A-31R-A115-07","TCGA-AQ-A0Y5-01A-11R-A14M-07")
query <- GDCquery(project = "TCGA-BRCA",
                  data.category = "Transcriptome Profiling",
                  platform = ,
                  data.type = "Gene Expression Quantification",
                  experimental.strategy = "RNA-Seq",
                  #file.type = "results",
                  barcode = listSamples,
                  #legacy = TRUE,
                  access="open",
                  workflow.type="STAR - Counts")

GDCdownload(
  query,
  token.file,
  method = "api",
  directory = "GDCdata",
  files.per.chunk = NULL
)


############################################读入分析
#GDC函数下载的所有文件在一个一个的文件夹中，非常不便，所以在读入之前要将其放到同一个文件夹下
dir.create('TCGA-BRCA')

# 获取子文件夹列表
subfolders <- dir("GDCdata\\TCGA-BRCA\\Transcriptome_Profiling\\Gene_Expression_Quantification\\.", full.names = TRUE)

# 遍历每个子文件夹
for (dirname in subfolders) {
  # 使用 list.files 函数找到子文件夹下的 .tsv 文件
  files <- list.files(dirname, pattern = "*.tsv", full.names = TRUE)
  
  # 复制每个 .tsv 文件到 TCGA-BRCA 文件夹
  for (file in files) {
    dest_path <- file.path("TCGA-BRCA", basename(file))
    file.copy(file, dest_path)
  }
}

##开始读入
#install.packages("readr")
library(readr)
dir("C:\\Users\\lenovo\\Desktop\\Survival analysis\\TCGA-BRCA",full.names = T)
fs=list.files("C:\\Users\\lenovo\\Desktop\\Survival analysis\\TCGA-BRCA")
column_names <- c("gene_id",
                  "gene_name",
                  "gene_type",
                  "unstranded",
                  "stranded_first",
                  "stranded_second",
                  "tpm_unstranded",
                  "fpkm_unstranded",
                  "fpkm_uq_unstranded"
)
TGCA_files_list=lapply(fs, function(x){
  a=read_tsv(file.path("C:\\Users\\lenovo\\Desktop\\Survival analysis\\TCGA-BRCA",x),
             skip = 6,
             col_names = column_names)
})
class(TGCA_files_list)

###合并读取结果
TGCA_files_count=do.call(rbind,TGCA_files_list)
dim(TGCA_files_count)
###do.call可以批量合并列表里的元素，cbind函数是对数据框按列合并，前提是行的数量相同，cbind在这是do.call函数里的一个参数。
gene=rbind(TGCA_files_list[[1]],
           TGCA_files_list[[2]],
           TGCA_files_list[[3]],
           TGCA_files_list[[4]],
           TGCA_files_list[[5]],
           TGCA_files_list[[6]],
           TGCA_files_list[[4]])

#读取json文件（fromJSON函数来自R自带R包jsonlite）
#TCGA的JSON文件通常用于存储元数据（metadata），即关于数据的数据。
#这些元数据提供了有关TCGA数据集中各个文件的重要信息，帮助研究人员理解数据的来源、类型以及与之关联的临床信息等。
mtdt <- jsonlite::fromJSON("biospecimen.project-tcga-brca.2024-03-14.json")
TCGAids <- mtdt$submitter_id
TCGAids[[1]]#查看submitter_id中的第一行，即第一个样本

library(dplyr)
gene_clinical <- inner_join(mtdt, TCGA_clinical, by = c("submitter_id" = "X"))

library(purrr)
library(dplyr)
#若有必要，先确保 submitter_id 列为字符类型
gene_clinical$samples <-
  map(gene_clinical$samples, ~ .x %>%
        mutate(submitter_id = as.character(submitter_id)) %>%
        mutate(submitter_id = substr(submitter_id, 1, nchar(submitter_id) - 1)))

######提取全基因全样本表达矩阵
options(stringsAsFactors = F)
library(WGCNA)

######提取TNBC的表达矩阵
counts <- read.table("C:\\Users\\lenovo\\Desktop\\Survival analysis\\TNBC_expr_log_ALL.txt", header = T, check.names = F)
# exprData <- log2(edgeR::cpm(counts+1))

# 先计算CPM
cpm_counts <- edgeR::cpm(counts + 1)
# 接着进行z-score标准化
exprData <- scale(cpm_counts)

# 这样就在每个数据框的第一行添加了一个新的变量 "clinical_id"，其值为原数据框的列名。

#####################################################从表达矩阵中查找目的基因
# 确保这两行存在于 exprData 中
# 此处选择了三个基因同时绘制
if ("PGLS" %in% rownames(exprData) && "PFKFB1" %in% rownames(exprData) 
    && "CBR4" %in% rownames(exprData)) {
  
  # 提取并转置行数据
  PGLS_values <- exprData[rownames(exprData) == "PGLS", ]
  PFKFB1_values <- exprData[rownames(exprData) == "PFKFB1", ]
  CBR4_values <- exprData[rownames(exprData) == "CBR4", ]
  
  # 创建新数据框，列名是基因名，行名是原表达矩阵的列名
  SP <- data.frame(PGLS = as.vector(PGLS_values),
                   PFKFB1 = as.vector(PFKFB1_values),
                   CBR4 = as.vector(CBR4_values))
  
  # 设置行名
  rownames(SP) <- colnames(exprData)
  
} else {
  print("Either 'PGLS','CBR4' or 'PFKFB1' does not exist in the expression matrix.")
}

# 未返回else内容则证明这三个基因均在exprData中

#####################################################将SP中的基因表达数据与临床数据进行匹配
library(dplyr)
library(tidyr)

# 将 SP 数据框的行名添加为新列
SP_with_submitter_id <- cbind(submitter_id = rownames(SP), SP)
rownames(SP_with_submitter_id) <- NULL

# 初始化空的数据框用于最后合并
gene_clinical_SP <- data.frame()

# 遍历 gene_clinical 的所有子数据框并执行匹配合并操作
for (i in seq_along(gene_clinical[[4]])) {
  # 获取当前子数据框
  current_df <- gene_clinical[[4]][[i]]

  # 只对包含 submitter_id 列的子数据框执行操作
  if ("submitter_id" %in% names(current_df)) {
    # 执行匹配合并操作并添加至最终的合并数据框中
    merged_data <- left_join(current_df, SP_with_submitter_id, by = "submitter_id")

    # 合并到最终结果中
    gene_clinical_SP <- bind_rows(gene_clinical_SP, merged_data)
  }
}

# 最终得到的是名为gene_clinical_SP的数据框，
# 其中包含了与SP匹配合并的所有行，并且不包含那些不存在于子数据框或与SP不匹配的行。

gene_clinical_SP <- gene_clinical_SP[!is.na(gene_clinical_SP$PGLS), ]

###############################################将临床数据与基因数据合并
# 确保submitter_id列都是字符型
gene_clinical_SP$submitter_id <- as.character(gene_clinical_SP$submitter_id)
gene_clinical$submitter_id <- as.character(gene_clinical$submitter_id)

# 提取两个数据框中submitter_id前10位字符
gene_clinical_SP$id_match <- substr(gene_clinical_SP$submitter_id, 1, 10)
gene_clinical$id_match <- substr(gene_clinical$submitter_id, 1, 10)

# 进行匹配合并
Gene_clinical_SP <- merge(gene_clinical_SP, gene_clinical, by = "id_match", all = FALSE)

# 重命名gene_clinical中的submitter_id列
names(Gene_clinical_SP)[names(Gene_clinical_SP) == "submitter_id.y"] <- "submitter_id2"

# 删除临时的id_match列
Gene_clinical_SP$id_match <- NULL

# 由于我们设置了all = FALSE，所以未进行匹配操作的行已经被自动删除

################将基因表达数据打上标签
# 对基因按照表达量分为High、Low和NS三组，以±0.5为界限
Gene_clinical_SP$PGLS <- ifelse(Gene_clinical_SP$PGLS > 0.5, "High", 
                                  ifelse(Gene_clinical_SP$PGLS < -0.5, "Low", "NS"))

Gene_clinical_SP$PFKFB1 <- ifelse(Gene_clinical_SP$PFKFB1 > 0.5, "High", 
                                   ifelse(Gene_clinical_SP$PFKFB1 < -0.5, "Low", "NS"))

Gene_clinical_SP$CBR4 <- ifelse(Gene_clinical_SP$CBR4 > 0.5, "High", 
                                   ifelse(Gene_clinical_SP$CBR4 < -0.5, "Low","NS"))
# 假设 PGLS 已经被转换成了字符向量，现在将其转化为有序因子
Gene_clinical_SP$PGLS <- factor(Gene_clinical_SP$PGLS, levels = c("NS", "Low", "High"), ordered = TRUE)
Gene_clinical_SP$PFKFB1 <- factor(Gene_clinical_SP$PFKFB1, levels = c("NS", "Low", "High"), ordered = TRUE)
Gene_clinical_SP$CBR4 <- factor(Gene_clinical_SP$CBR4, levels = c( "NS","Low", "High"), ordered = TRUE)
#Gene_clinical_SP$CBR4 <- replace(Gene_clinical_SP$CBR4, Gene_clinical_SP$CBR4 == "NS", "High")
Gene_clinical_SP$CBR4[Gene_clinical_SP$CBR4 == "NS"] <- "Low"

##########在得到Gene_clinical_SP注意在里面查看所选择的基因的分组信息具体是什么

##################尝试绘制生存曲线图
library(survival)
library(survminer)
library(dplyr)
Gene_clinical_SP$time <- as.numeric(Gene_clinical_SP$time) #将时间转化为数值类型
time <- Gene_clinical_SP$time
#将生存转化为TRUE, FALSE
Gene_clinical_SP$status <- ifelse(Gene_clinical_SP$status == "Alive", FALSE , TRUE)  
status <- Gene_clinical_SP$status
PGLS <- Gene_clinical_SP$PGLS
PFKFB1 <- Gene_clinical_SP$PFKFB1
CBR4 <- Gene_clinical_SP$CBR4
# 构建Surv对象
surv_object <- Surv(time = Gene_clinical_SP$time, event = Gene_clinical_SP$status)
# 创建一个cox比例风险回归模型
library(survival)
cox_modelF <- coxph(Surv(time, status) ~ PGLS, data = Gene_clinical_SP)
cox_modelS <- coxph(Surv(time, status) ~ PFKFB1, data = Gene_clinical_SP)
cox_modelS <- coxph(Surv(time, status) ~ CBR4, data = Gene_clinical_SP)

summary(cox_modelF) # 查看模型摘要以了解统计显著性和HR值
summary(cox_modelS)

library(survminer)
F=ggsurvplot(survfit(Surv(time, status) ~ PGLS, data = Gene_clinical_SP))
S=ggsurvplot(survfit(Surv(time, status) ~ PFKFB1, data = Gene_clinical_SP))
A=ggsurvplot(survfit(Surv(time, status) ~ CBR4, data = Gene_clinical_SP))
print(F)
print(S)
print(A)

######经过上述步骤，一个普通的生存图就绘制好了。但是还不够美观，因此接下来对其进行美化


####这里运行F_plot相关代码块之后会报错。这是一个很经典的错误。
###由于PGLS表达量与生存的相关性不显著，因此在分组的时候只被分为了一组，因此生存曲线只有一条线。
##如果一定要对这个基因作生存曲线的话就需要修改第297行代码。

print(F)
# 确保PGLS被正确地编码为因子
Gene_clinical_SP$PGLS <- factor(Gene_clinical_SP$PGLS, levels = c("High", "NS"))
# 使用survfit函数计算生存函数
surv_fit_PGLS <- survfit(Surv(time, status) ~ PGLS, data = Gene_clinical_SP)
F_plot <- ggsurvplot(surv_fit_PGLS,
                     pval = TRUE,  # 显示P值
                     conf.int = TRUE,  # 显示置信区间
                     risk.table = TRUE,  # 显示风险表
                     legend.labs = levels(Gene_clinical_SP$PGLS),  # 设置图例标签
                     palette = c("#E69F00", "#56B4E9"),  # 自定义颜色
                     ggtheme = theme_minimal(),  # 使用简洁主题
                     title = "Effect of PGLS status on time to live",  # 添加主标题
                     xlab = "Follow-High time (days)",  # X轴标签
                     ylab = "Survival function estimates",  # Y轴标签
                     font.main = 18,  # 主标题字体大小
                     font.subtitle = 14,  # 副标题字体大小
                     censor.shape = 15,  # 修改删失数据点的形状
                     censor.size = 3,  # 修改删失数据点的大小
                     surv.scale = "percent")  # 按百分比比例显示生存率
# 更新标题和图例的位置
F_plot$plot <- F_plot$plot +
  theme(
    plot.title = element_text(hjust = 0.5, size = 18), # 主标题居中
    legend.position = c(0.9, 0.9), # 根据旧版本的警告信息，这里可能需要调整为内部布局参数
    legend.justification = c("right", "top")
  )
print(F_plot)


#####如果一定要画，以下是F_plot正确画法。因为只有一个分组，因此将分组信息改为一个就好
print(F)
# 确保PGLS被正确地编码为因子
Gene_clinical_SP$PGLS <- factor(Gene_clinical_SP$PGLS, levels = c("High"))
# 使用survfit函数计算生存函数
surv_fit_PGLS <- survfit(Surv(time, status) ~ PGLS, data = Gene_clinical_SP)
F_plot <- ggsurvplot(surv_fit_PGLS,
                     pval = TRUE,  # 显示P值
                     conf.int = TRUE,  # 显示置信区间
                     risk.table = TRUE,  # 显示风险表
                     legend.labs = levels(Gene_clinical_SP$PGLS),  # 设置图例标签
                     palette = c("#56B4E9"),  # 自定义颜色
                     ggtheme = theme_minimal(),  # 使用简洁主题
                     title = "Effect of PGLS status on time to live",  # 添加主标题
                     xlab = "Follow-High time (days)",  # X轴标签
                     ylab = "Survival function estimates",  # Y轴标签
                     font.main = 18,  # 主标题字体大小
                     font.subtitle = 14,  # 副标题字体大小
                     censor.shape = 15,  # 修改删失数据点的形状
                     censor.size = 3,  # 修改删失数据点的大小
                     surv.scale = "percent")  # 按百分比比例显示生存率
# 更新标题和图例的位置
F_plot$plot <- F_plot$plot +
  theme(
    plot.title = element_text(hjust = 0.5, size = 18), # 主标题居中
    legend.position = c(0.9, 0.9), # 根据旧版本的警告信息，这里可能需要调整为内部布局参数
    legend.justification = c("right", "top")
  )
print(F_plot)



#############

print(S)
Gene_clinical_SP$PFKFB1 <- factor(Gene_clinical_SP$PFKFB1, levels = c("NS", "Low"))
# 使用survfit函数计算生存函数
surv_fit_PFKFB1 <- survfit(Surv(time, status) ~ PFKFB1, data = Gene_clinical_SP)
S_plot <- ggsurvplot(surv_fit_PFKFB1,
                     pval = TRUE,  # 显示P值
                     conf.int = TRUE,  # 显示置信区间
                     risk.table = TRUE,  # 显示风险表
                     legend.labs = levels(Gene_clinical_SP$PFKFB1),  # 设置图例标签
                     palette="lancet",  # 使用柳叶刀自动配色
                     ggtheme = theme_minimal(),  # 使用简洁主题
                     title = "Effect of PFKFB1 status on time to live",  # 添加主标题
                     xlab = "Follow-High time (days)",  # X轴标签
                     ylab = "Survival function estimates",  # Y轴标签
                     font.main = 18,  # 主标题字体大小
                     font.subtitle = 14,  # 副标题字体大小
                     censor.shape = 15,  # 修改删失数据点的形状
                     censor.size = 3,  # 修改删失数据点的大小
                     surv.scale = "percent")  # 按百分比比例显示生存率
print(S_plot)




# 使用survfit函数计算生存函数
surv_fit_CBR4 <- survfit(Surv(time, status) ~ CBR4, data = Gene_clinical_SP)

library(survminer) # 引入survminer包以便于绘制和美化生存曲线

print(A)
Gene_clinical_SP$CBR4
Gene_clinical_SP$CBR4 <- factor(Gene_clinical_SP$CBR4, levels = c("High","Low"))
# 绘制并美化生存曲线
A_plot <- ggsurvplot(surv_fit_CBR4,
                     pval = TRUE,  # 显示P值
                     conf.int = TRUE,  # 显示置信区间
                     risk.table = TRUE,  # 显示风险表
                     legend.labs = levels(Gene_clinical_SP$CBR4),  # 设置图例标签
                     palette = c("#FF0000","#56B4E9"),  # 自定义颜色，需要RGB颜色代码
                     ggtheme = theme_minimal(),  # 使用简洁主题
                     title = "Effect of CBR4 status on time to live",  # 添加主标题
                     xlab = "Follow-High time (days)",  # X轴标签
                     ylab = "Survival function estimates",  # Y轴标签
                     font.main = 18,  # 主标题字体大小
                     font.subtitle = 14,  # 副标题字体大小
                     censor.shape = 15,  # 修改删失数据点的形状
                     censor.size = 3,  # 修改删失数据点的大小
                     surv.scale = "percent")  # 按百分比比例显示生存率

# 更新标题和图例的位置
A_plot$plot <- A_plot$plot +
  theme(
    plot.title = element_text(hjust = 0.5, size = 18), # 主标题居中
    legend.position = c(0.9, 0.9), # 根据旧版本的警告信息，这里可能需要调整为内部布局参数
    legend.justification = c("right", "top")
  )

print(A_plot)









