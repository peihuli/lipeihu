############PCA+sva对于data1和data2

####
#install.packages("limma")
options(BioC_mirror="https://mirrors.tuna.tsinghua.edu.cn/bioconductor")
local({r <- getOption("repos")  
r["CRAN"] <- "http://mirrors.tuna.tsinghua.edu.cn/CRAN/"   
options(repos=r)}) 
#BiocManager::install("GEOquery")
#BiocManager::install("bladderbatch")
library(sva)
library(genefilter)
library(bladderbatch)
library('Biobase')
library('GEOquery')

###### 加载需要用到的R包
library(psych)
library(reshape2)
library(ggplot2)
library(factoextra)
library(sva)

AML_Sample_info<-read.delim("E:\\yan\\APL\\data\\0rawdata\\AML_Sample_info_sur.txt",stringsAsFactors = F)

GSE122505_Dataset_1<-read.delim("E:\\yan\\APL\\data\\0rawdata\\1.0-数据预处理\\GSE122505_Dataset_1.txt",
                                check.names = F,
                                stringsAsFactors = F)
GSE122511_Dataset_2<-read.delim("E:\\yan\\APL\\data\\0rawdata\\1.0-数据预处理\\GSE122511_Dataset_2.txt",
                                check.names = F,
                                stringsAsFactors = F)

length(intersect(rownames(GSE122505_Dataset_1),rownames(GSE122511_Dataset_2)))

data1_2_inter_gene<-intersect(rownames(GSE122505_Dataset_1),rownames(GSE122511_Dataset_2))

data_1_2<-cbind(GSE122505_Dataset_1[match(data1_2_inter_gene,rownames(GSE122505_Dataset_1)),],
                GSE122511_Dataset_2[match(data1_2_inter_gene,rownames(GSE122511_Dataset_2)),])

gsm_type<-c(paste("data1",AML_Sample_info$Condition[match(colnames(GSE122505_Dataset_1),AML_Sample_info$GSM)],sep = "_"),
            paste("data2",AML_Sample_info$Condition[match(colnames(GSE122511_Dataset_2),AML_Sample_info$GSM)],sep = "_"))

data_1_2[1:5,1:5]
batch_info<-data.frame(GSM=colnames(data_1_2),
                       batch=c(rep(1,ncol(GSE122505_Dataset_1)),rep(2,ncol(GSE122511_Dataset_2))),
                       Cancer=AML_Sample_info$Condition[match(colnames(data_1_2),AML_Sample_info$GSM)]
                       )
table(batch_info$Cancer)
table(batch_info$batch)





######批处理前PCA分析及可视化
data_1_2.pca<-prcomp(t(data_1_2[1:5000,]),
                     scale=T,
                     rank=10,###要使用的主成分的最大数目。
                     retx=T,##是否应返回已旋转的变量。
                     )

pdf("E:\\yan\\APL\\sva_batch\\picture\\before_batch.pdf",height = 8,width = 10)
fviz_pca_ind(data_1_2.pca, 
             label="none", #隐藏个人标签
             addEllipses = TRUE, ####浓度椭圆
             #ellipse.type="norm",
             ellipse.level=0.9,##椭圆大小
             habillage = gsm_type,
             palette = c("#00afbb","#e7b800","#fc4e07","#FF9966"),###填充颜色
             #mean.point=F
             )
dev.off()


#######使用sva包进行批次处理
batch_info$hasCancer <- as.numeric(batch_info$Cancer == "CASE") 
model <- model.matrix(~hasCancer, data = batch_info)
combat_edata <- ComBat(dat = data_1_2, batch = batch_info$batch, mod = model)
dim(combat_edata)
combat_edata[1:5,1:5]


######批处理后PCA分析及可视化
combat_edata.pca<-prcomp(t(combat_edata[1:5000,]),
                         scale=T,
                         rank=10,###要使用的主成分的最大数目。
                         retx=T,##是否应返回已旋转的变量。
                        )

pdf("E:\\yan\\APL\\sva_batch\\picture\\after_batch.pdf",height = 8,width = 10)
fviz_pca_ind(combat_edata.pca, 
             label="none", #隐藏个人标签
             addEllipses = TRUE, ####浓度椭圆
             #ellipse.type="norm",
             ellipse.level=0.9,##椭圆大小
             habillage = gsm_type,
             palette = c("#00afbb","#e7b800","#fc4e07","#FF9966"),###填充颜色
             #mean.point=F
            )
dev.off()

write.table(combat_edata,"E:\\yan\\APL\\data\\0rawdata\\1.0-数据预处理\\combat_edata.txt",quote=F,sep="\t")

combat_edata<-read.delim("E:\\yan\\APL\\data\\0rawdata\\1.0-数据预处理\\combat_edata.txt",
                         check.names = F,
                         stringsAsFactors = F)

AML_Sample_info_combat<-AML_Sample_info[match(colnames(combat_edata),AML_Sample_info$GSM),]

##########################################################################################################
#########整合数据集，求差异基因
health_gsm<-AML_Sample_info_combat$GSM[which(AML_Sample_info_combat$Disease=="healthy")]
####842个整合数据集中健康样本的GSM

case_gsm<-AML_Sample_info_combat$GSM[which(AML_Sample_info_combat$FAB=="M3")]
###117个整合数据集的M3型的GSM

case_health_DEG<-combat_edata[,match(c(case_gsm,health_gsm),colnames(combat_edata))]
case_health_DEG<-as.data.frame(case_health_DEG)
#case_health_DEG[1:5,1:5]
######从整合数据集中提取相关数据并整理

# library(limma)
# library(dplyr)
#####使用limma计算DEG
group<-data.frame(GSM=c(case_gsm,health_gsm),
                  type=c(rep("M3",length(case_gsm)),rep("healthy",length(health_gsm))))
group$PML<-as.vector(t(case_health_DEG[which(rownames(case_health_DEG)=="PML"),
                                             match(group$GSM,colnames(case_health_DEG))]))
group$RARA<-as.vector(t(case_health_DEG[which(rownames(case_health_DEG)=="RARA"),
                                              match(group$GSM,colnames(case_health_DEG))]))
####生成分类注释文件

group_list<-group$type
design <- model.matrix(~0+factor(group_list))
colnames(design)=levels(factor(group_list))
###M型相对于healthy的差异基因
contrast.matrix<-makeContrasts("M3-healthy",
                               levels = design)

fit <- lmFit(case_health_DEG,design)
fit2 <- contrasts.fit(fit, contrast.matrix) 
fit2 <- eBayes(fit2) 
DEG<-topTable(fit2, coef=1, n=Inf) %>% na.omit()  ## coef比较分组 n基因数
######求差异基因
head(DEG)
threshold <- as.factor(ifelse(DEG$adj.P.Val < 0.05 &abs(DEG$logFC) >= 1 ,
                              ifelse(DEG$logFC >= 1 ,'Up','Down'),'Not'))
#######阈值为p.adj<0.05,logFC>1
table(threshold)
names(threshold)<-rownames(DEG)

DEG_PML_RARA<-DEG[match(c("PML","RARA"),rownames(DEG)),c(1,4,5)]
#########查看PML和RARA基因

p_value_pearson<-cor.test(as.numeric(case_health_DEG[match("PML",rownames(case_health_DEG)),1:length(case_gsm)]),
                          as.numeric(case_health_DEG[match("RARA",rownames(case_health_DEG)),1:length(case_gsm)]),
                          method = c("pearson"))
# p_value_pearson$p.value
# p_value_pearson$estimate

p_value_spearman<-cor.test(as.numeric(case_health_DEG[match("PML",rownames(case_health_DEG)),1:length(case_gsm)]),
                           as.numeric(case_health_DEG[match("RARA",rownames(case_health_DEG)),1:length(case_gsm)]),
                           method = c("spearman"))

mean_h_PML<-mean(as.numeric(case_health_DEG[match("PML",rownames(case_health_DEG)),
                                                  (length(case_gsm)+1):ncol(case_health_DEG)]))
mean_h_RARA<-mean(as.numeric(case_health_DEG[match("RARA",rownames(case_health_DEG)),
                                                   (length(case_gsm)+1):ncol(case_health_DEG)]))

mean_M_PML<-mean(as.numeric(case_health_DEG[match("PML",rownames(case_health_DEG)),
                                                  1:length(case_gsm)]))
mean_M_RARA<-mean(as.numeric(case_health_DEG[match("RARA",rownames(case_health_DEG)),
                                                   1:length(case_gsm)]))

PML_cha<-mean_M_PML-mean_h_PML
RARA_cha<-mean_M_RARA-mean_h_RARA


DEG_PML_RARA$M_mean<-c(mean_M_PML,mean_M_RARA)
DEG_PML_RARA$health_mean<-c(mean_h_PML,mean_h_RARA)
DEG_PML_RARA$mean_diff<-c(PML_cha,RARA_cha)
DEG_PML_RARA$FC<-c(2^PML_cha,2^RARA_cha)

DEG_PML_RARA$pearson_cor<-c(p_value_pearson$estimate,p_value_pearson$estimate)
DEG_PML_RARA$pearson_pvalue<-c(p_value_pearson$p.value,p_value_pearson$p.value)

DEG_PML_RARA$spearman_cor<-c(p_value_spearman$estimate,p_value_spearman$estimate)
DEG_PML_RARA$spearman_pvalue<-c(p_value_spearman$p.value,p_value_spearman$p.value)


#####绘制火山图
library(ggplot2)
pdf("E:\\yan\\APL\\sva_batch\\picture\\batch_M3_DEG_volcano.pdf",height = 8,width = 10)
ggplot(DEG,aes(x=logFC,y=-log10(adj.P.Val),colour=threshold)) +
  
  xlab("log2(Fold Change)")+ylab("-log10(pvalue)") +
  
  geom_point(size = 2,alpha=1) +
  
  ylim(0,400) + xlim(-8,8) +
  
  scale_color_manual(values=c("blue","grey", "red"))+
  geom_vline(xintercept = c(-1, 1), lty = 2,colour="#000000")+ #增加虚线
  geom_hline(yintercept = c(1), lty = 2,colour="#000000")+
  theme(
    axis.text=element_text(size=10),
    axis.title=element_text(size=10)
  )
dev.off()

######绘制热图
DEG_gene<-rownames(DEG)[which(threshold!="Not")]
#####提取差异基因数据
M3_DEG<-case_health_DEG[match(DEG_gene,rownames(case_health_DEG)),]

annotation_col <- data.frame(
  type = factor(group$type,levels=c("M3","healthy")),
  PML=group$PML,
  RARA=group$RARA
)
rownames(annotation_col) = group$GSM
head(annotation_col)

ann_colors = list(
  #FAB = c(AML="#e3e59a", AMKL="#559a50"),
  type = c(M3="#8c90c2",healthy="#ec8970")
)

range(scale(M3_DEG))
bk <- seq(-3,3,0.01)

library(pheatmap)

pdf("E:\\yan\\APL\\sva_batch\\picture\\batch_M3_DEG_heatmap.pdf",height = 8,width = 10)
pheatmap(M3_DEG,#数据 
         #cluster_col =F,
         scale = "row",#行归一化
         border=F, #去每个热图格子边框线
         clustering_distance_cols = "euclidean",#列进行欧式距离聚类
         clustering_method = "ward.D",
         #cluster_row = FALSE,#不对行进行聚类
         #color = colorRampPalette(c("#4e8db6", "white", "#c65955"))(100),#自定义颜色
         color = c(colorRampPalette(colors = c("#4e8db6","white"))(length(bk)/2),
                   colorRampPalette(colors = c("white","#c65955"))(length(bk)/2)),
         show_colnames=F,#不显示列名
         show_rownames = F,
         main="batch_M3_healthy_DEG",#主标题
         annotation_col = annotation_col,#构建列注释信息
         annotation_colors = ann_colors,
         breaks = bk
)
dev.off()

DEG_downgene<-rownames(DEG)[which(threshold=="Down")]
DEG_upgene<-rownames(DEG)[which(threshold=="Up")]

###############与单个数据集的DEG的交集情况
combat_data1_down<-intersect(DEG_downgene_data1,DEG_downgene)
combat_data1_up<-intersect(DEG_upgene_data1,DEG_upgene)

combat_data2_down<-intersect(DEG_downgene_data2,DEG_downgene)
combat_data2_up<-intersect(DEG_upgene_data2,DEG_upgene)

data1_2_downgene<-intersect(DEG_downgene_data1,DEG_downgene_data2)
data1_2_upgene<-intersect(DEG_upgene_data1,DEG_upgene_data2)

combat_12_down<-intersect(DEG_downgene,data1_2_downgene)
combat_12_up<-intersect(DEG_upgene,data1_2_upgene)

write.table()
#BiocManager::install("UpSetR")
library(UpSetR)


aaa<-upset(fromList(list(data1_downgene=DEG_downgene_data1,
                         data1_upgene=DEG_upgene_data1,
                         data2_downgene=DEG_downgene_data2,
                         data2_upgene=DEG_upgene_data2,
                         combat_downgene=DEG_downgene,
                         combat_upgene=DEG_upgene)),
           nsets = 6,
           order.by = "freq",
           keep.order = F,
           mb.ratio = c(0.55,0.45),
           text.scale = rep(1.5,6),
           number.angles = 0,
           point.size = 3,
           line.size = 1,
           mainbar.y.label = "Intersection Size",
           sets.x.label = "Set Size"
)
print(aaa)



#####提取差异基因数据
########
M0_gsm<-AML_Sample_info_combat$GSM[which(AML_Sample_info_combat$FAB=="M0")]
M1_gsm<-AML_Sample_info_combat$GSM[which(AML_Sample_info_combat$FAB=="M1")]
M2_gsm<-AML_Sample_info_combat$GSM[which(AML_Sample_info_combat$FAB=="M2")]
M3_gsm<-AML_Sample_info_combat$GSM[which(AML_Sample_info_combat$FAB=="M3")]
M4_gsm<-AML_Sample_info_combat$GSM[which(AML_Sample_info_combat$FAB=="M4")]
M5_gsm<-AML_Sample_info_combat$GSM[which(AML_Sample_info_combat$FAB=="M5")]
M6_gsm<-AML_Sample_info_combat$GSM[which(AML_Sample_info_combat$FAB=="M6")]
M7_gsm<-AML_Sample_info_combat$GSM[which(AML_Sample_info_combat$FAB=="M7")]
length(c(M0_gsm,M1_gsm,M2_gsm,M3_gsm,M4_gsm,M5_gsm,M6_gsm,M7_gsm))

################################################################################################
###########上下调基因
data2_interDEG<-GSE122511_Dataset_2[match(c(down_inter,up_inter),rownames(GSE122511_Dataset_2)),
                                    c(match(c(M0_gsm,M1_gsm,M2_gsm,M3_gsm,
                                              M4_gsm,M5_gsm,M6_gsm,M7_gsm),
                                            colnames(GSE122511_Dataset_2)))]
data2_interDEG[1:5,1:5]

library(NMF)
#####使用M3与健康的DEG及所有的M型做NMF
set.seed(1234)
estim_r<-nmf(data2_interDEG,2:6,nrun=10)
plot(estim_r)#####cophenetic值随K变化的最大变动的前一个点
plot(2:6,estim_r$measures$cophenetic,type="b")

set.seed(1234)
res <- nmf(data2_interDEG, 5, nrun = 10)
coefmap(res)#####展示从最佳拟合结果中获得的簇（聚类数）和一致性矩阵的层次聚类。
coefmap(minfit(res))
basismap(res)#########W矩阵分解图

consensusmap(res, annColors=list(c='blue')
             , labCol='sample ', main='Cluster stability'
             , sub='Consensus matrix and all covariates')

group <- predict(res)
group <- as.data.frame(group)
group$group <- paste0('Cluster',group$group)
group$sample <- rownames(group)

sum(table(group$group))

library(pheatmap)

group$gsm<-colnames(data2_interDEG)
group$FAB<-AML_Sample_info$FAB[match(group$gsm,AML_Sample_info$GSM)]
group$PML<-as.vector(t(GSE122511_Dataset_2[which(rownames(GSE122511_Dataset_2)=="PML"),
                                           match(group$gsm,colnames(GSE122511_Dataset_2))]))
group$RARA<-as.vector(t(GSE122511_Dataset_2[which(rownames(GSE122511_Dataset_2)=="RARA"),
                                            match(group$gsm,colnames(GSE122511_Dataset_2))]))
group <- group[order(group$group),]

table(group$FAB)
data2_interDEG<-data2_interDEG[,match(group$gsm,colnames(data2_interDEG))]


annotation_col <- data.frame(
  FAB = factor(group$FAB,levels=c("M0","M1","M2","M3","M4","M5","M6","M7")),
  NMF_cluster = factor(group$group,levels=c("Cluster1","Cluster2","Cluster3","Cluster4","Cluster5")),
  PML=group$PML,
  RARA=group$RARA
)

rownames(annotation_col) = rownames(group)
head(annotation_col)
table(annotation_col$FAB,annotation_col$NMF_cluster)
ann_colors = list(
  FAB = c(M0="#99CC33",M1="#FF6666",M2="#336699",M3="#FF9900",M4="#33CC99",M5="#CC6600",M6="#CCCC00",M7="#99CCCC"),
  NMF_cluster = c(Cluster1="#8c90c2",Cluster2="#ec8970",Cluster3="#CCCC00",Cluster4="#0099CC",Cluster5="#009966")
  # PML=colorRampPalette(c(min(group$PML),median(group$PML),max(group$PML)),c("blue","white","red"),space = "RGB")
)

range(data2_interDEG)
pdf("E:\\yan\\APL\\NMF图\\M3_health_DEG_NMF.pdf",width=15,height=10)
bk <- seq(3,14,0.01)
#pheatmap(data2_allM_M3_DEG)
pheatmap(data2_interDEG,#数据 
         cluster_col = F,
         #cluster_row = F,#不对行进行聚类
         #scale = "row",#行归一化
         border=F, #去每个热图格子边框线
         clustering_distance_cols = "euclidean",#列进行欧式距离聚类
         clustering_method = "ward.D",
         color = colorRampPalette(c("#4e8db6", "white", "#c65955"))(100),#自定义颜色
         # color = c(colorRampPalette(colors = c("#4e8db6","white"))(length(bk)/2),
         #           colorRampPalette(colors = c("white","#c65955"))(length(bk)/2)),
         show_colnames = F,#不显示列名
         show_rownames = F,
         main="data2_M3_health_DEG_NMF",#主标题
         annotation_col = annotation_col,#构建列注释信息
         annotation_colors = ann_colors,
         #annotation_legend = FALSE,####去除图例
         # breaks = bk
)
dev.off()













