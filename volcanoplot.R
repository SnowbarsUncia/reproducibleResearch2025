setwd ("C:\\Users\\User\\Desktop\\VosprIssled") # Установка рабочей директории

# Загрузка необходимых пакетов:
if (!("openxlsx" %in% installed.packages())) install.packages("openxlsx")
library (openxlsx)
if (!("ggplot2" %in% installed.packages())) install.packages("ggplot2")
library (ggplot2)
if (!("ggpubr" %in% installed.packages())) install.packages("ggpubr")
library(ggpubr)
if (!("scales" %in% installed.packages())) install.packages("scales")
library(scales)
if (!("BiocManager" %in% installed.packages())) install.packages("BiocManager")
library(BiocManager)
BiocManager::install("EnhancedVolcano")
BiocManager::install("DESeq2")
library(EnhancedVolcano)
library(DESeq2)

#Анализ дифференциальной экспрессии
count_table <- read.delim("allSamples.featureCounts.txt", skip=1, row.names="Geneid")
sample_table <- data.frame(condition=c("DL", "DL", "DL", "control", "control", "control"))

ddsFullCountTable <- DESeqDataSetFromMatrix(countData = count_table[,6:11], colData = sample_table, design = ~ condition)

dds <- DESeq(ddsFullCountTable)
res <- results(dds)

#Визуализация данных
EnhancedVolcano(res, lab = rownames(res),
                x = 'log2FoldChange', y = 'pvalue',
                pCutoff=0.05, pCutoffCol = 'padj', FCcutoff = 1,
                title="S. cerevisiae", subtitle="5 mM D-lactate vs control",
                col = c("grey30", "grey30", "grey30", "red2"),
                xlab="", ylab = bquote(~-Log[10] ~ italic(p)),
                caption="", selectLab = "", legendPosition = 'none')

DEGs <- res[abs(res$log2FoldChange) > 1 & res$padj < 0.05 & complete.cases(res$padj), ]
DEGs <- DEGs[order(DEGs$log2FoldChange), ]
DEGs$Transcript <- row.names(DEGs)
library(openxlsx)
write.xlsx(x = DEGs, file = "DEGs_yeast.xlsx", rowNames = TRUE)

EnhancedVolcano(res, lab = rownames(res),
                x = 'log2FoldChange', y = 'pvalue',
                pCutoff=0.05, pCutoffCol = 'padj', FCcutoff = 1,
                title="S. cerevisiae", subtitle="5 mM D-lactate vs control",
                col = c("grey30", "grey30", "grey30", "red2"),
                xlab="", ylab = bquote(~-Log[10] ~ italic(p)),
                caption="", selectLab = rownames(DEGs), legendPosition = 'none')
