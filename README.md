# TCGA-LUAD-DESeq2-Analysis
Differential gene expression analysis of TCGA-LUAD RNA-Seq data using DESeq2 and Bioconductor packages.
DESeq2 ve Bioconductor paketleri kullanılarak TCGA-LUAD RNA-Seq verilerinin diferansiyel gen ekspresyon analizi.
 #1. Paket Kurulumu
 # ========================================
BiocManager::install(c(
  "TCGAbiolinks",
  "DESeq2",
  "SummarizedExperiment",
  "EnhancedVolcano",
  "pheatmap",
  "org.Hs.eg.db"))
# 2. Paketleri Yükle 
# ==========================================
library(TCGAbiolinks)
library(DESeq2)
library(SummarizedExperiment)
library(EnhancedVolcano)
library(pheatmap)
library(org.Hs.eg.db)
# 3. TCGA-LUAD RNA-seq VERİSİNİ SORGULA
# =========================================
query <- GDCquery(
  project       = "TCGA-LUAD",
  data.category = "Transcriptome Profiling",
  data.type     = "Gene Expression Quantification",
  workflow.type = "STAR - Counts")
  # 4. VERD0YD0 D0NDD0R VE HAZIRLA
# NOT: TCGA-LUAD ~600 C6rnek iC'erir, indirme uzun sC<rebilir.
# ==========================================
GDCdownload(query, files.per.chunk = 10)
data <- GDCprepare(query)
# 5. METAVERİLERİNİ HAZIRLA
# ==========================================
metadata <- colData(data)

# Yalnızca Tümör ve Normal doku örneklerini seç
metadata <- metadata[metadata$sample_type %in% c("Primary Tumor", "Solid Tissue Normal"), ]

# Count matrisini ve koşul vektörünü oluştur
counts    <- assay(data)[, rownames(metadata)]
condition <- factor(metadata$sample_type)
# 6. DESeq2 NESNESİNİ OLUŞTUR VE ANALİZİ ÇALIŞTIR
# ============================================
dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData   = data.frame(condition),
  design    = ~ condition)
#  DESeq2 sayımı genleri filtrele (toplam okuma sayısı > 10)
dds <- dds[rowSums(counts(dds)) > 10, ]
# DESeq2 analizini çalıştır
dds <- DESeq(dds)
# 7. DEG SONUÇLARINI ÇIKAR
# ============================================
# AyarlanmD1E p-deDerine gC6re sD1rala
res <- res[order(res$padj), ]

# Anlamlı genleri filtrele (padj < 0.05 ve |log2FC| > 1)
sig_genes <- subset(res, padj < 0.05 & abs(log2FoldChange) > 1)

# SonuC'ları CSV olarak kaydet
write.csv(
  as.data.frame(sig_genes),
  file         = "TCGA_LUAD_DEG_results.csv",
  fileEncoding = "UTF-8")
# 8. VOLCANO PLOT
# ==========================================
EnhancedVolcano(
  res,
  lab      = rownames(res),
  x        = "log2FoldChange",
  y        = "padj",
  pCutoff  = 0.05,
  FCcutoff = 1,
  title    = "TCGA-LUAD: Tumor vs Normal",
  subtitle = "DESeq2 Differential Expression")
  # 9. EN ÖNEMLİ 30 GENİN ISI HARİTASI
# ========================================
top30 <- head(rownames(sig_genes), 30)

# Varyans stabilizasyon (blind=FALSE: DEG bilgisi kullanılır)
vsd <- vst(dds, blind = FALSE)

# Z-skoru normalizasyonu için satır ortalaması çıkar
mat <- assay(vsd)[top30, ]
mat <- mat - rowMeans(mat)

# DÜZELTME: annotation_col iC'in row.names süttun isimleriyle eşleştirilmeli
pheatmap(
  mat,
  annotation_col = data.frame(
    condition = condition,
    row.names = colnames(mat)),
  show_rownames = TRUE,
  cluster_cols  = TRUE,
  main          = "Top 30 Differentially Expressed Genes")

