#library
library(clusterProfiler)
library(org.Drenale.eg.db)
library(tidyr)
library(dplyr)
library(cowplot)
library(ggplot2)
library(forcats)
library(enrichplot)
library(stringr)
library(rtracklayer)
#GUTS
drenale_proteome <- read.delim("D.renale functional_annotation.tsv", header = TRUE, sep = '\t')
degs <- read.csv("expression_matrix")

deseqgutsvshead = read.delim("deseq_guts_head.txt", header = TRUE, sep = "\t")
deseqgutsvsphar = read.delim("deseq_guts_pharynx.txt", header = TRUE, sep = "\t")

genesME <- read.delim("merged_turquoise.txt", header = FALSE)
names(genesME)[1] <- "GID"
Ortho<- read.table("Drenale_orthology_cov50_Id50_e-05_sim0.80.proteinortho.2.csv", header = TRUE, sep = "\t")
tx2 <- read.delim("Tx2gene", header = TRUE, sep = ",")


#GENES UP or DOWN
apply_regulation <- function(df) {
  df$Regulation <- ifelse(df$padj <= 0.005 & abs(df$log2FoldChange) >= log2(4), 
                          ifelse(df$log2FoldChange > 0, "UP", "DOWN"), 
                          "No DE")
  df$GID <- row.names(df)
  return(df)
}
#Guts
deseqgutsvshead <- apply_regulation(deseqgutsvshead)
deseqgutsvsphar <- apply_regulation(deseqgutsvsphar)
deseqgutsvshead <- deseqgutsvshead[, c("GID", "Regulation")]
deseqgutsvsphar <- deseqgutsvsphar[, c("GID", "Regulation")]
degs <- merge(deseqgutsvsphar, deseqgutsvshead, by = "GID")
names(degs)[grep("^Regulation", names(degs))] <- c("Regulation_deseqgutsvshead", "Regulation_deseqgutsvspharynx")
ME <- degs
ME <- merge(degs, genesME, by = "GID")
ME <- ME[rowSums(ME == "No DE") < (ncol(ME) - 1), ]
ME<- ME[rowSums(is.na(ME)) < (ncol(ME) - 1), ]

#GO 
names(drenale_proteome)[1] <- "GID"
drenale_proteome$GID <- sub("\\.t[1-9]$", "", drenale_proteome$GID)
protME <- merge(genesME, drenale_proteome, by = "GID")

#Clusterprofiler
geneME <- ME$GID

#GO
go_mf           <-                                 enrichGO(gene = geneME,
                                                    OrgDb         = org.Drenale.eg.db,
                                                    ont           = "MF",
                                                    pAdjustMethod = "BH",
                                                    pvalueCutoff  = 0.05,
                                                    keyType       = "GID",
                                                    qvalueCutoff  = 0.2,
                                                    readable      = TRUE)
#Simplify
go_mf <- simplify(go_mf1, cutoff=0.7, by="p.adjust", select_fun=min)

#barplot 
bp <- barplot(go_bp, showCategory = 10, title = "Biological Process")
mf <- barplot(go_mf, showCategory = 10, title = "Molecular Function")
cc <- barplot(go_cc, showCategory = 10, title = "Cellular Component")
plot_grid(bp, mf, cc, ncol = 1)

#KEGG
DKEGG <- read.delim("DKEGG.tsv", header = TRUE, sep = " ")
DKEGG <- DKEGG[, c("KO", "GID")]
DPATHWAYS <- read.delim("DKPATHWAYS", header = TRUE, sep = " ")
DPATHWAYS <- DPATHWAYS[, c("PATHWAY", "GID")]
DPATHWAYS <- DPATHWAYS[grep("^map", DPATHWAYS$PATHWAY), ]
DPATHWAYS <- merge(DPATHWAYS, tx2, by = "GID", all.x = "TRUE")

x <- enricher(geneME, pvalueCutoff = 0.05, qvalueCutoff = 0.02, pAdjustMethod = "BH", universe = universe, TERM2GENE = DPATHWAYS)
filtered_x <- as.data.frame(x)

descriptions <- c(
  "Glycerolipid metabolism",
  "Fat digestion and absorption",
  "Sulfur metabolism",
  "Metabolism of xenobiotics by cytochrome P450",
  "Drug metabolism - cytochrome P450",
  "Arachidonic acid metabolism",
  "Chemical carcinogenesis - DNA adducts",
  "Biosynthesis of secondary metabolites",
  "Glycine, serine and threonine metabolism",
  "Sulfur relay system")

filtered_x <- mutate(filtered_x, descriptions = descriptions)
filtered_x <- mutate(filtered_x, richFactor = Count / as.numeric(sub("/\\d+", "", BgRatio)))

#ggplot
ggplot(filtered_x, showCategory = 3, 
       aes(richFactor, fct_reorder(descriptions, richFactor))) + 
  geom_segment(aes(xend=0, yend = descriptions)) +
  geom_point(aes(color=p.adjust, size = Count)) +
  scale_color_viridis_c(guide=guide_colorbar(reverse=TRUE)) +
  scale_size_continuous(range=c(2, 10)) +
  theme_minimal() + 
  xlab("Rich factor") +
  ylab(NULL) + 
  ggtitle("KEGG Enriched")

#Ort 
Ortho <- Ortho[, c("Drenale.fasta", "Rculicivorax.fasta", "Sbaturini.fasta", "Tspiralis.fasta", "Tmuris.fasta", "Celegans.fasta")]
Ortho <- Ortho %>%
  separate_rows(Drenale.fasta, sep = ",")
Ortho$Drenale.fasta <- sub("Drenale\\|", "", Ortho$Drenale.fasta)
names(Ortho)[1] <- "GID"
protME <- merge(tx2, genesME, by ="GID")
protME$GID <- NULL
names(protME)[1] <- "TXNAME"
Orthome <- merge(Ortho, protME, by = "TXNAME")
Orthome$Counts <- apply(Orthome, 1, function(row) sum(grepl("\\*", row)))
Orthome <- subset(Orthome, Counts != 5)
names(Ortho)[25] <- "TXNAME"
#Filtrar por Orth
#C.elegans
Ortho <- Ortho[, c("Drenale.fasta", "Celegans.fasta")]
Ortho <- Ortho %>%
  separate_rows(Celegans.fasta, sep = ",")
Ortho$Celegans.fasta <- sub("Celegans\\|", "", Ortho$Celegans.fasta)
Ortho$Drenale.fasta <- sub("Drenale\\|", "", Ortho$Drenale.fasta)
names(Ortho)[1] <- "TXNAME"

library(org.Ce.eg.db)
x <- org.Ce.egENSEMBLPROT2EG
x <- as.data.frame(x)
x$prot_id <- sub("\\.[1-9]$", "", x$prot_id)
names(x)[2] <- "Celegans.fasta"
Orthome <- merge(Orthome, x, by = "Celegans.fasta")
y <- org.Ce.egGENENAME
y <- as.data.frame(y)
Orthome <- merge(Orthome, y, by = "gene_id")
Orthome$GID <- sub("\\.t[1-9]$", "", Orthome$TXNAME)
Orthome <- merge(Orthome, ME, by = "GID", all.x = TRUE)
Orthome$GID <- NULL
Orthome <- unique(Orthome)
Orthome1 <- Orthome

