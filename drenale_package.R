#drenale_db
DSym <- read.csv("DSym.csv", header = TRUE, sep = ",")
DGO <- read.csv("DGO.csv", header = TRUE, sep = ",")
DChr <- read.csv("DChr.csv", header = TRUE, sep = ",")
DKegg <- read.table("DKEGG.tsv", header = TRUE, sep = " ")
DKPathway <- read.table("DKPATHWAYS", header = TRUE, sep = " ")
tx2gene <- read.csv("Tx2gene", header = TRUE, sep = ",")

#library
library(dplyr)
library(tidyr)
library(AnnotationForge)
library(rtracklayer)
library(stringr)

#DataGO
drenale_original <- read.delim("D.renale functional_annotation.tsv", header = TRUE, sep = "\t")
deep <- read.delim("Drenale Deepest GO CLEAN", header = TRUE, sep = "\t")
goterms <- read.delim("GOterms.tsv")

#DataGFF3
gff3_data <- import("Drenale_annotation.gff3", format = "GFF")
gff3_data = as.data.frame(gff3_data)

#create dsym, dgo, dchr
#dsym
DSym <- drenale_original[,c("Sequence.Name")]
DSym <- as.data.frame(DSym)
DSym$GID <- sub("\\.t[1-9]$", "", DSym$DSym)
DSym = unique(DSym)
DSym$SYMBOL <- DSym$DSym
DSym <- DSym[!is.na(DSym$ANNOTA), ]
colnames(DSym) <- c("GID", "SYMBOL")

#dchr
DChr <- gff3_data[,c("ID", "seqnames", "type")]
DChr <- DChr[DChr[, "type"] == "gene",]
DChr$type <- NULL
colnames(DChr) <- c("GID", "CHROMOSOME")
DChr <- DChr[!is.na(DChr$GID), ]
DChr <- DChr[!is.na(DChr$CHROMOSOME), ]
DChr <- DSym

#dgo ALLGOS
drenale_only_GO <- drenale_original[, c("Sequence.Name", "Mapping.GO.Term", "GOs")]
caracteres_deseados <- c("IEA", "IDA", "IBA", "TAS", "ISS", "IMP", "IPI", "HDA", "NAS", "IGI", "ND", "ISA", "IC", "IEP", "EXP", "HMP", "RCA", "HEP", "ISO", "ISM", "IGC", "ADA")
drenale_only_GO$Mapping.GO.Term <- gsub(paste0(".*\\b(", paste(caracteres_deseados, collapse="|"), ")\\b.*"), "\\1", drenale_only_GO$Mapping.GO.Term)
drenale_only_GO$Sequence.Name <- sub("\\.t[1-9]$", "", drenale_only_GO$Sequence.Name)
drenale_only_GO = unique(drenale_only_GO)
drenale_only_GO <- drenale_only_GO[!is.na(drenale_only_GO$GOs), ]
drenale_only_GO <- drenale_only_GO[drenale_only_GO[,3]!="-",]
DGO <- drenale_only_GO
colnames(DGO) <- c("GID", "GO", "EVIDENCE")
DGO <- DGO[!is.na(DGO$GID), ]
DGO <- DGO[!is.na(DGO$GO), ]
DGO <- DGO[!is.na(DGO$EVIDENCE), ]
colnames(DGO) <- c("GID", "EVIDENCE", "GO")
DGO <- DGO[,c("GID", "GO", "EVIDENCE")]
DGO$EVIDENCE <- ifelse(is.na(DGO$EVIDENCE) | DGO$EVIDENCE == " ", "ND", DGO$EVIDENCE)
DGO <- DGO %>%
  separate_rows(GO, sep = ",")

filas_obsoletas <- grepl("obsolete", goterms$Name, ignore.case = TRUE)
DGO_filtrado <- DGO[!(DGO$GO %in% goterms$ID[filas_obsoletas]), , drop = FALSE]
DGO <- DGO_filtrado 

library(tidyr)

#dgo ONLY DEEPS
deep$deep = paste(deep$deep_bp, deep$deep_mf, deep$deep_cc)
deep$Sequence.Name <- sub("\\.t[1-9]$", "", deep$Sequence.Name)
deep <- deep[, c("Sequence.Name","deep")]
deep <- as.data.frame(deep)
deep$deep <- sapply(strsplit(deep$deep, " "), function(x) paste(x, collapse = ", "))
deep$deep <- sapply(strsplit(deep$deep, ",\\s*(,\\s*)*"), function(x) paste(x, collapse = ", "))
DGO <- deep %>%
  separate_rows(deep, sep = ",")
DGO$deep <- str_trim(DGO$deep)
DGO$deep[DGO$deep == "NA"] <- NA
DGO <- DGO[!is.na(DGO$deep), ]
colnames(DGO) <- c("GID", "GO")
merged_df <- merge(DGO, drenale_only_GO, by = "GID", all.x = TRUE)
merged_df$GO.y <- NULL
DGO = merged_df
colnames(DGO) <- c("GID", "GO", "EVIDENCE")
DGO$EVIDENCE <- ifelse(is.na(DGO$EVIDENCE) | DGO$EVIDENCE == "", "ND", DGO$EVIDENCE)
DGO <- unique(DGO)

#dkegg
drenale_original$Sequence.Name <- sub("\\.t[1-9]$", "", drenale_original$Sequence.Name)
dkegg <- drenale_original[, c("Sequence.Name", "KEGG_ko")]
dkegg <- dkegg[dkegg[,2]!="-",]
dkegg <- dkegg[!is.na(dkegg$KEGG_ko), ]
dkegg <- dkegg %>%
  separate_rows(KEGG_ko, sep = ",")
dkegg <- unique(dkegg)
colnames(dkegg) <- c("GID", "KO")
#--Borrar KO:--#
dkegg <- dkegg %>%
  mutate(KO = str_remove(KO, "ko:"))


#dkpathways
dkpathways <- drenale_original[, c("Sequence.Name", "KEGG_Pathway")]
dkpathways <- dkpathways[dkpathways[,2]!="-",]
dkpathways <- dkpathways[!is.na(dkpathways$KEGG_Pathway), ]
dkpathways <- dkpathways %>%
  separate_rows(KEGG_Pathway, sep = ",")
dkpathways <- unique(dkpathways)
colnames(dkpathways) <- c("GID", "PATHWAY")

#proteins

Ortho<- read.table("Drenale_orthology_cov50_Id50_e-05_sim0.80.proteinortho.2.csv", header = TRUE, sep = "\t")
Ortho <- Ortho[, c("Drenale.fasta")]
Ortho <- as.data.frame(Ortho)
Ortho <- Ortho %>%
  separate_rows(Ortho, sep = ",")
Ortho$Ortho <- sub("Drenale\\|", "", Ortho$Ortho)
Ortho <- subset(Ortho, Ortho != "*")
DProt <- DSym
DProt <- DProt[, c("GID", "DSym")]
names(DProt)[2] <- "PROTEINS"

#Control

comunDSymvsDGO <- merge(drenale_only_GO, DSym, by = "GID")
only_DGO <- anti_join(drenale_only_GO, DChr, by = "GID")
only_DGOdepvsDGO <- anti_join(DChr, drenale_only_GO, by = "GID")


DChr <- bind_rows(DChr, anti_join(DSym, DChr, by = "GID"))
DChr$SYMBOL <- NULL
DChr$CHROMOSOME <- as.character(DChr$CHROMOSOME)
DChr$CHROMOSOME <- ifelse(is.na(DChr$CHROMOSOME), "NOCHR1", DChr$CHROMOSOME)

#"Unique"
DSym <- unique(DSym)
DChr <- unique(DChr)
DGO <- unique(DGO)
dkegg <- unique(dkegg)
dkpathways <- unique(dkpathways)

write.csv(DSym, "DSym.csv", row.names = FALSE)
write.csv(DChr, "DChr.csv", row.names = FALSE)
write.csv(DGO, "DGO.csv", row.names = FALSE)
write.table(dkegg, "DKEGG.tsv", row.names = FALSE)
write.table(dkpathways, "DKPATHWAYS", row.names = FALSE)
write.table(DProt, "DProt", row.names =  FALSE)

#ensamblar .org.db

makeOrgPackage(gene_info=DSym, chromosome=DChr, go=DGO, ko=DKegg, pathway=DKPathway, transcripID=tx2gene,
               version="0.2",
               maintainer="davidcalupina08@gmail.com",
               author="David Calupiña",
               outputDir = "",
               tax_id="513045",
               genus="Dioctophyme",
               species="renale",
               goTable="go")

install.packages("org.Drenale.eg.db", repos = NULL)


                                                                                                                                                                                                    