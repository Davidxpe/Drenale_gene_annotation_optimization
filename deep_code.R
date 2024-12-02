#Scrip para dividir en tres columnas BP, MF, CC
library(tidyr)
library(dplyr)
drenale_original <- read.delim("D.renale functional_annotation.tsv")
#go.obo en df
GOterms <- read.delim(file = "GOterms.tsv", header = TRUE, sep = "\t")

#Quedarse solo con el ID, el término GO para dividir por categoría 
drenale_GO <-drenale_original[,c("Sequence.Name", "GOs")]
drenale_GO <- drenale_GO %>%
  separate_rows(GOs, sep = ",")
drenale_GO <- drenale_GO %>%
  left_join(GOterms[, c("ID", "Namespace")], by = c("GOs" = "ID"))
concatenate_values <- function(x) {
  if (length(x) > 1) {
    return(paste(x, collapse = ", "))
  } else {
    return(x)
  }
}
drenale_GO <- drenale_GO %>%
  pivot_wider(names_from = Namespace, values_from = GOs,
              values_fn = list(GOs = concatenate_values),
              names_sep = "_", names_repair = "minimal")


#Scrips para quedarse con los más "profundos", ejemplo con "cellular_component

drenale_only_GO$cellular_component <- ifelse(is.na(drenale_only_GO$cellular_component), NA, strsplit(drenale_only_GO$cellular_component, ", "))
library(ontologyIndex)
#deep función
get_filtered_descendants <- function(terms, cell_contents) {
  terms_without_filtered_descendants <- character(0)
  for (term in unlist(terms)) {
    if (!is.na(term)) {
      descendants <- get_descendants(term, ontology = go_ontology)
      filtered_descendants <- setdiff(intersect(descendants, unlist(cell_contents)), term)
      if (length(filtered_descendants) == 0) {
        terms_without_filtered_descendants <- c(terms_without_filtered_descendants, term)
      }
    }
  }
  
  return(terms_without_filtered_descendants)
}

#columna donde se almacena deep
drenale_only_GO$deep_cc <- NA
for (i in 1:nrow(drenale_only_GO)) {
  terms <- unlist(drenale_only_GO$cellular_component[i])
  terms_filtered <- get_filtered_descendants(terms, drenale_only_GO$cellular_component[i])
  terms_combined <- paste(terms_filtered, collapse = ", ")
  drenale_only_GO$deep_cc[i] <- terms_combined
}
