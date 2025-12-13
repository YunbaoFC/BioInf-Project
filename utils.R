################################################################################
# UTILS
################################################################################

getwd()
# setwd("Bureau/Sapienza/BioInf/BioInf-Project") # <- Modify this
getwd() # <- To check you're in the right directory

###############################################################################
# Function to convert from Genes Ensembl -> Genes Symbol
###############################################################################
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl") # Connect to Ensembl
translate_genes <- function(ids, map) {
  sapply(ids, function(x) {
    sym <- map$hgnc_symbol[map$ensembl_gene_id == x]
    if(length(sym) == 0 || sym == "") x else sym
  })
}
translate_peptide <- function(ids, map) {
  sapply(ids, function(x) {
    sym <- map$hgnc_symbol[map$ensembl_peptide_id == x]
    if(length(sym) == 0 || sym == "") x else sym
  })
}

################################################################################
# 1. Filter the BIOGRID-ORGANISM-Homo_sapiens file
################################################################################

# Read the original file
biogrid <- read.delim("Files/Original_files/BIOGRID-ORGANISM-Homo_sapiens-5.0.252.tab3.txt")
# Filter what we want to get
human_interactions <- subset(
  biogrid,
  Organism.ID.Interactor.A == 9606 &
  Organism.ID.Interactor.B == 9606 &
  Experimental.System.Type == "physical"  
)
# Only keep the 2 columns Official.Symbol.Interactor A and B
ppi_edges <- human_interactions[, c("Official.Symbol.Interactor.A", "Official.Symbol.Interactor.B")]
head(ppi_edges)

################################################################################
# 2. Translate Huri Genes Ensemble ID into Genes Symbols
################################################################################
#BiocManager::install("biomaRt") # If package not installed yet
library(biomaRt)
huri <- read.delim("Files/Original_files/HI-union.tsv", header = FALSE, stringsAsFactors = FALSE) # Load original file
colnames(huri) <- c("GeneA", "GeneB")
all_genes <- unique(c(huri$GeneA, huri$GeneB)) # Get IDs
# biomaRt request
gene_map <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters = "ensembl_gene_id",
  values = all_genes,
  mart = ensembl
)
# Replace in table
huri$GeneA <- translate_genes(huri$GeneA, gene_map)
huri$GeneB <- translate_genes(huri$GeneB, gene_map)
head(huri)

################################################################################
# 3. Filter String + Translate
################################################################################
library(tidyr)
String <- read.delim("Files/Original_files/9606.protein.links.v12.0.onlyAB.txt")
# Split in 3 distinct columns
String_sep <- separate(String, col = "protein1.protein2.combined_score", into = c("ProteinA", "ProteinB", "Score"), sep = " ")
# Clean IDs
String_sep$ProteinA <- sub("^9606\\.", "", String_sep$ProteinA)
String_sep$ProteinB <- sub("^9606\\.", "", String_sep$ProteinB)
# Unique vectors
all_genes_string <- unique(c(String_sep$ProteinA, String_sep$ProteinB))
# biomaRt request
gene_map <- getBM(
  attributes = c("ensembl_peptide_id", "hgnc_symbol"),
  filters = "ensembl_peptide_id",
  values = all_genes_string,
  mart = ensembl
)
# Translate in table
String_sep$ProteinA <- translate_peptide(String_sep$ProteinA, gene_map)
String_sep$ProteinB <- translate_peptide(String_sep$ProteinB, gene_map)
head(String_sep)

################################################################################
# 3. Filter Reactome
################################################################################
Reactome <- read.delim("Files/Original_files/FIsInGene_04142025_with_annotations.txt")
reactome_genes <- subset(Reactome, select = c(Gene1, Gene2))

################################################################################
# 4. Clean redundant interactions + self loops + isolate LCC
################################################################################
library(igraph)
clean_interactome <- function(df,
                              colA = 1,
                              colB = 2,
                              score_col = NULL) {
  
  # Rename columns
  df <- df[, c(colA, colB, score_col), drop = FALSE]
  colnames(df)[1:2] <- c("GeneA", "GeneB")
  #if (!is.null(score_col)) colnames(df)[3] <- "Score"
  
  # Surp les self-loops
  df <- df[df$GeneA != df$GeneB, ]
  # Supr doubles
  df$pair_id <- apply(df[, c("GeneA", "GeneB")], 1, function(x)
    paste(sort(x), collapse = "_"))
  if (!is.null(score_col)) {
    df <- aggregate(Score ~ pair_id + GeneA + GeneB,
                    data = df,
                    FUN = max)
  } else {
    df <- df[!duplicated(df$pair_id), ]
  }
  df$pair_id <- NULL
  # Get graph
  g <- graph_from_data_frame(df, directed = FALSE)
  # Get only LCC
  comp <- components(g)
  lcc_id <- which.max(comp$csize)
  g_lcc <- induced_subgraph(
    g,
    vids = V(g)[comp$membership == lcc_id]
  )
  # Filter Table for LCC
  lcc_nodes <- V(g_lcc)$name
  df_lcc <- df[
    df$GeneA %in% lcc_nodes &
      df$GeneB %in% lcc_nodes, ]
  
    return(list(
    table_full_clean = df,
    table_LCC = df_lcc,
    graph_LCC = g_lcc
  ))
}

biogrid <- clean_interactome(ppi_edges[, c("Official.Symbol.Interactor.A", "Official.Symbol.Interactor.B")])
huri <- clean_interactome(huri[, c("V1", "V2")])
string <- clean_interactome(String_sep[, c("ProteinA", "ProteinB", "Score")])
reactome <- clean_interactome(reactome_genes[, c("Gene1", "Gene2")])

# Write new files
write.table(biogrid$table_LCC,file = "Files/Biogrid.txt",sep = "\t",row.names = FALSE,quote = FALSE)
write.table(huri$table_LCC, "Files/Huri.txt", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(string$table_LCC, "Files/String.txt", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(reactome$table_LCC, "Files/Reactome.txt", sep = "\t", row.names = FALSE, quote = FALSE)

################################################################################
# 5. Filter Cardiomyopathy file
################################################################################
Cardiomyopathy <- read.delim("Files/Original_files/diseases_.tsv")
Cardiomyopathy <- subset(Cardiomyopathy, select = c(Associated.genes))
colnames(Cardiomyopathy) <- c("Genes")
head(Cardiomyopathy)
write.table(Cardiomyopathy,file = "Files/Cardiomyopathy.txt",sep = "\t",row.names = FALSE,quote = FALSE)
