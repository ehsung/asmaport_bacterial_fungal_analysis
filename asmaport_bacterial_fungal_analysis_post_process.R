# Install Packages
# library(devtools)
# install_github("zdk123/SpiecEasi")
# install.packages("microeco")
#
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("phyloseq")
#
# install.packages("ggraph")

# Load Packages
library(SpiecEasi)
library(microeco)
library(phyloseq)
library(dplyr)
library(tibble)
library(igraph)
library(ggraph)


# Running this part again because I didn't save the phyloseq data ----
# Load in mainData
ASV16S_raw <- read.csv("ASV16S_mainData.csv")
ASVITS_raw <- read.csv("ASVITS_mainData.csv")
metadata_raw <- read.csv("metadata_cleaned.csv")


# Function to subset dataframe into diseases
subset_func <- function(rawASV, rawITS, rawMeta, disease_input) {

  # Subset samples to input disease
  disease_samplelist <- rawMeta %>%
    filter(disease == disease_input) %>%
    pull(Sample_Name)

  # Subset ASV down to input disease samplelist
  disease_ASV_df <- rawASV %>%
    select(all_of(c("X", disease_samplelist, "Taxonomy")))

  # Subset ITS down to input disease samplelist
  disease_ITS_df <- rawITS %>%
    select(all_of(c("X", disease_samplelist, "Taxonomy")))

  # Remove 0 and 1 read rows
  non_zero_singleton_disease_ASV_df <- disease_ASV_df %>%
    filter(rowSums(select(., where(is.numeric))) > 1)

  non_zero_singleton_disease_ITS_df <- disease_ITS_df %>%
    filter(rowSums(select(., where(is.numeric))) > 1)

  # Set Threshold to 2 -- for ASV needs to be in at least 2 samples
  threshold <- 2

  non_rare_ASV_df <- non_zero_singleton_disease_ASV_df %>%
    filter(rowSums(select(., all_of(disease_samplelist)) > 1) >= threshold)

  non_rare_ITS_df <- non_zero_singleton_disease_ITS_df %>%
    filter(rowSums(select(., all_of(disease_samplelist)) > 1) >= threshold)


  # Keep >100 read rows --- This can be used to filter ASVs below a certain read
  # Currently set to 0 and is not used. Is kept just in case.
  cleaned_disease_ASV_df <- non_rare_ASV_df %>%
    filter(rowSums(select(., where(is.numeric))) >= 0)

  cleaned_disease_ITS_df <- non_rare_ITS_df %>%
    filter(rowSums(select(., where(is.numeric))) >= 0)


  # Split into OTUs and TAXs dataframes for both AVS and ITS
  ASV_OTU <- cleaned_disease_ASV_df %>%
    select(-Taxonomy) %>%
    column_to_rownames(var = "X")

  ASV_TAX <- cleaned_disease_ASV_df %>%
    select(X, Taxonomy) %>%
    column_to_rownames(var = "X")

  ITS_OTU <- cleaned_disease_ITS_df %>%
    select(-Taxonomy) %>%
    column_to_rownames(var = "X")

  ITS_TAX <- cleaned_disease_ITS_df %>%
    select(X, Taxonomy) %>%
    column_to_rownames(var = "X")

  return(list(ASV_OTU, ASV_TAX, ITS_OTU, ITS_TAX))
}


# Generate separate disease dataframes
disease_HC_data <- subset_func(ASV16S_raw, ASVITS_raw, metadata_raw, "HC")
disease_AR_data <- subset_func(ASV16S_raw, ASVITS_raw, metadata_raw, "AR")
disease_ARAS_data <- subset_func(ASV16S_raw, ASVITS_raw, metadata_raw, "ARAS")
disease_AS_data <- subset_func(ASV16S_raw, ASVITS_raw, metadata_raw, "AS")



# Print dim of dataframes
dim_print_func <- function(df, disease_input) {
  print(disease_input)
  print(dim(df[[1]]))
  print(dim(df[[2]]))
  print(dim(df[[3]]))
  print(dim(df[[4]]))
  cat("\n")
}

print("Disease_Group")
print("ASV_OTU dimension")
print("ASV_TAX dimension")
print("ITS_OTU dimension")
print("ITS_TAX dimension")
cat("\n")

dim_print_func(disease_HC_data, "HC")
dim_print_func(disease_AR_data, "AR")
dim_print_func(disease_ARAS_data, "ARAS")
dim_print_func(disease_AS_data, "AS")

# Create phyloseq objects
create_phyloseq_object_function <- function(df) {
  ASV16S_otu_matrix <- otu_table(as.matrix(df[[1]]), taxa_are_rows = TRUE)
  ASV16S_tax_matrix <- tax_table(as.matrix(df[[2]]))
  ASVITS_otu_matrix <- otu_table(as.matrix(df[[3]]), taxa_are_rows = TRUE)
  ASVITS_tax_matrix <- tax_table(as.matrix(df[[4]]))

  ASV16S_phyloseq <- phyloseq(ASV16S_otu_matrix, ASV16S_tax_matrix)
  ASVITS_phyloseq <- phyloseq(ASVITS_otu_matrix, ASVITS_tax_matrix)

  return(list(ASV16S_phyloseq, ASVITS_phyloseq))
}

HC_ASV_ITS_phyloseq <- create_phyloseq_object_function(disease_HC_data)
AR_ASV_ITS_phyloseq <- create_phyloseq_object_function(disease_AR_data)
ARAS_ASV_ITS_phyloseq <- create_phyloseq_object_function(disease_ARAS_data)
AS_ASV_ITS_phyloseq <- create_phyloseq_object_function(disease_AS_data)






# Start of post processing script to generate gml files ----
# Load in mainData rds files
TwoSamReadThreshold_hc_spiec_easi <- readRDS("TwoSamReadThreshold_hc_spiec_easi.rds")
TwoSamReadThreshold_as_spiec_easi <- readRDS("TwoSamReadThreshold_as_spiec_easi.rds")
TwoSamReadThreshold_ar_spiec_easi <- readRDS("TwoSamReadThreshold_ar_spiec_easi.rds")
TwoSamReadThreshold_aras_spiec_easi <- readRDS("TwoSamReadThreshold_aras_spiec_easi.rds")


obtain_network_information <- function(input_name, speci_easi_obj, combined_phyloseq_obj) {
  # Step 1: Build igraph object from SpiecEasi result
  net <- adj2igraph(getRefit(speci_easi_obj))

  # Step 2: Assign seqType (ASV or ITS)
  seqType <- c(rep("ASV", ntaxa(combined_phyloseq_obj[[1]])), rep("ITS", ntaxa(combined_phyloseq_obj[[2]])))
  V(net)$seqType <- seqType

  # Step 3: Assign Node seqName
  seqName <- c(taxa_names(combined_phyloseq_obj[[1]]), taxa_names(combined_phyloseq_obj[[2]]))
  V(net)$seqName <- seqName

  # Step 4: Assign Node taxaName
  taxName <- c(as.character(as.data.frame(tax_table(combined_phyloseq_obj[[1]]))$Taxonomy), as.character(as.data.frame(tax_table(combined_phyloseq_obj[[2]]))$Taxonomy))
  V(net)$taxName <- taxName

  # Step 5: Get edges and weights
  beta_mat <- getOptBeta(speci_easi_obj)
  edge_vals <- beta_mat[lower.tri(beta_mat)][which(getRefit(speci_easi_obj)[lower.tri(getRefit(speci_easi_obj))] == 1)]

  # Step 6: Assign edge attributes (weights and signs)
  E(net)$weight <- abs(edge_vals)
  E(net)$sign <- ifelse(edge_vals > 0, "positive", "negative")

  # Step 7: Export to GML (Gephi input)
  write_graph(net, file = paste0(input_name, "_gephi_network.gml"), format = "gml")
}



# Obtain gml for gephi
obtain_network_information("TwoSamReadThreshold_hc_spiec_easi", TwoSamReadThreshold_hc_spiec_easi, HC_ASV_ITS_phyloseq)
obtain_network_information("TwoSamReadThreshold_as_spiec_easi", TwoSamReadThreshold_as_spiec_easi, AS_ASV_ITS_phyloseq)
obtain_network_information("TwoSamReadThreshold_ar_spiec_easi", TwoSamReadThreshold_ar_spiec_easi, AR_ASV_ITS_phyloseq)
obtain_network_information("TwoSamReadThreshold_aras_spiec_easi", TwoSamReadThreshold_aras_spiec_easi, ARAS_ASV_ITS_phyloseq)
