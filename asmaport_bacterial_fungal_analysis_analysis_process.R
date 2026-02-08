# Install Packages
install.packages("openxlsx")

# Load Packages
library(readxl)
library(dplyr)
library(tidyr)
library(openxlsx)

# Load in data
ar_data_df <- read_excel("TwoSamReadThreshold_ar_spiec_easi_gephi_network_graphdata.xlsx", sheet = "Analysis")		# ar_spiec-easi_network.xlsx
aras_data_df <- read_excel("TwoSamReadThreshold_aras_spiec_easi_gephi_network_graphdata.xlsx", sheet = "Analysis")	# aras_spiec-easi_network.xlsx
as_data_df <- read_excel("TwoSamReadThreshold_as_spiec_easi_gephi_network_graphdata.xlsx", sheet = "Analysis")		# as_spiec-easi_network.xlsx
hc_data_df <- read_excel("TwoSamReadThreshold_hc_spiec_easi_gephi_network_graphdata.xlsx", sheet = "Analysis")		# hc_spiec-easi_network.xlsx

# Select only the meaningful columns
ar_data_df_filtered <- ar_data_df %>%
  select(Source_ID, Target_ID, Source_Taxa, Target_Taxa, Weight, sign, ConnectionType)

aras_data_df_filtered <- aras_data_df %>%
  select(Source_ID, Target_ID, Source_Taxa, Target_Taxa, Weight, sign, ConnectionType)

as_data_df_filtered <- as_data_df %>%
  select(Source_ID, Target_ID, Source_Taxa, Target_Taxa, Weight, sign, ConnectionType)

hc_data_df_filtered <- hc_data_df %>%
  select(Source_ID, Target_ID, Source_Taxa, Target_Taxa, Weight, sign, ConnectionType)




"""# ConnectionType Counter"""

# Functions
Connection_Type_Counter <- function(df, weight_threshold) {
  "
  Filters down by a weight_threshold input and then counts the remaining connection types.
  "
  df_threshold <- df %>%
    filter(Weight > weight_threshold) %>% # Cut number of connections at weight_threshold
    count(ConnectionType) %>% # Count number of connection types
    mutate(Total_Connections = sum(n),
           ConnectionType_Percent = 100 * n / sum(n),
           Weight_Threshold = weight_threshold)


  return(df_threshold)
}

Weight_Threshold_Looper <- function(df) {
  "
  Apply the Connection_Type_Counter function and loop through thresholds to create a dataframe of connection type counts
  "
  thresholds <- c(0, 0.000001, 0.000005, 0.00001, 0.00005, 0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)

  # Apply function over thresholds and combine
  result_list <- lapply(thresholds, function(th) Connection_Type_Counter(df, th))
  result_df <- bind_rows(result_list)

  # Reshape dataframe
  final_table <- result_df %>%
    select(Weight_Threshold, Total_Connections, ConnectionType, ConnectionType_Percent) %>%
    pivot_wider(names_from = ConnectionType, values_from = ConnectionType_Percent)

  return(final_table)
}

# Obtain weight_thresholds dataframes
ar_weighted_cuts <- Weight_Threshold_Looper(ar_data_df_filtered)
aras_weighted_cuts <- Weight_Threshold_Looper(aras_data_df_filtered)
as_weighted_cuts <- Weight_Threshold_Looper(as_data_df_filtered)
hc_weighted_cuts <- Weight_Threshold_Looper(hc_data_df_filtered)

# Export
wb <- createWorkbook()

# Add sheets
addWorksheet(wb, "ar_weighted_cuts")
writeData(wb, "ar_weighted_cuts", ar_weighted_cuts)

addWorksheet(wb, "aras_weighted_cuts")
writeData(wb, "aras_weighted_cuts", aras_weighted_cuts)

addWorksheet(wb, "as_weighted_cuts")
writeData(wb, "as_weighted_cuts", as_weighted_cuts)

addWorksheet(wb, "hc_weighted_cuts")
writeData(wb, "hc_weighted_cuts", hc_weighted_cuts)

# Save the workbook
saveWorkbook(wb, "weighted_thresholds_ConnectionType_summary.xlsx", overwrite = TRUE)




"""# Taxonomy Compare"""

# Functions
Taxonomy_Compare_Counter <- function(df, weight_threshold) {
  "
  Filters down by a weight_threshold input and then counts the remaining connection types.
  "
  df_threshold <- df %>%
    filter(Weight > weight_threshold) # Cut number of connections at weight_threshold

  ASV_df_threshold <- df_threshold %>%
    filter(ConnectionType == "ASV")

  ITS_df_threshold <- df_threshold %>%
    filter(ConnectionType == "ITS")

  Mix_df_threshold <- df_threshold %>%
    filter(ConnectionType == "Mix")


  # Create lists of unique ASV/ITS
  ASV_list <- c(unique(c(ASV_df_threshold$Source_ID, ASV_df_threshold$Target_ID)))
  ITS_list <- c(unique(c(ITS_df_threshold$Source_ID, ITS_df_threshold$Target_ID)))
  Mix_list <- c(unique(c(Mix_df_threshold$Source_ID, Mix_df_threshold$Target_ID)))

  # Count the subsets
  ASV_list_count <- length(ASV_list)
  ITS_list_count <- length(ITS_list)
  Mix_list_count <- length(Mix_list)
  Mix_list_ITS_count <- sum(startsWith(Mix_list, "ASVITS"))
  Mix_list_ASV_count <- Mix_list_count - Mix_list_ITS_count

  # Count the unique ASV/ITS not found in the other connections
  ASV_unique_count <- length(setdiff(ASV_list, union(ITS_list, Mix_list)))
  ITS_unique_count <- length(setdiff(ITS_list, union(ASV_list, Mix_list)))
  Mix_unique_count <- length(setdiff(Mix_list, union(ASV_list, ITS_list)))

  Mix_unique <- setdiff(Mix_list, union(ASV_list, ITS_list))
  Mix_unique_ITS_count <- sum(startsWith(Mix_unique, "ASVITS"))
  Mix_unique_ASV_count <- Mix_unique_count - Mix_unique_ITS_count

  return(c(ASV_list_count, ITS_list_count, Mix_list_count, Mix_list_ASV_count, Mix_list_ITS_count, ASV_unique_count, ITS_unique_count, Mix_unique_count, Mix_unique_ASV_count, Mix_unique_ITS_count))
}


Taxonomy_Compare_Weight_Threshold_Looper <- function(df) {
  "
  Apply the Taxonomy_Compare_Counter function and loop through thresholds to create a dataframe of ASV, ITS, Mix counts
  "
  thresholds <- c(0, 0.000001, 0.000005, 0.00001, 0.00005, 0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)

  # Initialize an empty list to store results
  results_list <- list()

  # Loop through each threshold and apply your function
  for (th in thresholds) {
    counts <- Taxonomy_Compare_Counter(df, weight_threshold = th)

    # Store counts with the threshold labeled
    results_list[[as.character(th)]] <- c(Threshold = th,
                                          ASV_count = counts[1],
                                          ITS_count = counts[2],
                                          Mix_count = counts[3],
                                          Mix_ASV_count = counts[4],
                                          Mix_ITS_count = counts[5],
                                          ASV_unique_count = counts[6],
                                          ITS_unique_count = counts[7],
                                          Mix_unique_count = counts[8],
                                          Mix_unique_ASV_count = counts[9],
                                          Mix_unique_ITS_ITS = counts[10])
  }

  # Combine into a single data frame
  results_df <- do.call(rbind, results_list) %>% as.data.frame()

  # Clean rownames
  rownames(results_df) <- NULL

  # return result
  return(results_df)
}

# Obtain taxonomy list weight_thresholds dataframes
ar_taxonomy_weighted_cuts <- Taxonomy_Compare_Weight_Threshold_Looper(ar_data_df_filtered)
aras_taxonomy_weighted_cuts <- Taxonomy_Compare_Weight_Threshold_Looper(aras_data_df_filtered)
as_taxonomy_weighted_cuts <- Taxonomy_Compare_Weight_Threshold_Looper(as_data_df_filtered)
hc_taxonomy_weighted_cuts <- Taxonomy_Compare_Weight_Threshold_Looper(hc_data_df_filtered)

# Export
wb <- createWorkbook()

# Add sheets
addWorksheet(wb, "ar_taxonomy_weighted_cuts")
writeData(wb, "ar_taxonomy_weighted_cuts", ar_taxonomy_weighted_cuts)

addWorksheet(wb, "aras_taxonomy_weighted_cuts")
writeData(wb, "aras_taxonomy_weighted_cuts", aras_taxonomy_weighted_cuts)

addWorksheet(wb, "as_taxonomy_weighted_cuts")
writeData(wb, "as_taxonomy_weighted_cuts", as_taxonomy_weighted_cuts)

addWorksheet(wb, "hc_taxonomy_weighted_cuts")
writeData(wb, "hc_taxonomy_weighted_cuts", hc_taxonomy_weighted_cuts)

# Save the workbook
saveWorkbook(wb, "weighted_thresholds_taxonomy_sequence_summary.xlsx", overwrite = TRUE)





"""# Taxonomy Aggregate
Focus will be on the Mix ConnectionType, where we are examining between kingdoms. What we are seeing here are ASVs with the most connections/degrees and then top taxonomy with the most connections/degrees
"""

install.packages("openxlsx")
install.packages("UpSetR")

# Load Packages
library(readxl)
library(tidyverse)
library(dplyr)
library(tidyr)
library(openxlsx)
library(ggplot2)
library(readr)
library(UpSetR)

# Load in data
ar_data_df <- read_excel("TwoSamReadThreshold_ar_spiec_easi_gephi_network_graphdata_SeqTypeBodysiteUpdate.xlsx", sheet = "Node_Analysis")	# ar_spiec-easi_network.xlsx
aras_data_df <- read_excel("TwoSamReadThreshold_aras_spiec_easi_gephi_network_graphdata_SeqTypeBodysiteUpdate.xlsx", sheet = "Node_Analysis")	# aras_spiec-easi_network.xlsx
as_data_df <- read_excel("TwoSamReadThreshold_as_spiec_easi_gephi_network_graphdata_SeqTypeBodysiteUpdate.xlsx", sheet = "Node_Analysis")	# as_spiec-easi_network.xlsx
hc_data_df <- read_excel("TwoSamReadThreshold_hc_spiec_easi_gephi_network_graphdata_SeqTypeBodysiteUpdate.xlsx", sheet = "Node_Analysis")	# hc_spiec-easi_network.xlsx

# Select only the meaningful columns
ar_data_df_filtered <- ar_data_df %>%
  select(Node, ConnectionType, Taxa, bodysite, seqtype_bodysite, Count) %>%
  filter(ConnectionType == "Mix")

aras_data_df_filtered <- aras_data_df %>%
  select(Node, ConnectionType, Taxa, bodysite, seqtype_bodysite, Count) %>%
  filter(ConnectionType == "Mix")

as_data_df_filtered <- as_data_df %>%
  select(Node, ConnectionType, Taxa, bodysite, seqtype_bodysite, Count) %>%
  filter(ConnectionType == "Mix")

hc_data_df_filtered <- hc_data_df %>%
  select(Node, ConnectionType, Taxa, bodysite, seqtype_bodysite, Count) %>%
  filter(ConnectionType == "Mix")

# Create mapping to re-map the associated taxonomy
ar_data_df_taxonomy_map <- ar_data_df_filtered %>%
  select(Node, Taxa) %>%
  distinct()

aras_data_df_taxonomy_map <- aras_data_df_filtered %>%
  select(Node, Taxa) %>%
  distinct()

as_data_df_taxonomy_map <- as_data_df_filtered %>%
  select(Node, Taxa) %>%
  distinct()

hc_data_df_taxonomy_map <- hc_data_df_filtered %>%
  select(Node, Taxa) %>%
  distinct()

all_data_df_taxonomy_map <- rbind(ar_data_df_taxonomy_map, aras_data_df_taxonomy_map, as_data_df_taxonomy_map, hc_data_df_taxonomy_map) %>%
  unique()

print(nrow(ar_data_df_taxonomy_map))
print(nrow(aras_data_df_taxonomy_map))
print(nrow(as_data_df_taxonomy_map))
print(nrow(hc_data_df_taxonomy_map))

print(nrow(all_data_df_taxonomy_map))

# Aggregate and find highest degree ASVs
ar_data_df_ASV_aggregate <- ar_data_df_filtered %>%
  group_by(Node) %>%
  summarise(TotalCount = sum(Count), .groups = "drop") %>%
  arrange(desc(TotalCount)) %>%
  left_join(ar_data_df_taxonomy_map, by = "Node")

aras_data_df_ASV_aggregate <- aras_data_df_filtered %>%
  group_by(Node) %>%
  summarise(TotalCount = sum(Count), .groups = "drop") %>%
  arrange(desc(TotalCount)) %>%
  left_join(aras_data_df_taxonomy_map, by = "Node")

as_data_df_ASV_aggregate <- as_data_df_filtered %>%
  group_by(Node) %>%
  summarise(TotalCount = sum(Count), .groups = "drop") %>%
  arrange(desc(TotalCount)) %>%
  left_join(as_data_df_taxonomy_map, by = "Node")

hc_data_df_ASV_aggregate <- hc_data_df_filtered %>%
  group_by(Node) %>%
  summarise(TotalCount = sum(Count), .groups = "drop") %>%
  arrange(desc(TotalCount)) %>%
  left_join(hc_data_df_taxonomy_map, by = "Node")

# Plots the highest degree decending
degree_desc_plot <- function(df_aggregate, df_disease) {
  df_aggregate %>%
    arrange(desc(TotalCount)) %>%
    ggplot(aes(x = reorder(Node, -TotalCount), y = TotalCount)) +
    geom_point(size = 3, color = "steelblue") +
    labs(title = paste(df_disease, "Degrees/Connection Counting for each ASV"),
  x = "Node",
  y = "Sum Total Degrees") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

# Plot degree
ar_data_df_ASV_aggregate_degree_plot <- degree_desc_plot(ar_data_df_ASV_aggregate, "ar")
aras_data_df_ASV_aggregate_degree_plot <- degree_desc_plot(aras_data_df_ASV_aggregate, "aras")
as_data_df_ASV_aggregate_degree_plot <- degree_desc_plot(as_data_df_ASV_aggregate, "as")
hc_data_df_ASV_aggregate_degree_plot <- degree_desc_plot(hc_data_df_ASV_aggregate, "hc")

ar_data_df_ASV_aggregate_degree_plot

aras_data_df_ASV_aggregate_degree_plot

as_data_df_ASV_aggregate_degree_plot

hc_data_df_ASV_aggregate_degree_plot

# Find matching and mismatching ASVs between disease groups for the entire dataset
ar_data_df_ASV_aggregate_full <- ar_data_df_ASV_aggregate %>%
  pull(Node)

aras_data_df_ASV_aggregate_full <- aras_data_df_ASV_aggregate %>%
  pull(Node)


as_data_df_ASV_aggregate_full <- as_data_df_ASV_aggregate %>%
  pull(Node)


hc_data_df_ASV_aggregate_full <- hc_data_df_ASV_aggregate %>%
  pull(Node)


# Create a presence absence matrix with the full set
all_Nodes_full <- union(ar_data_df_ASV_aggregate_full, union(aras_data_df_ASV_aggregate_full, union(as_data_df_ASV_aggregate_full, hc_data_df_ASV_aggregate_full)))

presence_df_full <- tibble(Node = all_Nodes_full) %>%
  mutate(
    ar_list = as.integer(Node %in% ar_data_df_ASV_aggregate_full),
    aras_list = as.integer(Node %in% aras_data_df_ASV_aggregate_full),
    as_list = as.integer(Node %in% as_data_df_ASV_aggregate_full),
    hc_list = as.integer(Node %in% hc_data_df_ASV_aggregate_full)
  ) %>%
    mutate(AppearCount = ar_list + aras_list + as_list + hc_list) %>%
    left_join(all_data_df_taxonomy_map, by = "Node") # Add back in the taxonomy

write_csv(presence_df_full, "mixed_connection_ASV_Unique_presence_absence_matrix.csv")




"""# Upset Chart"""

# Prepare dataframe for upset
presence_absence_df <- presence_df_full %>%
  select(Node, ar_list, aras_list, as_list, hc_list) %>%
  as.data.frame()
rownames(presence_absence_df) <- presence_absence_df$Node
presence_absence_df$Node <- NULL


# Generate upset
options(repr.plot.width = 12, repr.plot.height = 12)
upset(
  presence_absence_df,
  nsets = ncol(presence_absence_df),
  nintersects = NA,
  order.by = "degree",
  text.scale = c(3, 3, 3, 2, 3, 4)
)




"""# Relative Mean Abundance Taxonomy Prevelance"""

# Packages
library(readxl)
library(tidyverse)
library(ggplot2)
library(dplyr)
library(tidyr)
library(rlang)
library(RColorBrewer)

# Load in data
merged_maindata_df <- read_excel("16S_ITS_Merged_mainData.xlsx", sheet = "Merge_mainData")
mixed_asv_only_df <- read_excel("16S_ITS_Merged_mainData.xlsx", sheet = "Mixed_ASV_Only")
metadata_cleaned_df <- read.csv("metadata_cleaned.csv")

abundance_stacked_waterfall_plot <- function(reads_data, meta_data, taxa_input, n_top) {
  id_cols_all <- c("ASV", "Taxonomy", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  id_cols <- intersect(names(reads_data), id_cols_all)
  sample_cols <- setdiff(names(reads_data), id_cols)

  # Convert to long format
  long_df <- reads_data %>%
    select(c(ASV, sample_cols)) %>%
    tidyr::pivot_longer(all_of(sample_cols), names_to = "Sample_Name", values_to = "Reads") %>%
    left_join(meta_data %>% select(Sample_Name, disease), by = "Sample_Name") %>%
    left_join(reads_data %>% select(ASV, all_of(taxa_input)), by = "ASV") %>%
    mutate(!!taxa_input := tidyr::replace_na(!!sym(taxa_input), "Unassigned"))

  # per-sample relative abundance
  long_df <- long_df %>%
    group_by(Sample_Name) %>%
    mutate(TotalReads = sum(Reads, na.rm = TRUE), RelAbund = ifelse(TotalReads > 0, Reads / TotalReads, 0)) %>%
    ungroup()

  # global mean per taxon to choose top N
  global_means <- long_df %>%
    group_by(!!sym(taxa_input)) %>%
    summarise(global_mean_rel = mean(RelAbund, na.rm = TRUE), .groups = "drop")

  top_taxa <- global_means %>%
    arrange(desc(global_mean_rel)) %>%
    slice_head(n = n_top) %>%
    pull(!!sym(taxa_input))

  # mean per disease x taxon
  mean_by_tax <- long_df %>%
    group_by(disease, !!sym(taxa_input)) %>%
    summarise(mean_rel = mean(RelAbund, na.rm = TRUE), .groups = "drop")

  # split into top vs non-top
  top_part <- mean_by_tax %>%
    filter((!!sym(taxa_input)) %in% top_taxa) %>%
    mutate(Taxon_Group = !!sym(taxa_input)) %>%
    select(disease, Taxon_Group, mean_rel)

  # compute "Other" strictly as residual so the stack sums to 100
  other_part <- mean_by_tax %>%
    filter(!(!!sym(taxa_input)) %in% top_taxa) %>%
    group_by(disease) %>%
    summarise(other_mean = sum(mean_rel, na.rm = TRUE), .groups = "drop") %>%
    mutate(Taxon_Group = "Other") %>%
    select(disease, Taxon_Group, mean_rel = other_mean)

  out <- bind_rows(top_part, other_part) %>%
    group_by(disease) %>%
    mutate(percent_abund = 100 * mean_rel / sum(mean_rel)) %>%
    ungroup()

  # order taxa by global abundance with "Other" last
  tax_order <- c(setdiff(top_taxa, "Other"), "Other")
  out$Taxon_Group <- factor(out$Taxon_Group, levels = tax_order)

  return(out)
}

mean_rel_df <- abundance_stacked_waterfall_plot(mixed_asv_only_df, metadata_cleaned_df, "Species", n_top=10)
options(repr.plot.width = 12, repr.plot.height = 8)

# Base plot
base_palette <- brewer.pal(12, "Paired")

my_colors <- colorRampPalette(base_palette)(length(unique(mean_rel_df$Taxon_Group)))

# Replace the color for "Other" with gray
names(my_colors) <- unique(mean_rel_df$Taxon_Group)
my_colors["Other"] <- "gray70"

ggplot(mean_rel_df, aes(x = disease, y = percent_abund, fill = Taxon_Group)) +
  geom_bar(stat = "identity", width = 0.80, color = "black") +
  scale_x_discrete(expand = c(0, 0)) +
  scale_fill_manual(values = my_colors) +
  labs(
    title = "Species",
    x = "Disease Group",
    y = "Mean Relative Abundance (%)",
    fill = "Taxon"
  ) +
  theme_minimal(base_size = 20) +
  theme(
    plot.title = element_text(
    face = "bold",
    size = 35,
    hjust = 0.5),
    legend.position = "bottom",
    legend.box = "horizontal",
    panel.grid.major.x = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(size = 15, face = "bold"),
    axis.title.x = element_blank(),
    axis.text.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold")
  )

# Generate a table with relative abundances
mean_rel_df %>%
  pivot_wider(
    names_from = disease,
    values_from = c(mean_rel, percent_abund)) %>%
      select(-c(mean_rel_AR, mean_rel_ARAS, mean_rel_AS, mean_rel_HC)) %>%
      rename_with(~ gsub("percent_abund_", "", .x))

mean_rel_df <- abundance_stacked_waterfall_plot(mixed_asv_only_df, metadata_cleaned_df, "Genus", n_top=10)
options(repr.plot.width = 12, repr.plot.height = 8)

# Base plot
base_palette <- brewer.pal(12, "Paired")

my_colors <- colorRampPalette(base_palette)(length(unique(mean_rel_df$Taxon_Group)))

# Replace the color for "Other" with gray
names(my_colors) <- unique(mean_rel_df$Taxon_Group)
my_colors["Other"] <- "gray70"

ggplot(mean_rel_df, aes(x = disease, y = percent_abund, fill = Taxon_Group)) +
  geom_bar(stat = "identity", width = 0.80, color = "black") +
  scale_x_discrete(expand = c(0, 0)) +
  scale_fill_manual(values = my_colors) +
  labs(
    title = "Genera",
    x = "Disease Group",
    y = "Mean Relative Abundance (%)",
    fill = "Taxon"
  ) +
  theme_minimal(base_size = 20) +
  theme(
    plot.title = element_text(
    face = "bold",
    size = 35,
    hjust = 0.5),
    legend.position = "bottom",
    legend.box = "horizontal",
    panel.grid.major.x = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(size = 15, face = "bold"),
    axis.title.x = element_blank(),
    axis.text.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold")
  )

# Generate a table with relative abundances
mean_rel_df %>%
  pivot_wider(
    names_from = disease,
    values_from = c(mean_rel, percent_abund)) %>%
      select(-c(mean_rel_AR, mean_rel_ARAS, mean_rel_AS, mean_rel_HC)) %>%
      rename_with(~ gsub("percent_abund_", "", .x))

mean_rel_df <- abundance_stacked_waterfall_plot(mixed_asv_only_df, metadata_cleaned_df, "Phylum", n_top=8)
options(repr.plot.width = 12, repr.plot.height = 8)

# Base plot
base_palette <- brewer.pal(12, "Paired")

my_colors <- colorRampPalette(base_palette)(length(unique(mean_rel_df$Taxon_Group)))

# Replace the color for "Other" with gray
names(my_colors) <- unique(mean_rel_df$Taxon_Group)
my_colors["Other"] <- "gray70"

ggplot(mean_rel_df, aes(x = disease, y = percent_abund, fill = Taxon_Group)) +
  geom_bar(stat = "identity", width = 0.80, color = "black") +
  scale_x_discrete(expand = c(0, 0)) +
  scale_fill_manual(values = my_colors) +
  labs(
    title = "Phyla",
    x = "Disease Group",
    y = "Mean Relative Abundance (%)",
    fill = "Taxon"
  ) +
  theme_minimal(base_size = 20) +
  theme(
    plot.title = element_text(
    face = "bold",
    size = 35,
    hjust = 0.5),
    legend.position = "bottom",
    legend.box = "horizontal",
    panel.grid.major.x = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(size = 14, face = "bold"),
    axis.title.x = element_blank(),
    axis.text.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold")
  )

# Generate a table with relative abundances
mean_rel_df %>%
  pivot_wider(
    names_from = disease,
    values_from = c(mean_rel, percent_abund)) %>%
      select(-c(mean_rel_AR, mean_rel_ARAS, mean_rel_AS, mean_rel_HC)) %>%
      rename_with(~ gsub("percent_abund_", "", .x))




"""# Statistial Difference between Disease and Control for High Abundance Taxonomies"""

library(tidyverse)
library(readxl)
library(dplyr)
library(tidyr)
library(purrr)
library(stringr)

taxa_stat_test <- function(merged_maindata_df,
                           metadata_cleaned_df,
                           taxa_vec,
                           level,
                           group_col = "disease",
                           sample_col = "Sample_Name") {

  clean_taxa_vec <- setdiff(unique(taxa_vec), "Other") # Remove the Other


  id_cols <- c("ASV", "Taxonomy", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  read_cols <- setdiff(names(merged_maindata_df), id_cols)

  # Long Format -- Relative Proportional Abundance
  df_long_rel <- merged_maindata_df %>%
    pivot_longer(
      cols      = all_of(read_cols),
      names_to  = sample_col,
      values_to = "Count"
    ) %>%
    left_join(metadata_cleaned_df, by = sample_col) %>%
    group_by(.data[[sample_col]]) %>%
    mutate(
      SampleTotal = sum(Count, na.rm = TRUE),
      RelAbund    = if_else(SampleTotal > 0, Count / SampleTotal, NA_real_)
    ) %>%
    ungroup()

  # Aggregate at taxa level
  taxa_abund <- df_long_rel %>%
    filter(.data[[level]] %in% clean_taxa_vec) %>%
    group_by(taxon = .data[[level]], .data[[sample_col]], .data[[group_col]]) %>%
    summarise(
      RelAbund_taxon = sum(RelAbund, na.rm = TRUE),
      .groups = "drop"
    )

  # Median and IQR Calculation, to identify direction of abundance
  disease_stats_long <- taxa_abund %>%
    group_by(taxon, .data[[group_col]]) %>%
    summarise(
      median = median(RelAbund_taxon, na.rm = TRUE),
      IQR    = IQR(RelAbund_taxon,    na.rm = TRUE),
      .groups = "drop"
    ) %>%
    pivot_longer(
      cols      = c(median, IQR),
      names_to  = "stat",
      values_to = "value"
    ) %>%
    mutate(stat = paste0(stat, "_", .data[[group_col]])) %>%
    select(taxon, stat, value)

  # Kruskal–Wallis
  kruskal_tbl <- taxa_abund %>%
    group_by(taxon) %>%
    summarise(
      kruskal_p   = kruskal.test(RelAbund_taxon ~ .data[[group_col]])$p.value,
      .groups     = "drop"
    ) %>%
    mutate(
      kruskal_adj = p.adjust(kruskal_p, method = "BH") # Benjamini-Hochberg FDR correction
    ) %>%
    select(taxon, kruskal_adj) %>%
    mutate(
      stat  = "Kruskal_adj",
      value = kruskal_adj
    ) %>%
    select(taxon, stat, value)

  # Pairwise Wilcoxon per taxon
  pairwise_long <- map_dfr(unique(taxa_abund$taxon), function(tax) {
    tmp <- filter(taxa_abund, taxon == !!tax)

    pw <- pairwise.wilcox.test(
      tmp$RelAbund_taxon,
      tmp[[group_col]],
      p.adjust.method = "BH" # Benjamini-Hochberg FDR correction
    )$p.value

    if (is.null(pw)) return(NULL)

    pw_df <- as.data.frame(as.table(pw), stringsAsFactors = FALSE)

    pw_df %>%
      filter(!is.na(Freq)) %>%
      rename(group1 = Var1, group2 = Var2, p_adj = Freq) %>%
      mutate(
        taxon = tax,
        stat  = paste0("Wilcoxon_", group1, "_vs_", group2),
        value = p_adj
      ) %>%
      select(taxon, stat, value)
  })

  # Export as Table
  all_stats_long <- bind_rows(
    disease_stats_long,
    kruskal_tbl,
    pairwise_long
  )

  summary_tab <- all_stats_long %>%
    pivot_wider(
      names_from  = stat,
      values_from = value
    ) %>%
    arrange(taxon)

  list(
    level       = level,
    taxa_abund  = taxa_abund,     # sample-level relabund per taxon
    summary_tab = summary_tab     # rows = taxa, cols = stats
  )
}

# Load in data
merged_maindata_df <- read_excel("/content/16S_ITS_Merged_mainData.xlsx", sheet = "Merge_mainData")
mixed_asv_only_df <- read_excel("/content/16S_ITS_Merged_mainData.xlsx", sheet = "Mixed_ASV_Only")
metadata_cleaned_df <- read.csv("/content/metadata_cleaned.csv")

# Run the Species Top Relative Abundance Code first to get mean_rel_df table list

print(unique(mean_rel_df$Taxon_Group))

res <- taxa_stat_test(
  merged_maindata_df,
  metadata_cleaned_df,
  taxa_vec = unique(mean_rel_df$Taxon_Group),
  level    = "Species"
)

res$summary_tab

write.csv(
  res$summary_tab,
  file = "Species_taxa_stat_results.csv",
  row.names = FALSE
)

# Run the Genus Top Relative Abundance Code first to get mean_rel_df table list

print(unique(mean_rel_df$Taxon_Group))

res <- taxa_stat_test(
  merged_maindata_df,
  metadata_cleaned_df,
  taxa_vec = unique(mean_rel_df$Taxon_Group),
  level    = "Genus"
)

res$summary_tab

write.csv(
  res$summary_tab,
  file = "Genus_taxa_stat_results.csv",
  row.names = FALSE
)

# Run the Phylum Top Relative Abundance Code first to get mean_rel_df table list

print(unique(mean_rel_df$Taxon_Group))

res <- taxa_stat_test(
  merged_maindata_df,
  metadata_cleaned_df,
  taxa_vec = unique(mean_rel_df$Taxon_Group),
  level    = "Phylum"
)

res$summary_tab

write.csv(
  res$summary_tab,
  file = "Phylum_taxa_stat_results.csv",
  row.names = FALSE
)