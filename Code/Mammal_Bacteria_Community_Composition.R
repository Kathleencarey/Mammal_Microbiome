# Required libraries
library(phyloseq)
library(ggplot2)
library(dplyr)
library(vegan)
library(viridis)
library(DESeq2)
library(microViz)
library(scales)
library(ggtext)

# Load data

# Load bacteria (16S) counts 
bacteria <- read.csv("Bacteria_Reads.csv", check.names = FALSE)

# Load metadata 
bacteria_metadata <- read.csv("Mammal_Capture_Data.csv")

# Create row names of unique samples
rownames(bacteria_metadata) <- bacteria_metadata$Sample_Label

# Make these factors for DESeq2
bacteria_metadata$pair_label <- as.factor(bacteria_metadata$pair_label)
bacteria_metadata$Sample_Type <- as.factor(bacteria_metadata$Sample_Type)

# Remove taxonomy information and create OTU matrix
otu_mat <- as.matrix(bacteria[, 9:ncol(bacteria)])

# Create taxa matrix
taxa_mat <- as.matrix(bacteria[, 1:8])

# Create phyloseq objects
OTU <- otu_table(otu_mat, taxa_are_rows = TRUE)
TAX <- tax_table(taxa_mat)

# Convert metadata into a sample metadata object
sampledata <- sample_data(bacteria_metadata)

# Create phlyoseq objust
physeq <- phyloseq(OTU, TAX, sampledata)

# Set seed for analysis reproducibility
set.seed(123) 

# Rarify OTU table
OTU_rar = rarefy_even_depth(OTU, rngseed = TRUE, replace = FALSE) 

# Convert OTU table into an R data frame
phylo_bacteria = data.frame(otu_table(OTU_rar))

# Create phyloseq object
phylo_bacteria <- phyloseq(OTU_rar, TAX, sampledata) 

# Standardize taxonomic classifications by removing undetected taxonomy
phylo_bacteria <- phylo_bacteria %>% tax_fix() %>% 
  phyloseq_validate(remove_undetected = TRUE)

# Calculate Bray-Curtis distance
rodent_bc <- phyloseq::distance(phylo_bacteria, method = "bray") 

# Extract the sample metadata from the phyloseq object + convert it into a R data frame
rodent_meta <- data.frame(sample_data(phylo_bacteria))

# Set seed for reproducibility
set.seed(123) 

# Run perMANOVA analysis with strata for each pair (rodent individual)
adonis2(rodent_bc ~ Sample_Type, data = bacteria_metadata, strata = bacteria_metadata$pair_label, permutations = 999)

# Run test for multivariate dispersion
bd <- betadisper(rodent_bc, bacteria_metadata$Sample_Type)

# View results
print(bd)

# Run permutation test
permutest(bd, permutations = 999)

# Plot NMDS

# Run NMDS on the distance matrix
nmds_result_bacteria <- metaMDS(rodent_bc, k = 2, trymax = 100)

# Check stress level
nmds_result_bacteria$stress

# Extract NMDS scores
nmds_scores_bacteria <- as.data.frame(scores(nmds_result_bacteria, display = "sites"))

# Merge NMDS scores with metadata
nmds_scores_bacteria$SampleID <- rownames(nmds_scores_bacteria)
bacteria_metadata$SampleID <- rownames(bacteria_metadata)

# Merge NMDS scores with metadata
merged_scores_bacteria <- left_join(nmds_scores_bacteria, bacteria_metadata, by = "SampleID")

# Ellipse function
veganCovEllipse <- function(cov, center = c(0, 0), scale = 1, npoints = 100) {
  theta <- seq(0, 2 * pi, length.out = npoints)
  circle <- cbind(cos(theta), sin(theta))
  t(center + scale * t(circle %*% chol(cov)))
}

# # Create a scale factor of 95% confidence ellipses
scale_factor <- sqrt(qchisq(0.95, df = 2))  

# Create ellipses
df_ell_bacteria <- data.frame()
for (g in unique(merged_scores_bacteria$Sample_Type)) {
  group_data <- merged_scores_bacteria[merged_scores_bacteria$Sample_Type == g, c("NMDS1", "NMDS2")]
  if (nrow(group_data) >= 3) {
    cov_matrix <- cov(group_data)
    center <- colMeans(group_data)
    ellipse_coords <- veganCovEllipse(cov = cov_matrix, center = center, scale = scale_factor)
    df_ell_bacteria <- rbind(df_ell_bacteria, data.frame(NMDS1 = ellipse_coords[,1],
                                                         NMDS2 = ellipse_coords[,2],
                                                         Sample_Type = as.factor(g)))
  }
}

# Create custom colors and fill for plotting
viridis_colors <- viridis(100)
dark_blue <- viridis_colors[40]  
dark_yellow <- viridis_colors[100] 

light_blue <- alpha(dark_blue, 0.6)
light_yellow <- alpha(dark_yellow, 0.6)


custom_colors <- c("SWAB" = dark_blue, 
                   "SCAT" = dark_yellow)

custom_fills <- c("SWAB" = light_blue, 
                  "SCAT" = light_yellow)    

# Final NMDS Plot

Method_nmds_bacteria <- ggplot(merged_scores_bacteria, aes(x = NMDS1, y = NMDS2, color = Sample_Type, fill = Sample_Type)) +
  geom_polygon(data = df_ell_bacteria, 
               aes(x = NMDS1, y = NMDS2, group = Sample_Type, fill = Sample_Type, color = Sample_Type),
               alpha = 0.1, linewidth = 1.5, show.legend = FALSE) +
  geom_path(data = df_ell_bacteria,
            aes(x = NMDS1, y = NMDS2, group = Sample_Type, color = Sample_Type),
            linetype = 2, linewidth = 1.2, show.legend = FALSE) +
  geom_line(aes(group = pair_label),
            linewidth = 0.9, alpha = 1.0, color = "black") +
  geom_point(aes(fill = Sample_Type, shape = Species), size = 5, color = "black", stroke = 1, alpha = 1) +
  geom_point(data = subset(merged_scores_bacteria, Ear_Tag_ID == 11 & Species == "PE-PO"),
    aes(x = NMDS1, y = NMDS2, shape = Species), size = 5, fill = NA, color = "red", stroke = 1.5) +
  theme_classic() +
  theme(
    axis.text.y.left = element_text(size = 20, colour = "grey35"),
    axis.text.x = element_text(size = 20, colour = "grey35"),
    axis.line = element_line(colour = "grey35"),
    axis.title.x = element_text(size = 20, margin = margin(t = 15)),
    axis.title.y = element_text(size = 20, margin = margin(r = 15)),
    legend.position = c(0.05, 1.0),
    legend.justification = c("left", "top"), 
    legend.direction = "horizontal",
    legend.text = ggtext::element_markdown(size = 12),
    legend.title = element_blank(),
    legend.background = element_rect(fill = "white"),
    legend.key.spacing.y = unit(0.5, "cm")) +
  labs(x = "NMDS1", y = "NMDS2") +
  scale_y_continuous(limits = c(-0.6, 1.0), breaks = c(-1.0, -0.5, 0, 0.5, 1.0)) +
  scale_color_manual(values = custom_colors, labels = c("SCAT" = "Gut", "SWAB" = "Skin")) +
  scale_fill_manual(values = custom_fills, labels = c("SCAT" = "Gut", "SWAB" = "Skin")) +
  scale_shape_manual(values = c("RE-HU" = 21, "NE-FL" = 22, "PE-PO" = 24),
    labels = c( "RE-HU" = "<i>R. humulis</i>", "NE-FL" = "<i>N. floridana</i>", "PE-PO" = "<i>P. polionotus</i>")) + 
  guides(color = "none", fill = guide_legend(title = "Sample Type", order = 1,
      override.aes = list(shape = 21, color = "black", size = 5)),
    shape = guide_legend(title = "Species", order = 2, override.aes = list(fill = "white", color = "black", size = 5)))

# View Plot
Method_nmds_bacteria

#### Differential Abundance Analysis ####

# Remove extremely rare taxa (keep taxa present in at least 10% of samples)
phylo_bacteria_filtered <- prune_taxa(
  rowSums(otu_table(phylo_bacteria) == 0) < ncol(otu_table(phylo_bacteria)) * 0.9, 
  phylo_bacteria)

# Confirm that all taxa with zero counts are removed
zero_count_taxa <- rowSums(otu_table(phylo_bacteria_filtered) == 0) == ncol(otu_table(phylo_bacteria_filtered))

sum(zero_count_taxa) 

# Verify sample_type and pair_label are in phylo
colnames(sample_data(phylo_bacteria_filtered))

#add a plus for pair label
rodent_deseq2_obj <- phyloseq_to_deseq2(phylo_bacteria_filtered, ~ Sample_Type + pair_label)

# Run Wald Test
rodent_deseq2 <- DESeq(rodent_deseq2_obj, test = "Wald", fitType = "parametric", sfType = "poscounts")

# Define alpha for significance
alpha <- 0.05

#Make sample_type a factor for analysis
bacteria_metadata$Sample_Type <- factor(bacteria_metadata$Sample_Type)

# Create comparison list
contrast_list <- list(
  "SCAT vs SWAB" = c("Sample_Type", "SCAT", "SWAB"))


# Initialize a list to store results for each contrast
results_list <- list()

# Loop through each contrast, extract results, and add comparison information
for (contrast_name in names(contrast_list)) {
  contrast <- contrast_list[[contrast_name]]
  
  # Extract results for the specific contrast
  deseq_res <- results(rodent_deseq2, contrast = contrast, cooksCutoff = FALSE)
  
  # Filter significant results
  sigtab <- deseq_res[which(deseq_res$pvalue < alpha), ]
  
  # Convert DESeq2 results to data frame
  sigtab_df <- as.data.frame(sigtab)
  
  # Add comparison information
  sigtab_df$comparison <- contrast_name
  
  # Extract taxonomy table and ensure it matches with the results
  tax_data <- tax_table(phylo_bacteria_filtered)
  tax_data_df <- as.data.frame(tax_data)
  
  # Add identifier column to merge taxonomy data with DESeq2 results
  sigtab_df$id <- rownames(sigtab_df)
  tax_data_df$id <- rownames(tax_data_df)
  
  # Combine results with taxonomy data
  combined_df <- left_join(sigtab_df, tax_data_df, by = "id")
  
  # Add to results list
  results_list[[contrast_name]] <- combined_df
}

# Combine all results into one data frame
Bacteria_Differential_Abundance <- bind_rows(results_list)

# Inspect the combined data
head(Bacteria_Differential_Abundance)

# filter by significant padj (adjusted p-value)
significant_bacteria <- Bacteria_Differential_Abundance %>%
  filter(padj < 0.05) %>% 
  arrange(desc(log2FoldChange))

# View results
print(significant_bacteria)

#### Review Serratia found in OTUs ####

Taxa_dataframe <- as.data.frame(tax_table(phylo_bacteria))
Taxa_dataframe$TaxonID <- rownames(Taxa_dataframe)
Serratia_taxa <- Taxa_dataframe[Taxa_dataframe$g == "Serratia", ]
print(Serratia_taxa)

Serratia_ids <- Serratia_taxa$TaxonID

Serratia_otu <- as.data.frame(
  otu_table(phylo_bacteria)[Serratia_ids, , drop = FALSE]
)

