# Required libraries

library(phyloseq)
library(ggplot2)
library(dplyr)
library(vegan)
library(viridis)
library(microViz)
library(stringr)
library(forcats)

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

#### Observed Richness and Diversity ####

# Get OTU table
otu <- as(otu_table(phylo_bacteria), "matrix")

# Count taxa present in each sample
if (taxa_are_rows(phylo_bacteria)) {
  observed <- colSums(otu > 0)
} else {
  observed <- rowSums(otu > 0)
}

# Create richness data frame
richness <- data.frame(
  SampleID = names(observed),
  Observed = as.numeric(observed))

# Combine richness data with metadata
richness <- richness %>%
  left_join(
    bacteria_metadata %>%
      mutate(SampleID = rownames(bacteria_metadata)) %>%
      select(SampleID, Species, Ear_Tag_ID, Sample_Type, pair_label),
    by = "SampleID"
  )

# View richness data frame
richness

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


# Plot Bacterial Richness
bacteria_richness <- 
  ggplot(richness, aes(x = Sample_Type, y = Observed, group = pair_label)) +
  geom_line(linewidth = 1, color = "black") +
  geom_point(aes(fill = Sample_Type, shape = Species), size = 5, color = "black", stroke = 1, alpha = 1) + 
  geom_point(data = subset(richness, Ear_Tag_ID == 11 & Species == "PE-PO"),
             aes(x = Sample_Type, y = Observed, shape = Species), size = 5, fill = NA, color = "red", stroke = 1.5) +
  labs(x = NULL, y = "Bacterial Richness") +
  scale_x_discrete(limits = c("SCAT", "SWAB"), labels = c("SCAT" = "Gut","SWAB" = "Skin")) +
  scale_fill_manual(values = custom_fills, labels = c("SCAT" = "Gut", "SWAB" = "Skin")) + 
  scale_shape_manual(values = c("RE-HU" = 21, "NE-FL" = 22, "PE-PO" = 24),
    labels = c( "RE-HU" = "<i>R. humulis</i>","NE-FL" = "<i>N. floridana</i>", "PE-PO" = "<i>P. polionotus</i>")) + 
  theme_classic() +
  theme(
    axis.title.y = element_text(size = 20, margin = margin(r = 15)),
    axis.text.y = element_text(size = 20, color = "grey35"),
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 20),
    axis.line = element_line(color = "grey35"),
    plot.title = element_text(size = 20, hjust = 0.5),
    legend.position = c(0.05, 1.0),
    legend.justification = c("left", "top"), 
    legend.direction = "horizontal",
    legend.text = ggtext::element_markdown(size = 12),
    legend.title = element_blank(),
    legend.background = element_rect(fill = "white"),
    legend.key.spacing.y = unit(0.5, "cm")) +
  guides(
    color = "none", fill = "none", shape = guide_legend(title = NULL, order = 2, ncol = 1,
      override.aes = list(fill = "white", color = "black", size = 5)))

# View plot
print(bacteria_richness)

# Calculate Bacterial Shannon Diversity
shannon <- diversity(t(OTU), index = "shannon")

# View shannon data frame
shannon

# Add Shannon index to metadata 
bacteria_metadata$Shannon <- shannon

# Plot shannon
bacteria_shannon <- 
  ggplot(
    bacteria_metadata, 
    aes(
      x = Sample_Type, 
      y = Shannon, 
      group = pair_label)) +
  geom_line(linewidth = 1, color = "black") +
  geom_point(aes(fill = Sample_Type, shape = Species), size = 5, color = "black", stroke = 1, alpha = 1) + 
  geom_point(data = subset(bacteria_metadata, Ear_Tag_ID == 11 & Species == "PE-PO"), 
             aes(x = Sample_Type, y = Shannon, shape = Species), size = 5, fill = NA, color = "red", stroke = 1.5) +
  labs(x = NULL, y = "Bacterial Diversity") +
  scale_x_discrete(limits = c("SCAT", "SWAB"), labels = c("SCAT" = "Gut","SWAB" = "Skin")) +
  scale_fill_manual(values = custom_fills) +
  scale_shape_manual(values = c("RE-HU" = 21,"NE-FL" = 22,"PE-PO" = 24)) + 
  scale_y_continuous(limits = c(1, 5)) + 
  theme_classic() +
  theme(
    axis.title.y = element_text(size = 20, margin = margin(r = 15)),
    axis.text.y = element_text(size = 20, color = "grey35"),
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 20),
    axis.line = element_line(color = "grey35"),
    plot.title = element_text(size = 20, hjust = 0.5),
    legend.position = "none")

# View plot
print(bacteria_shannon)


#### Relative Abundance by Family ####

# Step 1: Melt phyloseq object
melt_bacteria <- psmelt(phylo_bacteria) %>%
  mutate(Sample_Type = as.factor(Sample_Type))

# Step 2: Identify top 10 taxa per group
top_taxa_bacteria <- melt_bacteria %>%
  group_by(Sample_Type, f) %>%
  summarise(mean_abundance = mean(Abundance), .groups = "drop") %>%
  group_by(Sample_Type) %>%
  slice_max(mean_abundance, n = 10, with_ties = FALSE) %>%
  ungroup()

# Step 3: Label taxa as "Other" if not in top 10 for that group
melt_bacteria <- melt_bacteria %>%
  left_join(
    top_taxa_bacteria,
    by = c("Sample_Type", "f")
  ) %>%
  mutate(
    f = ifelse(is.na(mean_abundance), "Other", f)
  ) %>%
  dplyr::select(-mean_abundance)

# Step 4: Summarize + calculate relative abundance
relative_abundance_bacteria <- melt_bacteria %>%
  group_by(Sample_Type, f) %>%
  summarise(Abundance = sum(Abundance), .groups = "drop") %>%
  group_by(Sample_Type) %>%
  mutate(Rel_Abundance = Abundance / sum(Abundance)) %>%
  ungroup()

# Create a map from abbreviations to the full taxonomic ranks
rank_map <- c(
  "k" = "kingdom",
  "p" = "phylum",
  "c" = "class",
  "o" = "order",
  "f" = "family",
  "g" = "genus",
  "s" = "species"
)

# Create a clean taxonomic name 
relative_abundance_bacteria <- relative_abundance_bacteria %>%
  mutate(
    rank_letter = str_extract(f, " ([a-z])$") %>% str_trim(),
    taxon_name = str_remove(f, " ([a-z])$"),
    rank_full = rank_map[rank_letter],
    f_clean = case_when(
      f == "Other" ~ "Other",
      !is.na(rank_full) ~ paste0(taxon_name, " (", rank_full, ")"),
      TRUE ~ f
    )
  )

# Add cleaned taxonomic name and relabel columns
Final_relative_abundance_bacteria <- relative_abundance_bacteria %>%
  arrange(
    factor(Sample_Type, levels = c("SWAB", "SCAT")),
    desc(Rel_Abundance)
  ) %>%
  transmute(
    Microbiome = recode(Sample_Type, SWAB = "Skin", SCAT = "Gut"),
    Kingdom = "Bacteria",
    `Taxonomic Classification` = f_clean,
    `Relative Abundance` = round(Rel_Abundance, 4)
  )


# Save as table
write.csv(Final_relative_abundance_bacteria, "Bacteria_Relative_Abundance.csv", row.names = FALSE)
