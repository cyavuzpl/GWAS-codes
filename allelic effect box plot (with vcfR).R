# Load necessary libraries
library(vcfR)  # For handling VCF (Variant Call Format) files
library(dplyr) # For data manipulation tasks like filtering and summarizing
library(ggplot2) # For creating visualizations, like boxplots

# Read the VCF file containing genetic data
vcf <- read.vcfR("xxxxx.vcf")  # Replace "xxxxx.vcf" with your file path

# Define the chromosome and position of interest in the dataset
chr <- "5"  # Specify the chromosome (ensure this matches the format in your VCF file, e.g., "chr5" vs "5")
pos <- 169387367  # Specify the position of the SNP of interest

# Extract data for the SNP at the defined chromosome and position
vcf_filtered <- vcf[vcf@fix[, "CHROM"] == chr & as.numeric(vcf@fix[, "POS"]) == pos, ]

# Convert the filtered VCF data to a tidy format, focusing on genotype information
snp_tidy <- vcfR2tidy(vcf_filtered)$gt

# View the structure of the SNP data
head(snp_tidy)

# Load phenotype data (e.g., trait measurements like oil content)
phenotypes <- read.csv("xxxxx.phe", header = TRUE, sep = "")  # Ensure correct delimiter (e.g., ',' or '\t')
colnames(phenotypes)[1] <- "Indiv"  # Rename the first column to "Indiv" to match with genotype data

# Merge the genotype data (snp_tidy) with the phenotype data (phenotypes) by individual ID
gwas_data <- merge(snp_tidy, phenotypes, by = "Indiv")

# Clean and modify the genotype data
gwas_data <- gwas_data %>%
  mutate(
    gt_GT_alleles = case_when(
      gt_GT_alleles == "A|A" ~ "A",  # If genotype is homozygous "A|A", label as "A"
      gt_GT_alleles == "G|G" ~ "G",  # If genotype is homozygous "G|G", label as "G"
      TRUE ~ gt_GT_alleles  # Keep other genotypes unchanged
    )
  ) %>%
  filter(!gt_GT_alleles %in% c("A|G", "G|A"))  # Remove heterozygous genotypes ("A|G" and "G|A")

# Check the cleaned genotype-phenotype data
head(gwas_data)

# Calculate the frequency of each allele (e.g., "A" and "G")
counts <- gwas_data %>%
  group_by(gt_GT_alleles) %>%
  summarise(n = n())

# Prepare data for visualization: select relevant columns and filter out unwanted genotypes
datanew <- gwas_data %>%
  filter(!gt_GT_alleles %in% c("G|A", "A|G", "G/A", "A/G")) %>%  # Remove heterozygous genotypes
  select(gt_GT_alleles, Oil) %>%  # Keep genotype and phenotype (Oil) columns
  rename(Allele = gt_GT_alleles, Value = Oil) %>%  # Rename columns for clarity
  mutate(Trait = "Oil")  # Add a "Trait" column to indicate what trait is being analyzed

# Calculate the average value of the trait (e.g., oil content) for each allele
avg_oil_content <- datanew %>%
  group_by(Allele) %>%
  summarise(avg_value = mean(Value, na.rm = TRUE))  # Compute the average while ignoring missing values

# Calculate the percent change of each allele's trait value compared to the reference allele (assumed to be "A/A")
percent_change <- avg_oil_content %>%
  mutate(percent_change = (avg_value / avg_value[Allele == "A/A"] - 1) * 100)  # Percent change relative to "A/A"
percent_change$percent_change[percent_change$Allele == "A/A"] <- NA  # Set percent change for "A/A" to NA

# Define y-axis positioning for text labels to avoid overlap with the plot
y_max <- max(datanew$Value, na.rm = TRUE)  # Get the highest value in the trait data
y_offset <- 0.1 * y_max  # Set a small offset to place labels slightly above the highest data point
y_offset_percent <- 0.1 * y_max  # Additional offset for percent change labels

# Calculate the count of individuals for each allele
count_data <- datanew %>%
  group_by(Allele) %>%
  summarise(n = n())

# Create a visualization: a boxplot with jittered points and text labels for counts and percent change
ggplot(datanew, aes(x = Allele, y = Value, fill = Allele)) +  # Plot trait values by allele
  geom_boxplot(trim = FALSE, alpha = 0.7, outlier.shape = NA) +  # Create boxplot (hide whiskers and outliers)
  geom_jitter(aes(color = Allele), width = 0.1, size = 3, shape = 16, alpha = 0.6) +  # Add jittered points for each individual
  geom_text(data = count_data, 
            aes(x = Allele, y = y_max + y_offset, label = paste("N =", n)), 
            inherit.aes = FALSE, color = "black", size = 4, vjust = 0) +  # Add count labels above each boxplot
  geom_text(data = percent_change, 
            aes(x = Allele, y = y_max + y_offset + y_offset_percent, 
                label = ifelse(!is.na(percent_change), paste(round(percent_change, 1), "%"), "")), 
            color = "blue", size = 4, vjust = 0) +  # Add percent change labels (excluding reference "A/A")
  scale_fill_manual(values = c("brown", "grey", "blue", "purple")) +  # Specify custom colors for alleles
  theme_minimal(base_size = 15) +  # Use a minimal theme for the plot
  labs(y = "Trait (%)") +  # Label the y-axis (adjusted for your trait)
  theme(
    panel.background = element_blank(),  # Remove background grid
    panel.grid = element_blank(),  # Remove gridlines
    axis.line = element_line(color = "black")  # Add a border around the axes
  ) +
  xlab("xxxxxxx")  # Label for the x-axis (replace "xxxxxxx" with your desired label)
