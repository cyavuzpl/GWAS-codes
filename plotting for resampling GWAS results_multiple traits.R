library(tidyverse)
library(data.table)
library(RColorBrewer)
library(ggrepel)

gwas.dat <- read.csv("combined_signals.csv")

# Define a custom palette with brown, blue, red, and similar colors
custom_palette <- c("brown", "blue", "red", "green", "purple")

# Ensure there are enough colors for all unique traits
traits <- unique(gwas.dat$basename)
num_traits <- length(traits)

# Print number of unique traits
cat("Number of unique traits: ", num_traits, "\n")

# Adjust the color palette to match the number of traits
if (num_traits > length(custom_palette)) {
  custom_palette <- rep(custom_palette, length.out = num_traits)  # Repeat the colors if needed
}

# Assign unique colors to each trait
gwas.dat$color <- factor(gwas.dat$basename, levels = traits)
gwas.dat$color <- custom_palette[as.numeric(gwas.dat$color)]

# Ensure that gwas.dat$color has the correct length
cat("Length of color vector: ", length(gwas.dat$color), "\n")

colnames(gwas.dat)[2] <- "CHROM"
gwas.dat <- as.data.frame(gwas.dat)

nCHR <- length(unique(gwas.dat$CHROM))
gwas.dat$BPcum <- NA
s <- 0
nbp <- numeric(nCHR)  # Initialize with proper length

# Ensure CHROM is treated as numeric, if possible
gwas.dat$CHROM <- as.numeric(as.character(gwas.dat$CHROM))

# Remove any rows with NA in CHROM or POS before processing
gwas.dat <- gwas.dat[!is.na(gwas.dat$CHROM) & !is.na(gwas.dat$POS), ]

# Add a large space between chromosomes by adjusting the cumulative position
for (i in sort(unique(gwas.dat$CHROM))) {
  chr_subset <- gwas.dat[gwas.dat$CHROM == i, ]
  
  if (nrow(chr_subset) > 0) {
    nbp[i] <- max(chr_subset$POS, na.rm = TRUE)
    gwas.dat[gwas.dat$CHROM == i, "BPcum"] <- chr_subset$POS + s
    s <- s + nbp[i] + 50000000  # Increase the gap to 50,000,000 between chromosomes
  }
}

axis.set <- gwas.dat %>% 
  group_by(CHROM) %>% 
  summarize(center = (max(BPcum) + min(BPcum)) / 2, maxBP=max(BPcum))

# Calculate limits for both x and y axes to ensure they are the same
x_lim <- range(gwas.dat$BPcum, na.rm = TRUE)
y_lim <- c(0, 1)  # Define the y axis limit



F1 <- ggplot(gwas.dat, aes(x = BPcum, y = support)) +
  geom_point(alpha = 0.75, size = 2, shape = 21, aes(fill = basename)) +  # Use fill for the data points
  geom_label_repel(
    data = subset(gwas.dat, support > 0.1),
    aes(label = SNP),  # Label with SNP names
    nudge_x = -100, 
    nudge_y = 0,    
    size = 2, 
    box.padding = 0.5, 
    point.padding = 1.0,
    force = 90,  
    segment.size = 0.2,  
    segment.color = "grey60",  
    segment.linetype = "dashed",  
    direction = "y",
    max.overlaps = 25 
  ) +
  geom_hline(yintercept = c(0.05, 0.1), color = c("grey40", "darkred"), linetype = "dashed") +
  scale_x_continuous(
    breaks = axis.set$center,
    labels = paste("chr", axis.set$CHROM, sep = ""),
    expand = c(0.05, 0.05),
    limits = c(0, max(gwas.dat$BPcum) + 50000000)
  ) +
  scale_y_continuous(expand = c(0, 0), limits = y_lim) +
  labs(x = NULL, y = "RMIP") +
  scale_fill_manual(values = custom_palette, name = "Traits") +  # Use fill for the legend
  scale_color_manual(values = custom_palette, guide = "none") +  # Suppress color legend
  theme_classic(base_size = 15, base_family = "NimbusSan") +
  theme(
    legend.position = "top",  # Adjust to your preference
    panel.border = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text.x = element_text(size = 12, vjust = 0.5, family = "NimbusSan"),
    axis.text.y = element_text(size = 12, family = "NimbusSan")
  )

F1

ggsave("resamplingGWAS.png", width = 15,   height = 8)

