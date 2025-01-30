### Rachel Gray
### January 2025
### Visualing Nucleotide Diversity and Heterozygosity Results
#############################################################
png("PIXY_hybrid_barplot.png", width = 1754, height = 1240, res = 300)
windows(width = 8.27, height = 11.69)  # test A4 sizing 
par(mfrow = c(2,2), mar = c(4, 4, 4,4), oma = c(3, 3, 3, 3)) # 4 panel figure, adjusting margins
### PIXY (Nucloetide diversity)
## hybrids only 
dat <- read.csv("PIXY_hybrid_barplot.csv", header=FALSE)
# check
head(dat)
# give headings
colnames(dat) <- c("Sample", "Category", "Color", "Value")
# check
head(dat)
# make factors so the order stays same
dat$Sample <- factor(dat$Sample, levels = unique(dat$Sample))
# turn off scientific notation
options(scipen = 999) 
barplot(
  height = dat$Value,               # Heights of the bars
  names.arg = dat$Sample,           # Labels for the x-axis
  col = dat$Color,                  # Colours for the bars
  las = 2,# Rotate x-axis labels to vertical
  xaxt = "n",# turn off the xaxis labels 
  axes =TRUE, # but keep the line
  cex.names = 0.5,    # Adjust size of x-axis labels
  cex=0.75, # Adjust size of y-axis labels
  main = "",  # Title of the plot
  xlab = "",               # Label for the x-axis
  ylab = "Nucleotide Diversity", # Label for the y-axis
  border = NA,   # no border for each bar
  ylim = c(0, 0.0012) # limit for y axis
)
abline(h = 0, col = "black", lwd = 1.5)  # Add a horizontal axis line at y = 0
abline(h = 0.000620379, col = "black",lwd=1.5,lty = 2)  ##Floreana Mean
abline(h = 0.000303346, col = "black",lwd=1.5,lty = 2)  ##Espanola Mean
abline(h = 0.000454094, col = "black",lwd=1.5,lty = 2)  ##VW Mean
dev.off()
### VCFTOOLS het
## hybrids only 
png("VCFTOOLS_het_hybrid_BARPLOTdat.png", width = 1754, height = 1240, res = 300)
at <- read.csv("VCFTOOLS_het_hybrid_BARPLOTdat.csv", header=FALSE)
# check
head(dat)
# give headings
colnames(dat) <- c("Sample", "Category", "Color", "Value")
# check
head(dat)
# make factors so the order stays same
dat$Sample <- factor(dat$Sample, levels = unique(dat$Sample))
# make plot
barplot(
  height = dat$Value,               # Heights of the bars
  names.arg = dat$Sample,           # Labels for the x-axis
  col = dat$Color,                  # Colours for the bars
  las = 2,                       # Rotate x-axis labels to vertical
  cex.names = 0.5,    # Adjust size of x-axis labels
  cex=0.75, # Adjust size of y-axis labels
  xaxt = "n",# turn off the xaxis labels 
  axes =TRUE, # but keep the line
  main = "",  # Title of the plot
  xlab = "",               # Label for the x-axis
  ylab = "Heterozygosity", # Label for the y-axis
  border = NA,   # no border for each bar
  ylim = c(0, 0.0012) # limit for y axis
)
abline(h = 0, col = "black", lwd = 1.5)  # Add a horizontal axis line at y = 0
abline(h = 0.000619916, col = "black",lwd=1.5,lty = 2)  ##Floreana Mean
abline(h = 0.000267633, col = "black",lwd=1.5,lty = 2)  ##Espanola Mean
abline(h = 0.000391828, col = "black",lwd=1.5,lty = 2)  ##VW Mean
dev.off()
## PIXY All lineages box plot
png("PIXY_all_species_BOXPLOTdat.png", width = 1754, height = 1240, res = 300)
dat <- read.csv("PIXY_all_species_BOXPLOTdat.csv")
head(dat)

split_data <- split(dat$PIXY, dat$Group) 
# Calculate mean for each group 
mean_het <- sapply(split_data, mean) 
# Create a data frame for the results 
mean_het_df <- data.frame(Group = names(mean_het), Mean_Het = mean_het) 
# Print the results 
print(mean_het_df)

# Create the boxplot 
boxplot(PIXY ~ Group, data=dat, main = "", xlab = "", ylab = "Nucleotide Diversity", xaxt = "n",  ylim=c(0, 0.001), cex.axis =0.75,las = 1)
# Adjust x-axis labels 
groups <- levels(factor(dat$Group))
par(cex.axis=0.8) 
# Reduce text size 
axis(side = 1, at = 1:length(groups), labels = groups, las = 2)
sample_sizes <- table(dat$Group)

boxplot_stats <- boxplot(PIXY ~ Group, data = dat, plot = FALSE) 

# Get stats without plotting again 
for (i in 1:length(groups)) { 
  group_name <- groups[i] 
  sample_size <- sample_sizes[group_name] # Position text above the boxplot whiskers 
  text(x = i, y = boxplot_stats$stats[5, i] + 0.0001, labels = paste0("n=", sample_size), col = "darkred", cex = 0.8) }
dev.off()

## VCFTOOLS het all lineages box plot
png("VCFTOOLS_het_all_species_BOXPLOTdat.png", width = 1754, height = 1240, res = 300)
dat <- read.csv("VCFTOOLS_het_all_species_BOXPLOTdat.csv")
head(dat)

split_data <- split(dat$VCFTOOLS, dat$Group) 
# Calculate mean for each group 
mean_het <- sapply(split_data, mean) 
# Create a data frame for the results 
mean_het_df <- data.frame(Group = names(mean_het), Mean_Het = mean_het) 
# Print the results 
print(mean_het_df)

# Create the boxplot 
boxplot(VCFTOOLS ~ Group, data=dat, main = "", xlab = "", ylab = "Heterozygosity", xaxt = "n", ylim=c(0, 0.001),cex.axis =0.75, las = 1)
# Adjust x-axis labels 
groups <- levels(factor(dat$Group))
par(cex.axis=0.8) 
# Reduce text size 
axis(side = 1, at = 1:length(groups), labels = groups, las = 2)
sample_sizes <- table(dat$Group)

boxplot_stats <- boxplot(VCFTOOLS ~ Group, data = dat, plot = FALSE) 

# Get stats without plotting again 
for (i in 1:length(groups)) { 
  group_name <- groups[i] 
  sample_size <- sample_sizes[group_name] # Position text above the boxplot whiskers 
  text(x = i, y = boxplot_stats$stats[5, i] + 0.0001, labels = paste0("n=", sample_size), col = "darkred", cex = 0.8) }
dev.off()
