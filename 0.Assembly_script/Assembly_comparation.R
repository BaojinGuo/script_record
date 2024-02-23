### Previous step: Align assembly genomes
## minimap2 -t 128 -x asm5 /scratch/pawsey0399/bguo1/0.assembly/01.hifi_assembly/S2_HIFI_RESULT/S2_hifi.asm.bp.p_ctg.fa /scratch/pawsey0399/bguo1/0.assembly/01.hifi_assembly/S1_HIFI_RESULT/S1_hifi.asm.bp.p_ctg.fa >S1_TO_S2.paf

# Load necessary libraries
library(pafr)      # For handling PAF (Pairwise mApping Format) files
library(patchwork) # For creating composite plots
library(tidyverse) # For data manipulation and visualization

# Read the PAF file into a data frame
df <- read_paf("S1_TO_S2.paf")

# Display the first few rows of the data frame
df %>% as.data.frame() %>% head()

# Create a dotplot visualization from the PAF data
dotplot(df, label_seqs = TRUE)

# Create a synteny plot for a specific chromosome pair
# (Adjust 'q_chrom' and 't_chrom' to match the desired chromosomes)
plot_synteny(df, q_chrom = "S1_1", t_chrom = "S2_1", centre = TRUE)

# Create a coverage plot, filling based on the query sequences
plot_coverage(df, fill = "qname")
