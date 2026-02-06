# ============================================
# Heritability Analysis (Concise Version)
# ============================================

library(lme4)
library(lmerTest)
library(tidyr)
library(dplyr)

# 1. Read and prepare data
dat_wide <- read.csv(
  "trait.txt",
  header = TRUE,
  stringsAsFactors = FALSE,
  na.strings = c("NaN", "NA", "")
)

dat_long <- dat_wide %>%
  pivot_longer(
    cols = -`Genotype.code`,
    names_to = "Environment",
    values_to = "Trait"
  ) %>%
  filter(!is.na(Trait)) %>%
  rename(Genotype = `Genotype.code`) %>%
  mutate(
    Genotype = as.factor(Genotype),
    Environment = as.factor(Environment)
  )

# Data overview
cat("=== DATA OVERVIEW ===\n")
cat(sprintf("Observations: %d\n", nrow(dat_long)))
cat(sprintf("Genotypes: %d\n", n_distinct(dat_long$Genotype)))
cat(sprintf("Environments/Replicates: %d\n", n_distinct(dat_long$Environment)))
cat(sprintf("Missing rate: %.1f%%\n", 
            100 * (1 - nrow(dat_long) / (nrow(dat_wide) * (ncol(dat_wide) - 1)))))

# 2. Fit mixed model
m <- lmer(
  Trait ~ (1 | Genotype) + (1 | Environment),
  data = dat_long,
  REML = TRUE
)

# 3. Model summary
cat("\n=== MODEL SUMMARY ===\n")
print(summary(m))

# 4. Extract variance components
vc <- as.data.frame(VarCorr(m))
sigma_G <- vc %>% filter(grp == "Genotype") %>% pull(vcov)      # Genotypic variance
sigma_E <- vc %>% filter(grp == "Environment") %>% pull(vcov)   # Environmental variance
sigma_R <- vc %>% filter(grp == "Residual") %>% pull(vcov)      # Residual variance

# 5. Extract BLUPs and PEVs
cat("\n=== BLUP AND PEV EXTRACTION ===\n")
re <- ranef(m, condVar = TRUE)$Genotype
pev <- attr(re, "postVar")

# Handle different PEV matrix structures
if (length(dim(pev)) == 3) {
  pev_vector <- pev[1, 1, ]  # 3D array
} else if (is.matrix(pev)) {
  pev_vector <- diag(pev)    # 2D matrix
} else {
  pev_vector <- as.numeric(pev)  # Vector
}

cat(sprintf("Number of genotypes with PEV: %d\n", length(pev_vector)))
cat(sprintf("Mean PEV: %.6f\n", mean(pev_vector)))
cat(sprintf("SD of PEV: %.6f\n", sd(pev_vector)))
cat(sprintf("Range of PEV: [%.6f, %.6f]\n", min(pev_vector), max(pev_vector)))

# 6. Significance tests
cat("\n=== SIGNIFICANCE TESTS ===\n")

# Fixed effects (intercept)
anova_table <- anova(m)
print(anova_table)

# Random effects (likelihood ratio tests)
cat("\nRandom effects likelihood ratio tests:\n")

# Genotype effect
m_no_geno <- lmer(Trait ~ (1 | Environment), data = dat_long, REML = FALSE)
m_full <- lmer(Trait ~ (1 | Genotype) + (1 | Environment), data = dat_long, REML = FALSE)
lrt_geno <- anova(m_no_geno, m_full)
print(lrt_geno)
p_geno <- lrt_geno$`Pr(>Chisq)`[2]

# Environment effect
m_no_env <- lmer(Trait ~ (1 | Genotype), data = dat_long, REML = FALSE)
lrt_env <- anova(m_no_env, m_full)
print(lrt_env)
p_env <- lrt_env$`Pr(>Chisq)`[2]

# 7. Calculate heritability
# Cullis broad-sense heritability
H2_cullis <- 1 - mean(pev_vector) / (2 * sigma_G)

# Traditional broad-sense heritability
H2_trad <- sigma_G / (sigma_G + sigma_R)  # Single plant level

# 8. Create comprehensive results table
results <- data.frame(
  Parameter = c("Genotypic variance (σ²g)", 
                "Environmental variance (σ²e)", 
                "Residual variance (σ²r)",
                "Total variance",
                "Proportion of genotypic variance",
                "Proportion of environmental variance",
                "Proportion of residual variance",
                "Mean PEV",
                "SD of PEV",
                "Genotype effect p-value",
                "Environment effect p-value",
                "Cullis broad-sense heritability (H²c)",
                "Traditional broad-sense heritability (H²)"),
  Value = c(
    sprintf("%.6f", sigma_G),
    sprintf("%.6f", sigma_E),
    sprintf("%.6f", sigma_R),
    sprintf("%.6f", sigma_G + sigma_E + sigma_R),
    sprintf("%.1f%%", sigma_G/(sigma_G + sigma_E + sigma_R)*100),
    sprintf("%.1f%%", sigma_E/(sigma_G + sigma_E + sigma_R)*100),
    sprintf("%.1f%%", sigma_R/(sigma_G + sigma_E + sigma_R)*100),
    sprintf("%.6f", mean(pev_vector)),
    sprintf("%.6f", sd(pev_vector)),
    ifelse(p_geno < 0.001, "<0.001", sprintf("%.6f", p_geno)),
    ifelse(p_env < 0.001, "<0.001", sprintf("%.6f", p_env)),
    sprintf("%.6f", H2_cullis),
    sprintf("%.6f", H2_trad)
  ),
  Significance = c(
    "", "", "", "", "", "", "", "", "",
    ifelse(p_geno < 0.001, "***",
           ifelse(p_geno < 0.01, "**",
                  ifelse(p_geno < 0.05, "*", "ns"))),
    ifelse(p_env < 0.001, "***",
           ifelse(p_env < 0.01, "**",
                  ifelse(p_env < 0.05, "*", "ns"))),
    "", ""
  )
)

cat("\n=== FINAL RESULTS ===\n")
print(results, row.names = FALSE)

# 9. Save detailed BLUP and PEV table
blup_pev_df <- data.frame(
  Genotype = rownames(re),
  BLUP = re[, 1],
  PEV = pev_vector,
  Reliability = 1 - (pev_vector / sigma_G)
) %>% 
  arrange(desc(BLUP)) %>%
  mutate(
    Rank = 1:n(),
    Reliability = round(Reliability, 4)
  )

# 10. Save all results to files
write.csv(results, "Heritability_Analysis_Results.csv", row.names = FALSE)
write.csv(blup_pev_df, "Genotype_BLUP_PEV_Detailed.csv", row.names = FALSE)

# 11. Create summary statistics for PEV
pev_summary <- data.frame(
  Statistic = c("Mean", "SD", "Minimum", "25th Percentile", 
                "Median", "75th Percentile", "Maximum", "CV (%)"),
  Value = c(
    mean(pev_vector),
    sd(pev_vector),
    min(pev_vector),
    quantile(pev_vector, 0.25),
    median(pev_vector),
    quantile(pev_vector, 0.75),
    max(pev_vector),
    sd(pev_vector)/mean(pev_vector)*100
  )
)

write.csv(pev_summary, "PEV_Statistics.csv", row.names = FALSE)

# 12. Output recommendations and summary
cat("\n=== ANALYSIS COMPLETED ===\n")
cat(sprintf("Recommended heritability: Cullis H² = %.4f\n", H2_cullis))
cat(sprintf("Mean PEV: %.6f\n", mean(pev_vector)))
cat(sprintf("Average reliability: %.4f\n", mean(blup_pev_df$Reliability)))

cat("\n=== FILES GENERATED ===\n")
cat("1. Heritability_Analysis_Results.csv - All statistical results\n")
cat("2. Genotype_BLUP_PEV_Detailed.csv    - BLUP, PEV, and reliability for each genotype\n")
cat("3. PEV_Statistics.csv                - Summary statistics of PEV\n")
