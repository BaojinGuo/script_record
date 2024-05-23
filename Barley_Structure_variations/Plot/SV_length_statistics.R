library(ggplot2)
library(dplyr)





sv_data <- read.table("4.combineSV_Clipper_sniffle_het_to_miss80_maf25.SV100k.length", header = FALSE, col.names = c("Chromosome", "Position", "SVType", "SVLength"))
sv_data$SVLength <- as.numeric(sv_data$SVLength)
sv_stats <- sv_data %>%
    group_by(SVType) %>%
    summarize(
        Count = n(),
        MeanLength = mean(SVLength, na.rm = TRUE),
        MedianLength = median(SVLength, na.rm = TRUE),
        MinLength = min(SVLength, na.rm = TRUE),
        MaxLength = max(SVLength, na.rm = TRUE),
        TotalLength = sum(SVLength, na.rm = TRUE)
    )
write.csv(sv_stats,"sv_stats.csv",row.names = F)
p1<-ggplot(sv_data %>% filter(SVType %in% c("INS", "DEL")), aes(x = SVLength, fill = SVType)) +
    geom_histogram(position = "dodge", bins = 50) +
    scale_fill_manual(values = c("INS" = "#E41A1C", "DEL" = "#377EB8")) +
    labs(title = "Length Distribution of INS and DEL",
         x = "Length",
         y = "Count") +
    theme_minimal() +
    theme(
        text = element_text(size = 14),
        plot.title = element_text(size = 16, hjust = 0.5)
    )
ggsave("ins_del_length_distribution.tiff", plot = p1, width = 8, height = 6, units = "in", bg = "white")
p2 <- ggplot(sv_data %>% filter(SVType == "DUP"), aes(x = SVLength, fill = SVType)) +
    geom_histogram(fill = "#4DAF4A", bins = 50) +
    labs(title = "Length Distribution of DUP",
         x = "Length",
         y = "Count") +
    theme_minimal() +
    theme(
        text = element_text(size = 14),
        plot.title = element_text(size = 16, hjust = 0.5)
    )
ggsave("dup_length_distribution.tiff", plot = p2, width = 8, height = 6, units = "in", bg = "white")
p3 <- ggplot(sv_data %>% filter(SVType == "INV"), aes(x = SVLength, fill = SVType)) +
    geom_histogram(fill = "#984EA3", bins = 50) +
    labs(title = "Length Distribution of INV",
         x = "Length",
         y = "Count") +
    theme_minimal() +
    theme(
        text = element_text(size = 14),
        plot.title = element_text(size = 16, hjust = 0.5)
    )
ggsave("inv_length_distribution.tiff", plot = p3, width = 8, height = 6, units = "in", bg = "white")
sv_data <- sv_data %>%
    mutate(SVLength = ifelse(SVType == "DEL", abs(SVLength), SVLength))
sv_data$SVType <- factor(sv_data$SVType, levels = c("DEL", "INS", "DUP", "INV"))
y_limits <- list(
    DEL = c(0, max(sv_data$SVLength[sv_data$SVType == "DEL"], na.rm = TRUE)),
    INS = c(0, max(sv_data$SVLength[sv_data$SVType == "INS"], na.rm = TRUE)),
    DUP = c(0, max(sv_data$SVLength[sv_data$SVType == "DUP"], na.rm = TRUE)),
    INV = c(0, max(sv_data$SVLength[sv_data$SVType == "INV"], na.rm = TRUE))
)
ggplot(sv_data, aes(x = SVType, y = SVLength, fill = SVType)) +
    geom_boxplot() +
    labs(title = "Distribution of Structural Variant Lengths",
         x = "",
         y = "Length") +
    theme_minimal() +
    facet_wrap(~ SVType, scales = "free_x", nrow = 1) +  # Adjusted scales to "free_x"
    theme(
        axis.text.x = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 14, face = "bold"),
        plot.title = element_text(hjust = 0.5)
    ) +
    scale_fill_manual(values = c("DEL" = "#E41A1C", "INS" = "#377EB8", "DUP" = "#4DAF4A", "INV" = "#984EA3"))+theme(
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        plot.title = element_text(size = 16, hjust = 0.5),
        text = element_text(family = "Times New Roman", face = "bold", size = 20),
        legend.text = element_text(size = 16),
        panel.background = element_rect(fill = "white"))+theme(legend.position="none")
