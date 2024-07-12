
library(patchwork)
library(ggplot2)


data<-read.csv("103sample.AUS_special_SV.PLOT.csv",header = T)
data$Chromosome<-as.factor(data$Chromosome)
data$POS<-as.numeric(data$POS)
data$TYPE<-as.factor(data$TYPE)
desired_order <-c("INS","DEL","DUP","INV")
data$TYPE<-factor(data$TYPE,levels=desired_order)
custom_colors <- c("DEL" = "#66c2a5", "INS" = "#fc8d62")
p<-ggplot(data, aes(x = POS, y = Chromosome, color = TYPE)) +
    geom_point(shape = '|', size = 5, alpha = 0.7) +
    labs(title = "", x = "", y = "Chromosome") +
    theme_minimal() +
    theme(
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        plot.title = element_text(size = 16, hjust = 0.5)
    ) +
    scale_color_manual(values = custom_colors)+ scale_x_continuous(breaks = c(0, 100000000, 200000000, 300000000, 400000000, 500000000, 600000000, 700000000, 800000000), labels = c("0Mbp", "", "200Mbp", "", "400Mbp", "", "600Mbp", "", "800Mbp"))+
    theme(text = element_text(family = "Times New Roman", face = "bold", size = 20))+
    theme(axis.text = element_text(size = 16,face = "bold"))+
    theme(legend.text = element_text(size = 16))+theme(legend.position = "bottom",  # 将图例放置在底部
                                                       legend.title = element_blank())




dat<-read.table("103sample.AUS_special_SV.density.bed",header = FALSE, sep = "\t", stringsAsFactors = FALSE)
colnames(dat) <- c('chrom', 'start', 'end', 'count')
dat$start_mb <- dat$start / 1e6
chrom_order <- c("7H", "6H", "5H", "4H", "3H", "2H", "1H")
dat$chrom <- factor(dat$chrom, levels = chrom_order)
density_plot <- ggplot(dat, aes(x = start_mb, y = count)) +
    geom_line() +
    labs(title = '',
         x = '',
         y = 'SV Count') +
    theme_minimal() +
    theme(legend.title = element_blank()) +
    facet_wrap(~chrom, ncol = 1)+
    theme(
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        plot.title = element_text(size = 16, hjust = 0.5))+scale_x_continuous(breaks = seq(0, 800, 200), labels = c("0Mbp", "200Mbp", "400Mbp", "600Mbp", "800Mbp"))+theme(text = element_text(family = "Times New Roman", face = "bold", size = 20)) +
    theme(axis.text = element_text(size = 16, face = "bold"))+geom_hline(yintercept = 5, linetype = "dashed", color = "red")
combined_plot <- (density_plot / p) + plot_layout(ncol = 1)
ggsave("103sample_AUS_special_SV_density.tiff",density_plot,width = 12, height = 10, dpi = 300, bg = "white")
ggsave("103sample_AUS_special_SV_distribution.tiff",p, width = 12, height = 2, dpi = 300, bg = "white")
ggsave("103sample_AUS_special_SV_distribution&density.tiff",combined_plot,width = 12, height = 16, dpi = 300, bg = "white")
