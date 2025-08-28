library(pegas)
library(ape)

# 找到所有 aln fasta 文件
files <- list.files(pattern = "line_.*\\.aln\\.fa$")

# 存放结果
hap_summary <- data.frame(File = character(),
                          HapCount = integer(),
                          stringsAsFactors = FALSE)

all_assign <- data.frame(File = character(),
                         Sequence = character(),
                         Haplotype = integer(),
                         stringsAsFactors = FALSE)

all_freq <- list()

for (f in files) {
    cat("Processing:", f, "\n")
    
    # 读入序列
    data <- read.dna(f, format = "fasta")
    
    # 计算单倍型
    h <- haplotype(data)
    
    # 单倍型频率
    hf <- haploFreq(data, what = 1, haplo = h)
    
    # 输出单倍型序列
    hap_fasta <- paste0(f, ".hap.fa")
    write.dna(h, file = hap_fasta, format = "fasta", nbcol = -1, colsep = "")
    
    # 单倍型数
    hap_count <- nrow(h)
    writeLines(paste("Number of haplotypes:", hap_count), paste0(f, ".hap.count.txt"))
    
    # 关键修改：展开 haplotype index
    hap_index <- attr(h, "index")
    assign_table <- do.call(rbind, lapply(seq_along(hap_index), function(i) {
        data.frame(
            File = f,
            Sequence = hap_index[[i]],  # 序列名
            Haplotype = paste0("Hap", i)  # 给 hap 起编号
        )
    }))
    write.table(assign_table, file = paste0(f, ".hap.assign.txt"),
                sep = "\t", quote = FALSE, row.names = FALSE)
    
    # 输出 freq 表
    hf_df <- as.data.frame(hf)
    hf_df$Haplotype <- paste0("Hap", seq_len(nrow(hf_df)))
    hf_df$File <- f
    write.table(hf_df, file = paste0(f, ".hap.freq.txt"),
                sep = "\t", quote = FALSE, row.names = FALSE)
    
    # 收集 summary 信息
    hap_summary <- rbind(hap_summary, data.frame(File = f, HapCount = hap_count))
    all_assign <- rbind(all_assign, assign_table)
    all_freq[[f]] <- hf_df
}
