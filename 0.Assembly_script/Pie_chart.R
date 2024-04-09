library(ggplot2)

data <- read.csv("expan_constract.csv",header = T)

# Function to create pie chart for each row and save as PDF
create_pie_chart <- function(row) {
    # Create data frame for current row
    row_data <- data.frame(
        category = c("Increase", "Decrease"),
        value = c(row$Increase, row$Decrease)
    )
    
    # Create pie chart
    pie_chart <- ggplot(row_data, aes(x = "", y = value, fill = category)) +
        geom_bar(stat = "identity", width = 1) +
        coord_polar("y", start = 0) +
        labs(title = paste("Taxon_ID:", row$Taxon_ID), fill = "Category") +
        geom_text(aes(label = value), position = position_stack(vjust = 0.5)) +  # Add labels for values
        scale_fill_manual(values = c("Increase" = "red", "Decrease" = "green")) +  # Set fill colors
        theme_void()  # Remove axes and background
    
    # Save pie chart as PDF
    pdf_file <- paste("pie_chart_", gsub("[^[:alnum:]]", "_", row$Taxon_ID), ".pdf", sep = "")
    ggsave(pdf_file, plot = pie_chart, device = "pdf")
}

# Apply function to create pie chart for each row and save as PDF
lapply(1:nrow(data), function(i) create_pie_chart(data[i,]))
