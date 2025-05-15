# Volcano plot generation function
generate_volcano_plot <- function(df, filename) {
  df <- as.data.frame(df)
  x_1 <- df$log2FoldChange
  y_1 <- -log10(df$padj)
  
  # Set colors based on significance
  color <- ifelse(x_1 > 1 & y_1 > -log10(0.05), "red", 
                  ifelse(x_1 < -1 & y_1 > -log10(0.05), "blue", "gray"))
  
  # Generate volcano plot
  png(filename = paste0("volcano_plot/",filename, ".jpg"), width = 2400, height = 1800, res = 300)
  plot(x_1, y_1, pch = 20, col = color, 
       main = paste("volcano_plot/",filename, "Volcano Plot"), 
       xlab = "log2FC", ylab = "-log10 Padj")
  dev.off()
}
