MvAplot <- function(exprData, pdffilename, pcolor="black", lcolor="red") {
  
  # Start the PDF device to save the plots
  pdf(pdffilename)
  
  # Loop over each sample (column) in exprData, except for the first column (reference sample)
  for (i in 2:ncol(exprData)) {
    
    # Calculate M and A values
    M <- log2(exprData[,i] + 1) - log2(exprData[,1] + 1)  # log-ratio between sample i and sample 1
    A <- (log2(exprData[,i] + 1) + log2(exprData[,1] + 1)) / 2  # average log expression
    
    # Plot M vs. A
    plot(A, M,
         col = pcolor,
         pch = 20,          # Set point character
         cex = 0.5,         # Set point size
         xlab = "Average log expression (A)",    # x-axis label
         ylab = "Log ratio (M)",                 # y-axis label
         main = paste("Sample", i, "vs Sample 1")  # plot title
    )
    
    # Add a horizontal line at y = 0 with the specified line color
    abline(h = 0, col = lcolor, lwd = 1.5)
  }
  
  # Close the PDF device
  dev.off()
}

exprData <- as.matrix(read.table("raw_count.txt", header=TRUE, row.names=1, sep="\t"))
MvAplot(exprData, "MvA_plots.pdf")


