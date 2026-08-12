MvAplot <- function(exprData, pdffilename, pcolor = "black", lcolor = "red") {
  # ensure data is numeric
  exprData <- as.data.frame(lapply(exprData, as.numeric))
  
  # extract column names
  column_names <- colnames(exprData)
  
  # first sample and its log-transformed values
  S_1 <- as.numeric(exprData[[column_names[1]]])
  log2_S1 <- log2(S_1 + 1)  # with offset to avoid log(0)
  
  pdf(pdffilename)
  
  # iterate over the columns (from the second one)
  for (colname in column_names[-1]) {
    col <- as.numeric(exprData[[colname]])
    log2_col <- log2(col + 1)  # with offset
    
    # M and A values
    M <- log2_S1 - log2_col
    A <- (log2_S1 + log2_col) / 2
    
    # plot for each sample vs. sample 1
    plot(A, M, 
         main = paste("MvA plot for", colname), 
         xlab = "A", 
         ylab = "M", 
         col = pcolor,        # points color
         pch = 16,            # solid circle
         cex = 0.6)           # point size 
    
    # add a horizontal line at y = 0
    abline(h = 0, col = lcolor)  
    
  }
  
  # close the PDF device
  dev.off()
}

# first function usage
# MvAplot(exprData = exprData, pdffilename = "plots_output.pdf", pcolor = "blue", lcolor = "green")

# exprData <- read.table('raw_count.txt', header=TRUE, sep='\t', row.names=1)
# annot <- read.table('gene_annot.txt', header=T, sep='\t', row.names='Symbol', quote="") 

TMMnorm <- function(exprData, annot, Mtrim=0.02, Atrim = c(0,8)) {

  # calculate sequencing depth (sum of each column)
  sequencing_depth <- colSums(exprData)
  
  # scale data: divide by sequencing depth and multiply by 10^6
  scaled_data <- t(t(exprData) / sequencing_depth) * 10^6
  
  # get column names
  column_names <- colnames(exprData)
  
  # extract and log-transform first sample
  S_1 <- as.numeric(scaled_data[, 1])
  log2_S1 <- log2(S_1 + 1)  # avoid log(0)
  
  # initialize lists for M and A values
  M_values <- list()
  A_values <- list()
  
  # iterate over columns (starting from the second)
  for (colname in column_names[-1]) {
    col <- as.numeric(scaled_data[, colname])
    log2_col <- log2(col + 1)  # avoid log(0)
    
    # calculate M and A values
    M <- log2_S1 - log2_col
    A <- (log2_S1 + log2_col) / 2
    
    # store M and A values
    M_values[[colname]] <- M
    A_values[[colname]] <- A
  }
  
  # SF with respect to S_1
  SF <- numeric(length(column_names) - 1) # init
  
  # SF for each sample
  for (i in 2:length(column_names)) {
    A <- A_values[[column_names[i]]]
    M <- M_values[[column_names[i]]]
    
    # trim M values corresponding to extreme A values
    valid_indices <- which(A >= Atrim[1] & A <= Atrim[2])
    valid_M <- M[valid_indices]
    
    # order M values and trim the most extreme values
    trimmed_M <- sort(valid_M)[(floor(length(valid_M) * Mtrim) + 1):(length(valid_M) - floor(length(valid_M) * Mtrim))]
    
    # trimmed mean of M
    SF[i - 1] <- mean(trimmed_M)
  }
  
  # names for SF
  names(SF) <- column_names[-1]
  
  # data normalization with computed scaling factors
  data_norm <- scaled_data[, -1]  # exclude the first column (subject 1)
  for (i in 1:nrow(data_norm)) {
    data_norm[i, ] <- data_norm[i, ] * (2^SF)  # row-wise multiplication
  }
  
  # append the first column of scaled data (subject 1) to the normalized data
  data_norm <- cbind(scaled_data[, 1], data_norm)
  
  # scale normalized data by gene length
  genes_lengths <- annot[, 3] 
  data_final <- data_norm
  for (i in 1:ncol(data_norm)) {
    data_final[, i] <- data_norm[, i] / genes_lengths
  }
  
  data_final <- data_final * 1e3
  
  # output data matrix
  normalized_matrix <- data_final  # normalized data in original scale
  
  # assign original row/column names to output data matrix
  colnames(normalized_matrix) <- colnames(exprData)
  rownames(normalized_matrix) <- rownames(exprData)
  
  # SF in original scale
  scaling_factors <- 2^SF
  scaling_factors <- c(1, scaling_factors)  # Append first factor equal to 1
  names(scaling_factors)[1] <- colnames(exprData)[1]  # add subject 1 name
  
  # list to store normalized data and scaling factors
  result <- list(
    normalized_matrix = normalized_matrix,
    scaling_factors = scaling_factors
  )
  return(result)
}

# second function usage
# exprData <- read.table('raw_count.txt', header=TRUE, sep='\t', row.names=1)
# annot <- read.table('gene_annot.txt', header=T, sep='\t', row.names='Symbol', quote="") 
# results <- TMMnorm(exprData, annot)

#DEbyEdgeR <- function(rawdata, groups, alpha = 0.05) {}






