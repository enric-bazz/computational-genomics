#raw_count.txt is a tab delimited file that contains raw data 
#and patients are divided in 2 groups: 1) asymptomatic, 2) symptomatic

#gene_annot.txt provides additional info on different genes, including their length

## FIRST FUNCTION##

#input: #exprData = numeric data matrix w\ same format as raw_count (gene IDs = rows, subject IDs = columns)
        #pdffilename = name of .pdf file where to save plots
        #pcolor = color used to draw points
        #lcolor = color used to draw line
#output: .pdf file w\ MvA plots of each sample VS. sample 1 (one plot per page). For each plot:
         #provide meaningful names for x-axis and y-axis and title
         #draw horizontal line corresponding to x-axis (line color must be lcolor)
         #tune the size of the plot and size/type of the points to make plot readable

MvAplot <- function ( exprData, pdffilename, pcolor = 'black', lcolor = 'red'){
  
  # Read exprData as a data frame
  exprData <- read.table(exprData, header = TRUE, row.names = 1, sep = "\t") #specifies that the first column represents the row names
  
  # Open a PDFto save plots
  pdf(pdffilename) 
  
  samples <- colnames(exprData) #number of columns in exprData matrix
  
  #first sample and its log-transformed values
  S_1 <- as.numeric(exprData[,1])
  log2_S1 <- log2(S_1 + 1) 
  
  # Loop through each sample to create MA plots
  for (sample in samples[-1]) { #start from sample 2 since sample 1 is the reference sample 
    
    # Calculate M and A values
    M <- log2_S1 - log2(as.numeric(exprData[, sample]) + 1) # M = log2 mean #add 1 in both terms in order to exclude cases of log(0)
    A <- (log2(as.numeric(exprData[, sample]) + 1) + log2_S1 ) / 2  #  A = log2 average intensity #add 1 in both terms in order to exclude cases of log(0)
    
    # Print each plot on a new page in the PDF
    plot(A, M, pch = 16, cex = 0.6, col = pcolor, xlab = 'A = (log2(r*1) + log2(r*i))/2', ylab ='M = log2(r*1) - log2(r*i)', main = paste("MvA Plot - Sample", colnames(exprData)[sample], "vs Sample", colnames(exprData)[1]))
    
    # Horizontal line corresponding to x-axis
    abline(h = 0, col = lcolor, lwd = 1.5, lty = 2)
  }
  
  # Close the PDF device
  dev.off()
  
  
}

#test
#MvAplot('/Users/aurorativeron/Desktop/COMPUTER ENGINEERING/SECOND YEAR/FIRST SEMESTER/COMPUTATIONAL GENOMICS/homework 2/raw_count.txt','plot_MVA.pdf',pcolor='black',lcolor='red')

##SECOND FUNCTION##

#input:#exprData = numeric data matrix w\ same format as raw_count (gene IDs = rows, subject IDs = columns)
       #annot = data frame that provides additional info on different genes, including their length
       #Mtrim = number between 0 and 1 indicating fraction of observations to be deleted from each end (positive and negative values) of the sorted vector M before calculating the mean (Ms are the log ratios defined as in the MvA plot)
       #Atrim = vector of 2 elements indicating lower and upper thresholds to trim the most extreme values of A ( A is the average in log2 scale)
# what it does: #scales data by their sequencing depth and multiply by 10^6
                #calculates scaling factors SF (w\ respect to sample 1) by trimming the most extreme values of A and taking the trimmed means of M values ( use R function mean)
                #normalize data by their scaling factor w\ respect to sample 1
                #scales genes by their length and multiply by 10^3
#output: #list of 2 elements: i) normalized matrix (in original scale, not log) , ii) vector of scaling factors (in original scale, not log)

TMMnorm <- function (exprData, annot, Mtrim =0.2, Atrim = c(0,8)){
  
  # Read exprData as data frame
  exprData <- read.table(exprData, header = TRUE, row.names = 1, sep = "\t")
  
  # Read annot as data frame 
  annot <- read.table(annot, header = TRUE, row.names = 1, sep = "\t")
  
  gene_lengths <- annot[ ,ncol(annot)] #column of gene lengths in annot is always the last column
  
  # Scale data by sequencing depth and multiply by 10^6 ( = calculate CPM )
  scaled_data <- t(t(exprData)/colSums(exprData))*(10^6)
  
  samples <- colnames(scaled_data)
  
  # Initialize storage for M and A values
  M_values <- list()
  A_values <- list()
  
  # First sample and its log-transformed values
  S_1 <- as.numeric(scaled_data[,1])
  log2_S1 <- log2(S_1 + 1) 
  
  # Loop through each sample to calculate M and A values
  for (sample in samples[-1]) { #start from sample 2 since sample 1 is the reference sample 
    
    M <- (log2_S1 - log2(as.numeric(scaled_data[, sample]) +1)) # M = log2 mean #add 1 in both terms in order to exclude cases of log(0)
    A <- (log2(as.numeric(scaled_data[, sample]) +1) + log2_S1) / 2  #  A = log2 average intensity #add 1 in both terms in order to exclude cases of log(0)
    
    M_values[[sample]] <- M
    A_values[[sample]] <- A
  }
  
  # Initialize vector to store SFs
  SF <- numeric(length(samples)-1)
  
  # Loop through eac sample to calculate the scaling factors
  for(i in 2:length(samples)){
    A <-A_values[[samples[i]]]
    M <-M_values[[samples[i]]]
    
    # Trim A based on threshold
    valid_As <- which(A>=Atrim[1] & A<=Atrim[2])
    
    # Filter M indexes based on A values
    valid_Ms <- M[valid_As]
    
    # Trim M based on threshold
    trimmed_M <- sort(valid_Ms)[floor(length(valid_Ms)*Mtrim + 1):(length(valid_Ms)-floor(length(valid_Ms)*Mtrim))]
    
    # Trimmed mean of M and store in SF vector
    SF[i-1] <- mean(trimmed_M, na.rm=TRUE)
  }
  
  # SFs in original scale
  scaling_factors <- c(1,2^SF)
  names(scaling_factors) <- samples
  
  # Normalize data by SFs (excluding first sample)
  normalized_data <- scaled_data[,-1] #exclude subject 1
  
  for (i in 1:nrow(normalized_data)){
    normalized_data[i,] <- normalized_data[i,]*(2^SF)
  }
  
  # Append first column of scaled data (subject 1) to the normalized data
  normalized_data <- cbind(scaled_data[,1],normalized_data)
  
  # Scale by gene length and multiply by 10^3
  normalized_final_matrix <- apply(normalized_data, 2, function(x) (x / gene_lengths) * 10^3)
  
  # Assign original row and column names to the output matrix
  colnames(normalized_final_matrix) <- samples
  rownames(normalized_final_matrix) <- rownames(exprData)
  
  # Return normalized data and scaling factors (both in original scale)
  results <- list(normalized_matrix = normalized_final_matrix, scaling_factors = scaling_factors)
  
  return(results)
}
#test
TMMnorm('/Users/aurorativeron/Desktop/COMPUTER ENGINEERING/SECOND YEAR/FIRST SEMESTER/COMPUTATIONAL GENOMICS/homework 2/raw_count.txt','/Users/aurorativeron/Desktop/COMPUTER ENGINEERING/SECOND YEAR/FIRST SEMESTER/COMPUTATIONAL GENOMICS/homework 2/gene_annot.txt',0.2,c(2,8))

##THRID FUNCTION##

#input: #rawdata = numeric data matrix w\ same format as raw count
        #groups = vector of labels corresponding to rawdata column labels
        #alpha = significance level
#output: #vector of 4 elements w\ i) # of selected genes in correspondence to significance level alpha
         #                       ii) corresponding estimate of expected # of FP considering G0 =0.8G ( G = tot # of analyzed genes)
         #                      iii) corresponding estimate of expected # of FN
         #                       iv) corresponding FDR
         # a matrix w\ p_values, q_values (obtained w\ Benjamini Hockberg procedure) and log fold change of each gene in columns. row names mus correspond to the gene names

DEbyEdgeR <- function(rawdata, groups, alpha = 0.05) {
  
  # Read rawdata as a data frame
  if (is.character(rawdata)) {
    rawdata <- read.table(rawdata, header = TRUE, row.names = 1, sep = "\t")
  }
  
  # Load edgeR library
  library(edgeR)
  
  # (1) Prepare the data
  DGEObject <- DGEList(counts = rawdata) # Initialize object with counts
  DGEObject <- calcNormFactors(DGEObject) # TMM normalization
  
  # (2) Create design matrix
  sampleGroupFactor <- factor(groups, levels = c("Patient", "CTRL")) # Convert group labels to factors
  designMatrix <- model.matrix(~ 0 + sampleGroupFactor) # Create design matrix without intercept
  colnames(designMatrix) <- levels(sampleGroupFactor) # Assign group names to the columns of the design matrix
  
  # (3) Dispersion estimation
  DGEObject <- estimateGLMCommonDisp(DGEObject, designMatrix) # Common dispersion
  DGEObject <- estimateGLMTrendedDisp(DGEObject, designMatrix) # Trended dispersion
  DGEObject <- estimateGLMTagwiseDisp(DGEObject, designMatrix) # Tagwise dispersion
  
  # (4) Fit GLM model
  GLMFitResult <- glmFit(DGEObject, designMatrix)
  
  # (5) Perform hypothesis testing
  contrast <- makeContrasts(Patient - CTRL, levels = designMatrix) # Make contrast
  GLMTestResult <- glmLRT(GLMFitResult, contrast = contrast)
  
  # (6) Extract p_values, q_values, and log fold changes
  pValues <- as.vector(GLMTestResult$table$PValue) # Ensure pValues is a numeric vector
  
  # qValue computation
  m <- length(pValues) # Total number of tests
  sorted_pValues <- sort(pValues) # Ordered p_values
  rank <- rank(sorted_pValues, ties.method = "min") # Use rank to create rank vector
  qValues <- sorted_pValues * m / rank # Compute q_values with usual formula
  qValues <- qValues[order(order(pValues))] # Match q-values to the original p-values order
  
  logFoldChanges <- GLMTestResult$table$logFC # Extract log fold changes
  
  # (7) Summarize results
  G <- nrow(rawdata) # Total number of genes
  numSelected <- sum(pValues < alpha) # Number of selected genes
  G0 <- 0.8 * G # Estimated number of non-DE genes
  FPs <- min(numSelected, G0 * alpha) # Expected number of false positives
  FNs <- G - FPs - max(0, numSelected - FPs) - min(G - numSelected, G0 - FPs) # Expected number of false negatives
  fdr <- ifelse(numSelected > 0, FPs / numSelected, NA) # False discovery rate
  
  # (8) Create output
  resultStats <- c(
    NumSelected = numSelected,
    FalsePositives = FPs,
    FalseNegatives = FNs,
    FDR = fdr
  )
  resultMatrix <- cbind(Pvalue = pValues, Qvalue = qValues, LogFC = logFoldChanges)
  rownames(resultMatrix) <- rownames(rawdata)
  
  # Return results
  return(list(Statistics = resultStats, ResultsMatrix = resultMatrix))
}


rawdata <- read.table('/Users/aurorativeron/Desktop/COMPUTER ENGINEERING/SECOND YEAR/FIRST SEMESTER/COMPUTATIONAL GENOMICS/homework 2/raw_count.txt', header = TRUE, row.names = 1, sep = "\t")
indCTRL <- grep("CTRL",colnames(rawdata))
nCTRL <- length(indCTRL) #number of control samples
ncols <- ncol(rawdata) #total number of samples corresponds to the number of columns in rawdata
#assign groups: "CTRL" for control samples, "Patient" for the rest
groups <- c(rep("Patient", (ncols - nCTRL)), rep("CTRL", nCTRL))
DEbyEdgeR('/Users/aurorativeron/Desktop/COMPUTER ENGINEERING/SECOND YEAR/FIRST SEMESTER/COMPUTATIONAL GENOMICS/homework 2/raw_count.txt',groups,0.5)
