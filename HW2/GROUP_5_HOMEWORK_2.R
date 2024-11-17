### Function 1: MvAplot()
## INPUT
# - exprData    : a numeric data matrix that is supposed to have the same format of the raw count data
#                 provided in stem (i.e. genes’IDs on rows and subjects’IDs on columns)

# - pdffilename : the name of the .pdf file where to save the plots

# - pcolor      : the color to be used to draw the points

# - lcolor      : the color to be used to draw the line

## OUTPUT
# - Generates a .pdf file with MvA plots of each sample vs. sample 1 (one plot per page). For each plot:
#   • Provide meaningful names for the x-axis, y-axis and plot title.
#   • Draw a horizontal line corresponding to the x-axis (i.e. y = 0); the line must have as color the one
#     specified in lcolor.
#   • Tune the size of the plot and the size/type of the points to make the plot readable

# CODE
MvAplot <- function(exprData, pdffilename, pcolor, lcolor) {
  
  #Read input data as data frame
  exprData <- read.table(exprData, header=TRUE, sep='\t', row.names=1)
  
  # Log-Ratios of subject 1 against all other subjects
  # Compute M values of each gene for each pair of (subject 1, subject i) by:
  # - bringing all the data except the first column in log2 scale
  # - subtracting it column-wise to the first subject raw counts in log2 scale
  M <- log2(exprData[,1]+1)-apply(exprData[,-1]+1,2,log2)
  
  # Average expression of each subject and subject 1
  # Compute A values of each gene for each pair of (subject 1, subject i) by:
  # - bringing all the data except the first column in log2 scale
  # - adding column-wise the first subject raw counts in log2 scale
  # - divide the sum by a 2 factor --> average of raw counts in log2 scale
  A <- (apply(exprData[,-1]+1,2,log2)+log2(exprData[,1]+1))/2
  
  #NOTE: In both computations of A and M a value of 1 is added to all data in order to avoid a zero input
  # in the log2() function that would return a -Inf value
  
  # Save subjects names in a vector of length=ncol(exprData)
  subjects_names <- colnames(exprData)
  
  # Open pdf file for plots
  pdf(pdffilename)
  
  for(i in 1:(ncol(exprData)-1)){ #cycles through each subject (total number of columns = ncol(exprData)-1 )
   
    # Scatter plot of (A,M) values of each subject with graphics::plot() function
    graphics::plot(A[,i], M[,i], #defines points coordinates in MvA plot
                   type = 'p', #set plot type to points
                   col = pcolor, #set color of points
                   pch = 20, #set point-symbol: full dots
                   cex = 0.6, #set dot size scaling-factor: makes them smaller
                   main = paste('MvA plot of subject',subjects_names[1],'vs. subject',subjects_names[i+1],sep=' '), #set plot title
                   xlab = 'Average: A = (log2(r*1) + log2(r*i))/2', #set x-axis label
                   ylab = 'LogRatio: M = log2(r*1) - log2(r*i)', #set y-axes label
                   ann = T, axes = T, asp = 1) #set further graphic settings
    
    # Horizontal line for M=0 (coincides with x-axis)
    abline(h=0, col=lcolor)
  }
  dev.off() #close device for plotting
  
}  




### Function 2: TMMnorm()
## INPUT
# - exprData : a numeric data matrix that is supposed to have the same format of the raw count
#              data provided in stem (i.e. genes’IDs on rows and subjects’IDs on columns)

# - annot    : a data frame that is supposed to have the same format of the annotation data provided in stem

# - Mtrim    : a number between 0 and 1 indicating the fraction of observations to be deleted from each
#              end (positive and negative values) of the sorted vector M before calculating the mean (where Ms are the log ratios defined as in the MvA plot)

# - Atrim    : a vector of 2 elements indicating the lower and upper thresholds to trim the most
#              extreme values of A (where A is the average in log2 scale) 

## OUTPUT
# - out : list of 2 elements: 
#         i) the normalized matrix (in the original scale, not in log scale); 
#         ii) the vector of the scaling factors (in the original scale, not in log scale)

# CODE
TMMnorm <- function(exprData, annot, Mtrim, Atrim) {
  
  #Read input files as data frames
  exprData <- read.table(exprData, header=TRUE, sep='\t', row.names=1)
  annot <- read.table(annot, header=T, sep='\t', row.names='Symbol', quote="") #deactivate quote to continue file reading
  
  #Save vector of genes lengths from annot matrix: size=[#genes=nrow(exprData),1]
  genes_lengths <- annot[,3] 
  
  #Compute sequencing depths vector = total number of reads for each subject: size=[#subjects=ncol(exprData),1]
  sequencing_depth <- apply(exprData,2,sum) 
  
  #Scale data by sequencing depth: divide data by number of reads of each subject
  data.scaled <- apply(exprData, 1, function(row) {row.scaled <- row/sequencing_depth
                                                   return(row.scaled)})
 
  #Multiply data by 10^6 factor
  data.scaled <- t(data.scaled)*1e6 
  #NOTE: Whenever I use apply() on rows in a matrix object I need to transpose the output to maintain the shape
  
  ## Calculate scaling factor SF
  #Compute A values of scaled data
  A <- (apply(data.scaled[,-1]+1,2,log2)+log2(data.scaled[,1]+1))/2
  
  #Trim A values according to Atrim thresholds
  Atrim <- sort(Atrim, decreasing = F) #sort to guarantee correct order of values
  A_filter <- (A>Atrim[1] & A<Atrim[2])
  #returns a boolean matrix with 0 (FALSE) on elements trimmed (discarded) according to Atrim and
  # 1 (TRUE) on elements we want to keep
  
  #NOTE: Conditions on A trimming
  # A=0 if and only if both raw counts in subject 1 and the other subject are 0 (this gives (log(0+1)+log(0+1))/2=0 with the correction)
  # That also means M=0 for such subjects (log(0+1)-log(0+1)=0), thus excluding (strict inequality A>0) or including (A>=0) the A=Atrim[1]=0 value 
  # in the trimming procedure of A does affect the final values of the TMM (values of 0 are either included or not in the mean computation).
  
  
  #Compute M values
  M <- log2(data.scaled[,1]+1)-apply(data.scaled[,-1]+1,2,log2)
  
  #NOTE: in both computations of A and M a value of 1 is added to all data in order to avoid a zero input
  # in the log2() function that would return a -Inf value
  
  #Filter M values according to Atrim
  M_filter <- M #assing each element to the corresponding M value
  M_filter[!(A_filter)] <- NaN #set to NaN elements discarded with A trim (A_filter=0 and !(A_filter)=1)
  #returns a matrix with NaN on elements trimmed according to Atrim and original M values on other elements (the shape and size are maintained)

  #Compute scaling factor for each subject pair: size=[#subjects-1=ncol(exprData)-1,1]
  SF <- apply(M_filter, 2, mean, trim = Mtrim, na.rm = T)
  # Computes an algebraic mean of M values where the elements trimmed based on Atrim are set to NaN (do not get counted)
  # and elements kept based on Atrim have their original value
  # 'trim' argument defines the fraction (0 to 0.5) of observations to be trimmed from each end of M before the mean is computed
  # 'na.rm = T' ensures that NaN values are excluded from the computation
  
  #NOTE: On the original paper a weighted mean is used, with weights = 1/(i-th gene variance)
  #NOTE: In the instructions Mtrim can vary in the interval [0,1], the R mean() function only accepts values in [0,0.5]
  # values above 0.5 would not make sense since they would exceed the maximum acceptable value: removing a fraction > 0.5 from both ends means removing all the elements
  
  #Data normalization with computed scaling factors
  #Scale each column of scaled data (except for column 1) by multiplying each value for 2^SF of each subject
  data.norm <- apply(data.scaled[,-1], 1, function(row){row.normalized <- row*(2^SF) #row-wise multiplication
                                                        return(row.normalized)})
  data.norm <- t(data.norm)
  #NOTE: Whenever I use apply() on rows in a matrix object I need to transpose the output to maintain the shape
  
  #Append first column of scaled data (subject 1, not normalized) to the left of normalized data
  data.norm <- cbind(data.scaled[,1],data.norm) 
  
  #Scale normalized data 
  #Scale genes by gene length 
  data.final <- apply(data.norm, 2, function(column){column.scaled <- column/genes_lengths #column-wise multiplication
                                                     return(column.scaled)})
  #NOTE: with this procedure it is assumed that lengths in the 'gene_annot.txt' file are in the same order as
  # genes in the 'raw_count.txt' file (the gene in each row coincides in the two dataframes), true for the given files
  
  #Multiply data by 10^3
  data.final <- data.final*1e3 
  
  #Build output data matrix
  normalized_matrix <- data.final #normalized data in original scale
  
  #Assign original row/column names to output data matrix
  colnames(normalized_matrix) <- colnames(exprData)
  rownames(normalized_matrix) <- rownames(exprData)
  
  #Build vector of Scaling Factors in original scale
  scaling_factors <- 2^SF
  #Append first factor equal to 1 in original scale (2^0)
  scaling_factors <- c(1,scaling_factors)
  names(scaling_factors)[1] <- colnames(exprData)[1] #add subject 1 name
  
  out <- list(Normalized_exprData = normalized_matrix, ScalingFactors = scaling_factors) #output list
  return(out)
} 




### Function 3: DEbyEdgeR()
## INPUT
# - rawdata : a numeric data matrix that is supposed to have the same format of the raw count data
#             provided in stem (i.e. genes’IDs on rows and subjects’IDs on columns)

# - groups  : a vector of labels corresponding to the rawdata columns labels

# - alpha   : the significance level (default 0.05)

## OUTPUT
# - out : list of two elements
#      • A vector of 4 elements with: i) the number of selected genes in correspondence to the significance level alpha; 
#                                     ii) the corresponding estimate of the expected number of false positives (FP) considering G0=0.8*G (G is the total number of analyzed genes); 
#                                     iii) the corresponding estimate of the expected number of false negatives (FN); 
#                                     iv) the corresponding FDR.
#      • A matrix with p values, q values (in this case obtained with the Benjamini-Hockberg procedure) and log
#           fold change of each gene in the columns. The matrix must have row names corresponding to the gene
#           names.

# CODE
DEbyEdgeR <- function(rawdata, groups, alpha) {
  
  #NOTE: the input variable 'groups' is expected to be a character vector of length = ncol(rawdata)
  # with values corresponding to colnames(rawdata), each element identifying a specific subject/sample
  
  #Read input data as data frame
  rawdata <- read.table(rawdata, header=TRUE, sep='\t', row.names=1)
  
  # Identify cases and controls by recognizing expression in character string
  ind_patients <- grep("Group",groups) # get samples related to patients
  ind_CTRL <- grep("CTRL",groups) # get samples related to control
  
  # Load library
  library(edgeR)
  
  # Create the DGE object
  DGEObject <- DGEList(counts = rawdata)   
  DGEObject <- calcNormFactors(DGEObject)   #calculates SF using TMM normalization

  # Turn groups vector (a char) into a factor with 2 levels: 'Case' for cases, 'CTRL' for controls
  groups[ind_patients] <- 'Case' #unify patients char
  groups[ind_CTRL] <- 'CTRL' #unify controls char 
  # Factor
  Group <- factor(groups, levels = c('Case','CTRL')) #turn into factor

  # DESIGN matrix 
  design <- model.matrix(~ 0 + Group) # model with no intercept
  rownames(design) <- colnames(rawdata)  

  # Fit values of phi - estimate dispersions
  DGEObject <- estimateGLMCommonDisp(DGEObject,design, verbose=F) # Common dispersion
  DGEObject <- estimateGLMTrendedDisp(DGEObject,design) # Trended dispersion
  DGEObject <- estimateGLMTagwiseDisp(DGEObject,design) # Tagwise dispersion 
  
  # GLM fit
  fit <- glmFit(DGEObject, design) #model fit (that accounts for raw NB data and scaling factors and seq. depth) 
  
  # Decide which factors/conditions I want to test
  comp <- makeContrasts(Subject = GroupCase - GroupCTRL, levels=design)
  # log ratio test (LRT)
  RES <- glmLRT(fit, contrast = comp[,"Subject"])

  
  ## Expected outcomes counts
  G <- nrow(rawdata) #number of tests = number of genes
  G0 <- round(0.8*G, digits = 0) #expected number of negatives in the data set (as integer value)
  
  # Count selected genes for given alpha (no multiple tests correction)
  selected_genes <- sum(RES$table$PValue<alpha) #number of selected genes
  
  #Compute estimated FP = min(expected number of selected genes for given alpha, observed number of selected genes)
  FP <- min(round(G0*alpha, digits = 0), selected_genes) #expected false positives
  
  #Compute estimated FN = total number of tests - expected FP - expected TP - expected TN
  # Corrected FN calculation to ensure non-negative values
  TP <- max(0, selected_genes-FP)  # Expected true positives
  TN <- min(G-selected_genes, G0-FP)  # Expected true negatives
  FN <- G - FP - TP - TN  # Expected false negatives
  
  #Compute FDR = expected FP / observed number of selected genes (if the selected number is zero, then FDR=0)
  FDR <- ifelse(selected_genes==0,0,FP/selected_genes) #False Discovery Rate
  
  
  
  ## Benjamini-Hochberg correction for multiple tests
  # Compute ranks: tie breaking assigning maximum
  k <- rank(RES$table$PValue, ties.method = 'max')
  
  # Compute q-values according to B-H
  q_values <- RES$table$PValue*G/k
  
  #NOTE: p-values are not sorted since the rank() function does not require sorted input 
  
  #NOTE: A further step could be implemented to ensure q-values are monotonic with the following line
  #q_values <- pmin(1, cummin(q_values[order(RES$table$PValue)])[rank(RES$table$PValue)])
  
  # This would allow to select genes controlling for the FDR by simply selecting rows according to q_values<alpha
  
  
  # Compute log Fold Changes 
  # Retrieve logFC values directly from the LRT output
  logFC <- RES$table$logFC

  #NOTE: Alternatively   
  # - Compute normalized expression data with our TMMnorm() function as
  #rawdata_norm <- TMMnorm(rawdata, 'gene_annot.txt', Mtrim = 0.02, Atrim = c(0,8))$Normalized_exprData 
  
  # - Compute log Fold Change by using the definition: logFC = log2(Normalized_data_average_cases / Normalized_data_average_controls)
  #logFC <- log2(apply(rawdata_norm[,ind_patients]+1, 1, mean) / apply(rawdata_norm[,ind_CTRL]+1, 1, mean))
  
  # - A term equal to 1 is added to all data to prevent having a zero input in the log2() function
  # - The values are slightly different probably due to different normalization procedure 
  
  
  
  # Create output vector of expected outcomes
  expected_rates <- c(SelectedGenes = selected_genes, FalsePositives = FP, FalseNegatives = FN, FalseDiscoveryRate = FDR)
  
  # Create output matrix of Differential Expression test results
  DE_results <- cbind(RES$table$PValue,q_values,logFC)
  # Assign row/column names
  rownames(DE_results) <- rownames(rawdata)
  colnames(DE_results) <- c('p-value','q-value','logFC')
  
  ##OUTPUT
  out <- list(DEtestOutcomes = expected_rates, DEtestResults = DE_results)
  return(out)
}

