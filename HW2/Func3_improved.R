
#FN Calculation:
  
#Modified to ensure non-negative values by properly calculating TP and TN first

#Q-values Calculation:
  
#Added monotonicity enforcement to q-values 
#This ensures q-values are properly ordered with respect to p-values

DEbyEdgeR <- function(rawdata, groups, alpha = 0.05) {
  # Identify cases and controls by recognizing expression in character string
  ind_patients <- grep("Group",groups) # get samples related to patients
  ind_CTRL <- grep("CTRL",groups) # get samples related to control
  
  # Load library
  library(edgeR)
  
  # Create the DGE object
  y <- DGEList(counts = rawdata)   
  y <- calcNormFactors(y)   # This calculates the SF using the TMM normalization
  
  # Turn groups vector (a char) into a factor with 2 levels: cases, controls
  groups[ind_patients] <- 'Case' #unify patients char
  groups[ind_CTRL] <- 'CTRL' #unify controls char 
  # Factor
  Group <- factor(groups, levels = c('Case','CTRL')) #turn into factor
  
  # DESIGN matrix 
  design <- model.matrix(~ 0 + Group) # model with no intercept
  rownames(design) <- colnames(rawdata)  
  
  # Fit values of phi
  y <- estimateGLMCommonDisp(y,design, verbose=F) #phi common to the entire dataset
  y <- estimateGLMTrendedDisp(y,design) #phi depends on mu
  y <- estimateGLMTagwiseDisp(y,design) #phi is gene specific
  
  # GLM fit
  fit <- glmFit(y, design) 
  
  # Decide which factors/conditions I want to test
  comp <- makeContrasts(Subject = GroupCase - GroupCTRL, levels=design)
  # log ratio test (LRT)
  RES <- glmLRT(fit, contrast = comp[,"Subject"])
  
  #Expected outcomes counts
  G <- nrow(rawdata) #number of tests = number of genes
  G0 <- round(0.8*G) #expected number of negatives in the data set
  
  # Count selected genes for given alpha (no multiple tests correction)
  selected_genes <- sum(RES$table$PValue < alpha) #number of selected genes
  
  #Compute estimated FP = min(expected number of selected genes for given alpha, observed number of selected genes)
  FP <- min(round(G0*alpha), selected_genes) #expected false positives
  
  #Compute estimated FN = total number of tests - expected FP - expected TP - expected TN
  # Corrected FN calculation to ensure non-negative values
  TP <- selected_genes - FP  # Expected true positives
  TN <- min(G-selected_genes, G0-FP)  # Expected true negatives
  FN <- G - FP - TP - TN  # Expected false negatives
  
  #Compute FDR = expected FP / observed number of selected genes (if the selected number is zero, then FDR=0)
  FDR <- ifelse(selected_genes==0, 0, FP/selected_genes) #False Discovery Rate
  
  ## Benjamini-Hochberg correction for multiple tests
  #Compute ranks: tie breaking assigning maximum
  k <- rank(RES$table$PValue, ties.method = 'max')
  
  #Compute q-values according to B-H
  q_values <- RES$table$PValue*G/k
  # Ensure q-values are monotonic
  q_values <- pmin(1, cummin(q_values[order(RES$table$PValue)])[rank(RES$table$PValue)])
  
  ## Normalized expression data (compute with provided function)
  rawdata_norm <- TMMnorm(rawdata, 'gene_annot.txt', Mtrim = 0.02, Atrim = c(0,8))[[1]]
  
  #log Fold Change, use definition:
  #logFC = log2(Normalized_data_average_cases / Normalized_data_average_controls)
  logFC <- log2(apply(rawdata_norm[,ind_patients]+1, 1, mean) / 
                apply(rawdata_norm[,ind_CTRL]+1, 1, mean))
  
  # Create output vector of expected outcomes
  expected_rates <- c(SelectedGenes = selected_genes, 
                     FalsePositives = FP, 
                     FalseNegatives = FN, 
                     FalseDiscoveryRate = FDR)
  
  # Create output matrix of DE test results
  DE_results <- cbind(RES$table$PValue, q_values, logFC)
  # Assign row/column names
  rownames(DE_results) <- rownames(rawdata)
  colnames(DE_results) <- c('p-value','q-value','logFC')
  
  ##OUTPUT
  out <- list(expected_rates, DE_results)
  return(out)
}
