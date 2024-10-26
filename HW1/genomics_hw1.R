# Author: Kristina Svetova 
# load SNP data 
# file_path <- 'SNPdata.txt'
# SNPdata <- read.table(file_path, header = TRUE, sep = "\t", row.names = 1)

# function to compute MAF
compute_MAF <- function(SNPdata) {
  n_snps <- nrow(SNPdata) # number of SNPs
  n_subjects <- ncol(SNPdata) # number of subjects
  mafs <- numeric(n_snps)
  mafs_names <- rownames(SNPdata)
  
  # iterate through each SNP
  for (i in 1:n_snps) {
    snp_row <- SNPdata[i,]
    count_0 <- sum(snp_row == 0) # AA
    count_1 <- sum(snp_row == 1) # Aa
    count_2 <- sum(snp_row == 2) # aa
    
    N <- count_0 + count_1 + count_2
    
    p <- (2 * count_0 + count_1) / (2 * N)
    q <- 1 - p
    
    mafs[i] <- q
  }
  
  names(mafs) <- mafs_names
  print(paste("Number of SNPs:", n_snps))
  print(paste("Number of subjects:", n_subjects))
  
  return(mafs)
}

# function to perform HWE test for controls
HWE_test <- function(SNPdata, pdf_file_name) {
  n_snps <- nrow(SNPdata)
  p_values <- numeric(n_snps)
  p_values_names <- rownames(SNPdata)
  
  # iterate through each SNP
  for (i in 1:n_snps) {
    snp_row <- SNPdata[i,]
    count_0 <- sum(snp_row == 0) # AA
    count_1 <- sum(snp_row == 1) # Aa
    count_2 <- sum(snp_row == 2) # aa
    
    N <- count_0 + count_1 + count_2
    
    p <- (2 * count_0 + count_1) / (2 * N)
    q <- 1 - p
    
    expected_0 <- N * p^2
    expected_1 <- 2 * N * p * q
    expected_2 <- N * q^2
    
    obs_value <- ((count_0 - expected_0)^2 / expected_0) + 
      ((count_1 - expected_1)^2 / expected_1) + 
      ((count_2 - expected_2)^2 / expected_2)
    
    p_value <- pchisq(obs_value, df = 1, lower.tail = FALSE)
    p_values[i] <- p_value
  }
  
  names(p_values) <- p_values_names
  p_values_sorted <- sort(p_values)
  lowest_p_values_names <- names(head(p_values_sorted, 5))
  print(lowest_p_values_names)
  
  pdf(pdf_file_name)
  boxplot(p_values, main = "Boxplot of computed p-values", 
          xlab = "SNPs", 
          ylab = "p-values", 
          las = 2)
  dev.off()
  
  return(p_values)
}

# function to perform SNP association test
SNP_association_test <- function(filepath, indCTRL, MAFth=0.01, HWEalpha=0.01) {
  SNPdata <- read.table(filepath, header = TRUE, sep = "\t", row.names = 1)
  
  # compute MAF and HWE p-values
  maf_values <- compute_MAF(SNPdata)
  hwe_p_values <- HWE_test(SNPdata[, indCTRL], "test_boxplot.pdf")
  
  # filter based on MAF and HWE p-values
  filtered_snps <- which(!(maf_values < MAFth | hwe_p_values < HWEalpha))
  
  genotype_counts_ctrl <- apply(SNPdata[, indCTRL], 1, function(row){
    AA <- sum(row == 0)
    Aa <- sum(row == 1)
    aa <- sum(row == 2)
    return(c(AA, Aa, aa))
  })
  
  genotype_counts_case <- apply(SNPdata[, -indCTRL], 1, function(row){
    AA <- sum(row == 0)
    Aa <- sum(row == 1)
    aa <- sum(row == 2)
    return(c(AA, Aa, aa))
  })
  
  genotype_counts_ctrl <- t(genotype_counts_ctrl)
  genotype_counts_case <- t(genotype_counts_case)
  
  colnames(genotype_counts_ctrl) <- c("AA_ctrl", "Aa_ctrl", "aa_ctrl")
  colnames(genotype_counts_case) <- c("AA_case", "Aa_case", "aa_case")
  
  results <- data.frame(
    AA_ctrl = genotype_counts_ctrl[filtered_snps, "AA_ctrl"],
    Aa_ctrl = genotype_counts_ctrl[filtered_snps, "Aa_ctrl"],
    aa_ctrl = genotype_counts_ctrl[filtered_snps, "aa_ctrl"],
    AA_case = genotype_counts_case[filtered_snps, "AA_case"],
    Aa_case = genotype_counts_case[filtered_snps, "Aa_case"],
    aa_case = genotype_counts_case[filtered_snps, "aa_case"],
    row.names = rownames(SNPdata)[filtered_snps]
  )
  
  # chi-squared computation
  chi_squared <- apply(results, 1, function(row){
    N_ctrl <- sum(row[1:3])
    N_case <- sum(row[4:6])
    N <- N_ctrl + N_case
    
    p_ctrl <- (2 * row['AA_ctrl'] + row['Aa_ctrl']) / (2 * N_ctrl)
    p_case <- (2 * row['AA_case'] + row['Aa_case']) / (2 * N_case)
    q_ctrl <- 1 - p_ctrl
    q_case <- 1 - p_case
    
    expected_ctrl <- c(N_ctrl * p_ctrl^2, 2 * N_ctrl * p_ctrl * q_ctrl, N_ctrl * q_ctrl^2)
    expected_case <- c(N_case * p_case^2, 2 * N_case * p_case * q_case, N_case * q_case^2)
    
    chi_squared_ctrl <- sum((row[1:3] - expected_ctrl)^2 / expected_ctrl)
    chi_squared_case <- sum((row[4:6] - expected_case)^2 / expected_case)
    
    return(chi_squared_ctrl + chi_squared_case)
  })
  
  results$chi_squared_obs <- chi_squared
  
  # p-values computation
  p_values <- pchisq(results$chi_squared_obs, df = 1, lower.tail = FALSE)
  results$pval <- p_values
  
  # q-values computation using Benjamini-Hochberg procedure
  m <- length(p_values)
  p_values_sorted <- sort(p_values)
  ranks <- rank(p_values_sorted)
  
  q_values <- p_values_sorted * m / ranks
  results$qval <- q_values[match(results$pval, p_values_sorted)]
  
  # final results matrix
  final_results <- results[, c("AA_ctrl", "Aa_ctrl", "aa_ctrl", "AA_case", "Aa_case", "aa_case", "pval", "qval")]
  
  return(final_results)
}

# test
association_results <- SNP_association_test('SNPdata.txt', 1201:2000)
association_results
