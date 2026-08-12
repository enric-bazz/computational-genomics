# SNP.data.txt has chromosome_position as row names and subject_ID as columns
#1200 patients (cases) in column 1:1200 and 800 control subjects (controls) in column 1201:2000
# 0 = homozygous genotype AA (0 mutated alleles)
# 1 = heterozygous genotype Aa (1 mutated allele)
# 2 = homozygous genotype aa (2 mutated alleles)
# Submit a single R script which implements 3 functions 
# ! use requested names of the functions and names of the input parameters
# ! R script must contain only the 3 functions (no other code)
# ! R script can have any name but must be an .R file
# ! comment code and use meaningful variable names in the body of the functions 


## FIRST FUNCTION ##

#Takes as input a numeric data matrix supposed to have the same format of the genetic data provided 
#Prints on the terminal the numbers of SNPs and subjects (computing them from input data matrix)
#Calculates the minor allele frequency (MAF) for each SNP
#Returns the MAF of each SNP as a vector of numeric values w\ names corresponding to the SNP IDs ( =chromosome_position) with
#same order they had in input matrix

compute_MAF <-function (SNPdata){
  #Get number of SNPs and subjects from the input matrix
  num_SNPs <- nrow(SNPdata)     #scans rows of SNPdata matrix to obtain number of SNPs
  num_subjects <- ncol(SNPdata) #scans columns of SNPdata matrix to obtain number of subjects
  
  #Print number of SNPs and subjects
  cat('Number of SNPs:', num_SNPs,'\n') #using 'cat' rather than 'print' because it concatenates the outputs
  cat('Number of subjects:', num_subjects,'\n')
  
  #Store MAF values for each obtained SNP
  MAF_values <- numeric(num_SNPs)
  
  #Assign names to the MAF vector corresponding to the SNP IDs
  MAF_names <- rownames(SNPdata)
  
  #Iterate through each SNP
  for (i in 1:num_SNPs){
    
    #Get the SNP row (which corresponds to the alleles for all subjects)
    SNP<- SNPdata[i,]
    
    #Count occurrences of alleles 0 (AA), 1 (Aa), 2(aa)
    allele_0<- sum(SNP==0, na.rm=TRUE)    #AA #na-rm removes possible missing Nan values
    allele_1 <- sum(SNP==1, na.rm =TRUE)  #Aa
    allele_2 <- sum(SNP==2, na.rm=TRUE)   #aa
    
    #Total number of alleles
    total_alleles <- allele_0 + allele_1 + allele_2
    
    #Calculate allele frequencies 
    p<- (2*allele_0 + allele_1)/(2*total_alleles) #formula is (2*AA + Aa)/(total number of alleles)
    q<- 1- p                                      #this is the MAF value
    
    #Store results in the empy matrix initialized before
    MAF_values[i] <- q
  }
  names(MAF_values)<-MAF_names
  
  #Return vector of MAF values
  return (MAF_values)
}


## SECOND FUNCTION ##

#Takes as input a numeric data matrix supposed to have the same format of genetic data provided 
#Computes a HWE test for each SNP given as input by :
#Calculating the chi^2_obs of the HWE test from data
#Computing the p_value using the function pchisq()
#Prints on the terminal the names of SNPs having lowest p_value
#Saves in a pdf file a boxplot of the computed values (w\ xlabel,ylabel,title); the filename/path is specified in the
#input parameter pdf_file_name
#Returns the HWE test p_values for each SNP as a vector of numeric values w\ names corresponding to 
#the SPD IDs ( =chromosome_position) w\ same order they had in input matrix

HWE_test <- function (SNPdata, pdf_file_name){
  num_SNPs <- nrow(SNPdata)           #number of SNPs
  p_val <-numeric(num_SNPs)           #empty numeric vector of length equal to the number of rows in SNPdata amtrix
  p_values_names <- rownames(SNPdata)
  
  #Loop over each SNP to calcualte the HWE test
  for (i in 1:num_SNPs){
    
    #Extract data of current SNP
    snp_data <-SNPdata[i,]
    
    #Calculation of expected counts and occurrences of alleles 0,1,2
    allele_0 <- sum(snp_data==0,na.rm=TRUE) #AA
    allele_1 <- sum(snp_data==1,na.rm=TRUE) #Aa
    allele_2 <- sum(snp_data==2, na.rm=TRUE) #aa
    
    #Total number of alleles
    total_alleles <- allele_0 + allele_1 + allele_2
    
    #Allele frequency parameters
    p<-(2*allele_0 + allele_1)/(2*total_alleles) #two alleles per subject
    q<- 1-p
    
    #Expected frequencies (not required, but implemented them anyway)
    exp_freq_all_0 <- p^2    #AA
    exp_freq_all_1 <- 2*p*q  #Aa
    exp_freq_all_2 <- q^2    #aa
    
    #Expected genotype counts obtained by multiplying each expected genotype frequency by the total number of alleles
    expected_count_allele_0 <-(p^2)*total_alleles   #AA 
    expected_count_allele_1 <-(2*p*q)*total_alleles #Aa
    expected_count_allele_2 <-(q^2)*total_alleles   #aa
    
    #Create observed values with chi-square statistics
    observed_values<-((allele_0-expected_count_allele_0)^2/expected_count_allele_0)+((allele_1-expected_count_allele_1)^2/expected_count_allele_1)+((allele_2 - expected_count_allele_2)^2/expected_count_allele_2)
    
    #Compute p_values
    p_value_statistics <- pchisq(observed_values, df=1, lower.tail=FALSE) #'df=1' because there's 1 degree of freedom, 'lower.tail=FALSE' sets probabilities P(X>chi-square)
    
    p_val[i] <- p_value_statistics
  }
  
  names(p_val)<-p_values_names
  sorted_pval <-sort(p_val)
  lowest_pval_SNPs <- names(head(sorted_pval,5))
  cat('SNPs with the lowest p-values:', lowest_pval_SNPs,'\n')
  
  #Save boxplot to PDF
  pdf(pdf_file_name)
  boxplot(p_val, main = 'Boxplot of HWE p_values for SNPs in data matrix', xlab = 'SNPs', ylab='p_values', las=2)
  dev.off() #shuts down the specified device to provide better control
  
  #Return vector of p_values
  return(p_val)
}

## THIRD FUNCTION ##

#Takes as input:
#A file name (entire path) in the filepath input parameter; this file is supposed to have the same format of the
#genetic data provided
#A vector of indexes indicating in which columns the input file has data from control subjects
#A threshold to filter SNPs w\ lower MAF (default 0.01) in input parameter MAFth
#A significance level alpha to be used to filter SNPs w\ lower p_values because possibly not in HWE (default 0.01) in the
#input parameter HWEalpha
#Reads the file and analyses the different SNPs:
#Filtering out SNPs w\ MAF<MAFth or HWE - p-value < HWEalpha (use previously developed functions to compute this)
#Calculating the chi^2_obs of the association test from the data
#Computing the p_value using the function pchisq()
#Computes the q_value for each SNP using the Benjamini-Hochberg procedure (don't use R built in fx)
#Returns a matrix w\ chromosome_position as row names and 8 columns w\ names c('AA_ctrl', 'Aa_ctrl', 'aa_ctrl', 'AA_case', 
#'Aa_case', 'aa_case','pval', 'qval'); for each row the matrix reports the number of occurences of each genotype for controls 
#and cases respectively, the p_value and q_value

SNP_association_test <-function(filepath,indCTRL,MAFth=0.01,HWEalpha=0.01){
  
  #Read data from the file
  SNPdata <- read.table(filepath, header = TRUE, sep='\t', row.names=1)   
  
  #Compute MAF and HWE p_values for each SNP using the previously computed functions
  MAF_values <- compute_MAF(SNPdata)
  HWE_p_values<-HWE_test(SNPdata[,indCTRL],'test_boxplot.pdf')
  
  #Filter out unwanted SNP values
  filtered_SNPs <- which(!(MAF_values < MAFth & HWE_p_values < HWEalpha))
  
  #Count genotypes AA, Aa, aa in control group
  genotype_count_ctrl <- apply(SNPdata[,indCTRL],1,function(row){
    AA <- sum (row == 0)
    Aa <- sum(row ==1)
    aa <- sum(row==2)
    return(c(AA,Aa,aa))
  })
  
  #Count genotypes AA, Aa, aa in case group
  genotype_count_cases <- apply(SNPdata[,-indCTRL],1,function(row){
    AA<-sum(row==0)
    Aa <- sum(row==1)
    aa<- sum(row==2)
    return(c(AA,Aa,aa))
  })
  
  #Transpose the obtained matrices for an easier indexing 
  genotype_count_ctrl <- t(genotype_count_ctrl)
  genotype_count_cases<- t(genotype_count_cases)
  
  colnames(genotype_count_ctrl) <- c('AA_ctrl', 'Aa_ctrl','aa_ctrl')
  colnames(genotype_count_cases) <-c('AA_case','Aa_case','aa_case')
  
  #Create a results data frame with genotype counts for filtered SNPs
  results <- data.frame(AA_ctrl = genotype_count_ctrl[filtered_SNPs,'AA_ctrl'], Aa_ctrl = genotype_count_ctrl[filtered_SNPs,'Aa_ctrl'],aa_ctrl = genotype_count_ctrl[filtered_SNPs,'aa_ctrl'],AA_case=genotype_count_cases[filtered_SNPs,'AA_case'],Aa_case=genotype_count_cases[filtered_SNPs,'Aa_case'],aa_case=genotype_count_cases[filtered_SNPs,'aa_case'],row.names=rownames(SNPdata)[filtered_SNPs])
  
  #Chi_square computation for each SNP
  chi_squared <- apply (results,1, function(row){
    num_ctrl <- sum(row[1:3])
    num_case <- sum(row[4:6])
    total <- num_ctrl + num_case
    
    #Xalculate allele frequencies in control and case groups 
    pval_ctrl <- (2*row['AA_ctrl'] + row['Aa_ctrl'])/(2*num_ctrl)
    pval_case <- (2*row['AA_case'] + row['Aa_case'])/(2*num_case)
    qval_ctrl <- 1 - pval_ctrl
    qval_case <- 1- pval_case
    
    #Calculate expected genotype counts 
    expected_ctrl <- c(num_ctrl*pval_ctrl^2,2*num_ctrl*pval_ctrl*qval_ctrl, num_ctrl*qval_ctrl^2)
    expected_case <- c(num_case*pval_case^2,2*num_case*pval_case*qval_case, num_case*qval_case^2)
    
    #Compute chi_square statistics
    chi_squared_ctrl <- sum((row[1:3] - expected_ctrl)^2/expected_ctrl)
    chi_squared_case <- sum((row[4:6]-expected_case)^2/expected_case)
    
    #Return the sum of chi_square statistics for combined control and case groups
    return(chi_squared_ctrl + chi_squared_case)
  })
  
  results$chi_squared_observed <-chi_squared #assigns contents of chi_squared to a new column named chi_squared_observed in dataframe results
  
  #p_value computation for chi_square test 
  p_values <- pchisq(results$chi_squared_observed, df=2,lower.tail=FALSE)
  results$pval <- p_values # assigns contents of p_values to new column named pval in dataframe results
  
  #q_value computation with Benjamini - Hochberg method
  k<- length(p_values)
  p_values_sorted<-sort(p_values)
  ranks<-rank(p_values_sorted)
  
  q_values<-p_values_sorted*k/ranks
  results$qval <-q_values[match(results$pval, p_values_sorted)]
  
  #Final result matrix
  final_results <- results[, c('AA_ctrl','Aa_ctrl','aa_ctrl','AA_case','Aa_case','aa_case','pval','qval')]
  
  #Return final results
  return(final_results)
  
  # print(final_results) #un-comment to test the function 
  
}

