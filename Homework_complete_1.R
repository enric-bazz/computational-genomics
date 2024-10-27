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

compute_MAF<-function(SNPdata){
  
  #Get the number of SNPs and subjects from the  input matrix
  num_SNPs <-nrow(SNPdata)    #scans row of SNPdata matrix 
  num_subjects<-ncol(SNPdata) #scans columns of SNPdata matrix
  
  #Print number of SNPs and subjects
  cat('Number of SNPs:', num_SNPs,'\n')          #using 'cat' rather than 'print' because it concatenates the outputs
  cat('Number of subjects:', num_subjects,'\n')
  
  #Initialize an empty vector to store MAF values obtained for each SNP
  MAF_values <- numeric(num_SNPs) #'numeric' creates a numeric vector of length num_SNPs
  
  #For each SNP (so for each row)
  for (i in 1 : num_SNPs){
    
    #Get the SNP row (allele for all subjects)
    SNP<-SNPdata[i, ]
    
    #Count occurences of alleles 0 (AA) ,1 (Aa) ,2 (aa)
    all_0 <- sum(SNP==0,na.rm = TRUE)  #homozugous major alleles  #na.rm removes possible missing NaN values (as specified)
    all_1 <- sum(SNP==1,na.rm = TRUE)  #heterozygous alleles
    all_2 <- sum(SNP==2,na.rm = TRUE)  #homozygous minor alleles
    
    #Total number of alleles (two per subject)
    total_all < -2*sum(!is.na(SNP)) #is.na(SNP) returns TRUE for every element that's NA (missing), and FALSE 
                                    #otherwise; by negating it w\ ! we get that it returns TRUE for every non missing element
    
    #Calculate allele frequencies
    minor_all_freq <- ( 2 * all_2 + all_1 )/total_all # theory formula (2*aa + AA)/(total # alleles)
    
    #Minor allele frequency is the smallest between the major and minor frequency found
    MAF<- min(minor_all_freq, 1 - minor_all_freq) #respectively q and p =1-q 
    
    #Store results in teh empty matrix initialized before
    MAF_values[i]<-MAF
  }
  
  #Assign names to the MAF vector corresponding to the SNP IDs (=row names)
  names(MAF_values)<- rownames(SNPdata)
  
  #Return the name vector of MAF values
  return(MAF_values)
  # print(MAF_values) #un-comment to test function and see output 
}



## SECOND FUNCTION ##

#Takes as input a numeric data matrix supposed to have the same format of genetic data provided 
#Computes a HWE test for each SNP given as input by :
  #Calculating the chi^2_obs of the HWE test from data
  #Computing the p-value using the function pchisq()
#Prints on the terminal the names of SNPs having lowest p_value
#Saves in a pdf file a boxplot of the computed values (w\ xlabel,ylabel,title); the filename/path is specified in the
#input parameter pdf_file_name
#Returns the HWE test p_values for each SNP as a vector of numeric values w\ names corresponding to 
#the SPD IDs ( =chromosome_position) w\ same order they had in input matrix

HWE_test<- function (SNPdata, pdf_file_name){
  
  #Initialize an empty vector to store p_values for each SNP
  p_val<- numeric(nrow(SNPdata)) #empty numeric vector of length equal to # of rows in SNP data matrix
  
  #Loop over each SNP row to calculate the HWE test
  for (i in 1 : nrow(SNPdata)){
    
    #Extract data of current SNP
    snp_data<-table(SNPdata[i,])
    
    #Calculation of expected counts in order to compute the chi^2_obs statistic 
    #I implemented the same code used in the 1st function to obtain allele counts
    num_SNPs <-nrow(SNPdata) #scans rows of SNPdata matrix 
    
    #Count occurences of alleles 0,1,2
    all_0 <- sum(snp_data==0,na.rm=TRUE)  #homozugous major alleles
    all_1 <- sum(snp_data==1,na.rm=TRUE)  #heterozygous alleles
    all_2 <- sum(snp_data==2,na.rm=TRUE)  #homozygous minor alleles
    
    #Total number of alleles (two per subject)
    total_all<- 2 * sum(!is.na(snp_data))
    
    
    #Allele frequency parameters
    p <-(2 * all_0 + all_1)/total_all
    q <- 1 - p
    
    #Expected frequencies (not required, but implemented them anyway)
    exp_freq_all_0 <- p^2    #AA
    exp_freq_all_1 <- 2*p*q  #Aa
    exp_freq_all_2 <- q^2    #aa
    
    #Expected genotype counts obtained by multiplying each expected genotype frequency by the total numer of individuals
    expected_count_all_0 = (p^2)*sum(!is.na(SNP))   #AA
    expected_count_all_1 = (2*p*q)*sum(!is.na(SNP)) #Aa
    expected_count_all_2 = (q^2)*sum(!is.na(SNP))   #aa
    
    #Create vector to store observed and expected counts
    observed_counts<- c(all_0, all_1, all_2)
    expected_counts<- c(expected_count_all_0, expected_count_all_1, expected_count_all_2)
    
    #Calculate chi^2_obs statistics   
    chi_square_val <- sum((observed_counts - expected_counts))^2/(expected_counts )
    
    #Compute p_values
    p_val[i]<- pchisq(chi_square_val, df=1, lower.tail=FALSE) #df = 1: there's 1 degree of freedom 
                                                              #lower.tail = FALSE sets probabilites P[X>x^2_obs]
    
  }
  
  #Print SNPs with lowest p_values
  lowest_pval_SNP <- names(SNPdata)[order(p_val)[1:10]] #prints top 10 w\ lowest p-values # 'order' does the same as 'sort'
  cat('SNPs with lowest p-values:', lowest_pval_SNP,'\n')
  
  #Save boxplot to PDF
  pdf(pdf_file_name)
  boxplot(p_val, main ='HWE p_values for SNPs in data matrix ', xlab = 'SNPs', ylab = 'p-values')
  dev.off() #shuts down the specified device to provide better control 
  
  #Return vector of p_values
  return(p_val)
  #print(p_val) #un-comment to test function asn see output 
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


SNP_association_test <- function (filepath, indCTRL, MAFth=0.01, HWEalpha=0.01){
  
  #Line to read the filepath data to be sure it's read in the rigth way
  filepath<- read.table(filepath, header = TRUE,  sep = "\t", row.names = 1) 
  
  #Filter out SNPs with MAF < MAFth or (HWEp_value) < HWEalpha
  #Recall 1st and 2nd function to do so
  values_MAF <- compute_MAF(filepath[, indCTRL]) 
  pvalues_HWE <- HWE_test(filepath[ , indCTRL], 'test_boxplot.pdf') 
  
  #Filter out SNPs with MAF < MAFth or HWEp_value < HWEalpha 
  #Create vector that stores indexes of acceptable SNP IDs
  ind_filtered_snps <- names(values_MAF)[values_MAF >= MAFth & pvalues_HWE >= HWEalpha] 
  
  #Create dataset containing valid SNPs
  data <- data[ ,ind_filtered_snps, drop=FALSE] #columns are indexes of acceptable SNP IDs
                                                #'drop = FALSE' prevent matrix to drop to a vector 
  
  #Initialize empty vector to store p_values fo chi^2_obs test
  p_values <- numeric (length(ind_filtered_snps))
  names(p_values)<- ind_filtered_snps #assigns to every p_value the corresponding SNP ID
  
  #Get cases and controls data based on indexes in indCTRL
  control_data <- data[ ,indCTRL,drop=FALSE]  #columns specified by indCTRL
  case_data <- data[,-indCTRL,drop = FALSE]   #columns not specified by indCTRL
  
  #Loop over each SNP to calculate the chi^2_obs test for association
  for (snp in ind_filtered_snps){
    
    #Extract genotype data for cases and controls
    case_genotypes <- case_data[ ,snp]
    control_genotypes <- control_data[ ,snp]
    
    #Calculate observed frequencies
    obs_freq <- table(factor(case_genotypes,levels = c('AA_case','Aa_case','aa_case')), factor(control_genotypes, levels = c('AA_ctrl','Aa_ctrl','aa_ctrl')))
    
    #Compute row and column sums
    row_tot <-rowSums(obs_freq)
    col_tot <-colSums(obs_freq)
    grand_tot<-sum(obs_freq)
    
    #Calculate expected frequencies under the null hypothesis
    exp_freq <- outer(row_tot,col_tot)/grand_tot
    
    #Calculate chi^2_obs statistics
    chi_square_stat <- sum((obs_freq - exp_freq)^2/exp_freq , na.rm = TRUE)
    
    #Calculate p_value from chi^2_obs distribution w\ 1 degree of freedom
    p_values[snp]<- 1 -pchisq(chi_square_stat, df=1)
  }
  
  #Compute q-value for each SNP using the Benjamini - Hochberg procedure
  k <- length(p_values)
  p_values_sorted <- sort (p_values)           #sorts p:values and obtains indeced
  p_values_ranked <- p_values[p_values_sorted] #p values in ascending order
  q_val <- numeric (k)                         #initialize q values
  
  #Calculate q values w\ BH method
  for ( i in seq_along(p_values_ranked)){
    rank <- i
    q_val[p_values_sorted[i]] <- p_values_ranked[i]*k/rank
  }
  
  #We want non-decreasing q values
  for ( i in (k-1):1){
    q_val[p_values_sorted[i]] <- min(q_val[p_values_sorted[i]],q_val[p_values_sorted[i+1]])
  }
  
  #Return matrix with chromosome_position as row names and 8 columns w\ names c('AA_ctrl', 'Aa_ctrl', 'aa_ctrl', 'AA_case', 'Aa_case', 'aa_case','pval', 'qval'), for each row the matrix reports the # of occurences of each genotype for controls and cases respectively, the p-value and q-value
  results <- cbind (MAF = values_MAF[ind_filtered_snps], HWE_pvalues =pvalues_HWE[ind_filtered_snps], pval = p_values, qval = q_val)
  
  return(results)
}


