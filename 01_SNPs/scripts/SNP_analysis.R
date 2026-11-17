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
#Returns the MAF of each SNP as a vector of numeric values with names corresponding to the SNP IDs
#(=chromosome_position) with the same order they had in the input matrix

compute_MAF <-function (SNPdata){
    #Coerce to a matrix: row-wise indexing (SNPdata[i,]) inside a loop is very
    #slow on a data.frame (each call rebuilds a one-row data.frame across all
    #columns) but fast on a matrix, so this keeps the function usable on large inputs
    SNPdata <- as.matrix(SNPdata)
    
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
        SNP <- SNPdata[i,]
        
        #Count occurrences of genotypes 0 (AA), 1 (Aa), 2 (aa)
        allele_0 <- sum(SNP==0, na.rm=TRUE)   #AA (na.rm=TRUE ignores any missing/NA values)
        allele_1 <- sum(SNP==1, na.rm=TRUE)   #Aa
        allele_2 <- sum(SNP==2, na.rm=TRUE)   #aa
        
        #Total number of subjects with a non-missing genotype for this SNP
        total_subjects <- allele_0 + allele_1 + allele_2
        
        #Calculate allele frequencies
        p <- (2*allele_0 + allele_1)/(2*total_subjects) #freq. of allele A: (2*AA + Aa)/(total number of alleles)
        q <- 1 - p                                      #freq. of allele a: this is the MAF value
        
        #Store result in the vector initialized before
        MAF_values[i] <- q
    }
    names(MAF_values) <- MAF_names
    
    #Return vector of MAF values
    return (MAF_values)
}


## SECOND FUNCTION ##

#Takes as input a numeric data matrix supposed to have the same format of the genetic data provided
#Computes a HWE test for each SNP given as input by:
#Calculating the chi^2_obs of the HWE test from the data
#Computing the p_value using the function pchisq()
#Prints on the terminal the name(s) of the SNP(s) having the lowest p_value
#Saves in a pdf file a boxplot of the computed values (with xlabel, ylabel, title); the filename/path
#is specified in the input parameter pdf_file_name
#Returns the HWE test p_values for each SNP as a vector of numeric values with names corresponding to
#the SNP IDs (=chromosome_position) with the same order they had in the input matrix

HWE_test <- function (SNPdata, pdf_file_name){
    #Coerce to a matrix for fast row-wise indexing in the loop below (see note
    #in compute_MAF())
    SNPdata <- as.matrix(SNPdata)
    
    num_SNPs <- nrow(SNPdata)           #number of SNPs
    p_val <- numeric(num_SNPs)          #empty numeric vector of length equal to the number of rows in SNPdata matrix
    p_values_names <- rownames(SNPdata)
    
    #Loop over each SNP to calculate the HWE test
    for (i in 1:num_SNPs){
        
        #Extract data of current SNP
        snp_data <- SNPdata[i,]
        
        #Count occurrences of genotypes 0 (AA), 1 (Aa), 2 (aa)
        allele_0 <- sum(snp_data==0, na.rm=TRUE) #AA
        allele_1 <- sum(snp_data==1, na.rm=TRUE) #Aa
        allele_2 <- sum(snp_data==2, na.rm=TRUE) #aa
        
        #Total number of subjects with a non-missing genotype for this SNP
        total_subjects <- allele_0 + allele_1 + allele_2
        
        #Allele frequency parameters
        p <- (2*allele_0 + allele_1)/(2*total_subjects) #freq. of allele A (two alleles per subject)
        q <- 1 - p                                      #freq. of allele a
        
        #Expected genotype counts under Hardy-Weinberg equilibrium, obtained by
        #multiplying each expected genotype frequency (p^2, 2pq, q^2) by the
        #total number of subjects
        expected_count_allele_0 <- (p^2)*total_subjects     #AA
        expected_count_allele_1 <- (2*p*q)*total_subjects   #Aa
        expected_count_allele_2 <- (q^2)*total_subjects     #aa
        
        #Chi-square statistic comparing observed vs expected genotype counts.
        #A genotype with expected count 0 (only possible for a monomorphic SNP,
        #i.e. p=0 or q=0) also has observed count 0, so it contributes 0 to the
        #statistic instead of the undefined 0/0
        chi_sq_term <- function(obs, exp) ifelse(exp == 0, 0, (obs-exp)^2/exp)
        observed_values <- chi_sq_term(allele_0, expected_count_allele_0) +
            chi_sq_term(allele_1, expected_count_allele_1) +
            chi_sq_term(allele_2, expected_count_allele_2)
        
        #Compute p_value. df=1 because one parameter (allele frequency) is
        #estimated from the data; lower.tail=FALSE gives P(X > chi-square)
        p_value_statistics <- pchisq(observed_values, df=1, lower.tail=FALSE)
        
        p_val[i] <- p_value_statistics
    }
    
    names(p_val) <- p_values_names
    
    #Identify the SNP(s) with the lowest p-value (handles ties correctly)
    min_pval <- min(p_val, na.rm = TRUE)
    lowest_pval_SNPs <- names(p_val)[!is.na(p_val) & p_val == min_pval]
    cat('SNP(s) with the lowest HWE p-value:', lowest_pval_SNPs,'\n')
    
    #Save boxplot to PDF
    pdf(pdf_file_name)
    boxplot(p_val, main = 'Boxplot of HWE p_values for SNPs in data matrix', xlab = 'SNPs', ylab = 'p_values', las = 2)
    dev.off() #shuts down the specified device to provide better control
    
    #Return vector of p_values
    return(p_val)
}

## THIRD FUNCTION ##

#Takes as input:
#A file name (entire path) in the filepath input parameter; this file is supposed to have the same
#format of the genetic data provided
#A vector of indexes indicating in which columns the input file has data from control subjects
#A threshold to filter SNPs with lower MAF (default 0.01) in input parameter MAFth
#A significance level alpha to be used to filter SNPs with lower p_values because possibly not in
#HWE (default 0.01) in the input parameter HWEalpha
#Reads the file and analyses the different SNPs:
#Filtering out SNPs with MAF < MAFth or HWE p-value < HWEalpha (using the previously developed
#functions to compute this)
#Calculating the chi^2_obs of the association test from the data
#Computing the p_value using the function pchisq()
#Computes the q_value for each SNP using the Benjamini-Hochberg procedure (not using an R built-in fx)
#Returns a matrix with chromosome_position as row names and 8 columns with names
#c('AA_ctrl', 'Aa_ctrl', 'aa_ctrl', 'AA_case', 'Aa_case', 'aa_case', 'pval', 'qval'); for each row the
#matrix reports the number of occurrences of each genotype for controls and cases respectively, the
#p_value and q_value

SNP_association_test <-function(filepath,indCTRL,MAFth=0.01,HWEalpha=0.01){
    
    #Read data from the file
    SNPdata <- read.table(filepath, header = TRUE, sep='\t', row.names=1)   
    
    #Compute MAF (on all subjects) and HWE p_values (on controls only) using the
    #previously defined functions
    MAF_values <- compute_MAF(SNPdata)
    #drop=FALSE keeps a matrix-like object (with SNP row names) even when a
    #single control column is selected
    HWE_p_values <- HWE_test(SNPdata[,indCTRL,drop=FALSE],'test_boxplot.pdf')
    
    #Filter out SNPs with (MAF < MAFth) OR (HWE p-value < HWEalpha):
    #keep only SNPs that pass BOTH thresholds
    filtered_SNPs <- which(MAF_values >= MAFth & HWE_p_values >= HWEalpha)
    
    #Count genotypes AA, Aa, aa in control group
    genotype_count_ctrl <- apply(SNPdata[,indCTRL,drop=FALSE],1,function(row){
        AA <- sum(row == 0, na.rm = TRUE)
        Aa <- sum(row == 1, na.rm = TRUE)
        aa <- sum(row == 2, na.rm = TRUE)
        return(c(AA,Aa,aa))
    })
    
    #Count genotypes AA, Aa, aa in case group
    genotype_count_cases <- apply(SNPdata[,-indCTRL,drop=FALSE],1,function(row){
        AA <- sum(row == 0, na.rm = TRUE)
        Aa <- sum(row == 1, na.rm = TRUE)
        aa <- sum(row == 2, na.rm = TRUE)
        return(c(AA,Aa,aa))
    })
    
    #Transpose the obtained matrices for easier indexing
    genotype_count_ctrl <- t(genotype_count_ctrl)
    genotype_count_cases <- t(genotype_count_cases)
    
    colnames(genotype_count_ctrl) <- c('AA_ctrl', 'Aa_ctrl','aa_ctrl')
    colnames(genotype_count_cases) <- c('AA_case','Aa_case','aa_case')
    
    #Create a results data frame with genotype counts for the filtered SNPs only
    results <- data.frame(
        AA_ctrl = genotype_count_ctrl[filtered_SNPs,'AA_ctrl'],
        Aa_ctrl = genotype_count_ctrl[filtered_SNPs,'Aa_ctrl'],
        aa_ctrl = genotype_count_ctrl[filtered_SNPs,'aa_ctrl'],
        AA_case = genotype_count_cases[filtered_SNPs,'AA_case'],
        Aa_case = genotype_count_cases[filtered_SNPs,'Aa_case'],
        aa_case = genotype_count_cases[filtered_SNPs,'aa_case'],
        row.names = rownames(SNPdata)[filtered_SNPs]
    )
    
    #Chi-square computation for each SNP: test of independence between genotype
    #and case/control status (3 genotypes x 2 groups contingency table)
    chi_squared <- apply(results,1, function(row){
        num_ctrl <- sum(row[1:3])
        num_case <- sum(row[4:6])
        total <- num_ctrl + num_case
        
        #Row totals (cases + controls combined) for each genotype category
        genotype_totals <- c(row['AA_ctrl'] + row['AA_case'],
                             row['Aa_ctrl'] + row['Aa_case'],
                             row['aa_ctrl'] + row['aa_case'])
        
        #Expected counts under the null hypothesis of no association between
        #genotype and case/control status
        expected_ctrl <- genotype_totals * num_ctrl/total
        expected_case <- genotype_totals * num_case/total
        
        #Compute chi-square statistic comparing observed vs expected counts.
        #A genotype category with expected count 0 necessarily also has observed
        #count 0 (nobody carries that genotype in either group), so it
        #contributes 0 to the statistic instead of the undefined 0/0
        chi_sq_term <- function(obs, exp) ifelse(exp == 0, 0, (obs-exp)^2/exp)
        chi_squared_ctrl <- sum(chi_sq_term(row[1:3], expected_ctrl))
        chi_squared_case <- sum(chi_sq_term(row[4:6], expected_case))
        
        #Return the combined chi-square statistic for the SNP
        return(chi_squared_ctrl + chi_squared_case)
    })
    
    results$chi_squared_observed <- chi_squared #store the observed chi-square statistic for each SNP
    
    #p_value computation for the chi-square association test. df=2 because the
    #contingency table has (3 genotypes - 1) x (2 groups - 1) = 2 degrees of freedom;
    #lower.tail=FALSE gives P(X > chi-square)
    p_values <- pchisq(results$chi_squared_observed, df=2, lower.tail=FALSE)
    results$pval <- p_values #store the p-value for each SNP
    
    #q_value computation with Benjamini-Hochberg method (step-down procedure)
    k <- length(p_values)
    p_order <- order(p_values)          #indices that sort p-values in ascending order
    p_sorted <- p_values[p_order]
    i <- seq_len(k)                     #rank of each sorted p-value
    raw_q <- p_sorted*k/i
    
    #Enforce monotonicity: going from the largest to the smallest p-value, each
    #q-value cannot be smaller than the q-value of the next larger p-value
    q_sorted <- rev(cummin(rev(raw_q)))
    q_sorted <- pmin(q_sorted, 1)        #q-values cannot exceed 1
    
    #Restore the original SNP order
    q_values <- numeric(k)
    q_values[p_order] <- q_sorted
    results$qval <- q_values
    
    #Final result matrix, with only the columns requested. The instructions ask
    #for a matrix: as.matrix() keeps the SNP IDs as row names and gives a numeric
    #matrix since all columns are numeric
    final_results <- as.matrix(results[, c('AA_ctrl','Aa_ctrl','aa_ctrl','AA_case','Aa_case','aa_case','pval','qval')])
    
    #Return final results
    return(final_results)
}

