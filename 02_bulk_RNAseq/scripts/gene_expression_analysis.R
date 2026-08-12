### Function 1: MvAplot()

## INPUT

# - exprData    : a numeric data matrix with the same format as the raw count data
#                 provided in stem (genes' IDs on rows and subjects' IDs on columns)
#
# - pdffilename : the name of the .pdf file where the plots will be saved
#
# - pcolor      : color used to draw the points (default: "black")
#
# - lcolor      : color used to draw the horizontal reference line (default: "red")

## OUTPUT

# - Generates a .pdf file containing MvA plots comparing each sample with sample 1,
#   with one plot per page. Each plot should:
#
#   • Have meaningful x-axis, y-axis, and title labels.
#   • Include a horizontal reference line at y = 0 using the color specified by lcolor.
#   • Use appropriate point size and type to ensure the plot is readable.

MvAplot <- function(exprData, pdffilename, pcolor = "black", lcolor = "red") {
    
    # Compute M values for each gene and each comparison with sample 1
    M <- log2(exprData[, 1] + 1) -
        apply(exprData[, -1] + 1, 2, log2)
    
    # Compute A values for each gene and each comparison with sample 1
    A <- (apply(exprData[, -1] + 1, 2, log2) +
              log2(exprData[, 1] + 1)) / 2
    
    # Sample names
    subjects_names <- colnames(exprData)
    
    # Open PDF file for plots
    pdf(pdffilename)
    
    for (i in 1:(ncol(exprData) - 1)) {
        
        graphics::plot(
            A[, i], M[, i],
            type = "p",
            col = pcolor,
            pch = 20,
            cex = 0.6,
            main = paste(
                "MvA plot of subject",
                subjects_names[1],
                "vs. subject",
                subjects_names[i + 1],
                sep = " "
            ),
            xlab = "Average: A = (log2(r*1) + log2(r*i))/2",
            ylab = "LogRatio: M = log2(r*1) - log2(r*i)",
            ann = TRUE,
            axes = TRUE,
            asp = 1
        )
        
        # Add horizontal reference line at M = 0
        abline(h = 0, col = lcolor)
    }
    
    dev.off()
}




### Function 2: TMMnorm()

## INPUT

# - exprData : numeric data matrix with the same format as the raw count data
#              provided in stem (genes' IDs on rows and subjects' IDs on columns)
#
# - annot    : data frame with the same format as the annotation data provided in stem
#
# - Mtrim    : fraction of observations to trim from each end of the sorted M
#              values before calculating their mean
#
# - Atrim    : vector of two elements specifying the lower and upper thresholds
#              used to trim extreme A values

## OUTPUT

# - out : list containing:
#
#   i) the normalized expression matrix on the original, non-logarithmic scale
#   ii) the vector of scaling factors on the original, non-logarithmic scale

TMMnorm <- function(exprData, annot, Mtrim = 0.02, Atrim = c(0, 8)) {

    # Extract gene lengths and sequencing depths
    genes_lengths <- annot[,3]
    sequencing_depth <- apply(exprData,2,sum)

    # Normalize counts by sequencing depth and express them as counts per million
    data.scaled <- apply(exprData, 1, function(row) {row.scaled <- row/sequencing_depth
                                                   return(row.scaled)})
    
    # Multiply data by 10^6 factor
    data.scaled <- t(data.scaled)*1e6  # transpose to maintain shape
    
    # Compute A values between sample 1 and each other sample
    A <- (apply(data.scaled[,-1] + 1, 2, log2) +
              log2(data.scaled[,1] + 1)) / 2
    
    # Keep genes whose A values fall within the specified trimming interval
    Atrim <- sort(Atrim)
    A_filter <- A > Atrim[1] & A < Atrim[2]
    # Strict inequalities exclude observations exactly at the trimming boundaries
    
    # Compute M values relative to sample 1; +1 prevents log2(0)
    M <- log2(data.scaled[,1] + 1) -
        apply(data.scaled[,-1] + 1, 2, log2)
    
    # Apply A filtering by excluding observations outside the specified range
    M_filter <- M
    M_filter[!A_filter] <- NaN
    
    # Compute trimmed mean of M values after excluding genes removed by A trimming
    SF <- apply(M_filter, 2, mean, trim = Mtrim, na.rm = TRUE)
    
    # Apply the TMM scaling factors to samples relative to sample 1
    data.norm <- apply(data.scaled[,-1], 1, function(row) {
        row * (2^SF)
    })
    data.norm <- t(data.norm)
    
    # Reattach the reference sample, whose scaling factor is 1
    data.norm <- cbind(data.scaled[,1], data.norm)
    
    # Normalize for gene length and express values per kilobase
    data.final <- apply(data.norm, 2, function(column) {
        column / genes_lengths
    })
    data.final <- data.final * 1e3
    
    # Restore original dimensions and names
    normalized_matrix <- data.final
    colnames(normalized_matrix) <- colnames(exprData)
    rownames(normalized_matrix) <- rownames(exprData)
    
    # Convert log2 scaling factors back to their original scale
    scaling_factors <- 2^SF
    scaling_factors <- c(1, scaling_factors)
    names(scaling_factors) <- colnames(exprData)
    
    out <- list(
        Normalized_exprData = normalized_matrix,
        ScalingFactors = scaling_factors
    )
    
    return(out)
} 




### Function 3: DEbyEdgeR()

## INPUT

# - rawdata : numeric data matrix with the same format as the raw count data
#             provided in stem (genes' IDs on rows and subjects' IDs on columns)
#
# - groups  : vector of group labels corresponding to the columns of rawdata
#
# - alpha   : significance level (default: 0.05)

## OUTPUT

# - out : list containing:
#
#   i) a vector with the number of selected genes, estimated false positives,
#      estimated false negatives, and estimated false discovery rate (FDR);
#
#   ii) a matrix containing the p-value, Benjamini-Hochberg adjusted q-value,
#       and log2 fold change for each gene, with gene names as row names.

DEbyEdgeR <- function(rawdata, groups, alpha = 0.05) {
    
    # Identify cases and controls by recognizing expression in character string
    ind_patients <- grep("Group", groups)
    ind_CTRL <- grep("CTRL", groups)
    
    # Load library
    library(edgeR)
    
    # Create the DGE object
    DGEObject <- DGEList(counts = rawdata)
    DGEObject <- calcNormFactors(DGEObject)
    
    # Turn groups vector into a factor with 2 levels: 'Case' for cases and 'CTRL' for controls
    groups[ind_patients] <- 'Case'
    groups[ind_CTRL] <- 'CTRL'
    
    Group <- factor(groups, levels = c('Case','CTRL'))
    
    # Design matrix
    design <- model.matrix(~ 0 + Group)
    rownames(design) <- colnames(rawdata)
    
    # Estimate dispersions
    DGEObject <- estimateGLMCommonDisp(DGEObject, design, verbose = FALSE)
    DGEObject <- estimateGLMTrendedDisp(DGEObject, design)
    DGEObject <- estimateGLMTagwiseDisp(DGEObject, design)
    
    # Fit negative binomial GLM
    fit <- glmFit(DGEObject, design)
    
    # Test Case vs CTRL
    comp <- makeContrasts(Subject = GroupCase - GroupCTRL, levels = design)
    RES <- glmLRT(fit, contrast = comp[, "Subject"])
    
    
    ## Expected outcomes counts
    
    G <- nrow(rawdata)
    G0 <- round(0.8 * G, digits = 0)
    
    # Number of genes selected at the specified significance level
    selected_genes <- sum(RES$table$PValue < alpha)
    
    # Estimated false positives
    FP <- min(round(G0 * alpha, digits = 0), selected_genes)
    
    # Estimated true positives and true negatives
    TP <- max(0, selected_genes - FP)
    TN <- min(G - selected_genes, G0 - FP)
    
    # Estimated false negatives
    FN <- G - FP - TP - TN
    
    # Estimated false discovery rate
    FDR <- ifelse(selected_genes == 0, 0, FP / selected_genes)
    
    
    ## Benjamini-Hochberg correction
    
    # Sort p-values and calculate their ranks
    p_order <- order(RES$table$PValue)
    p_sorted <- RES$table$PValue[p_order]
    ranks <- seq_along(p_sorted)
    
    # Calculate the raw BH-adjusted values
    q_sorted <- p_sorted * G / ranks
    
    # Enforce monotonicity of the adjusted p-values
    q_sorted <- rev(cummin(rev(q_sorted)))
    
    # q-values cannot exceed 1
    q_sorted <- pmin(q_sorted, 1)
    
    # Restore the original gene order
    q_values <- numeric(G)
    q_values[p_order] <- q_sorted
    
    
    # Log fold changes from the edgeR LRT output
    logFC <- RES$table$logFC
    
    
    # Create output vector of expected outcomes
    expected_rates <- c(
        SelectedGenes = selected_genes,
        FalsePositives = FP,
        FalseNegatives = FN,
        FalseDiscoveryRate = FDR
    )
    
    # Create output matrix of differential expression results
    DE_results <- cbind(
        RES$table$PValue,
        q_values,
        logFC
    )
    
    rownames(DE_results) <- rownames(rawdata)
    colnames(DE_results) <- c('p-value', 'q-value', 'logFC')
    
    
    ## OUTPUT
    
    out <- list(
        DEtestOutcomes = expected_rates,
        DEtestResults = DE_results
    )
    
    return(out)
}

