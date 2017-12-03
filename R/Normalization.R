#' Normalize count matrix 
#'
#' Normalize count matrix
#'
#' @param count count matrix
#' @param method normalization method
#' @references
#' Ahmed Metwally (ametwa2@uic.edu)
#' @export
normalize = function(count, method = "css"){
 # col.data=0 ## this line is for CRAN package
  if(method == "css")
  {
    cat("Normalization using CSS method \n")
    otu = metagenomeSeq::newMRexperiment(count)
    p.1 = metagenomeSeq::cumNormStatFast(otu, pFlag = TRUE)
    otu.2 = metagenomeSeq::cumNorm(otu, p = p.1)
    count.normalized = metagenomeSeq::MRcounts(otu.2, norm = TRUE)
  }
  else if(method == "tmm")
  {
    cat("Normalization using TMM method \n")
    factors = edgeR::calcNormFactors(count, method="TMM")
    eff.lib.size = colSums(count) * factors
    ref.lib.size = mean(eff.lib.size) #Use the mean of the effective library sizes as a reference library size
    count.normalized = sweep(count, MARGIN = 2, eff.lib.size, "/") * ref.lib.size 
  }
  else if(method == "ra")
  {
    cat("Normalization using Relative Abundance (RA) method \n")
    count.normalized  = apply(count, 2, function(x) (x/sum(x)))
  }
  else if(method == "log10")
  {
    cat("Normalization using log10 of the RA method \n")
    count.normalized  = apply(count, 2, function(x) log10(x/sum(x) + 1)) 
  }
  else if(method == "median_ratio")
  {
    cat("Normalization using Median-Ratio method \n")
    col.data = as.data.frame(cbind(colnames(count), rep(1:2, length.out= ncol(count)), rep(1:2, length.out= ncol(count))))
    rownames(col.data) = col.data[,1]
    col.data =  col.data[,-1]
    colnames(col.data) = c("test","condition")
    data.deseq = DESeq2::DESeqDataSetFromMatrix(countData = count, colData = col.data, ~ condition)
    cds = DESeq2::estimateSizeFactors( data.deseq )
    count.normalized = t( t(DESeq2::counts(cds)) / DESeq2::sizeFactors(cds) )
  }
  return(count.normalized)
}
