#' Permute group labels
#'
#' Permutes the group label of the samples in order to construct the AR empirical distibution
#'
#' @param perm.dat dataframe has the Count, Group, ID, Time
#' @param n.perm number of permutations
#' @param method The fitting method (negative binomial, LOWESS)
#' @param points The points at which the prediction should happen
#' @param lev the two level's name
#' @return returns the fitted model for all the permutations
#' @import plyr
#' @import utils
#' @references
#' Ahmed Metwally (ametwa2@uic.edu)
#' @examples 
#' data(metalonda_test_data)
#' n.sample = 5
#' n.timepoints = 10
#' n.perm = 3
#' n.group = 2
#' Group = factor(c(rep(0, n.sample*n.timepoints), rep(1, n.sample*n.timepoints)))
#' Time = rep(rep(1:n.timepoints, times = n.sample), 2)
#' ID = factor(rep(1:(2*n.sample), each = n.timepoints))
#' points = seq(1, 10, length.out = 10)
#' aggregate.df = data.frame(Count = metalonda_test_data[1,], Time = Time, Group = Group, ID = ID)
#' prm = permutation(aggregate.df, n.perm = 3, method = "nbinomial", points)
#' @export
permutation = function(perm.dat, n.perm = 500, method = "nbinomial", points, lev){

  n.subjects = length(unique(perm.dat$ID))
  cat("# of Subjects = ", n.subjects, "\n")
  
  ## Start permutation
  pp = list() 
  perm = 0 # to be able to store the value
  pp = llply(1:n.perm, function(j){
    cat("Permutation #", j,  "\n")
    for( i in levels(perm.dat$ID)){
      perm.uniq.len = 1
      m=0
      
      while(perm.uniq.len == 1){
        # cat("In While m = ", m, "\n")
        m=m+1
        perm.dat[which(perm.dat$ID == i),]$Group = rep(sample(c(0,1),1), sum(perm.dat$ID == i))# ,  replace = TRUE) #rep(sample(1:2), each = time.point)
        perm.uniq.len = length(unique(perm.dat$Group))
      }
    }
    
    
    g.0 = perm.dat[perm.dat$Group == 0, ]
    g.1 = perm.dat[perm.dat$Group == 1, ]
    g.min = max(sort(g.0$Time)[1], sort(g.1$Time)[1])
    g.max = min(sort(g.0$Time)[length(g.0$Time)], sort(g.1$Time)[length(g.1$Time)])
    
    if(g.min > min(points) | g.max < max(points))
    {
      cat("Special Case: generated permutation is out of range \n")
      assign(paste("Model", j, sep = "_"), NULL)
    } 
    else if (length(which(sum(g.0$Count) == 0 | sum(g.1$Count)==0)))
    {
      cat("Special Case: zero for all variable of one group \n")
      assign(paste("Model", j, sep = "_"), NULL)
    }
    else
    {
      perm = curveFitting(df = perm.dat, method = method, points)
      assign(paste("Model", j, sep = "_"), perm)
    }
  }, .parallel = FALSE, .progress = "none")
    
  
  pp[sapply(pp, is.null)] = NULL
  return(pp)
}  
