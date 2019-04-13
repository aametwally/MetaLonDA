#' Fit longitudinal data
#'
#' Fits longitudinal samples from the same group using negative binomial smoothing splines or LOWESS
#' 
#' @param df dataframe has the Count, Group, ID, Time
#' @param method fitting method (nbinomial, lowess)
#' @param points points at which the prediction should happen
#' @return returns the fitted model
#' @import gss
#' @import stats
#' @references
#' Ahmed Metwally (ametwall@stanford.edu)
#' @examples 
#' data(metalonda_test_data)
#' n.sample = 5 
#' n.timepoints = 10 
#' n.group = 2 
#' Group = factor(c(rep(0, n.sample*n.timepoints), rep(1, n.sample*n.timepoints)))
#' Time = rep(rep(1:n.timepoints, times = n.sample), 2)
#' ID = factor(rep(1:(2*n.sample), each = n.timepoints))
#' points = seq(1, 10, length.out = 10)
#' aggretage.df = data.frame(Count = metalonda_test_data[1,], Time = Time, Group = Group, ID = ID)
#' cf = curveFitting(df = aggretage.df, method= "nbinomial", points)
#' @export
curveFitting = function(df, method = "nbinomial", points){
  
  ## Seprate the two groups
  group.null = df
  group.0 = df[df$Group==0, ]
  group.1 = df[df$Group==1, ]
  
  ## Fitting 
  if(method == "nbinomial"){
    
    ## null model
    mod.null = gssanova(Count ~ Time, data = group.null, family = "nbinomial", skip.iter=TRUE)
    mod.null.nbinomial.project = project(mod.null, c("Time"))
    
    ## full model
    mod.0 = gssanova(Count ~ Time, data = group.0, family = "nbinomial", skip.iter=TRUE)
    mod.1 = gssanova(Count ~ Time, data = group.1, family = "nbinomial", skip.iter=TRUE)
    mod.0.nbinomial.project = project(mod.0, c("Time"))
    mod.1.nbinomial.project = project(mod.1, c("Time"))
  }
  else if(method == "lowess"){
    
    ## null model
    mod.null = loess(Count ~ Time, data = group.null)
    
    ## full model
    mod.0 = loess(Count ~ Time, data = group.0)
    mod.1 = loess(Count ~ Time, data = group.1)
  }
  
  
  
  ## Calculate goodness of fit F-statistic for the non nbinomial models
  if(method !="nbinomial")
  {
    rss.null = summary(mod.null)$rss
    rss.full = summary(mod.0)$rss+summary(mod.1)$rss
    f.stat = (rss.null - rss.full)/rss.null
  }
  
  
  ## Estimate values at the provided time points
  est.null = predict(mod.null, data.frame(Time = points), se = TRUE)
  est.0 = predict(mod.0, data.frame(Time = points), se = TRUE)
  est.1 = predict(mod.1, data.frame(Time = points), se = TRUE)

  ## prepare dataframe for plotting
  if(method != "nbinomial")
  {
    ## Curve dataframe
    dd.null = data.frame(Time = points, Count = est.null$fit, Group = "NULL", ID = "NULL")
    dd.0 = data.frame(Time = points, Count = est.0$fit, Group = "fit.0", ID = "fit.0")
    dd.1 = data.frame(Time = points, Count = est.1$fit, Group = "fit.1", ID = "fit.1")
    
    ## Confidence interval dataframe
    dd.null.u95 = data.frame(Time = points, Count = (est.null$fit + 1.96*est.null$se), Group = "null.u", ID = "null.u")
    dd.null.l95 = data.frame(Time = points, Count = (est.null$fit - 1.96*est.null$se), Group = "null.l", ID = "null.l")
    dd.0.u95 = data.frame(Time = points, Count = (est.0$fit + 1.96*est.0$se), Group = "fit.0.u", ID = "fit.0.u")
    dd.0.l95 = data.frame(Time = points, Count = (est.0$fit - 1.96*est.0$se), Group = "fit.0.l", ID = "fit.0.l")
    dd.1.u95 = data.frame(Time = points, Count = (est.1$fit + 1.96*est.1$se), Group = "fit.1.u", ID = "fit.1.u")
    dd.1.l95 = data.frame(Time = points, Count = (est.1$fit - 1.96*est.1$se), Group = "fit.1.l", ID = "fit.1.l")
  } else{
    
    ## Curve dataframe
    dd.null = data.frame(Time = points, Count = mod.null$nu/exp(est.null$fit), Group = "null", ID = "null")
    dd.0 = data.frame(Time = points, Count = mod.0$nu/exp(est.0$fit), Group = "fit.0", ID = "fit.0")
    dd.1 = data.frame(Time = points, Count = mod.1$nu/exp(est.1$fit), Group = "fit.1", ID = "fit.1")
    
    ## Confidence interval dataframe
    dd.null.u95 = data.frame(Time = points, Count = mod.null$nu/exp(est.null$fit +1.96*est.null$se), Group = "null.u", ID = "null.u")
    dd.null.l95 = data.frame(Time = points, Count = mod.null$nu/exp(est.null$fit -1.96*est.null$se), Group = "null.l", ID = "null.l")
    dd.0.u95 = data.frame(Time = points, Count = mod.0$nu/exp(est.0$fit + 1.96*est.0$se), Group = "fit.0.u", ID = "fit.0.u")
    dd.0.l95 = data.frame(Time = points, Count = mod.0$nu/exp(est.0$fit - 1.96*est.0$se), Group = "fit.0.l", ID = "fit.0.l")
    dd.1.u95 = data.frame(Time = points, Count = mod.1$nu/exp(est.1$fit + 1.96*est.1$se), Group = "fit.1.u", ID = "fit.1.u")
    dd.1.l95 = data.frame(Time = points, Count = mod.1$nu/exp(est.1$fit - 1.96*est.1$se), Group = "fit.1.l", ID = "fit.1.l")
  }
  
  ## Return the results
  if(method != "nbinomial")
  {
    output = list(f.stat = f.stat, rss.null = rss.null, rss.full = rss.full, 
                  dd.null = dd.null, dd.0 = dd.0, dd.1 = dd.1, mod.null = mod.null, 
                  mod.0 = mod.0, mod.1 = mod.1, dd.null.u95 = dd.null.u95, dd.null.l95 = dd.null.l95,
                  dd.0.u95 = dd.0.u95, dd.0.l95 = dd.0.l95, dd.1.u95 = dd.1.u95, dd.1.l95= dd.1.l95)
  }
  else
  {
    output = list(dd.null = dd.null, dd.0 = dd.0, dd.1 = dd.1, mod.null = mod.null, mod.0 = mod.0, 
                  mod.1 = mod.1, mod.0.nbinomial.project = mod.0.nbinomial.project, 
                  mod.1.nbinomial.project = mod.1.nbinomial.project, 
                  mod.null.nbinomial.project = mod.null.nbinomial.project,
                  dd.null.u95 = dd.null.u95, dd.null.l95= dd.null.l95,
                  dd.0.u95 = dd.0.u95, dd.0.l95= dd.0.l95, dd.1.u95 = dd.1.u95, dd.1.l95= dd.1.l95)
  }
  return(output)
}


#' Calculate Area Ratio (AR) of each feature's time interval
#'
#' Calculate Area Ratio (AR) of each feature's time interval
#' 
#' @param curve.fit.df gss data object of the fitted spline
#' @return returns the area ratio for all time intervals
<<<<<<< HEAD
#' @import caTools
#' @references
#' Ahmed Metwally (ametwa2@uic.edu)
=======
#' @import pracma
#' @references
#' Ahmed Metwally (ametwall@stanford.edu)
>>>>>>> v1.1.5
#' @examples 
#' data(metalonda_test_data)
#' n.sample = 5
#' n.timepoints = 10
#' n.group= 2
#' Group = factor(c(rep(0,n.sample*n.timepoints), rep(1,n.sample*n.timepoints)))
#' Time = rep(rep(1:n.timepoints, times = n.sample), 2)
#' ID = factor(rep(1:(2*n.sample), each = n.timepoints))
#' points = seq(1, 10, length.out = 10)
#' aggregate.df = data.frame(Count = metalonda_test_data[1,], Time = Time, Group = Group, ID = ID)
#' model = curveFitting(df = aggregate.df, method= "nbinomial", points)
#' intervalArea(model)
#' @export
intervalArea = function(curve.fit.df){
  size = length(curve.fit.df$dd.null$Time)
  ar = numeric(size - 1)
  
  ## Calculate the absoulte and the sign of each interval area
  ar.abs = numeric(size - 1)
  ar.sign = numeric(size - 1)
  for(i in 1:(size - 1)){
    area.0 = trapz(curve.fit.df$dd.0$Time[i:(i+1)], curve.fit.df$dd.0$Count[i:(i+1)])
    area.1 = trapz(curve.fit.df$dd.1$Time[i:(i+1)], curve.fit.df$dd.1$Count[i:(i+1)])
    area.null = trapz(curve.fit.df$dd.null$Time[i:(i+1)], curve.fit.df$dd.null$Count[i:(i+1)])
    
    ar[i] = (abs(area.0) - abs(area.1)) / max(abs(area.0), abs(area.1))
    if(is.na(ar[i])){
     ar[i]=1
    }
    ar.abs[i] = abs(ar[i])
    ar.sign[i] = ar[i]/abs(ar[i])
  }
  
  return(list(ar = ar, ar.abs = ar.abs, ar.sign = ar.sign))
}



#' Find significant time intervals
#'
#' Identify significant time intervals
#' 
#' @param adjusted.pvalue vector of the adjusted p-value
#' @param threshold p-value cut off
#' @param sign vector hold area sign of each time interval 
#' @return returns a list of the start and end points of all significant time intervals
#' @references
#' Ahmed Metwally (ametwall@stanford.edu)
#' @examples 
#' p = c(0.04, 0.01, 0.02, 0.04, 0.06, 0.2, 0.06, 0.04)
#' sign = c(1, 1, 1, 1, -1, -1, 1, 1)
#' findSigInterval(p, threshold = 0.05, sign)
#' @export
findSigInterval = function(adjusted.pvalue, threshold = 0.05, sign)
{
  sig = which(adjusted.pvalue < threshold)
  sign = sign[sig]
  padj = adjusted.pvalue[sig]
  start = numeric()
  end = numeric()
  p = numeric()
  dom = numeric()
  
  if(length(sig) == 0)
  {
    cat("No Significant Intevals Found \n")
  }
  else if(length(sig) == 1)
  {
    start = sig[1]
    end = sig [1]
    p = padj[1]
    dom = sign[1]
  }
  else
  {
    start = sig[1]
    
    if((sig[2] - sig[1]) != 1 | sign[2] != sign[1])
    {
      end = c(end, sig[1])
      dom = c(dom, sign[1])
      p = c(p, padj[1])
    }
    
    for(i in 2:length(sig))
    {
      if(i != length(sig))
      {
        if((sig[i] - sig[i-1]) > 1 | sign[i] != sign[i-1])
        {
          start= c(start, sig[i])
        }
        
        if((sig[i+1] - sig[i]) != 1 | sign[i+1] != sign[i])
        {
          end = c(end, sig[i])
          dom = c(dom, sign[i])
          p = c(p, mean(adjusted.pvalue[start[length(start)] : end[length(end)]]))
        }
      }
      else
      {
        if((sig[i]-sig[i-1]) > 1 | sign[i] != sign[i-1])
        {
          start= c(start, sig[i])
        }
        end= c(end, sig[i])
        dom = c(dom, sign[i])
        p = c(p, mean(adjusted.pvalue[start[length(start)] : end[length(end)]]))
      }
    }
  }
  
  return(list(start = start, end = end, pvalue = p, dominant = dom))
}




#' Calculate Area Ratio (AR) of each feature's time interval for all permutations
#'
#' Fits longitudinal samples from the same group using negative binomial or LOWESS for all permutations
#' 
#' @param perm list has all the permutated models
#' @return returns a list of all permutation area ratio
#' @references
#' Ahmed Metwally (ametwall@stanford.edu)
#' @examples 
#' data(metalonda_test_data)
#' n.sample = 5 # sample size;
#' n.timepoints = 10 # time point;
#' n.perm = 3
#' n.group= 2 # number of group;
#' Group = factor(c(rep(0,n.sample*n.timepoints), rep(1,n.sample*n.timepoints)))
#' Time = rep(rep(1:n.timepoints, times = n.sample), 2)
#' ID = factor(rep(1:(2*n.sample), each = n.timepoints))
#' points = seq(1, 10, length.out = 10)
#' aggretage.df = data.frame(Count = metalonda_test_data[1,], Time = Time, Group = Group, ID = ID)
#' perm = permutation(aggretage.df, n.perm = 3, method = "nbinomial", points)
#' areaPermutation(perm)
#' @export
areaPermutation = function(perm)
{
  ar.list = list()
  list.len = length(perm)
  for (j in 1:list.len)
  {
    ar.list[[j]] = intervalArea(perm[[j]])
  }
  
  return(ar.list)
}
