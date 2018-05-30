#' Metagenomic Longitudinal Differential Abundance Analysis for one feature
#'
#' Find significant time intervals of the one feature
#' 
#' @param Count matrix has the number of reads that mapped to each feature in each sample.
#' @param Time vector of the time label of each sample.
#' @param Group vector of the group label of each sample.
#' @param ID vector of the subject ID label of each sample.
#' @param n.perm number of permutations.
#' @param fit.method fitting method (nbinomial, lowess).
#' @param points points at which the prediction should happen.
#' @param text Feature's name.
#' @param parall boolean to indicate whether to use multicore.
#' @param pvalue.threshold p-value threshold cutoff for identifing significant time intervals.
#' @param adjust.method multiple testing correction method.
#' @param time.unit time unit used in the Time vector (hours, days, weeks, months, etc.)
#' @return returns a list of the significant time intervals for the tested feature.
#' @import parallel
#' @import doParallel
#' @import stats
#' @references
#' Ahmed Metwally (ametwa2@uic.edu)
#' @examples 
#' data(metalonda_test_data)
#' n.sample = 5
#' n.timepoints = 10
#' n.group = 2
#' Group = factor(c(rep(0, n.sample*n.timepoints), rep(1,n.sample*n.timepoints)))
#' Time = rep(rep(1:n.timepoints, times = n.sample), 2)
#' ID = factor(rep(1:(2*n.sample), each = n.timepoints))
#' points = seq(1, 10, length.out = 10)
#' \dontrun{
#' output.nbinomial = metalonda(Count = metalonda_test_data[1,], Time = Time, Group = Group,
#' ID = ID, fit.method =  "nbinomial", n.perm = 10, points = points,
#' text = rownames(metalonda_test_data)[1], parall = FALSE, pvalue.threshold = 0.05,
#' adjust.method = "BH")
#' }
#' @export
metalonda = function(Count, Time, Group, ID, n.perm = 500, fit.method = "nbinomial", 
                      points, text = 0, parall = FALSE, pvalue.threshold = 0.05, 
                      adjust.method = "BH", time.unit = "days")
{
  cat("Start MetaLonDA \n")

  # Extract groups
  Group = as.character(Group)
  group.levels = sort(unique(Group))

  
  
  if(length(group.levels) > 2){
    stop("You have more than two phenotypes.")
  } else if(length(group.levels) < 2){
    stop("You have less than two phenotypes.")
  }

  
  gr.1 = as.character(group.levels[1])
  gr.2 = as.character(group.levels[2])

  Group[which(Group == gr.1)] = 0
  Group[which(Group == gr.2)] = 1

  
  
  ## Preprocessing: add pseudo counts for zero abudnance features
  Count = Count + 1e-8

  
  ## Form MetaLonDA dataframe
  aggregate.df = data.frame(Count = Count, Time = Time, Group = Group, ID = ID)


  ## Visualize feature's abundance accross different time points  
  visualizeFeature(aggregate.df, text, group.levels, unit = time.unit)


  
  group.0 = aggregate.df[aggregate.df$Group == 0, ]
  group.1 = aggregate.df[aggregate.df$Group == 1, ]
  points.min = max(sort(group.0$Time)[1], sort(group.1$Time)[1])
  points.max = min(sort(group.0$Time)[length(group.0$Time)], sort(group.1$Time)[length(group.1$Time)])
  points = points[which(points >= points.min & points <= points.max)]
  


  
  cat("Start Curve Fitting \n") 
  #if (fit.method == "ss")
  #{
  #  cat("Fitting: Gaussian SS \n") 
  #  model= curveFitting(df = aggregate.df, method= "ss", points)
  #}  else if (fit.method == "lowess")
  if (fit.method == "lowess")
  {
    cat("Fitting: LOWESS \n")
    model = curveFitting(df = aggregate.df, method= "lowess", points)
  } else if (fit.method == "nbinomial")
  {
    cat("Fitting: NB SS \n")
    model = tryCatch({
      curveFitting(df = aggregate.df, method= "nbinomial", points)
      },  error = function(err) {
        print(paste("ERROR in gss = ", err, sep="")); 
        return("ERROR")
        })
  } else {
    cat("You have entered unsupported fitting method\n")
    quit()
    
  }
  
  ## Visualize feature's trajectories spline
  visualizeFeatureSpline(aggregate.df, model, fit.method, text, group.levels, unit = time.unit)
 
   
  ## Calculate area under the fitted curve for each time interval
  cat("Calculate Area Under the Fitted Curves \n")
  area = intervalArea(model)
  
  
  ## Run in Parallel
  if(parall == TRUE) {
    max.cores = detectCores()
    desired.cores = max.cores - 1		
    cl = makeCluster(desired.cores)
    registerDoParallel(cl)
  } 

  
  ## Permutation 
  cat("Start Permutation \n")
  perm  = permutation(aggregate.df, n.perm, fit.method, points, lev = group.levels)
  
  ## Area p-value per unit interval
  area.perm = areaPermutation(perm)


  if(parall == TRUE) {
    stopCluster(cl)
  } 
  
  a1 = do.call(rbind, area.perm)
  a2 = do.call(rbind, a1[,2])

  ##  Visualize AR empirical distribution
  #visualizeARHistogram(a2, text, fit.method)
  
  ## Calculate AR p-value 
  pvalue.area = sapply(1:(length(points)-1), function(i){
    sum(a2[,i] >= area$ar.abs[i])/length(a2[,i])
  } )



  ## Identify significant time inetrval based on the adjusted p-value 
  cat("p-value Adjustment Method = ", adjust.method, "\n")
  adjusted.pvalue = p.adjust(pvalue.area, method = adjust.method)
  interval = findSigInterval(adjusted.pvalue, threshold = pvalue.threshold, sign = area$ar.sign)
  st = points[interval$start]
  en = points[interval$end + 1]
  
  
  if(length(st) > 0)
  {
    ## Visualize sigificant area
    visualizeArea(aggregate.df, model, fit.method, st, en, text, group.levels, unit = time.unit)
  }
  
  cat("\n\n")
  
  output.details = list(feature = text, significant.interval = cbind(start = st, end = en), 
       intervals.pvalue = pvalue.area, adjusted.pvalue = adjusted.pvalue, area.sign = area$ar.sign)
  output.summary = data.frame(feature = rep(text, length(interval$start)), start = st, end = en,
                 dominant = interval$dominant, pvalue = interval$pvalue)
  
  return(list(detailed = output.details, summary = output.summary))
}


#' Metagenomic Longitudinal Differential Abundance Analysis for all Features
#'
#' Identify significant features and their significant time interval
#' 
#' @param Count Count matrix of all features
#' @param Time Time label of all samples
#' @param Group Group label of all samples
#' @param ID individual ID label for samples
#' @param n.perm number of permutations
#' @param fit.method The fitting method (nbinomial, lowess)
#' @param num.intervals The number of time intervals at which metalonda test differential abundance 
#' @param parall logic to indicate whether to use multicore
#' @param pvalue.threshold p-value threshold cutoff
#' @param adjust.method Multiple testing correction methods
#' @param time.unit time unit used in the Time vector (hours, days, weeks, months, etc.)
#' @param norm.method normalization method to be used to normalize count matrix (css, tmm, ra, log10, median_ratio) 
#' @param prefix prefix for the output figure
#' @return Returns a list of the significant features a long with their significant time intervals
#' @references
#' Ahmed Metwally (ametwa2@uic.edu)
#' @examples 
#' \dontrun{
#' data(metalonda_test_data)
#' n.sample = 5
#' n.timepoints = 10
#' n.group = 2
#' Group = factor(c(rep(0, n.sample*n.timepoints), rep(1,n.sample*n.timepoints)))
#' Time = rep(rep(1:n.timepoints, times = n.sample), 2)
#' ID = factor(rep(1:(2*n.sample), each = n.timepoints))
#' points = seq(1, 10, length.out = 10)
#' output.nbinomial = metalondaAll(Count = metalonda_test_data, Time = Time, Group = Group,
#' ID = ID, n.perm = 10, fit.method =  "nbinomial", num.intervals = 100, 
#' parall = FALSE, pvalue.threshold = 0.05, adjust.method = "BH", 
#' time.unit = "hours", norm.method = "none", prefix = "Test")
#' }
#' @export
metalondaAll = function(Count, Time, Group, ID, n.perm = 500,
                         fit.method = "nbinomial", num.intervals = 100, parall = FALSE, 
                         pvalue.threshold = 0.05, adjust.method = "BH", time.unit = "days", 
                         norm.method = "none", prefix = "Output")
{
  ## Check the dimentions of the annotation vectors and count matrix
  if(length(Time) == length(ID))
  {
    if(ncol(Count) == length(Group))
    {
      if(length(Time) == length(Group))
      {
        cat("Dimensionality check passed\n")
      } else
      {
        stop("The length of the annotation vectors don't match or not consistent 
             with the number of columns of the count matrix.")
      }
    } else
    {
      stop("The length of the annotation vectors don't match or not consistent 
           with the number of columns of the count matrix.")
    }
  } else
  {
    stop("The length of the annotation vectors don't match or not consistent 
         with the number of columns of the count matrix.")
  }
  
  
  ## Normalization
  if(norm.method != "none")
  {
    if(norm.method == "css" | norm.method == "tmm" | norm.method == "ra" | norm.method == "log10" | norm.method == "median_ratio")
    {
      cat("Normalizaton method = ", norm.method, "\n")
      Count = normalize(Count, method = norm.method)
    } else{
      stop("You have entered a wrong normalization method")
    }
  }
  
  ## Specify the test/prediction timepoints for metalonda
  if(num.intervals == "none")
    #points = floor(seq(min(Time), max(Time)))
    points = seq(min(Time), max(Time))
  else
    #points = floor(seq(min(Time), max(Time), length.out = num.intervals + 1))
    points = seq(min(Time), max(Time), length.out = num.intervals + 1)
  
  cat("Prediction Points = ")
  print(points)
  cat("\n")
  
  ## Filter out the taxa that always have zero of one/both group
  Group = as.character(Group)
  group.levels = sort(unique(Group))
  if(length(group.levels) > 2){
    stop("You have more than two phenotypes.")
  }
  gr.1 = group.levels[1]
  gr.2 = group.levels[2]
  
  
  
  # ## Preprocessing: Remove features that have zeros in all time points of one group 
  # q = Count[,which(ID %in% ID[which(Group == gr.1)])]
  # w = Count[,which(ID %in% ID[which(Group == gr.2)])]
  # rm = which(apply(q, 1, sum) == 0 | apply(w, 1, sum) == 0)
  # 
  # if(length(rm) == 0)
  # {
  #   data.count.filt = Count
  # } else
  # {
  #   data.count.filt = Count[-rm, ]
  # }
  # 
  # data.count.filt = as.matrix(data.count.filt)
  
  data.count.filt = as.matrix(Count)
  
  
  
  ## Apply metalonda for each feature
  n.features = nrow(data.count.filt)
  detailed = list()
  summary = list()
  for (i in 1:n.features)
  {
    cat ("Feature  = ", rownames(data.count.filt)[i], "\n")
    out = metalonda(Count = data.count.filt[i,], Time = Time, Group = Group, ID = ID,
                          fit.method = fit.method, n.perm = n.perm, points = points, 
                          text=rownames(data.count.filt)[i], parall = parall, pvalue.threshold, adjust.method, time.unit)
    
    detailed[[i]] = out$detailed
    summary[[i]] = out$summary
  }
  
  summary.tmp = do.call(rbind, summary)
  summary.tmp$dominant[which(summary.tmp$dominant == 1)] = gr.1
  summary.tmp$dominant[which(summary.tmp$dominant == -1)] = gr.2
  
  ## Output table and figure that summarize the significant time intervals
  write.csv(summary.tmp, file = sprintf("%s_MetaLonDA_TimeIntervals.csv", prefix), row.names = FALSE)
  visualizeTimeIntervals(interval.details = summary.tmp, prefix, unit = time.unit)
  
  return(list(output.detail = detailed, output.summary = summary.tmp))
}
