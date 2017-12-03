#' Visualize Longitudinal Feature
#'
#' Visualize Longitudinal Feature
#'
#' @param df dataframe has the Count, Group, ID, Time
#' @param text feature name
#' @param group.levels The two level's name
#' @param unit time interval unit
#' @import ggplot2
#' @import grDevices
#' @import graphics
#' @references
#' Ahmed Metwally (ametwa2@uic.edu)
#' @examples 
#' data(metalonda_test_data)
#' n.sample = 5 # sample size;
#' n.timepoints = 10 # time point;
#' n.group= 2 # number of group;
#' Group = factor(c(rep(0,n.sample*n.timepoints), rep(1,n.sample*n.timepoints)))
#' Time = rep(rep(1:n.timepoints, times = n.sample), 2)
#' ID = factor(rep(1:(2*n.sample), each = n.timepoints))
#' points = seq(1, 10, length.out = 10)
#' aggregate.df = data.frame(Count = metalonda_test_data[1,], Time = Time, Group = Group, ID = ID)
#' visualizeFeature(aggregate.df, text = rownames(metalonda_test_data)[1], Group)
#' @export
visualizeFeature = function (df, text, group.levels, unit = "days")
{
  cat("Visualizing Feature = ", text, "\n")
  Count=0; Time=0; ID=0; Group=0 ## This line is just to pass the CRAN checks for the aes in ggplot2
  
  p = ggplot(df, aes(Time, Count, colour = Group, group = interaction(Group, ID)))
  p = p + geom_point(size=1, alpha=0.5) + geom_line(size=1, alpha=0.7) +  theme_bw() +
    ggtitle(paste("Feature = ", text, sep = "")) + labs(y = "Normalized Count", x = sprintf("Time (%s)", unit)) +
    scale_colour_manual(values = c("skyblue", "pink"), breaks = c("0", "1"),
                        labels = c(group.levels[1], group.levels[2])) +
    theme(axis.text.x = element_text(colour="black", size=12, angle=0, hjust=0.5, vjust=0.5, face="bold"),
          axis.text.y = element_text(colour="black", size=12, angle=0, hjust=0.5, vjust=0.5, face="bold"),
          axis.title.x = element_text(colour="black", size=15, angle=0, hjust=.5, vjust=0.5, face="bold"),
          axis.title.y = element_text(colour="black", size=15, angle=90, hjust=.5, vjust=.5, face="bold"),
          legend.text=element_text(size=15, face="plain"), legend.title = element_blank()) +
    theme(legend.position="top") + scale_x_continuous(breaks = waiver())

  ggsave(filename=paste("Feature_", text, ".tiff", sep=""), dpi = 300, height = 10, width = 15, units = 'cm')
}



#' Visualize the feature trajectory with the fitted Splines
#'
#' Plot the longitudinal features along with the fitted splines
#'
#' @param df dataframe has the Count , Group, ID, Time
#' @param model the fitted model
#' @param method The fitting method (negative binomial, LOWESS)
#' @param group.levels The two level's name
#' @param text feature name
#' @param unit time unit used in the Time vector (days, weeks, months, etc.)
#' @import ggplot2
#' @import grDevices
#' @import graphics
#' @references
#' Ahmed Metwally (ametwa2@uic.edu)
#' @export
visualizeFeatureSpline = function (df, model, method, text, group.levels, unit = "days")
{ 
  cat("Visualizing the Splines of Feature =  ", text, "\n")
    
  Count=0;Time=0;ID=0;Group=0;lnn=0 ## This line is just to pass the CRAN checks for the aes in ggplot2
  dd.null = model$dd.null
  dd.0 = model$dd.0
  dd.1 = model$dd.1
  
  ln = factor(c(rep("longdash", nrow(df)), rep("longdash", nrow(dd.0)), rep("longdash", nrow(dd.1))))
  size = c(rep(1, nrow(df)), rep(1, nrow(dd.0)), rep(1, nrow(dd.1)))
  dm = rbind(df[,c("Time", "Count", "Group", "ID")], dd.0, dd.1)
  dm$lnn=ln
  dm$sz= size
  
  p = ggplot(dm, aes(Time, Count, colour = Group, group = interaction(Group, ID)))
  p = p + theme_bw() + geom_point(size=1, alpha=0.5) + geom_line(aes(linetype=lnn), size=1, alpha=0.5) + 
    ggtitle(paste("Feature = ", text, " '", method, "'", sep = "")) + labs(y = "Normalized Count", x = sprintf("Time (%s)", unit)) +
    scale_colour_manual(values = c("skyblue", "pink", "blue", "firebrick",
                                   "blue",  "blue",  "firebrick", "firebrick"), 
                        breaks = c("0", "1", "fit.0", "fit.1"),
                        labels = c(group.levels[1], group.levels[2], paste(group.levels[1], ".fit", sep=""), paste(group.levels[2], ".fit", sep="")))+
    theme(axis.text.x = element_text(colour="black", size=12, angle=0, hjust=0.5, vjust=0.5, face="bold"),
          axis.text.y = element_text(colour="black", size=12, angle=0, hjust=0.5, vjust=0.5, face="bold"),
          axis.title.x = element_text(colour="black", size=15, angle=0, hjust=.5, vjust=0.5, face="bold"),
          axis.title.y = element_text(colour="black", size=15, angle=90, hjust=.5, vjust=.5, face="bold"), 
          legend.text=element_text(size=15, face="plain"), legend.title = element_blank()) +
    theme(legend.position="top") + scale_x_continuous(breaks = waiver()) + guides(linetype=FALSE, size =FALSE)
  
  ggsave(filename=paste("Feature_", text, "_CurveFitting_", method, ".tiff", sep=""), dpi = 300, height = 10, width = 15, units = 'cm')
}


#' Visualize Area Ratio (AR) empirical distribution
#'
#' Visualize Area Ratio (AR) empirical distribution for each time interval
#'
#' @param permuted Permutation of the permuted data
#' @param text Feature name
#' @param method fitting method
#' @import ggplot2
#' @import grDevices
#' @import graphics
#' @references
#' Ahmed Metwally (ametwa2@uic.edu)
#' @export
visualizeARHistogram = function(permuted, text, method){
  cat("Visualizing AR distribution of each time interval for feature = ", text, "\n")
  n = ncol(permuted)
  r = ceiling(sqrt(n))
  c = ceiling(sqrt(n))
	tiff(paste("Feature_", text, "_AR_distribution_", method, ".tiff", sep = ""), res = 300, height = r*5, width = c*5, units = 'cm')
  par(mfrow=c(r,c))
  
  for( i in 1:ncol(permuted)){
    hist(permuted[,i], xlab = "AR Ratio", ylab = "Frequency", 
         breaks = 10, col = "yellow", border = "red", 
         main = paste("Interval # ", i, sep=""), xlim = c(0,1))
  }
  dev.off()
}



#' Visualize significant time interval
#'
#' Visualize significant time interval
#'
#' @param aggregate.df Dataframe has the Count, Group, ID, Time
#' @param model.ss The fitted model
#' @param method Fitting method (negative binomial or LOWESS)
#' @param start Vector of the start points of the time intervals
#' @param end Vector of the end points of the time intervals
#' @param text Feature name
#' @param group.levels Level's name
#' @param unit time unit used in the Time vector (days, weeks, months, etc.)
#' @import ggplot2
#' @import grDevices
#' @import graphics
#' @references
#' Ahmed Metwally (ametwa2@uic.edu)
#' @export
visualizeArea = function(aggregate.df, model.ss, method, start, end, text, group.levels, unit = "days")
{
  cat("Visualizing significant intervals of feature = ", text, "\n")
  Time = 0 ## This line is just to pass the CRAN checks for the aes in ggplot2
  sub.11 = list()
  sub.10 = list()
  xx = NULL
  for(i in 1:length(start))
  {
    sub.11[[i]] = subset(model.ss$dd.1, Time >= start[i] & Time <= end[i])  
    sub.10[[i]] = subset(model.ss$dd.0, Time >= start[i] & Time <= end[i])
    cmd = sprintf('geom_ribbon(data=sub.10[[%d]], aes(ymin = sub.11[[%d]]$Count, ymax = Count), colour= "grey3", fill="grey69", 
                  alpha = "0.6")', i, i)
    if (i != 1)
    {
      xx = paste(xx, cmd, sep = "+")
    } else
    {
      xx = cmd
    }
  }
  
  # ddNULL = model_ss$ddNULL
  dd.0 = model.ss$dd.0
  dd.1 = model.ss$dd.1
  
  dm = rbind(dd.0, dd.1)
  p1 = 'ggplot(dm, aes(Time, Count, colour = Group, group = interaction(Group, ID))) + 
  theme_bw() + geom_point(size = 1, alpha = 0.5) + geom_line(size = 1, alpha = 0.5) + 
  ggtitle(paste("Feature = ", text, sep = "")) + labs(y = "Normalized Count", x = sprintf("Time (%s)", unit)) +
  scale_colour_manual(values = c("blue", "firebrick"), 
  breaks = c("fit.0", "fit.1"),
  labels = c(paste(group.levels[1], ".fit", sep = ""), paste(group.levels[2], ".fit", sep = ""))) +
  theme(axis.text.x = element_text(colour = "black", size = 12, angle = 0, hjust = 0.5, vjust = 0.5, face = "bold"),
  axis.text.y = element_text(colour = "black", size = 12, angle = 0, hjust = 0.5, vjust = 0.5, face = "bold"),
  axis.title.x = element_text(colour = "black", size = 15, angle = 0, hjust = 0.5, vjust = 0.5, face = "bold"),
  axis.title.y = element_text(colour = "black", size = 15, angle = 90, hjust = 0.5, vjust = 0.5, face = "bold"), 
  legend.text = element_text(size = 15, face="plain"), legend.title = element_blank()) +
  theme(legend.position = "top") + scale_x_continuous(breaks = waiver())' 
  p2 = xx  
  p3 = paste(p1, p2, sep="+")
  p = eval(parse(text = p3))
  ggsave(filename=paste("Feature_", text, "_SignificantInterval_", method, ".tiff", sep=""), dpi = 300, height = 10, width = 15, units = 'cm')
}



#' Visualize all significant time intervals for all tested features
#'
#' Visualize all significant time intervals for all tested features
#'
#' @param interval.details Dataframe has infomation about significant interval (feature name, start, end, dominant, p-value)
#' @param prefix prefix for the output figure
#' @import ggplot2
#' @import grDevices
#' @import graphics
#' @references
#' Ahmed Metwally (ametwa2@uic.edu)
#' @export
visualizeTimeIntervals = function(interval.details, prefix = "MetaLonDA_timeline")
{
  feature=0;dominant=0;ID=0;Group=0;lnn=0 ## This line is just to pass the CRAN checks for the aes in ggplot2
  interval.details$dominant = as.factor(interval.details$dominant)
  interval.details$pvalue = as.numeric((interval.details$pvalue))
  interval.details = interval.details[order(interval.details$feature), ]
  
  
  ggplot(interval.details, aes(ymin = start , ymax = end, x = feature, xend = feature)) + 
    geom_linerange(aes(color = dominant), size = 1) + 
    coord_flip() +  scale_colour_manual(values=c("firebrick", "blue")) +
    labs(x = "Feature", y = "Time (Days)", colour="Dominant") + 
     theme(axis.text.x = element_text(colour = "black", size = 10, angle = 0, hjust = 0.5, vjust = 0.5, face = "bold"),
           axis.text.y = element_text(colour = "black", size = 8, angle = 0, vjust = 0.5, face = "bold"),
           axis.title.x = element_text(colour = "black", size = 15, angle = 0, hjust = 0.5, vjust = 0.5, face = "bold"),
           axis.title.y = element_text(colour = "black", size = 15, angle = 90, hjust = 0.5, vjust = 0.5, face = "bold"),
           legend.text = element_text(size = 15, face = "plain")) + 
    theme(panel.grid.minor =   element_blank(),
          panel.grid.major.y = element_line(colour = "white", size = 6),
          panel.grid.major.x = element_line(colour = "white",size = 0.75)) +
    theme(legend.position="top", panel.border = element_rect(colour = "black", fill = NA, size = 2))
  ggsave(filename = paste(prefix, "_TimeIntervals_", ".tiff", sep=""), dpi = 300, height = 30, width = 20, units = 'cm')
}
