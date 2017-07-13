#' @title barplot for ASDA objects
#'
#' @description This is a function to visualize the discriminant vector from the ASDA
#'     method. The plot is constructed as a ggplot barplot and the main purpose of it is
#'     to visually inspect the sparsity of the discriminant vectors. The main things to
#'     look for are how many parameters are non-zero and if there is any structure in
#'     the ones that are non-zero, but the structure is dependent on the order you specify
#'     your variables. For time-series data, this could mean that a chunk of variables are
#'     non-zero that are close in time, meaning that there is some particular event that is
#'     best for discriminating between the classes that you have.
#'
#' @param asdaObj Object from the \code{ASDA} function.
#' @param numDVs Number of discriminant vectors (DVs) to plot. This is limited by the
#'        number of DVs outputted from the \code{ASDA} function or k-1 DVs where k
#'        is the number of classes. The first 1 to numDVs are plotted.
#' @param xlab Label to put under every plot
#' @param ylab Vector of y-axis labels for each plot, e.g. if there are three DVs, then
#'        \code{ylab = c('Discriminant Vector 1', 'Discriminant Vector 2', 'Discriminant Vector 3')}
#'        is a valid option.
#' @param getList Logical value indicating whether the output should be a list of the plots
#'        or the plots stacked in one plot using the gridExtra package. By default the function
#'        produces a single plot combining all plots of the DVs.
#' @param main Main title for the plots, this is not used if getList is set to \code{TRUE}.
#'
#' @param ... Extra arguments to \code{grid.arrange}.
#'
#' @return \code{barplot.ASDA} returns either a single combined plot or a list of
#'         individual ggplot objects.
#'
#' @seealso \code{\link{ASDA}}
#' @note This function is used as a quick diagnostics tool for the output from the ASDA function.
#'       Feel free to look at the code to customize the plots in any way you like.
#' @examples
#'   \dontrun{
#'     # Generate and ASDA object with your data, e.g.
#'     # Prepare training and test set
#'     # This is a very small data set, I advise you to try it on something with more
#'     # variables, e.g. something from this source: http://www.cs.ucr.edu/~eamonn/time_series_data/
#'     train <- c(1:40,51:90,101:140)
#'     Xtrain <- iris[train,1:4]
#'     nX <- normalize(Xtrain)
#'     Xtrain <- nX$Xc
#'     Ytrain <- iris[train,5]
#'     Xtest <- iris[-train,1:4]
#'     Xtest <- normalizetest(Xtest,nX)
#'     Ytest <- iris[-train,5]
#'     # Run the method
#'     resIris <- ASDA(Xtrain,Ytrain)
#'
#'     # Look at the barplots of the DVs
#'     barplot(resIris)
#'   }
#' @rdname barplot.ASDA
#' @export barplot.ASDA
barplot.ASDA <- function(asdaObj, numDVs = 1, xlabel, ylabel, getList = FALSE, main, ...){
  # Found on: https://rpubs.com/Koundy/71792
  theme_Publication <- function(base_size=14, base_family="helvetica") {
    (theme_foundation(base_size=base_size, base_family=base_family)
     + theme(plot.title = element_text(face = "bold",
                                       size = rel(1.2), hjust = 0.5),
             text = element_text(),
             panel.background = element_rect(colour = NA),
             plot.background = element_rect(colour = NA),
             panel.border = element_rect(colour = NA),
             axis.title = element_text(face = "bold",size = rel(1)),
             axis.title.y = element_text(angle=90,vjust =2),
             axis.title.x = element_text(vjust = -0.2),
             axis.text = element_text(),
             axis.line = element_line(colour="black"),
             axis.ticks = element_line(),
             panel.grid.major = element_line(colour="#f0f0f0"),
             panel.grid.minor = element_blank(),
             legend.key = element_rect(colour = NA),
             legend.position = "bottom",
             legend.direction = "horizontal",
             legend.key.size= unit(0.2, "cm"),
             legend.spacing = unit(0, "cm"),
             legend.title = element_text(face="italic"),
             plot.margin=unit(c(10,5,5,5),"mm"),
             strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
             strip.text = element_text(face="bold")
     ))

  }
  DVs <- as.data.frame(asdaObj$beta)
  if(numDVs)
    if(missing(ylabel) || length(ylab)!=numDVs){
      colnames(DVs) <- paste('DV',1:numDVs, sep='')
      ylabel <- colnames(DVs)
    }else{
      colnames(DVs) <- ylab
    }
  if(missing(xlabel)){
    xLab <- "Feature number"
  }else{
    xLab <- xlabel
  }
  DVs$featureNum <- 1:dim(DVs)[1]
  plotList <- list()
  for(i in 1:numDVs){
    ggObj <- ggplot(data = DVs, aes_string(x = 'featureNum', y=colnames(DVs)[i]))+
      geom_bar(stat='identity',position='identity') +
      ylab(ylabel[i]) +
      theme_Publication() +
      geom_hline(aes(yintercept=0))
    if(i == numDVs || getList){
      ggObj <- ggObj + xlab(xLab)
    }else{
      ggObj <- ggObj + xlab('')
    }
    plotList[[i]] <- ggObj
  }
  if(missing(main)){
    # Calculate the proportion of nonzero features
    rS <- rowSums(abs(DVs[,-dim(DVs)[2]]))
    propNZ <- sum(rS!=0)/dim(DVs)[1]
    main = textGrob(paste('Total sparsity:', round(propNZ,3)), gp=gpar(fontfamily = "helvetica", fontface = "bold", fontsize = 16))
  }else{
    main = textGrob(main, gp=gpar(fontfamily = "helvetica", fontface = "bold", fontsize = 16))
  }
  if(getList){
    return(plotList)
  }
  ml <- marrangeGrob(plotList, nrow=numDVs, ncol=1, top=NULL)
  return(grid.arrange(grobs = ml, top = main))
}
