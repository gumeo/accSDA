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
#' @param xlabel Label to put under every plot
#' @param ylabel Vector of y-axis labels for each plot, e.g. if there are three DVs, then
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
#'     # Generate and ASDA object with your data, e.g.
#'     # Prepare training and test set
#'     # This is a very small data set, I advise you to try it on something with more
#'     # variables, e.g. something from this source: http://www.cs.ucr.edu/~eamonn/time_series_data/
#'     # or possibly run this on the Gaussian data example from the ASDA function
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
#'     ASDABarPlot(resIris)
#' @rdname ASDABarPlot
#' @export ASDABarPlot
ASDABarPlot <- function(asdaObj, numDVs = 1, xlabel, ylabel, getList = FALSE, main, ...){
  # Found on: https://rpubs.com/Koundy/71792
  theme_Publication <- function(base_size=14, base_family="Helvetica") {
    (ggthemes::theme_foundation(base_size=base_size, base_family=base_family)
     + ggplot2::theme(plot.title = ggplot2::element_text(face = "bold",
                                       size = ggplot2::rel(1.2), hjust = 0.5),
             text = ggplot2::element_text(),
             panel.background = ggplot2::element_rect(colour = NA),
             plot.background = ggplot2::element_rect(colour = NA),
             panel.border = ggplot2::element_rect(colour = NA),
             axis.title = ggplot2::element_text(face = "bold",size = ggplot2::rel(1)),
             axis.title.y = ggplot2::element_text(angle=90,vjust =2),
             axis.title.x = ggplot2::element_text(vjust = -0.2),
             axis.text = ggplot2::element_text(),
             axis.line = ggplot2::element_line(colour="black"),
             axis.ticks = ggplot2::element_line(),
             panel.grid.major = ggplot2::element_line(colour="#f0f0f0"),
             panel.grid.minor = ggplot2::element_blank(),
             legend.key = ggplot2::element_rect(colour = NA),
             legend.position = "bottom",
             legend.direction = "horizontal",
             legend.key.size= ggplot2::unit(0.2, "cm"),
             legend.spacing = ggplot2::unit(0, "cm"),
             legend.title = ggplot2::element_text(face="italic"),
             plot.margin=ggplot2::unit(c(10,5,5,5),"mm"),
             strip.background=ggplot2::element_rect(colour="#f0f0f0",fill="#f0f0f0"),
             strip.text = ggplot2::element_text(face="bold")
     ))

  }
  DVs <- as.data.frame(asdaObj$beta)
  if(numDVs)
    if(missing(ylabel) || length(ylabel)!=numDVs){
      colnames(DVs) <- paste('DV',1:numDVs, sep='')
      ylabel <- colnames(DVs)
    }else{
      colnames(DVs) <- ylabel
    }
  if(missing(xlabel)){
    xLab <- "Feature number"
  }else{
    xLab <- xlabel
  }
  DVs$featureNum <- 1:dim(DVs)[1]
  plotList <- list()
  for(i in 1:numDVs){
    ggObj <- ggplot2::ggplot(data = DVs, ggplot2::aes_string(x = 'featureNum', y=colnames(DVs)[i]))+
      ggplot2::geom_bar(stat='identity',position='identity') +
      ggplot2::ylab(ylabel[i]) +
      theme_Publication() +
      ggplot2::geom_hline(ggplot2::aes(yintercept=0))
    if(i == numDVs || getList){
      ggObj <- ggObj + ggplot2::xlab(xLab)
    }else{
      ggObj <- ggObj + ggplot2::xlab('')
    }
    plotList[[i]] <- ggObj
  }
  if(missing(main)){
    # Calculate the proportion of nonzero features
    rS <- base::rowSums(abs(DVs[,-base::dim(DVs)[2]]))
    propNZ <- base::sum(rS!=0)/base::dim(DVs)[1]
    main = grid::textGrob(paste('Total sparsity:', round(propNZ,3)), gp=grid::gpar(fontfamily = "Helvetica", fontface = "bold", fontsize = 16))
  }else{
    main = grid::textGrob(main, gp=grid::gpar(fontfamily = "Helvetica", fontface = "bold", fontsize = 16))
  }
  if(getList){
    return(plotList)
  }
  ml <- gridExtra::marrangeGrob(plotList, nrow=numDVs, ncol=1, top=NULL)
  return(gridExtra::grid.arrange(grobs = ml, top = main))
}
