#' Plot Pulses
#'
#' Plots pulses in 2D, 3D waterfall or sequence based on users choice.
#'
#' This function combines pulses in a data frame that can be used graphically.
#' Several inputs are based on the user, for example, normalization of the pulses, color ramping and axis labels.
#' Options on the type of plot include "2D", "3D", "Sequence" and "Heatmap"
#'
#'
#' @param xPulses A list of data frames containing the pulses to be plotted.
#' @param xTime A list of vectors of time values for xPulses.
#' @param yPulses An optional data frame for values to be plotted against xPulses.
#' @param yTime An optional time for yPulses.  This must be present if yPulses are used.
#' @param normalization A string that indicates the normalization type.  Options: "none", "height", "area".
#' @param newXlab The x axis label.
#' @param newYlab The y axis label.
#' @param colRamp A boolean value if the pulses are colored by a single color or matlab like gradient.
#' @param showLegend A boolean value to display the legend.
#' @param plotType A character describing the plot type.  Options include: "2D", "3D", "Sequence" and "Heatmap"
#' @return A plot of the pulses.
#' @examples
#' data("pulseData")
#'
#' xPulses = pulseData$`32AMU`$pulses
#' xTime = pulseData$`32AMU`$time
#' plotPulse(xPulses, xTime)
#'
#'
#'
#'
#'
#'
#' @seealso
#' \code{\link{pulseNorm}}
#' @export plotPulse

plotPulse = function(xPulses, xTime, yPulses = NULL, yTime = NULL, newXlab = "", newYlab = "", colRamp = F, showLegend = F, plotType = "2D"){

  lineWidth = 2
  cexSize = 2
  pulseNames = unlist(lapply(xPulses,colnames),F,F)



  # if(normalization == "height"){
  #   for(i in 1:length(xPulses)){
  #     tempPulses = xPulses[[i]]
  #     if(is.null(dim(tempPulses))){
  #       tempPulses = tempPulses / max(tempPulses, na.rm = T)
  #     }else{
  #       for(j in 1:dim(tempPulses)[2]){
  #         tempPulses[,j] = tempPulses[,j] / max(tempPulses[,j], na.rm = T)
  #       }
  #     }
  #     xPulses[[i]] = tempPulses
  #   }
  # }else if(normalization == "area"){
  #   for(i in 1:length(xPulses)){
  #     tempPulses = xPulses[[i]]
  #     if(is.null(dim(tempPulses))){
  #       tempPulses = tempPulses / TAPcode::moments(tempPulses, na.rm = T)
  #     }else{
  #       for(j in 1:dim(tempPulses)[2]){
  #         tempPulses[,j] = tempPulses[,j] / max(tempPulses[,j], na.rm = T)
  #       }
  #     }
  #     xPulses[[i]] = tempPulses
  #   }
  # }

#
#   if(!(normalization == "none")){
#     for(i in 1:length(xPulses)){
#       xPulses[[i]] = TAPcodeV2::pulseNorm(xPulses[[i]], xTime[[i]], normalization)
#     }
#     if(!is.null(yPulses))
#       yPulses = TAPcodeV2::pulseNorm(yPulses,yTime,normalization)
#   }


  # check for vectors to matrix
  for(i in 1:length(xPulses)){
    if(is.null(dim(xPulses[[i]]))){
      xPulses[[i]] = as.matrix(xPulses[[i]])
    }
  }

  if(plotType == "Heatmap"){
    xPulsesDim = lapply(xPulses,dim)
    minDim = rep(0,length(xPulses))
    for(i in 1:length(xPulses)){
      minDim[i] = xPulsesDim[[i]][1]
    }
    minDim = min(minDim)
    catMat = NULL
    for(i in 1:length(xPulses)){
      catMat = rbind(catMat, t(as.matrix(xPulses[[i]][1:minDim,])))
    }
    heatmap(catMat, Rowv = NA, Colv = NA, scale = "none", xlab = newXlab, ylab = newYlab, col = colorRamps::matlab.like(256))
  }else{
    if(colRamp){
      pal = grDevices::rainbow(sum(unlist(lapply(xPulses,ncol)),F,F) + 1)
    }else{
      palfun1 = colorRampPalette(c("royalblue", "paleturquoise"))
      palfun2 = colorRampPalette(c("seagreen", "springgreen"))
      palfun3 = colorRampPalette(c("gold4","yellow"))
      palfun4 = colorRampPalette(c("darkred", "pink"))
      palTracker = 1
      pal = "0"
      for(i in 1:length(xPulses)){
        tempPalFun = switch(palTracker,
                            palfun1,
                            palfun2,
                            palfun3,
                            palfun4)
        pal = c(pal,tempPalFun(dim(xPulses[[i]])[2] ))

        if(palTracker == 4){
          palTracker = 1
        }else{
          palTracker = palTracker + 1
        }
      }
      pal = pal[-1]
    }
    tempColor = pal


    if(is.null(yPulses)){
      minXaxis = min(unlist(lapply(xTime,min, na.rm = T),F,F), na.rm = T)
      maxXaxis = max(unlist(lapply(xTime,max, na.rm = T),F,F), na.rm = T)
      minYaxis = min(unlist(lapply(xPulses,min, na.rm = T),F,F), na.rm = T)
      maxYaxis = max(unlist(lapply(xPulses,max, na.rm = T),F,F), na.rm = T)
    }else{
      minXaxis = min(yPulses, na.rm = T)
      maxXaxis = max(yPulses, na.rm = T)
      minYaxis = min(unlist(lapply(xPulses,min),F,F), na.rm = T)
      maxYaxis = max(unlist(lapply(xPulses,max),F,F), na.rm = T)
    }

    if(is.null(yPulses)){
      numYpoints = sum(unlist(lapply(xPulses,ncol))  * unlist(lapply(xPulses,nrow),F,F))
    }else{
      numYpoints = sum(apply(yPulses,2,max))*length(xPulses)
    }



    if(showLegend){
      legScale = .15
    }else{
      legScale = 0
    }

    if(plotType != "Sequence"){
      plot(c(minXaxis,maxXaxis + (maxXaxis - minXaxis)*legScale), c(minYaxis,maxYaxis),
           type = "n", xlab = newXlab, ylab = newYlab, cex.lab = cexSize, cex.axis = cexSize, cex.main = cexSize, cex.sub = cexSize)
    }else{
      plot(c(1,numYpoints +  numYpoints * legScale), c(minYaxis,maxYaxis),
           type = "n", xlab = newXlab, ylab = newYlab, cex.lab = cexSize, cex.axis = cexSize, cex.main = cexSize, cex.sub = cexSize)
    }


    if(is.null(yPulses)){
      palFlag = 1
      timeFlag = 1
      if(plotType == "3D"){
        for(i in 1:length(xPulses)){
          for(j in 1:dim(xPulses[[i]])[2]){
            subPal = colorRampPalette(c(pal[palFlag], pal[palFlag + 1]))
            lines( (xTime[[i]]/2 + (palFlag)/length(pal) * maxXaxis/2  ),
                   (xPulses[[i]][,j]/2  + (palFlag)/length(pal) * maxYaxis/2 ),
                   col =   subPal(length(xTime[[i]])), lwd = lineWidth)
            palFlag = palFlag + 1
          }
        }
      }else if(plotType == "Sequence"){
        for(i in 1:length(xPulses)){
          for(j in 1:dim(xPulses[[i]])[2]){
            lines( seq(timeFlag,(timeFlag + length(xTime[[i]])-1)), xPulses[[i]][,j], col = pal[palFlag], lwd = lineWidth)
            timeFlag = timeFlag + length(xTime[[i]])
            palFlag = palFlag + 1
          }
        }
      }else{
        for(i in 1:length(xPulses)){
          for(j in 1:dim(xPulses[[i]])[2]){
            lines(xTime[[i]], xPulses[[i]][,j], col = pal[palFlag], lwd = lineWidth)
            palFlag = palFlag + 1
          }
        }
      }

    }else{
      palFlag = 1
      timeFlag = 0
      if(plotType == "3D"){
        for(i in 1:length(xPulses)){
          for(j in 1:dim(xPulses[[i]])[2]){
            minLen = min(c(length(yPulses[,j]), length(xPulses[[i]][,j])))
            lines( (yPulses[1:minLen,j]/2 + (palFlag)/length(pal) * maxXaxis/2  ),
                   (xPulses[[i]][1:minLen,j]/2  + (palFlag)/length(pal) * maxYaxis/2 ), col = pal[palFlag], lwd = lineWidth)
            palFlag = palFlag + 1
          }
        }
      }else if(plotType == "Sequence"){
        for(i in 1:length(xPulses)){
          for(j in 1:dim(xPulses[[i]])[2]){
            minLen = min(c(length(yPulses[,j]), length(xPulses[[i]][,j])))
            lines(yPulses[1:minLen,j] + timeFlag  , xPulses[[i]][1:minLen,j], col = pal[palFlag], lwd = lineWidth)
            timeFlag = timeFlag + max(yPulses[,j])
            palFlag = palFlag + 1
          }
        }
      }else{
        for(i in 1:length(xPulses)){
          for(j in 1:dim(xPulses[[i]])[2]){
            minLen = min(c(length(yPulses[,j]), length(xPulses[[i]][,j])))
            lines(yPulses[1:minLen,j], xPulses[[i]][1:minLen,j], col = pal[palFlag], lwd = lineWidth)
            palFlag = palFlag + 1
          }
        }
      }
    }

    if(plotType == "seq" ){
      legMax = numYpoints
    }else{
      legMax = maxXaxis
    }

    if(showLegend){
      legend("topright",
             pulseNames,
             lty=rep(1,length(pulseNames)),
             lwd=rep(2.5,length(pulseNames)),
             col=tempColor)
    }
  }




}
