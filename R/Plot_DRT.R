PlotDRT <- function (name, reps, burnin = 0.1) {
  
 #library(TeachingDemos)
 #library(DescTools) 
 
  opar <-par(no.readonly=TRUE)
  on.exit(par(opar))

  ## Read the log
  readFileName <- paste0(name, ".log")
  comp0 <- read.table(readFileName, header = TRUE, sep = '\t')
  comp0a = c()
  while(length(comp0a) < 1 ) {value <- readline(
    "enter parameter, type VAR to list variable names, or hit enter to cancel: ")
  comp0a <- as.data.frame(cbind(comp0[which(names(comp0) == value)]))
  if (value == "var") {print(colnames(comp0, prefix = "col"))}
  if (value == "VAR") {print(colnames(comp0, prefix = "col"))}
  if (value == "") {break} # breaks when hit enter
  }
  
  comp0a$calibr <- "0"
  comp0b <- tail(comp0a, round(dim(comp0a)[1] - (dim(comp0a)[1] * burnin)))
  newdat <- comp0b

  ## Extract values from replicates
  for (i in 1 : reps){
    
    ## Read replicates log files
    fileRep <- paste0(name, ".Rep", i, ".log")  
    compName <- paste0("comp", i)
    compName <- paste0("comp", i, "a")
    temp <- read.table(fileRep, header = TRUE, sep = '\t')
    tempa <- as.data.frame(cbind(temp[which(names(temp) == value)]))
    names (tempa) <- c(value)
    tempa$calibr <- as.character(i)
    
    ##  remove burnin
    tempb <- tail(tempa, round(dim(tempa)[1] - (dim(tempa)[1] * burnin)))
    assign(compName, tempb)
    newdat <- rbind(newdat, tempb)
    print (paste("log of replicate", i, "processed"))
  }
  newdat$calibr <- as.factor(newdat$calibr)
  
  ##  MAKE GRAPHS
  
  # compute a matrix with mean, median, lower and higher 95 %HPD on height stat
  stat = matrix(NA, length(unique(newdat$calibr)), 4)
  colnames(stat) = c("mean","median","lowerHPD","HigherHPD")
  calibr = seq(0, length(unique(newdat$calibr))-1, 1)
  stats = cbind(calibr, stat)
  for (i in 1: dim(stats)[1]){
    stats[i, 2] = mean(newdat[newdat$calibr == i-1, 1])
    stats[i, 3] = median(newdat[newdat$calibr == i-1, 1])
    stats[i, 4] = emp.hpd(newdat[newdat$calibr == i-1, 1])[1]
    stats[i, 5] = emp.hpd(newdat[newdat$calibr == i-1, 1])[2]
  }
  
  # Export the data into a text table that summarizes the results
  write.table (stats, file = paste0(value,".stats.csv"), sep=",")

  min = min(stats[, 4])
  max = max(stats[, 5])
  max_with_margin_for_legend = max + (max * 0.1)
  par(mfrow=c(2, 1))
  par(mar = rep(2, 4))
  plot (stats[1, 1], stats[1, 3], xlim = c(0, dim(stats)[1]), 
        ylim = c(min, max_with_margin_for_legend), 
        pch=16,col="red", xaxt="n", xlab="", ylab = value, 
        main = paste0("Date-randomization test performed on ", value))
  points (stats[1, 1], stats[1, 4], pch = '-', col = "red")
  points (stats[1, 1], stats[1, 5], pch = '-', col = "red")
  segments (stats[1, 1], stats[1, 4], stats[1, 1], stats[1, 5], col = "red")
  for (i in 2: dim(stats)[1]){
    points(stats[i, 1], stats[i, 3], pch = 16, col = "black")
    points(stats[i, 1], stats[i, 4], pch = '-', col = "black")
    points(stats[i, 1], stats[i, 5], pch = '-', col = "black")
    segments(stats[i, 1], stats[i, 4], stats[i, 1], stats[i, 5], col = "black")
  }
  abline(h = stats[1,4], col = "lightcoral", lty = 5,lwd = 1.5)
  abline(h = stats[1,5], col = "lightcoral", lty = 5,lwd = 1.5)
  
  legend("top",c ("Real", "Randomized"),col = c("red", "black"), lty = 1, 
         pch = 16, cex = 1, horiz = TRUE, bty = "n")
  
  # compute statistic
  
  minmax = vector("list", dim(stats)[1])
  for (i in 1: dim(stats)[1]){
    minmax [[i]] = c(stats[i, 4], stats[i, 5])
  }
  # Test each interval against true dataset and create a vector 
  # with a value per randomized dataset filled with 0 if no overlap and a number
  # (the width of overlap) if overlap
  
  is.overlap = NULL
  for (i in 2: dim(stats)[1]){
    new <- Overlap(minmax[[1]], minmax[[i]])
    is.overlap <- append(is.overlap, new)
  }
  DRT_result <- all(is.overlap == 0)
  overlapping_with_true_dataset <- which(is.overlap != 0)
  # Print result under the plot
  if (DRT_result == TRUE) {
    mtext("No Overlapping: DRT SUCCESSFULLY PASSED !!!", side = 1,
        line = 1, col = "red")
  } else {
    mtext(paste("Overlapping: DRT FAILED !!!: date randomized dataset(s)",
                      paste(overlapping_with_true_dataset, collapse = ", "), 
                      "is/are overlapping with true dataset"), side = 1,
                line = 1, col = "red")
  }

# plot DRT log10 rate
  plot (stats[1, 1], log10(stats[1, 3]), xlim = c(0, dim(stats)[1]), 
        ylim = c(log10(min), log10(max_with_margin_for_legend)),
        pch = 16, col = "red", xaxt = "n", xlab = "",
        ylab = paste("Log10", value),
        main = paste0("Date-randomization test performed on ", value))
  points (stats[1, 1], log10(stats[1, 4]), pch = '-', col = "red")
  points (stats[1, 1], log10(stats[1, 5]), pch = '-', col = "red")
  segments (stats[1, 1], log10(stats[1, 4]), stats[1, 1], log10(stats[1, 5]),
            col = "red")
  
  for (i in 2: dim(stats)[1]){
    points(stats[i, 1], log10(stats[i, 3]), pch = 16, col = "black")
    points(stats[i, 1], log10(stats[i, 4]), pch = '-', col = "black")
    points(stats[i, 1], log10(stats[i, 5]), pch = '-', col = "black")
    segments(stats[i, 1], log10(stats[i, 4]), stats[i, 1], log10(stats[i, 5]),
             col = "black")
  }
  abline(h = log10(stats[1,4]), col = "lightcoral", lty = 5,lwd = 1.5)
  abline(h = log10(stats[1,5]), col = "lightcoral", lty = 5,lwd = 1.5)
  
  legend("top", c ("Real", "Randomized"),col = c("red", "black"), lty = 1, 
         pch = 16, cex = 1, horiz = TRUE, bty = "n")
  DRT_result <- all(is.overlap == 0)
  overlapping_with_true_dataset <- which(is.overlap != 0)
  # Print result under the plot
  if (DRT_result == TRUE) {
    mtext("No Overlapping: DRT SUCCESSFULLY PASSED !!!", side = 1, line = 1,
          col = "red")
  } else {
    mtext(paste("Overlapping: DRT FAILED !!!: date randomized dataset(s)", 
                      paste(overlapping_with_true_dataset, collapse = ", "), 
                      "is/are overlapping with true dataset"), side = 1,
                line = 1, col = "red")
  }
  
  pdf(paste0("Fig_DRT_", value, ".pdf"))

 tidy = F
  
  plot (stats[1, 1], stats[1, 3], xlim = c(0, dim(stats)[1]), 
        ylim = c(min, max_with_margin_for_legend), 
        pch=16,col="red", xaxt="n", xlab="", ylab = value, 
        main = paste0("Date-randomization test performed on ", value))
  points (stats[1, 1], stats[1, 4], pch = '-', col = "red")
  points (stats[1, 1], stats[1, 5], pch = '-', col = "red")
  segments (stats[1, 1], stats[1, 4], stats[1, 1], stats[1, 5], col = "red")
  
  for (i in 2: dim(stats)[1]){
    points(stats[i, 1], stats[i, 3], pch = 16, col = "black")
    points(stats[i, 1], stats[i, 4], pch = '-', col = "black")
    points(stats[i, 1], stats[i, 5], pch = '-', col = "black")
    segments(stats[i, 1], stats[i, 4], stats[i, 1], stats[i, 5], col = "black")
  }
  abline(h = stats[1,4], col = "lightcoral", lty = 5,lwd = 1.5)
  abline(h = stats[1,5], col = "lightcoral", lty = 5,lwd = 1.5)
  legend("top",c ("Real", "Randomized"),col = c("red", "black"), lty = 1, 
         pch = 16, cex = 1, horiz = TRUE, bty = "n")

  minmax = vector("list", dim(stats)[1])
  for (i in 1: dim(stats)[1]){
    minmax [[i]] = c(stats[i, 4], stats[i, 5])
  }

  is.overlap = NULL
  for (i in 2: dim(stats)[1]){
    new <- Overlap(minmax[[1]], minmax[[i]])
    is.overlap <- append(is.overlap, new)
  }
  DRT_result <- all(is.overlap == 0)
  overlapping_with_true_dataset <- which(is.overlap != 0)
  # Print result under the plot
  if (DRT_result == TRUE) {
    mtext("No Overlapping: DRT SUCCESSFULLY PASSED !!!", side = 1, line = 1,
          col = "red")
  } else {
    mtext(paste("Overlapping: DRT FAILED !!!: date randomized dataset(s)", 
                      paste(overlapping_with_true_dataset, collapse = ", "),
                      "is/are overlapping with true dataset"), side = 1,
                line = 1, col = "red")
  }

  dev.off()
  
  pdf(paste0("Fig_DRT_log_", value, ".pdf"))

 tidy = F

# plot DRT log10 rate

  plot (stats[1, 1], log10(stats[1, 3]), xlim = c(0, dim(stats)[1]), 
        ylim = c(log10(min), log10(max_with_margin_for_legend)),
        pch = 16, col = "red", xaxt = "n", xlab = "",
        ylab = paste("Log10", value),
        main = paste0("Date-randomization test performed on ", value))
  points (stats[1, 1], log10(stats[1, 4]), pch = '-', col = "red")
  points (stats[1, 1], log10(stats[1, 5]), pch = '-', col = "red")
  segments (stats[1, 1], log10(stats[1, 4]), stats[1, 1], log10(stats[1, 5]),
            col = "red")
  
  for (i in 2: dim(stats)[1]){
    points(stats[i, 1], log10(stats[i, 3]), pch = 16, col = "black")
    points(stats[i, 1], log10(stats[i, 4]), pch = '-', col = "black")
    points(stats[i, 1], log10(stats[i, 5]), pch = '-', col = "black")
    segments(stats[i, 1], log10(stats[i, 4]), stats[i, 1], log10(stats[i, 5]),
             col = "black")
  }
  abline(h = log10(stats[1,4]), col = "lightcoral", lty = 5,lwd = 1.5)
  abline(h = log10(stats[1,5]), col = "lightcoral", lty = 5,lwd = 1.5)
  
  legend("top", c ("Real", "Randomized"),col = c("red", "black"), lty = 1, 
         pch = 16, cex = 1, horiz = TRUE, bty = "n")
  
  DRT_result <- all(is.overlap == 0)
  overlapping_with_true_dataset <- which(is.overlap != 0)
  if (DRT_result == TRUE) {
    mtext("No Overlapping: DRT SUCCESSFULLY PASSED !!!", side = 1, line = 1,
          col = "red")
  } else {
    mtext(paste("Overlapping: DRT FAILED !!!: date randomized dataset(s)", 
                      paste(overlapping_with_true_dataset, collapse = ", "),
                      "is/are overlapping with true dataset"), side = 1,
                line = 1, col = "red")
  }
  
  dev.off()
  
}
