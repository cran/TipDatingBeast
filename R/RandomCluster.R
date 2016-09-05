RandomCluster <- function(name, reps = 20, loadCluster = T, writeTrees = T) {
  
  #library(mclust)
  
  options (warn = -1)
  
  inFileName <- paste0(name, ".xml")
  inFile <- readLines(inFileName)
  ver2 <- grep(pattern = "version=\"2.0\"", x = inFile, value = F)  
  ver1 <- grep(pattern = "version=\"1.0\"", x = inFile, value = F)  
  ver <- length(ver1) + length(ver2)
  
  if (ver == 1) {

  matchLines <- grep(pattern = "date value=", x = inFile, value = T)
  nameLinesRaw <- grep(pattern = "<taxon id=", x = inFile, value = T)
  nameLines = substring(nameLinesRaw, 14, nchar(nameLinesRaw))
  nameLines = substring(nameLines, 1, nchar(nameLines)-2)  
  if (length(matchLines) == 0) {stop(
    "No dates found, check BEAST input file")}
  min <- 1
  tips <- length(matchLines)
  
  if (loadCluster == F){
    dates <- as.numeric(regmatches(matchLines[], 
                                   gregexpr('\\(?[0-9,.]+', matchLines[])))
    d_clust <- Mclust(as.matrix(as.numeric(dates[])), G = min:tips, modelNames = "E")
    n <- dim(d_clust$z)[2]
    if (n == 1) {stop(
      "One cluster found, input clusters manually or use RandomDates")}
    cat("Number of clusters found:", n, "\n")
    clusterOut <- cbind(nameLines, d_clust$classification)
    nameCluster <- paste0(name, ".clusters.csv")
    write.table(clusterOut, file = nameCluster, 
                row.names=F, col.names=F, sep=",")
    } else {
    nameCluster <- paste0(name, ".clusters.csv")
    clusterOut <- read.table(nameCluster, sep=",")}
  matchLinesPosition <- grep(pattern = "date value=", x = inFile, value = F)
  matchFileName <- grep(pattern = "fileName", x = inFile, value = T)
  matchFileNamePosition <- grep(pattern = "fileName", x = inFile, value = F)
  matchLinesCluster <- cbind(matchLines, clusterOut[,2])
  cluster <- split(matchLinesCluster[,1], matchLinesCluster[,2])
  nCluster <- length(cluster)
  clusterList <- matchLinesCluster[, 2]
  if (loadCluster == T){
    if (nCluster == 1) {stop("Error, only One cluster found, check input")}
    cat("Number of clusters loaded:", nCluster, "\n")}
  # Loops
  for (i in 1 : reps){
    newFile <- inFile

    for (j in 1 : tips){
      clusterL <- clusterList
      clusterId <- clusterOut[j,2]
      if (clusterId != 0){
        clusterL <- clusterL[clusterL != clusterId]
        if ("0" %in% clusterL) {
          clusterL <- clusterL[clusterL != "0"]}
        if ("0" %in% clusterList && nCluster == 2){
          clusterL <- clusterList[clusterList != "0"]}
        as.numeric(clusterL)
        pool <- matchLinesCluster[matchLinesCluster[,2] %in% clusterL,1]
        if (length(pool) == 1) {stop(
          "Single element in cluster, Check input or use RandomDates")}
        newFile[matchLinesPosition[j]] <- sample(pool, size=1)
      }
    }
    
    log = paste0("\\.log")
    matchLog <- grep(pattern = log, x = inFile, value = T)
    matchLogPosition <- grep(pattern = log, x = inFile, value = F)
    logRep <- paste0("\\.Rep", i, log)
    if (length(matchLogPosition) != 0) {
      newFile [matchLogPosition] <- gsub(log, logRep, matchLog)}
    
    trees = paste0("\\.trees")
    matchTrees <- grep(pattern = trees, x = inFile, value = T)
    matchTreesPosition <- grep(pattern = trees, x = inFile, value = F)
    treesRep <- paste0("\\.Rep", i, trees)
    if (length(matchTreesPosition) != 0) {
      newFile [matchTreesPosition] <- gsub(trees, treesRep, matchTrees)}
    
    csv = paste0("\\.csv")
    matchCsv <- grep(pattern = csv, x = inFile, value = T)
    matchCsvPosition <- grep(pattern = csv, x = inFile, value = F) 
    csvRep <- paste0("\\.Rep", i, csv)
    if (length(matchCsvPosition) != 0) {
      newFile [matchCsvPosition] <- gsub(csv, csvRep, matchCsv)}
    
    ops = paste0("\\.ops")
    matchOps <- grep(pattern = ops, x = inFile, value = T)
    matchOpsPosition <- grep(pattern = ops, x = inFile, value = F) 
    opsRep <- paste0("\\.Rep", i, ops)
    if (length(matchOpsPosition) != 0) {
      newFile [matchOpsPosition] <- gsub(ops, opsRep, matchOps)}
    
    if (writeTrees == F) {
      logA <- grep(pattern = "<logTree id=", x = newFile, value = T)
      logAn <- grep(pattern = "<logTree id=", x = newFile, value = F)
      newFile [logAn] <- paste0("\t\t<!-- \n", logA)
      logB <- grep(pattern = "</logTree>", x = newFile, value = T)
      logBn <- grep(pattern = "</logTree>", x = newFile, value = F)
      newFile [logBn] <- paste0(logB, "\n", " \t\t -->")
      }
    
    out <- paste0(name, ".Rep", i, ".xml")
    cat (newFile, file = out, sep = "\n")
  }
  
  cat ("Replicates done:", i,"\n")
}

if (ver == 2) {

  numberTaxa <- length(grep("taxon=", inFile))
  line <- grep(pattern = "traitname=\"date|traitname=\'date", x = inFile)
  line <- line + 1
  if (length(line) == 0) {stop(
    "No date info found, check Beast input file")}
  datePositions = c()
  
  repeat {
    if (length(grep("value=", inFile[line])) > 0) line <- line + 1
    if (length(grep("alignment", inFile[line])) > 0) break
    if (length(grep("=", inFile[line])) > 0)
      {datePositions <- c(datePositions, line)}
    line <- line + 1
  }
  
  numberDates <- length(datePositions)
  dateLines <- inFile[datePositions]
  dateLines <- trimws(dateLines)
  date <- unlist(strsplit(dateLines, "="))
  dateHap <- date[c(T, F)]
  dateHap <- dateHap[1: numberDates]
  dateValues <- date[c(F, T)]
  lastLine <- length(grep("<taxa", dateValues))
  
  if (lastLine == 1){
    lastDate <- tail(dateValues, 2)
    lastDate <- unlist(strsplit(lastDate, " "))
    lastDate <- head(lastDate, 1)
    dateValues <- head(dateValues, numberTaxa-1)
    dateValues <- c(dateValues, lastDate)
  }
  
  dateValues <- gsub(",$", "", dateValues)
  
  if (loadCluster == F){
    min <- 1
    d_clust <- Mclust(as.matrix(as.numeric(dateValues)), G = min:numberDates)
    n <- dim(d_clust$z)[2]
    if (n == 1) {stop(
      "Single cluster found, input clusters manually or use RandomDates")}
    cat("Number of clusters found:", n, "\n")
    clusterOut <- cbind(dateHap, d_clust$classification)
    nameCluster <- paste0(name, ".clusters.csv")
    write.table(clusterOut, file = nameCluster,
                row.names = F, col.names = F, sep = ",")
    } else {
    nameCluster <- paste0(name, ".clusters.csv")
    clusterOut <- read.table(nameCluster, sep = ",")
    }

  matchLinesCluster <- cbind(dateValues, clusterOut[,2])
  cluster <- split(matchLinesCluster[,1], matchLinesCluster[,2])
  n <- length(cluster)
  check <- lapply(cluster, length)
  clusterList <- matchLinesCluster[, 2]
  if (loadCluster == T){
    if (n == 1) {stop("Error, only One cluster found, check input")}
    cat("Number of clusters loaded:", n, "\n")
    }
  # Loop
  for(i in 1 :reps) {
    newFile <- inFile 
    for (j in 1 : numberTaxa){
      clusterL <- clusterList
      clusterId <- as.numeric(as.character(matchLinesCluster[j, 2]))
      if (clusterId != 0){
        clusterL <- clusterL[clusterL != clusterId]
        if ("0" %in% clusterL) {
          clusterL <- clusterL[clusterL != "0"]}
        if ("0" %in% clusterList && n == 2){
          clusterL <- clusterList[clusterList != "0"]}
        as.numeric(clusterL)
        pool <- matchLinesCluster[matchLinesCluster[, 2] %in% clusterL, 1]
        if (length(pool) == 1) {stop(
          "Single element in cluster, check input or use RandomDates")}
        dateValues[j] <- sample(pool, size = 1)
      }  else {
        dateValues[j] <- matchLinesCluster[j, 1]}
    }
    newDate <- paste0("\t\t\t", dateHap, "=", dateValues)
    newFile[datePositions] <- paste0(newDate, ",")
    if(lastLine == 1){newFile[(datePositions[numberDates])] <-
      paste0 (newDate[numberDates], "\t\t\t\t<taxa id=",
              date[numberDates* 2 + 1], "=", date[numberDates * 2 + 2])}
    
    log = paste0("\\.log")
    matchLog <- grep(pattern = log, x = inFile, value = T)
    matchLogPosition <- grep(pattern = log, x = inFile, value = F)
    logRep <- paste0("\\.Rep", i, log)
    if (length(matchLogPosition) != 0) {
      newFile [matchLogPosition] <- gsub(log, logRep, matchLog)}
    
    trees = paste0("\\.trees")
    matchTrees <- grep(pattern = trees, x = inFile, value = T)
    matchTreesPosition <- grep(pattern = trees, x = inFile, value = F)
    treesRep <- paste0("\\.Rep", i, trees)
    if (length(matchTreesPosition) != 0) {
      newFile [matchTreesPosition] <- gsub(trees, treesRep, matchTrees)}
    
    csv = paste0("\\.csv")
    matchCsv <- grep(pattern = csv, x = inFile, value = T)
    matchCsvPosition <- grep(pattern = csv, x = inFile, value = F) 
    csvRep <- paste0("\\.Rep", i, csv)
    if (length(matchCsvPosition) != 0) {
      newFile [matchCsvPosition] <- gsub(csv, csvRep, matchCsv)}

    if (writeTrees == F) {
      logA <- grep(pattern = "\\.trees", x = newFile, value = T)
      logAn <- grep(pattern = "\\.trees", x = newFile, value = F)
      newFile [logAn] <- paste0("\t<!-- \n ", logA)
      ctr <- 0
      
      repeat {
        ctr <- ctr + 1
        ctr2 <- 0
        ctr2 <- logAn + ctr
        temp <- grep(pattern = "</logger>", newFile[ctr2], value = F)
        temp <- length(temp)
        if (temp != 0) break
        if (ctr == 100) stop("Error, check files, no tree block found")
        }
      newFile [ctr2] <- paste0("\t</logger>", "\n", "\t-->")
    }

    out <- paste0(name, ".Rep", i, ".xml")
    cat (newFile, file = out, sep = "\n")
  }
  cat ("Replicates done:", i,"\n")
}

if (ver != 1 & ver != 2) {stop("Error, check BEAST input file")}
  
}