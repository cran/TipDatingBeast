RandomDates <- function(name, reps = 20, writeTrees = T) {
 
  opar <-par(no.readonly=TRUE)
  on.exit(par(opar))
  
  inFileName <- paste0(name, '.xml')
  inFile <- readLines(inFileName)
  versionLine <- readLines(inFileName,n=1)
  ver <- grepl(pattern = 'version=\"1.0', x = versionLine)
  ver2 <- grepl(pattern = 'version=\"2.', x = versionLine)  
  ver1 = F
  if (ver == T & ver2 ==  F) {ver1 = T}

  if (ver1 == T) {
    
    matchLines <- grep(pattern = "date value=", x = inFile, value = T)
    if (length(matchLines) == 0) {stop(
      "No dates found, check BEAST input files")}  
    matchLinesPosition <- grep(pattern = "date value=", x = inFile, value = F)
    matchFileName <- grep(pattern = "fileName", x = inFile, value = T)
    matchFileNamePosition <- grep(pattern = "fileName", x = inFile, value = F)
    
    # Loop
    
    for (i in 1 : reps){
      
      newFile <- inFile
      randLines <- sample(matchLines)
      newFile[matchLinesPosition] <- randLines
      
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
  
  if (ver2 == T) {
    
    linearDates=F
    
    numberTaxa <- length(grep('taxon=', inFile))
    line <- grep(pattern = 'traitname=\"date|traitname=\'date', x = inFile)
    
    linearDates <- grepl(pattern = 'value=', inFile[line])
    
    if (linearDates == T) {
      
      dateLine <- inFile[line]
      step1 <- gsub('\">','',strsplit(dateLine, 'value=\"')[[1]][2])
      step2 <- unlist(strsplit(step1, ","))
      numberDates <- length(step2)
      
      date <- unlist(strsplit(step2, "="))
      dateHap <- date[c(T, F)]
      dateHap <- dateHap[1: numberDates]
      dateValues <- date[c(F, T)]
      
      dateValues <- gsub(",$", "", dateValues) 
      
      for(i in 1 :reps) {
        newFile <- inFile 
        dateValues <- sample(dateValues)
        
        newDate <- paste0(dateHap, "=", dateValues)
        
        newDateLine <- paste0(strsplit(dateLine, 'value=\"')[[1]][1],'value=\"',paste0(newDate,collapse=','),'\">') 
        
        newFile[line] <- newDateLine
        
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
    }
    if (linearDates == F) {
      
      line <- line + 1
      if (length(line) == 0) {stop(
        "No date info found, check BEAST input file")}  
      datePositions = c()
      
      repeat {
        if (length(grep("value=", inFile[line])) > 0) line <- line + 1
        if (length(grep("alignment", inFile[line])) > 0) break
        if (length(grep("=", inFile[line])) > 0) {datePositions <- c(datePositions, line)}
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
      
      # Loop
      
      for(i in 1 :reps) {
        newFile <- inFile 
        dateValues <- sample(dateValues)
        newDate <- paste0("\t\t\t", dateHap, "=", dateValues)
        newFile[datePositions] <- paste0(newDate, ",")
        datePositions[numberDates]
        if(lastLine == 1){newFile[(datePositions[numberDates])] <-
          paste0 (newDate[numberDates], "\t\t\t\t<taxa id=",
                  date[numberDates* 2 + 1],
                  "=", date[numberDates * 2 + 2])}
        
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
    }
    cat ("Replicates done:", i,"\n")
  }
  if (ver1 == F & ver2 == F) {stop("Error, check BEAST input file -version not recognized")}
}

