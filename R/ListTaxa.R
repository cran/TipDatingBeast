ListTaxa <- function(name) {

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

    endTaxaLine <- grep(pattern = "</taxa>", x = inFile, value = F)
    taxaLine <- grep(pattern = "<taxon id=", x = inFile, value = T)
    taxaLinePosition <- grep(pattern = "<taxon id=", x = inFile, value = F)
    taxaLine <- unlist(strsplit(taxaLine, "\""))
    taxa <- taxaLine[c (F, T, F)]
    numberTaxa <- length(taxa)
    if (length(taxa) == 0) {stop(
      "No date info found, check Beast input file")}  
    matchFileName <- grep(pattern = "fileName", x = inFile, value = T)
    matchFileNamePosition <- grep(pattern = "fileName", x = inFile, value = F)
    
    
  # Loop
  for (i in 1 : numberTaxa) { 
    cat ("Taxon", i, "is", taxa[i], "\n")
  }
}

  if (ver2 == T) {
    
    linearDates=F
    
    numberTaxa <- length(grep('taxon=', inFile))

    line <- grep("sequence id=", inFile)
    line <- line[1]

  datePositions = c()
  
  repeat {
    if (length(grep("</data>", inFile[line])) > 0) break
    if (length(grep("taxon=", inFile[line])) > 0) {
      datePositions <- c(datePositions, line)}
    line <- line + 1
  }
  
  numberDates <- length(datePositions)
  dateLines <- inFile[datePositions]
  dateLines <- trimws(dateLines)
  date <- unlist(strsplit(dateLines, " totalcount"))
  dateHap <- date[c(T, F)]
  dateHap <- unlist(strsplit(dateHap, "taxon="))  
  dateHap <- dateHap[c(F, T)] 
  dateHap <- gsub("\"","",dateHap)
  
  dateHap <- dateHap[1: numberDates]
  dateValues <- date[c(F, T)]
  lastLine <- length(grep("<taxa", dateValues))
  
  dateValues <- gsub(",$", "", dateValues)
  
  lineTrees <- grep(pattern ="@Tree.t:", x = inFile)
  lineTree <- tail(lineTrees, 1) 
  treeLine <- inFile[lineTree]
  treePart <- tail(unlist(strsplit(treeLine, "@Tree.t:")), 1)
  treeName <- head(unlist(strsplit(treePart, "\"")), 1)
  # Loop
  for (i in 1 : numberTaxa) { 
    cat ("Taxon", i, "is", dateHap[i], "\n")
  }
}

  if (ver1 == F & ver2 == F) {stop("Error, check BEAST input file -version not recognized")}
}
