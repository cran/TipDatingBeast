ListTaxa <- function(name) {

options (warn = -1)
  
  inFileName <- paste0(name, ".xml")
  inFile <- readLines(inFileName)
  ver2 <- grep(pattern = "version=\"2.0\"", x = inFile, value = F)  
  ver1 <- grep(pattern = "version=\"1.0\"", x = inFile, value = F)  
  ver <- length(ver1) + length(ver2)
  
  if (ver == 1) {
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
    if (length(grep("=", inFile[line])) > 0) {
      datePositions <- c(datePositions, line)}
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

if (ver != 1 & ver != 2) {stop("Error, check BEAST input file")}
}
