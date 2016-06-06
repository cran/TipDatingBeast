PlotLOOCV <- function(name, burnin = 0.1) {

options (warn = -1)

# library (TeachingDemos)

# read original file
  
inFileName <- paste0(name, ".xml")
inFile <- readLines(inFileName)
ver2 <- grep(pattern = "version=\"2.0\"", x = inFile, value = F)  
ver1 <- grep(pattern = "version=\"1.0\"", x = inFile, value = F)  
ver <- length(ver1) + length(ver2)
  
if (ver == 1) {
matchLines <- grep(pattern = "date value=", x = inFile, value = T)
values <- na.omit(as.numeric (gsub("[^0-9].", "", matchLines)))
taxa = length(values)
for (i in 1 : taxa){
  ## Read replicates log files
  fileRep <- paste0(name, ".Taxon", i, ".log")  
  temp1 <- read.table(fileRep, header = TRUE, sep = '\t')
  bn <- dim(temp1)[1] - round(dim(temp1)[1]*(burnin), 0)
  temp2 <- tail(temp1, bn)
  coltemp <- colnames(temp2[2])
  col <- unlist(strsplit(colnames(temp2[2]), split = "age.")) [2]
  
  if (i == 1) {data <- cbind(temp2[, 2])} else {
    data <- cbind(data, temp2[, 2])}
  if (i == 1) {colnames = col} else {
    colnames = c(colnames, col)}
  
    print (paste("log of taxon", i, "processed"))
}
}

if (ver == 2) {

  taxa <- length(grep("taxon=", inFile))
  line <- grep(pattern = "traitname=\"date|traitname=\'date", x = inFile)
  line <- line + 1
  if (length(line) == 0) {stop(
    "No date info found, check BEAST input file")}
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
  values <- date[c(F, T)]
  values <- na.omit(as.numeric (gsub("[^\\d]+", "", values, perl = T)))
  
  for (i in 1 : taxa){
    ## Read replicates log files
    fileRep <- paste0(name, ".Taxon", i, ".log")  
    temp1 <- read.table(fileRep, header = TRUE, sep = '\t')
    bn <- dim(temp1)[1] - round(dim(temp1)[1]*(burnin), 0)
    temp2 <- tail(temp1, bn)
    coltemp <- colnames(temp2[3])
    col <- unlist(strsplit(colnames(temp2[3]), split = "height.")) [2]
    
    if (i == 1) {data <- cbind(temp2[, 3])} else {
      data <- cbind(data, temp2[, 3])}
    if (i == 1) {colnames = col} else {
      colnames = c(colnames, col)}
    print (paste("log of taxon", i, "processed"))
  }
}
colnames(data) = colnames
### Convert in calendar years, last sampling date
if (ver == 1) {
data = ceiling(max(values)) - data
}

stats = matrix('NA',3,dim(data)[2])
colnames(stats) = colnames(data)
rownames(stats) = c("median", "min", "max")

for (i in 1:dim(stats)[2]){
  stats[1,i] = median(data[, i])  
  stats[2,i] = emp.hpd(data[, i], conf = 0.95)[1]
  stats[3,i] = emp.hpd(data[, i], conf = 0.95)[2]
}

### Export the data into a text table that summarizes the results
write.table (stats, paste0(name, "_leave_one_out.txt"))
mindata = floor(as.numeric(min(stats[2, ])))
maxdata = ceiling(as.numeric(max(stats[3, ])))
xsp <- maxdata - mindata
xlim <- c(mindata - xsp * 0.35, maxdata)
ylim <- c(-2*taxa, 0)

fail = NULL
plot(1, xlim = xlim, ylim = ylim, axes = F, xlab = "",
     ylab = "", type = "n", ann = F)
for (i in 1 : taxa){
  median = (as.numeric(stats[1, i]))
  min = round((as.numeric(stats[2, i])), 3)
  max = round((as.numeric(stats[3, i])), 3)
  label = substr(colnames[i], 1, 20)
  arrows (min, -2*i, max, -2*i, angle = 90, code = 3, length = 0.08, lwd = 1)
  points (median,  -2*i, cex = 1, pch = 3)
  if (min <= values[i] & values[i] <= max) {pt = 1; col = "black"}
  else {pt = 20; col = "red"; fail <- append(fail, i)}
  points (values[i],  -2*i, cex = 1, pch = pt, col = col, bg = col)
  text ((mindata - xsp * 0.20), -2*i, labels = label, cex = 0.5)
  } 

#Scale in bottom
axis(1, at = seq(mindata, maxdata), labels = seq(mindata, maxdata),
     cex.axis = 0.7, lwd = 1, line = 0)
LOOCV_result <- all(fail == 0)

# Print result over the plot
if (LOOCV_result == TRUE) {
  mtext("Pass!!! Age estimation for all taxa inside the expected 95% HPD",
        side = 3, line = 1, col = "red")
} else {mtext(paste("Age estimation for taxon/taxa",
                    paste(fail, collapse = ", "), 
                    "is/are not overlapping with expected 95% HPD"), side = 3,
					line = 2, col = "red");
  mtext("Attention !!! check LOOCV report file", side = 3, line = 1, 
        col = "red");
  write.table (colnames[fail], row.names = fail, sep = ",",
  col.names = "Taxon not overlapping with estimated 95% HPD (position, name)",
  paste0(name, "_LOOCV_report.txt"))
}

pdf(paste0("Fig_LOOCV_", name, ".pdf"))

plot(1, xlim = xlim, ylim = ylim, axes = F, xlab = "",
     ylab = "", type = "n", ann = F)
for (i in 1 : taxa){
  median = (as.numeric(stats[1, i]))
  min = round((as.numeric(stats[2, i])), 3)
  max = round((as.numeric(stats[3, i])), 3)
  label = substr(colnames[i], 1, 20)
  arrows (min, -2*i, max, -2*i, angle = 90, code = 3, length = 0.08, lwd = 1)
  points (median,  -2*i, cex = 1, pch = 3)
  if (min <= values[i] & values[i] <= max) {pt = 1; col = "black"}
  else {pt = 20; col = "red"; fail}
  points (values[i],  -2*i, cex = 1, pch = pt, col = col, bg = col)
  text ((mindata - xsp * 0.20), -2*i, labels = label, cex = 0.5)
} 
#Scale in bottom
axis(1, at = seq(mindata, maxdata), labels = seq(mindata, maxdata),
     cex.axis = 0.7, lwd = 1, line = 0)
LOOCV_result <- all(fail == 0)

# Print result over the plot
if (LOOCV_result == TRUE) {
  mtext("Pass!!! Age estimation for all taxa inside the expected 95% HPD",
        side = 3, line = 1, col = "red")
} else {mtext(paste("Age estimation for taxon/taxa",
                    paste(fail, collapse = ", "), 
                    "is/are not overlapping with expected 95% HPD"), side = 3,
					line = 2, col = "red");
  mtext("Attention !!! check LOOCV report file", side = 3, line = 1, 
        col = "red")
}

dev.off()
}
