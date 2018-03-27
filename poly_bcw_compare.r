# This program compares BCW and Polysolver results and prints the results of Polysolver
# and any inconsistencies between the files.

# clears the console in RStudio
cat("\014")

# clears environment
rm(list = ls())

# Install packages if needed:
list.of.packages <-
  c("plyr", "stringi", "optparse")
new.packages <-
  list.of.packages[!(list.of.packages %in% installed.packages()[, "Package"])]
if (length(new.packages) > 0) {
  install.packages(new.packages, dependencies = TRUE, repos = 'http://cran.us.r-project.org')
}

# Required packages
suppressMessages(require(plyr))
suppressMessages(require(stringi))
suppressMessages(require(optparse))

# This passes the command line arguments to the variables
option_list = list(
  make_option(
    c("-s", "--winnersFile"),
    type = "character",
    default = NULL,
    help = "winners.hla.txt",
    metavar = "character"
  ),
  make_option(
    c("-p", "--bcwFile"),
    type = "character",
    default = NULL,
    help = "neon_study_output.csv",
    metavar = "character"
  )
)

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)


if (is.null(opt$winnersFile)) {
  print_help(opt_parser)
  stop("At least one argument must be supplied input file seqzFile",
       call. = FALSE)
}

if (is.null(opt$bcwFile)) {
  print_help(opt_parser)
  stop("At least one argument must be supplied for the BCW file",
       call. = FALSE)
}

# Name input files
winnersFile <- opt$winnersFile
bcwFile <- opt$bcwFile
# for troubleshooting
winnersFile <- "C:/Users/jconte/Documents/R/winners/example/example_winners.hla"
bcwFile <- "C:/Users/jconte/Documents/R/winners/example/example_output.csv"
setwd("C:/Users/jconte/Documents/R/winners/example/")

# Get project name
name <- tail(strsplit(winnersFile, split = '[/]')[[1]], 1)
projectName <- head(strsplit(name, split = '[.]')[[1]], 1)
projectName <- head(strsplit(name, split = '[_]')[[1]], 1)

# Load Neon Data
winnersVCF <-
  read.table(
    winnersFile,
    sep = "\t",
    stringsAsFactors = FALSE,
    row.names = NULL,
    header = FALSE
  )

# Load BCW Data
suppressWarnings (bcwVCF <-
                    read.csv(
                      bcwFile,
                      stringsAsFactors = FALSE,
                      row.names = NULL,
                      header = TRUE
                    ))

#put the winners in a format easy to compare
HLA_TYPE <- c(
  winnersVCF$V2[1],
  winnersVCF$V3[1],
  winnersVCF$V2[2],
  winnersVCF$V3[2],
  winnersVCF$V2[3],
  winnersVCF$V3[3]
)

HLA_TYPE <- gsub("hla_", "", HLA_TYPE)
HLA_TYPE <- sub("_", "*", HLA_TYPE)
HLA_TYPE <- gsub("_", ":", HLA_TYPE)
HLA_TYPE <- casefold(HLA_TYPE, upper = TRUE)

winnersdf <- data.frame(HLA_TYPE, stringsAsFactors = FALSE)

# # remove duplicates and name rows
winnersdf <-
  data.frame(unique(winnersdf[, 1]), stringsAsFactors = FALSE)
colnames(winnersdf) <- "HLA_TYPE"

# Length of winnersdf
winLength <- nrow(winnersdf)

# Get BCW file info
# Make sure columns are named properly
if ("HLA.TYPE" %in% names(bcwVCF) == TRUE) {
  bcwVCF$HLA_TYPE <- bcwVCF$HLA.TYPE
}

if ("CASE.ID" %in% names(bcwVCF) == TRUE) {
  bcwVCF$CASE_ID <- bcwVCF$CASE.ID
}

# make sure case id and file match
bcwVCF <- bcwVCF[with(bcwVCF,  grepl(projectName, CASE_ID)), ]

# make sure the HLA TYPE is an allele by checking if it has a '*'
bcwVCF <- bcwVCF[with(bcwVCF,  grepl("\\*", HLA_TYPE)), ]


# Reformat BCW results
HLA_TYPE <- bcwVCF$HLA_TYPE


# Loop removes letter from HLA calls from bcw
for (i in 1:length(HLA_TYPE))
{
  # function to check the last charachter
  substrRight <- function(x, n) {
    substr(x, nchar(x) - n + 1, nchar(x))
  }
  
  # Check if ends with a letter
  letter1 <- substrRight(HLA_TYPE[i], 1)
  
  # Check if it is a letter
  is.letter <- function(x)
    grepl("[[:alpha:]]", x)
  
  # if it is a letter remove it
  if (is.letter(letter1) == TRUE) {
    HLA_TYPE[i] <- substr(HLA_TYPE[i], 1, nchar(HLA_TYPE[i]) - 1)
    HLA_TYPE[i]
  }
  
}

# Get new dataframe with BCW results
bcwdf <- data.frame(HLA_TYPE, stringsAsFactors = FALSE)

# for troubleshooting
# bcwdf$HLA_TYPE[1]<-'A*24:03'
# bcwdf$HLA_TYPE[2]<-'A*24:03'

# # Remove duplicates and name column
bcwdf <- data.frame(unique(bcwdf[, 1]), stringsAsFactors = FALSE)
colnames(bcwdf) <- "HLA_TYPE"


# Length of BCW
bcwLength <- length(bcwdf$HLA_TYPE)



# Compare BCW and Poly and print response
# Initialize dataframe for different HLA aleles
dat <-
  data.frame(character(), character(), stringsAsFactors = FALSE)
hlaWin  <- character()

# Define variable
test <- character()

for (i in 1:bcwLength) {
  # See if they match
  test <- grepl(bcwdf$HLA_TYPE[i], winnersdf$HLA_TYPE[1:winLength])
  
  # if there are no matches, put it into a df
  if (any(test) == FALSE) {
    dat1 <-
      data.frame(bcwdf$HLA_TYPE[i], winnersdf$HLA_TYPE[i], stringsAsFactors = FALSE)
    colnames(dat1) <- c("bcwdf.HLA_TYPE", "winnersdf.HLA_TYPE")
    dat <- rbind(dat, dat1)
  }
  
  # This is a complicated way to check if BCW is in HLA and which is longest and adds it to a new vector,
  # if it does not exist, an NA is added to the vector.
  if (stri_detect_fixed(winnersdf$HLA_TYPE[1], bcwdf$HLA_TYPE[i]) == TRUE) {
    if (length(bcwdf$HLA_TYPE[i]) > winnersdf$HLA_TYPE[1]) {
      hlaWin[i] <- bcwdf$HLA_TYPE[i]
    }
    else {
      hlaWin[i] <- winnersdf$HLA_TYPE[1]
    }
  }
  
  if (stri_detect_fixed(winnersdf$HLA_TYPE[2], bcwdf$HLA_TYPE[i]) == TRUE) {
    if (length(bcwdf$HLA_TYPE[i]) > winnersdf$HLA_TYPE[2]) {
      hlaWin[i] <- bcwdf$HLA_TYPE[i]
    }
    else {
      hlaWin[i] <- winnersdf$HLA_TYPE[2]
    }
  }
  
  if ((winnersdf$HLA_TYPE[3] %in% winnersdf$HLA_TYPE) == TRUE) {
    if (stri_detect_fixed(winnersdf$HLA_TYPE[3], bcwdf$HLA_TYPE[i]) == TRUE) {
      if (length(bcwdf$HLA_TYPE[i]) > winnersdf$HLA_TYPE[3]) {
        hlaWin[i] <- bcwdf$HLA_TYPE[i]
      }
      else {
        hlaWin[i] <- winnersdf$HLA_TYPE[3]
      }
    }
  }
  if ((winnersdf$HLA_TYPE[4] %in% winnersdf$HLA_TYPE) == TRUE) {
    if (stri_detect_fixed(winnersdf$HLA_TYPE[4], bcwdf$HLA_TYPE[i]) == TRUE) {
      if (length(bcwdf$HLA_TYPE[i]) > winnersdf$HLA_TYPE[4]) {
        hlaWin[i] <- bcwdf$HLA_TYPE[i]
      }
      else {
        hlaWin[i] <- winnersdf$HLA_TYPE[4]
      }
    }
  }
}

# Make a new data frame
dattest <- dat

# This will replace the first na with the first value from dattest
if (any(is.na(hlaWin)) == TRUE) {
  for (i in 1:6) {
    if (is.na(hlaWin[i] == TRUE)) {
      hlaWin[i] <- dattest$bcwdf.HLA_TYPE[1]
      dattest <- dattest[-1,]
    }
  }
}

# Remove NA's
hlaWin <- hlaWin[!is.na(hlaWin)]

# this checks to make sure there are hla a alleles and hla b alleles.
# if it is missing one it adds it.
a <- grep("A", hlaWin, perl = TRUE, value = TRUE)
a_bcw <- grep("A", bcwdf$HLA_TYPE, perl = TRUE, value = TRUE)
if (length(a) != 0) {
  if (length(a) != length(a_bcw)) {
    if (grepl(a_bcw[1], a) != TRUE) {
      hlaWin <- c(hlaWin, a_bcw[1])
    }
    if (grepl(a_bcw[2], a) != TRUE) {
      hlaWin <- c(hlaWin, a_bcw[2])
    }
  }
  if (length(a) == 0) {
    hlaWin <- c(hlaWin, a_bcw)
  }
}


b <- grep("B", hlaWin, perl = TRUE, value = TRUE)
b_bcw <- grep("B", bcwdf$HLA_TYPE, perl = TRUE, value = TRUE)
if (length(b) != 0) {
  if (length(b) != length(b_bcw)) {
    if (grepl(b_bcw[1], b) != TRUE) {
      hlaWin <- c(hlaWin, b_bcw[1])
    }
    if (grepl(b_bcw[2], b) != TRUE) {
      hlaWin <- c(hlaWin, b_bcw[2])
    }
  }
}
if (length(b) == 0) {
  hlaWin <- c(hlaWin, b_bcw)
}



# If hla_c is not in hlaWin, it gets added here
test <- grepl("C", hlaWin)


if (any(test) != TRUE)
{
  c <- winnersdf[with(winnersdf,  grepl("C", HLA_TYPE)), ]
  hlaWin <- c(hlaWin, c)
}


# Make like original winners file
hlaWin <- paste0("hla_", hlaWin)
hlaWin <- sub("\\*", "_", hlaWin)
hlaWin <- gsub(":", "_", hlaWin)
hlaWin <- casefold(hlaWin, upper = FALSE)

# get all of the alleles
a <- grep("hla_a", hlaWin, perl = TRUE, value = TRUE)
b <- grep("hla_b", hlaWin, perl = TRUE, value = TRUE)
c <- grep("hla_c", hlaWin, perl = TRUE, value = TRUE)

# define V1
V1 <- winnersVCF$V1

# Define V2
V2 <- c(a[1], b[1], c[1])

# Define V3
if (a[2] %in% a == TRUE) {
  V3 <- a[2]
} else{
  V3 <- ""
}

if (b[2] %in% b == TRUE) {
  V3[2] <- b[2]
} else{
  V3[2] <- ""
}

if (c[2] %in% c == TRUE) {
  V3[3] <- c[2]
} else{
  V3[3] <- ""
}

# Get new data frame
winnersVCFnew <- data.frame(V1, V2, V3, stringsAsFactors = FALSE)

# name new file
winName <- gsub("[.][\\s\\S]*$", "", name, perl = T)


# Polysolver HLA_A and HLA_B alleles
polyHLAAB <-
  grep("^A+|^B+",
       winnersdf$HLA_TYPE,
       perl = true,
       value = TRUE)
lengthPolyHLAAB <- length(polyHLAAB)

# Logical if poly and bcw are the same length
lengthBP <- grepl(lengthPolyHLAAB, bcwLength)

## This section formats everything for printing
# Message function to print the information
message <- function(...)
  write(paste(...), file = stderr(), append = TRUE)

# Get current time
current_time <-
  paste(strsplit(as.character(Sys.time()), split = " +")[[1]], collapse =
          " ")

# If everything is the same
if (nrow(dat) == 0 && lengthBP == TRUE) {
  info <- c(
    projectName,
    "The HLA genotypes from Polysolver are:",
    winnersdf$HLA_TYPE,
    "\n",
    "BCW and Polysolver are consistent."
  )
  
  message(current_time, "\n")
  message(info)
  
  write.csv(
    info,
    file = paste0(projectName, ".consistent.BCW.Polysolver.csv"),
    row.names = FALSE
  )
  
  # If the alleles are inconsistent
} else if (nrow(dat) != 0 && lengthBP == TRUE) {
  info <- c(
    projectName,
    "The HLA genotypes from Polysolver are:",
    winnersdf$HLA_TYPE,
    "\n",
    "ALERT: BCW and Polysolver are inconsistent:",
    "The inconsistent HLA alleles from BCW are:",
    dat$bcwdf.HLA_TYPE,
    "The Polysolver file will be updated to match BCW"
  )
  message(current_time, "\n")
  message(info)
  
  write.csv(
    info,
    file = paste0(projectName, ".not.consistent.BCW.Polysolver.csv"),
    row.names = FALSE
  )
  
  # Write new info to hla file
  write.table(
    winnersVCFnew,
    file = paste0(winName, ".winners.hla.txt"),
    row.names = FALSE,
    sep = "\t",
    quote = FALSE,
    #row.names = FALSE,
    col.names = FALSE
  )
  
  # write old info to new file
  write.table(
    winnersVCF,
    file = paste0(winName, ".winners.hla.original.txt"),
    row.names = FALSE,
    sep = "\t",
    quote = FALSE,
    #row.names = FALSE,
    col.names = FALSE
  )
  
  # If the number of alleles do not match between Poly and BCW
} else if (nrow(dat) == 0 && lengthBP == FALSE) {
  info <- c(
    projectName,
    "The HLA genotypes from Polysolver are:",
    winnersdf$HLA_TYPE,
    "\n",
    "ALERT: BCW and Polysolver are inconsistent:",
    "The number of alleles do not match.",
    "Number of HLA A and B alleles for BCW:",
    bcwLength,
    "Number of HLA A and B alleles for Polysolver:",
    lengthPolyHLAAB,
    "the BCW results are:",
    bcwdf$HLA_TYPE,
    "The Polysolver file will be updated to match BCW"
  )
  message(current_time, "\n")
  message(info)
  
  write.csv(
    info,
    file = paste0(projectName, ".not.consistent.BCW.Polysolver.csv"),
    row.names = FALSE
  )
  
  # Write new info to hla file
  write.table(
    winnersVCFnew,
    file = paste0(winName, ".winners.hla.txt"),
    row.names = FALSE,
    sep = "\t",
    quote = FALSE,
    #row.names = FALSE,
    col.names = FALSE
  )
  
  # write old info to new file
  write.table(
    winnersVCF,
    file = paste0(winName, ".winners.hla.original.txt"),
    row.names = FALSE,
    sep = "\t",
    quote = FALSE,
    #row.names = FALSE,
    col.names = FALSE
  )
  # If the number of alleles do not match between Poly and BCW and the alleles are inconsistent
} else if (nrow(dat) != 0 && lengthBP == FALSE) {
  info <- c(
    projectName,
    "The HLA genotypes from Polysolver are:",
    winnersdf$HLA_TYPE,
    "\n",
    "ALERT: BCW and Polysolver are inconsistent:",
    "The inconsistent HLA alleles from BCW are:",
    dat$bcwdf.HLA_TYPE,
    "The number of alleles do not match.",
    "Number of HLA A and B alleles for BCW:",
    bcwLength,
    "Number of HLA A and B alleles for Polysolver:",
    lengthPolyHLAAB,
    "the BCW results are:",
    bcwdf$HLA_TYPE,
    "The Polysolver file will be updated to match BCW"
  )
  message(current_time, "\n")
  message(info)
  
  write.csv(
    info,
    file = paste0(projectName, ".not.consistent.BCW.Polysolver.csv"),
    row.names = FALSE
  )
  
  # Write new info to hla file
  write.table(
    winnersVCFnew,
    file = paste0(winName, ".winners.hla.txt"),
    row.names = FALSE,
    sep = "\t",
    quote = FALSE,
    #row.names = FALSE,
    col.names = FALSE
  )
  
  # write old info to new file
  write.table(
    winnersVCF,
    file = paste0(winName, ".winners.hla.original.txt"),
    row.names = FALSE,
    sep = "\t",
    quote = FALSE,
    #row.names = FALSE,
    col.names = FALSE
  )
}
