# Create Index file for demultiplexing SHARE-seq runs with deML

library(yaml)
args <- commandArgs(trailingOnly = TRUE)

#stop the script if no command line argument
if(length(args)==0){
  print("Please include the location of the .yaml file.")
  stop("Requires command line argument.")
}

print("Argument provided, started creating the index file.")

# 1. Fist provide all Round1 indices:
R1_rev.comp <- c("ATCACGTT", "CGATGTTT", "TTAGGCAT", "TGACCACT", "ACAGTGGT", "GCCAATGT", "CAGATCTG", "ACTTGATG", "GATCAGCG", "TAGCTTGT",
                 "GGCTACAG", "CTTGTACT", "TGGTTGTT", "TCTCGGTT", "TAAGCGTT", "TCCGTCTT", "TGTACCTT", "TTCTGTGT", "TCTGCTGT", "TTGGAGGT",
                 "TCGAGCGT", "TGATACGT", "TGCATAGT", "TTGACTCT", "TGCGATCT", "TTCCTGCT", "TAGTGACT", "TACAGGAT", "TCCTCAAT", "TGTGGTTG",
                 "TACTAGTC", "TTCCATTG", "TCGAAGTG", "TAACGCTG", "TTGGTATG", "TGAACTGG", "TACTTCGG", "TCTCACGG", "TCAGGAGG", "TAAGTTCG",
                 "TCCAGTCG", "TGTATGCG", "TCATTGAG", "TGGCTCAG", "TATGCCAG", "TCAGATTC", "TAGTCTTG", "TTCAGCTC", "TGTCTATC", "TATGTGGC",
                 "TTACTCGC", "TCGTTAGC", "TACCGAGC", "TGTTCTCC", "TTCGCACC", "TTGCGTAC", "TCTACGAC", "TGACAGAC", "TAGAACAC", "TCATCCTA",
                 "TGCTGATA", "TAGACGGA", "TGTGAAGA", "TCTCTTCA", "TTGTTCCA", "TGAAGCCA", "TACCACCA", "TGCGTGAA", "GGTGAGTT", "GATCTCTT",
                 "GTGTCCTT", "GACGGATT", "GCAACATT", "GGTCGTGT", "GAATCTGT", "GTACATCT", "GAGGTGCT", "GCATGGCT", "GTTAGCCT", "GTCGCTAT",
                 "GGAATGAT", "GAGCCAAT", "GCTCCTTG", "GTAAGGTG", "GAGGATGG", "GTTGTCGG", "GGATTAGG", "GATAGAGG", "GTGTGTCG", "GCAATCCG",
                 "GACCTTAG", "GCCTGTTC", "GCACTGTC", "GCTAACTC", "GATTCATC", "GTCTTGGC")
## add the constant regions
R1_reverse <- c()
for (i in 1:length(R1_rev.comp)) {
  R1_reverse[i] <- paste0("TCGGACGATCATGGG", R1_rev.comp[i], "CAAGTATGCAGCGCGCTCAAGCACGTGGAT")
}

# 2. Second, provide all Ad1.X (P5) indices
P5_reverse <- c("GCGATCTA", "ATAGAGAG", "AGAGGATA", "TCTACTCT", "CTCCTTAC", "TATGCAGT", "TACTCCTT", "AGGCTTAG", "GATTTCCA", "ATCATGTT",
                "TTTCATCA", "AGTCCGAC", "GCTAGAAA", "CTTGGTTA", "CGATACAC", "TTGATGGA", "TGCACGAA", "GGCAACCT", "ACATAAGG", "CGTTGCTG",
                "ATTGAACC", "ACGAATGT", "TGGGAATC", "GCAGTCCG", "GAACGGCT", "GACCCAAT", "AGTATGCA", "CCAAGCCC", "GCCACGTC", "AAATTTGC",
                "GAGGCTGC", "AACTCGGA", "CTTAATGC", "GTTATCGT", "CCCGCAGG", "AACAATCA", "TCCGTGCC", "GAATGATC", "ATGACCAT", "TTGGTACG",
                "TAAACTGG", "GGGCCGGT", "ACTTCTAG", "ATCTGGCG", "CCATGTGA", "TCGAGTTC", "AACGGTGG", "GTAACTTA", "CACGTCTC", "TTAGGCAA",
                "CAAGTTAA", "TGTTAAAG", "GGTCTACG", "CGCAAATA", "TCCTGGAT", "CAGGAACA", "CTGCGCGT", "TCGCCAGA", "TGTAGATT", "GGTCAGTA",
                "CCCTATCG", "TTCTAAGT", "AGATCTCT", "CCTTCACC", "CATTCGAT", "GCTCTTGA", "ACGTGGGC", "ACCGCCCA", "TCCAAGGG", "ACGGTAAT",
                "CTCGGACT", "CAACAAGT", "TGTATTAC", "TAGACGCC", "AGCAGCGC", "AATGGCAC", "CATACCTA", "TAGGTGTT", "GTTCGGAG", "TGCCGTTG",
                "CTACATTG", "GGGTAGCC", "CGGACTTT", "CCGCGGAA", "AAGTGCCT", "CACTGAAG", "CTACCGGC", "GGATTGAA", "GTGTGTGG", "GATAATAT",
                "TGCTTCGG", "ACCGATAC")

##################################
library(yaml)
# Provide a .yaml file containing the project names, used P5s and the first and last well for Round1 barcoding plates

## Read in data ----
yaml<-read_yaml(args[1])

# Check the structure of the yaml file
print(paste("The number of projects is", length(yaml), "."))
# yaml[[1]]$Name
# yaml[[1]]$Primer
# yaml[[1]]$Round1

# Get the necessary indices for the Round 1 barcodes for each sample

fun.extract.number <- function(x) {
  # Split the string by dot
  parts <- strsplit(x, "\\.")[[1]]
  # Extract the part after the dot and convert to numeric
  num <- as.numeric(parts[2])
  # If the number is less than 10, convert it to a single digit
  if (num < 10) {
    num <- as.numeric(substr(parts[2], 2, 2))
  }
  return(num)
}

# Return a list with all Round1 barcodes
Index1 <- list()
for (i in 1:length(yaml)) {
  Round1_indices <- sapply(yaml[[i]]$Round1, fun.extract.number)
  if (length(yaml[[i]]$Primer) == 1) {
    Index1[[i]] <- R1_reverse[Round1_indices[1]:Round1_indices[2]]
    } else {
    Index1[[i]] <- rep(R1_reverse[Round1_indices[1]:Round1_indices[2]], length(yaml[[i]]$Primer))
  }
}

# Return a list with all P5 (reverse chemistry)
Index2 <- list()
for (i in 1:length(yaml)){
  Round1_indices <- sapply(yaml[[i]]$Round1, fun.extract.number)
  P5_indices <- sapply(yaml[[i]]$Primer, fun.extract.number)
  if (length(yaml[[i]]$Primer) == 1) {
    Index2[[i]] <- rep(P5_reverse[P5_indices], length(R1_reverse[Round1_indices[1]:Round1_indices[2]]))
  } else {
    if (length(yaml[[i]]$Primer) == 2) {
      x <- rep(P5_reverse[P5_indices[1]], length(R1_reverse[Round1_indices[1]:Round1_indices[2]]))
      y <- rep(P5_reverse[P5_indices[2]], length(R1_reverse[Round1_indices[1]:Round1_indices[2]]))
      Index2[[i]] <- c(x, y)
    } else {
      if (length(yaml[[i]]$Primer) == 3) {
        x <- rep(P5_reverse[P5_indices[1]], length(R1_reverse[Round1_indices[1]:Round1_indices[2]]))
        y <- rep(P5_reverse[P5_indices[2]], length(R1_reverse[Round1_indices[1]:Round1_indices[2]]))
        z <- rep(P5_reverse[P5_indices[3]], length(R1_reverse[Round1_indices[1]:Round1_indices[2]]))
        Index2[[i]] <- c(x, y, z)
      } else {
        print("No more than 3 P5 indices supported per project.")
      }
    }
  }
}

Name <- list()
for (i in 1:length(yaml)){
  Name[[i]] <- rep(yaml[[i]]$Name, length(Index1[[i]]))
}

# now write this to a dataframe
deml.df <- data.frame(unlist(Index1), unlist(Index2), unlist(Name))
names(deml.df) <- c("#Index1", "Index2", "Name")
# head(deml.df)

# write the index file as .txt
write.table(deml.df, file = "index_SHARE.txt",
            sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

print("Successfully written index_SHARE.txt to working directory.")
