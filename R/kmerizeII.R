#kmerizeII v 0.1.0 Mar 2023
#'Filters NetCleave results to likely to be generated epitopes for Class II HLA alleles
#'
#'Given a FASTA file, NetCleave generates all possible epitopes (13-17mers) for 
#'Class II HLA alleles and scores the likelihood of generation. Peptides with 
#'scores greater than the specified threshold of epitope generation are
#'filtered to be included in the output file, where the 'relaxed' threshold is > 0.5, and the
#''stringent' threshold is > 0.6. The filtered file is output
#'to the temp directory, named, named 'ClassII_Peptides_<file>_.txt', 
#'where <file> is the name of the FASTA file. 
#'
#'@param fasta FASTA file for a given protein. Default 
#'value is set to the inst/extdata/ref_FASTA folder, which contains FASTA files
#'for the SARS-CoV2 reference proteins. If the user is using other FASTA files,
#'the full path to the file must be provided. FASTA files should only contain
#'ONE protein sequence per file.
#'@param threshold the threshold to use for filtering NetCleave predictions.
#'A relaxed threshold filters peptides that are likely to be generated,
#'whereas a stringent threshold filters peptides that are highly likely to be generated.
#'Options are 'relaxed' or 'stringent'. Default is set to 'relaxed'.
#'
#'@importFrom dplyr %>% filter select
#'@importFrom lgr lgr
#'@importFrom utils write.table
#'
#'@return The names of the FASTA files with generated peptides. Otherwise, "No peptides" is returned.
#'
#'@note The user's working directory MUST be set to where the input FASTA files are located if 
#'the bundled SARS-CoV2 FASTA files are not being used. 
#'@export

kmerizeII<-function(fasta = list.files(system.file("extdata/ref_FASTA", package = "COBRA"), full.names = TRUE), threshold = "relaxed"){
  
  protfile <- gsub(".faa", "", basename(fasta))
  
  ncfile <- paste(tempdir(), "/", protfile, "_NetCleave.csv", sep="")
  
  netFiles<-sapply(protfile, function(x) NULL)
  
  lgr$info(paste(threshold, " threshold selected for filtering NetCleave predictions", sep = ""))
  
  
  for(a in 1:length(ncfile)){
    
    netFiles[[a]]<-read.csv(ncfile[[a]], stringsAsFactors = F, colClasses = c("NULL", NA, "NULL", NA, NA))
    
    if(!nrow(netFiles[[a]])==0){
      if(threshold == 'relaxed'){
        t <- 0.5
      } else{
        t <- 0.6
      }
      
      netFiles[[a]]<- netFiles[[a]] %>%
        filter(prediction > t) %>%
        select('epitope')
    }
  }
  netOut<-netFiles[lapply(netFiles, nrow) > 0]
  netNoPep<-names(netFiles[lapply(netFiles, nrow) == 0])
  
    if(length(netOut) != 0){
      for(i in 1:length(netOut)){
        write.table(netOut[[i]], file = paste(tempdir(), "/ClassII_Peptides_", gsub(".faa", "", basename(names(netOut)[[i]])), ".txt", sep=""), row.names = F, quote = F, col.names = F)
        lgr$info(paste("Class II generated peptides for", gsub("_NetCleave.csv", "", basename(names(netOut)[[i]])), "successfully written to the temp directory"))
      }
      if(length(netNoPep) != 0){
        lgr$info(paste("There were no likely to be generated peptides for", netNoPep, "based on NetCleave evaluation and threshold selection"))
      }
      return(gsub("_NetCleave.csv", "", basename(names(netOut))))
    } else{
      return("No peptides")
    }
}

