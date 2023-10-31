#kmerizeI v 0.1.0 Mar2023
#'Generate possible generated peptides during peptide processing for HLA Class I alleles
#'
#'Generates a list of estimated 8-12 residue long peptides based on potential
#'split site information obtained from NetChop for Class I HLA alleles. Possible
#'ending amino acids are indicated by potential split sites, and possible beginning 
#'amino acid positions are those directly after a split position. The list is 
#'output to the temp directory. If 'netMHCpan' is specified as the binding affinity
#'predictor, the output will be named 'ClassI_Peptides_<file>_.txt', 
#'where <file> is the name of the FASTA file. If 'mhcflurry' is selected, 
#'the output will be named '<file>_pep.txt'.
#'
#'@param FASTAfile FASTA file for a given protein or group of proteins. Default 
#'value is set to the inst/extdata/ref_FASTA folder, which contains FASTA files
#'for the SARS-CoV2 reference proteins 
#'@param ba_predict Binding Affinity Prediction method used for Class I. Accepted values
#'are 'mhcflurry' or 'netmhcpan'
#'
#'@importFrom dplyr %>%
#'@importFrom lgr lgr
#'@importFrom utils read.table tail write.table
#'
#'@return The names of the FASTA files with generated peptides. Otherwise, "No peptides" is returned.
#'
#'@note The user's working directory MUST be set to where the input FASTA files are located if 
#'the bundled SARS-CoV2 FASTA files are not being used. 
#'@export

kmerizeI<-function(FASTAfile = list.files(system.file("extdata/ref_FASTA", package = "COBRA"), full.names = TRUE), ba_predict){
  
  netchopFile <- paste(tempdir(), "/", gsub(".faa", "", basename(FASTAfile)), "_netchopfinal.txt", sep="")
  
  netchopOutput<-NULL
  
  peptides <- sapply(netchopFile, function(x) NULL)
  
  for(j in 1:length(netchopFile)){
    
    #read in netchop output, skip first row, header = F
    #cannot just use header = T, because then first column is not read in as a column,
    #but as row names
    netchopOutput[[j]]<-read.table(netchopFile[[j]], header = F, stringsAsFactors = F, skip=1)
    
    #add column names
    colnames(netchopOutput[[j]]) <- c("pos", "AA", "C", "score", "ident", "rownum")
    
    prot_split<-split(netchopOutput[[j]], cumsum(netchopOutput[[j]]$pos == 1))
    
    for(k in 1:length(prot_split)){
      
      #potential stop amino acids
      possEndSeq<-prot_split[[k]] %>%
        filter(C %in% "S")
      
      #possible start amino acids are those directly after a split position 
      possBeginSeq <- prot_split[[k]] %>%
        filter(rownum %in% c(possEndSeq$rownum+1))
      
      #find first residues for all proteins
      firstResidue<-prot_split[[k]] %>%
        filter(pos == 1)
      
      #use setdiff to find which first residues are not already present in possible
      #beginning sequences and rbind 
      possBeginSeq <- rbind(firstResidue, possBeginSeq)
      
      #for each row in possBeginSeq, examine if the tail end of 8-12mers is in 
      #possEndSeq
      #if it is, generate a peptide from the possBeginSeq to possEndSeq
      for(a in 1:nrow(possBeginSeq)){
        for(b in 8:12){
          if(tail(seq(possBeginSeq[a,]$pos, possBeginSeq[a,]$pos + b-1), 1) %in% possEndSeq$pos){
            peptides[[j]]<-append(peptides[[j]], paste(prot_split[[k]][possBeginSeq[a,]$pos:tail(seq(possBeginSeq[a,]$pos, possBeginSeq[a,]$pos + b-1), 1),]$AA, collapse = ""))
          }
        }
      }
    }
    peptides[[j]] <- data.frame(Peptides = peptides[[j]], stringsAsFactors = F)
  }
  
  peptidefiles <- peptides[lapply(peptides, nrow) > 0]
  no_peptides <- names(peptides[lapply(peptides, nrow) == 0])

  if (length(peptidefiles) != 0) {
    for (k in 1:length(peptidefiles)) {
      if (ba_predict == "mhcflurry") {
        write(noquote(paste(peptides[[k]]$Peptides, collapse = " ")),
              file = paste(tempdir(), "/", gsub("_netchopfinal.txt", "", basename(names(peptidefiles)[[k]])),"_pep.txt",sep = ""))
      }
      else{
        #write generated peptide file to user's working directory
        write(noquote(as.vector(peptides[[k]])$Peptides),file = paste(tempdir(),"/ClassI_Peptides_",gsub("_netchopfinal.txt", "", basename(names(peptidefiles)[[k]])),".txt",sep = ""))
      }
      lgr$info(paste("Class I generated peptides for", gsub("_netchopfinal.txt", "", basename(names(peptidefiles)[[k]])),"successfully written to the temp directory"))
    }
    if(length(no_peptides) != 0){
      lgr$info(paste("There were no likely to be generated peptides for", no_peptides, "based on NetChop evaluation and threshold selection"))
    }
    return(gsub("_netchopfinal.txt", "", basename(names(peptidefiles))))
  } else{
      return("No peptides")
  }
}
