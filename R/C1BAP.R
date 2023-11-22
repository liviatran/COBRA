#C1BAP v 0.2.0.9000 10/31/23
#'Class I Binding Affinity Predictor
#'
#'Compiles a comprehensive list of All Binders, Strong Binders, Weak Binders, and 
#'Weak and Strong Binders for the genotypes present in a BIGDAWG formatted dataset
#'for Class I HLA alleles.
#'
#'@param filename A BIGDAWG formatted dataset, which should have an identifying
#'subject column, a phenotype, where control = 0 and case = 1, and two columns
#'per Class I HLA locus. Please see the vignette for an example of a BIGDAWG
#'formatted dataset.
#'are 'mhcflurry' or 'netmhcpan'
#'@param lookup_file Unique name for the look up table generated from LUG()
#'
#'@return A dataframe containing peptide binding affinity predictions for Class I genotypes present in the provided BIGDAWG formatted dataset.
#'If the lookup table provided does not exist, '' is returned.
#'
#'@importFrom BIGDAWG GetField
#'@importFrom dplyr %>% filter
#'@importFrom utils read.table
#'
#'@export

C1BAP<-function(filename, lookup_file){
  
  if(!file.exists(lookup_file)){
    return('')
  }
  else{
    lookuptable<-read.csv(lookup_file, stringsAsFactors = F)
  }
  
  dataset<-read.table(filename, sep="\t", header=T, check.names = F, stringsAsFactors = F)
  
  #shorten allele resolution to 2 field
  for(i in 1:nrow(dataset)){
    for(k in 3:ncol(dataset)){
      dataset[i,k]<-GetField(dataset[i,k], 2)
    }
  }
  
  ScoreOutput<-data.frame(ID = dataset[,1], 
                          Disease = dataset[,2],
                          "AllBinders" = NA,
                          "StrongBinders" = NA, 
                          "WeakBinders" = NA,
                          "WeakandStrongBinders" = NA,
                          stringsAsFactors = F)
  #subset to class I loci
  dataset<-dataset[,colnames(dataset) %in% c("A", "B", "C")]
  colnames(dataset)<-gsub(".1", "", colnames(dataset), fixed = TRUE)
  
  #paste HLA and locus to allele
  for(i in 1:length(colnames(dataset))){
      dataset[,i]<-paste(paste("HLA-", colnames(dataset)[[i]], sep=""), dataset[,i], sep="*")
  }
  
  #Generate scores based on user input dataset
  for(z in 1:nrow(dataset)){
    temp<-lookuptable %>%
      filter(Alleles %in% c(as.vector(unlist(dataset[z,]))))
    
    #raw counts across all Class I alleles
    ScoreOutput[z,]$AllBinders<-sum(temp$AllBinders)
    ScoreOutput[z,]$StrongBinders<-sum(temp$StrongBinders)
    ScoreOutput[z,]$WeakBinders<-sum(temp$WeakBinders)
    ScoreOutput[z,]$WeakandStrongBinders<-sum(temp$WeakandStrongBinders)
  }
  
  return(ScoreOutput)
}
