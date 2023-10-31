#C2BAP v 0.1.0 Mar2023
#'Class II Binding Affinity Predictor
#'
#'Compiles a comprehensive list of All Binders, Strong Binders, Weak Binders, and 
#'Weak and Strong Binders for the genotypes present in a BIGDAWG formatted dataset
#'for Class II HLA alleles. For DQ and DP, all possible heterodimer combinations
#'are generated from the BIGDAWG dataset, since DP and DQ haplotypes cannot be 
#'inferred. 
#'Forbidden DQ and DP heterodimers can be excluded during score generation. 
#'Forbidden DQ heterodimer rules are insert stuff here: <citation>
#'Forbidden DP heterodimer rules are based on peptide positions 31 and 50 (S.J. Mack, Personal Communication).
#'Users must select which DPA rule they would like implemented. 
#'Permissible heterodimers under the '31' rule include: DPA1 = Q, DPB1 = E OR if DPA1 = M, DPB1 = G.
#'Permissible heterodimers under the '50' rule include: DPA1 = R, DPB1 = E OR if DPA1 = Q, DPB1 = G. 
#'
#'@param filename A BIGDAWG formatted dataset, which should have an identifying
#'subject column, a phenotype, where control = 0 and case = 1, and two columns
#'per Class II HLA locus, except in the case of DRB3/4/5. Please see the vignette
#'for an example of a BIGDAWG formatted dataset.
#'@param excludeDQforbidden A Boolean value for whether forbidden DQ heterodimers should
#'be excluded from score evaluation. Default is set to TRUE. 
#'@param excludeDPforbidden A Boolean value for whether forbidden DP heterodimers should
#'be excluded from score evaluation. Default is set to TRUE. 
#'@param dpa_rule Dimerization rule to implement. Options include "31" or "50".
#'@param lookup_file Unique name for the look up table generated from BALUT()
#'
#'@importFrom BIGDAWG GetField
#'@importFrom dplyr %>% filter all_of
#'@importFrom SSHAARP BLAASD
#'@importFrom lgr lgr
#'@importFrom gtools mixedsort
#'@importFrom utils read.table
#'@importFrom stats na.omit
#'
#'@return A dataframe containing peptide binding affinity predictions for Class II genotypes present in the provided BIGDAWG formatted dataset
#''If the lookup table provided does not exist, '' is returned.
#'@references Position 31 rules: Hollenbach, J. et.al (2012). A combined DPA1~DPB1 amino acid epitope is the primary unit of selection on the HLA-DP heterodimer. Immunogenetics, 64(8), 559–569. https://doi.org/10.1007/s00251-012-0615-3
#'@references Position 50 rules:  Mack, S.J., Personal Communication
#'@export

C2BAP <- function(filename, excludeDQforbidden = TRUE, excludeDPforbidden = TRUE, dpa_rule = "50", lookup_file){

  if(!file.exists(lookup_file)){
    return("")
  } else{
    lookuptable<-read.csv(lookup_file, stringsAsFactors = F)
  }
  
  lgr$info(paste("Excluding forbidden DQ heterodimers?", excludeDQforbidden))
  lgr$info(paste("Excluding forbidden DP heterodimers?", excludeDPforbidden))
  
  if(excludeDPforbidden == TRUE){
    lgr$info(paste("Dimerization rule selected for excluding forbidden DP heterodimers:", dpa_rule))
  }
  
  dataset<-read.table(filename, sep="\t", header=T, check.names = F, stringsAsFactors = F, na.strings = c(" ", NA))
  
  subjs<-dataset[,c(1,2)]
  
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
  
  
  classIIloci<-colnames(dataset)[c(grepl("DRB|DP|DQ", colnames(dataset)))]
  
  #remove DQB1 and DPB1 if alpha subunits for DQ and DP are not present
  if("DQB1" %in% classIIloci & !"DQA1" %in% classIIloci){
    
    classIIloci<-classIIloci[!classIIloci %in% "DQB1"]
    
  }
  
  if("DPB1" %in% classIIloci & !"DPA1" %in% classIIloci){
    
    classIIloci<-classIIloci[!classIIloci %in% "DPB1"]
    
  }
  
  #subset dataset to only class II loci
  dataset<- dataset[, colnames(dataset) %in% classIIloci]
  
  #sort column names in alphabetical order, and bind with subject columns
  dataset<-cbind(subjs, dataset[mixedsort(colnames(dataset))]) 
  colnames(dataset)<-gsub(".1", "", colnames(dataset), fixed=TRUE)
  
  for (i in 3:ncol(dataset)){
    
    dataset[[i]]<- gsub(":", "", ifelse(is.na(dataset[,i])==FALSE, paste(colnames(dataset[i]),dataset[,i],sep="_"), NA))
    
  }
  
  #set up 'dimers' for for possible dimer combinations 
  dimers<-sapply(dataset[,1], function(x) list())
  dimers<-lapply(dimers, function(x) list(DP = NULL, DQ = NULL))
  
  for(d in 1:length(dimers)){
    
    dimers[[d]][["DRB1"]] <- c(t(dataset[d,colnames(dataset) == "DRB1"]))
    dimers[[d]][["DRB3"]] <- c(t(dataset[d,colnames(dataset) == "DRB3"]))
    dimers[[d]][["DRB4"]] <- c(t(dataset[d,colnames(dataset) == "DRB4"]))
    dimers[[d]][["DRB5"]] <- c(t(dataset[d,colnames(dataset) == "DRB5"]))
    
  }
  
  
  #if DP subunits are in column names, find all possible dimers 
  if(any(grep("DP", colnames(dataset)))){
    
    DPalleles<-dataset[,c(grepl("DP", colnames(dataset)))]
    colnames(DPalleles)<-gsub(".1", "", colnames(DPalleles), fixed = TRUE)
    
    for(i in 1:nrow(DPalleles)){
      for(j in 1:2){
        for(k in 3:4){
          if(excludeDPforbidden == TRUE){
            
            dpa<-BLAASD("DPA1")[[1]]
            dpb<-BLAASD("DPB1")[[1]]
            
            #ALLOWED - DPA1 position 31 = Q, DPB1 position 85 = E
            #DPA1 position 31 = M, DPB1 position 85 = G
            #DPA1 position 50 = R, DPB1 position 85 = E
            #DPA1 position 50 = Q, DPB1 position 85 = G
            
            dpa$trimmed_allele<-gsub(":", "", gsub("\\*", "_", dpa$trimmed_allele))
            dpb$trimmed_allele<-gsub(":", "", gsub("\\*", "_", dpb$trimmed_allele))
            
            dpa_allele<-(dpa %>%
                          filter(trimmed_allele == DPalleles[i,j]) %>%
                          select(all_of(dpa_rule)) %>%
                          distinct())[[1]]
          
            dpb_allele<-(dpb %>%
                          filter(trimmed_allele == DPalleles[i,k]) %>%
                          select(`85`) %>%
                          distinct())$`85` 
          
            if(dpa_rule == "31"){
              if((dpa_allele == "Q" & dpb_allele == "E") | (dpa_allele == "M" & dpb_allele == "G")){
                dimers[[i]][["DP"]]<-append(dimers[[i]][["DP"]], gsub("_", "", paste(DPalleles[i,j], DPalleles[i,k], sep = "-")))
              }
            } else{
              if((dpa_allele == "R" & dpb_allele == "E") | (dpa_allele == "Q" & dpb_allele == "G")){
                dimers[[i]][["DP"]]<-append(dimers[[i]][["DP"]], gsub("_", "", paste(DPalleles[i,j], DPalleles[i,k], sep = "-")))
              }
            }
          }  else{
              dimers[[i]][["DP"]]<-append(dimers[[i]][["DP"]], gsub("_", "", paste(DPalleles[i,j], DPalleles[i,k], sep = "-")))
        }
      }
    }
  }
}
  
  #if DQ subunits are in column names, find all possible dimers 
  if(any(grep("DQ", colnames(dataset)))){
    
    DQalleles<-dataset[,c(grepl("DQ", colnames(dataset)))]
    colnames(DQalleles)<-gsub(".1", "", colnames(DQalleles), fixed = TRUE)
    
    for(i in 1:nrow(DQalleles)){
      for(j in 1:2){
        for(k in 3:4){
          
          #if forbiddenDQ is true, assess if to be formed dimers meet the criteria 
          #for being considered forbidden; if false, append the dimer 
          if(excludeDQforbidden == TRUE){
            if((substr(DQalleles[i,j], 6, 7) == "01" & substr(DQalleles[i,k], 6, 7) %in% c("01", "05", "06")) |
               (substr(DQalleles[i,j], 6, 7) %in% c("02","03","04","05","06") & substr(DQalleles[i,k], 6, 7) %in% c("02","03","04"))){
              dimers[[i]][["DQ"]]<- append(dimers[[i]][["DQ"]], gsub("_", "", paste(DQalleles[i,j], DQalleles[i,k], sep = "-")))
            }
          }else{
            dimers[[i]][["DQ"]]<- append(dimers[[i]][["DQ"]], gsub("_", "", paste(DQalleles[i,j], DQalleles[i,k], sep = "-")))
            
          }
        }
      }
    }
  }
  
  #Generate scores based on user input dataset
  for(z in 1:length(dimers)){
    
    temp<-lookuptable %>%
      filter(Alleles %in% c(unique(paste('HLA-', na.omit(gsub("_(\\d{2})", '*\\1:', gsub("(D[PDQ]A1)(\\d{2})(\\d{2})-(D[PDQ]B1)(\\d{2})(\\d{2})", '\\1*\\2:\\3-\\4*\\5:\\6', as.vector(unlist(dimers[[z]]))))), sep = ""))))
    
    #raw counts across all Class II alleles
    ScoreOutput[z,]$AllBinders<-sum(temp$AllBinders)
    ScoreOutput[z,]$StrongBinders<-sum(temp$StrongBinders)
    ScoreOutput[z,]$WeakBinders<-sum(temp$WeakBinders)
    ScoreOutput[z,]$WeakandStrongBinders<-sum(temp$WeakandStrongBinders)
  }
  
  return(ScoreOutput)
}
