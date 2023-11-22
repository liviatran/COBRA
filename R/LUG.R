#LUG v 0.2.0.9000 10/31/23
#'Creates lookup tables for HLA Class I allele-peptide pairs based on binding affinity strength
#'
#'Creates a lookup table for the total number of binders, strong binders, weak binders,
#'and weak and strong binders for Class I CIWD alleles based on binding affinity
#'results from run_netmhcpan() or run_mhcflurry(). Strong Binders are categorized as 
#'peptide-allele pairs being less than 0.5 Rank, and Weak binders as less than or 
#'equal to 2.0 Rank, where Rank is the rank of the predicted binding score 
#'compared to a set of random natural peptides. Binding strength Rank value thresholds were
#'obtained from Reynisson et.al "NetMHCpan-4.1 and NetMHCIIpan-4.0: improved 
#'predictions of MHC antigen presentation by concurrent motif deconvolution and 
#'integration of MS MHC eluted ligand data" 2020, W449– W454m https://doi.org/10.1093/nar/gkaa379.
#'Only peptide-allele pairs with
#'at least loose binding peptide-allele pairs (less than 500 nm) will be evaluated. 
#'The final output is output to  the user's working directory, named 
#'LookUpTable_ClassI_<variant>_<predict-method>_<random alphanumeric sequence>.csv', where <variant> is the 
#'SARS-CoV2 variant, and <predict_method> is the binding affinity prediction 
#'method used. If the bundled SARS-CoV-2 FASTA files are used, the variant will be 'ref'. 
#'
#'@param wdc1 path to directory with Class I peptide binding affinity results. 
#'Default is set to the temp directory. 
#'@param ba_predict Binding Affinity Prediction method used for Class I. Accepted values
#'are 'mhcflurry' or 'netmhcpan'
#'
#'@importFrom dplyr %>% mutate distinct filter case_when
#'@importFrom data.table rbindlist fread
#'@importFrom lgr lgr
#'@importFrom utils read.table
#'
#'@return The full path for the Class I look up table
#'
#'
#'@export

LUG <- function(wdc1 = tempdir(), ba_predict){
  
  if(ba_predict == "mhcflurry"){
    variant <- strsplit(list.files(path = wdc1, pattern="_pep_flur")[1], "_")[[1]][4]
    BAfiles_list<-lapply(list.files(path = wdc1, pattern="_pep_flur", full.names = T), function(x) fread(x, select =c("peptide", "allele", "prediction_percentile")))    
    columns <-  c("Peptide", "Allele", "Rank")
  } else{
    variant <- strsplit(list.files(path = wdc1, pattern="_I_BindAff")[1], "_")[[1]][4]
    #read in Allele, Peptide, Rank, and BA columns from netMHCpan output 
    BAfiles_list<-lapply(lapply(list.files(path = wdc1, pattern="_I_BindAff", full.names = T), function(x) read.table(x, skip = 50, fill = TRUE)), "[", c(2,3,15))
    columns <-  c("Allele", "Peptide", "Rank")
  }
  
  BAfiles<-rbindlist(BAfiles_list)
  
  colnames(BAfiles) <- columns
  
  BAfiles<-BAfiles %>%	
    distinct(Allele, Peptide, .keep_all = T)
  
  if(ba_predict == "netmhcpan"){	
    BAfiles<-BAfiles %>%	
      filter(grepl("HLA", Allele))	
  }
  
  BAfiles$Rank<-as.numeric(BAfiles$Rank)
  
  BAfiles<-BAfiles %>%
    mutate(Rank = as.numeric(Rank)) %>%
      mutate(BindStrength = case_when(Rank < 0.5 ~ "SB",
                                      Rank > 2 ~ "VWB")) 
  
  BAfiles$BindStrength[is.na(BAfiles$BindStrength)]<-"WB"
  
  
  if(ba_predict == "mhcflurry"){	
    BAfiles$Allele <- gsub("(HLA-[A-Z])(\\d{2}:\\d{2})", '\\1*\\2', BAfiles$Allele)
  }
  
  BindingAffinityScores<-data.frame("Alleles" = sort(unique(BAfiles$Allele)),
                                    "AllBinders" = NA,
                                    "StrongBinders" = NA,
                                    "WeakBinders" = NA,
                                    "WeakandStrongBinders" = NA,
                                    stringsAsFactors = F)
  
  
  for(i in 1:nrow(BindingAffinityScores)){
    temp<-BAfiles %>%
      filter(Allele == BindingAffinityScores$Alleles[i])
    
    BindingAffinityScores[i,]$AllBinders<-nrow(temp)
    BindingAffinityScores[i,]$StrongBinders<-nrow(temp %>%
                                                    filter(BindStrength == "SB"))
    BindingAffinityScores[i,]$WeakBinders<-nrow(temp %>%
                                                  filter(BindStrength == "WB"))
    BindingAffinityScores[i,]$WeakandStrongBinders<-BindingAffinityScores[i,]$StrongBinders + BindingAffinityScores[i,]$WeakBinders
  }
  
  out_file <- tempfile(paste("LookUpTable_ClassI_", variant, "_", ba_predict, "_", sep = ""), tmpdir = getwd(), fileext = '.csv')
  
  write.csv(BindingAffinityScores, file = out_file, row.names = F)
  
  lgr$info(paste("Class I Lookup Table sucessfully output to: ", out_file, sep = ""))
  
  return(out_file)
}

