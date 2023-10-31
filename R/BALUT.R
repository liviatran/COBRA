#BALUT v 0.1.0 March2023
#'Creates lookup tables for HLA Class II allele-peptide pairs based on binding affinity strength
#'
#'Creates a lookup table for the total number of binders, strong binders, weak binders,
#'and weak and strong binders for Class II CIWD alleles based on binding affinity
#'results from run_netmhc2pan(). Strong Binders are categorized as 
#'peptide-allele pairs being less than 2.0 Rank, and Weak binders as less than or 
#'equal to 10.0 Rank, where Rank is the rank of the predicted binding score 
#'compared to a set of random natural peptides. Binding strength Rank value thresholds were
#'obtained from Reynisson et.al "NetMHCpan-4.1 and NetMHCIIpan-4.0: improved 
#'predictions of MHC antigen presentation by concurrent motif deconvolution and 
#'integration of MS MHC eluted ligand data" 2020, W449– W454m https://doi.org/10.1093/nar/gkaa379.
#'The final output is output to  the user's working directory, named 
#'LookUpTable_ClassII_<variant>_<threshold>_<random alphanumeric sequence>.csv', where <variant> is the 
#'SARS-CoV2 variant and <threshold> is the threshold selected. If the bundled SARS-CoV-2 FASTA files are used, the variant
#' will be 'ref'. 
#' 
#'@param wdc2 path to directory with Class II peptide binding affinity results. 
#'Default is set to the temp directory. 
#'@param threshold the threshold to use for filtering NetCleave predictions.
#'A relaxed threshold filters peptides that are likely to be generated,
#'whereas a stringent threshold filters peptides that are highly likely to be generated.
#'Options are 'relaxed' or 'stringent'. Default is set to 'relaxed'.
#'
#'@importFrom data.table rbindlist
#'@importFrom dplyr %>% distinct filter case_when
#'@importFrom lgr lgr
#'@importFrom utils read.csv write.csv read.table
#'
#'@return The full path for the Class II look up table
#'
#'@export

BALUT<- function(wdc2 = tempdir(), threshold = 'relaxed'){
  
  variant <- strsplit(list.files(path = wdc2, pattern="_II_BindAff")[1], "_")[[1]][4]
  
  #read in Allele, Peptide, Rank, and BA columns from netMHCpan output 
  BAfiles_list<-lapply(lapply(list.files(path = wdc2, pattern="_II_BindAff", full.names = T), function(x) read.table(x, skip = 11, fill = TRUE, row.names = NULL)), "[", c(2,3,13))
  
  BAfiles<-rbindlist(BAfiles_list)
  
  colnames(BAfiles) <- c("Allele", "Peptide", "Rank")
  
  BAfiles<-BAfiles %>%	
    distinct(Allele, Peptide, .keep_all = T)  %>%
    filter(grepl("DRB|DP|DQ", Allele))
  
  BAfiles<-BAfiles %>%
    mutate(Rank = as.numeric(BAfiles$Rank)) %>%
    mutate(BindStrength = case_when(Rank < 1 ~ "SB",
                                    Rank > 5 ~ "VWB")) 
  
  BAfiles$BindStrength[is.na(BAfiles$BindStrength)]<-"WB"
  
  BAfiles<-BAfiles %>%
    mutate(Allele = case_when(!grepl("HLA-", BAfiles$Allele) ~ paste("HLA-", BAfiles$Allele, sep = ""),
                              .default = Allele))
  
  BAfiles$Allele<-gsub("_(\\d{2})", '*\\1:', gsub("(D[PDQ]A1)(\\d{2})(\\d{2})-(D[PDQ]B1)(\\d{2})(\\d{2})", '\\1*\\2:\\3-\\4*\\5:\\6', BAfiles$Allele))
  
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
  
  out_file <- tempfile(paste("LookUpTable_ClassII_", variant, '_', threshold, "_", sep = ""), tmpdir = getwd(), fileext = '.csv')
  write.csv(BindingAffinityScores, file = out_file , row.names = F)
  
  lgr$info(paste("Class II Lookup Table sucessfully output to: ", out_file, sep = ""))
  
  return(out_file)
}