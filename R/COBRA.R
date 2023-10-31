#COBRA v 0.1.0 Mar2023
#'COVID-19 Binding Affinity Repertoire Assessment
#'
#'Generates peptide binding affinity scores for Class I and/or Class II HLA alleles
#'and the reference SARS-CoV-2 strain. Protein sequences for the reference SARS-COV-2
#'strain are bundled in the package. Users can also use COBRA for other virus
#'sequences.
#'
#'@param generatefiles a Boolean value for whether the lookup table needs to be generated. If TRUE,
#'the entire COBRA pipeline will be executed, and all parameters must be provided. If FALSE, 
#'only ds_filename, res_sel, and lookuptable_c1 and/or lookuptable_c2 (depending on res_sel chosen)
#'parameters need to be populated. Default value is FALSE. Note: we do not recommend generating
#'your own look up tables if they are already bundled in the package. Running the entire pipeline
#'for the SARS-CoV-2 sequence is computationally expensive. Please see more details in the vignette.
#'@param ds_filename a BIGDAWG formatted HLA genotype dataset, where there is 
#'an identifying subject column (SID), a phenotype, where control = 0 and 
#'case = 1 (Disease), and two columns per HLA locus. Please see the vignette
#'for an example of a BIGDAWG formatted dataset. 
#'@param res_sel pipeline selection for whether results should contain HLA Class I
#'evaluation, HLA Class II evaluation, or both. Options include "c1", "c2", or "both". 
#'@param c1_predict_method Binding Affinity Prediction method used for Class I. Accepted values
#'are 'mhcflurry' or 'netmhcpan'. Default is set to 'netmhcpan'. 
#'@param fasta_files FASTA file for a given protein. Default 
#'value is set to the inst/extdata/ref_FASTA folder, which contains FASTA files
#'for the SARS-CoV2 reference proteins. If the user is using other FASTA files,
#'the full path to the file must be provided. FASTA files should only contain
#'ONE protein sequence per file if Class II evaluation is desired. The FASTA file
#'should be follow this naming convention: '<protein>_<identifier>.faa', where 
#'<protein> is the protein name for the given sequence, and <identifier> is the 
#'the viral strain or organism identifier. EX: 'S_ref.faa', where 'S' is the Spike
#'protein, and 'ref' is the SARS-CoV2 reference strain. If the user wishes to regenerate
#'look up tables for the Beta or Omicron BA.1 sequences, 
#'list.files(system.file("extdata/<folder>", package = "COBRA"), full.names = TRUE) 
#'should be used, where the name of the <folder> is beta_FASTA or omiba1_FASTA. 
#'the vignette for more details. 
#'@param chop_path The path to the user's installation of NetChop. Required for Class I processing.
#'@param mhcpan_path The path to the user's installation of NetMHCpan. Required for Class I processing.
#'@param netmhcIIpanpath The path to the user's installation of NetMHCIIpan. Required for Class II processing.
#'@param exclude_forbidDQ A Boolean value for whether forbidden DQ heterodimers should
#'be excluded from score evaluation. Default is set to TRUE.
#'@param exclude_forbidDP A Boolean value for whether forbidden DP heterodimers should
#'be excluded from score evaluation. Default is set to TRUE. 
#'@param dpa_rule Dimerization rule to implement if excludeforbidDP is TRUE. 
#'Options include "31" or "50". Default is set to "50".
#'@param nc_path The path to the cloned NetCleave repo. Required for Class II processing.
#'@param threshold the threshold to use for filtering NetCleave predictions.
#'A relaxed threshold filters peptides that are likely to be generated,
#'whereas a stringent threshold filters peptides that are highly likely to be generated.
#'Options are 'relaxed' or 'stringent'. Default is set to 'relaxed'.
#'@param dev A Boolean value for whether the function is being executed for development
#'and testing purposes. Default is set to FALSE. This value should NEVER be altered
#'by the user.
#'@param lookuptable_c1 The path to the look up table generated for Class I. This parameter ONLY needs to be 
#'populated if generatefiles = FALSE and res_sel = 'c1' or 'both'. If generatefiles = TRUE, the generated look up table
#'name will be passed to the downstream processing functions. If generatefiles = FALSE, the
#'path to the look up table name must be provided.
#'@param lookuptable_c2 The path to the look up table generated for Class II. This parameter ONLY needs to be 
#'populated if generatefiles = FALSE and res_sel = 'c2' or 'both'. If generatefiles = TRUE, the generated look up table
#'name will be passed to the downstream processing functions. If generatefiles = FALSE, the
#'path to the look up table name must be provided. Class II peptide-allele mapping to the
#'look up table will use defualt values for exclude_forbidDP, exclude_forbidDQ, and dpa_rule, unless otherwise specified.
#'
#'
#'@importFrom lgr add_appender AppenderFile remove_appender
#'@importFrom future %<-% plan multisession
#'@importFrom dplyr group_by %>% summarise
#'@importFrom seqinr read.fasta
#'
#'@note The user's working directory MUST be set to where the input FASTA files are located if 
#'the bundled SARS-CoV2 FASTA files are not being used. 
#'@note If an error occurs in the pipeline AFTER the look up table has been generated,
#'generatefiles should be set to FALSE.
#'@note If generatefiles = TRUE, COBRA will take a long time to run. If Class I 
#'or Class II pipelines are chosen, the log will output to the console to update the
#'user on which step the execution is on. If both pipelines are chosen, the log
#'in the console will not update as the job executes, since they are separate processes. 
#'To view what step each pipeline is on, copy the 'Class I Log File' and 'Class II Log File' file
#'paths from the beginning of the log console. In the Terminal, run 'more <copied_path>', where
#'copied_path is the path to the Class I or Class II log. 
#'@note FASTA files for the Beta and Omicron BA.1 SARS-CoV2 Variants of Concern (VOCs) are included
#'in COBRA. These sequences were scraped from user submissions to the NCBI website and are not officially
#'published by NCBI, so they may not be  entirely accurate. The other VOCs were not included, 
#'due to an issue with formatting or completeness with the user submissions. 
#'Look up tables for the SARS-CoV-2 reference, Beta, and Omicron BA.1 proteomes are bundled in COBRA.
#'Link to Beta NCBI submission: https://www.ncbi.nlm.nih.gov/nuccore/MW598419.1
#'Link to Omicron BA.1 NCBI submission: https://www.ncbi.nlm.nih.gov/nuccore/OL672836
#'
#'@examples
#'#COBRA usage for generating look up tables for ONLY Class I HLA alleles
#'\donttest{COBRA(generatefiles = TRUE, res_sel = "c1", c1_predict_method = "netmhcpan", chop_path = "pathtoyournetchop/netchop-3.1/netchop", mhcpan_path = "pathtoyour/netMHCpan-4.1/netMHCpan", ds_filename = "your_bigdawg_dataset.txt")}
#'#COBRA usage for generating look up tables for ONLY Class II HLA alleles
#'\donttest{COBRA(generatefiles = TRUE, res_sel = "c2", nc_path = 'pathtoyour/NetCleave', threshold = "relaxed", netmhcIIpanpath = "pathtoyour/netMHCIIpan-4.1/netMHCIIpan", ds_filename = "your_bigdawg_dataset.txt")}
#'#COBRA usage for generating look up tables for Class I AND Class II HLA alleles
#'\donttest{COBRA(generatefiles = TRUE, res_sel = "both", c1_predict_method = "mhcflurry", exclude_forbidDQ = TRUE, exclude_forbidDP = TRUE, dpa_rule = "50", chop_path = "pathtoyour/netchop-3.1/netchop", mhcpan_path = "pathtoyour/netMHCpan-4.1/netMHCpan", netmhcIIpanpath = "pathtoyour/netMHCIIpan-4.1/netMHCIIpan",ds_filename = "your_bigdawg_dataset.txt", threshold = "relaxed", nc_path = 'pathtoyour/NetCleave')}
#'#COBRA usage for generating scores for an existing look up table for class I
#'\donttest{COBRA(generatefiles = FALSE, res_sel = 'c1', ds_filename = "yourbigdawg_dataset.txt")}
#'
#'@return A look up table/s for the pipeline selection, a COBRA log, and the scores generated for the input dataset,
#'output to the user's working directory. If 'both' was selected, a log for Class I and Class II will be output to 
#'the user's working directory. The COBRA log will not contain Class I and Class II processing logging information, 
#'since they are running in parallel.
#'The Class I look up table will be named LookUpTable_ClassI_<variant>_<predict-method>_<random alphanumeric sequence>.csv', 
#'where <variant> is the SARS-CoV2 variant, and <predict_method> is the binding affinity prediction method used. 
#'If the bundled SARS-CoV-2 FASTA files are used, the variant will be 'ref'. 
#'The Class II look up table will be named 'LookUpTable_ClassII_<variant>_<threshold>_<random alphanumeric sequence>.csv', 
#'where <variant> is the SARS-CoV2 variant and <threshold> is the threshold selected. If the bundled SARS-CoV-2 FASTA files are used, the variant will be 'ref'. 
#'The COBRA log will be named 'COBRA_<random alphanumeric sequence>.log'. 
#'The Class and Class II logs will be named 'c1_<random alphanumeric sequence>.log' and 'c2_<random alphanumeric sequence>.log', respectively. 
#'Generated score files will be named '<dataset>_ScoreOutput_<res>_<random alphanumeric sequence>', where <dataset> is the name of the dataset and <res> is the pipeline selection. 
#'
#'@export

COBRA <- function(generatefiles = FALSE, ds_filename = NULL, res_sel = NULL, lookuptable_c1 = NULL, lookuptable_c2 = NULL, c1_predict_method = "netmhcpan", fasta_files = list.files(system.file("extdata/ref_FASTA", package = "COBRA"), full.names = TRUE), exclude_forbidDQ = TRUE, exclude_forbidDP = TRUE, dpa_rule = "50", threshold = "relaxed", dev = FALSE, chop_path = NULL, mhcpan_path = NULL, netmhcIIpanpath = NULL, nc_path = NULL){
  
  log_file <- tempfile(tmpdir = getwd(), pattern = "COBRA_", fileext = ".log")
  lgr$add_appender(AppenderFile$new(log_file), name = "COBRA_log")
  lgr$info(paste("COBRA Log File: ", log_file, sep = ""))
  
  #parameter checks
  if(is.null(res_sel) | !res_sel %in% c('c1','c2','both')){
    return(lgr$fatal("No input was found for the res_sel parameter, or an invalid input was provided. Please provide a valid input for the res_sel parameter. Valid inputs include 'c1' for Class I, 'c2' for Class II, or 'both' for Class I and Class II."))
  }
  
  if(is.null(ds_filename)){
    return(lgr$fatal("No input was found for the ds_filename parameter. Please enter the full path to a BIGDAWG formatted dataset."))
  }
  
  if(!file.exists(ds_filename)){
    return(lgr$fatal("The input for ds_filename does not exist. Please make sure the input file exists."))
  }
  
  if(exclude_forbidDP == TRUE){
    if(!dpa_rule %in% c('31', '50')){
      return(lgr$fatal("DP rule entered for excluding forbidden DP heterodimers is not valid. Only '31' and '50' are valid."))
    }
  }
  
  lgr$info(paste("Entering COBRA pipeline... '", res_sel, "' selected as pipeline execution", sep = ""))

  if(generatefiles == TRUE){
    
    lgr$info(paste("FASTA files provided as input:", paste(basename(fasta_files), collapse = ", ")))
    
    lgr$info(paste("Temporary directory path:", tempdir()))
    
    lgr$info(paste(paste(rep("*", 20), collapse = ""), "PARAMETER CHECKS", paste(rep("*", 20), collapse = "")))
    
    #check FASTA files
    for(k in 1:length(fasta_files)){
      
      #nonexistent FASTA file
      if(!file.exists(fasta_files[[k]])){
        return(lgr$fatal(paste(fasta_files[[k]], "does not exist in the provided directory. Please make sure input FASTA files are spelled correctly, and that they are present in the working directory.")))
      }
      
      #improper FASTA formatting
      if(length(read.fasta(fasta_files[[k]])[[1]])==0){
        return(lgr$fatal(paste("No FASTA entries were detected in ", fasta_files[[k]], ". Please make sure the FASTA file is formatted properly and has a peptide sequence.", sep = "")))
      }
      
      if('>' %in% read.fasta(fasta_files[[k]])[[1]]){
        return(lgr$fatal(paste("There was a '>' detected in the peptide sequence for ", fasta_files[[k]], ". Please make sure the starting amino acid starts on a new line after the identifier in the FASTA file.", sep = "")))
      }
      
      if(' ' %in% read.fasta(fasta_files[[k]])[[1]]){
        return(lgr$fatal(paste("Unnecessary white space was detected in ", fasta_files[[k]], ". Please make sure there is no extraneous white space or empty lines in the FASTA file.", sep = "")))
      }
    }
    
    #c1 parameter checks
    if(res_sel == "c1" | res_sel == "both"){
      #check if prediction method for Class I has been specified if res is c1 or both 
      if(!c1_predict_method %in% c("netmhcpan", "mhcflurry")){
        return(lgr$fatal("The Class I peptide binding affinity prediction method is either missing, or is not supported. Only 'netmhcpan' and 'mhcflurry' are currently supported."))
      }
      if(is.null(chop_path)){
        return(lgr$fatal("A NetChop path was not provided. Please make sure to provide a path your NetChop executable. Please see choppR's documentation on how to download NetChop."))
      }
      if(!file.exists(chop_path)){
        return(lgr$fatal("The path provided for NetChop is not valid. Please make sure the path to NetChop is correct."))
      }
      if(is.null(mhcpan_path)){
        return(lgr$fatal("A NetMHCpan path was not provided. Please make sure to provide a path your NetMHCpan executable. Please see run_netmhcpan's documentation on how to download NetMHCpan."))
      }
      if(!file.exists(mhcpan_path)){
        return(lgr$fatal("The path provided for NetMHCpan is not valid. Please make sure the path to NetMHCpan is correct."))
      }
    }
    
    #c2 parameter checks
    if(res_sel == "c2" | res_sel == "both"){
      if(!dir.exists(nc_path)){
        return(lgr$fatal("The NetCleave path provided does not exist. Please double check the path and make sure it is correct."))
      }
      if(!file.exists(paste(nc_path, "NetCleave.py", sep = "/"))){
        return(lgr$fatal("NetCleave.py was not detected in the provided path. Please make sure the NetCleave repository was cloned properly. See the vignette for more details."))
      }
      if(is.null(nc_path)){
        return(lgr$fatal("A NetCleave path was not provided. Please provide a path to your local clone of the NetCleave repository. Please see cleaveR's documentation on how to clone the NetCleave repository."))
      }
      if(is.null(netmhcIIpanpath)){
        return(lgr$fatal("A NetMHCIIpan path was not provided. Please make sure to provide a path your NetMHCIIpan executable. Please see run_netmhc2pan's documentation on how to download NetMHCIIpan."))
      }
      if(!file.exists(netmhcIIpanpath)){
        return(lgr$fatal("The path provided for NetMHCIIpan is not valid. Please make sure the path to NetMHCIIpan is correct."))
      }
      if(!threshold %in% c("relaxed", "stringent")){
        return(lgr$fatal("The threshold entered is invalid. Please make sure the threshold is 'relaxed' or 'stringent'."))
      }
    }
    
    lgr$info(paste(paste(rep("*", 20), collapse = ""), "PARAMETER CHECKS COMPLETE", paste(rep("*", 20), collapse = "")))
    
    if(res_sel == "both"){
      out_dir = getwd()
    } else{
      out_dir = tempdir()
    }
    
    #CLASS I ONLY
    if(res_sel == "c1"){
      out1<-C1(c1_predict_method, fasta_files, chop_path, mhcpan_path, dev, ds_filename, "none")
      
      if(!grepl("LookUpTable_ClassI_", out1)){
        return(lgr$fatal("There was an error in generating the lookup table for Class I. Exiting job..."))
      }
    }
    
    #CLASS II ONLY
    if(res_sel == "c2"){
      out2<-C2(fasta_files, nc_path, threshold, netmhcIIpanpath, dev, ds_filename, "none")
      
      if(!grepl("LookUpTable_ClassII_", out2)){
        return(lgr$fatal("There was an error in generating the lookup table for Class II. Exiting job..."))
      }
    }
    
    #CLASS I AND CLASS II
    if(res_sel == "both") {
      
      plan(multisession)
      
      options(future.rng.onMisuse="ignore")
      
      c1 <- tempfile(pattern = "c1_", fileext = ".log", tmpdir = out_dir)
      lgr$info(paste("Class I Log File: ", c1, sep = ""))
      
      runC1 %<-% {
        C1(c1_predict_method, fasta_files, chop_path, mhcpan_path, dev, ds_filename, c1)
      }
      
      c2 <- tempfile(pattern = "c2_", fileext = ".log", tmpdir = out_dir)
      lgr$info(paste("Class II Log File: ", c2, sep = ""))
      
      runC2  %<-% {
        C2(fasta_files, nc_path, threshold, netmhcIIpanpath, dev, ds_filename, c2)
      }
      
      out1<-runC1
      out2<-runC2
    }
  } 
  
  cap1<-cap2<-NULL

  lgr$info(paste(paste(rep("*", 20), collapse = ""), "SCORE GENERATION", paste(rep("*", 20), collapse = "")))
  
  #generate binding affinity predictions for Class I
  if (res_sel == "c1" | res_sel == "both") {
    if (generatefiles == FALSE) {
      if (is.null(lookuptable_c1)) {
        return(lgr$fatal("A path to the look up table for Class I was not provided. Please provide a valid file path to the Class I look up table."))
      }
        if (file.exists(lookuptable_c1)) {
          out1 <- lookuptable_c1
        } else{
          return(lgr$fatal("The look up table provided for Class I does not exist. Please provide a valid file path to the Class I look up table."))
        }
      }
    lgr$info("Executing C1BAP...")
    cap1 <- C1BAP(ds_filename, out1)
    if(!is.data.frame(cap1)){
      lgr$info("There was a fatal error that prevented the Class I Pipeline from fully processing. Class I results will not be included in the final score output.")
      cap1 <- NULL
    } else{
    lgr$info("Class I Binding Affinity Predictions for dataset successfully generated")
    }
  }
  
  #generate binding affinity predictions for Class II
  if(res_sel == "c2" | res_sel == "both"){
    if(generatefiles == FALSE){
      if (is.null(lookuptable_c2)) {
        return(lgr$fatal("A path to the look up table for Class II was not provided. Please provide a valid file path to the Class II look up table."))
      }
    if(file.exists(lookuptable_c2)){
      out2<-lookuptable_c2
    } else{
      return(lgr$fatal("The look up table provided for Class II does not exist. Please provide a valid file path to the Class II look up table."))
    }
  }
    lgr$info("Executing C2BAP...")
    cap2 <- C2BAP(ds_filename, exclude_forbidDQ, exclude_forbidDP, dpa_rule, out2)
    if(!is.data.frame(cap2)){
      lgr$info("There was a fatal error that prevented the Class II Pipeline from fully processing. Class II results will not be included in the final score output.")
      cap2 <- NULL
    } else{
    lgr$info("Class II Binding Affinity Predictions for dataset successfully generated")
    }
  }
  
  if(res_sel == "both"){
    if(!is.data.frame(cap1) & !is.data.frame(cap2)){
      return(lgr$fatal("Both Class I and Class II Pipelines did not run to completion due to an error. Please check the log to troubleshoot where the error occurred."))
    }
  }

    final_out<-rbind(cap1[,1:6], cap2[,1:6])
      
    final_out<- final_out %>%
      group_by(ID, Disease) %>% 
      summarise(AllBinders = sum(AllBinders),
                StrongBinders = sum(StrongBinders),
                WeakBinders = sum(WeakBinders), 
                WeakandStrongBinders = sum(WeakandStrongBinders), .groups = "keep")
    
    out_file = tempfile(paste(gsub(".txt", "", basename(ds_filename)), "_ScoreOutput_", res_sel, "_", sep = ""), tmpdir = getwd(), fileext = ".csv")
    
    write.csv(final_out, file = out_file, row.names = F)
    
    lgr$info(paste("Score output for ", basename(ds_filename), " has been written to ", out_file, sep = ""))
    
    remove_appender("COBRA_log")
    
}
