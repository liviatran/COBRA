#C1 v 0.1.0 Mar2023
#'Class I Pipeline Wrapper Function
#'
#'Generates peptide binding affinity scores for Class I and the reference SARS-CoV-2 strain.
#'Protein sequences for the reference SARS-COV-2 strain are bundled in the package. 
#'Users can also use COBRA for other virus sequences. 
#'
#'@param c1_predict_method Binding Affinity Prediction method used for Class I. Accepted values
#'are 'mhcflurry' or 'netmhcpan'. Default is set to 'netmhcpan'. 
#'@param fasta_files FASTA file for a given protein. Default 
#'value is set to the inst/extdata/ref_FASTA folder, which contains FASTA files
#'for the SARS-CoV2 reference proteins. If the user is using other FASTA files,
#'the full path to the file must be provided. FASTA files should only contain
#'ONE protein sequence per file if Class II evaluation is desired. Please read
#'the vignette for more details. 
#'@param chop_path The path to the user's installation of NetChop
#'@param mhcpan_path The path to the user's installation of NetMHCpan
#'@param dev A Boolean value for whether the function is being executed for development
#'and testing purposes. Default is set to FALSE. This value should NEVER be altered
#'by the user.
#'@param ds a BIGDAWG formatted HLA genotype dataset, where there is 
#'an identifying subject column (SID), a phenotype, where control = 0 and 
#'case = 1 (Disease), and two columns per HLA locus. Please see the vignette
#'for an example of a BIGDAWG formatted dataset.
#'@param c1_file Temp file for Class I pipeline log
#'
#'@importFrom lgr lgr AppenderFile
#'
#'@return  The full path for the Class I look up table. Otherwise, a fatal logging message is returned.
#'
#'@note The user's working directory MUST be set to where the input FASTA files are located if 
#'the bundled SARS-CoV2 FASTA files are not being used. 
#'@export

C1<-function(c1_predict_method, fasta_files, chop_path, mhcpan_path, dev, ds, c1_file){
  
  if(c1_file != "none"){
    lgr$add_appender(AppenderFile$new(c1_file), name = "c1_log")
  }
  
  lgr$info(paste(paste(rep("*", 20), collapse = ""), "CLASS I PROCESSING", paste(rep("*", 20), collapse = "")))
  lgr$info(paste("Class I Pipeline Temp Directory:", tempdir(), sep = ""))
  lgr$info(paste(c1_predict_method, "chosen as Class I Peptide Binding Affinity Predictor", sep = " "))
  
  #run choppR
  lgr$info("Executing choppR...")
  chop_catch <- choppR(fasta_files, chop_path)
  if (!is.null(chop_catch)) {
    return(lgr$fatal(chop_catch))
  }
  
  #generate potential peptides for Class I 
  lgr$info("Executing kmerizeI...")
  new_fasta <- kmerizeI(fasta_files, c1_predict_method)
  
  if(length(new_fasta) == 1){
    if(new_fasta == "No peptides"){
        return(lgr$fatal("There were no likely to be generated peptides from any of the FASTA files provided, based on NetChop and kmerizeI evaluation. Exiting Class I pipeline..."))
    }
  }
  
    #execute Class I Binding Affinity Predictor chosen
    if(c1_predict_method == "netmhcpan"){
      lgr$info("Executing NetMHCpan...")
      ########CHANGE THIS BACK BEFORE SUBMISSION#################
      mhcpan_catch <- run_netmhcpan(new_fasta, mhcpan_path, dev)
      if (!is.null(mhcpan_catch)) {
          return(lgr$fatal(mhcpan_catch))
      } else{
        lgr$info(paste("NetMHCPan results successfuly written to the temp directory for", paste(gsub(".faa", "", new_fasta), collapse = ", ")))
      }
    } else{
      lgr$info("Checking Python configurations...")
      cp_check <- configure_python()
      if (!cp_check %in% c("Previously configured virtual environment activated successfully!", "The venv's Python interpreter was successfully activated, and all modules were downloaded successfully!")){
        return(lgr$fatal(cp_check))
      } else{
        lgr$info("Executing mhcflurry...")
        mf_check <- run_mhcflurry(new_fasta, dev)
        if(mf_check != TRUE){
          return(lgr$fatal(mf_check))
        }
        lgr$info(paste("mhcflurry results successfully written to the temp directory for", paste(gsub(".faa", "", basename(new_fasta)), collapse = ", ")))
      }
    }
  
    #generate look up table for Class I
    lgr$info("Executing LUG...")
    lu_file<-LUG(ba_predict = c1_predict_method)
    
    return(lu_file)
}
  
