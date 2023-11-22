#C2 v 0.2.0.9000 10/31/23
#'Class II Pipeline Wrapper Function
#'
#'Generates peptide binding affinity scores for Class II and the reference SARS-CoV-2 strain.
#'Protein sequences for the reference SARS-COV-2 strain are bundled in the package. 
#'Users can also use COBRA for other virus sequences. 
#'
#'@param fasta_files FASTA file for a given protein. Default 
#'value is set to the inst/extdata/ref_FASTA folder, which contains FASTA files
#'for the SARS-CoV2 reference proteins. If the user is using other FASTA files,
#'the full path to the file must be provided. FASTA files should only contain
#'ONE protein sequence per file if Class II evaluation is desired. Please read
#'the vignette for more details. 
#'@param netmhcIIpanpath The path to the user's installation of NetMHCIIpan
#'@param netcleave_path The path to the cloned NetCleave repo.
#'@param threshold the threshold to use for filtering NetCleave predictions.
#'A relaxed threshold filters peptides that are likely to be generated,
#'whereas a stringent threshold filters peptides that are highly likely to be generated.
#'Options are 'relaxed' or 'stringent'. Default is set to 'relaxed'.
#'@param dev A Boolean value for whether the function is being executed for development
#'and testing purposes. Default is set to FALSE. This value should NEVER be altered
#'by the user.
#'@param ds_filename a BIGDAWG formatted HLA genotype dataset, where there is 
#'an identifying subject column (SID), a phenotype, where control = 0 and 
#'case = 1 (Disease), and two columns per HLA locus. Please see the vignette
#'for an example of a BIGDAWG formatted dataset.
#'@param c2_file Temp file for Class II pipeline log
#'
#'@importFrom lgr lgr AppenderFile
#'
#'@return The full path for the Class II look up table. Otherwise, a fatal logging message is returned.
#'
#'@note The user's working directory MUST be set to where the input FASTA files are located if 
#'the bundled SARS-CoV2 FASTA files are not being used. 
#'@export
#'

C2<-function(fasta_files, netcleave_path, threshold, netmhcIIpanpath, dev, ds_filename, c2_file){
  
  if(c2_file != "none"){
    lgr$add_appender(AppenderFile$new(c2_file), name = "c2_log")
  }
  
  lgr$info(paste(paste(rep("*", 20), collapse = ""), "CLASS II PROCESSING", paste(rep("*", 20), collapse = "")))
  lgr$info(paste("Class II Pipeline Temp Directory:", tempdir(), sep = ""))
  
  lgr$info("Checking Python configurations...")
  cp_check <- configure_python()
  if (!cp_check %in% c("Previously configured virtual environment activated successfully!","The venv's Python interpreter was successfully activated, and all modules were downloaded successfully!")){
    return(lgr$fatal(cp_check))
  }
  
  #run NetCleave
  lgr$info("Executing cleaveR...")
  netcleave_catch <- cleaveR(fasta_files, netcleave_path)
  if(netcleave_catch != TRUE){
    return(lgr$fatal(netcleave_catch))
  }
  
  #generate potential peptides for Class II
  lgr$info("Executing kmerizeII...")
  new_fasta <- kmerizeII(fasta_files, threshold)
  
  if(length(new_fasta) == 1){
    if(new_fasta == "No peptides"){
      return(lgr$fatal("There were no likely to be generated peptides from any of the FASTA files provided, based on NetCleave and kmerizeII evaluation. Exiting Class II pipeline..."))
    }
  }
  
  #run NetMHCIIpan
  lgr$info("Executing netMHCIIpan")
  mhcIIpan_catch <- run_netmhc2pan(new_fasta, netmhcIIpanpath, dev)
  
  if (!is.null(mhcIIpan_catch)) {
    return(lgr$fatal(mhcIIpan_catch))
  }
  
  #generate look up table for Class II
  lgr$info("Executing BALUT")
  lu_file<-BALUT(threshold = threshold)
  
  return(lu_file)
  
}
