#run_mhcflurryI v 0.1.0 Mar2023
#'Runs mhcflurry on allele-peptide data 
#'
#'Runs mhcflurry to obtain peptide binding affinity for CIWD HLA Class I alleles and possible generated peptides.
#'CIWD files formatted specifically for mhcflurry are bundled in the package. mhcflurry
#'will run in the COBRA venv if venv was set to TRUE in configure_python_mhcflurry(). 
#'mhcflurry results are output to the temporary directory, and are named 'ClassICIWD_<n>_<file>_pep_flur.csv',
#'where <n> is number for the corresponding CIWD file and <file> is the name of the FASTA file.
#'The jobs are run in parallel, with the optimal number of jobs determined based on the user's
#'machine's processing power.
#'
#'@param FASTAfile FASTA file for a given protein or group of proteins. Default 
#'value is set to the inst/extdata/ref_FASTA folder, which contains FASTA files
#'for the SARS-CoV2 reference proteins 
#'@param dev A Boolean value for whether the function is being executed for development
#'and testing purposes. Default is set to FALSE. This value should NEVER be altered
#'by the user.
#'
#'@importFrom reticulate import
#'@importFrom lgr lgr
#'
#'@return TRUE. Otherwise, a string detailing the error that occurred. 
#'
#'@note The user's working directory MUST be set to where the input FASTA files are located if 
#'the bundled SARS-CoV2 FASTA files are not being used. 
#'@note GNU Parallel MUST be installed in order for parallel computation to 
#'occur. It can be installed from https://www.gnu.org/software/parallel/, or via
#'homebrew with 'brew install parallel'.
#'
#'@source mhcflurry Github Repo: https://github.com/openvax/mhcflurry
#'@source T. J. O’Donnell, et al. “MHCflurry 2.0: Improved pan-allele prediction of MHC I-presented peptides by incorporating antigen processing,” Cell Systems, 2020. https://doi.org/10.1016/j.cels.2020.06.010
#'@export

run_mhcflurry<-function(FASTAfile = list.files(system.file("extdata/ref_FASTA", package = "COBRA"), full.names = TRUE), dev = FALSE){
  
  lib_check<- system("python3 -c 'import mhcflurry'", ignore.stdout  = T, ignore.stderr = T)
    
    if(lib_check == 1){
      return("There was an issue with importing the mhcflurry module. Please see the vignette for troubleshooting tips. Exiting job...")
    }

  lgr$info("Successfully imported mhcflurry!")
  
  if(dev == FALSE){
    files <- list.files(system.file("extdata/CIWD_I_mhcflurry", package = "COBRA"), full.names = TRUE)
    ciwdfiles <- paste(files[!grepl("ClassICIWD_test.txt", files)], collapse= " ")
  } else{
    ciwdfiles<-"/Library/Frameworks/R.framework/Versions/4.2/Resources/library/COBRA/extdata/CIWD_I_mhcflurry/ClassICIWD_test.txt"
  }
  
  prot <- gsub(".faa", "", basename(FASTAfile))
  
  all_pep  <- paste(tempdir(), "/", prot, "_pep.txt", sep="")
  
  no_pep<-sapply(all_pep, function(x) NULL)
  
  #check if there are any peptide files with zero entries 
  for(z in 1:length(all_pep)){
    no_pep[[z]]<-readLines(all_pep[[z]])
  }
  
  peptidefiles<-paste(names(no_pep), collapse = " ")
  
  lgr$info(paste("Running mhcflurry in parallel for", gsub(".faa", "", paste(basename(FASTAfile), collapse = ", "))))
  
  system(paste("parallel -j 100 'python3 /Library/Frameworks/R.framework/Versions/4.2/Resources/library/COBRA/python/run_mhcflurry.py -a {1} -p {2}' ::: ", ciwdfiles,  " ::: ", peptidefiles, sep = ""))

  return(TRUE)
}
