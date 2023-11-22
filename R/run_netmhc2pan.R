#runNetMHCIIpan v 0.3.0.9000 11/20/23
#'Run NetMHCpan on Class II HLA allele-peptide data
#'
#'Runs NetMHCIIpan to obtain peptide binding affinity for CIWD/CWD HLA Class II alleles and possible generated peptides.
#'CIWD files formatted specifically for NetMHCIIpan are bundled in the package. NetMHCpan results are output to the temporary
#'directory, and are named 'ClassI_Peptides_<file>_ClassICIWD_<n>_BindAff.txt', 
#'where <file> is the name of the FASTA file, and <n> is number for the corresponding 
#'CIWD file. 53 jobs are executed in parallel for faster computation.
#'
#'@param faaname FASTA file for a given protein or group of proteins. Default 
#'value is set to the inst/extdata/ref_FASTA folder, which contains FASTA files
#'for the SARS-CoV2 reference proteins 
#'@param netMHC2pan_path The path to the user's installation of NetMHCIIpan
#'@param dev A Boolean value for whether the function is being executed for development
#'and testing purposes. Default is set to FALSE. This value should NEVER be altered
#'by the user.
#'
#'@importFrom lgr lgr
#'@importFrom utils read.table
#'@importFrom tools file_path_sans_ext
#'
#'@return NULL. Otherwise, a string detailing the error that occurred. 
#'
#'@note The user's working directory MUST be set to where the input FASTA files are located if 
#'the bundled SARS-CoV2 FASTA files are not being used. 
#'@note NetMHCIIpan MUST be downloaded and the path to its executable MUST be set up properly in order
#'for this function to work properly. NetMHCIIpan v 4.1 was used at the time of 
#'COBRA development. In order to set up the full path to NetMHCIIpan,
#'open the executable in a text editor, and enter the path to your NetMHCIIpan download
#'where it says '# full path to the NetMHCIIpan 4.1 directory (mandatory)'. 
#'@note GNU Parallel MUST be installed in order for parallel computation to 
#'occur. It can be installed from https://www.gnu.org/software/parallel/, or via
#'homebrew with 'brew install parallel'.
#'
#'@source NetMHCIIpan 4.1, downloaded from https://services.healthtech.dtu.dk/services/NetMHCIIpan-4.0/
#'@source Reynisson B, Barra C, Kaabinejadian S, Hildebrand WH, Peters B, Nielsen M. Improved Prediction of MHC II Antigen Presentation through Integration and Motif Deconvolution of Mass Spectrometry MHC Eluted Ligand Data. J Proteome Res. 2020 Jun 5;19(6):2304-2315. doi: 10.1021/acs.jproteome.9b00874. Epub 2020 Apr 30. PMID: 32308001.
#'@export

run_netmhc2pan<-function(faaname = list.files(system.file("extdata/ref_FASTA", package = "COBRA"), full.names = TRUE), netMHC2pan_path, dev = FALSE){
  
  prot <- file_path_sans_ext(basename(faaname))
  
  all_pep  <- paste(tempdir(), "/ClassII_Peptides_", prot, ".txt", sep="")
  
  no_pep<-sapply(all_pep, function(x) NULL)
  
  #check if there are any peptide files with zero entries 
  for(z in 1:length(all_pep)){
    
    no_pep[[z]]<-readLines(all_pep[[z]])
  }
  
  #filter out any peptide files with no generated peptides
  peptidefiles<-paste(names(no_pep[lengths(no_pep) > 0]), collapse = " ")
  
  if(dev == FALSE){
    files <- list.files(system.file("extdata/CIWD_C2", package = "COBRA"), full.names = TRUE)
    ciwdfiles <- paste(files[!grepl("ClassIICIWD_test.txt", files)], collapse= " ")
    ciwd_1 <- list.files(system.file("extdata/CIWD_C2", package = "COBRA"), full.names = TRUE)[[1]]
  }else{
      ciwdfiles<-ciwd_1 <- "/Library/Frameworks/R.framework/Versions/4.2/Resources/library/COBRA/extdata/CIWD_C2/ClassIICIWD_test.txt"
    }

  lgr$info(paste("Running NetMHCIIpan in parallel for", paste(prot, collapse = ", ")))
  
  system(paste("parallel -j 100 '", netMHC2pan_path, " -inptype 1 -f `echo {1}` -BA -a `cat {2}`' '>'", tempdir(), "/`echo {1/.}`_`echo {2/.}`_II_BindAff.txt ::: ", peptidefiles, " ::: " , ciwdfiles, sep=""))

  #check one Generated Peptides file to see if netMHCpan was executed 
  netmhcIIpan_check<-tryCatch({read.table(paste(gsub(".txt", "", strsplit(peptidefiles, " ")[[1]][1]), "_", gsub(".txt", "", basename(ciwd_1)), "_II_BindAff.txt", sep = ""), skip = 11, fill = TRUE)}, error = function(error) {return("error")})
  
  if(!is.data.frame(netmhcIIpan_check)){
    return("The full path to the netMHCIIpan directory was not set up correctly. Please make sure the environment for netMHCIIpan is set correctly. This can be done by using a text editor to open the netMHCIIpan executable, and providing the full path for netMHCIIpan in the setenv argument.")
  } 
}
