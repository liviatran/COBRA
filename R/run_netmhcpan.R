#runNetMHCpan v 0.3.0.9000 11/20/23
#'Run NetMHCpan on Class I HLA allele-peptide data
#'
#'Runs NetMHCpan to obtain peptide binding affinity for CIWD HLA Class I alleles and possible generated peptides.
#'CIWD files formatted specifically for NetMHCpan are bundled in the package. NetMHCpan results are output to the temporary
#'directory, and are named 'ClassI_Peptides_<file>_ClassICIWD_<n>_BindAff.txt', 
#'where <file> is the name of the FASTA file, and <n> is number for the corresponding 
#'CIWD file. 46 jobs are executed in parallel for faster computation.
#'
#'@param fastafilename FASTA file for a given protein or group of proteins. Default 
#'value is set to the inst/extdata/ref_FASTA folder, which contains FASTA files
#'for the SARS-CoV2 reference proteins 
#'@param netmhcpan_path The path to the user's installation of NetMHCpan
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
#'@note NetMHCpan MUST be downloaded and the path to its executable MUST be set up properly in order
#'for this function to work properly. NetMHCpan 4.0 was used at the time of COBRA
#'development. In order to set up the full path to NetMHCpan, open the executable
#'in a text editor, and enter the path to your NetMHCpan download where it says 
#'# full path to the NetMHCpan 4.0 directory (mandatory)'. 
#'@note GNU Parallel MUST be installed in order for parallel computation to 
#'occur. It can be installed from https://www.gnu.org/software/parallel/, or via
#'homebrew with 'brew install parallel'.
#'
#'@source NetMHCpan 4.1, downloaded from https://services.healthtech.dtu.dk/services/NetMHCpan-4.1/
#'@source Reynisson B, Alvarez B, Paul S, Peters B, Nielsen M. NetMHCpan-4.1 and NetMHCIIpan-4.0: improved predictions of MHC antigen presentation by concurrent motif deconvolution and integration of MS MHC eluted ligand data, Nucleic Acids Research, Volume 48, Issue W1, 02 July 2020, Pages W449– W454, https://doi.org/10.1093/nar/gkaa379
#'@export

run_netmhcpan<-function(fastafilename = list.files(system.file("extdata/ref_FASTA", package = "COBRA"), full.names = TRUE), netmhcpan_path, dev = FALSE){
  
  fileBase <- file_path_sans_ext(basename(fastafilename))
  
  all_pep  <- paste(tempdir(), "/ClassI_Peptides_", fileBase, ".txt", sep="")
  
  no_pep<-sapply(all_pep, function(x) NULL)
  
  #check if there are any peptide files with zero entries 
  for(z in 1:length(all_pep)){
    
    no_pep[[z]]<-readLines(all_pep[[z]])
  }
  
  peptidefiles<-paste(names(no_pep), collapse = " ")
  
  if(dev == FALSE){
    files <- list.files(system.file("extdata/CIWD_C1_netmhcpan", package = "COBRA"), full.names = TRUE)
    ciwdfiles <- paste(files[!grepl("ClassICIWD_test.txt", files)], collapse= " ")
    ciwd_1 <- paste(list.files(system.file("extdata/CIWD_C1_netmhcpan", package = "COBRA"), full.names = TRUE)[[1]], collapse= " ")
  }else{
    ciwdfiles<-ciwd_1 <- "/Library/Frameworks/R.framework/Versions/4.2/Resources/library/COBRA/extdata/CIWD_C1_netmhcpan/ClassICIWD_test.txt"
  }
  
  lgr$info(paste("Running NetMHCPan in parallel for", paste(fileBase, collapse = ", ")))
  
  system(paste("parallel -j 100 '", netmhcpan_path, " -p `echo {1}` -BA -a `cat {2}` -l 8,9,10,11,12' '>'", tempdir(), "/`echo {1/.}`_`echo {2/.}`_I_BindAff.txt ::: ", peptidefiles, " ::: " , ciwdfiles, sep=""), intern = TRUE)

  #check one Generated Peptides file to see if netMHCpan was executed 
  netmhcpan_check<-tryCatch({read.table(paste(gsub(".txt", "", strsplit(peptidefiles, " ")[[1]][1]), "_", gsub(".txt", "", basename(ciwd_1)), "_I_BindAff.txt", sep = ""), skip = 50, fill = TRUE)}, error = function(error) {return("error")})
  
  #if netmhcpan_check is not a dataframe ('error' was returned), then netMHCpan
  #was not set up correctly
  if(!is.data.frame(netmhcpan_check)){
    return("The full path to the netMHCpan directory was not set up correctly. Please make sure the environment for netMHCpan is set correctly. This can be done by using a text editor to open the netMHCpan executable, and providing the full path for netMHCpan in the setenv argument.")
  }
}
