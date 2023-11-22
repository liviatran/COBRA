#cleaveR v 0.3.0.9000 11/20/23
#'Runs NetCleave on input FASTA files
#'
#'Runs NetCleave on input FASTA files to predict likelihood of
#'lysosomal cleavage at the fourth amino acid position. The NetCleave output,
#'named '<file>_NetCleave.csv', where <file> is the name of the input FASTA file, is
#'output to the temp directory.
#'
#'@param fastafiles FASTA file for a given protein. Default 
#'value is set to the inst/extdata/ref_FASTA folder, which contains FASTA files
#'for the SARS-CoV2 reference proteins. If the user is using other FASTA files,
#'the full path to the file must be provided. FASTA files should only contain
#'ONE protein sequence per file.
#'@param netcleave_path The path to the cloned NetCleave repo.
#'
#'@importFrom lgr lgr
#'@importFrom utils capture.output
#'@importFrom seqinr read.fasta
#'@importFrom tools file_path_sans_ext
#'
#'@return TRUE. Otherwise, a string detailing the error that occurred. 
#'
#'@note NetCleave is not available via pip installation. Python scripts for NetCleave
#'were forked from https://github.com/BSC-CNS-EAPM/. Some scripts (NetCleave.py and
#'and cleavage_site_generator.py) were modified to accept an output path (--td) as 
#'an argument. The temp directory path is provided for the NetCleave predictions 
#'to be output to. Users MUST clone the Github repo with the modified 
#'scripts. The repo can be cloned by executing 'git clone https://github.com/liviatran/NetCleave.git'.
#'Link to the forked repo with modified scripts: https://github.com/liviatran/NetCleave
#'
#'@references Amengual-Rigo, P., Guallar, V. NetCleave: an open-source algorithm for predicting C-terminal antigen processing for MHC-I and MHC-II. Sci Rep 11, 13126 (2021). https://doi.org/10.1038/s41598-021-92632-y
#'@export

cleaveR<-function(fastafiles=list.files(system.file("extdata/ref_FASTA", package = "COBRA"), full.names = TRUE), netcleave_path){
  
  pwd <- getwd()
  
  on.exit(setwd(pwd))
  
  setwd(netcleave_path)
  
  fileBase <- file_path_sans_ext(basename(fastafiles))
    
  lgr$info(paste("Running NetCleave in parallel for", paste(fileBase, collapse = ", ")))
  
  for(i in 1:length(fastafiles)){
    f<-read.fasta(fastafiles[[i]])
    if(length(f[[1]])<16){
      lgr$info(paste(names(f[1]), "does not have a long enough peptide sequence to be evaluated properly by NetCleave. The minimum peptide length for NetCleave evaluation is 16 peptides. Final scores will not contain peptide binding affinity scores for this protein."))
    }
  }
  
  system(paste("parallel -j 100 'python3 NetCleave.py --predict {1} --pred_input 1 --mhc_class II --mhc_allele HLA --data_path ",  netcleave_path, "/data/training_data/II_mass-spectrometry_HLA --td ", tempdir(), "/' ::: ", paste(fastafiles, collapse = " "), sep=""), intern = T)
  
  netcleave_out <- paste(tempdir(), "/",  fileBase, "_NetCleave.csv", sep="")

  #if any anticipated NetCleave outputs are missing, stop function
  if(any(!file.exists(netcleave_out))){
    return("NetCleave failed. If non-bundled FASTA files are being used, please check that the provided FASTA files contain the full path, and only contain one protein sequence per file. 
           Please read the NetCleave output for further troubleshooting.")
  } else{
    lgr$info(paste("NetCleave successfully processed for ", paste(basename(fastafiles), collapse = " ,"), sep = ""))
  }
  return(TRUE)
}

