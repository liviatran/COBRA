#configure_python v 0.2.0.9000 10/31/23
#'Configure system with required modules to run Python code, either via a virtual environment or the user's machine
#'
#'A properly configured Python environment is required to run COBRA for Class II HLA alleles, 
#'or if the user wishes to use mhcflurry as the peptide binding affinity predictor. 
#'A virtual environment named "COBRA" will be created, and modules required to run
#'mhcflurry, or the Class II pipeline, will be installed to the virtual environment. Once successful, the project's
#'Python interpreter will be set to the virtual environment's Python interpreter.
#'Troubleshooting advice for virtual environment activation and installation can be found in the vignette.
#'
#'@note Python 3 must be installed on the user's system.
#'@note For internal use only.
#
#'@importFrom reticulate py_discover_config use_virtualenv virtualenv_list virtualenv_create virtualenv_root py_install
#'@importFrom lgr lgr
#'
#'@return ""Previously configured virtual environment activated successfully!" or a string detailing the error that occurred.
#'
#'@export

configure_python <- function(){
  
    if(!"COBRA" %in% virtualenv_list()){
      
      lgr$info(paste("Creating virtual environment named 'COBRA' and installing required modules to run mhcflurry or the Class II pipeline. The venv can be found in ", virtualenv_root(), ".", sep =""))
      
      virtualenv_create(envname = "COBRA")
      py_install(c("pandas", "argparse", "numpy", "matplotlib", "pathlib", "scikit-learn", "keras", "tensorflow", "mhcflurry", "biopython", "silence_tensorflow"), envname = "COBRA")
      use_virtualenv("COBRA")
      check<-tryCatch({a<-import('mhcflurry')}, error=function(cond){return("error")})
      
      if(check != "error"){
        lgr$info("The venv's Python interpreter was successfully activated, and all modules were downloaded successfully!")
      } else{
        return("There was an issue in setting up the venv's Python interpreter. Please see the vignette for troubleshooting tips. Exiting COBRA...")
      }
    } 
    
    else{
      #check current python interpreter
      use_virtualenv("COBRA")
      
      #have to activate venv by importing a module contained in it
      pkg_check<-tryCatch({a<-import('mhcflurry')}, error=function(cond){return("error")})
      #check current python interpreter
      env_check<-system("which python3", intern = T)
      
      if(grepl("COBRA", env_check) & pkg_check != "error"){
        lgr$info("Previously configured virtual environment activated successfully!")
      } else{
        return("The expected venv, 'COBRA', was detected in the list of virtual environment's, but is not currently active. Please see the vignette on how to activate the venv. Exiting COBRA...")
        }
    }
}
