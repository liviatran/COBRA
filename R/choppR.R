#choppR v 0.1.0 Mar2023
#'Runs NetChop on input FASTA files
#'
#'Runs NetChop on input FASTA files to predict likelihood of cleavage at a given
#'amino acid position. The NetChop output, named '<file>_chop.txt', where <file>
#'is the name of the input FASTA file, is output to a temp directory. The file is 
#'further processed, where extraneous text and white space are removed, and is output
#'to the temp directory as '<file>_chopadj.txt'. Finally, the processed file is
#'read into R, where column headers for multi-protein files are removed, row 
#'numbers are added as a column, and output to the temp directory as '<file>_netchopfinal.txt'.
#'NetChop parameters: threshold for specifying weak and strong binders was set to 0.1, and the
#'C-terminus prediction method was used.
#'
#'@param fastafile FASTA file for a given protein or group of proteins. Default 
#'value is set to the inst/extdata/ref_FASTA folder, which contains FASTA files
#'for the SARS-CoV2 reference proteins 
#'@param netchop_path The path to the user's installation of NetChop
#'
#'@importFrom dplyr %>%
#'@importFrom tibble add_column
#'@importFrom lgr lgr
#'@importFrom utils read.table write.table
#'
#'@return NULL. Otherwise, a string detailing the error that occurred. 
#'
#'@note The user's working directory MUST be set to where the input FASTA files are located if 
#'the bundled SARS-CoV2 FASTA files are not being used. 
#'@note NetChop MUST be downloaded and the path to its executable MUST be set up properly in order
#'for this function to work properly. In order to set up the full path to NetChop,
#'open the executable in a text editor, and enter the path to your NetChop download
#'where it says '# full path to the NetMHCpan 4.0 directory (mandatory)'. 
#'
#'@export
#'
#'@examples
#'#Example to run NetChop on a FASTA file named 'EMN_ref.faa'. Path should end in 
#'#/netchop-3.1/netchop'
#'\donttest{choppR("EMN_ref.faa", "path/to/NetChop/executable")}
#'
#'@source NetChop v 3.1, downloaded from https://services.healthtech.dtu.dk/services/NetMHCpan-4.1/
#'@references Nielsen M, Lundegaard C, Lund O, and Kesmir C. The role of the proteasome in generating cytotoxic T cell epitopes: Insights obtained from improved predictions of proteasomal cleavage. Immunogenetics., 57(1-2):33-41, 2005.

choppR<-function(fastafile = list.files(system.file("extdata/ref_FASTA", package = "COBRA"), full.names = TRUE), netchop_path){
    
    fileout<-fileout_adj<-netfile<-NULL
    
    for(k in 1:length(fastafile)){
      
      fileout[[k]]<-paste(tempdir(), "/", sapply(strsplit(basename(fastafile[[k]]), ".", fixed=TRUE), "[", 1), "_chop.txt", sep="")
      
      lgr$info(paste("Running NetChop on", gsub(".faa", "", basename(fastafile[[k]]))))
      
      system(paste(netchop_path, " ", fastafile[[k]], " -t 0.1 > ", fileout[[k]], sep=""))
      
      fileout_adj[[k]]<-paste(tempdir(), "/", sapply(strsplit(basename(fastafile[[k]]), ".", fixed=TRUE), "[", 1), "_chopadj.txt", sep="")
      
      system(paste("sed '/#/d; /-/d; /Number/d; /^$/d;' ", fileout[[k]], " > ", fileout_adj[[k]], sep=""))
      
      #read in netchop results
      #if path to netchop is not properly set up, error will be thrown since there
      #are no lines to read
      netfile[[k]]<-tryCatch({read.table(fileout_adj[[k]], header=TRUE, stringsAsFactors = F)}, error = function(error) {return("error")})
      
      #if netfile[[k]] is not a dataframe ('error' was returned), then NetChop
      #was not set up correctly
      if(!is.data.frame(netfile[[k]])){
        return("The full path to the NetChop directory was not set up correctly. Please make sure the environment for NetChop is set correctly. This can be done by using a text editor to open the NetChop executable, and providing the full path for NetChop in the setenv argument.")
      }
      
      netfile[[k]]<-netfile[[k]] %>%
        filter(!pos %in% "pos")
      
      netfile[[k]]<-netfile[[k]] %>%
        add_column(rownum = row.names(netfile[[k]]))
      
      lgr$info(paste("Outputting NetChop results for", gsub(".faa", "", basename(fastafile[[k]])), "to the temp directory"))
               
      #write final adjusted netchop file to temp directory 
      write.table(netfile[[k]], file = paste(tempdir(), "/", paste(gsub(".faa", "", basename(fastafile[[k]])), "_netchopfinal.txt", sep=""), sep =""), quote = F, row.names = F)
    }
}

