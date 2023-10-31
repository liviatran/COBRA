##CIWD_C1_netmhcpan
#'HLA Class I CIWD Allele files
#'
#'A folder containing all Common, Intermediate, Well-Documented (CIWD) v 3.0.0 Class I
#'HLA Alleles. These files will be used as input if NetMHCpan is selected as the
#'Binding Affinity Predictor for Class I alleles. The source file was parsed
#'to only include Class I CIWD alleles, to only alleles allowed by netMHCpan, 
#'and non-expression alleles (N, S, C, A and Q) were filtered out. The comprehensive
#'list of Class I CIWD HLA alleles was broken up to generate files of 45 alleles per file,
#'due to netMHCpan's maximum allele evaluation capacity. 
#' @docType data
#' @name CIWD_C1_netmhcpan
#' @format a folder with 46 .txt files containing HLA CLass I CIWD alelles
#' @source CIWD_CLASS_Palleles-all_loci-2020320-ihws-website.xlsx from https://www.ihiw18.org/component-immunogenetics/download-common-and-well-documented-alleles-3-0/
#' @references Hurley et. al. "Common, intermediate and well-documented HLA alleles in world populations: CIWD version 3.0.0". HLA (2020), Volume 95, Issue 6, pgs 503-637
NULL

##CIWD_C1_mhcflurry
#'HLA Class I CIWD Allele files
#'
#'A folder containing all Common, Intermediate, Well-Documented (CIWD) v 3.0.0 Class I
#'HLA Alleles. These files will be used as input if mhcflurry is selected as the
#'Binding Affinity Predictor for Class I alleles. The source file was parsed
#'to only include CIWD alleles, to only alleles allowed by mhcflurry,
#'and non-expression alleles (N, S, C, A and Q) were filtered out. The comprehensive
#'list of Class I CIWD HLA alleles was broken up to generate files of 45 alleles per file.
#'These files are formatted differently than CIWD_C1_netmhcpan files in that they have spaces 
#'instead of commas as the delimiter, and do not contain HLA-C12:139, as this allele is not supported by NetCleave.
#' @docType data
#' @name CIWD_C1_mhcflurry
#' @format a folder with 46 .txt files containing HLA Class I CIWD alelles
#' @source CIWD_CLASSI_Palleles-all_loci-2020320-ihws-website.xlsx from https://www.ihiw18.org/component-immunogenetics/download-common-and-well-documented-alleles-3-0/
#' @references Hurley et. al. "Common, intermediate and well-documented HLA alleles in world populations: CIWD version 3.0.0". HLA (2020), Volume 95, Issue 6, pgs 503-637
NULL

##CIWD_C2
#'HLA Class II CWD/CIWD alleles
#'
#'A folder containing all CIWD v 3.0.0 DRB1, DQB1, and DPB1 alleles, and CWD v 2.0.0
#'CWD v 2.0.0 DPA1 and DQA1 alleles. The source files were parsed
#'to only include Class II CIWD/CWD alleles, to only alleles allowed by netMHCIIpan, 
#'and non-expression alleles (N, S, C, A and Q) were filtered out. The 
#'comprehensive list of CIWD/CWD Class II HLA alleles was broken up to generate 
#'20 alleles per file to accommodate netMHCIIPan’s maximum allele evaluation capacity.  
#'The source files were also manipulated to create all possible heterodimers for DP and DQ, 
#'including forbidden heterodimers. 
#' @docType data
#' @name CIWD_C2
#' @format a folder with 105 .txt files containing HLA Class II CIWD/CWD alleles
#' @source CIWD_CLASS_Palleles-all_loci-2020320-ihws-website.xlsx from https://www.ihiw18.org/component-immunogenetics/download-common-and-well-documented-alleles-3-0/
#' @source cwd200.xls from https://www.uvm.edu/~igdawg/cwd.html
#' @references Hurley et. al. "Common, intermediate and well-documented HLA alleles in world populations: CIWD version 3.0.0". HLA (2020), Volume 95, Issue 6, pages 503-637
#' @references Mack et. al. "Common and well-documented HLA alleles: 2012 update to the CWD catalogue". Tissue Antigens (2013), Volume 8, Issue 4 pages 194-203
NULL

##ref_FASTA
#'SARS-CoV-2 Reference Sequence FASTA Files
#'
#'A folder containing FASTA files the SARS-CoV-2 reference, or Wuhan-Hu-1 proteome, which includes the
#'Envelope, Membrane, Nucleocapsid, Spike, nsp1-nsp11, nsp13-16, ORF3A, ORF6, ORF8, 
#'ORF10, ORF7a, ORF7b, and RdRp proteins. A custom script was used to scrape the 
#'FASTA sequences from https://www.ncbi.nlm.nih.gov/nuccore/NC_045512
#' @docType data
#' @name ref_FASTA
#' @format a folder with 26 FASTA files containing FASTA sequences for the SARS-CoV-2
#' reference proteome. 
#' @source https://www.ncbi.nlm.nih.gov/nuccore/NC_045512
NULL

##beta_FASTA
#'SARS-CoV-2 Beta Sequence FASTA Files
#'
#'A folder containing FASTA files the SARS-CoV-2 Beta proteome, which includes the
#'Envelope, Membrane, Nucleocapsid, Spike, nsp1-nsp11, nsp13-16, ORF3A, ORF6, ORF8, 
#'ORF10, ORF7a, ORF7b, and RdRp proteins. A custom script was used to scrape the 
#'FASTA sequences from https://www.ncbi.nlm.nih.gov/nuccore/MW598419.1. Note these
#'are user submitted sequences to NCBI, and are not officially published by NCBI
#' @docType data
#' @name beta_FASTA
#' @format a folder with 26 FASTA files containing FASTA sequences for the SARS-CoV-2
#' Beta proteome. 
#' @source https://www.ncbi.nlm.nih.gov/nuccore/MW598419.1
NULL

##omiba1_FASTA
#'SARS-CoV-2 Omicron BA.1 Sequence FASTA Files
#'
#'A folder containing FASTA files the SARS-CoV-2 Omicron BA.1 proteome, which includes the
#'Envelope, Membrane, Nucleocapsid, Spike, nsp1-nsp11, nsp13-16, ORF3A, ORF6, ORF8, 
#'ORF10, ORF7a, ORF7b, and RdRp proteins. A custom script was used to scrape the 
#'FASTA sequences from https://www.ncbi.nlm.nih.gov/nuccore/OL672836. Note these
#'are user submitted sequences to NCBI, and are not officially published by NCBI
#' @docType data
#' @name omniba1_FASTA
#' @format a folder with 26 FASTA files containing FASTA sequences for the SARS-CoV-2
#' Omicron BA.1 proteome. 
#' @source https://www.ncbi.nlm.nih.gov/nuccore/OL672836
NULL

##Lookup_Tables
#'Lookup Tables
#'
#'A folder containing look up tables generated for the SARS-CoV-2 reference, Beta, and Omicron BA.1 proteomes
#' @docType data
#' @name Lookup_Tables
#' @format a folder with 3 look up tables generated for the SARS-CoV-2 reference, Beta, and Omicron BA.1 proteomes

