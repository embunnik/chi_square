## Probability of observed spectral counts in capture experiments (Bunnik et al.) based on cellular protein abundance (Oehring et al.)

#### Data input is of the format:

## column_1		column_2		column_3		column_4
## spectra_capture	spectra_Oehring	total_spectra_capture	total_spectra_Oehring
## spectra_capture	spectra_Oehring	total_spectra_capture	total_spectra_Oehring
## spectra_capture	spectra_Oehring	total_spectra_capture	total_spectra_Oehring
## spectra_capture	spectra_Oehring	total_spectra_capture	total_spectra_Oehring
## spectra_capture	spectra_Oehring	total_spectra_capture	total_spectra_Oehring
## ....			....			....			.... 

## where: 
## spectra_capture is the number of spectra for a particular protein in all capture experiments combined 
## spectra_Oehring is the number of spectra for a particular protein in all fractions analyzed by Oehring et al. 
## total_spectra_capture is the total number of spectra for all proteins in all capture experiments combined
## total_spectra_Oehring is the total number of spectra for all proteins in all fractions analyzed by Oehring et al. 
##
## 
## store the function for the chi-square test 
## read in csv file containing data input, in above format
## perform the chi-square test on each row of the dataframe
## obtain a list of raw p-values
## adjust raw p-value using Benjami-Hochberg false discovery rate correction
## write csv file with the results, including the counts 

get_chisq <- function(data_input){
  mat <- matrix(as.numeric(data_input[c(1:4)]), ncol=2)
  f <- chisq.test(as.table(mat))
  return(c(data_input[1], f$p.value))
}

data_input <- read.csv("data_input_file.csv")
chisqs <- apply(data_input, 1, get_chisq)
p_values <- as.numeric(chisqs[2,])
BH_adjusted_p_values <- p.adjust(p_values,method="BH")
write.csv(cbind(data_input,BH_adjusted_p_values), file="chisq_result.csv", quote=FALSE)
