library(dplyr)
filename1 <- "./Téléchargements/Galaxy221-[msPurity.frag4feature_on_data_216_and_data_30__tsv].tsv"
filename2 <- "./Téléchargements/Galaxy278-[W4M_f4f_on_data_28,_data_216,_and_data_30__tsv].tsv"
file1 <- read.csv(filename1, sep="\t", header= TRUE)
file2 <- read.csv(filename2, sep="\t", header = TRUE)
colnames <- c("mzmed","mzmin","mzmax","rtmed","rtmin","rtmax","npeaks")
file1notin2 <- anti_join(file1,file2,by=colnames)
file2notin1 <- anti_join(file2,file1,by=colnames)
same <- semi_join(file1,file2,by=colnames)