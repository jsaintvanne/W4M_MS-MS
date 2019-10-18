#!/usr/bin/env Rscript
#W4M_run_f4f.R version 0.0.10
#created by Julien Saint-Vanne
xset_pa_filename_fix <- function(opt, pa, xset=NULL){
  if (!is.null(opt$mzML_files) && !is.null(opt$galaxy_names)){
    # NOTE: Relies on the pa@fileList having the names of files given as 'names' of the variables
    # needs to be done due to Galaxy moving the files around and screwing up any links to files

    filepaths <- trimws(strsplit(opt$mzML_files, ',')[[1]])
    filepaths <- filepaths[filepaths != ""]
    new_names <- basename(filepaths)

    galaxy_names <- trimws(strsplit(opt$galaxy_names, ',')[[1]])
    galaxy_names <- galaxy_names[galaxy_names != ""]

    nsave <- names(pa@fileList)
    old_filenames  <- basename(pa@fileList)

    pa@fileList <- filepaths[match(names(pa@fileList), galaxy_names)]
    names(pa@fileList) <- nsave

    pa@puritydf$filename <- basename(pa@fileList[match(pa@puritydf$filename, old_filenames)])
    pa@grped_df$filename <- basename(pa@fileList[match(pa@grped_df$filename, old_filenames)])
  }

  if(!is.null(xset)){
    if(!all(basename(pa@fileList)==basename(xset@filepaths))){
       if(!all(names(pa@fileList)==basename(xset@filepaths))){
          print('FILELISTS DO NOT MATCH')
          message('FILELISTS DO NOT MATCH')
          quit(status = 1)
       }else{
          xset@filepaths <- unname(pa@fileList)
       }
   }
 }

  return(list(pa, xset))
}

# ----- LOG FILE -----
log_file <- file("log.txt", open = "wt")
sink(log_file)
sink(log_file, type = "output")


# ----- PACKAGE -----
cat("\tSESSION INFO\n")

#Import the different functions
source_local <- function(fname) {
    argv <- commandArgs(trailingOnly = FALSE)
    base_dir <- dirname(substring(argv[grep("--file=", argv)], 8))
    source(paste(base_dir, fname, sep="/"))
}
source_local("lib.R")

#@author G. Le Corguille
# This function will
# - load the packages
# - display the sessionInfo
loadAndDisplayPackages <- function(pkgs) {
    for(pkg in pkgs) suppressPackageStartupMessages( stopifnot( library(pkg, quietly=TRUE, logical.return=TRUE, character.only=TRUE)))

    sessioninfo = sessionInfo()
    cat(sessioninfo$R.version$version.string,"\n")
    cat("Main packages:\n")
    for (pkg in names(sessioninfo$otherPkgs)) { cat(paste(pkg,packageVersion(pkg)),"\t") }; cat("\n")
    cat("Other loaded packages:\n")
    for (pkg in names(sessioninfo$loadedOnly)) { cat(paste(pkg,packageVersion(pkg)),"\t") }; cat("\n")
}

pkgs <- c("R.utils","batch") #"batch" necessary for parseCommandArgs function
loadAndDisplayPackages(pkgs)

#To have my own msPurity package
sourceDirectory("/home/jsaintvanne/W4M/msPurity/R")
#source("/home/jsaintvanne/W4M/msPurityW4M/R/all-generics.R")
#source("/home/jsaintvanne/W4M/msPurityW4M/R/averaging.R")
#source("/home/jsaintvanne/W4M/msPurityW4M/R/combine_annotations.R")
#source("/home/jsaintvanne/W4M/msPurityW4M/R/create_database.R")
#source("/home/jsaintvanne/W4M/msPurityW4M/R/flag-filter-remove.R")
#source("/home/jsaintvanne/W4M/msPurityW4M/R/iw-norm.R")
#source("/home/jsaintvanne/W4M/msPurityW4M/R/pcalc.R")
#source("/home/jsaintvanne/W4M/msPurityW4M/R/purityA-0-class.R")
#source("/home/jsaintvanne/W4M/msPurityW4M/R/purityA-av-spectra.R")
#source("/home/jsaintvanne/W4M/msPurityW4M/R/purityA-constructor.R")
#source("/home/jsaintvanne/W4M/msPurityW4M/R/purityA-create-msp.R")
#source("/home/jsaintvanne/W4M/msPurityW4M/R/purityA-filter-frag-spectra.R")
#source("/home/jsaintvanne/W4M/msPurityW4M/R/purityA-frag4feature.R")
#source("/home/jsaintvanne/W4M/msPurityW4M/R/spectral-complexity.R")
#source("/home/jsaintvanne/W4M/msPurityW4M/R/spectral_matching.R")
#source("/home/jsaintvanne/W4M/msPurityW4M/R/splinepurity.R")

modNamC <- "W4M_f4f" ## module name

cat("\nStart of the '", modNamC, "' Galaxy module call: ", format(Sys.time(), "%a %d %b %Y %X"), "\n\n", sep="")


# ----- ARGUMENTS -----
cat("\tARGUMENTS INFO\n\n")
args = parseCommandArgs(evaluate=FALSE) #interpretation of arguments given in command line as an R list of objects
#write.table(as.matrix(args), col.names=F, quote=F, sep='\t\t')
print(cbind(value = unlist(args)))


# ----- PROCESSING INFILE -----
cat("\n\n\tARGUMENTS PROCESSING INFO\n\n")

loadRData <- function(rdata_path, name){
#loads an RData file, and returns the named xset object if it is there
    load(rdata_path)
    return(get(ls()[ls() %in% name]))
}

# Saving the specific parameters
#pa parameter
#purityA object
if(args$pa != "NULL"){
	pa <- loadRData(args$pa, "pa")
}else{
	error_message <- "No pa object in your RData path enterred as pa parameter\n"
	cat(error_message)
	stop(error_message)
}

#xset parameter
#xcmsSet object derived from the same files as those used to create the purityA object
if(args$xset != "NULL"){
	load(args$xset)
  #xset <- loadRData(args$xset, c("xset","xdata"))
	#Transform XCMS object if needed
  if(!exists("xset")) {
    if(exists("xdata")) {
      xset<-getxcmsSetObject(xdata)
    } else {
      error_message="no xset and no xdata... Probably a problem"
      print(error_message)
      stop(error_message)
    }
  }
}else{
	error_message <- "No xset or xdata object in your RData path enterred as xset parameter\n"
	cat(error_message)
	stop(error_message)
}

#sampleMetadata parameter
#
if(!(is.null(args$sampleMetadata))){
  sampleMetadata <- args$sampleMetadata
}else{
  sampleMetadata <- NULL
}
#ppm parameter
#ppm tolerance between precursor mz and XCMS feature mz
if(args$ppm != "NULL"){
	ppm <- args$ppm
}else{
	ppm <- 5
}
#purity threshold parameter
#minimum purity of precursor to be included
if(args$plim != "NULL"){
	plim <- args$plim
}else{
	plim <- 0
}
#mostIntense parameter
#If TRUE the most intense precursor will be used.
#If FALSE the precursor closest to the center of the isolation window will be used
if(args$mostIntense == "true"){
	mostIntense <- TRUE
}else{
	mostIntense <- FALSE
}
#convert2RawRT parameter
#If retention time correction has been used in XCMS set this to TRUE
convert2RawRT = FALSE
#useGroup parameter
#Ignore individual peaks and just find matching fragmentation spectra within the (full) rtmin rtmax of each grouped feature
if(args$useGroup == "false"){
  useGroup = FALSE
}else{
  useGroup = TRUE
}
#createDB parameter
#(Deprecated, to be removed - use createDatabase function) SQLite database will be created of the results
createDB = FALSE
#out_dir parameter
#(Deprecated, to be removed - use createDatabase function) Path where database will be created
out_dir = "."
#db_name parameter
#(Deprecated, to be removed - use createDatabase function) If create_db is TRUE, a custom database name can be used, default is a time stamp
db_name = NA
#grp_peaklist parameter
#(Deprecated, to be removed - use createDatabase function) Can use any peak dataframe to add to databse. Still needs to be derived from the xset object though
grp_peaklist = NA
#use_group parameter
#(Deprecated, to be removed - replaced with useGroup argument for style consistency)
use_group = useGroup

cat("\n\n")


# ----- INFILE PROCESSING -----
cat("\tINFILE PROCESSING INFO\n\n")

#Correction of filenames (if Galaxy uses)
if(!useGroup){
  fix <- xset_pa_filename_fix(args, pa, xset)
  pa <- fix[[1]]
  xset <- fix[[2]]
  useGroup=FALSE
}else{
  # if are only aligning to the group not each file we do not need to align the files between the xsettempo and pa object
  fix <- xset_pa_filename_fix(args, pa)
  pa <- fix[[1]]
  useGroup=TRUE
}
print(pa@fileList)
print(xset@filepaths)


# ----- MAIN PROCESSING INFO -----
cat("\n\n\tMAIN PROCESSING INFO\n\n")
# Handle infiles
if (!exists("singlefile")) singlefile <- NULL
if (!exists("zipfile")) zipfile <- NULL
rawFilePath <- getRawfilePathFromArguments(singlefile, zipfile, args)
zipfile <- rawFilePath$zipfile
singlefile <- rawFilePath$singlefile
directory <- retrieveRawfileInTheWorkingDirectory(singlefile, zipfile)

# Check some character issues
md5sumList <- list("origin" = getMd5sum(directory))
checkXmlStructure(directory)
checkFilesCompatibilityWithXcms(directory)

cat("\n\n")

cat("\t\tCOMPUTE\n\n")

source("/home/jsaintvanne/W4M/W4M_MS-MS/W4M_run_f4f.R")
final_pa <- W4M_frag4feature(pa = pa, xset = xset, sampleMetadata = sampleMetadata, ppm = ppm, plim = plim, intense = mostIntense, 
                      convert2RawRT = convert2RawRT, useGroup = useGroup, create_db = FALSE, out_dir = '.', db_name = NA, 
                      grp_peaklist = NA, use_group = NA)

print(final_pa)


dir.create("tsv",showWarnings=FALSE)
dir.create("rdata",showWarnings=FALSE)
for(o in 1:length(final_pa)){
  cat("Saving object for ",names(final_pa)[o]," class...\n")
  pa <- final_pa[[o]]
  outputframe <- pa@grped_df
  #Order by mzmed
  outputframe <- outputframe[order(outputframe[,1]),]
  #if(use_group){
  #if(pa@f4f_link_type == 'group'){
  #    cols.dont.want <- NULL
    #cols.dont.want <- c("pid", "precurMtchID", "fileid") # if you want to remove multiple columns
  #    outputdata <- outputdata[, ! names(outputdata) %in% cols.dont.want, drop = F]
  #}else{
  #    print("la")
  #    cols.dont.want <- c("sample", "is_filled", "cid", "pid", "precurMtchID", "fileid") # if you want to remove multiple columns
  #    outputdata <- outputdata[, ! names(outputdata) %in% cols.dont.want, drop = F]
  #}
  file_path <- paste0('tsv/f4fresult_for_',names(final_pa)[o],"_class.tsv")
  write.table(outputframe, file=paste0(args$out_dir, "/",file_path), row.names=FALSE, sep='\t')

  
  save(pa,file=paste0(args$out_dir, "/rdata/f4fresult_for_",names(final_pa)[o],"_class.Rdata"))
}
#Build a comprehensive table for users
#outputdata <- pa@grped_df



cat("\nEnd of the '", modNamC, "' Galaxy module call: ", format(Sys.time(), "%a %d %b %Y %X"), "\n\n", sep="")

cat("\n\t\tDONE\n\n")


