
###################
#### FUNCTIONS ####
###################

# This function retrieve a xset like object
#@author Gildas Le Corguille lecorguille@sb-roscoff.fr
getxcmsSetObject <- function(xobject) {
    # XCMS 1.x
    if (class(xobject) == "xcmsSet")
        return (xobject)
    # XCMS 3.x
    if (class(xobject) == "XCMSnExp") {
        # Get the legacy xcmsSet object
        suppressWarnings(xset <- as(xobject, 'xcmsSet'))
        if (!is.null(xset@phenoData$sample_group))
            sampclass(xset) <- xset@phenoData$sample_group
        else
            sampclass(xset) <- "."
        return (xset)
    }
}

# This function get the raw file path from the arguments
#@author Gildas Le Corguille lecorguille@sb-roscoff.fr
getRawfilePathFromArguments <- function(singlefile, zipfile, args, prefix="") {
  if (!(prefix %in% c("","Positive","Negative","MS1","MS2"))) stop("prefix must be either '', 'Positive', 'Negative', 'MS1' or 'MS2'")

  if (!is.null(args[[paste0("zipfile",prefix)]])) zipfile <- args[[paste0("zipfile",prefix)]]

  if (!is.null(args[[paste0("singlefile_galaxyPath",prefix)]])) {
    singlefile_galaxyPaths <- args[[paste0("singlefile_galaxyPath",prefix)]]
    singlefile_sampleNames <- args[[paste0("singlefile_sampleName",prefix)]]
  }
  if (exists("singlefile_galaxyPaths")){
    singlefile_galaxyPaths <- unlist(strsplit(singlefile_galaxyPaths,"\\|"))
    singlefile_sampleNames <- unlist(strsplit(singlefile_sampleNames,"\\|"))

    singlefile <- NULL
    for (singlefile_galaxyPath_i in seq(1:length(singlefile_galaxyPaths))) {
      singlefile_galaxyPath <- singlefile_galaxyPaths[singlefile_galaxyPath_i]
      singlefile_sampleName <- singlefile_sampleNames[singlefile_galaxyPath_i]
      # In case, an url is used to import data within Galaxy
      singlefile_sampleName <- tail(unlist(strsplit(singlefile_sampleName,"/")), n=1)
      singlefile[[singlefile_sampleName]] <- singlefile_galaxyPath
    }
  }
  return(list(zipfile=zipfile, singlefile=singlefile))
}

# This function retrieve the raw file in the working directory
#   - if zipfile: unzip the file with its directory tree
#   - if singlefiles: set symlink with the good filename
#@author Gildas Le Corguille lecorguille@sb-roscoff.fr
retrieveRawfileInTheWorkingDirectory <- function(singlefile, zipfile) {
    if(!is.null(singlefile) && (length("singlefile")>0)) {
        for (singlefile_sampleName in names(singlefile)) {
            singlefile_galaxyPath <- singlefile[[singlefile_sampleName]]
            if(!file.exists(singlefile_galaxyPath)){
                error_message <- paste("Cannot access the sample:",singlefile_sampleName,"located:",singlefile_galaxyPath,". Please, contact your administrator ... if you have one!")
                print(error_message); stop(error_message)
            }

            if (!suppressWarnings( try (file.link(singlefile_galaxyPath, singlefile_sampleName), silent=T)))
                file.copy(singlefile_galaxyPath, singlefile_sampleName)

        }
        directory <- "."

    }
    if(!is.null(zipfile) && (zipfile != "")) {
        if(!file.exists(zipfile)){
            error_message <- paste("Cannot access the Zip file:",zipfile,". Please, contact your administrator ... if you have one!")
            print(error_message)
            stop(error_message)
        }

        #list all file in the zip file
        #zip_files <- unzip(zipfile,list=T)[,"Name"]

        #unzip
        suppressWarnings(unzip(zipfile, unzip="unzip"))

        #get the directory name
        suppressWarnings(filesInZip <- unzip(zipfile, list=T))
        directories <- unique(unlist(lapply(strsplit(filesInZip$Name,"/"), function(x) x[1])))
        directories <- directories[!(directories %in% c("__MACOSX")) & file.info(directories)$isdir]
        directory <- "."
        if (length(directories) == 1) directory <- directories

        cat("files_root_directory\t",directory,"\n")

    }
    return (directory)
}

# This function check if XML contains special caracters. It also checks integrity and completness.
#@author Misharl Monsoor misharl.monsoor@sb-roscoff.fr ABiMS TEAM
checkXmlStructure <- function (directory) {
    cat("Checking XML structure...\n")

    cmd <- paste0("IFS=$'\n'; for xml in $(find '",directory,"' -not -name '\\.*' -not -path '*conda-env*' -type f -iname '*.*ml*'); do if [ $(xmllint --nonet --noout \"$xml\" 2> /dev/null; echo $?) -gt 0 ]; then echo $xml;fi; done;")
    capture <- system(cmd, intern=TRUE)

    if (length(capture)>0){
        #message=paste("The following mzXML or mzML file is incorrect, please check these files first:",capture)
        write("\n\nERROR: The following mzXML or mzML file(s) are incorrect, please check these files first:", stderr())
        write(capture, stderr())
        stop("ERROR: xcmsSet cannot continue with incorrect mzXML or mzML files")
    }

}


# This function check if XML contain special characters
#@author Misharl Monsoor misharl.monsoor@sb-roscoff.fr ABiMS TEAM
deleteXmlBadCharacters<- function (directory) {
    cat("Checking Non ASCII characters in the XML...\n")

    processed <- F
    l <- system( paste0("find '",directory, "' -not -name '\\.*' -not -path '*conda-env*' -type f -iname '*.*ml*'"), intern=TRUE)
    for (i in l){
        cmd <- paste("LC_ALL=C grep '[^ -~]' \"", i, "\"", sep="")
        capture <- suppressWarnings(system(cmd, intern=TRUE))
        if (length(capture)>0){
            cmd <- paste("perl -i -pe 's/[^[:ascii:]]//g;'",i)
            print( paste("WARNING: Non ASCII characters have been removed from the ",i,"file") )
            c <- system(cmd, intern=TRUE)
            capture <- ""
            processed <- T
        }
    }
    if (processed) cat("\n\n")
    return(processed)
}


# This function will compute MD5 checksum to check the data integrity
#@author Gildas Le Corguille lecorguille@sb-roscoff.fr
getMd5sum <- function (directory) {
    cat("Compute md5 checksum...\n")
    # WHAT XCMS WILL FIND
    filepattern <- c("[Cc][Dd][Ff]", "[Nn][Cc]", "([Mm][Zz])?[Xx][Mm][Ll]","[Mm][Zz][Dd][Aa][Tt][Aa]", "[Mm][Zz][Mm][Ll]")
    filepattern <- paste(paste("\\.", filepattern, "$", sep=""),collapse="|")
    info <- file.info(directory)
    listed <- list.files(directory[info$isdir], pattern=filepattern, recursive=TRUE, full.names=TRUE)
    files <- c(directory[!info$isdir], listed)
    exists <- file.exists(files)
    files <- files[exists]

    library(tools)

    #cat("\n\n")

    return(as.matrix(md5sum(files)))
}

# This function check if xcms will found all the files
#@author Gildas Le Corguille lecorguille@sb-roscoff.fr ABiMS TEAM
checkFilesCompatibilityWithXcms <- function(directory) {
    cat("Checking files filenames compatibilities with xmcs...\n")
    # WHAT XCMS WILL FIND
    filepattern <- c("[Cc][Dd][Ff]", "[Nn][Cc]", "([Mm][Zz])?[Xx][Mm][Ll]","[Mm][Zz][Dd][Aa][Tt][Aa]", "[Mm][Zz][Mm][Ll]")
    filepattern <- paste(paste("\\.", filepattern, "$", sep=""),collapse="|")
    info <- file.info(directory)
    listed <- list.files(directory[info$isdir], pattern=filepattern, recursive=TRUE, full.names=TRUE)
    files <- c(directory[!info$isdir], listed)
    files_abs <- file.path(getwd(), files)
    exists <- file.exists(files_abs)
    files[exists] <- files_abs[exists]
    files[exists] <- sub("//","/",files[exists])

    # WHAT IS ON THE FILESYSTEM
    filesystem_filepaths <- system(paste0("find \"",getwd(),"/",directory,"\" -not -name '\\.*' -not -path '*conda-env*' -type f -name \"*\""), intern=T)
    filesystem_filepaths <- filesystem_filepaths[grep(filepattern, filesystem_filepaths, perl=T)]

    # COMPARISON
    if (!is.na(table(filesystem_filepaths %in% files)["FALSE"])) {
        write("\n\nERROR: List of the files which will not be imported by xcmsSet",stderr())
        write(filesystem_filepaths[!(filesystem_filepaths %in% files)],stderr())
        stop("\n\nERROR: One or more of your files will not be import by xcmsSet. It may due to bad characters in their filenames.")
    }
}