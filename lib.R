
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

# 2 functions to have a sampleMetadata looks like this 
#(with no extension for files and just basename) :
# MSMS 	   class
#file1    class1
#file2    class2
#file3    class1
#First from xcmsSet object
buildSamplemetadataFromXCMS <- function(xset){
	if(!(is.null(xset))){
  		#Create sampleMetadata for MSMS files with each one and its class (only MSMS)
  		MSMS <- NULL
  		classMSMS <- NULL
  		filepattern <- c("[Cc][Dd][Ff]", "[Nn][Cc]", "([Mm][Zz])?[Xx][Mm][Ll]", 
                         "[Mm][Zz][Dd][Aa][Tt][Aa]", "[Mm][Zz][Mm][Ll]")
        filepattern <- paste(paste("\\.", filepattern, "$", sep = ""), collapse = "|")
  		for(i in 1:length(xset@filepaths)) {
    		#Select only MSMS files
    		if(2 %in% unique(readMSData(xset@filepaths[i],mode="onDisk")@featureData@data$msLevel)) {
            	sampname <- gsub(filepattern, "",basename(xset@filepaths[i])) 
        		MSMS <- c(MSMS, sampname)
        		classMSMS <- c(classMSMS,as.character(xset@phenoData[i,"class"]))
    		}
  		}
  		MSMSclassfile <- data.frame(MSMS,classMSMS)
  		return(MSMSclassfile)
	}else{
		error_message <- "No xset available in buildSamplemetadataFromXCMS\n"
		cat(error_message)
		stop(error_message)
		return(NULL)
  	}	
}
#Second from file
buildSamplemetadataFromFile <- function(MSMSclassfile,xset){
	finaleMSMSclasses <- NULL
	if(!(is.null(MSMSclassfile))){
  		#Read the MSMSclasses
  		MSMSclasses <- read.csv(file=MSMSclassfile, sep="\t", header=TRUE)
  		#Verify the colnames
  		if(length(colnames(MSMSclasses)) == 2){
  			print("bon nombre de col")
  			if(colnames(MSMSclasses)[1] == "MSMS" && colnames(MSMSclasses)[2] == "class"){
  				print("nom de col good")
  				finaleMSMSclasses <- MSMSclasses
  			}else{
  				colnames(MSMSclasses) <- c("MSMS","class")
  				print(MSMSclasses)
  				finaleMSMSclasses <- MSMSclasses
  			}
  		}else{
  			print("mauvais nbre de col")
  			#Verify if MSMS col exists and keep it
  			if("MSMS" %in% colnames(MSMSclasses)){
  				MSMS <- MSMSclasses[,which(colnames(MSMSclasses) == "MSMS")]
  			}else{
  				error_message <- "No MSMS column in your sampleMetadata"
  				cat(error_message)
				stop(error_message)
				return(NULL)
  			}
  			#Verify if class col exists and keep it
  			if("class" %in% colnames(MSMSclasses)){
  				class <- MSMSclasses[,which(colnames(MSMSclasses) == "class")]
  			}else{
  				error_message <- "No class column in your sampleMetadata"
  				cat(error_message)
				stop(error_message)
				return(NULL)
  			}
  			MSMSclasses <- cbind(MSMS,class)
  			print(MSMSclasses)
  		}
  		#Keep only MSMS files
  		for(i in 1:length(xset@filepaths)) {
  			#remove .cdf, .mzXML filepattern
            filepattern <- c("[Cc][Dd][Ff]", "[Nn][Cc]", "([Mm][Zz])?[Xx][Mm][Ll]", 
                                    "[Mm][Zz][Dd][Aa][Tt][Aa]", "[Mm][Zz][Mm][Ll]")
            filepattern <- paste(paste("\\.", filepattern, "$", sep = ""), collapse = "|")
            sampname<-gsub(filepattern, "",basename(xset@filepaths[i]))  
    		if(2 %in% unique(readMSData(xset@filepaths[i],mode="onDisk")@featureData@data$msLevel)) {
        		addMSMS <- MSMSclasses[which(MSMSclasses$MSMS == sampname),]
        		print(addMSMS)
        		finaleMSMSclasses <- rbind(finaleMSMSclasses,addMSMS)
    		}
  		}
  		return(finaleMSMSclasses)
	}else{
  		error_message <- "No sampleMetadata file enter\n"
  		cat(error_message)
		stop(error_message)
		return(NULL)
	}
}
#Thirs from table to do ????


#This function create a tempo xcmsSet object with only peaks from class of MSMSfile we are looking at
modifyXsetObject <- function(xset, MSMSfile, class){
	cat(paste("--------------- Modify xcmsSet object for",MSMSfile,"from",class,"class ---------------\n"))
    cat("STEP 1 : find good groups from our class\n")
    #STEP 1 : select groups in xset@groups for the good class and save ids of groups deleted
    #Select groups which contain peaks in the same class as file class
    xset@phenoData <- xset@phenoData[which(gsub(filepattern,"",xset@phenoData[,"sample_name"]) == MSMSfile),]
    xset@filepaths <- xset@filepaths[which(basename(xset@filepaths) == rownames(xset@phenoData))]
    #Save line where we have one peak in one file of my class
    allsavedgroups <- xset@groups[which(xset@groups[,class] > 0),]
    #Save line number where 0 peaks is in my files of my class
    alldeletedgroups <- which(xset@groups[,class] == 0)
    cat("\tDeleting",length(alldeletedgroups),"group(s) from groups tabledata and stock",nrow(allsavedgroups),"group(s) from class",class,"from xset\n")

    cat("STEP 2 : add grpid for each peak and delete when they are not in the good group\n")
    #STEP 2 : Delete peaks from the wrong groups (groups to delete)   
    #Add new col for grpid number for each peak
    grpidcol <- rep(NA, length(xset@peaks))
    xset@peaks <- cbind(xset@peaks,grpidcol)
    #Complete grpid for each peak
    for(x in 1:length(xset@groupidx)){
        for(y in 1:length(xset@groupidx[[x]])){
            xset@peaks[xset@groupidx[[x]][y],"grpidcol"] <- x
        }
    }
    #Delete peaks where groupidx have to be delete
    peaktodelete <- NULL
    peaktodelete <- which(xset@peaks[,"grpidcol"] %in% alldeletedgroups)
    if(length(peaktodelete) > 0){
        xset@peaks <- xset@peaks[-peaktodelete,]
    }
    cat("\tWe now have",nrow(xset@peaks),"peaks and deleting",length(peaktodelete),"peaks which were not in groups where we have peak froms class \"",class,"\" !\n")
    
    cat("STEP 3 : Rebuild groupidx and groups\n")
    #STEP 3 : Rebuild groupidx with new row for groups and peaks (usefull for groupval after...)
    newgrpidx <- xset@groupidx
    newgrpidx <- lapply(newgrpidx, function(a) a <- NULL)
    for(p in 1:nrow(xset@peaks)){
        if(!(is.na(xset@peaks[p,"grpidcol"]))){
        	newgrpidx[[xset@peaks[p,"grpidcol"]]] <- c(newgrpidx[[xset@peaks[p,"grpidcol"]]],p)
        }else{
        	next
        }
        
    }
    if(length(alldeletedgroups) > 0){
        newgrpidx <- newgrpidx[-alldeletedgroups]
    }
    xset@groupidx <- newgrpidx
    cat("\tRebuild groupidx with",length(xset@groupidx),"groupidx ! \n")
    xset@groups <- allsavedgroups
    cat("\tRebuild groups with",nrow(xset@groups),"groups ! \n")

    #Delete line grpidcol after order them 
    xset@peaks <- subset(xset@peaks, select = -(which(colnames(xset@peaks)=="grpidcol")))

    return(xset)
}

#This function keep MSMS informations only from our MSMSfile we are working with
modifyPaObject <- function(pa, MSMSfile){
	cat(paste("--------------- Modify pa object for",MSMSfile,"---------------\n"))
	patempo <- pa
	cat("STEP 1 : modify pa@fileList\n")
	#patempo@fileList <- pa@fileList[which(gsub(filepattern,"",names(pa@fileList)) == MSMSfile)]

	cat("STEP 2 : modify pa@puritydf\n")
	fileIndex <- which(gsub(filepattern,"",names(pa@fileList)) == MSMSfile)
	patempo@puritydf <- pa@puritydf[which(pa@puritydf[,"fileid"] == fileIndex),]

	return(patempo)
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