#' @title Using a purityA object, link MS/MS datas to XCMS features
#'
#' @description
#'
#' **General**
#'
#' Assign fragmentation spectra (MS/MS) stored within a purityA class object to grouped features within an XCMS xset object.
#'
#' XCMS calculates individual chromatographic peaks for each mzML file (saved in xset@@peaks), these are then grouped together
#' (using xcms.group). Ideally the mzML files that contain the MS/MS spectra also contain sufficient MS1 scans for XCMS to detect
#' MS1 chromatographic features. If this is the case, to determine if a MS2 spectra is to be linked to an XCMS grouped feature,
#' the associated acquisition time of the MS/MS event has to be within the retention time window defined for the individual peaks
#' associated for each file. The precursor m/z value also has to be within the user ppm tolerance to XCMS feature.
#'
#' See below for representation of the linking (the \*------\* represent a many-to-many relationship) e.g. 1 or more MS/MS events can be
#' linked to 1 or more individual feature and an individual XCMS feature can be linked to 1 or more grouped XCMS features
#'
#' * \[grouped XCMS feature - across files\] \*------\*  \[individual XCMS feature - per file\] \*------\*  \[MS/MS spectra\]
#'
#' Alternatively, if the "useGroup" argument is set to TRUE, the full width of the grouped peak (determined as the minimum rtmin
#' and maximum rtmax of the all associated individual peaks) will be used. This option should be used if the mzML file with
#' MS/MS has very limited MS1 data and so individual chromatographic peaks might not be detected with the mzML files containing the
#' MS/MS data. However, it should be noted this may lead to potential inaccurate linking.
#'
#' * \[grouped XCMS peaks\] \*------\* \[MS/MS spectra\]
#'
#'
#' **Example LC-MS/MS processing workflow**
#'
#' The purityA object can be used for further processing including linking the fragmentation spectra to XCMS features, averaging fragmentation, database creation and spectral matching (from the created database). See below for an example workflow
#'
#'  * Purity assessments
#'    +  (mzML files) -> purityA -> (pa)
#'  * XCMS processing
#'    +  (mzML files) -> xcms.xcmsSet -> xcms.merge -> xcms.group -> xcms.retcor -> xcms.group -> (xset)
#'  * Fragmentation processing
#'    + (xset, pa) -> **frag4feature** -> filterFragSpectra -> averageAllFragSpectra -> createDatabase -> spectralMatching -> (sqlite spectral database)
#'
#' **Additional notes**
#'
#' * If using only a single file, then grouping still needs to be performed within XCMS before frag4feature can be used.
#' * Fragmentation spectra below a certain precursor ion purity can be be removed (see plim argument).
#' * A SQLite database can be created directly here but the functionality has been deprecated and the createDatabase function should now be used
#' * Can experience some problems when using XCMS version < 3 and obiwarp retention time correction.
#'
#'
#' @aliases frag4feature
#'
#' @param pa object; purityA object
#' @param xset object; xcmsSet object derived from the same files as those used to create the purityA object
#' @param ppm numeric; ppm tolerance between precursor mz and XCMS feature mz
#' @param plim numeric; minimum purity of precursor to be included
#' @param intense boolean; If TRUE the most intense precursor will be used. If FALSE the precursor closest to the center of the isolation window will be used
#' @param useGroup boolean; Ignore individual peaks and just find matching fragmentation spectra within the (full) rtmin rtmax of each grouped feature
#' @param convert2RawRT boolean; If retention time correction has been used in XCMS set this to TRUE
#' @param create_db boolean; (Deprecated, to be removed - use createDatabase function) SQLite database will be created of the results
#' @param db_name character; (Deprecated, to be removed - use createDatabase function) If create_db is TRUE, a custom database name can be used, default is a time stamp
#' @param out_dir character; (Deprecated, to be removed - use createDatabase function) Path where database will be created
#' @param grp_peaklist dataframe; (Deprecated, to be removed - use createDatabase function) Can use any peak dataframe to add to databse. Still needs to be derived from the xset object though
#' @param use_group boolean; (Deprecated, to be removed - replaced with useGroup argument for style consistency)
#' @return Returns a purityA object (pa) with the following slots populated:
#'
#' * pa@@grped_df: A dataframe of the grouped XCMS features linked to the associated fragmentation spectra precursor details is recorded here
#' * pa@@grped_ms2: A list of fragmentation spectra associated with each grouped XCMS feature is recorded here
#' * pa@@f4f_link_type: The linking method is recorded here (e.g. individual peaks or grouped - "useGroup=TRUE")
#'
#'
#' @examples
#'
#' sampleMetadata <- read.csv("my_sampleMetadata_file.csv", sep="\t", header=TRUE)
#' msmsPths <- list.files(system.file("extdata", "lcms",
#'                        "mzML", package="msPurityData"), full.names = TRUE,
#'                        pattern = "MSMS")
#' xset <- xcms::xcmsSet(msmsPths, nSlaves = 1)
#' xset <- xcms::group(xset)
#' xset <- xcms::retcor(xset)
#' xset <- xcms::group(xset)
#'
#' pa  <- purityA(msmsPths)
#' pa <- W4M_frag4feature(pa, xset,sampleMetadata)
#' @md
#' @export

W4M_frag4feature <- function(pa, xset, sampleMetadata=NULL, ppm=5, plim=NA, intense=TRUE, convert2RawRT=TRUE, useGroup=FALSE, create_db=FALSE,
                                out_dir='.', db_name=NA, grp_peaklist=NA, use_group=NA){
  
pkgs <- c("R.utils","batch") #"batch" necessary for parseCommandArgs function
loadAndDisplayPackages(pkgs)

#Need it to run my own msPurity (I changed just 1 thing I have to discuss !)
sourceDirectory("/home/jsaintvanne/W4M/msPurity/R")

xset <- getxcmsSetObject(xset)

#Plutot passer directement un tableau qu'un fichier ??
if(is.null(sampleMetadata)){
	sampleMetadataMSMS <- buildSamplemetadataFromXCMS(xset)
}else{
  #Verify if we have a file
  if(class(sampleMetadata) == "character"){
    print("sampleMetadata est un character")
    sampleMetadataMSMS <- buildSamplemetadataFromFile(sampleMetadata,xset)
  }else{
    sampleMetadataMSMS <- buildSamplemetadataFromTable(sampleMetadata,xset)
  }
}

#Verify if sampleMetadataMSMS is null !
if(nrow(sampleMetadataMSMS) == 0){
	error_message <- "No file has msLevel = 2 in your sampleMetadata...\n"
	cat(error_message)
	stop(error_message)
}
print(sampleMetadataMSMS)
for(line in 1:nrow(sampleMetadataMSMS)){
	
    use_group=TRUE
	convert2RawRT = FALSE
print(sampleMetadataMSMS[line,])
    #Save class and filename
    MSMSfile <- as.character(sampleMetadataMSMS[line,"MSMS"])
    class <- as.character(sampleMetadataMSMS[line,"class"])

    #Modify the XcmsSet object to have only peaks and groups from the good class
    xsettempo <- modifyXsetObject(xset,MSMSfile,class)

    #Modify the pa object to have only the file we are working on!
    patempo <- modifyPaObject(pa,MSMSfile)

    #Find if retention time have been processed for this file
    #remove .cdf, .mzXML filepattern and find the file index number
    filepattern <- c("[Cc][Dd][Ff]", "[Nn][Cc]", "([Mm][Zz])?[Xx][Mm][Ll]", 
                     "[Mm][Zz][Dd][Aa][Tt][Aa]", "[Mm][Zz][Mm][Ll]")
    filepattern <- paste(paste("\\.", filepattern, "$", sep = ""), collapse = "|")
    for(f in 1:length(xset@filepaths)){
    	sampname <- gsub(filepattern, "",basename(xset@filepaths[f]))
    	if(sampname == MSMSfile){
    		fileIndexMSMS <- f
    		break
    	}
    }
    for(n in 1:length(xset@.processHistory)){
        if(fileIndexMSMS %in% fileIndex(xset@.processHistory[[n]])){
            if(xset@.processHistory[[n]]@type == "Retention time correction"){
                cat("Retention time used for",MSMSfile,"\n")
                convert2RawRT= TRUE
            }
        }
    }

sourceDirectory("/home/jsaintvanne/W4M/msPurity/R")
    #STEP 4 : Run f4f function
    cat("Run frag4feature in msPurity package\n")
    patempo <- frag4feature(pa=patempo, 
                       xset=xsettempo, 
                       ppm=ppm, 
                       plim=plim,
                       intense=intense,
                       convert2RawRT=convert2RawRT,
                       useGroup=useGroup,
                       create_db=createDB,
                       out_dir=out_dir,
                       db_name='alldata.sqlite',
                       grp_peaklist=grp_peaklist,
                       use_group=useGroup)

    #STEP 5 : save modified things from pa object
    #1-Save filename and precFilename
    patempo@grped_df <- tibble::add_column(patempo@grped_df,filename=basename(xset@filepaths[f]),.after=length(patempo@grped_df))
    #2-Save grped_df
    cat(paste("From",nrow(pa@grped_df),"to "))
    pa@grped_df <- rbind(pa@grped_df,patempo@grped_df)
    cat(paste(nrow(pa@grped_df),"rows in grped_df\n"))
    #3-Save grped_ms2
    cat(paste("From",length(pa@grped_ms2),"to "))
    pa@grped_ms2 <- c(pa@grped_ms2,patempo@grped_ms2)
    cat(paste(length(pa@grped_ms2),"rows in grped_ms2\n"))

}
return(pa)
}
###################
#### FUNCTIONS ####
###################

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

# 3 functions to have a sampleMetadata looks like this 
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
  		class <- NULL
  		filepattern <- c("[Cc][Dd][Ff]", "[Nn][Cc]", "([Mm][Zz])?[Xx][Mm][Ll]", 
                       "[Mm][Zz][Dd][Aa][Tt][Aa]", "[Mm][Zz][Mm][Ll]")
        filepattern <- paste(paste("\\.", filepattern, "$", sep = ""), collapse = "|")
  		for(i in 1:length(xset@filepaths)) {
    		#Select only MSMS files
    		if(2 %in% unique(readMSData(xset@filepaths[i],mode="onDisk")@featureData@data$msLevel)) {
            	sampname <- gsub(filepattern, "",basename(xset@filepaths[i])) 
        		MSMS <- c(MSMS, sampname)
        		class <- c(class,as.character(xset@phenoData[i,"class"]))
    		}
  		}
  		MSMSclassfile <- data.frame(MSMS,class)
  		return(MSMSclassfile)
	}else{
		error_message <- "No xset available in buildSamplemetadataFromXCMS\n"
		cat(error_message)
		stop(error_message)
		return(NULL)
    }	
}
#Second from file
buildSamplemetadataFromFile <- function(sampleMetadata,xset){
	finaleMSMSclasses <- NULL
	if(!(is.null(sampleMetadata))){
  		#Read the MSMSclasses
  		MSMSclasses <- read.csv(file=sampleMetadata, sep="\t", header=TRUE)
  		finaleMSMSclasses <- buildSamplemetadataFromTable(MSMSclasses, xset)
        return(finaleMSMSclasses)
	}else{
  		error_message <- "No sampleMetadata file enter\n"
  		cat(error_message)
		stop(error_message)
		return(NULL)
	}
}
#Third from table
buildSamplemetadataFromTable <- function(sampleMetadata, xset){
    finaleMSMSclasses <- NULL
    #Verify the colnames
    if(length(colnames(sampleMetadata)) == 2){
        print("bon nombre de col")
        if(colnames(sampleMetadata)[1] == "MSMS" && colnames(sampleMetadata)[2] == "class"){
            print("nom de col good")
        }else{
            colnames(sampleMetadata) <- c("MSMS","class")
        }
    }else{
        print("mauvais nbre de col")
        #Verify if MSMS col exists and keep it
        if("MSMS" %in% colnames(sampleMetadata)){
            MSMS <- sampleMetadata[,which(colnames(sampleMetadata) == "MSMS")]
        }else{
            error_message <- "No MSMS column in your sampleMetadata"
            cat(error_message)
            stop(error_message)
            return(NULL)
        }
        #Verify if class col exists and keep it
        if("class" %in% colnames(sampleMetadata)){
            class <- sampleMetadata[,which(colnames(sampleMetadata) == "class")]
        }else{
            error_message <- "No class column in your sampleMetadata"
            cat(error_message)
            stop(error_message)
            return(NULL)
        }
        sampleMetadata <- cbind(MSMS,class)
    }
    #Keep only MSMS files
    for(i in 1:length(xset@filepaths)){
        #remove .cdf, .mzXML filepattern
        filepattern <- c("[Cc][Dd][Ff]", "[Nn][Cc]", "([Mm][Zz])?[Xx][Mm][Ll]", 
                         "[Mm][Zz][Dd][Aa][Tt][Aa]", "[Mm][Zz][Mm][Ll]")
        filepattern <- paste(paste("\\.", filepattern, "$", sep = ""), collapse = "|")
        sampname<-gsub(filepattern, "",basename(xset@filepaths[i]))  
        if(2 %in% unique(readMSData(xset@filepaths[i],mode="onDisk")@featureData@data$msLevel)){
            print("ajout")
            addMSMS <- sampleMetadata[which(sampleMetadata$MSMS == sampname),]
            print(addMSMS)
            finaleMSMSclasses <- rbind(finaleMSMSclasses,addMSMS)
        }
    }
    return(finaleMSMSclasses)
}

#This function create a tempo xcmsSet object with only peaks from class of MSMSfile we are looking at
modifyXsetObject <- function(xset, MSMSfile, class){
    cat(paste("--------------- Modify xcmsSet object for",MSMSfile,"from",class,"class ---------------\n"))
    cat("STEP 1 : find good groups from our class\n")
    #STEP 1 : select groups in xset@groups for the good class and save ids of groups deleted
    #Select groups which contain peaks in the same class as file class
    filepattern <- c("[Cc][Dd][Ff]", "[Nn][Cc]", "([Mm][Zz])?[Xx][Mm][Ll]", 
                     "[Mm][Zz][Dd][Aa][Tt][Aa]", "[Mm][Zz][Mm][Ll]")
    filepattern <- paste(paste("\\.", filepattern, "$", sep = ""), collapse = "|")
    xset@phenoData <- xset@phenoData[which(gsub(filepattern,"",basename(xset@filepaths)) == MSMSfile),]
    #xset@filepaths <- xset@filepaths[which(basename(xset@filepaths) == rownames(xset@phenoData))]
    #Save line where we have one peak in one file of my class
    print(head(xset@groups))
    allsavedgroups <- xset@groups[which(xset@groups[,class] > 0),]
    print(head(allsavedgroups))
    #Save line number where 0 peaks is in my files of my class
    print(class)
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
    cat("\tWe now have",nrow(xset@peaks),"peaks and deleting",length(peaktodelete),"peaks which were not in groups where we have peak from class \"",class,"\" !\n")
    
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
print(head(pa@puritydf))
	cat("STEP 2 : modify pa@puritydf\n")
    filepattern <- c("[Cc][Dd][Ff]", "[Nn][Cc]", "([Mm][Zz])?[Xx][Mm][Ll]", 
                     "[Mm][Zz][Dd][Aa][Tt][Aa]", "[Mm][Zz][Mm][Ll]")
    filepattern <- paste(paste("\\.", filepattern, "$", sep = ""), collapse = "|")
	fileIndex <- which(gsub(filepattern,"",names(pa@fileList)) == MSMSfile)
print(MSMSfile)
print(pa@fileList)
    print(fileIndex)
	patempo@puritydf <- pa@puritydf[which(pa@puritydf[,"fileid"] == fileIndex),]
print(head(patempo@puritydf))
	return(patempo)
}