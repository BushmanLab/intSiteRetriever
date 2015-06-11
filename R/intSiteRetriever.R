.connectToDB <- function(dbConn){
  if(is.null(dbConn)){
    library("RMySQL") #also loads DBI
    all_cons <- dbListConnections(MySQL())
    for (con in all_cons) {
      discCon <- dbDisconnect(con)
    }
    dbConn <- dbConnect(MySQL(), group="intSitesDev") #~/.my.cnf must be present
    dbConn
  }else{
    dbConn
  }
}

#this is a hack-around version of dbQuoteString that does the exact same thing
#but does not require an instance of DBIConnection to be passed along.  This
#gets around the issue where we have to SQL quote before a DB connection is
#established
.quoteForMySQL <- function(x){
  x <- gsub("'", "''", x, fixed=TRUE)
  str <- paste("'", encodeString(x), "'", sep="")
  str[is.na(x)] <- "NULL"
  SQL(str)
}

.disconnectFromDB <- function(dbConn, conn){ #conn is passed in from user, dbConn is the actual connection
  if(is.null(conn)){ #we made a temp connection which needs to be closed
    dbDisconnect(dbConn)
  }
}

#parse a vector of strings into a piped together list ready for SQL REGEXP
#is tolerant of MySQL's '%' wildcard which is unfortunately used pretty extensively in legacy code
#allows single distinct queries
.parseSetNames <- function(setName){
  .quoteForMySQL(paste0("^", setName, "$", collapse="|"))
}

.intSiteRetrieverQuery <- function(command, conn){
  dbConn <- .connectToDB(conn)
  res <- suppressWarnings(dbGetQuery(dbConn, command))
  .disconnectFromDB(dbConn, conn)
  res
}

getUniqueSites <- function(setName, conn=NULL){
    if (is.list(conn) && conn$sitesFromFiles == TRUE) {
        return(get_unique_sites_from_files(setName, conn))
    }
  .intSiteRetrieverQuery(paste0("SELECT sites.siteID,
                                        sites.chr,
                                        sites.strand,
                                        sites.position,
                                        samples.sampleName
                                 FROM sites, samples
                                 WHERE sites.sampleID = samples.sampleID
                                 AND samples.sampleName REGEXP ",
                                 .parseSetNames(setName),
                                "AND sites.multihitClusterID IS NULL;"), conn)
}

getMRCs <- function(setName, numberOfMRCs=3, conn=NULL){
    if (is.list(conn) && conn$sitesFromFiles == TRUE) {
        sites.metadata <- get_sites_metadata_from_files(setName, conn)
    } else { # from DB
      sites.metadata <- .intSiteRetrieverQuery(paste0("SELECT sites.siteID,
                                                             samples.refGenome,
                                                             samples.gender,
                                                             samples.sampleName
                                                      FROM sites, samples
                                                      WHERE sites.sampleID = samples.sampleID
                                                      AND samples.sampleName REGEXP ",
                                                      .parseSetNames(setName),
                                                      "AND sites.multihitClusterID IS NULL;"), conn)
    }

  sites_meta <- data.frame("siteID"=sites.metadata$siteID,
                           "gender"=tolower(sites.metadata$gender))

    stopifnot(length(unique(sites.metadata$refGenome)) == 1)
    ref_genome <- sites.metadata$refGenome[1] # all the same
  
  mrcs <- get_N_MRCs(sites_meta, get_reference_genome(ref_genome), numberOfMRCs)

  #keep output consistant across functions
  mrcs$siteID <- as.numeric(levels(mrcs$siteID))[mrcs$siteID]
  mrcs$position <- as.numeric(levels(mrcs$position))[mrcs$position]
  mrcs$strand <- as.character(mrcs$strand)
  mrcs$chr <- as.character(mrcs$chr)
  
  merge(mrcs, sites.metadata[c("siteID", "sampleName")])
}

getMultihitPositions <- function(setName, conn=NULL){
  .intSiteRetrieverQuery(paste0("SELECT sites.siteID,
                                        sites.chr,
                                        sites.strand,
                                        sites.position,
                                        samples.sampleName
                                 FROM sites, samples
                                 WHERE sites.sampleID = samples.sampleID
                                 AND samples.sampleName REGEXP ",
                                .parseSetNames(setName),
                                "AND sites.multihitClusterID IS NOT NULL;"), conn)
}

getMultihitLengths <- function(setName, conn=NULL){
  .intSiteRetrieverQuery(paste0("SELECT DISTINCT multihitlengths.multihitClusterID,
                                                 multihitlengths.length,
                                                 multihitlengths.count
                                        FROM sites, samples, multihitlengths
                                 WHERE multihitlengths.multihitClusterID = sites.multihitClusterID
                                 AND sites.sampleID = samples.sampleID
                                 AND samples.sampleName REGEXP ",
                                .parseSetNames(setName),
                                "AND sites.multihitClusterID IS NOT NULL;"), conn)
}

getUniquePCRbreaks <- function(setName, conn=NULL){
  .intSiteRetrieverQuery(paste0("SELECT pcrbreakpoints.breakpoint,
                                        pcrbreakpoints.count,
                                        sites.position AS integration,
                                        sites.siteID,
                                        sites.chr,
                                        sites.strand,
                                        samples.sampleName
                                 FROM sites, samples, pcrbreakpoints
                                 WHERE (sites.sampleID = samples.sampleID AND
                                        pcrbreakpoints.siteID = sites.siteID)
                                 AND samples.sampleName REGEXP ",
                                 .parseSetNames(setName), 
                                "AND sites.multihitClusterID IS NULL;"), conn)
}

setNameExists <- function(setName, conn=NULL){
    if (is.list(conn) && conn$sitesFromFiles == TRUE) {
        return(get_existing_sample_name_from_files(setName, conn))
    }
  res <- .intSiteRetrieverQuery(paste0("SELECT DISTINCT sampleName
                                                    FROM samples
                                                    WHERE sampleName REGEXP ", .parseSetNames(setName), ";"), conn)

  setName %in% res$sampleName
}

#only function that treats % as a wildcard rather than a literal
getSampleNamesLike <- function(setName, conn=NULL){
    if (is.list(conn) && conn$sitesFromFiles == TRUE) {
        return(get_sample_names_like_from_files(setName, conn))
    }
  parsedSetNames <- paste0(gsub("%", "(.*)", paste0("^", setName, "$")), collapse="|")

  res <- .intSiteRetrieverQuery(paste0("SELECT DISTINCT sampleName
                                                    FROM samples
                                                    WHERE sampleName REGEXP ", .quoteForMySQL(parsedSetNames), ";"), conn)

  sampleNames <- lapply(strsplit(parsedSetNames, "\\|")[[1]], function(x){
    res$sampleName[grepl(x, res$sampleName)]
  })

  data.frame("sampleNames"=unlist(sampleNames, use.names=F),
             "originalNames"=rep(setName, sapply(sampleNames, length)))
}

getRefGenome <- function(setName, conn=NULL){
    if (is.list(conn) && conn$sitesFromFiles == TRUE) {
        return(get_ref_genome_from_files(setName, conn))
    }
  .intSiteRetrieverQuery(paste0("SELECT samples.sampleName,
                                        samples.refGenome
                                 FROM samples
                                 WHERE samples.sampleName REGEXP ", .parseSetNames(setName), ";"), conn)
}

getUniqueSiteReadCounts <- function(setName, conn=NULL){
  .intSiteRetrieverQuery(paste0("SELECT samples.sampleName,
                                        SUM(pcrbreakpoints.count) AS readCount
                                 FROM sites, samples, pcrbreakpoints
                                 WHERE (sites.sampleID = samples.sampleID AND
                                        pcrbreakpoints.siteID = sites.siteID)
                                 AND samples.sampleName REGEXP ", .parseSetNames(setName),
                                "AND sites.multihitClusterID IS NULL
                                GROUP BY sites.sampleID;"), conn)
}

getUniqueSiteCounts <- function(setName, conn=NULL){
  .intSiteRetrieverQuery(paste0("SELECT samples.sampleName,
                                        COUNT(*) AS uniqueSites
                                 FROM sites, samples
                                 WHERE sites.sampleID = samples.sampleID
                                 AND samples.sampleName REGEXP ", .parseSetNames(setName),
                                "AND sites.multihitClusterID IS NULL
                                 GROUP BY sites.sampleID;"), conn)
}
