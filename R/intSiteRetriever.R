#parse a vector of strings into a piped together list ready for SQL REGEXP
#is tolerant of MySQL's '%' wildcard which is unfortunately used pretty extensively in legacy code
#allows single distinct queries
.parseSetNames <- function(setName){
    stop("obsolete")
}

.intSiteRetrieverQuery <- function(command, conn){
    stop("obsolete")
}

#' for a given sample names get sites.
#' @param setName vector of sample names
#' @export
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
                                ";"), conn)
}

#' creates match random controls.
#' @param setName vector of sample names
#' @param numberOfMRCs how many controls for each site
#' @param conn connection: DB or File connection
#' @export
getMRCs <- function(setName, numberOfMRCs=3, conn=NULL){
    if (is.list(conn) && conn$sitesFromFiles == TRUE) {
        sites.metadata <- get_sites_metadata_from_files(setName, conn)
        sites.metadata
    } else { # from DB
        sites.metadata <- .intSiteRetrieverQuery(paste0("SELECT sites.siteID,
                                                             samples.refGenome,
                                                             samples.gender,
                                                             samples.sampleName
                                                      FROM sites, samples
                                                      WHERE sites.sampleID = samples.sampleID
                                                      AND samples.sampleName REGEXP ",
                                                      .parseSetNames(setName),
                                                      ";"), conn)
    }

    sites_meta <- data.frame("siteID"=sites.metadata$siteID,
                           "gender"=tolower(sites.metadata$gender))

    stopifnot(length(unique(sites.metadata$refGenome)) == 1)
    ref_genome <- sites.metadata$refGenome[1] # all the same
  
    mrcs <- get_N_MRCs(sites_meta, get_reference_genome(ref_genome), numberOfMRCs)

    #keep output consistant across functions
    mrcs$siteID <- as.numeric(mrcs$siteID)
    mrcs$position <- as.numeric(mrcs$position)
    mrcs$strand <- as.character(mrcs$strand)
    mrcs$chr <- as.character(mrcs$chr)
  
    merge(mrcs, sites.metadata[c("siteID", "sampleName")])
}

#' multihits
#'
#' @param setName vector of sample names
#' @param conn connection: DB or File connection
#' @export
getMultihitPositions <- function(setName, conn=NULL){
    stop()
    .intSiteRetrieverQuery(paste0("SELECT sites.siteID,
                                        sites.chr,
                                        sites.strand,
                                        sites.position,
                                        samples.sampleName
                                 FROM sites, samples
                                 WHERE sites.sampleID = samples.sampleID
                                 AND samples.sampleName REGEXP ",
                                .parseSetNames(setName),
                                ";"), conn)
}

#' from vector of sample names to sql string: ( sample, ..., sampleX )
.sampleName_sql <- function(sampleName, conn) {
    samples_sql <- dbQuoteString(conn, sampleName)
    samples_sql <- paste(samples_sql, collapse=", ")
    paste("(", samples_sql, ")")
}

#' lengths distributions for multihits
#'
#' @param sampleName vector of sample names
#' @param conn connection: DB or File connection
#' @return df with 3 cols: sampleName, multihitID, length
#' @export
getMultihitLengths <- function(sampleName, conn=NULL){
    query <- paste( 
        "SELECT DISTINCT samples.sampleName, 
                multihitlengths.multihitID,
                multihitlengths.length 
        FROM samples JOIN multihitpositions
        ON samples.sampleID = multihitpositions.sampleID
        JOIN multihitlengths
        ON multihitpositions.multihitID = multihitlengths.multihitID ",
        "WHERE sampleName in ", 
        .sampleName_sql(sampleName, conn), 
        ";"
    )
    dbGetQuery(conn, query)
}

#' breakpoints
#'
#' @param setName vector of sample names
#' @param conn connection: DB or File connection
#' @export
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
                                ";"), conn)
}

.check_has_sample_ref_cols <- function(sample_ref) {
    return (all(c("sampleName", "refGenome") %in% names(sample_ref)))
}

.get_sample_table <- function(conn) {
    samples_in_db <- tbl(conn, "samples") 
    select(samples_in_db, sampleID, sampleName, refGenome)
}

.get_sample_ref_in_db <- function(sample_ref, conn) {
    samples_in_db <- tbl(conn, "samples") 
    samples_in_db <- select(samples_in_db, sampleID, sampleName, refGenome)
    inner_join(samples_in_db, sample_ref, by=c('sampleName', 'refGenome'), copy=TRUE)
}

#' do we have sample names for a given connection
#'
#' @param sample_ref: df with 2 cols: sampleName, refGenome
#' @param conn connection: DB or File connection
#' @return vector of TRUE/FALSE 
#' @export
setNameExists <- function(sample_ref, conn) {
    sample_ref
    if (is.list(conn) && "sitesFromFiles" %in% names(conn) && conn$sitesFromFiles == TRUE) {
        refGenome <- unique(sample_ref$refGenome)
        stopifnot(length(refGenome) == 1) # file connection is for 1 genome at present
        stopifnot(refGenome == conn$ref_genome)
        return(get_existing_sample_name_from_files(sample_ref$sampleName, conn))
    }
    stopifnot(.check_has_sample_ref_cols(sample_ref))
    
    sample_ref_in_db <- collect(.get_sample_table(conn))
    if (nrow(sample_ref_in_db) == 0) { # nothing is in db
        return(rep(FALSE, nrow(sample_ref))) 
    }
    sample_ref_in_db$in_db <- TRUE # to distuinguish after merging
    in_db <- merge(sample_ref, sample_ref_in_db, all.x=TRUE, sort=FALSE)$in_db
    in_db[is.na(in_db)] <- FALSE
    in_db
}

#' find replicates
#'
#' @param setName vector of sample names
#' @param conn connection: DB or File connection
#' @note only function that treats % as a wildcard rather than a literal
#' @export
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

#' name of the reference genome used
#'
#' @param setName vector of sample names
#' @param conn connection: DB or File connection
#' return  df with 2 cols: sampleName and refGenome
#' @export
getRefGenome <- function(setName, conn=NULL){
    if (is.list(conn) && conn$sitesFromFiles == TRUE) {
        return(get_ref_genome_from_files(setName, conn))
    }
  .intSiteRetrieverQuery(paste0("SELECT samples.sampleName,
                                        samples.refGenome
                                 FROM samples
                                 WHERE samples.sampleName REGEXP ", .parseSetNames(setName), ";"), conn)
}

.get_unique_sites <- function(sample_ref, conn) {
    sample_ref_in_db <- .get_sample_ref_in_db(sample_ref, conn)
    sites <- tbl(conn, "sites") 
    inner_join(sites, sample_ref_in_db)
}

#' counts
#'
#' @param setName vector of sample names
#' @param conn connection: DB or File connection
#' @export
getUniqueSiteReadCounts <- function(sample_ref, conn) {
    stopifnot(.check_has_sample_ref_cols(sample_ref))
    sample_ref_sites <- .get_unique_sites(sample_ref, conn)
    breakpoints <- tbl(conn, "pcrbreakpoints") 
    sample_ref_sites_breakpoints <- inner_join(sample_ref_sites, breakpoints) 
    sample_ref_sites_breakpoints_grouped <- group_by(
        sample_ref_sites_breakpoints, sampleName, refGenome)
    collect(summarize(sample_ref_sites_breakpoints_grouped, readCount=sum(count)))
}

#' unique counts for integration sites for a given sample(with fixed genome)
#'
#' @param sample_ref df with 2 cols: sampleName, refGenome
#' @param conn connection: DB or File connection
#' @export
getUniqueSiteCounts <- function(sample_ref, conn) {
    stopifnot(.check_has_sample_ref_cols(sample_ref))
    sample_ref_sites <- .get_unique_sites(sample_ref, conn)
    sample_ref_sites_grouped <- group_by(sample_ref_sites, sampleName, refGenome)
    collect(summarize(sample_ref_sites_grouped, uniqueSites=n()))
}
