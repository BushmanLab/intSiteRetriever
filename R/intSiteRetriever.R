.get_unique_sites <- function(sample_ref, conn) {
    sample_ref_in_db <- .get_sample_ref_in_db(sample_ref, conn)
    sites <- tbl(conn, "sites") 
    inner_join(sites, sample_ref_in_db)
}

#' for a given sample names get sites.
#' @param setName vector of sample names
#' @export
getUniqueSites <- function(sample_ref, conn) {
    if (is.list(conn) && "sitesFromFiles" %in% names(conn) && conn$sitesFromFiles == TRUE) {
        refGenome <- unique(sample_ref$refGenome)
        stopifnot(length(refGenome) == 1) # file connection is for 1 genome at present
        stopifnot(refGenome == conn$ref_genome)
        return(get_unique_sites_from_files(sample_ref$sampleName, conn))
    }
    stopifnot(.check_has_sample_ref_cols(sample_ref))
    sites <- .get_unique_sites(sample_ref, conn)
    collect( select(sites, 
        siteID, chr, strand, position, sampleName, refGenome)
    )
}

.get_multihitpositions <- function(sample_ref, conn) {
    sample_ref_in_db <- .get_sample_ref_in_db(sample_ref, conn)
    multihitpositions <- tbl(conn, "multihitpositions") 
    inner_join(multihitpositions, sample_ref_in_db, by="sampleID")
}

#' lengths distributions for multihits
#'
#' @param sampleName vector of sample names
#' @param conn connection: DB or File connection
#' @return df with 3 cols: sampleName, refGenome, multihitID, length
#' @export
getMultihitLengths <- function(sample_ref, conn) {
    if (is.list(conn) && "sitesFromFiles" %in% names(conn) && conn$sitesFromFiles == TRUE) {
        stop("getMultihitLengths is not implemented for file connection.")
    }
    stopifnot(.check_has_sample_ref_cols(sample_ref))
    samples_multihitpositions <- .get_multihitpositions(sample_ref, conn)
    multihit_lengths <- tbl(conn, "multihitlengths")
    collect(distinct(select(inner_join(samples_multihitpositions, multihit_lengths, by="multihitID"),
        sampleName, refGenome, multihitID, length)))
}

.get_breakpoints <- function(sample_ref, conn) {
    sample_ref_sites <- .get_unique_sites(sample_ref, conn)
    breakpoints <- tbl(conn, "pcrbreakpoints") 
    inner_join(sample_ref_sites, breakpoints) 
}

#' breakpoints
#'
#' @param setName vector of sample names
#' @param conn connection: DB or File connection
#' @export
getUniquePCRbreaks <- function(setName, conn) {
    if (is.list(conn) && "sitesFromFiles" %in% names(conn) && conn$sitesFromFiles == TRUE) {
        stop("getUniquePCRbreaks is not implemented for file connection.")
    }
    breakpoints <- .get_breakpoints(sample_ref, conn)
    collect(select(breakpoints,
        breakpoint, count, position, siteID, chr, strand, sampleName, refGenome)
    )
# column named kept as in DB ...sites.position AS integration...
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

#' counts
#'
#' @param setName vector of sample names
#' @param conn connection: DB or File connection
#' @export
getUniqueSiteReadCounts <- function(sample_ref, conn) {
    if (is.list(conn) && "sitesFromFiles" %in% names(conn) && conn$sitesFromFiles == TRUE) {
        stop("getUniqueSiteReadCounts does not implemented for file connection")
    }
    stopifnot(.check_has_sample_ref_cols(sample_ref))
    sample_ref_sites_breakpoints <- .get_breakpoints(sample_ref, conn) 
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

# TODO: functions below are not adjusted to PK as sampleName,refGenome

#' creates match random controls.
#' @param setName vector of sample names
#' @param numberOfMRCs how many controls for each site
#' @param conn connection: DB or File connection
#'
getMRCs <- function(setName, numberOfMRCs=3, conn=NULL){
    if (is.list(conn) && conn$sitesFromFiles == TRUE) {
        sites.metadata <- get_sites_metadata_from_files(setName, conn)
        sites.metadata
    } else { # from DB
        stop("broken")
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

#' find replicates
#'
#' @param setName vector of sample names
#' @param conn connection: DB or File connection
#' @note only function that treats % as a wildcard rather than a literal
#'
getSampleNamesLike <- function(setName, conn=NULL){
    if (is.list(conn) && conn$sitesFromFiles == TRUE) {
        return(get_sample_names_like_from_files(setName, conn))
    }
    stop("broken")
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
#'
getRefGenome <- function(setName, conn=NULL){
    if (is.list(conn) && conn$sitesFromFiles == TRUE) {
        return(get_ref_genome_from_files(setName, conn))
    }
    stop("broken")
  .intSiteRetrieverQuery(paste0("SELECT samples.sampleName,
                                        samples.refGenome
                                 FROM samples
                                 WHERE samples.sampleName REGEXP ", .parseSetNames(setName), ";"), conn)
}
