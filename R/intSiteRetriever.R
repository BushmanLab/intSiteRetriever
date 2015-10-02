.get_unique_sites <- function(sample_ref, conn) {
    sample_ref_in_db <- .get_sample_ref_in_db(sample_ref, conn)
    sites <- tbl(conn, "sites") 
    inner_join(sites, sample_ref_in_db)
}

#' for a given sample names get sites.
#' @param setName vector of sample names
#' @export
getUniqueSites <- function(sample_ref, conn) {
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
#' @param conn DB connection
#' @return df with 3 cols: sampleName, refGenome, multihitID, length
#' @export
getMultihitLengths <- function(sample_ref, conn) {
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
#' @param sample_ref df with 2 cols: sampleName and refGenome
#' @param conn DB connection
#' @export
getUniquePCRbreaks <- function(sample_ref, conn) {
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
    samples_in_db <- select(samples_in_db, sampleID, sampleName, refGenome, gender)
    inner_join(samples_in_db, sample_ref, by=c('sampleName', 'refGenome'), copy=TRUE)
}

#' do we have sample names for a given connection
#'
#' @param sample_ref: df with 2 cols: sampleName, refGenome
#' @param conn DB connection
#' @return vector of TRUE/FALSE 
#' @export
setNameExists <- function(sample_ref, conn) {
    stopifnot(.check_has_sample_ref_cols(sample_ref))
    
    sample_ref_in_db <- collect(.get_sample_table(conn))
    if (nrow(sample_ref_in_db) == 0) { # nothing is in db
        return(rep(FALSE, nrow(sample_ref))) 
    }
    ids <- paste0(sample_ref$sampleName, sample_ref$refGenome)
    ids_DB <- paste0(sample_ref_in_db$sampleName, sample_ref_in_db$refGenome)
    ids %in% ids_DB
}

#' counts
#'
#' @param setName vector of sample names
#' @param conn DB connection
#' @export
getUniqueSiteReadCounts <- function(sample_ref, conn) {
    stopifnot(.check_has_sample_ref_cols(sample_ref))
    sample_ref_sites_breakpoints <- .get_breakpoints(sample_ref, conn) 
    sample_ref_sites_breakpoints_grouped <- group_by(
        sample_ref_sites_breakpoints, sampleName, refGenome)
    collect(summarize(sample_ref_sites_breakpoints_grouped, readCount=sum(count)))
}

#' unique counts for integration sites for a given sample(with fixed genome)
#'
#' @param sample_ref df with 2 cols: sampleName, refGenome
#' @param conn DB connection
#' @export
getUniqueSiteCounts <- function(sample_ref, conn) {
    stopifnot(.check_has_sample_ref_cols(sample_ref))
    sample_ref_sites <- .get_unique_sites(sample_ref, conn)
    sample_ref_sites_grouped <- group_by(sample_ref_sites, sampleName, refGenome)
    collect(summarize(sample_ref_sites_grouped, uniqueSites=n()))
}


#' creates match random controls.
#'
#' @param sampe_ref  df with 2 cols: sampleName, refGenome 
#' @param numberOfMRCs how many controls for each site
#' @param DB connection
#' @return df with cols: siteID, position, strand, chr, sampleName, refGenome
#' @export
#'
getMRCs <- function(sample_ref, conn, numberOfMRCs=3) {
    stopifnot(.check_has_sample_ref_cols(sample_ref))
    sites <- .get_unique_sites(sample_ref, conn) 
    sites.metadata <- collect(select(sites, 
        siteID, gender, sampleName, refGenome))

    sites_meta <- data.frame("siteID"=sites.metadata$siteID,
                           "gender"=tolower(sites.metadata$gender))

    stopifnot(length(unique(sites.metadata$refGenome)) == 1)
    ref_genome <- sites.metadata$refGenome[1] # all the same
  
    mrcs <- get_N_MRCs(sites_meta, get_reference_genome(ref_genome), numberOfMRCs)
  
    merge(mrcs, sites.metadata[c("siteID", "sampleName", "refGenome")])
}
