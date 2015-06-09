library(dplyr)

#' creates connection from RData files
#'
#' @param sites_final_path vector of file paths to sites.final.RData
#' @param sampleInfo file path to sample information
#' required columns are alias(synonym of setName, sampleNames) and gender
#' @return connection (represented as a list)
#' connection has: sites df, sample_sex df, sitesFromFiles members.
#'
#' @note list has member sitesFromFiles that is TRUE
#' @note connection is used instead of DB connection
create_connection_from_files <- function(sampleInfo, sites_final_path) {
    sample_sex <- .read_sampleInfo(sampleInfo)
    sites <- .read_sites(sites_final_path)
    connection <- list(sites=sites, sitesFromFiles=TRUE, sample_sex=sample_sex)
    connection
}

get_unique_sites_from_files <- function(sampleNames, connection) {
    stopifnot(connection$sitesFromFiles == TRUE)
    stopifnot(all(sampleNames %in% (connection$sites)$sampleName))
    filter(connection$sites, sampleName %in% sampleNames)
}

get_sample_names_like_from_files <- function(sample_name_prefix, connection) {
    known_sample_names <- unique((connection$sample_sex)$alias)
    sample_name_regex <- gsub("%", "(.*)", sample_name_prefix)

    matches <- .sample_name_prefix_match(known_sample_names, sample_name_regex)
    matches <- gsub("(.*)", "%", matches, fixed=TRUE) # back to input format

    replicates <- data.frame(sampleNames=known_sample_names,
        originalNames=matches)
    filter(replicates,  ! is.na(originalNames))
}

#' for a vector of sample names find matches to regex vector.
#'
#' @param sample_name vector of full sample names(replicates)
#' @sample_prefix regex for sample prefixes
#' @return vector of the same length as sample_name
#' the value is either sample_prefix or NA if no hit is found.
.sample_name_prefix_match <- function(sample_name, sample_prefix) {
    sapply(sample_name, function(known_sample) {
        .find_match(known_sample, sample_prefix)
    }) 
}


#' for one sample and all possible prefixes find if
#' there is a match
.find_match <- function(sample, possible_prefix) {
    stopifnot(length(sample) != 1)
    matches <- sapply(possible_prefix, function(prefix) {
        if(grepl(prefix, sample)) {
            prefix
        } else {
            NA
        }
    })
    match_index <- which( ! is.na(matches))
    # replicate should only match one prefix
    stopifnot(length(match_index) <= 1)
    if (length(match_index) == 0) {
        as.character(NA) # no hit
    } else {
        matches[match_index] # single hit
    }
}

.read_sites <- function(sites_final_path) {
    do.call(rbind, lapply(sites_final_path, .get_sites_from_rdata))
}

.get_sites_from_rdata <- function(sites_final_path) {
    if ( ! file.exists(sites_final_path)) {
        stop(paste("Cannot find sites.final.RData file at:",  
            sites_final_path))
    }
    load(sites_final_path)
    sites <- data.frame(
        chr=as.character(seqnames(sites.final)),
        strand=as.character(strand(sites.final)),
        position=start(flank(sites.final, -1, start=TRUE)),
        sampleName=mcols(sites.final)$sampleName,
        stringsAsFactors=FALSE
    )
}

#TODO: extract sample reading and checking into mini-package
#TODO: see also intSiteCaller read_sample_files.R
.read_sampleInfo <- function(sampleInfo_file) {
    sample_sex <- read.delim(sampleInfo_file, stringsAsFactors=FALSE)
    .check_sex_is_correct(sample_sex$gender)
    .check_sampleName_is_correct(sample_sex$alias)
    sample_sex
}

.valid_sex <- c('m', 'f')

.check_sex_is_correct <- function(gender) {
    stopifnot(all(gender %in% .valid_sex))
}

.check_sampleName_is_correct <- function(samples) {
    stopifnot( ! (is.null(samples)))
    stopifnot(length(samples) == length(unique(samples)))
}
