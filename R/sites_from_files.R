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
