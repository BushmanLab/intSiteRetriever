#' creates connection from RData files and sample information
#'
#' @param sites_final_path vector of file paths to sites.final.RData
#' @param sampleInfo file path to sample information(tab separated file)
#' @param ref_genome name of genome, at present the same 
#' genome is used for all the samples
#'
#' @return connection (represented as a list)
#' connection has: sites df, sample_sex df, sitesFromFiles members
#'
#' @note stop if sample names from RData do not present in sample information file
#' @note list has member sitesFromFiles that is TRUE
#' @note connection is used instead of DB connection
#' @note required columns are alias(synonym of setName, sampleNames) and gender
#' @export
create_connection_from_files <- function(sampleInfo, sites_final_path, ref_genome="hg18") {
    sample_sex <- .read_sampleInfo(sampleInfo)
    sites <- .read_sites(sites_final_path)
    sites$siteID <- 1:nrow(sites)
    .check_rdata_sample_info_consistent(sample_sex, sites)
    connection <- list(sites=sites, sitesFromFiles=TRUE, 
        sample_sex=sample_sex, ref_genome=ref_genome)
    connection
}

get_unique_sites_from_files <- function(sampleNames, connection) {
    stopifnot(connection$sitesFromFiles == TRUE)
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

get_existing_sample_name_from_files <- function(sample_names, connection) {
    known_sample_names <- unique((connection$sample_sex)$alias)
    sample_names %in% known_sample_names
}

get_ref_genome_from_files <- function(sample_names, connection) {
    sample_exist <- get_existing_sample_name_from_files(
        sample_names, connection)
    ref_genomes <- rep(NA, length(sample_exist))
    ref_genomes[sample_exist] <- connection$ref_genome
    data.frame(sampleName=sample_names, refGenome=ref_genomes, stringsAsFactors=FALSE)
}

get_sites_metadata_from_files <- function(sample_names, connection) {
    current_samples <- filter(connection$sites, sampleName %in% sample_names)
    metadata <- dplyr::select(current_samples, siteID, sampleName)
    metadata$refGenome <- connection$ref_genome
    metadata <- merge(metadata, connection$sample_sex, by.x="sampleName", by.y="alias")
    metadata <- dplyr::select(metadata, siteID, refGenome, gender, sampleName)
    metadata
}

#' for a vector of sample names find matches to regex vector.
#'
#' @param sample_name vector of full sample names(replicates)
#' @param sample_prefix regex for sample prefixes
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
    stopifnot(length(sample) == 1)
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
    stopifnot("gender" %in% names(sample_sex))
    stopifnot("alias" %in% names(sample_sex))
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

.check_rdata_sample_info_consistent <- function(samples, sites) {
    stopifnot(length(setdiff(sites$sampleName, samples$alias)) == 0)
}
