context("for a given prefix name finds all the replicates")

sites_final_path <- Sys.glob("./data/*/sites.final.RData")
sampleInfo_path <- "./data/sampleInfo.tsv"

connection <- create_connection_from_files(
    sampleInfo_path, sites_final_path)

setName_prefix <- c("HIV_CTRL_noLig%", "pool3%")

replicates <- getSampleNamesLike(setName_prefix, connection)

test_that("returns correct datastructure", {
    expect_is(replicates, "data.frame")
    expect_named(replicates, c("sampleNames", "originalNames"), ignore.order=TRUE)
})

test_that("find all replicate sampleNames counts for prefix", {
    N_REPLICATES <- 4
    expect_equal(nrow(replicates), N_REPLICATES*length(setName_prefix))
    expect_equal(sum(
        replicates$originalNames == setName_prefix[1]), N_REPLICATES)
    expect_equal(sum(
        replicates$originalNames == setName_prefix[2]), N_REPLICATES)
})

test_that("find all replicate sampleNames values for prefix", {
    N_REPLICATES <- 4
    expect_sampleNames <- c("pool3-1", "pool3-2", "pool3-3", "pool3-4")
    pool <- filter(replicates, originalNames == "pool3%")
    expect_equal(as.character(pool$sampleNames), expect_sampleNames)
})
