context("integration sites from files not DB")

sites_final_path <- Sys.glob("./data/*/sites.final.RData")
sampleInfo_path <- "./data/sampleInfo.tsv"

connection <- create_connection_from_files(
    sampleInfo_path, sites_final_path)

test_that("can create file connection, represented as a list", {
    expect_true(is.list(connection))
})

test_that("fails if cannot parse alias and gender", {
    expect_error(create_connection_from_files(
        "./data/sampleInfo_csv_but_need_tsv.csv",
        sites_final_path)
    )
})

test_that("RData samples that are not in sampleInfo are invalid", {
    sampleInfo_path_inconsistent <- "./data/smaller_sampleInfo.tsv"
    expect_error(create_connection_from_files(info, sampleInfo_path_inconsistent)) 
})

test_that("can not create file connection if files are not found", {
    expect_error(create_connection_from_files("this does not exist", sites_final_path))
    expect_error(create_connection_from_files(sampleInfo_path, c("this does not exist")))
})

test_that("file connection have sites element and confirmation", {
    expect_true("sites" %in% names(connection))
    expect_true("sitesFromFiles" %in% names(connection))
    expect_true("sample_sex" %in% names(connection))
    expect_true("ref_genome" %in% names(connection))
})

test_that("sites element has 5 columns", {
    expected <- c("siteID", "chr", "strand", "position", "sampleName")
    actual <- names(connection$sites) 
    expect_equal(length(setdiff(expected, actual)), 0)
})

create_sample_ref <- function(sampleName) {
    data_frame(
        sampleName=sampleName,
        refGenome=rep("hg18", length(sampleName))
    )
}

test_that("can get sites that are present in files", {
    expect_equal(nrow(getUniqueSites(create_sample_ref("pool1-1"), connection)), 4)
})

test_that("can get sites that are present in different files", {
    expect_equal(nrow(getUniqueSites(
        create_sample_ref(c("pool1-1", "clone1-1")), connection)), 4 + 1)
    expect_equal(nrow(getUniqueSites(
        create_sample_ref(c("pool1-1", "clone1-1", "clone4-3")), connection)), 4 + 1 + 2)
})

test_that("if sampleName is not found it is ignored", {
    expect_equal(nrow(getUniqueSites(
        create_sample_ref(c("pool1-1", "sample that does not exist")), connection)), 4
    )
})

context("MRC sites from DB")

source("database.R") # provide db_name
db_conn <- src_sqlite(db_name)
sample_ref <- data_frame(
    sampleName=c("sample1", "sample2"),
    refGenome=c("hg18", "hg18")
)

test_that("correct number of MRCs", {
    # we don't have human genome on travis becouse of memory restriction
    skip_if_not_installed("BSgenome.Hsapiens.UCSC.hg18")
    expect_equal(
        nrow(getUniqueSites(sample_ref, db_conn)) * 3,
        nrow(getMRCs(sample_ref, db_conn))
    )
})

test_that("correct columns of MRCs", {
    # we don't have human genome on travis becouse of memory restriction
    skip_if_not_installed("BSgenome.Hsapiens.UCSC.hg18")
    mrcs <- getMRCs(sample_ref, db_conn)
    expect_named(mrcs, c("siteID", "position", "strand", "chr", "sampleName", "refGenome"),
        ignore.order=TRUE)
})

