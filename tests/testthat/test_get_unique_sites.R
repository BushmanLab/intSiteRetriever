context("integration sites")
source("database.R") # provide db_name

read_conn <- src_sqlite(db_name)

sample_ref <- data_frame(
    sampleName=c("sample1", "sample2", "sample2", "sample3"),
    refGenome=c("hg18", "hg18", "hgXXX", "hgYYY")
)

test_that("sites element has 5 columns", {
    expected <- c("siteID", "chr", "strand", 
          "position", "sampleName", 'refGenome')
    actual <- getUniqueSites(sample_ref, read_conn)
    expect_named(actual, expected, ignore.order=TRUE)
})

sample_ref <- data_frame(
    sampleName=c("sample1"),
    refGenome=c("hg18")
)

test_that("can get sites that are present in files", {
    expect_equal(nrow(getUniqueSites(sample_ref, read_conn)), 3)
})

sample_ref <- data_frame(
    sampleName=c("sample1", "XXX"),
    refGenome=c("hg18", "YYY")
)

test_that("if sampleName is not found it is ignored", {
    expect_equal(nrow(getUniqueSites(sample_ref, read_conn)), 3)
})

context("MRC sites")

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

