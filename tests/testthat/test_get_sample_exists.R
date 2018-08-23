context("check if sites for samples exists")
source("database.R") # provide db_name

read_conn <- dplyr::src_sqlite(db_name)

sample_ref <- dplyr::data_frame(
    sampleName=c("sample1", "sample2", "NOT_THERE", "sample2", "sample3"),
    refGenome=c("hg18", "hg18", "UNKNOWN_GENOME", "hgXXX", "hgYYY")
)

test_that("check that samples are in", {
    is_exist <- setNameExists(sample_ref, read_conn)
    expect_equal(is_exist, c(TRUE, TRUE, FALSE, TRUE, TRUE))
})
