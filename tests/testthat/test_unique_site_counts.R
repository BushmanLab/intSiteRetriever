source("database.R") # provide db_name

read_conn <- src_sqlite(db_name)

context("Unique sites counts")

sample_ref <- data_frame(
    sampleName=c("sample1", "sample2", "sample2", "sample3"),
    refGenome=c("hg18", "hg18", "hgXXX", "hgYYY")
)

test_that("for samples not in db return nothing", {
    sample_ref <- data_frame(
        sampleName=c("NOT_IN_DB_sample1", "NOT_IN_DB_sample2"),
        refGenome=c("NOT_IN_DB_hg18", "NOT_IN_DB_hg18")
    )
    expect_equal(nrow(getUniqueSiteCounts(sample_ref,  read_conn)), 0)
})

test_that("for samples in db with sites return counts", {
    sample_count <- getUniqueSiteCounts(sample_ref, read_conn)
    expect_equal(filter(sample_count, sampleName=="sample1", refGenome=="hg18")$uniqueSites,
        3
    )
    expect_equal(filter(sample_count, sampleName=="sample2", refGenome=="hg18")$uniqueSites,
        2
    )
})

rm(read_conn)
