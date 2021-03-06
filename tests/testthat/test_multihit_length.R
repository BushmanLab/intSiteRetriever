source("database.R") # provide db_name

read_conn <- src_sqlite(db_name)

context("setNameExists")

test_that("can find sample that is in db", {
    sample_ref <- data_frame(
        sampleName=c("sample1"),
        refGenome=c("hg18")
    )
    expect_true(setNameExists(sample_ref, read_conn))
})

test_that("can NOT find sample that is NOT in db", {
    sample_ref <- data_frame(
        sampleName=c("sample1_NOT_THERE"),
        refGenome=c("hg18")
    )
    expect_false(setNameExists(sample_ref, read_conn))
})

test_that("can present/absent samples", {
    sample_ref <- data_frame(
        sampleName=c("sample3", "sample1_NOT_THERE"),
        refGenome=c("hgYYY", "hg18")
    )
    expect_equal(setNameExists(sample_ref, read_conn), c(TRUE, FALSE))
})

test_that("the same sample but different genomes", {
    sample_ref <- data_frame(
        sampleName=c("sample2", "sample2", "sample2"),
        refGenome=c("hg18", "hgXXX", "NOT_IN_DB")
    )
    expect_equal(setNameExists(sample_ref, read_conn), c(TRUE, TRUE, FALSE))
})

context("check multihit length")

samples <- data_frame(
    sampleName=c("sample1", "sample2"),
    refGenome=c("hg18", "hg18")
)

result <- getMultihitLengths(samples, read_conn)

test_that("all tables exists in db", {
    expect_equal(dbListTables(dbConn), 
        c("multihitlengths", "multihitpositions", "pcrbreakpoints", "samples", "sites"))
})


test_that("return dataframe with 3 columns", {
    expect_is(result, "data.frame")
    expect_named(result, c("sampleName", "refGenome", "multihitID", "length"))
})

test_that("has all samples", {
    expect_equal(unique(result$sampleName), samples$sampleName)
})

test_that("has 4 lengths for sample1", {
    expect_equal(nrow(filter(result, sampleName =="sample1", refGenome=="hg18")), 4)
})

test_that("values are correct for sample1", {
    expected <- data.frame(
        sampleName=c("sample1", "sample1", "sample1", "sample1"),
        refGenome=c("hg18", "hg18", "hg18", "hg18"),
        multihitID=c(1L, 1L, 2L, 2L),
        length=c(100L, 120L, 100L, 300L),
        stringsAsFactors=FALSE
    )
    expect_equivalent(filter(result, sampleName =="sample1", refGenome=="hg18"), expected)
})

test_that("has 1 lengths for sample2", {
    expect_equal(nrow(filter(result, sampleName =="sample2", refGenome=="hg18")), 1)
})

test_that("return nothing for non-existing sample", {
    expect_equal(nrow(getMultihitLengths(data_frame(sampleName="it does not exist",
        refGenome="hg18"), read_conn)), 0)
})
