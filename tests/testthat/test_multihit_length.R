sql_schema <- "../../integration_site_schema.sql"
db_name <- "test_database"
system(paste('sqlite3', db_name, '<', sql_schema))

dbConn <- dbConnect(RSQLite::SQLite(), db_name)

dbSendQuery(dbConn, 'INSERT INTO samples VALUES (1, "sample1", "hg18", "m", "asdfas");' )
dbSendQuery(dbConn, 'INSERT INTO samples VALUES (2, "sample2", "hg18", "m", "asdfas");' )
dbSendQuery(dbConn, 'INSERT INTO samples VALUES (3, "sample2", "hgXXX", "m", "asdfas");' )
dbSendQuery(dbConn, 'INSERT INTO samples VALUES (4, "sample3", "hgYYY", "m", "asdfas");' )

dbSendQuery(dbConn, 'INSERT INTO multihitpositions VALUES (1, 1, 1234, "chr1", "+");' )
dbSendQuery(dbConn, 'INSERT INTO multihitpositions VALUES (1, 1, 1234, "chr10", "+");' )
dbSendQuery(dbConn, 'INSERT INTO multihitpositions VALUES (2, 1, 9999, "chr1", "+");' )
dbSendQuery(dbConn, 'INSERT INTO multihitpositions VALUES (3, 2, 1111, "chr10", "+");' )

dbSendQuery(dbConn, 'INSERT INTO multihitlengths VALUES (1, 100, 123)')
dbSendQuery(dbConn, 'INSERT INTO multihitlengths VALUES (1, 120, 2)')
dbSendQuery(dbConn, 'INSERT INTO multihitlengths VALUES (2, 300, 345)')
dbSendQuery(dbConn, 'INSERT INTO multihitlengths VALUES (2, 100, 5)')
dbSendQuery(dbConn, 'INSERT INTO multihitlengths VALUES (3, 333, 3)')

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

samples <- c("sample1", "sample2")

result <- getMultihitLengths(samples, dbConn)

test_that("all tables exists in db", {
    expect_equal(dbListTables(dbConn), 
        c("multihitlengths", "multihitpositions", "pcrbreakpoints", "samples", "sites"))
})


test_that("return dataframe with 3 columns", {
    expect_is(result, "data.frame")
    expect_named(result, c("sampleName", "multihitID", "length"))
})

test_that("has all samples", {
    expect_equal(unique(result$sampleName), samples)
})

test_that("has 4 lengths for sample1", {
    expect_equal(nrow(filter(result, sampleName =="sample1")), 4)
})

test_that("values are correct for sample1", {
    expected <- data.frame(
        sampleName=c("sample1", "sample1", "sample1", "sample1"),
        multihitID=c(1, 1, 2, 2),
        length=c(100, 120, 100, 300),
        stringsAsFactors=FALSE
    )
    expect_equal(filter(result, sampleName =="sample1"), expected)
})

test_that("has 1 lengths for sample2", {
    expect_equal(nrow(filter(result, sampleName =="sample2")), 1)
})

test_that("return nothing for non-existing sample", {
    expect_equal(nrow(getMultihitLengths("it does not exist", dbConn)), 0)
})
