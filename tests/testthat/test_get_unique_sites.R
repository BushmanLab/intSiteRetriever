context("integration sites from files not DB")

sites_final_path <- Sys.glob("../../data/*/sites.final.RData")
sampleInfo_path <- "../../data/sampleInfo.tsv"

connection <- create_connection_from_files(
    sampleInfo_path, sites_final_path)

test_that("can create file connection, represented as a list", {
    expect_true(is.list(connection))
})

test_that("RData samples that are not in sampleInfo are invalid", {
    sampleInfo_path_inconsistent <- "../../data/smaller_sampleInfo.tsv"
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
})

test_that("can get sites that are present in files", {
    expect_equal(nrow(getUniqueSites("pool1-1", connection)), 4)
})

test_that("can get sites that are present in different files", {
    expect_equal(nrow(getUniqueSites(
        c("pool1-1", "clone1-1"), connection)), 4 + 1)
    expect_equal(nrow(getUniqueSites(
        c("pool1-1", "clone1-1", "clone4-3"), connection)), 4 + 1 + 2)
})

test_that("fail if sampleName is not found", {
    expect_error(getUniqueSites(
        c("pool1-1", "sample that does not exist"), connection))
})
