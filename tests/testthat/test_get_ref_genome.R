context("can get reference genome info")

sites_final_path <- Sys.glob("../../data/*/sites.final.RData")
sampleInfo_path <- "../../data/sampleInfo.tsv"

connection <- create_connection_from_files(
    sampleInfo_path, sites_final_path)

sample_names_existing <- c("pool2-4", "HIV_CTRL_noLig-1", "clone2-3")
sample_names_non_existing <- c("pool42", "beVeryDRY", "SOLID")
sample_names <- c(sample_names_existing, sample_names_non_existing)

test_that("ref genomes are correct or NA", {
    is_exist <- getRefGenome(sample_names, connection)
    expect_equal(is_exist, c(
        rep("hg18", length(sample_names_existing)),
        rep(NA, length(sample_names_non_existing))
    ))
})
