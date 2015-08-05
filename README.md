[![Travis-CI Build Status](https://travis-ci.org/BushmanLab/intSiteRetriever.svg?branch=master)](https://travis-ci.org/BushmanLab/intSiteRetriever)

![codecov.io](http://codecov.io/github/BushmanLab/intSiteRetriever/branch.svg?branch=master)

# Get integration sites from Database or from Rdata files

For a given sample name gets integration site and matched random controls(MRCs).



# Interface

At present we provide functions that work by loading data from database
or from RData files produced by intSiteCaller.

We need to create connection for both cases. 
If we have RData files and sampleInfo.tsv connection can be established by:

```r
sites_final_path <- Sys.glob("../../data/*/sites.final.RData")
sampleInfo_path <- "../../data/sampleInfo.tsv"

connection <- create_connection_from_files(
    sampleInfo_path, sites_final_path)
```

After connection is established we can get sites by:

```r
getUniqueSites(c("pool1-1", "clone1-1"), connection)
```

Functions:

* getUniqueSites
* getMRCs
* setNameExists
* getReferenceGenome

## Generation of random positions on genome

* get_reference_genome finds genome object for UCSC name
* get_random_positions generate gender-specific positions for a genome
* get_N_MRCs generate MRCs for sites that are not in DB


# Testing of the library components

Run in the R console:

```bash
library(devtools)
devtools::test()
```

Database schema is from file `intsitesdev.sqlite`.


# Installation

From the R terminal run:
```
devtools::install_github("BushmanLab/intSiteRetriever")
```

If source is required, first get source from github, from bash:

```
git clone https://github.com/BushmanLab/intSiteRetriever.git
```

and from R terminal:

```
devtools::document()
devtools::install()
library(intSiteRetriever)
```

After that public API can be accessed.

# continuous integration 
