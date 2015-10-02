[![Travis-CI Build Status](https://travis-ci.org/BushmanLab/intSiteRetriever.svg?branch=master)](https://travis-ci.org/BushmanLab/intSiteRetriever)

[![codecov.io](http://codecov.io/github/BushmanLab/intSiteRetriever/coverage.svg?branch=master)](http://codecov.io/github/BushmanLab/intSiteRetriever?branch=master)

# Get integration sites from Database

For a given sample name gets integration site and matched random controls(MRCs).

# Installation

From the R terminal  or RStudio console run:
```
devtools::install_github("BushmanLab/intSiteRetriever")
```

If source is required, first get source from github(from bash):
```
git clone https://github.com/BushmanLab/intSiteRetriever.git
```
and from R terminal:
```{r}
devtools::document()
devtools::install()
library(intSiteRetriever)
```

# Connection to Database

For MySQL database config file should be in home directory and called .my.cnf.

.my.cnf file format is:
```
[GROUP_NAME]
user=YYYYYYY
password=XXXXXX
host=microbYYYY.med.upenn.edu
port=3309
database=intsites_miseq
```

We can establish connection in R:
```{r}
dbConn <- dbConnect(MySQL(), group='GROUP_NAME')
info <- dbGetInfo(dbConn)
dbConn <- src_sql("mysql", dbConn, info = info)
```
The group should be the same as group in `.my.cnf` file.
Note that by default connection will expire in about 5 minutes for MySQL.

#Interface

After connection is established we can get sites:

```{r}
sample_ref <- data_frame(
    sampleName=c("sample1", "sample2"),
    refGenome=c("hg18", "hg18")
)
sites <- getUniqueSites(sample_ref, dbConn)
```

Public functions are:

* setNameExists
* getUniqueSites
* getUniqueSiteCounts
* getUniqueSiteReadCounts
* getUniquePCRbreaks
* getMultihitLengths
* getMRCs
* get_N_MRCs
* get_random_positions
* get_reference_genome

Documentation for functions:
```{r}
?getUniqueSites
```

#Dependencies

Dependencies can be loaded with:
```
library(GenomicRanges)
library(BSgenome)
library(RMySQL)
library(dplyr)
```

# Testing of the library components

Run in the R console:
```bash
library(devtools)
devtools::test()
```

# Schema

Database schema is from file `integration_site_schema.sql`.

# Continuous Integration 

Travis is used for testing on each commit and Codecov used for code coverage metrics.
