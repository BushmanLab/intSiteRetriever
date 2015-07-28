sql_schema <- "../../integration_site_schema.sql"
db_name <- "test_database"
system(paste('sqlite3', db_name, '<', sql_schema))

dbConn <- dbConnect(RSQLite::SQLite(), db_name)

res <- dbSendQuery(dbConn, 'INSERT INTO samples VALUES (1, "sample1", "hg18", "m", "asdfas");' )
res <- dbSendQuery(dbConn, 'INSERT INTO samples VALUES (2, "sample2", "hg18", "m", "asdfas");' )
res <- dbSendQuery(dbConn, 'INSERT INTO samples VALUES (3, "sample2", "hgXXX", "m", "asdfas");' )
res <- dbSendQuery(dbConn, 'INSERT INTO samples VALUES (4, "sample3", "hgYYY", "m", "asdfas");' )

res <- dbSendQuery(dbConn, 'INSERT INTO sites VALUES (1, 1, 34, "chr1", "+");' )
res <- dbSendQuery(dbConn, 'INSERT INTO sites VALUES (2, 1, 34, "chr10", "+");' )
res <- dbSendQuery(dbConn, 'INSERT INTO sites VALUES (3, 1, 99, "chr1", "+");' )
res <- dbSendQuery(dbConn, 'INSERT INTO sites VALUES (4, 2, 11, "chr10", "+");' )
res <- dbSendQuery(dbConn, 'INSERT INTO sites VALUES (5, 2, 11, "chr10", "+");' )

res <- dbSendQuery(dbConn, 'INSERT INTO pcrbreakpoints VALUES (1, 100, 123);' )
res <- dbSendQuery(dbConn, 'INSERT INTO pcrbreakpoints VALUES (1, 50, 3);' )
res <- dbSendQuery(dbConn, 'INSERT INTO pcrbreakpoints VALUES (1, 45, 1);' )
res <- dbSendQuery(dbConn, 'INSERT INTO pcrbreakpoints VALUES (2, 150, 33);' )
res <- dbSendQuery(dbConn, 'INSERT INTO pcrbreakpoints VALUES (4, 175, 666);' )

res <- dbSendQuery(dbConn, 'INSERT INTO multihitpositions VALUES (1, 1, 1234, "chr1", "+");' )
res <- dbSendQuery(dbConn, 'INSERT INTO multihitpositions VALUES (1, 1, 1234, "chr10", "+");' )
res <- dbSendQuery(dbConn, 'INSERT INTO multihitpositions VALUES (2, 1, 9999, "chr1", "+");' )
res <- dbSendQuery(dbConn, 'INSERT INTO multihitpositions VALUES (3, 2, 1111, "chr10", "+");' )

res <- dbSendQuery(dbConn, 'INSERT INTO multihitlengths VALUES (1, 100, 123)')
res <- dbSendQuery(dbConn, 'INSERT INTO multihitlengths VALUES (1, 120, 2)')
res <- dbSendQuery(dbConn, 'INSERT INTO multihitlengths VALUES (2, 300, 345)')
res <- dbSendQuery(dbConn, 'INSERT INTO multihitlengths VALUES (2, 100, 5)')
res <- dbSendQuery(dbConn, 'INSERT INTO multihitlengths VALUES (3, 333, 3)')
dbClearResult(res)
