library(masstrixR)
library(tidyverse)

# get ymdb example
ymdb_example_file <- system.file("testdata", "ymdb_example.txt", package="masstrixR")

# read ymdb
ymdb_compound_list <- read_tsv(ymdb_example_file, col_types = cols(exactmass = col_double()))

# adducts used for DB generation
adducts <- c("[M+H]+", "[M+Na]+")

# create compound list for DB creation
ymdb_new_compound_list <- prepareCompoundList(ymdb_compound_list, adductList = adducts)

# check if compound list is valid and create SQLite DB
if(validateCompoundList(ymdb_new_compound_list)) {
  dbFileName <- createDb(ymdb_new_compound_list, "test")
}
