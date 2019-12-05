# we generate series.RData in shell.
load('series.RData')
head(series.accession)
head(series.description)
library(usethis)
use_data(series.accession)
