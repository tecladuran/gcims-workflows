# Load all custom GC-IMS functions from the R/ folder
files <- list.files('R', pattern = '\\.R$', full.names = TRUE)
sapply(files, source)
