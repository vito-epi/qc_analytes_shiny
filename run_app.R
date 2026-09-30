tmp <- tempfile(fileext = ".R")

download.file(
  "https://raw.githubusercontent.com/vito-epi/qc_analytes_shiny/main/app.R",
  destfile = tmp,
  mode = "wb"
)

source(tmp)