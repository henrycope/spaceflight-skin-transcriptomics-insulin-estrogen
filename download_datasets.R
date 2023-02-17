library(jsonlite)
library(purrr)

# Define API endpoint and datasets to be downloaded
api_url <- "https://osdr.nasa.gov/genelab/data/glds/files/"
dataset_ids <- c("238", "239", "240", "241", "254")

# Check if main data directory exists, and create it if it doesn't
if (!file.exists("data")) {
  dir.create("data")
}

# Check if each dataset subdirectory exist, and create them if none of them do
dataset_dirs <- file.path("data", dataset_ids)

if (all(file.exists(dataset_dirs))) {
  message("Data appears to already be downloaded, please delete existing files before redownloading.")
  q() # exit R
}

walk(dataset_dirs, dir.create, recursive = TRUE)

# Make API request
request_url <- paste0(api_url, paste0(dataset_ids, collapse = ","))
response <- fromJSON(request_url)

# Process API response, to get urls for raw counts and metadata files to download
file_urls <- map(response$studies, ~ .x$study_files$remote_url[grepl("*genes.results|*metadata_GLDS*", .x$study_files$remote_url)])

# Download files
walk2(file_urls, dataset_dirs, ~ walk2(.x, .y, ~ {
  file_local_path <- file.path(.y, sub(".*file=([^=]*$)", "\\1", .x)) # regex to get just file name
  
  download.file(file.path("https://osdr.nasa.gov/", .x), destfile = file_local_path)
}))