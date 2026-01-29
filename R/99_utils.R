print_tree <- function(path = ".", prefix = "") {
  files <- list.files(path, all.files = TRUE, full.names = TRUE)
  files <- files[!basename(files) %in% c(".", "..", ".git", ".Rproj.user")]
  
  for (i in seq_along(files)) {
    cat(
      prefix,
      if (i == length(files)) "└── " else "├── ",
      basename(files[i]),
      "\n",
      sep = ""
    )
    if (dir.exists(files[i])) {
      print_tree(
        files[i],
        paste0(prefix, if (i == length(files)) "    " else "│   ")
      )
    }
  }
}
