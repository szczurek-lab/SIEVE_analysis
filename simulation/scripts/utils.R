# Title     : Utils
# Objective : Provide utensil functions for R
# Created by: senbaikang
# Created on: 07.03.21

if(!"scales" %in% installed.packages()){
  install.packages("scales", repos = repository)
}
library(scales)

source(file = "magic_color.R")

load.data <- function(
  file,
  sep = "\t",
  additional_labels = NULL
) {
  data <- read.table(file = file, header = TRUE, sep = sep)

  if (!is.null(additional_labels)) {
    stopifnot(length(additional_labels[[1]]) == length(additional_labels[[2]]))

    for (i in seq_len(length(additional_labels[[1]]))) {
      data[[additional_labels[[1]][i]]] <- additional_labels[[2]][i]
    }
  }

  return(data)
}

load.rds <- function(
  file,
  additional_labels = NULL
) {
  data <- readRDS(file = file)

  if (!is.null(additional_labels)) {
    stopifnot(length(additional_labels[[1]]) == length(additional_labels[[2]]))

    for (i in seq_len(length(additional_labels[[1]]))) {
      data[[additional_labels[[1]][i]]] <- additional_labels[[2]][i]
    }
  }

  return(data)
}

scientific <- function(vals){
  sapply(
    vals,
    function(x) {
      if (x == 0)
        return("0")

      x <- ifelse(
        grepl("[eE]", x),
        x,
        label_scientific()(x)
      )

      if (grepl("^1[eE]", x)) {
        x <- gsub("^1[eE]", "10^", x)
      } else {
        x <- gsub("[eE]", "%*%10^", x)
      }

      parse(text = gsub(
        "[+]",
        "",
        x
      ))
    }
  )
}

define.legend <- function(
  use.common.legend,
  common.legend,
  original.legend
) {
  if (use.common.legend) {
    common.legend
  } else {
    original.legend
  }
}
