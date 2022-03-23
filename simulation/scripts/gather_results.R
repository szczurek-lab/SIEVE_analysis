# Title     : Gather results from simulations
# Created by: senbaikang
# Created on: 07.03.21

repository <- "https://stat.ethz.ch/CRAN/"

if(!"dplyr" %in% installed.packages()){
  install.packages("dplyr", repos = repository)
}
library(dplyr)

if(!"ggplot2" %in% installed.packages()){
  install.packages("ggplot2", repos = repository)
}
library(ggplot2)

if(!"lemon" %in% installed.packages()){
  install.packages("lemon", repos = repository)
}
library(lemon)

if(!"ggpubr" %in% installed.packages()){
  install.packages("ggpubr", repos = repository)
}
library(ggpubr)


source(file = "utils.R")


###################################################
#                    Constants                    #
###################################################

true.eff.seq.err.rate <- 1.999E-3
true.ado.rate <- 0.16334
true.wildtype.overdispersion <- 100.0
true.alt.overdispersion <- 2.5

isa.mu.rate.threshold <- 3.0E-6

dot.size <- 0.5 / .pt

use.common.legend <- TRUE

strip.text <- element_text(size = 6)

output.prefix <- "./"


###################################################
#                    Variables                    #
###################################################

tree.info <- data.frame()
param.info <- data.frame()
var.info <- data.frame()


###################################################
#            Import simulation results            #
###################################################

new.col.names <- c("simulated_data", "mutation_rate")

# simulation 2
tree.info.2 <- load.data(
  file = "../simulated_data_02/trees_info_updated.tsv",
  data.label.name = new.col.names[1],
  data.label = 2L,
  mu.rate.name = new.col.names[2],
  mu.rate = 1.0E-6
)
site.info.2 <- load.data(
  file = "../simulated_data_02/sieve_sites_info.tsv",
  data.label.name = new.col.names[1],
  data.label = 2L,
  mu.rate.name = new.col.names[2],
  mu.rate = 1.0E-6
)
param.info.2 <- load.data(
  file = "../simulated_data_02/params_info.tsv",
  data.label.name = new.col.names[1],
  data.label = 2L,
  mu.rate.name = new.col.names[2],
  mu.rate = 1.0E-6
)
var.info.2 <- load.data(
  file = "../simulated_data_02/variants_info.tsv",
  data.label.name = new.col.names[1],
  data.label = 2L,
  mu.rate.name = new.col.names[2],
  mu.rate = 1.0E-6
)
ado.info.2 <- load.data(
  file = "../simulated_data_02/ado_info.tsv",
  data.label.name = new.col.names[1],
  data.label = 2L,
  mu.rate.name = new.col.names[2],
  mu.rate = 1.0E-6
)
allelic.info.2 <- load.rds(
  file = "../simulated_data_02/allelic_info.rds",
  data.label.name = new.col.names[1],
  data.label = 2L,
  mu.rate.name = new.col.names[2],
  mu.rate = 1.0E-6
)
parallel.mut.2 <- load.data(
  file = "../simulated_data_02/summary.tsv",
  sep = " ",
  data.label.name = new.col.names[1],
  data.label = 2L,
  mu.rate.name = new.col.names[2],
  mu.rate = 1.0E-6
)

# simulation 3
tree.info.3 <- load.data(
  file = "../simulated_data_03/trees_info_updated.tsv",
  data.label.name = new.col.names[1],
  data.label = 3L,
  mu.rate.name = new.col.names[2],
  mu.rate = 3.0E-5
)
site.info.3 <- load.data(
  file = "../simulated_data_03/sieve_sites_info.tsv",
  data.label.name = new.col.names[1],
  data.label = 3L,
  mu.rate.name = new.col.names[2],
  mu.rate = 3.0E-5
)
param.info.3 <- load.data(
  file = "../simulated_data_03/params_info.tsv",
  data.label.name = new.col.names[1],
  data.label = 3L,
  mu.rate.name = new.col.names[2],
  mu.rate = 3.0E-5
)
var.info.3 <- load.data(
  file = "../simulated_data_03/variants_info.tsv",
  data.label.name = new.col.names[1],
  data.label = 3L,
  mu.rate.name = new.col.names[2],
  mu.rate = 3.0E-5
)
ado.info.3 <- load.data(
  file = "../simulated_data_03/ado_info.tsv",
  data.label.name = new.col.names[1],
  data.label = 3L,
  mu.rate.name = new.col.names[2],
  mu.rate = 3.0E-5
)
allelic.info.3 <- load.rds(
  file = "../simulated_data_03/allelic_info.rds",
  data.label.name = new.col.names[1],
  data.label = 3L,
  mu.rate.name = new.col.names[2],
  mu.rate = 3.0E-5
)
parallel.mut.3 <- load.data(
  file = "../simulated_data_03/summary.tsv",
  sep = " ",
  data.label.name = new.col.names[1],
  data.label = 3L,
  mu.rate.name = new.col.names[2],
  mu.rate = 3.0E-5
)

# simulation 4
tree.info.4 <- load.data(
  file = "../simulated_data_04/trees_info_updated.tsv",
  data.label.name = new.col.names[1],
  data.label = 4L,
  mu.rate.name = new.col.names[2],
  mu.rate = 8.0E-6
)
site.info.4 <- load.data(
  file = "../simulated_data_04/sieve_sites_info.tsv",
  data.label.name = new.col.names[1],
  data.label = 4L,
  mu.rate.name = new.col.names[2],
  mu.rate = 8.0E-6
)
param.info.4 <- load.data(
  file = "../simulated_data_04/params_info.tsv",
  data.label.name = new.col.names[1],
  data.label = 4L,
  mu.rate.name = new.col.names[2],
  mu.rate = 8.0E-6
)
var.info.4 <- load.data(
  file = "../simulated_data_04/variants_info.tsv",
  data.label.name = new.col.names[1],
  data.label = 4L,
  mu.rate.name = new.col.names[2],
  mu.rate = 8.0E-6
)
ado.info.4 <- load.data(
  file = "../simulated_data_04/ado_info.tsv",
  data.label.name = new.col.names[1],
  data.label = 4L,
  mu.rate.name = new.col.names[2],
  mu.rate = 8.0E-6
)
allelic.info.4 <- load.rds(
  file = "../simulated_data_04/allelic_info.rds",
  data.label.name = new.col.names[1],
  data.label = 4L,
  mu.rate.name = new.col.names[2],
  mu.rate = 8.0E-6
)
parallel.mut.4 <- load.data(
  file = "../simulated_data_04/summary.tsv",
  sep = " ",
  data.label.name = new.col.names[1],
  data.label = 4L,
  mu.rate.name = new.col.names[2],
  mu.rate = 8.0E-6
)

# simulation 11
tree.info.11 <- load.data(
  file = "../simulated_data_11/trees_info_updated.tsv",
  data.label.name = new.col.names[1],
  data.label = 11L,
  mu.rate.name = new.col.names[2],
  mu.rate = 1.0E-6
)
site.info.11 <- load.data(
  file = "../simulated_data_11/sieve_sites_info.tsv",
  data.label.name = new.col.names[1],
  data.label = 11L,
  mu.rate.name = new.col.names[2],
  mu.rate = 1.0E-6
)
param.info.11 <- load.data(
  file = "../simulated_data_11/params_info.tsv",
  data.label.name = new.col.names[1],
  data.label = 11L,
  mu.rate.name = new.col.names[2],
  mu.rate = 1.0E-6
)
var.info.11 <- load.data(
  file = "../simulated_data_11/variants_info.tsv",
  data.label.name = new.col.names[1],
  data.label = 11L,
  mu.rate.name = new.col.names[2],
  mu.rate = 1.0E-6
)
ado.info.11 <- load.data(
  file = "../simulated_data_11/ado_info.tsv",
  data.label.name = new.col.names[1],
  data.label = 11L,
  mu.rate.name = new.col.names[2],
  mu.rate = 1.0E-6
)
allelic.info.11 <- load.rds(
  file = "../simulated_data_11/allelic_info.rds",
  data.label.name = new.col.names[1],
  data.label = 11L,
  mu.rate.name = new.col.names[2],
  mu.rate = 1.0E-6
)
parallel.mut.11 <- load.data(
  file = "../simulated_data_11/summary.tsv",
  sep = " ",
  data.label.name = new.col.names[1],
  data.label = 11L,
  mu.rate.name = new.col.names[2],
  mu.rate = 1.0E-6
)

# simulation 12
tree.info.12 <- load.data(
  file = "../simulated_data_12/trees_info_updated.tsv",
  data.label.name = new.col.names[1],
  data.label = 12L,
  mu.rate.name = new.col.names[2],
  mu.rate = 8.0E-6
)
site.info.12 <- load.data(
  file = "../simulated_data_12/sieve_sites_info.tsv",
  data.label.name = new.col.names[1],
  data.label = 12L,
  mu.rate.name = new.col.names[2],
  mu.rate = 8.0E-6
)
param.info.12 <- load.data(
  file = "../simulated_data_12/params_info.tsv",
  data.label.name = new.col.names[1],
  data.label = 12L,
  mu.rate.name = new.col.names[2],
  mu.rate = 8.0E-6
)
var.info.12 <- load.data(
  file = "../simulated_data_12/variants_info.tsv",
  data.label.name = new.col.names[1],
  data.label = 12L,
  mu.rate.name = new.col.names[2],
  mu.rate = 8.0E-6
)
ado.info.12 <- load.data(
  file = "../simulated_data_12/ado_info.tsv",
  data.label.name = new.col.names[1],
  data.label = 12L,
  mu.rate.name = new.col.names[2],
  mu.rate = 8.0E-6
)
allelic.info.12 <- load.rds(
  file = "../simulated_data_12/allelic_info.rds",
  data.label.name = new.col.names[1],
  data.label = 12L,
  mu.rate.name = new.col.names[2],
  mu.rate = 8.0E-6
)
parallel.mut.12 <- load.data(
  file = "../simulated_data_12/summary.tsv",
  sep = " ",
  data.label.name = new.col.names[1],
  data.label = 12L,
  mu.rate.name = new.col.names[2],
  mu.rate = 8.0E-6
)

# simulation 13
tree.info.13 <- load.data(
  file = "../simulated_data_13/trees_info_updated.tsv",
  data.label.name = new.col.names[1],
  data.label = 13L,
  mu.rate.name = new.col.names[2],
  mu.rate = 3.0E-5
)
site.info.13 <- load.data(
  file = "../simulated_data_13/sieve_sites_info.tsv",
  data.label.name = new.col.names[1],
  data.label = 13L,
  mu.rate.name = new.col.names[2],
  mu.rate = 3.0E-5
)
param.info.13 <- load.data(
  file = "../simulated_data_13/params_info.tsv",
  data.label.name = new.col.names[1],
  data.label = 13L,
  mu.rate.name = new.col.names[2],
  mu.rate = 3.0E-5
)
var.info.13 <- load.data(
  file = "../simulated_data_13/variants_info.tsv",
  data.label.name = new.col.names[1],
  data.label = 13L,
  mu.rate.name = new.col.names[2],
  mu.rate = 3.0E-5
)
ado.info.13 <- load.data(
  file = "../simulated_data_13/ado_info.tsv",
  data.label.name = new.col.names[1],
  data.label = 13L,
  mu.rate.name = new.col.names[2],
  mu.rate = 3.0E-5
)
allelic.info.13 <- load.rds(
  file = "../simulated_data_13/allelic_info.rds",
  data.label.name = new.col.names[1],
  data.label = 13L,
  mu.rate.name = new.col.names[2],
  mu.rate = 3.0E-5
)
parallel.mut.13 <- load.data(
  file = "../simulated_data_13/summary.tsv",
  sep = " ",
  data.label.name = new.col.names[1],
  data.label = 13L,
  mu.rate.name = new.col.names[2],
  mu.rate = 3.0E-5
)

# simulation 21
tree.info.21 <- load.data(
  file = "../simulated_data_21/trees_info_updated.tsv",
  data.label.name = new.col.names[1],
  data.label = 21L,
  mu.rate.name = new.col.names[2],
  mu.rate = 1.0E-6
)
site.info.21 <- load.data(
  file = "../simulated_data_21/sieve_sites_info.tsv",
  data.label.name = new.col.names[1],
  data.label = 21L,
  mu.rate.name = new.col.names[2],
  mu.rate = 1.0E-6
)
param.info.21 <- load.data(
  file = "../simulated_data_21/params_info.tsv",
  data.label.name = new.col.names[1],
  data.label = 21L,
  mu.rate.name = new.col.names[2],
  mu.rate = 1.0E-6
)
var.info.21 <- load.data(
  file = "../simulated_data_21/variants_info.tsv",
  data.label.name = new.col.names[1],
  data.label = 21L,
  mu.rate.name = new.col.names[2],
  mu.rate = 1.0E-6
)
ado.info.21 <- load.data(
  file = "../simulated_data_21/ado_info.tsv",
  data.label.name = new.col.names[1],
  data.label = 21L,
  mu.rate.name = new.col.names[2],
  mu.rate = 1.0E-6
)
allelic.info.21 <- load.rds(
  file = "../simulated_data_21/allelic_info.rds",
  data.label.name = new.col.names[1],
  data.label = 21L,
  mu.rate.name = new.col.names[2],
  mu.rate = 1.0E-6
)
parallel.mut.21 <- load.data(
  file = "../simulated_data_21/summary.tsv",
  sep = " ",
  data.label.name = new.col.names[1],
  data.label = 21L,
  mu.rate.name = new.col.names[2],
  mu.rate = 1.0E-6
)

# simulation 22
tree.info.22 <- load.data(
  file = "../simulated_data_22/trees_info_updated.tsv",
  data.label.name = new.col.names[1],
  data.label = 22L,
  mu.rate.name = new.col.names[2],
  mu.rate = 8.0E-6
)
site.info.22 <- load.data(
  file = "../simulated_data_22/sieve_sites_info.tsv",
  data.label.name = new.col.names[1],
  data.label = 22L,
  mu.rate.name = new.col.names[2],
  mu.rate = 8.0E-6
)
param.info.22 <- load.data(
  file = "../simulated_data_22/params_info.tsv",
  data.label.name = new.col.names[1],
  data.label = 22L,
  mu.rate.name = new.col.names[2],
  mu.rate = 8.0E-6
)
var.info.22 <- load.data(
  file = "../simulated_data_22/variants_info.tsv",
  data.label.name = new.col.names[1],
  data.label = 22L,
  mu.rate.name = new.col.names[2],
  mu.rate = 8.0E-6
)
ado.info.22 <- load.data(
  file = "../simulated_data_22/ado_info.tsv",
  data.label.name = new.col.names[1],
  data.label = 22L,
  mu.rate.name = new.col.names[2],
  mu.rate = 8.0E-6
)
allelic.info.22 <- load.rds(
  file = "../simulated_data_22/allelic_info.rds",
  data.label.name = new.col.names[1],
  data.label = 22L,
  mu.rate.name = new.col.names[2],
  mu.rate = 8.0E-6
)
parallel.mut.22 <- load.data(
  file = "../simulated_data_22/summary.tsv",
  sep = " ",
  data.label.name = new.col.names[1],
  data.label = 22L,
  mu.rate.name = new.col.names[2],
  mu.rate = 8.0E-6
)

# simulation 23
tree.info.23 <- load.data(
  file = "../simulated_data_23/trees_info_updated.tsv",
  data.label.name = new.col.names[1],
  data.label = 23L,
  mu.rate.name = new.col.names[2],
  mu.rate = 3.0E-5
)
site.info.23 <- load.data(
  file = "../simulated_data_23/sieve_sites_info.tsv",
  data.label.name = new.col.names[1],
  data.label = 23L,
  mu.rate.name = new.col.names[2],
  mu.rate = 3.0E-5
)
param.info.23 <- load.data(
  file = "../simulated_data_23/params_info.tsv",
  data.label.name = new.col.names[1],
  data.label = 23L,
  mu.rate.name = new.col.names[2],
  mu.rate = 3.0E-5
)
var.info.23 <- load.data(
  file = "../simulated_data_23/variants_info.tsv",
  data.label.name = new.col.names[1],
  data.label = 23L,
  mu.rate.name = new.col.names[2],
  mu.rate = 3.0E-5
)
ado.info.23 <- load.data(
  file = "../simulated_data_23/ado_info.tsv",
  data.label.name = new.col.names[1],
  data.label = 23L,
  mu.rate.name = new.col.names[2],
  mu.rate = 3.0E-5
)
allelic.info.23 <- load.rds(
  file = "../simulated_data_23/allelic_info.rds",
  data.label.name = new.col.names[1],
  data.label = 23L,
  mu.rate.name = new.col.names[2],
  mu.rate = 3.0E-5
)
parallel.mut.23 <- load.data(
  file = "../simulated_data_23/summary.tsv",
  sep = " ",
  data.label.name = new.col.names[1],
  data.label = 23L,
  mu.rate.name = new.col.names[2],
  mu.rate = 3.0E-5
)

# simulation 31
tree.info.31 <- load.data(
  file = "../simulated_data_31/trees_info_updated.tsv",
  data.label.name = new.col.names[1],
  data.label = 31L,
  mu.rate.name = new.col.names[2],
  mu.rate = 1.0E-6
)
site.info.31 <- load.data(
  file = "../simulated_data_31/sieve_sites_info.tsv",
  data.label.name = new.col.names[1],
  data.label = 31L,
  mu.rate.name = new.col.names[2],
  mu.rate = 1.0E-6
)
param.info.31 <- load.data(
  file = "../simulated_data_31/params_info.tsv",
  data.label.name = new.col.names[1],
  data.label = 31L,
  mu.rate.name = new.col.names[2],
  mu.rate = 1.0E-6
)
var.info.31 <- load.data(
  file = "../simulated_data_31/variants_info.tsv",
  data.label.name = new.col.names[1],
  data.label = 31L,
  mu.rate.name = new.col.names[2],
  mu.rate = 1.0E-6
)
ado.info.31 <- load.data(
  file = "../simulated_data_31/ado_info.tsv",
  data.label.name = new.col.names[1],
  data.label = 31L,
  mu.rate.name = new.col.names[2],
  mu.rate = 1.0E-6
)
allelic.info.31 <- load.rds(
  file = "../simulated_data_31/allelic_info.rds",
  data.label.name = new.col.names[1],
  data.label = 31L,
  mu.rate.name = new.col.names[2],
  mu.rate = 1.0E-6
)
parallel.mut.31 <- load.data(
  file = "../simulated_data_31/summary.tsv",
  sep = " ",
  data.label.name = new.col.names[1],
  data.label = 31L,
  mu.rate.name = new.col.names[2],
  mu.rate = 1.0E-6
)

# simulation 32
tree.info.32 <- load.data(
  file = "../simulated_data_32/trees_info_updated.tsv",
  data.label.name = new.col.names[1],
  data.label = 32L,
  mu.rate.name = new.col.names[2],
  mu.rate = 8.0E-6
)
site.info.32 <- load.data(
  file = "../simulated_data_32/sieve_sites_info.tsv",
  data.label.name = new.col.names[1],
  data.label = 32L,
  mu.rate.name = new.col.names[2],
  mu.rate = 8.0E-6
)
param.info.32 <- load.data(
  file = "../simulated_data_32/params_info.tsv",
  data.label.name = new.col.names[1],
  data.label = 32L,
  mu.rate.name = new.col.names[2],
  mu.rate = 8.0E-6
)
var.info.32 <- load.data(
  file = "../simulated_data_32/variants_info.tsv",
  data.label.name = new.col.names[1],
  data.label = 32L,
  mu.rate.name = new.col.names[2],
  mu.rate = 8.0E-6
)
ado.info.32 <- load.data(
  file = "../simulated_data_32/ado_info.tsv",
  data.label.name = new.col.names[1],
  data.label = 32L,
  mu.rate.name = new.col.names[2],
  mu.rate = 8.0E-6
)
allelic.info.32 <- load.rds(
  file = "../simulated_data_32/allelic_info.rds",
  data.label.name = new.col.names[1],
  data.label = 32L,
  mu.rate.name = new.col.names[2],
  mu.rate = 8.0E-6
)
parallel.mut.32 <- load.data(
  file = "../simulated_data_32/summary.tsv",
  sep = " ",
  data.label.name = new.col.names[1],
  data.label = 32L,
  mu.rate.name = new.col.names[2],
  mu.rate = 8.0E-6
)

# simulation 33
tree.info.33 <- load.data(
  file = "../simulated_data_33/trees_info_updated.tsv",
  data.label.name = new.col.names[1],
  data.label = 33L,
  mu.rate.name = new.col.names[2],
  mu.rate = 3.0E-5
)
site.info.33 <- load.data(
  file = "../simulated_data_33/sieve_sites_info.tsv",
  data.label.name = new.col.names[1],
  data.label = 33L,
  mu.rate.name = new.col.names[2],
  mu.rate = 3.0E-5
)
param.info.33 <- load.data(
  file = "../simulated_data_33/params_info.tsv",
  data.label.name = new.col.names[1],
  data.label = 33L,
  mu.rate.name = new.col.names[2],
  mu.rate = 3.0E-5
)
var.info.33 <- load.data(
  file = "../simulated_data_33/variants_info.tsv",
  data.label.name = new.col.names[1],
  data.label = 33L,
  mu.rate.name = new.col.names[2],
  mu.rate = 3.0E-5
)
ado.info.33 <- load.data(
  file = "../simulated_data_33/ado_info.tsv",
  data.label.name = new.col.names[1],
  data.label = 33L,
  mu.rate.name = new.col.names[2],
  mu.rate = 3.0E-5
)
allelic.info.33 <- load.rds(
  file = "../simulated_data_33/allelic_info.rds",
  data.label.name = new.col.names[1],
  data.label = 33L,
  mu.rate.name = new.col.names[2],
  mu.rate = 3.0E-5
)
parallel.mut.33 <- load.data(
  file = "../simulated_data_33/summary.tsv",
  sep = " ",
  data.label.name = new.col.names[1],
  data.label = 33L,
  mu.rate.name = new.col.names[2],
  mu.rate = 3.0E-5
)

# simulation 41
tree.info.41 <- load.data(
  file = "../simulated_data_41/trees_info_updated.tsv",
  data.label.name = new.col.names[1],
  data.label = 41L,
  mu.rate.name = new.col.names[2],
  mu.rate = 1.0E-6
)
site.info.41 <- load.data(
  file = "../simulated_data_41/sieve_sites_info.tsv",
  data.label.name = new.col.names[1],
  data.label = 41L,
  mu.rate.name = new.col.names[2],
  mu.rate = 1.0E-6
)
param.info.41 <- load.data(
  file = "../simulated_data_41/params_info.tsv",
  data.label.name = new.col.names[1],
  data.label = 41L,
  mu.rate.name = new.col.names[2],
  mu.rate = 1.0E-6
)
var.info.41 <- load.data(
  file = "../simulated_data_41/variants_info.tsv",
  data.label.name = new.col.names[1],
  data.label = 41L,
  mu.rate.name = new.col.names[2],
  mu.rate = 1.0E-6
)
ado.info.41 <- load.data(
  file = "../simulated_data_41/ado_info.tsv",
  data.label.name = new.col.names[1],
  data.label = 41L,
  mu.rate.name = new.col.names[2],
  mu.rate = 1.0E-6
)
allelic.info.41 <- load.rds(
  file = "../simulated_data_41/allelic_info.rds",
  data.label.name = new.col.names[1],
  data.label = 41L,
  mu.rate.name = new.col.names[2],
  mu.rate = 1.0E-6
)
parallel.mut.41 <- load.data(
  file = "../simulated_data_41/summary.tsv",
  sep = " ",
  data.label.name = new.col.names[1],
  data.label = 41L,
  mu.rate.name = new.col.names[2],
  mu.rate = 1.0E-6
)

# simulation 42
tree.info.42 <- load.data(
  file = "../simulated_data_42/trees_info_updated.tsv",
  data.label.name = new.col.names[1],
  data.label = 42L,
  mu.rate.name = new.col.names[2],
  mu.rate = 8.0E-6
)
site.info.42 <- load.data(
  file = "../simulated_data_42/sieve_sites_info.tsv",
  data.label.name = new.col.names[1],
  data.label = 42L,
  mu.rate.name = new.col.names[2],
  mu.rate = 8.0E-6
)
param.info.42 <- load.data(
  file = "../simulated_data_42/params_info.tsv",
  data.label.name = new.col.names[1],
  data.label = 42L,
  mu.rate.name = new.col.names[2],
  mu.rate = 8.0E-6
)
var.info.42 <- load.data(
  file = "../simulated_data_42/variants_info.tsv",
  data.label.name = new.col.names[1],
  data.label = 42L,
  mu.rate.name = new.col.names[2],
  mu.rate = 8.0E-6
)
ado.info.42 <- load.data(
  file = "../simulated_data_42/ado_info.tsv",
  data.label.name = new.col.names[1],
  data.label = 42L,
  mu.rate.name = new.col.names[2],
  mu.rate = 8.0E-6
)
allelic.info.42 <- load.rds(
  file = "../simulated_data_42/allelic_info.rds",
  data.label.name = new.col.names[1],
  data.label = 42L,
  mu.rate.name = new.col.names[2],
  mu.rate = 8.0E-6
)
parallel.mut.42 <- load.data(
  file = "../simulated_data_42/summary.tsv",
  sep = " ",
  data.label.name = new.col.names[1],
  data.label = 42L,
  mu.rate.name = new.col.names[2],
  mu.rate = 8.0E-6
)

# simulation 43
tree.info.43 <- load.data(
  file = "../simulated_data_43/trees_info_updated.tsv",
  data.label.name = new.col.names[1],
  data.label = 43L,
  mu.rate.name = new.col.names[2],
  mu.rate = 3.0E-5
)
site.info.43 <- load.data(
  file = "../simulated_data_43/sieve_sites_info.tsv",
  data.label.name = new.col.names[1],
  data.label = 43L,
  mu.rate.name = new.col.names[2],
  mu.rate = 3.0E-5
)
param.info.43 <- load.data(
  file = "../simulated_data_43/params_info.tsv",
  data.label.name = new.col.names[1],
  data.label = 43L,
  mu.rate.name = new.col.names[2],
  mu.rate = 3.0E-5
)
var.info.43 <- load.data(
  file = "../simulated_data_43/variants_info.tsv",
  data.label.name = new.col.names[1],
  data.label = 43L,
  mu.rate.name = new.col.names[2],
  mu.rate = 3.0E-5
)
ado.info.43 <- load.data(
  file = "../simulated_data_43/ado_info.tsv",
  data.label.name = new.col.names[1],
  data.label = 43L,
  mu.rate.name = new.col.names[2],
  mu.rate = 3.0E-5
)
allelic.info.43 <- load.rds(
  file = "../simulated_data_43/allelic_info.rds",
  data.label.name = new.col.names[1],
  data.label = 43L,
  mu.rate.name = new.col.names[2],
  mu.rate = 3.0E-5
)
parallel.mut.43 <- load.data(
  file = "../simulated_data_43/summary.tsv",
  sep = " ",
  data.label.name = new.col.names[1],
  data.label = 43L,
  mu.rate.name = new.col.names[2],
  mu.rate = 3.0E-5
)

# simulation 51
tree.info.51 <- load.data(
  file = "../simulated_data_51/trees_info_updated.tsv",
  data.label.name = new.col.names[1],
  data.label = 51L,
  mu.rate.name = new.col.names[2],
  mu.rate = 1.0E-6
)
site.info.51 <- load.data(
  file = "../simulated_data_51/sieve_sites_info.tsv",
  data.label.name = new.col.names[1],
  data.label = 51L,
  mu.rate.name = new.col.names[2],
  mu.rate = 1.0E-6
)
param.info.51 <- load.data(
  file = "../simulated_data_51/params_info.tsv",
  data.label.name = new.col.names[1],
  data.label = 51L,
  mu.rate.name = new.col.names[2],
  mu.rate = 1.0E-6
)
var.info.51 <- load.data(
  file = "../simulated_data_51/variants_info.tsv",
  data.label.name = new.col.names[1],
  data.label = 51L,
  mu.rate.name = new.col.names[2],
  mu.rate = 1.0E-6
)
ado.info.51 <- load.data(
  file = "../simulated_data_51/ado_info.tsv",
  data.label.name = new.col.names[1],
  data.label = 51L,
  mu.rate.name = new.col.names[2],
  mu.rate = 1.0E-6
)
allelic.info.51 <- load.rds(
  file = "../simulated_data_51/allelic_info.rds",
  data.label.name = new.col.names[1],
  data.label = 51L,
  mu.rate.name = new.col.names[2],
  mu.rate = 1.0E-6
)
parallel.mut.51 <- load.data(
  file = "../simulated_data_51/summary.tsv",
  sep = " ",
  data.label.name = new.col.names[1],
  data.label = 51L,
  mu.rate.name = new.col.names[2],
  mu.rate = 1.0E-6
)

# simulation 52
tree.info.52 <- load.data(
  file = "../simulated_data_52/trees_info_updated.tsv",
  data.label.name = new.col.names[1],
  data.label = 52L,
  mu.rate.name = new.col.names[2],
  mu.rate = 8.0E-6
)
site.info.52 <- load.data(
  file = "../simulated_data_52/sieve_sites_info.tsv",
  data.label.name = new.col.names[1],
  data.label = 52L,
  mu.rate.name = new.col.names[2],
  mu.rate = 8.0E-6
)
param.info.52 <- load.data(
  file = "../simulated_data_52/params_info.tsv",
  data.label.name = new.col.names[1],
  data.label = 52L,
  mu.rate.name = new.col.names[2],
  mu.rate = 8.0E-6
)
var.info.52 <- load.data(
  file = "../simulated_data_52/variants_info.tsv",
  data.label.name = new.col.names[1],
  data.label = 52L,
  mu.rate.name = new.col.names[2],
  mu.rate = 8.0E-6
)
ado.info.52 <- load.data(
  file = "../simulated_data_52/ado_info.tsv",
  data.label.name = new.col.names[1],
  data.label = 52L,
  mu.rate.name = new.col.names[2],
  mu.rate = 8.0E-6
)
allelic.info.52 <- load.rds(
  file = "../simulated_data_52/allelic_info.rds",
  data.label.name = new.col.names[1],
  data.label = 52L,
  mu.rate.name = new.col.names[2],
  mu.rate = 8.0E-6
)
parallel.mut.52 <- load.data(
  file = "../simulated_data_52/summary.tsv",
  sep = " ",
  data.label.name = new.col.names[1],
  data.label = 52L,
  mu.rate.name = new.col.names[2],
  mu.rate = 8.0E-6
)

# simulation 53
tree.info.53 <- load.data(
  file = "../simulated_data_53/trees_info_updated.tsv",
  data.label.name = new.col.names[1],
  data.label = 53L,
  mu.rate.name = new.col.names[2],
  mu.rate = 3.0E-5
)
site.info.53 <- load.data(
  file = "../simulated_data_53/sieve_sites_info.tsv",
  data.label.name = new.col.names[1],
  data.label = 53L,
  mu.rate.name = new.col.names[2],
  mu.rate = 3.0E-5
)
param.info.53 <- load.data(
  file = "../simulated_data_53/params_info.tsv",
  data.label.name = new.col.names[1],
  data.label = 53L,
  mu.rate.name = new.col.names[2],
  mu.rate = 3.0E-5
)
var.info.53 <- load.data(
  file = "../simulated_data_53/variants_info.tsv",
  data.label.name = new.col.names[1],
  data.label = 53L,
  mu.rate.name = new.col.names[2],
  mu.rate = 3.0E-5
)
ado.info.53 <- load.data(
  file = "../simulated_data_53/ado_info.tsv",
  data.label.name = new.col.names[1],
  data.label = 53L,
  mu.rate.name = new.col.names[2],
  mu.rate = 3.0E-5
)
allelic.info.53 <- load.rds(
  file = "../simulated_data_53/allelic_info.rds",
  data.label.name = new.col.names[1],
  data.label = 53L,
  mu.rate.name = new.col.names[2],
  mu.rate = 3.0E-5
)
parallel.mut.53 <- load.data(
  file = "../simulated_data_53/summary.tsv",
  sep = " ",
  data.label.name = new.col.names[1],
  data.label = 53L,
  mu.rate.name = new.col.names[2],
  mu.rate = 3.0E-5
)

# combine data
tree.info <- rbind(
  tree.info.2, tree.info.3, tree.info.4,
  tree.info.11, tree.info.12, tree.info.13,
  tree.info.21, tree.info.22, tree.info.23,
  tree.info.31, tree.info.32, tree.info.33,
  tree.info.41, tree.info.42, tree.info.43,
  tree.info.51, tree.info.52, tree.info.53
)
rm(
  tree.info.2, tree.info.3, tree.info.4,
  tree.info.11, tree.info.12, tree.info.13,
  tree.info.21, tree.info.22, tree.info.23,
  tree.info.31, tree.info.32, tree.info.33,
  tree.info.41, tree.info.42, tree.info.43,
  tree.info.51, tree.info.52, tree.info.53
)
site.info <- rbind(
  site.info.2, site.info.3, site.info.4,
  site.info.11, site.info.12, site.info.13,
  site.info.21, site.info.22, site.info.23,
  site.info.31, site.info.32, site.info.33,
  site.info.41, site.info.42, site.info.43,
  site.info.51, site.info.52, site.info.53
)
rm(
  site.info.2, site.info.3, site.info.4,
  site.info.11, site.info.12, site.info.13,
  site.info.21, site.info.22, site.info.23,
  site.info.31, site.info.32, site.info.33,
  site.info.41, site.info.42, site.info.43,
  site.info.51, site.info.52, site.info.53
)
param.info <- rbind(
  param.info.2, param.info.3, param.info.4,
  param.info.11, param.info.12, param.info.13,
  param.info.21, param.info.22, param.info.23,
  param.info.31, param.info.32, param.info.33,
  param.info.41, param.info.42, param.info.43,
  param.info.51, param.info.52, param.info.53
)
rm(
  param.info.2, param.info.3, param.info.4,
  param.info.11, param.info.12, param.info.13,
  param.info.21, param.info.22, param.info.23,
  param.info.31, param.info.32, param.info.33,
  param.info.41, param.info.42, param.info.43,
  param.info.51, param.info.52, param.info.53
)
var.info <- rbind(
  var.info.2, var.info.3, var.info.4,
  var.info.11, var.info.12, var.info.13,
  var.info.21, var.info.22, var.info.23,
  var.info.31, var.info.32, var.info.33,
  var.info.41, var.info.42, var.info.43,
  var.info.51, var.info.52, var.info.53
)
rm(
  var.info.2, var.info.3, var.info.4,
  var.info.11, var.info.12, var.info.13,
  var.info.21, var.info.22, var.info.23,
  var.info.31, var.info.32, var.info.33,
  var.info.41, var.info.42, var.info.43,
  var.info.51, var.info.52, var.info.53
)
ado.info <- rbind(
  ado.info.2, ado.info.3, ado.info.4,
  ado.info.11, ado.info.12, ado.info.13,
  ado.info.21, ado.info.22, ado.info.23,
  ado.info.31, ado.info.32, ado.info.33,
  ado.info.41, ado.info.42, ado.info.43,
  ado.info.51, ado.info.52, ado.info.53
)
rm(
  ado.info.2, ado.info.3, ado.info.4,
  ado.info.11, ado.info.12, ado.info.13,
  ado.info.21, ado.info.22, ado.info.23,
  ado.info.31, ado.info.32, ado.info.33,
  ado.info.41, ado.info.42, ado.info.43,
  ado.info.51, ado.info.52, ado.info.53
)
allelic.info <- rbind(
  allelic.info.2, allelic.info.3, allelic.info.4,
  allelic.info.11, allelic.info.12, allelic.info.13,
  allelic.info.21, allelic.info.22, allelic.info.23,
  allelic.info.31, allelic.info.32, allelic.info.33,
  allelic.info.41, allelic.info.42, allelic.info.43,
  allelic.info.51, allelic.info.52, allelic.info.53
)
rm(
  allelic.info.2, allelic.info.3, allelic.info.4,
  allelic.info.11, allelic.info.12, allelic.info.13,
  allelic.info.21, allelic.info.22, allelic.info.23,
  allelic.info.31, allelic.info.32, allelic.info.33,
  allelic.info.41, allelic.info.42, allelic.info.43,
  allelic.info.51, allelic.info.52, allelic.info.53
)
parallel.mut <- rbind(
  parallel.mut.2, parallel.mut.3, parallel.mut.4,
  parallel.mut.11, parallel.mut.12, parallel.mut.13,
  parallel.mut.21, parallel.mut.22, parallel.mut.23,
  parallel.mut.31, parallel.mut.32, parallel.mut.33,
  parallel.mut.41, parallel.mut.42, parallel.mut.43,
  parallel.mut.51, parallel.mut.52, parallel.mut.53
)
rm(
  parallel.mut.2, parallel.mut.3, parallel.mut.4,
  parallel.mut.11, parallel.mut.12, parallel.mut.13,
  parallel.mut.21, parallel.mut.22, parallel.mut.23,
  parallel.mut.31, parallel.mut.32, parallel.mut.33,
  parallel.mut.41, parallel.mut.42, parallel.mut.43,
  parallel.mut.51, parallel.mut.52, parallel.mut.53
)

# For facets
cell.num.labels <- c("40 cells", "100 cells")
names(cell.num.labels) <- c("40", "100")
allelic.raw.var <- c(
  "High mean\nLow variance",
  "High mean\nMedium variance",
  "Low mean\nHigh variance"
)
names(allelic.raw.var) <- c("2", "10", "20")

# process columns
tree.info$cell_num <- as.factor(tree.info$cell_num)
tree.info$coverage_mean <- as.factor(tree.info$coverage_mean)
tree.info$coverage_variance <- as.factor(tree.info$coverage_variance)
tree.info$dataset <- as.factor(tree.info$dataset)
tree.info$tool <- as.factor(tree.info$tool)
tree.info$snv_type <- as.factor(tree.info$snv_type)
tree.info$tool_setup <- as.factor(tree.info$tool_setup)
tree.info$max_clades <- as.numeric(tree.info$max_clades)
tree.info$RF_distance <- as.numeric(tree.info$RF_distance)
tree.info$normalized_RF_distance <- as.numeric(tree.info$normalized_RF_distance)
tree.info$weighted_RF_distance <- as.numeric(tree.info$weighted_RF_distance)
tree.info$rooted_branch_score_difference <- as.numeric(tree.info$rooted_branch_score_difference)
tree.info[[new.col.names[1]]] <- as.factor(tree.info[[new.col.names[1]]])
tree.info[[new.col.names[2]]] <- as.factor(tree.info[[new.col.names[2]]])

site.info$cell_num <- as.factor(site.info$cell_num)
site.info$coverage_mean <- as.factor(site.info$coverage_mean)
site.info$coverage_variance <- as.factor(site.info$coverage_variance)
site.info$dataset <- as.factor(site.info$dataset)
site.info$tool <- as.factor(site.info$tool)
site.info$snv_type <- as.factor(site.info$snv_type)
site.info$tool_setup <- as.factor(site.info$tool_setup)
site.info$max_clades <- as.numeric(site.info$max_clades)
site.info$RF_distance <- as.numeric(site.info$RF_distance)
site.info$normalized_RF_distance <- as.numeric(site.info$normalized_RF_distance)
site.info$weighted_RF_distance <- as.numeric(site.info$weighted_RF_distance)
site.info$rooted_branch_score_difference <- as.numeric(site.info$rooted_branch_score_difference)
site.info$num_candidate_mutated_sites <- as.integer(site.info$num_candidate_mutated_sites)
site.info$num_background_sites <- as.integer(site.info$num_background_sites)
site.info$log10_num_background_sites <- log10(site.info$num_background_sites)
site.info[[new.col.names[1]]] <- as.factor(site.info[[new.col.names[1]]])
site.info[[new.col.names[2]]] <- as.factor(site.info[[new.col.names[2]]])

param.info$cell_num <- as.factor(param.info$cell_num)
param.info$coverage_mean <- as.factor(param.info$coverage_mean)
param.info$coverage_variance <- as.factor(param.info$coverage_variance)
param.info$dataset <- as.factor(param.info$dataset)
param.info$tool <- as.factor(param.info$tool)
param.info$snv_type <- as.factor(param.info$snv_type)
param.info$tool_setup <- as.factor(param.info$tool_setup)
param.info$eff_seq_err_rate <- as.numeric(param.info$eff_seq_err_rate)
param.info$allelic_seq_cov <- as.numeric(param.info$allelic_seq_cov)
param.info$allelic_seq_cov_raw_var <- as.numeric(param.info$allelic_seq_cov_raw_var)
param.info$ado_rate <- as.numeric(param.info$ado_rate)
param.info$gamma_shape <- as.numeric(param.info$gamma_shape)
param.info$wild_overdispersion <- as.numeric(param.info$wild_overdispersion)
param.info$alternative_overdispersion <- as.numeric(param.info$alternative_overdispersion)
param.info$zygosity_rate <- as.numeric(param.info$zygosity_rate)
param.info$estimates_type <- as.factor(param.info$estimates_type)
param.info$deletion_rate <- as.numeric(param.info$deletion_rate)
param.info$insertion_rate <- as.numeric(param.info$insertion_rate)
param.info$population_size <- as.numeric(param.info$population_size)
param.info[[new.col.names[1]]] <- as.factor(param.info[[new.col.names[1]]])
param.info[[new.col.names[2]]] <- as.factor(param.info[[new.col.names[2]]])

var.info$cell_num <- as.factor(var.info$cell_num)
var.info$coverage_mean <- as.factor(var.info$coverage_mean)
var.info$coverage_variance <- as.factor(var.info$coverage_variance)
var.info$dataset <- as.factor(var.info$dataset)
var.info$tool <- as.factor(var.info$tool)
var.info$snv_type <- as.factor(var.info$snv_type)
var.info$tool_setup <- as.factor(var.info$tool_setup)
var.info$true_positive <- as.numeric(var.info$true_positive)
var.info$false_positive <- as.numeric(var.info$false_positive)
var.info$true_negative <- as.numeric(var.info$true_negative)
var.info$false_negative <- as.numeric(var.info$false_negative)
var.info$recall <- as.numeric(var.info$recall)
var.info$precision <- as.numeric(var.info$precision)
var.info$fall_out <- as.numeric(var.info$fall_out)
var.info$f1_score <- as.numeric(var.info$f1_score)
var.info$true_positive_hetero_mu <- as.numeric(var.info$true_positive_hetero_mu)
var.info$false_positive_hetero_mu <- as.numeric(var.info$false_positive_hetero_mu)
var.info$true_negative_hetero_mu <- as.numeric(var.info$true_negative_hetero_mu)
var.info$false_negative_hetero_mu <- as.numeric(var.info$false_negative_hetero_mu)
var.info$recall_hetero_mu <- as.numeric(var.info$recall_hetero_mu)
var.info$precision_hetero_mu <- as.numeric(var.info$precision_hetero_mu)
var.info$fall_out_hetero_mu <- as.numeric(var.info$fall_out_hetero_mu)
var.info$f1_score_hetero_mu <- as.numeric(var.info$f1_score_hetero_mu)
var.info$true_positive_homo_mu <- as.numeric(var.info$true_positive_homo_mu)
var.info$false_positive_homo_mu <- as.numeric(var.info$false_positive_homo_mu)
var.info$true_negative_homo_mu <- as.numeric(var.info$true_negative_homo_mu)
var.info$false_negative_homo_mu <- as.numeric(var.info$false_negative_homo_mu)
var.info$recall_homo_mu <- as.numeric(var.info$recall_homo_mu)
var.info$precision_homo_mu <- as.numeric(var.info$precision_homo_mu)
var.info$fall_out_homo_mu <- as.numeric(var.info$fall_out_homo_mu)
var.info$f1_score_homo_mu <- as.numeric(var.info$f1_score_homo_mu)
var.info[[new.col.names[1]]] <- as.factor(var.info[[new.col.names[1]]])
var.info <- var.info %>%
  mutate(prop_true_homo_ref_as_hetero_mu_in_called_pos = (false_positive_hetero_mu - true_homo_mu_as_hetero_mu) / (false_positive_hetero_mu + true_positive_hetero_mu)) %>%
  mutate(prop_true_homo_mu_as_hetero_mu_in_called_pos = true_homo_mu_as_hetero_mu / (false_positive_hetero_mu + true_positive_hetero_mu))

fsa.var.info <- var.info[var.info[[new.col.names[2]]] > isa.mu.rate.threshold,]
var.info[[new.col.names[2]]] <- as.factor(var.info[[new.col.names[2]]])
fsa.var.info[[new.col.names[2]]] <- as.factor(fsa.var.info[[new.col.names[2]]])

ado.info$cell_num <- as.factor(ado.info$cell_num)
ado.info$coverage_mean <- as.factor(ado.info$coverage_mean)
ado.info$coverage_variance <- as.factor(ado.info$coverage_variance)
ado.info$dataset <- as.factor(ado.info$dataset)
ado.info$tool <- as.factor(ado.info$tool)
ado.info$snv_type <- as.factor(ado.info$snv_type)
ado.info$tool_setup <- as.factor(ado.info$tool_setup)
ado.info$true_positive <- as.numeric(ado.info$true_positive)
ado.info$false_positive <- as.numeric(ado.info$false_positive)
ado.info$true_negative <- as.numeric(ado.info$true_negative)
ado.info$false_negative <- as.numeric(ado.info$false_negative)
ado.info$recall <- as.numeric(ado.info$recall)
ado.info$precision <- as.numeric(ado.info$precision)
ado.info$fall_out <- as.numeric(ado.info$fall_out)
ado.info$f1_score <- as.numeric(ado.info$f1_score)
ado.info$true_positive_one_ado <- as.numeric(ado.info$true_positive_one_ado)
ado.info$false_positive_one_ado <- as.numeric(ado.info$false_positive_one_ado)
ado.info$true_negative_one_ado <- as.numeric(ado.info$true_negative_one_ado)
ado.info$false_negative_one_ado <- as.numeric(ado.info$false_negative_one_ado)
ado.info$recall_one_ado <- as.numeric(ado.info$recall_one_ado)
ado.info$precision_one_ado <- as.numeric(ado.info$precision_one_ado)
ado.info$fall_out_one_ado <- as.numeric(ado.info$fall_out_one_ado)
ado.info$f1_score_one_ado <- as.numeric(ado.info$f1_score_one_ado)
ado.info[[new.col.names[1]]] <- as.factor(ado.info[[new.col.names[1]]])
ado.info[[new.col.names[2]]] <- as.factor(ado.info[[new.col.names[2]]])

allelic.info$cell_num <- as.factor(allelic.info$cell_num)
allelic.info$coverage_mean <- as.factor(allelic.info$coverage_mean)
allelic.info$coverage_variance <- as.factor(allelic.info$coverage_variance)
allelic.info$dataset <- as.factor(allelic.info$dataset)
allelic.info$tool <- as.factor(allelic.info$tool)
allelic.info$tool_setup <- as.factor(allelic.info$tool_setup)
allelic.info$cell_names <- as.factor(allelic.info$cell_names)
allelic.info[[new.col.names[1]]] <- as.factor(allelic.info[[new.col.names[1]]])
allelic.info[[new.col.names[2]]] <- as.factor(allelic.info[[new.col.names[2]]])

# Rename tools.
tree.info <- tree.info %>%
  mutate(tool = case_when(
    grepl("^cellphy", tool, ignore.case = TRUE, fixed = FALSE) ~ sub("^cellphy", "CellPhy", x = tool, fixed = FALSE),
    grepl("^sciphi", tool, ignore.case = TRUE, fixed = FALSE) ~ sub("^sciphi", "SCIPhI", tool, fixed = FALSE),
    grepl("^sieve", tool, ignore.case = TRUE, fixed = FALSE) ~ sub("^sieve", "SIEVE", tool, fixed = FALSE),
    grepl("^sifit", tool, ignore.case = TRUE, fixed = FALSE) ~ sub("^sifit", "SiFit", tool, fixed = FALSE)
  ))

site.info <- site.info %>%
  mutate(tool = case_when(
    grepl("^cellphy", tool, ignore.case = TRUE, fixed = FALSE) ~ sub("^cellphy", "CellPhy", x = tool, fixed = FALSE),
    grepl("^sciphi", tool, ignore.case = TRUE, fixed = FALSE) ~ sub("^sciphi", "SCIPhI", tool, fixed = FALSE),
    grepl("^sieve", tool, ignore.case = TRUE, fixed = FALSE) ~ sub("^sieve", "SIEVE", tool, fixed = FALSE),
    grepl("^sifit", tool, ignore.case = TRUE, fixed = FALSE) ~ sub("^sifit", "SiFit", tool, fixed = FALSE)
  ))

param.info <- param.info %>%
  mutate(tool = case_when(
    grepl("^cellphy", tool, ignore.case = TRUE, fixed = FALSE) ~ sub("^cellphy", "CellPhy", x = tool, fixed = FALSE),
    grepl("^sciphi", tool, ignore.case = TRUE, fixed = FALSE) ~ sub("^sciphi", "SCIPhI", tool, fixed = FALSE),
    grepl("^sieve", tool, ignore.case = TRUE, fixed = FALSE) ~ sub("^sieve", "SIEVE", tool, fixed = FALSE),
    grepl("^sifit", tool, ignore.case = TRUE, fixed = FALSE) ~ sub("^sifit", "SiFit", tool, fixed = FALSE)
  ))

var.info <- var.info %>%
  mutate(tool = case_when(
    grepl("^sciphi", tool, ignore.case = TRUE, fixed = FALSE) ~ sub("^sciphi", "SCIPhI", tool, fixed = FALSE),
    grepl("^sieve", tool, ignore.case = TRUE, fixed = FALSE) ~ sub("^sieve", "SIEVE", tool, fixed = FALSE),
    grepl("^monovar", tool, ignore.case = TRUE, fixed = FALSE) ~ sub("^monovar", "Monovar", tool, fixed = FALSE)
  ))

fsa.var.info <- fsa.var.info %>%
  mutate(tool = case_when(
    grepl("^sciphi", tool, ignore.case = TRUE, fixed = FALSE) ~ sub("^sciphi", "SCIPhI", tool, fixed = FALSE),
    grepl("^sieve", tool, ignore.case = TRUE, fixed = FALSE) ~ sub("^sieve", "SIEVE", tool, fixed = FALSE),
    grepl("^monovar", tool, ignore.case = TRUE, fixed = FALSE) ~ sub("^monovar", "Monovar", tool, fixed = FALSE)
  ))

ado.info <- ado.info %>%
  mutate(tool = case_when(
    grepl("^sieve", tool, ignore.case = TRUE, fixed = FALSE) ~ sub("^sieve", "SIEVE", tool, fixed = FALSE),
  ))

allelic.info <- allelic.info %>%
  mutate(tool = case_when(
    grepl("^sieve", tool, ignore.case = TRUE, fixed = FALSE) ~ sub("^sieve", "SIEVE", tool, fixed = FALSE),
  ))


# Filter out some tools.
tree.info <- tree.info %>%
  filter(tool == "CellPhy" | tool == "SCIPhI" | tool == "SiFit" | tool_setup != "none_univ")

param.info <- param.info %>%
  filter(tool == "CellPhy" | tool == "SCIPhI" | tool == "SiFit" | tool_setup == "none_univ")

var.info <- var.info %>%
  filter(tool == "SCIPhI" | tool == "MonoVar" | tool_setup != "none_univ")

fsa.var.info <- fsa.var.info %>%
  filter(tool == "SCIPhI" | tool == "MonoVar" | tool_setup != "none_univ")

ado.info <- ado.info %>%
  filter(tool_setup != "none_univ")

# process tool names

# rename legends of tools if necessary and replace 'tool' column
tree.info$tool <- rename.tools(tree.info[c("tool", "snv_type", "tool_setup")])
# set colors
tree.info$tool <- as.factor(tree.info$tool)
tree.info.pretiffied <- prettify.colors(levels(tree.info$tool))
tree.info$tool <- factor(tree.info$tool, levels = tree.info.pretiffied[["tool"]])

# rename legends of tools if necessary and replace 'tool' column
site.info$tool <- rename.tools(site.info[c("tool", "snv_type", "tool_setup")])
# set colors
site.info$tool <- as.factor(site.info$tool)
site.info.pretiffied <- prettify.colors(levels(site.info$tool))
site.info$tool <- factor(site.info$tool, levels = site.info.pretiffied[["tool"]])

# rename legends of tools if necessary and replace 'tool' column
param.info$tool <- rename.tools(param.info[c("tool", "snv_type", "tool_setup")])
# set colors
param.info$tool <- as.factor(param.info$tool)
param.info.pretiffied <- prettify.colors(levels(param.info$tool))
param.info$tool <- factor(param.info$tool, levels = param.info.pretiffied[["tool"]])

# rename legends of tools if necessary and replace 'tool' column
var.info$tool <- rename.tools(var.info[c("tool", "tool_setup")])
fsa.var.info$tool <- rename.tools(fsa.var.info[c("tool", "tool_setup")])
# set colors
var.info$tool <- as.factor(var.info$tool)
var.info.pretiffied <- prettify.colors(levels(var.info$tool))
var.info$tool <- factor(var.info$tool, levels = var.info.pretiffied[["tool"]])
fsa.var.info$tool <- as.factor(fsa.var.info$tool)
fsa.var.info.pretiffied <- prettify.colors(levels(fsa.var.info$tool))
fsa.var.info$tool <- factor(fsa.var.info$tool, levels = fsa.var.info.pretiffied[["tool"]])

homo.var.prop <- var.info %>%
  select(
    cell_num, simulated_data, mutation_rate, coverage_mean, coverage_variance, dataset, tool,
    true_positive_homo_mu, false_positive_homo_mu, true_negative_homo_mu, false_negative_homo_mu
  ) %>%
  mutate(
    true_homo = true_positive_homo_mu + false_negative_homo_mu,
    genotype_sum = true_positive_homo_mu + false_positive_homo_mu + true_negative_homo_mu + false_negative_homo_mu
  ) %>%
  mutate(
    true_homo_prop = (true_positive_homo_mu + false_negative_homo_mu) / genotype_sum
  )

homo.var <- lapply(
  split(homo.var.prop, homo.var.prop$simulated_data),
  function(x) {
    .x <- x %>% filter(tool != "Monovar")
    ret <- list()
    ret$info <- .x %>%
      distinct(cell_num, simulated_data, mutation_rate, coverage_mean, coverage_variance)

    t1 <- .x %>%
      select(dataset, tool, genotype_sum) %>%
      arrange(tool, dataset)
    ret$genotype_sum <- matrix(
      t1$genotype_sum,
      nrow = length(unique(t1$dataset)),
      dimnames = list(
        dataset = distinct(t1, dataset)$dataset,
        tool = distinct(t1, tool)$tool
      )
    )

    t2 <- .x %>%
      select(dataset, tool, true_homo_prop) %>%
      arrange(tool, dataset)
    ret$true_homo_prop <- matrix(
      t2$true_homo_prop,
      nrow = length(unique(t2$dataset)),
      dimnames = list(
        dataset = distinct(t2, dataset)$dataset,
        tool = distinct(t2, tool)$tool
      )
    )

    num_tool <- dim(ret$genotype_sum)[2L]
    if (num_tool > 1L) {
      ret$info$identical_genotype_sum <- all(
        vapply(
          2L:num_tool,
          function(y) identical(ret$genotype_sum[, 1L], ret$genotype_sum[, y]),
          FUN.VALUE = logical(1L),
          USE.NAMES = FALSE
        )
      )

      ret$info$identical_true_homo_prop <- all(
        vapply(
          2L:num_tool,
          function(y) identical(ret$true_homo_prop[, 1L], ret$true_homo_prop[, y]),
          FUN.VALUE = logical(1L),
          USE.NAMES = FALSE
        )
      )
    }

    ret$info$min_genotype_sum <- min(ret$genotype_sum[, 1L])
    ret$info$max_genotype_sum <- max(ret$genotype_sum[, 1L])

    ret$info$min_true_homo_prop <- min(ret$true_homo_prop[, 1L])
    ret$info$max_true_homo_prop <- max(ret$true_homo_prop[, 1L])

    return(ret)
  }
)

homo.var.info <- bind_rows(
  lapply(
    homo.var,
    function(x) x$info
  )
) %>%
  arrange(mutation_rate) %>%
  group_by(mutation_rate) %>%
  mutate(
    global_min_genotype_sum = min(min_genotype_sum),
    global_max_genotype_sum = max(max_genotype_sum),
    global_min_true_homo_prop = min(min_true_homo_prop),
    global_max_true_homo_prop = max(max_true_homo_prop)
  )

parallel.mut.prop <- parallel.mut %>%
  group_by(mutation_rate) %>%
  mutate(
    min_prop_parallel_mu_sites = min(prop_parallel_mu_sites),
    max_prop_parallel_mu_sites = max(prop_parallel_mu_sites),
    min_prop_parallel_mu_entries = min(prop_parallel_mu_entries),
    max_prop_parallel_mu_entries = max(prop_parallel_mu_entries)
  ) %>%
  distinct(
    mutation_rate,
    min_prop_parallel_mu_sites,
    max_prop_parallel_mu_sites,
    min_prop_parallel_mu_entries,
    max_prop_parallel_mu_entries
  ) %>%
  arrange(mutation_rate)

# rename legends of tools if necessary and replace 'tool' column
ado.info$tool <- rename.tools(ado.info[c("tool", "tool_setup")])
# set colors
ado.info$tool <- as.factor(ado.info$tool)
ado.info.pretiffied <- prettify.colors(levels(ado.info$tool))
ado.info$tool <- factor(ado.info$tool, levels = ado.info.pretiffied[["tool"]])

# rename legends of tools if necessary and replace 'tool' column
allelic.info$tool <- rename.tools(allelic.info[c("tool", "tool_setup")])
# set colors
allelic.info$tool <- as.factor(allelic.info$tool)
allelic.info.pretiffied <- prettify.colors(levels(allelic.info$tool))
allelic.info$tool <- factor(allelic.info$tool, levels = allelic.info.pretiffied[["tool"]])

common.legend <- list(
  tool = c(tree.info.pretiffied[[1]][1:2], var.info.pretiffied[[1]]),
  color = c(tree.info.pretiffied[[2]][1:2], var.info.pretiffied[[2]]),
  fill = c(tree.info.pretiffied[[3]][1:2], var.info.pretiffied[[3]])
)

data_quality_legend <- list(
  quality = c(
    "High mean\nLow variance",
    "High mean\nMedium variance",
    "Low mean\nHigh variance"
  ),
  color = common.legend$color[1L:3L],
  fill = common.legend$fill[1L:3L]
)
names(data_quality_legend$color) <- c("2", "10", "20")
names(data_quality_legend$fill) <- c("2", "10", "20")

###################################################
#         Generate tree comparison graphs         #
###################################################

# set no legends for single image
no.legend <- theme(legend.position = 'hidden')

# Tree: normalized RF distance
normalized.rf.dist.plot <-
  ggplot(
    tree.info,
    aes_(
      x = as.name(new.col.names[2]),
      y = ~normalized_RF_distance
    )
  ) +
  geom_boxplot(
    lwd = 1 / .pt,
    fatten = 1.2 / .pt,
    alpha = 1.0,
    outlier.alpha = 1.0,
    outlier.size = dot.size,
    outlier.shape = NA,
    aes_(
      fill = ~tool,
      color = ~tool
    )
  ) +
  scale_fill_manual(
    name = "Method",
    breaks = common.legend[[1]],
    values = common.legend[[3]]
  ) +
  scale_color_manual(
    name = "Method",
    breaks = common.legend[[1]],
    values = common.legend[[2]]
  ) +
  geom_point(
    position = position_jitterdodge(),
    aes_(color = ~tool),
    size = dot.size,
    alpha = 0.5
  ) +
  facet_grid(
    ~coverage_variance ~ ~cell_num,
    labeller = labeller(
      coverage_variance = allelic.raw.var,
      cell_num = cell.num.labels
    )
  ) +
  scale_x_discrete(
    breaks = as.numeric(levels(tree.info[[new.col.names[2]]])),
    label = scientific
  ) +
  ylim(0, NA) +
  geom_hline(
    yintercept = 0,
    linetype="dashed",
    color = "gray43",
    size = 1 / .pt
  ) +
  labs(
    y = "Normalised Robinson-Foulds distance",
    x = "Mutation rate"
  ) +
  theme_bw() +
  theme(
    text = element_text(size = 7),
    axis.text = element_text(size = 7),
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 6),
    legend.title.align = 0.5,
    legend.box.spacing = unit(0.3, "mm"),
    legend.position = "bottom",
    legend.key.size = unit(4, 'mm'),
    strip.text = strip.text,
    strip.background = element_blank(),
    panel.grid = element_blank()
  )

# Tree: rooted branch score difference
rooted.branch.score.diff.plot <-
  ggplot(
    tree.info[which(!is.na(tree.info$rooted_branch_score_difference)),],
    aes_(
      x = as.name(new.col.names[2]),
      y = ~rooted_branch_score_difference
    )
  ) +
  geom_boxplot(
    lwd = 1 / .pt,
    fatten = 1.2 / .pt,
    alpha = 1.0,
    outlier.alpha = 1.0,
    outlier.size = dot.size,
    outlier.shape = NA,
    aes_(
      fill = ~tool,
      color = ~tool
    )
  ) +
  scale_fill_manual(
    name = "Method",
    breaks = common.legend[[1]],
    values = common.legend[[3]]
  ) +
  scale_color_manual(
    name = "Method",
    breaks = common.legend[[1]],
    values = common.legend[[2]]
  ) +
  geom_point(
    position = position_jitterdodge(),
    aes_(color = ~tool),
    size = dot.size,
    alpha = 0.5
  ) +
  facet_grid(
    ~coverage_variance ~ ~cell_num,
    labeller = labeller(
      coverage_variance = allelic.raw.var,
      cell_num = cell.num.labels
    )
  ) +
  scale_x_discrete(
    breaks = as.numeric(levels(tree.info[[new.col.names[2]]])),
    label = scientific
  ) +
  ylim(0, NA) +
  geom_hline(
    yintercept = 0,
    linetype="dashed",
    color = "gray43",
    size = 1 / .pt
  ) +
  labs(
    y = "Rooted branch score distance",
    x = "Mutation rate"
  ) +
  theme_bw() +
  theme(
    text = element_text(size = 7),
    axis.text = element_text(size = 7),
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 6),
    legend.title.align = 0.5,
    legend.box.spacing = unit(0.3, "mm"),
    legend.position = "bottom",
    legend.key.size = unit(4, 'mm'),
    strip.text = strip.text,
    strip.text.y = element_blank(),
    strip.background = element_blank(),
    panel.grid = element_blank()
  )

# Correlation of the number of background sites and the rooted branch score difference.
cor.plot <- ggscatter(
    site.info,
    x = "log10_num_background_sites",
    y = "rooted_branch_score_difference",
    color = new.col.names[2],
    size = 2 / .pt,
    xlab = "Log10 of number of background sites",
    ylab = "Rooted branch score distance",
    facet.by = c(
      "coverage_variance",
      "cell_num"
      ),
    panel.labs = list(
      coverage_variance = allelic.raw.var,
      cell_num = cell.num.labels
    ),
    add = "reg.line",
    add.params = list(
      size = 1 / .pt,
      color = "red",
      fill = "gray"
      ),
    conf.int = TRUE
  ) +
  scale_color_discrete(
    name = "Mutation rate",
    breaks = as.numeric(levels(site.info[[new.col.names[2]]])),
    labels = scientific
  ) +
  scale_x_continuous(
    breaks = seq(3.5, 5.5, 1),
    limits = c(3.5, NA),
    expand = expansion(mult = c(0.1, 0.3))
  ) +
  stat_cor(
    method = "kendall",
    cor.coef.name = "tau",
    label.x.npc = 0.35,
    size = 6 / .pt
  ) +
  theme_bw() +
  theme(
    text = element_text(size = 7),
    axis.text = element_text(size = 7),
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 6),
    legend.title.align = 0.5,
    legend.box.spacing = unit(0.3, "mm"),
    legend.position = "bottom",
    legend.key.size = unit(4, 'mm'),
    strip.text = strip.text,
    strip.background = element_blank(),
    panel.grid = element_blank()
  )

# Proportions of homozygous mutations.
# prop.homo.mu.plot <-
#   ggplot(
#     homo.var.prop %>%
#     distinct(cell_num, simulated_data, mutation_rate, coverage_mean, coverage_variance, dataset, genotype_sum, true_homo_prop),
#     aes_(
#       x = as.name(new.col.names[2]),
#       y = ~true_homo_prop
#     )
#   ) +
#   geom_boxplot(
#     lwd = 1 / .pt,
#     fatten = 1.2 / .pt,
#     alpha = 1.0,
#     outlier.alpha = 1.0,
#     outlier.size = dot.size,
#     outlier.shape = NA,
#     aes_(
#       fill = ~coverage_variance,
#       color = ~coverage_variance
#     )
#   ) +
#   scale_fill_manual(
#     name = "Data quality",
#     breaks = names(data_quality_legend[[3]]),
#     labels = data_quality_legend[[1]],
#     values = data_quality_legend[[3]]
#   ) +
#   scale_color_manual(
#     name = "Data quality",
#     breaks = names(data_quality_legend[[2]]),
#     labels = data_quality_legend[[1]],
#     values = data_quality_legend[[2]]
#   ) +
#   geom_point(
#     position = position_jitterdodge(),
#     aes_(color = ~coverage_variance),
#     size = dot.size,
#     alpha = 0.5
#   ) +
#   facet_wrap(
#     ~cell_num,
#     labeller = labeller(
#       cell_num = cell.num.labels
#     )
#   ) +
#   scale_x_discrete(
#     breaks = as.numeric(levels(tree.info[[new.col.names[2]]])),
#     label = scientific
#   ) +
#   ylim(0, 0.015) +
#   labs(
#     y = "Porportion of true double mutant genotypes",
#     x = "Mutation rate"
#   ) +
#   theme_bw() +
#   theme(
#     text = element_text(size = 7),
#     axis.text = element_text(size = 7),
#     legend.text = element_text(size = 6),
#     legend.title = element_text(size = 6),
#     legend.title.align = 0.5,
#     legend.box.spacing = unit(0.3, "mm"),
#     legend.position = "bottom",
#     legend.key.size = unit(4, 'mm'),
#     strip.text = strip.text,
#     strip.background = element_blank(),
#     panel.grid = element_blank()
#   )

# Proportions of sites with parallel mutations.
# prop.parallel.mu.sites.plot <-
#   ggplot(
#     parallel.mut,
#     aes_(
#       x = as.name(new.col.names[2]),
#       y = ~prop_parallel_mu_sites
#     )
#   ) +
#   geom_boxplot(
#     lwd = 1 / .pt,
#     fatten = 1.2 / .pt,
#     alpha = 1.0,
#     outlier.alpha = 1.0,
#     outlier.size = dot.size,
#     outlier.shape = NA,
#     aes_(
#       fill = ~coverage_variance,
#       color = ~coverage_variance
#     )
#   ) +
#   scale_fill_manual(
#     name = "Data quality",
#     breaks = names(data_quality_legend[[3]]),
#     labels = data_quality_legend[[1]],
#     values = data_quality_legend[[3]]
#   ) +
#   scale_color_manual(
#     name = "Data quality",
#     breaks = names(data_quality_legend[[2]]),
#     labels = data_quality_legend[[1]],
#     values = data_quality_legend[[2]]
#   ) +
#   geom_point(
#     position = position_jitterdodge(),
#     aes_(color = ~coverage_variance),
#     size = dot.size,
#     alpha = 0.5
#   ) +
#   facet_wrap(
#     ~cell_num,
#     labeller = labeller(
#       # coverage_variance = allelic.raw.var,
#       cell_num = cell.num.labels
#     )
#   ) +
#   scale_x_discrete(
#     breaks = as.numeric(levels(tree.info[[new.col.names[2]]])),
#     label = scientific
#   ) +
#   ylim(0, 0.015) +
#   labs(
#     y = "Porportion of true double mutant genotypes",
#     x = "Mutation rate"
#   ) +
#   theme_bw() +
#   theme(
#     text = element_text(size = 7),
#     axis.text = element_text(size = 7),
#     legend.text = element_text(size = 6),
#     legend.title = element_text(size = 6),
#     legend.title.align = 0.5,
#     legend.box.spacing = unit(0.3, "mm"),
#     legend.position = "bottom",
#     legend.key.size = unit(4, 'mm'),
#     strip.text = strip.text,
#     strip.background = element_blank(),
#     panel.grid = element_blank()
#   )

# Parameters: ado rate
ado.rate.plot <- ggplot(
  param.info[which(!is.na(param.info$ado_rate)),],
  aes_(
    x = as.name(new.col.names[2]),
    y = ~ado_rate
  )
) +
  geom_boxplot(
    lwd = 1 / .pt,
    fatten = 1.2 / .pt,
    alpha = 1.0,
    outlier.alpha = 1.0,
    outlier.size = dot.size,
    outlier.shape = NA,
    aes_(
      fill = ~tool,
      color = ~tool
    )
  ) +
  scale_fill_manual(
    name = "Method",
    breaks = param.info.pretiffied[[1]],
    values = param.info.pretiffied[[3]]
  ) +
  scale_color_manual(
    name = "Method",
    breaks = param.info.pretiffied[[1]],
    values = param.info.pretiffied[[2]]
  ) +
  geom_point(
    position = position_jitterdodge(),
    aes_(color = ~tool),
    size = dot.size,
    alpha = 0.5
  ) +
  facet_grid(
    ~coverage_variance ~ ~cell_num,
    labeller = labeller(
      coverage_variance = allelic.raw.var,
      cell_num = cell.num.labels
    )
  ) +
  scale_x_discrete(
    breaks = as.numeric(levels(param.info[[new.col.names[2]]])),
    labels = scientific
  ) +
  ylim(0, 0.55) +
  geom_hline(
    yintercept = true.ado.rate,
    linetype="dashed",
    color = "gray43",
    size = 1 / .pt
  ) +
  labs(
    y = "ADO rate",
    x = "Mutation rate"
  ) +
  theme_bw() +
  theme(
    text = element_text(size = 7),
    axis.text = element_text(size = 7),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 6),
    legend.title.align = 0.5,
    legend.box.spacing = unit(0.3, "mm"),
    legend.position = "bottom",
    legend.key.size = unit(4, 'mm'),
    strip.text = strip.text,
    strip.background = element_blank(),
    panel.grid = element_blank()
  )

# Parameters: sequencing error rate
eff.seq.err.rate.plot <- ggplot(
  param.info[which(!is.na(param.info$eff_seq_err_rate)),],
  aes_(
    x = as.name(new.col.names[2]),
    y = ~eff_seq_err_rate
  )
) +
  geom_boxplot(
    lwd = 1 / .pt,
    fatten = 1.2 / .pt,
    alpha = 1.0,
    outlier.alpha = 1.0,
    outlier.size = dot.size,
    outlier.shape = NA,
    aes_(
      fill = ~tool,
      color = ~tool
    )
  ) +
  scale_fill_manual(
    name = "Method",
    breaks = param.info.pretiffied[[1]],
    values = param.info.pretiffied[[3]]
  ) +
  scale_color_manual(
    name = "Method",
    breaks = param.info.pretiffied[[1]],
    values = param.info.pretiffied[[2]]
  ) +
  geom_point(
    position = position_jitterdodge(),
    aes_(color = ~tool),
    size = dot.size,
    alpha = 0.5
  ) +
  facet_grid(
    ~coverage_variance ~ ~cell_num,
    labeller = labeller(
      coverage_variance = allelic.raw.var,
      cell_num = cell.num.labels
    )
  ) +
  scale_x_discrete(
    breaks = as.numeric(levels(param.info[[new.col.names[2]]])),
    labels = scientific
  ) +
  geom_hline(
    yintercept = true.eff.seq.err.rate,
    linetype="dashed",
    color = "gray43",
    size = 1 / .pt
  ) +
  labs(
    y = "Effective sequencing error rate",
    x = "Mutation rate"
  ) +
  theme_bw() +
  theme(
    text = element_text(size = 7),
    axis.text = element_text(size = 7),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 6),
    legend.title.align = 0.5,
    legend.box.spacing = unit(0.3, "mm"),
    legend.position = "bottom",
    legend.key.size = unit(4, 'mm'),
    strip.text = strip.text,
    strip.text.y = element_blank(),
    strip.background = element_blank(),
    panel.grid = element_blank()
  )

# Parameters: wildtype overdispersion
wildtype.overdispersion.plot <- ggplot(
  param.info[which(!is.na(param.info$wild_overdispersion)),],
  aes_(
    x = as.name(new.col.names[2]),
    y = ~wild_overdispersion
  )
) +
  geom_boxplot(
    lwd = 1 / .pt,
    fatten = 1.2 / .pt,
    alpha = 1.0,
    outlier.alpha = 1.0,
    outlier.size = dot.size,
    outlier.shape = NA,
    aes_(
      fill = ~tool,
      color = ~tool
    )
  ) +
  scale_fill_manual(
    name = "Method",
    breaks = param.info.pretiffied[[1]],
    values = param.info.pretiffied[[3]]
  ) +
  scale_color_manual(
    name = "Method",
    breaks = param.info.pretiffied[[1]],
    values = param.info.pretiffied[[2]]
  ) +
  geom_point(
    position = position_jitterdodge(),
    aes_(color = ~tool),
    size = dot.size,
    alpha = 0.5
  ) +
  facet_grid(
    ~coverage_variance ~ ~cell_num,
    labeller = labeller(
      coverage_variance = allelic.raw.var,
      cell_num = cell.num.labels
    )
  ) +
  scale_x_discrete(
    breaks = as.numeric(levels(param.info[[new.col.names[2]]])),
    labels = scientific
  ) +
  ylim(NA, 120) +
  geom_hline(
    yintercept = true.wildtype.overdispersion,
    linetype="dashed",
    color = "gray43",
    size = 1 / .pt
  ) +
  labs(
    y = "Wildtype overdispersion",
    x = "Mutation rate"
  ) +
  theme_bw() +
  theme(
    text = element_text(size = 7),
    axis.text = element_text(size = 7),
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 6),
    legend.title.align = 0.5,
    legend.box.spacing = unit(0.3, "mm"),
    legend.position = "bottom",
    legend.key.size = unit(4, 'mm'),
    strip.text = strip.text,
    strip.text.x = element_blank(),
    strip.text.y = element_blank(),
    strip.background = element_blank(),
    panel.grid = element_blank()
  )

# Parameters: alternative overdispersion
alt.overdispersion.plot <- ggplot(
  param.info[which(!is.na(param.info$alternative_overdispersion)),],
  aes_(
    x = as.name(new.col.names[2]),
    y = ~alternative_overdispersion
  )
) +
  geom_boxplot(
    lwd = 1 / .pt,
    fatten = 1.2 / .pt,
    alpha = 1.0,
    outlier.alpha = 1.0,
    outlier.size = dot.size,
    outlier.shape = NA,
    aes_(
      fill = ~tool,
      color = ~tool
    )
  ) +

  scale_fill_manual(
    name = "Method",
    breaks = param.info.pretiffied[[1]],
    values = param.info.pretiffied[[3]]
  ) +
  scale_color_manual(
    name = "Method",
    breaks = param.info.pretiffied[[1]],
    values = param.info.pretiffied[[2]]
  ) +
  geom_point(
    position = position_jitterdodge(),
    aes_(color = ~tool),
    size = dot.size,
    alpha = 0.5
  ) +
  facet_grid(
    ~coverage_variance ~ ~cell_num,
    labeller = labeller(
      coverage_variance = allelic.raw.var,
      cell_num = cell.num.labels
    )
  ) +
  scale_x_discrete(
    breaks = as.numeric(levels(param.info[[new.col.names[2]]])),
    labels = scientific
  ) +
  ylim(2.0, 2.8) +
  geom_hline(
    yintercept = true.alt.overdispersion,
    linetype="dashed",
    color = "gray43",
    size = 1 / .pt
  ) +
  labs(
    y = "Alternative overdispersion",
    x = "Mutation rate"
  ) +
  theme_bw() +
  theme(
    text = element_text(size = 7),
    axis.text = element_text(size = 7),
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 6),
    legend.title.align = 0.5,
    legend.box.spacing = unit(0.3, "mm"),
    legend.position = "bottom",
    legend.key.size = unit(4, 'mm'),
    strip.text = strip.text,
    strip.text.x = element_blank(),
    strip.background = element_blank(),
    panel.grid = element_blank()
  )


# Variant calling: mutation plot
# recall
recall.plot <- ggplot(
  var.info[which(!is.na(var.info$recall)),],
  aes_(
    x = as.name(new.col.names[2]),
    y = ~recall
  )
) +
  geom_boxplot(
    lwd = 1 / .pt,
    fatten = 1.2 / .pt,
    alpha = 1.0,
    outlier.alpha = 1.0,
    outlier.size = dot.size,
    outlier.shape = NA,
    aes_(
      fill = ~tool,
      color = ~tool
    )
  ) +
  scale_fill_manual(
    name = "Method",
    breaks = var.info.pretiffied[[1]],
    values = var.info.pretiffied[[3]]
  ) +
  scale_color_manual(
    name = "Method",
    breaks = var.info.pretiffied[[1]],
    values = var.info.pretiffied[[2]]
  ) +
  geom_point(
    position = position_jitterdodge(),
    aes_(color = ~tool),
    size = dot.size,
    alpha = 0.5
  ) +
  facet_grid(
    ~coverage_variance ~ ~cell_num,
    labeller = labeller(
      coverage_variance = allelic.raw.var,
      cell_num = cell.num.labels
    )
  ) +
  scale_x_discrete(
    breaks = as.numeric(levels(var.info[[new.col.names[2]]])),
    label = scientific
  ) +
  ylim(0.6, 1.0) +
  geom_hline(
    yintercept = 1,
    linetype="dashed",
    color = "gray43",
    size = 1 / .pt
  ) +
  labs(
    y = "Recall",
    x = "Mutation rate"
  ) +
  theme_bw() +
  theme(
    text = element_text(size = 7),
    axis.text = element_text(size = 7),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 6),
    legend.title.align = 0.5,
    legend.box.spacing = unit(0.3, "mm"),
    legend.position = "bottom",
    legend.key.size = unit(4, 'mm'),
    strip.text = strip.text,
    strip.text.y = element_blank(),
    strip.background = element_blank(),
    panel.grid = element_blank()
  )

# precision
precision.plot <- ggplot(
  var.info[which(!is.na(var.info$precision)),],
  aes_(
    x = as.name(new.col.names[2]),
    y = ~precision
  )
) +
  geom_boxplot(
    lwd = 1 / .pt,
    fatten = 1.2 / .pt,
    alpha = 1.0,
    outlier.alpha = 1.0,
    outlier.size = dot.size,
    outlier.shape = NA,
    aes_(
      fill = ~tool,
      color = ~tool
    )
  ) +
  scale_fill_manual(
    name = "Method",
    breaks = var.info.pretiffied[[1]],
    values = var.info.pretiffied[[3]]
  ) +
  scale_color_manual(
    name = "Method",
    breaks = var.info.pretiffied[[1]],
    values = var.info.pretiffied[[2]]
  ) +
  geom_point(
    position = position_jitterdodge(),
    aes_(color = ~tool),
    size = dot.size,
    alpha = 0.5
  ) +
  facet_grid(
    ~coverage_variance ~ ~cell_num,
    labeller = labeller(
      coverage_variance = allelic.raw.var,
      cell_num = cell.num.labels
    )
  ) +
  scale_x_discrete(
    breaks = as.numeric(levels(var.info[[new.col.names[2]]])),
    label = scientific
  ) +
  ylim(0.55, 1.0) +
  geom_hline(
    yintercept = 1,
    linetype="dashed",
    color = "gray43",
    size = 1 / .pt
  ) +
  labs(
    y = "Precision",
    x = "Mutation rate"
  ) +
  theme_bw() +
  theme(
    text = element_text(size = 7),
    axis.text = element_text(size = 7),
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 6),
    legend.title.align = 0.5,
    legend.box.spacing = unit(0.3, "mm"),
    legend.position = "bottom",
    legend.key.size = unit(4, 'mm'),
    strip.text = strip.text,
    strip.text.x = element_blank(),
    strip.text.y = element_blank(),
    strip.background = element_blank(),
    panel.grid = element_blank()
  )

# f1 score
f1.score.plot <- ggplot(
  var.info[which(!is.na(var.info$f1_score)),],
  aes_(
    x = as.name(new.col.names[2]),
    y = ~f1_score
  )
) +
  geom_boxplot(
    lwd = 1 / .pt,
    fatten = 1.2 / .pt,
    alpha = 1.0,
    outlier.alpha = 1.0,
    outlier.size = dot.size,
    outlier.shape = NA,
    aes_(
      fill = ~tool,
      color = ~tool
    )
  ) +
  scale_fill_manual(
    name = "Method",
    breaks = var.info.pretiffied[[1]],
    values = var.info.pretiffied[[3]]
  ) +
  scale_color_manual(
    name = "Method",
    breaks = var.info.pretiffied[[1]],
    values = var.info.pretiffied[[2]]
  ) +
  geom_point(
    position = position_jitterdodge(),
    aes_(color = ~tool),
    size = dot.size,
    alpha = 0.5
  ) +
  facet_grid(
    ~coverage_variance ~ ~cell_num,
    labeller = labeller(
      coverage_variance = allelic.raw.var,
      cell_num = cell.num.labels
    )
  ) +
  scale_x_discrete(
    breaks = as.numeric(levels(fsa.var.info[[new.col.names[2]]])),
    label = scientific
  ) +
  ylim(0.7, 1.0) +
  geom_hline(
    yintercept = 1,
    linetype="dashed",
    color = "gray43",
    size = 1 / .pt
  ) +
  labs(
    y = "F1 score",
    x = "Mutation rate"
  ) +
  theme_bw() +
  theme(
    text = element_text(size = 7),
    axis.text = element_text(size = 7),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 6),
    legend.title.align = 0.5,
    legend.box.spacing = unit(0.3, "mm"),
    legend.position = "bottom",
    legend.key.size = unit(4, 'mm'),
    strip.text = strip.text,
    strip.text.x = element_blank(),
    strip.background = element_blank(),
    panel.grid = element_blank()
  )

# fallout
fallout.plot <- ggplot(
  var.info[which(!is.na(var.info$fall_out)),],
  aes_(
    x = as.name(new.col.names[2]),
    y = ~fall_out
  )
) +
  geom_boxplot(
    lwd = 1 / .pt,
    fatten = 1.2 / .pt,
    alpha = 1.0,
    outlier.alpha = 1.0,
    outlier.size = dot.size,
    outlier.shape = NA,
    aes_(
      fill = ~tool,
      color = ~tool
    )
  ) +
  scale_fill_manual(
    name = "Method",
    breaks = var.info.pretiffied[[1]],
    values = var.info.pretiffied[[3]]
  ) +
  scale_color_manual(
    name = "Method",
    breaks = var.info.pretiffied[[1]],
    values = var.info.pretiffied[[2]]
  ) +
  geom_point(
    position = position_jitterdodge(),
    aes_(color = ~tool),
    size = dot.size,
    alpha = 0.5
  ) +
  facet_grid(
    ~coverage_variance ~ ~cell_num,
    labeller = labeller(
      coverage_variance = allelic.raw.var,
      cell_num = cell.num.labels
    )
  ) +
  scale_x_discrete(
    breaks = as.numeric(levels(fsa.var.info[[new.col.names[2]]])),
    label = scientific
  ) +
  ylim(0, NA) +
  geom_hline(
    yintercept = 0,
    linetype="dashed",
    color = "gray43",
    size = 1 / .pt
  ) +
  labs(
    y = "False positive rate",
    x = "Mutation rate"
  ) +
  theme_bw() +
  theme(
    text = element_text(size = 7),
    axis.text = element_text(size = 7),
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 6),
    legend.title.align = 0.5,
    legend.box.spacing = unit(0.3, "mm"),
    legend.position = "bottom",
    legend.key.size = unit(4, 'mm'),
    strip.text = strip.text,
    strip.text.x = element_blank(),
    strip.background = element_blank(),
    panel.grid = element_blank()
  )

# Variant calling: heterozygous mutation plot
# recall
hetero.recall.plot <- ggplot(
  var.info[which(!is.na(var.info$recall_hetero_mu)),],
  aes_(
    x = as.name(new.col.names[2]),
    y = ~recall_hetero_mu
  )
) +
  geom_boxplot(
    lwd = 1 / .pt,
    fatten = 1.2 / .pt,
    alpha = 1.0,
    outlier.alpha = 1.0,
    outlier.size = dot.size,
    outlier.shape = NA,
    aes_(
      fill = ~tool,
      color = ~tool
    )
  ) +
  scale_fill_manual(
    name = "Method",
    breaks = common.legend[[1]],
    values = common.legend[[3]]
  ) +
  scale_color_manual(
    name = "Method",
    breaks = common.legend[[1]],
    values = common.legend[[2]]
  ) +
  geom_point(
    position = position_jitterdodge(),
    aes_(color = ~tool),
    size = dot.size,
    alpha = 0.5
  ) +
  facet_grid(
    ~coverage_variance ~ ~cell_num,
    labeller = labeller(
      coverage_variance = allelic.raw.var,
      cell_num = cell.num.labels
    )
  ) +
  scale_x_discrete(
    breaks = as.numeric(levels(var.info[[new.col.names[2]]])),
    label = scientific
  ) +
  ylim(0.5, 1.0) +
  geom_hline(
    yintercept = 1,
    linetype="dashed",
    color = "gray43",
    size = 1 / .pt
  ) +
  labs(
    y = "Recall",
    x = "Mutation rate"
  ) +
  theme_bw() +
  theme(
    text = element_text(size = 7),
    axis.text = element_text(size = 7),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 6),
    legend.title.align = 0.5,
    legend.box.spacing = unit(0.3, "mm"),
    legend.position = "bottom",
    legend.key.size = unit(4, 'mm'),
    strip.text = strip.text,
    strip.text.x = element_blank(),
    strip.text.y = element_blank(),
    strip.background = element_blank(),
    panel.grid = element_blank()
  )

# precision
hetero.precision.plot <- ggplot(
  var.info[which(!is.na(var.info$precision_hetero_mu)),],
  aes_(
    x = as.name(new.col.names[2]),
    y = ~precision_hetero_mu
  )
) +
  geom_boxplot(
    lwd = 1 / .pt,
    fatten = 1.2 / .pt,
    alpha = 1.0,
    outlier.alpha = 1.0,
    outlier.size = dot.size,
    outlier.shape = NA,
    aes_(
      fill = ~tool,
      color = ~tool
    )
  ) +
  scale_fill_manual(
    name = "Method",
    breaks = common.legend[[1]],
    values = common.legend[[3]]
  ) +
  scale_color_manual(
    name = "Method",
    breaks = common.legend[[1]],
    values = common.legend[[2]]
  ) +
  geom_point(
    position = position_jitterdodge(),
    aes_(color = ~tool),
    size = dot.size,
    alpha = 0.5
  ) +
  facet_grid(
    ~coverage_variance ~ ~cell_num,
    labeller = labeller(
      coverage_variance = allelic.raw.var,
      cell_num = cell.num.labels
    )
  ) +
  scale_x_discrete(
    breaks = as.numeric(levels(var.info[[new.col.names[2]]])),
    label = scientific
  ) +
  ylim(0.55, 1.0) +
  geom_hline(
    yintercept = 1,
    linetype="dashed",
    color = "gray43",
    size = 1 / .pt
  ) +
  labs(
    y = "Precision",
    x = "Mutation rate"
  ) +
  theme_bw() +
  theme(
    text = element_text(size = 7),
    axis.text = element_text(size = 7),
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 6),
    legend.title.align = 0.5,
    legend.box.spacing = unit(0.3, "mm"),
    legend.position = "bottom",
    legend.key.size = unit(4, 'mm'),
    strip.text = strip.text,
    strip.text.x = element_blank(),
    strip.text.y = element_blank(),
    strip.background = element_blank(),
    panel.grid = element_blank(),
    plot.margin = unit(c(0, 0, 0, 0),"mm")
  )

# f1 score
hetero.f1.score.plot <- ggplot(
  var.info[which(!is.na(var.info$f1_score_hetero_mu)),],
  aes_(
    x = as.name(new.col.names[2]),
    y = ~f1_score_hetero_mu
  )
) +
  geom_boxplot(
    lwd = 1 / .pt,
    fatten = 1.2 / .pt,
    alpha = 1.0,
    outlier.alpha = 1.0,
    outlier.size = dot.size,
    outlier.shape = NA,
    aes_(
      fill = ~tool,
      color = ~tool
    )
  ) +
  scale_fill_manual(
    name = "Method",
    breaks = common.legend[[1]],
    values = common.legend[[3]]
  ) +
  scale_color_manual(
    name = "Method",
    breaks = common.legend[[1]],
    values = common.legend[[2]]
  ) +
  geom_point(
    position = position_jitterdodge(),
    aes_(color = ~tool),
    size = dot.size,
    alpha = 0.5
  ) +
  facet_grid(
    ~coverage_variance ~ ~cell_num,
    labeller = labeller(
      coverage_variance = allelic.raw.var,
      cell_num = cell.num.labels
    )
  ) +
  scale_x_discrete(
    breaks = as.numeric(levels(var.info[[new.col.names[2]]])),
    label = scientific
  ) +
  ylim(0.65, 1.0) +
  geom_hline(
    yintercept = 1,
    linetype="dashed",
    color = "gray43",
    size = 1 / .pt
  ) +
  labs(
    y = "F1 score",
    x = "Mutation rate"
  ) +
  theme_bw() +
  theme(
    text = element_text(size = 7),
    axis.text = element_text(size = 7),
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 6),
    legend.title.align = 0.5,
    legend.box.spacing = unit(0.3, "mm"),
    legend.key.size = unit(4, 'mm'),
    legend.position = "bottom",
    strip.text = strip.text,
    strip.text.x = element_blank(),
    strip.text.y = element_blank(),
    strip.background = element_blank(),
    panel.grid = element_blank()
  )

# fallout
hetero.fallout.plot <- ggplot(
  var.info[which(!is.na(var.info$fall_out_hetero_mu)),],
  aes_(
    x = as.name(new.col.names[2]),
    y = ~fall_out_hetero_mu
  )
) +
  geom_boxplot(
    lwd = 1 / .pt,
    fatten = 1.2 / .pt,
    alpha = 1.0,
    outlier.alpha = 1.0,
    outlier.size = dot.size,
    outlier.shape = NA,
    aes_(
      fill = ~tool,
      color = ~tool
    )
  ) +
  scale_fill_manual(
    name = "Method",
    breaks = var.info.pretiffied[[1]],
    values = var.info.pretiffied[[3]]
  ) +
  scale_color_manual(
    name = "Method",
    breaks = var.info.pretiffied[[1]],
    values = var.info.pretiffied[[2]]
  ) +
  geom_point(
    position = position_jitterdodge(),
    aes_(color = ~tool),
    size = dot.size,
    alpha = 0.5
  ) +
  facet_grid(
    ~coverage_variance ~ ~cell_num,
    labeller = labeller(
      coverage_variance = allelic.raw.var,
      cell_num = cell.num.labels
    )
  ) +
  scale_x_discrete(
    breaks = as.numeric(levels(var.info[[new.col.names[2]]])),
    label = scientific
  ) +
  ylim(0, 0.15) +
  geom_hline(
    yintercept = 0,
    linetype="dashed",
    color = "gray43",
    size = 1 / .pt
  ) +
  labs(
    y = "False positive rate",
    x = "Mutation rate"
  ) +
  theme_bw() +
  theme(
    text = element_text(size = 7),
    axis.text = element_text(size = 7),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 6),
    legend.title.align = 0.5,
    legend.box.spacing = unit(0.3, "mm"),
    legend.position = "bottom",
    legend.key.size = unit(4, 'mm'),
    strip.text = strip.text,
    strip.text.y = element_blank(),
    strip.background = element_blank(),
    panel.grid = element_blank()
  )

# Variant calling: homozygous mutation plot
# recall
homo.recall.plot <- ggplot(
  fsa.var.info[which(!is.na(fsa.var.info$recall_homo_mu) & fsa.var.info$tool != "SCIPhI"),],
  aes_(
    x = as.name(new.col.names[2]),
    y = ~recall_homo_mu
  )
) +
  geom_boxplot(
    lwd = 1 / .pt,
    fatten = 1.2 / .pt,
    alpha = 1.0,
    outlier.alpha = 1.0,
    outlier.size = dot.size,
    outlier.shape = NA,
    aes_(
      fill = ~tool,
      color = ~tool
    )
  ) +
  scale_fill_manual(
    name = "Method",
    breaks = common.legend[[1]],
    values = common.legend[[3]]
  ) +
  scale_color_manual(
    name = "Method",
    breaks = common.legend[[1]],
    values = common.legend[[2]]
  ) +
  geom_point(
    position = position_jitterdodge(),
    aes_(color = ~tool),
    size = dot.size,
    alpha = 0.5
  ) +
  facet_grid(
    ~coverage_variance ~ ~cell_num,
    labeller = labeller(
      coverage_variance = allelic.raw.var,
      cell_num = cell.num.labels
    )
  ) +
  scale_x_discrete(
    breaks = as.numeric(levels(fsa.var.info[[new.col.names[2]]])),
    label = scientific
  ) +
  ylim(0, 1) +
  geom_hline(
    yintercept = 1,
    linetype="dashed",
    color = "gray43",
    size = 1 / .pt
  ) +
  labs(
    y = "Recall",
    x = "Mutation rate"
  ) +
  theme_bw() +
  theme(
    text = element_text(size = 7),
    axis.text = element_text(size = 7),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 6),
    legend.title.align = 0.5,
    legend.box.spacing = unit(0.3, "mm"),
    legend.position = "bottom",
    legend.key.size = unit(4, 'mm'),
    strip.text = strip.text,
    strip.text.x = element_blank(),
    strip.background = element_blank(),
    panel.grid = element_blank()
  )

# precision
homo.precision.plot <- ggplot(
  fsa.var.info[which(!is.na(fsa.var.info$precision_homo_mu) & fsa.var.info$tool != "SCIPhI"),],
  aes_(
    x = as.name(new.col.names[2]),
    y = ~precision_homo_mu
  )
) +
  geom_boxplot(
    lwd = 1 / .pt,
    fatten = 1.2 / .pt,
    alpha = 1.0,
    outlier.alpha = 1.0,
    outlier.size = dot.size,
    outlier.shape = NA,
    aes_(
      fill = ~tool,
      color = ~tool
    )
  ) +
  scale_fill_manual(
    name = "Method",
    breaks = common.legend[[1]],
    values = common.legend[[3]]
  ) +
  scale_color_manual(
    name = "Method",
    breaks = common.legend[[1]],
    values = common.legend[[2]]
  ) +
  geom_point(
    position = position_jitterdodge(),
    aes_(color = ~tool),
    size = dot.size,
    alpha = 0.5
  ) +
  facet_grid(
    ~coverage_variance ~ ~cell_num,
    labeller = labeller(
      coverage_variance = allelic.raw.var,
      cell_num = cell.num.labels
    )
  ) +
  scale_x_discrete(
    breaks = as.numeric(levels(fsa.var.info[[new.col.names[2]]])),
    label = scientific
  ) +
  ylim(0, 1) +
  geom_hline(
    yintercept = 1,
    linetype="dashed",
    color = "gray43",
    size = 1 / .pt
  ) +
  labs(
    y = "Precision",
    x = "Mutation rate"
  ) +
  theme_bw() +
  theme(
    text = element_text(size = 7),
    axis.text = element_text(size = 7),
    axis.title.y = element_blank(),
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 6),
    legend.title.align = 0.5,
    legend.box.spacing = unit(0.3, "mm"),
    legend.position = "bottom",
    legend.key.size = unit(4, 'mm'),
    strip.text = strip.text,
    strip.text.x = element_blank(),
    strip.background = element_blank(),
    panel.grid = element_blank(),
    plot.margin = unit(c(0, 0, 0, 0),"mm")
  )

# f1 score
homo.f1.score.plot <- ggplot(
 fsa.var.info[which(!is.na(fsa.var.info$f1_score_homo_mu) & fsa.var.info$tool != "SCIPhI"),],
 aes_(
   x = as.name(new.col.names[2]),
   y = ~f1_score_homo_mu
 )
) +
  geom_boxplot(
    lwd = 1 / .pt,
    fatten = 1.2 / .pt,
    alpha = 1.0,
    outlier.alpha = 1.0,
    outlier.size = dot.size,
    outlier.shape = NA,
    aes_(
      fill = ~tool,
      color = ~tool
    )
  ) +
  scale_fill_manual(
    name = "Method",
    breaks = common.legend[[1]],
    values = common.legend[[3]]
  ) +
  scale_color_manual(
    name = "Method",
    breaks = common.legend[[1]],
    values = common.legend[[2]]
  ) +
  geom_point(
    position = position_jitterdodge(),
    aes_(color = ~tool),
    size = dot.size,
    alpha = 0.5
  ) +
  facet_grid(
    ~coverage_variance ~ ~cell_num,
    labeller = labeller(
      coverage_variance = allelic.raw.var,
      cell_num = cell.num.labels
    )
  ) +
  scale_x_discrete(
    breaks = as.numeric(levels(fsa.var.info[[new.col.names[2]]])),
    label = scientific
  ) +
  ylim(0, 1) +
  geom_hline(
    yintercept = 1,
    linetype="dashed",
    color = "gray43",
    size = 1 / .pt
  ) +
  labs(
    y = "F1 score",
    x = "Mutation rate"
  ) +
  theme_bw() +
  theme(
    text = element_text(size = 7),
    axis.text = element_text(size = 7),
    axis.title.y = element_blank(),
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 6),
    legend.title.align = 0.5,
    legend.box.spacing = unit(0.3, "mm"),
    legend.position = "none",
    legend.key.size = unit(4, 'mm'),
    strip.text = strip.text,
    strip.text.x = element_blank(),
    strip.background = element_blank(),
    panel.grid = element_blank()
  )

# fallout
homo.fallout.plot <- ggplot(
  fsa.var.info[which(!is.na(fsa.var.info$fall_out_homo_mu) & fsa.var.info$tool != "SCIPhI"),],
  aes_(
    x = as.name(new.col.names[2]),
    y = ~fall_out_homo_mu
  )
) +
  geom_boxplot(
    lwd = 1 / .pt,
    fatten = 1.2 / .pt,
    alpha = 1.0,
    outlier.alpha = 1.0,
    outlier.size = dot.size,
    outlier.shape = NA,
    aes_(
      fill = ~tool,
      color = ~tool
    )
  ) +
  scale_fill_manual(
    name = "Method",
    breaks = fsa.var.info.pretiffied[[1]],
    values = fsa.var.info.pretiffied[[3]]
  ) +
  scale_color_manual(
    name = "Method",
    breaks = fsa.var.info.pretiffied[[1]],
    values = fsa.var.info.pretiffied[[2]]
  ) +
  geom_point(
    position = position_jitterdodge(),
    aes_(color = ~tool),
    size = dot.size,
    alpha = 0.5
  ) +
  facet_grid(
    ~coverage_variance ~ ~cell_num,
    labeller = labeller(
      coverage_variance = allelic.raw.var,
      cell_num = cell.num.labels
    )
  ) +
  scale_x_discrete(
    breaks = as.numeric(levels(fsa.var.info[[new.col.names[2]]])),
    label = scientific
  ) +
  ylim(0, 0.05) +
  geom_hline(
    yintercept = 0,
    linetype="dashed",
    color = "gray43",
    size = 1 / .pt
  ) +
  labs(
    y = "False positive rate",
    x = "Mutation rate"
  ) +
  theme_bw() +
  theme(
    text = element_text(size = 7),
    axis.text = element_text(size = 7),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 6),
    legend.title.align = 0.5,
    legend.box.spacing = unit(0.3, "mm"),
    legend.position = "bottom",
    legend.key.size = unit(4, 'mm'),
    strip.text = strip.text,
    strip.background = element_blank(),
    panel.grid = element_blank()
  )

# Source analysis of false positives for heterozygous mutations: originally homozygous reference
true.homo.ref.as.hetero.mu.in.called.pos.plot <- ggplot(
  var.info,
  aes_(
    x = as.name(new.col.names[2]),
    y = ~prop_true_homo_ref_as_hetero_mu_in_called_pos
  )
) +
  geom_boxplot(
    lwd = 1 / .pt,
    fatten = 1.2 / .pt,
    alpha = 1.0,
    outlier.alpha = 1.0,
    outlier.size = dot.size,
    outlier.shape = NA,
    aes_(
      fill = ~tool,
      color = ~tool
    )
  ) +
  scale_fill_manual(
    name = "Method",
    breaks = var.info.pretiffied[[1]],
    values = var.info.pretiffied[[3]]
  ) +
  scale_color_manual(
    name = "Method",
    breaks = var.info.pretiffied[[1]],
    values = var.info.pretiffied[[2]]
  ) +
  geom_point(
    position = position_jitterdodge(),
    aes_(color = ~tool),
    size = dot.size,
    alpha = 0.5
  ) +
  facet_grid(
    ~coverage_variance ~ ~cell_num,
    labeller = labeller(
      coverage_variance = allelic.raw.var,
      cell_num = cell.num.labels
    )
  ) +
  scale_x_discrete(
    breaks = as.numeric(levels(var.info[[new.col.names[2]]])),
    label = scientific
  ) +
  geom_hline(
    yintercept = 0,
    linetype="dashed",
    color = "gray43",
    size = 1 / .pt
  ) +
  labs(
    y = "Proportion of true wildtype genotype in predicted single mutant genotype",
    x = "Mutation rate"
  ) +
  theme_bw() +
  theme(
    text = element_text(size = 7),
    axis.text = element_text(size = 7),
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 6),
    legend.title.align = 0.5,
    legend.box.spacing = unit(0.3, "mm"),
    legend.position = "bottom",
    legend.key.size = unit(4, 'mm'),
    strip.text = strip.text,
    strip.text.y = element_blank(),
    strip.background = element_blank(),
    panel.grid = element_blank()
  )

# Source analysis of false positives for heterozygous mutations: originally homozygous mutation
true.homo.mu.as.hetero.mu.in.called.pos.plot <- ggplot(
  var.info,
  aes_(
    x = as.name(new.col.names[2]),
    y = ~prop_true_homo_mu_as_hetero_mu_in_called_pos
  )
) +
  geom_boxplot(
    lwd = 1 / .pt,
    fatten = 1.2 / .pt,
    alpha = 1.0,
    outlier.alpha = 1.0,
    outlier.size = dot.size,
    outlier.shape = NA,
    aes_(
      fill = ~tool,
      color = ~tool
    )
  ) +
  scale_fill_manual(
    name = "Method",
    breaks = var.info.pretiffied[[1]],
    values = var.info.pretiffied[[3]]
  ) +
  scale_color_manual(
    name = "Method",
    breaks = var.info.pretiffied[[1]],
    values = var.info.pretiffied[[2]]
  ) +
  geom_point(
    position = position_jitterdodge(),
    aes_(color = ~tool),
    size = dot.size,
    alpha = 0.5
  ) +
  facet_grid(
    ~coverage_variance ~ ~cell_num,
    labeller = labeller(
      coverage_variance = allelic.raw.var,
      cell_num = cell.num.labels
    )
  ) +
  scale_x_discrete(
    breaks = as.numeric(levels(var.info[[new.col.names[2]]])),
    label = scientific
  ) +
  geom_hline(
    yintercept = 0,
    linetype="dashed",
    color = "gray43",
    size = 1 / .pt
  ) +
  labs(
    y = "Proportion of true double mutant genotypes in predicted single mutant genotype",
    x = "Mutation rate"
  ) +
  theme_bw() +
  theme(
    text = element_text(size = 7),
    axis.text = element_text(size = 7),
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 6),
    legend.title.align = 0.5,
    legend.box.spacing = unit(0.3, "mm"),
    legend.position = "bottom",
    legend.key.size = unit(4, 'mm'),
    strip.text = strip.text,
    strip.background = element_blank(),
    panel.grid = element_blank()
  )

# ADO calling: mutation plot
# recall
ado.recall.plot <- ggplot(
  ado.info[which(!is.na(ado.info$recall)),],
  aes_(
    x = as.name(new.col.names[2]),
    y = ~recall
  )
) +
  geom_boxplot(
    lwd = 1 / .pt,
    fatten = 1.2 / .pt,
    alpha = 1.0,
    outlier.alpha = 1.0,
    outlier.size = dot.size,
    outlier.shape = NA,
    aes_(
      fill = ~tool,
      color = ~tool
    )
  ) +
  scale_fill_manual(
    name = "Method",
    breaks = ado.info.pretiffied[[1]],
    values = ado.info.pretiffied[[3]]
  ) +
  scale_color_manual(
    name = "Method",
    breaks = ado.info.pretiffied[[1]],
    values = ado.info.pretiffied[[2]]
  ) +
  geom_point(
    position = position_jitterdodge(),
    aes_(color = ~tool),
    size = dot.size,
    alpha = 0.5
  ) +
  facet_grid(
    ~coverage_variance ~ ~cell_num,
    labeller = labeller(
      coverage_variance = allelic.raw.var,
      cell_num = cell.num.labels
    )
  ) +
  scale_x_discrete(
    breaks = as.numeric(levels(ado.info[[new.col.names[2]]])),
    label = scientific
  ) +
  ylim(0, 1.0) +
  geom_hline(
    yintercept = 1,
    linetype="dashed",
    color = "gray43",
    size = 1 / .pt
  ) +
  labs(
    y = "Recall",
    x = "Mutation rate"
  ) +
  theme_bw() +
  theme(
    text = element_text(size = 7),
    axis.text = element_text(size = 7),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 6),
    legend.title.align = 0.5,
    legend.box.spacing = unit(0.3, "mm"),
    legend.position = "bottom",
    legend.key.size = unit(4, 'mm'),
    strip.text = strip.text,
    strip.text.y = element_blank(),
    strip.background = element_blank(),
    panel.grid = element_blank()
  )

# precision
ado.precision.plot <- ggplot(
  ado.info[which(!is.na(ado.info$precision)),],
  aes_(
    x = as.name(new.col.names[2]),
    y = ~precision
  )
) +
  geom_boxplot(
    lwd = 1 / .pt,
    fatten = 1.2 / .pt,
    alpha = 1.0,
    outlier.alpha = 1.0,
    outlier.size = dot.size,
    outlier.shape = NA,
    aes_(
      fill = ~tool,
      color = ~tool
    )
  ) +
  scale_fill_manual(
    name = "Method",
    breaks = ado.info.pretiffied[[1]],
    values = ado.info.pretiffied[[3]]
  ) +
  scale_color_manual(
    name = "Method",
    breaks = ado.info.pretiffied[[1]],
    values = ado.info.pretiffied[[2]]
  ) +
  geom_point(
    position = position_jitterdodge(),
    aes_(color = ~tool),
    size = dot.size,
    alpha = 0.5
  ) +
  facet_grid(
    ~coverage_variance ~ ~cell_num,
    labeller = labeller(
      coverage_variance = allelic.raw.var,
      cell_num = cell.num.labels
    )
  ) +
  scale_x_discrete(
    breaks = as.numeric(levels(ado.info[[new.col.names[2]]])),
    label = scientific
  ) +
  ylim(0, 1.0) +
  geom_hline(
    yintercept = 1,
    linetype="dashed",
    color = "gray43",
    size = 1 / .pt
  ) +
  labs(
    y = "Precision",
    x = "Mutation rate"
  ) +
  theme_bw() +
  theme(
    text = element_text(size = 7),
    axis.text = element_text(size = 7),
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 6),
    legend.title.align = 0.5,
    legend.box.spacing = unit(0.3, "mm"),
    legend.position = "bottom",
    legend.key.size = unit(4, 'mm'),
    strip.text = strip.text,
    strip.text.x = element_blank(),
    strip.text.y = element_blank(),
    strip.background = element_blank(),
    panel.grid = element_blank()
  )

# f1 score
ado.f1.score.plot <- ggplot(
  ado.info[which(!is.na(ado.info$f1_score)),],
  aes_(
    x = as.name(new.col.names[2]),
    y = ~f1_score
  )
) +
  geom_boxplot(
    lwd = 1 / .pt,
    fatten = 1.2 / .pt,
    alpha = 1.0,
    outlier.alpha = 1.0,
    outlier.size = dot.size,
    outlier.shape = NA,
    aes_(
      fill = ~tool,
      color = ~tool
    )
  ) +
  scale_fill_manual(
    name = "Method",
    breaks = ado.info.pretiffied[[1]],
    values = ado.info.pretiffied[[3]]
  ) +
  scale_color_manual(
    name = "Method",
    breaks = ado.info.pretiffied[[1]],
    values = ado.info.pretiffied[[2]]
  ) +
  geom_point(
    position = position_jitterdodge(),
    aes_(color = ~tool),
    size = dot.size,
    alpha = 0.5
  ) +
  facet_grid(
    ~coverage_variance ~ ~cell_num,
    labeller = labeller(
      coverage_variance = allelic.raw.var,
      cell_num = cell.num.labels
    )
  ) +
  scale_x_discrete(
    breaks = as.numeric(levels(ado.info[[new.col.names[2]]])),
    label = scientific
  ) +
  ylim(0, 1.0) +
  geom_hline(
    yintercept = 1,
    linetype="dashed",
    color = "gray43",
    size = 1 / .pt
  ) +
  labs(
    y = "F1 score",
    x = "Mutation rate"
  ) +
  theme_bw() +
  theme(
    text = element_text(size = 7),
    axis.text = element_text(size = 7),
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 6),
    legend.title.align = 0.5,
    legend.box.spacing = unit(0.3, "mm"),
    legend.position = "bottom",
    legend.key.size = unit(4, 'mm'),
    strip.text = strip.text,
    strip.text.x = element_blank(),
    strip.background = element_blank(),
    panel.grid = element_blank()
  )

# fallout
ado.fallout.plot <- ggplot(
  ado.info[which(!is.na(ado.info$fall_out)),],
  aes_(
    x = as.name(new.col.names[2]),
    y = ~fall_out
  )
) +
  geom_boxplot(
    lwd = 1 / .pt,
    fatten = 1.2 / .pt,
    alpha = 1.0,
    outlier.alpha = 1.0,
    outlier.size = dot.size,
    outlier.shape = NA,
    aes_(
      fill = ~tool,
      color = ~tool
    )
  ) +
  scale_fill_manual(
    name = "Method",
    breaks = ado.info.pretiffied[[1]],
    values = ado.info.pretiffied[[3]]
  ) +
  scale_color_manual(
    name = "Method",
    breaks = ado.info.pretiffied[[1]],
    values = ado.info.pretiffied[[2]]
  ) +
  geom_point(
    position = position_jitterdodge(),
    aes_(color = ~tool),
    size = dot.size,
    alpha = 0.5
  ) +
  facet_grid(
    ~coverage_variance ~ ~cell_num,
    labeller = labeller(
      coverage_variance = allelic.raw.var,
      cell_num = cell.num.labels
    )
  ) +
  scale_x_discrete(
    breaks = as.numeric(levels(ado.info[[new.col.names[2]]])),
    label = scientific
  ) +
  ylim(0, 0.05) +
  geom_hline(
    yintercept = 0,
    linetype="dashed",
    color = "gray43",
    size = 1 / .pt
  ) +
  labs(
    y = "False positive rate",
    x = "Mutation rate"
  ) +
  theme_bw() +
  theme(
    text = element_text(size = 7),
    axis.text = element_text(size = 7),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 6),
    legend.title.align = 0.5,
    legend.box.spacing = unit(0.3, "mm"),
    legend.position = "bottom",
    legend.key.size = unit(4, 'mm'),
    strip.text = strip.text,
    strip.background = element_blank(),
    panel.grid = element_blank()
  )


##################################################
#                Generate figures                #
##################################################

# Figure 2
# tree structure
tree.plot <- ggarrange(
  rooted.branch.score.diff.plot,
  normalized.rf.dist.plot,
  ncol = 2L,
  align = "h",
  labels = c("a", "b"),
  font.label = list(face = "bold", size = 7),
  legend = "none"
  )

# recall and precision for heterozygous mutation genotype
hetero.plot <- ggarrange(
  hetero.recall.plot,
  hetero.precision.plot,
  nrow = 2L,
  align = "v",
  labels = c("c", "d"),
  vjust = c(1.5, 0.4),
  font.label = list(face = "bold", size = 7),
  legend = "none"
)

hetero.plot <- annotate_figure(
  hetero.plot,
  top = text_grob(
    label = "Single mutant genotype calling",
    size = 7,
    hjust = 0.4
  )
)

# recall and precision for homozygous mutation genotype
homo.plot <- ggarrange(
  homo.recall.plot,
  homo.precision.plot,
  nrow = 2L,
  align = "v",
  labels = c("e", "f"),
  vjust = c(1.5, 0.4),
  font.label = list(face = "bold", size = 7),
  legend = "none"
)

homo.plot <- annotate_figure(
  homo.plot,
  top = text_grob(
    label = "Double mutant genotype calling",
    size = 7
  )
)

# combine into figure 2
fig2.plot <- ggarrange(
  tree.plot,
  ggarrange(
    hetero.plot,
    homo.plot,
    ncol = 2L,
    align = "h"
  ),
  nrow = 2L,
  align = "h",
  heights = c(1, 2),
  common.legend = TRUE,
  legend = "bottom",
  legend.grob = get_legend(normalized.rf.dist.plot)
)

ggsave(
  paste0(output.prefix, "fig2.pdf"),
  plot = fig2.plot,
  device = "pdf",
  width = 183,
  height = 215,
  units = "mm",
  dpi = 300
)

# Figure S1
ggsave(
  paste0(output.prefix, "s1.pdf"),
  plot = cor.plot,
  device = "pdf",
  width = 89,
  height = 90,
  units = "mm",
  dpi = 300
)

# Figure S2
param.plot <- ggarrange(
  ggarrange(
    eff.seq.err.rate.plot,
    wildtype.overdispersion.plot,
    nrow = 2L,
    align = "v",
    labels = c("a", "c"),
    font.label = list(face = "bold", size = 7),
    legend = "none"
  ),
  ggarrange(
    ado.rate.plot,
    alt.overdispersion.plot,
    nrow = 2L,
    align = "v",
    labels = c("b", "d"),
    font.label = list(face = "bold", size = 7),
    legend = "none"
  ),
  ncol = 2L,
  align = "hv",
  common.legend = TRUE,
  legend = "bottom",
  legend.grob = get_legend(eff.seq.err.rate.plot)
)

ggsave(
  paste0(output.prefix, "s2.pdf"),
  plot = param.plot,
  device = "pdf",
  width = 183,
  height = 150,
  units = "mm",
  dpi = 300
)

# Figure S3
# F1 score and false positive rate for heterozygous mutation genotype
hetero.splot <- ggarrange(
  hetero.fallout.plot,
  hetero.f1.score.plot,
  nrow = 2L,
  align = "v",
  labels = c("a", "b"),
  hjust = 0,
  font.label = list(face = "bold", size = 7),
  legend = "none"
)

hetero.splot <- annotate_figure(
  hetero.splot,
  top = text_grob(
    label = "Single mutant genotype calling",
    size = 7
  )
)

# F1 score and false positive rate for homozygous mutation genotype
homo.splot <- ggarrange(
  homo.fallout.plot,
  homo.f1.score.plot,
  nrow = 2L,
  align = "v",
  labels = c("c", "d"),
  hjust = 0,
  font.label = list(face = "bold", size = 7),
  legend = "none"
)

homo.splot <- annotate_figure(
  homo.splot,
  top = text_grob(
    label = "Double mutant genotype calling",
    size = 7
  )
)

# combine into figure S3
figs3.plot <- ggarrange(
  hetero.splot,
  homo.splot,
  ncol = 2L,
  align = "h",
  common.legend = TRUE,
  legend = "bottom",
  legend.grob = get_legend(true.homo.ref.as.hetero.mu.in.called.pos.plot)
)

ggsave(
  paste0(output.prefix, "s3.pdf"),
  plot = figs3.plot,
  device = "pdf",
  width = 183,
  height = 150,
  units = "mm",
  dpi = 300
)

# Figure S4
figs4.plot <- ggarrange(
  true.homo.ref.as.hetero.mu.in.called.pos.plot,
  true.homo.mu.as.hetero.mu.in.called.pos.plot,
  ncol = 2L,
  align = "h",
  labels = c("a", "b"),
  font.label = list(face = "bold", size = 7),
  common.legend = TRUE,
  legend = "bottom",
  legend.grob = get_legend(true.homo.ref.as.hetero.mu.in.called.pos.plot)
)

ggsave(
  paste0(output.prefix, "s4.pdf"),
  plot = figs4.plot,
  device = "pdf",
  width = 183,
  height = 110,
  units = "mm",
  dpi = 300
)

# Figure S5: ADO calling
figs5.plot <- ggarrange(
  ggarrange(
    ado.recall.plot,
    ado.precision.plot,
    nrow = 2,
    align = "v",
    labels = c("a", "b"),
    font.label = list(face = "bold", size = 7),
    legend = "none"
  ),
  ggarrange(
    ado.fallout.plot,
    ado.f1.score.plot,
    nrow = 2,
    align = "v",
    labels = c("c", "d"),
    font.label = list(face = "bold", size = 7),
    legend = "none"
  ),
  ncol = 2,
  align = "h",
  common.legend = TRUE,
  legend = "bottom",
  legend.grob = get_legend(ado.recall.plot)
)

ggsave(
  paste0(output.prefix, "s5.pdf"),
  plot = figs5.plot,
  device = "pdf",
  width = 183,
  height = 150,
  units = "mm",
  dpi = 300
)
