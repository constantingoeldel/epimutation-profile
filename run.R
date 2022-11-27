if (!requireNamespace("BiocManager", quietly = TRUE ))
  install.install.packages("BiocManager")
BiocManager::install(version = "3.16")

BiocManager::install("AlphaBeta")
BiocManager::install("data.table")
BiocManager::install("dplyr")
remotes::install_github("extendr/rextendr")
library(rextendr)

source("./alphabeta.R")

window_size = 5
for (i in seq(0, 100, by=window_size)) {
  print(i)
  nodeFile = "nodelist.fn"
  edgeFile = "edgelist.fn"
  setwd(sprintf("/mnt/extStorage/constantin/windows/%s", i))
  name = sprintf("epimutation_rate_estimation_window_%s",i)
  directory = sprintf("/mnt/extStorage/constantin/windows/%s/", i)
  run.alphabeta.new(nodelist=nodeFile,
                    edelist=edgeFile,
                    name=name,
                    input.dir=directory,
                    output.dir=directory)
  print("Finished calculating epimutation rate for window")
}





