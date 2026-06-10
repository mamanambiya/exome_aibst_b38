#!/users/mamana/miniconda3/envs/r-upsetr/bin/Rscript
# Title     : TODO
# Objective : TODO
# Created by: mamana
# Created on: 2019-06-21

# library(rJava)
library(UpSetR)
# library(tidyverse)
# library(venneuler)
library(grid)
    
rawSets <- read.csv(
          file = "${csv_file}",
          header = TRUE, sep = ",", stringsAsFactors = FALSE
        )
# Replace the NA's
rawSets[is.na(rawSets)] <- 0

pdf("${upset_plot}")

upset(rawSets,
  nsets = 12, number.angles = 30, point.size = 3.5, line.size = 2
)

# upset(rawSets,
#   nsets = 10, number.angles = 30, point.size = 3.5, line.size = 2,
#   mainbar.y.label = "Private variants", sets.x.label = "Dataset"
# )

# grid.text(
#   "@littlemissdata",
#   x = 0.90,
#   y = 0.05,
#   gp = gpar(
#     fontsize = 10,
#     fontface = 3
#   )
# )

dev.off()
