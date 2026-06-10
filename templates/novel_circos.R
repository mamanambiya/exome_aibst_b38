#!/usr/bin/env Rscript

# Title     : TODO
# Objective : TODO
# Created by: mamana
# Created on: 2019-06-21

library(circlize)
#par(bg = "black") => set image background color is black
#plot circos and create a network:
mat <- read.table("${circos_csv}", header=TRUE)

pdf("${circos_plot}")
chordDiagram(
    mat, 
    annotationTrack = c("grid", "axis"), 
    directional = 1, 
    transparency = 0, 
    preAllocateTracks = list( track.height = uh(4, "mm"), track.margin = c(uh(4, "mm"), 0) )
)
#add a text in each sector 
circos.track(
    track.index = 2, 
    panel.fun = function(x, y) { 
        sector.index = get.cell.meta.data("sector.index") 
        xlim = get.cell.meta.data("xlim") 
        ylim = get.cell.meta.data("ylim") 
        circos.text(
            mean(xlim), 
            mean(ylim), 
            sector.index, 
            col = "white", 
            cex = 0.6, 
            niceFacing = TRUE
            )
    }
)
#add a distic sector in each sector
# highlight.sector(
#     rownames(mat), 
#     track.index = 1, 
#     col = "red", 
#     text = "Chitwan", 
#     cex = 0.8, 
#     text.col = "white", 
#     niceFacing = TRUE
# )
# highlight.sector(
#     colnames(mat), 
#     track.index = 1, 
#     col = "green", 
#     text = "Bardia", 
#     cex = 0.8, 
#     text.col = "white", 
#     niceFacing = TRUE
# )
# highlight.sector(
#     colnames(mat2), 
#     track.index = 1, 
#     col = "blue", 
#     text = "Khaptad", 
#     cex = 0.8, 
#     text.col = "white", 
#     niceFacing = TRUE
# )
# highlight.sector(
#     colnames(mat4), 
#     track.index = 1, 
#     col = "black", 
#     text = "Suklaphata", 
#     cex = 0.8, 
#     text.col = "white", 
#     niceFacing = TRUE
# )
circos.clear()
dev.off()