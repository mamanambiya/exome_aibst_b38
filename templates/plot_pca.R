#!/usr/bin/env Rscript
# Title     : TODO
# Objective : TODO
# Created by: mamana
# Created on: 2019-06-21

library(ggplot2)
library(reshape2)

shapes <- c(0:25)
evecDat = read.table("${input_evec}", col.names=c("Sample", "PC1", "PC2", "PC3", "PC4", "PC5", "Population", "Group", "Supergroup"))

# if (nlevels(evecDat\$Group) > 1) {
#     output_pdf <- "${output_pdf}_pc1_pc2.tiff"
#     p <- ggplot(data=evecDat, aes(x=PC1, y=PC2, colour=Supergroup, shape=Supergroup)) +
#     geom_point(size=2) +
#     scale_shape_manual(values=rep(shapes, times=4))
#     p + theme( panel.background= element_rect(fill="#E3E3EE"))
#     ggsave(output_pdf, units='in', dpi=120)

#     # scale_shape_identity(breaks = unique(evecDat\$Pop), guide = 'legend')
#     output_pdf <- "${output_pdf}_pc1_pc3.tiff"
#     p <- ggplot(data=evecDat, aes(x=PC1, y=PC3, colour=Supergroup, shape=Supergroup)) +
#         geom_point(size=2) +
#         scale_shape_manual(values=rep(shapes, times=4))
#     p + theme( panel.background= element_rect(fill="#E3E3EE"))
#     ggsave(output_pdf, units='in', dpi=120)
# } else {
    shapes <- c(3, 1, 3, 2, 5, 3, 5, 5, 6, 3, 5, 7, 9,10,11, 3, 3,14,15,16,17,18,19,20, 5,21,22,12, 3,0,4,8)
    output_pdf <- "${output_pdf}_pc1_pc2.tiff"
    p <- ggplot(data=evecDat, aes(x=PC1, y=PC2, colour=Population, shape=Population)) +
    geom_point(size=2) +
    scale_shape_manual(values=rep(shapes, times=4))
    p + theme( panel.background=element_rect(fill="#E3E3EE")) + theme(legend.position="right")
    ggsave(output_pdf, units='in', dpi=120, width=16, height=8.5)

    output_pdf <- "${output_pdf}_pc1_pc3.tiff"
    p <- ggplot(data=evecDat, aes(x=PC1, y=PC3, colour=Population, shape=Population)) +
    geom_point(size=2) +
    scale_shape_manual(values=rep(shapes, times=4))
    p + theme( panel.background=element_rect(fill="#E3E3EE"))
    ggsave(output_pdf, units='in', dpi=120, width=16, height=8.5)
# }
