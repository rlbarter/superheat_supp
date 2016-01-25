

library(scatterplot3d)
library(rgl)


load("v1_locations.RData")
load("voxel_clusters.RData")



################### ROTATING GIF ###################
voxel.locs <- data.frame(v1_locations)
voxel.locs <- voxel.locs[1:1294,]

voxel.locs$col <- rep("#f46d43", nrow(voxel.locs))
voxel.locs$col[membership == 2] <- "#999999"
width <- 700
open3d()
par3d(windowRect = 50 + c( 0, 0, width, width ) )

rgl.bg(color = "lightgray") # Setup the background color

rgl.spheres(voxel.locs$X1, voxel.locs$X2, voxel.locs$X3,
            color=voxel.locs$col, radius=0.3)

rgl.bbox(color=c("#333377","black"), emission="#333377",
         specular="#3333FF", shininess=5, alpha=0.8, 
         xlen = 0, ylen = 0, zlen = 0) 

movie3d(spin3d(axis = c(0, 0, 1)), duration = 12,
        dir = getwd(), movie="voxel_loc")


# tips: http://www.sthda.com/english/wiki/a-complete-guide-to-3d-visualization-device-system-in-r-r-software-and-data-visualization
