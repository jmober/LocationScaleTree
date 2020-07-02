##############################################################################
#####   TREE-STRUCTURED SCALE EFFECTS IN BINARY AND ORDINAL REGRESSION   #####
#####                     Electronic Supplement                          #####
##############################################################################
#####                    Author: Moritz Berger                           #####
##############################################################################

# Exemplary analysis fitting a tree-structured location-scale model 
# (2.3 Illustrative Example)



# load all necessary files 
load("./confidence.rda") # the original data is available from
                         # http://www.gesis.org/allbus
source("./functions.R")
source("./LocScaleTree.R")
source("./plot.LocScaleTree.R")

# response vector and design matrix 
y      <- data[,1]
DM_kov <- data[,-1]


# fit tree (note: might take some time, depending on the number of permutations)
mod <- LocScaleTree(y,
                    DM_kov,
                    alpha=0.05,
                    nperm=50,
                    trace=TRUE)

# load fit (model based on 1000 permutations (nperm=1000))
load("./mod_confidence.rda")

# fitted threshold coefficients
mod$theta_hat

# location term 
par(mar=c(0,0,2,0))
plot(mod,"beta", ellipse_a=2.4, ellipse_b=0.25, ellipse_x=1.5, ellipse_y=0, branch_adj=-0.4, title="location term", draw_numbers=T, 
     cex.branches=1.8, cex.coefs=1.8, cex.main=2.5, cex.number=1.3)

# variance term 
plot(mod,"gamma", ellipse_a=0.5, ellipse_b=0.2, ellipse_x=0.1, ellipse_y=0, branch_adj=-0.4, title="variance term", draw_numbers=T, 
     cex.branches=1.8, cex.coefs=1.8, cex.main=2.5, cex.number=1.3)

