# This creates the internal vector permorder, which is used to control 
#   permutation of COA columns to get multiple COAs.  
# For purposes of reproducibility, we don't want it to *have* to be random 
#  everytime.  So we use this one pre-established random order.

# DO NOT RERUN THIS CODE as it will replace the reproducible order!!

# Must be in package directory when this is done.  It creates a file
# ./R/sysdata.rda in the package.
setwd("ShapWork")
permorder <- order ( runif ( 1000 ) )
usethis::use_data ( permorder, internal=TRUE )
