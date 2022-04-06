reldata <- function(x, byrow = TRUE, bycol = FALSE, rowfirst=TRUE)
{
# Relativize x by row sum, column maximum or both.
# See also scale().
# offers option to process row or column first
# wisconsin standardization: rel(x, byrow = TRUE, bycol=TRUE, rowfirst=FALSE)
# SCG 2012-02-21
#

    if(rowfirst) {
        if(byrow) {
            rsum <- apply(x, 1, sum)
            rsum[rsum == 0] <- 1
            x <- sweep(x, 1, rsum, "/")
        }
        if(bycol) {
            cmax <- apply(x, 2, max)
            cmax[cmax == 0] <- 1
            x <- sweep(x, 2, cmax, "/")
        }
    } else {
        if(bycol) {
            cmax <- apply(x, 2, max)
            cmax[cmax == 0] <- 1
            x <- sweep(x, 2, cmax, "/")
        }
        if(byrow) {
            rsum <- apply(x, 1, sum)
            rsum[rsum == 0] <- 1
            x <- sweep(x, 1, rsum, "/")
        }
    }
    x
}

