library(lattice)
library(grid)

panel.violin.hack <-
function (x, y, box.ratio = 1, box.width = box.ratio/(1 + box.ratio),
     horizontal = TRUE, alpha = plot.polygon$alpha, border =  
plot.polygon$border,
     lty = plot.polygon$lty, lwd = plot.polygon$lwd, col = plot.polygon 
$col,
     varwidth = FALSE, bw = NULL, adjust = NULL, kernel = NULL,
     window = NULL, width = NULL, n = 50, from = NULL, to = NULL,
     cut = NULL, na.rm = TRUE, ...)
{
     if (all(is.na(x) | is.na(y)))
         return()
     x <- as.numeric(x)
     y <- as.numeric(y)
     plot.polygon <- trellis.par.get("plot.polygon")
     darg <- list()
     darg$bw <- bw
     darg$adjust <- adjust
     darg$kernel <- kernel
     darg$window <- window
     darg$width <- width
     darg$n <- n
     darg$from <- from
     darg$to <- to
     darg$cut <- cut
     darg$na.rm <- na.rm
     my.density <- function(x) {
         ans <- try(do.call("density", c(list(x = x), darg)),
             silent = TRUE)
         if (inherits(ans, "try-error"))
             list(x = rep(x[1], 3), y = c(0, 1, 0))
         else ans
     }
     numeric.list <- if (horizontal)
         split(x, factor(y))
     else split(y, factor(x))
     levels.fos <- as.numeric(names(numeric.list))
     d.list <- lapply(numeric.list, my.density)
     dx.list <- lapply(d.list, "[[", "x")
     dy.list <- lapply(d.list, "[[", "y")
     max.d <- sapply(dy.list, max)
     if (varwidth)
         max.d[] <- max(max.d)
     xscale <- current.panel.limits()$xlim
     yscale <- current.panel.limits()$ylim
     height <- box.width
     if (horizontal) {
         for (i in seq_along(levels.fos)) {
             if (is.finite(max.d[i])) {
                 pushViewport(viewport(y = unit(levels.fos[i],
                   "native"), height = unit(height, "native"),
                   yscale = c(max.d[i] * c(-1, 1)), xscale = xscale))
                 grid.polygon(x = c(dx.list[[i]], rev(dx.list[[i]])),
                   y = c(dy.list[[i]], -rev(dy.list[[i]])),  
default.units = "native",
# this is the point at which the index is added
                   gp = gpar(fill = col[i], col = border, lty = lty,
                     lwd = lwd, alpha = alpha))
                 popViewport()
             }
         }
     }
     else {
         for (i in seq_along(levels.fos)) {
             if (is.finite(max.d[i])) {
                 pushViewport(viewport(x = unit(levels.fos[i],
                   "native"), width = unit(height, "native"),
                   xscale = c(max.d[i] * c(-1, 1)), yscale = yscale))
                 grid.polygon(y = c(dx.list[[i]], rev(dx.list[[i]])),
                   x = c(dy.list[[i]], -rev(dy.list[[i]])),  
default.units = "native",
# this is the point at which the index is added
                   gp = gpar(fill = col[i], col = border, lty = lty,
                     lwd = lwd, alpha = alpha))
                 popViewport()
             }
         }
     }
     invisible()
}


load('~/Desktop/ZNF143_ChIP/actDownMatchedRepDown.Rdata')
load('~/Desktop/ZNF143_ChIP/actOnTssMatchedRep.Rdata')

actDownMatchedRepDown$status <- ifelse(actDownMatchedRepDown$status == "Activated", paste0(actDownMatchedRepDown$status, " (28)"), actDownMatchedRepDown$status)
actDownMatchedRepDown$status <- ifelse(actDownMatchedRepDown$status == "Matched to Activated", paste0(actDownMatchedRepDown$status, " (31)"), actDownMatchedRepDown$status)
actDownMatchedRepDown$status <- ifelse(actDownMatchedRepDown$status == "Repressed", paste0(actDownMatchedRepDown$status, " (30)"), actDownMatchedRepDown$status)

actOnTssMatchedRep$status <- ifelse(actOnTssMatchedRep$status == "Activated", paste0(actOnTssMatchedRep$status, " (5)"), actOnTssMatchedRep$status)
actOnTssMatchedRep$status <- ifelse(actOnTssMatchedRep$status == "Matched to Activated", paste0(actOnTssMatchedRep$status, " (27)"), actOnTssMatchedRep$status)
actOnTssMatchedRep$status <- ifelse(actOnTssMatchedRep$status == "Repressed", paste0(actOnTssMatchedRep$status, " (33)"), actOnTssMatchedRep$status)


pdf(file = "Violin_BWplot_downsteam.pdf" ,width=5,height=3.5)

trellis.par.set(box.umbrella = list(lty = 1, col="black", lwd=1),
                box.rectangle = list(col = 'black', lwd=1), plot.symbol = list(col='blue', lwd=1, pch ='.'))
print(
    bwplot(peakIntens ~ status,
         data=actDownMatchedRepDown, #subset=(abs(Lat)<60),
           xlab='Gene Response',
           ylab='ChIP Intensity',
         #main = 'Changes in Rates',
         ylim = c(0,10500),
           aspect = 1.0, 
           horizontal=FALSE,
           between=list(y=0.7, x=0.7),
           scales=list(x=list(alternating=c(1,1,1,0,0,0),rot=30, font = 1, cex=0.8),
               y=list(alternating=c(1,1))),
           panel = function(..., box.ratio) {
               
                panel.violin.hack(..., col = c("grey95", "grey95", "grey95"),
                            varwidth = FALSE, box.ratio = 10)
               panel.abline(h=1, lty =2)
               panel.bwplot(..., col='red',
                            cex=0.8, pch='|', fill='transparent', box.ratio = .2, do.out = FALSE)
	       panel.stripplot(..., col='#b3b3b3CC', do.out=FALSE, jitter.data=TRUE, amount = 0.2, pch = 16, cex = 0.3)
           },
           par.settings = list(box.rectangle=list(col='black'),
               strip.background=list(col="grey90"),
               plot.symbol = list(pch='.', cex = 0.1))
           )
           )
dev.off()


pdf(file = "Violin_BWplot_TSS.pdf" ,width=5,height=3.5)

trellis.par.set(box.umbrella = list(lty = 1, col="black", lwd=1),
                box.rectangle = list(col = 'black', lwd=1), plot.symbol = list(col='blue', lwd=1, pch ='.'))
print(
    bwplot(peakIntens ~ status,
         data=actOnTssMatchedRep, #subset=(abs(Lat)<60),
           xlab='Gene Response',
           ylab='ChIP Intensity',
         #main = 'Changes in Rates',
         ylim = c(0,10500),
           aspect = 1.0, 
           horizontal=FALSE,
           between=list(y=0.7, x=0.7),
           scales=list(x=list(alternating=c(1,1,1,0,0,0),rot=30, font = 1, cex=0.8),
               y=list(alternating=c(1,1))),
           panel = function(..., box.ratio) {
               
                panel.violin.hack(..., col = c("grey95", "grey95", "grey95"),
                            varwidth = FALSE, box.ratio = 10)
               panel.abline(h=1, lty =2)
               panel.bwplot(..., col='red',
                            cex=0.8, pch='|', fill='transparent', box.ratio = .2, do.out = FALSE)
	       panel.stripplot(..., col='#b3b3b3CC', do.out=FALSE, jitter.data=TRUE, amount = 0.2, pch = 16, cex = 0.3)
           },
           par.settings = list(box.rectangle=list(col='black'),
               strip.background=list(col="grey90"),
               plot.symbol = list(pch='.', cex = 0.1))
           )
           )
dev.off()

