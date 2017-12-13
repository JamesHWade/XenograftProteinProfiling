# Split data for heatmaps and correlation plots
DataSplitHeat <- function(){
        setwd("D:/Box Sync/XPP_Data/")
        library(tidyverse)
        library(reshape2)
        library(gplots) # heatmap.2
        library(Hmisc) # for rcorr function
        library(corrplot)
        
        dat <- read_csv("compiledNormalized.csv")
        
        dat.1 <- dat
        dat.6 <- filter(dat.1, CellLine == "GBM6")
        dat.26 <- filter(dat.1, CellLine == "GBM26")
        dat.1h <- filter(dat.1, TimePoint == "1h")
        dat.24h <- filter(dat.1, TimePoint == "24h")
        dat.6.1 <- filter(dat.1, TimePoint == "1h",
                          CellLine == "GBM6")
        dat.6.24 <- filter(dat.1, TimePoint == "24h",
                           CellLine == "GBM6")
        dat.26.1 <- filter(dat.1, TimePoint == "1h",
                           CellLine == "GBM26")
        dat.26.24 <- filter(dat.1, TimePoint == "24h",
                            CellLine == "GBM26")
        
        casting.1a <- dcast(dat.1, TimePoint + CellLine + Treatment + Target ~
                                    Replicate + n, 
                           value.var = "NormLog", fun.aggregate = mean)

        casting.1 <- dcast(dat.1, TimePoint + CellLine + Treatment ~ Target, 
                           value.var = "NormLog", fun.aggregate = mean)
        casting.6 <- dcast(dat.6, TimePoint + CellLine + Treatment ~ Target,
                           value.var = "NormLog", fun.aggregate = mean)
        casting.26 <- dcast(dat.26, TimePoint + CellLine + Treatment ~Target,
                            value.var = "NormLog", fun.aggregate = mean)
        casting.1h <- dcast(dat.1h, TimePoint + CellLine + Treatment ~ Target,
                            value.var = "NormLog", fun.aggregate = mean)
        casting.24h <- dcast(dat.24h, TimePoint + CellLine + Treatment ~ Target,
                             value.var = "NormLog", fun.aggregate = mean)
        casting.6.1 <- dcast(dat.6.1, TimePoint + CellLine + Treatment ~ Target,
                             value.var = "NormLog", fun.aggregate = mean)
        casting.6.24 <- dcast(dat.6.24, TimePoint + CellLine + Treatment ~
                                      Target,
                              value.var = "NormLog", fun.aggregate = mean)
        casting.26.1 <- dcast(dat.26.1, TimePoint + CellLine + Treatment ~
                                      Target,
                              value.var = "NormLog", fun.aggregate = mean)
        casting.26.24 <- dcast(dat.26.24, TimePoint + CellLine + Treatment ~ 
                                       Target,
                               value.var = "NormLog", fun.aggregate = mean)
        
        rownames(casting.1a) = paste(casting.1a$Treatment,
                                    casting.1a$CellLine,
                                    casting.1a$TimePoint,
                                    casting.1a$Target)
        
        rownames(casting.1) = paste(casting.1$Treatment,
                                    casting.1$CellLine,
                                    casting.1$TimePoint,
                                    sep = " ")
        rownames(casting.6) = paste(casting.6$Treatment,
                                    casting.6$CellLine,
                                    casting.6$TimePoint,
                                    sep = " ")
        rownames(casting.26) = paste(casting.26$Treatment,
                                     casting.26$CellLine, 
                                     casting.26$TimePoint,
                                     sep = " ")
        rownames(casting.1h) = paste(casting.1h$Treatment,
                                     casting.1h$CellLine,
                                     casting.1h$TimePoint,
                                     sep = " ")
        rownames(casting.24h) = paste(casting.24h$Treatment,
                                      casting.24h$CellLine,
                                      casting.24h$TimePoint,
                                      sep = " ")
        rownames(casting.6.1) = paste(casting.6.1$Treatment,
                                      casting.6.1$CellLine, 
                                      casting.6.1$TimePoint,
                                      sep = " ")
        rownames(casting.6.24) = paste(casting.6.24$Treatment,
                                       casting.6.24$CellLine,
                                       casting.6.24$TimePoint,
                                       sep = " ")
        rownames(casting.26.1) = paste(casting.26.1$Treatment,
                                       casting.26.1$CellLine,
                                       casting.26.1$TimePoint,
                                       sep = " ")
        rownames(casting.26.24) = paste(casting.26.24$Treatment,
                                        casting.26.24$CellLine,
                                        casting.6.24$TimePoint,
                                        sep = " ")

        casting.1a[,c(1,2,3,4)] <- NULL
        casting.1[,c(1,2,3)] <- NULL
        casting.6[,c(1,2,3)] <- NULL
        casting.26[,c(1,2,3)] <- NULL
        casting.1h[,c(1,2,3)] <- NULL
        casting.24h[,c(1,2,3)] <- NULL
        casting.6.1[,c(1,2,3)] <- NULL
        casting.6.24[,c(1,2,3)] <- NULL
        casting.26.1[,c(1,2,3)] <- NULL
        casting.26.24[,c(1,2,3)] <- NULL
        
        casting.1a.m <<-data.matrix(casting.1a)
        casting.1.m <<- data.matrix(casting.1)
        casting.6.m <<- data.matrix(casting.6)
        casting.26.m <<- data.matrix(casting.26)
        casting.1h.m <<- data.matrix(casting.1h)
        casting.24h.m <<- data.matrix(casting.24h)
        casting.6.1.m <<- data.matrix(casting.6.1)
        casting.6.24.m <<- data.matrix(casting.6.24)
        casting.26.1.m <<- data.matrix(casting.26.1)
        casting.26.24.m <<- data.matrix(casting.26.24)
        
        allThings <<- rcorr(t(casting.1a))
        both.both.tar <<- rcorr(casting.1.m)
        both.both.txt <<- rcorr(t(casting.1.m))
        GBM6.both.tar <<- rcorr(casting.6.m)
        GBM6.both.txt <<- rcorr(t(casting.6.m))
        GBM26.both.tar <<- rcorr(casting.26.m)
        GBM26.both.txt <<- rcorr(t(casting.26.m))
        both.1.tar <<- rcorr(casting.1h.m)
        both.1.txt <<- rcorr(t(casting.1h.m))
        both.24.tar <<- rcorr(casting.24h.m)
        both.24.txt <<- rcorr(t(casting.24h.m))
        GBM6.1.tar <<- rcorr(casting.6.1.m)
        GBM6.1.txt <<- rcorr(t(casting.6.1.m))
        GBM6.24.tar <<- rcorr(casting.6.24.m)
        GBM6.24.txt <<- rcorr(t(casting.6.24.m))
        GBM26.1.tar <<- rcorr(casting.26.1.m)
        GBM26.1.txt <<- rcorr(t(casting.26.1.m))
        GBM26.24.tar <<- rcorr(casting.26.24.m)
        GBM26.24.txt <<- rcorr(t(casting.26.24.m))
}

# Heatmaps
HeatPlots <- function(){
        # hmcol <- colorRampPalette(c("darkblue", "white", "darkred"))
        setwd("D:/Box Sync/XPP_Plots/Heatmaps")
        hmcol <- colorRampPalette(c("#3a5387", "#e8e8e8", "#b25752"))(256)
        
        png('Heatmap_All.png', width = 8, height = 8,
            units = "in", res = 200)
        heatmap(casting.1.m, scale = "none", col = hmcol, margins = c(10, 10),
                  main = "Everything: GBM 6 & 26 at 1 h & 24 h")
        dev.off()
        
        png('Heatmap_GBM6.png', width = 8, height = 8,
            units = "in", res = 200)
        heatmap(casting.6.m, scale = "none", col = hmcol, margins = c(10, 10),
                  main = "GBM 6 with 1h & 24h")
        dev.off()
        
        png('Heatmap_GBM26.png', width = 8, height = 8,
            units = "in", res = 200)
        heatmap(casting.26.m, scale = "none", col = hmcol, margins = c(10, 10),
                  main = "GBM 26 with 1h & 24h")
        dev.off()
        
        png('Heatmap_1h_both.png', width = 8, height = 8,
            units = "in", res = 200)
        heatmap(casting.1h.m, scale = "none", col = hmcol, margins = c(10, 10),
                  main = "GBM 6 & 26 with 1h")
        dev.off()
        
        png('Heatmap_24h_both.png', width = 8, height = 8,
            units = "in", res = 200)
        heatmap(casting.24h.m, scale = "none", col = hmcol, margins = c(10, 10),
                  main = "GBM 6 & 26 with 24h")
        dev.off()
        
        png('Heatmap_GBM6_1h.png', width = 8, height = 8,
            units = "in", res = 200)
        heatmap(casting.6.1.m, scale = "none", col = hmcol, margins = c(10, 10),
                  main = "GBM 6 with 1h")
        dev.off()
        
        png('Heatmap_GBM6_24h.png', width = 8, height = 8,
            units = "in", res = 200)
        heatmap(casting.6.24.m, scale = "none", col = hmcol, margins = c(10, 10),
                  main = "GBM 6 with 24h")
        dev.off()
        
        png('Heatmap_GBM26_1h.png', width = 8, height = 8,
            units = "in", res = 200)
        heatmap(casting.26.1.m, scale = "none", col = hmcol, margins = c(10, 10),
                  main = "GBM 26 with 1h")
        dev.off()
        
        png('Heatmap_GBM26_24h.png', width = 8, height = 8,
            units = "in", res = 200)
        heatmap(casting.26.24.m, scale = "none", col = hmcol, margins = c(10, 10),
                  main = "GBM 26 with 24h")
        dev.off()
}

# Corrplots
CorrPlots <- function(){
        setwd("D:/Box Sync/XPP_Plots/CorrPlots")
        
        png('Corrplot_All_The_Things.png', width = 8, height = 8,
            units = "in", res = 400)
        corrplot(allThings$r, method = "color", tl.pos = "n")
        dev.off()
        
        png('Corrplot_All_Targets.png', width = 8, height = 8,
            units = "in", res = 200)
        corrplot(both.both.tar$r,
                 title = "Cell Lines: GBM 6, 26 Time Points: 1 h, 24 h",
                 method = "color",
                 # diag = FALSE,
                 # type = "lower",
                 tl.col = "black", tl.cex = 1, cl.cex = 1,
                 mar = c(0, 0, 2, 0))
        dev.off()
        
        png('Corrplot_All_Treatments.png', width = 8, height = 8,
            units = "in", res = 200)
        corrplot(both.both.txt$r,
                 title = "Cell Lines: GBM 6, 26 Time Points: 1 h, 24 h",
                 method = "color",
                 # diag = FALSE,
                 # type = "lower",
                 tl.col = "black", tl.cex = 1, cl.cex = 1,
                 mar = c(0, 0, 2, 0))
        dev.off()
        
        # Time Points
        
        png('Corrplot_1h_Targets.png', width = 8, height = 8,
            units = "in", res = 200)
        corrplot(both.1.tar$r,
                 title = "Cell Lines: GBM 6, 26 Time Points: 1 h",
                 method = "color",
                 # diag = FALSE,
                 # type = "lower",
                 tl.col = "black", tl.cex = 1.25, cl.cex = 1,
                 mar = c(0, 0, 2, 0))
        dev.off()
        
        png('Corrplot_1h_Treatments.png', width = 8, height = 8,
            units = "in", res = 200)
        corrplot(both.1.txt$r,
                 title = "Cell Lines: GBM 6, 26 Time Points: 1 h",
                 method = "color",
                 # diag = FALSE,
                 # type = "lower",
                 tl.col = "black", tl.cex = 1.25, cl.cex = 1,
                 mar = c(0, 0, 2, 0))
        dev.off()
        
        png('Corrplot_24h_Targets.png', width = 8, height = 8,
            units = "in", res = 200)
        corrplot(both.24.tar$r,
                 title = "Cell Lines: GBM 6, 26 Time Points: 24 h",
                 method = "color",
                 # diag = FALSE,
                 # type = "lower",
                 tl.col = "black", tl.cex = 1.25, cl.cex = 1,
                 mar = c(0, 0, 2, 0))
        dev.off()
        
        png('Corrplot_24h_Treatments.png', width = 8, height = 8,
            units = "in", res = 200)
        corrplot(both.24.txt$r,
                 title = "Cell Lines: GBM 6, 26 Time Points: 24 h",
                 method = "color",
                 # diag = FALSE,
                 # type = "lower",
                 tl.col = "black", tl.cex = 1.25, cl.cex = 1,
                 mar = c(0, 0, 2, 0))
        dev.off()
        
        # GBM 6
        
        png('Corrplot_GBM6_Targets.png', width = 8, height = 8,
            units = "in", res = 200)
        corrplot(GBM6.both.tar$r,
                 title = "Cell Lines: GBM 6 Time Points: 1 h, 24 h",
                 method = "color",
                 # diag = FALSE,
                 # type = "lower",
                 tl.col = "black", tl.cex = 1.25, cl.cex = 1,
                 mar = c(0, 0, 2, 0))
        dev.off()
        
        
        png('Corrplot_GBM6_Treatments.png', width = 8, height = 8,
            units = "in", res = 200)
        corrplot(GBM6.both.txt$r,
                 title = "Cell Lines: GBM 6 Time Points: 1 h, 24 h",
                 method = "color",
                 # diag = FALSE,
                 # type = "lower",
                 tl.col = "black", tl.cex = 1.25, cl.cex = 1,
                 mar = c(0, 0, 2, 0))
        dev.off()
        
        png('Corrplot_GBM6_1h_Treatments.png', width = 8, height = 8,
            units = "in", res = 200)
        corrplot(GBM6.1.txt$r,
                 title = "Cell Lines: GBM 6 Time Points: 1 h",
                 method = "color",
                 # diag = FALSE,
                 # type = "lower",
                 tl.col = "black", tl.cex = 1.25, cl.cex = 1,
                 mar = c(0, 0, 2, 0))
        dev.off()
        
        png('Corrplot_GBM6_1h_Targets.png', width = 8, height = 8,
            units = "in", res = 200)
        corrplot(GBM6.1.tar$r,
                 title = "Cell Lines: GBM 6 Time Points: 1 h",
                 method = "color",
                 # diag = FALSE,
                 # type = "lower",
                 tl.col = "black", tl.cex = 1.25, cl.cex = 1,
                 mar = c(0, 0, 2, 0))
        dev.off()
        
        png('Corrplot_GBM6_24h_Targets.png', width = 8, height = 8,
            units = "in", res = 200)
        corrplot(GBM6.24.tar$r,
                 title = "Cell Lines: GBM 6 Time Points: 24 h",
                 method = "color",
                 # diag = FALSE,
                 # type = "lower",
                 tl.col = "black", tl.cex = 1.25, cl.cex = 1,
                 mar = c(0, 0, 2, 0))
        dev.off()
        
        png('Corrplot_GBM6_24h_Treatments.png', width = 8, height = 8,
            units = "in", res = 200)
        corrplot(GBM6.24.txt$r,
                 title = "Cell Lines: GBM 6 Time Points: 24 h",
                 method = "color",
                 # diag = FALSE,
                 # type = "lower",
                 tl.col = "black", tl.cex = 1.25, cl.cex = 1,
                 mar = c(0, 0, 2, 0))
        dev.off()
        
        # GBM 26
        
        png('Corrplot_GBM26_Targets.png', width = 8, height = 8,
            units = "in", res = 200)
        corrplot(GBM26.both.tar$r,
                 title = "Cell Lines: GBM 26 Time Points: 1 h, 24 h",
                 method = "color",
                 # diag = FALSE,
                 # type = "lower",
                 tl.col = "black", tl.cex = 1.25, cl.cex = 1,
                 mar = c(0, 0, 2, 0))
        dev.off()
        
        
        png('Corrplot_GBM26_Treatments.png', width = 8, height = 8,
            units = "in", res = 200)
        corrplot(GBM26.both.txt$r,
                 title = "Cell Lines: GBM 26 Time Points: 1 h, 24 h",
                 method = "color",
                 # diag = FALSE,
                 # type = "lower",
                 tl.col = "black", tl.cex = 1.25, cl.cex = 1,
                 mar = c(0, 0, 2, 0))
        dev.off()
        
        png('Corrplot_GBM26_1h_Treatments.png', width = 8, height = 8,
            units = "in", res = 200)
        corrplot(GBM26.1.txt$r,
                 title = "Cell Lines: GBM 26 Time Points: 1 h",
                 method = "color",
                 # diag = FALSE,
                 # type = "lower",
                 tl.col = "black", tl.cex = 1.25, cl.cex = 1,
                 mar = c(0, 0, 2, 0))
        dev.off()
        
        png('Corrplot_GBM26_1h_Targets.png', width = 8, height = 8,
            units = "in", res = 200)
        corrplot(GBM26.1.tar$r,
                 title = "Cell Lines: GBM 26 Time Points: 1 h",
                 method = "color",
                 # diag = FALSE,
                 # type = "lower",
                 tl.col = "black", tl.cex = 1.25, cl.cex = 1,
                 mar = c(0, 0, 2, 0))
        dev.off()
        
        png('Corrplot_GBM26_24h_Treatments.png', width = 8, height = 8,
            units = "in", res = 200)
        corrplot(GBM26.24.txt$r,
                 title = "Cell Lines: GBM 26 Time Points: 24 h",
                 method = "color",
                 # diag = FALSE,
                 # type = "lower",
                 tl.col = "black", tl.cex = 1.25, cl.cex = 1,
                 mar = c(0, 0, 2, 0))
        dev.off()
        
        png('Corrplot_GBM26_24h_Targets.png', width = 8, height = 8,
            units = "in", res = 200)
        corrplot(GBM26.24.tar$r,
                 title = "Cell Lines: GBM 26 Time Points: 24 h",
                 method = "color",
                 # diag = FALSE,
                 # type = "lower",
                 tl.col = "black", tl.cex = 1.25, cl.cex = 1,
                 mar = c(0, 0, 2, 0))
        dev.off()
}

GoGoHeatmaps <- function(i){
        DataSplitHeat()
        HeatPlots()
        CorrPlots()
}
