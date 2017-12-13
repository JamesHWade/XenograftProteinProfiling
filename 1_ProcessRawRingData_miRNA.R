

PlotCombineNetShifts <- function(loc = 'plots', ch){
        g1 <- ggplot(dat, aes(x = Cycle, y = NetShift, 
                              color = Target, group = Cycle)) + 
                geom_boxplot() + facet_wrap(~Target)
        
        g2 <- ggplot(dat, aes(x = Cycle, y = NetShift, 
                              color = factor(Ring), group = Ring)) +
                geom_point() +
                ggtitle(name) + facet_wrap(~Target)
        
        ggsave(plot = g2, filename = paste0(loc, "/", name,
                                            "_netShiftsCombined", 
                                            "_boxplot.png"),
               width = 8, height = 6)
        
        ggsave(plot = g2, filename = paste0(loc, "/", name,
                                            "_netShiftsCombined",  
                                            ".png"),
               width = 8, height = 6)
}

Fit <- function(loc = 'plots', ch){
        dat <- read_csv(paste0(paste0(loc, "/", name, "_netShiftsCombined_ch", 
                                      ch, ".csv")))
        
        dat$Replicate <- as.factor((dat$Ring - 1) %% 4 + 1)
        
        dat.fit <- select(dat, c(Target, Experiment, Cycle, NetShift,
                                 Replicate, Ring))
        
        dat.fit <- filter(dat.fit, Target != "Off target")
        
        fit <- list()
        targetList <- unique(dat.fit$Target)
        
        ggplot(dat, aes(x = Cycle, y = NetShift, color = Target)) + 
                geom_point() +
                geom_smooth(formula = y ~ x / (1 + x)) +
                facet_wrap(~Target)
        
        for(i in 1:length(targetList)) {
                tar <- targetList[i]
                print(tar)
                dat.tar <- filter(dat.fit, Target == tar)
                y <- dat.tar$newShift
                x <- dat.tar$Step
                #fit[[i]] <- nls(y ~ SSlogis(x, Asym, xmid, scal))
                fit[[i]] <- nls(formula = y ~ A.2 + (A.1-A.2)/(1 + (x/x.0)^p),
                                start = list(A.2 =max(y),
                                             A.1 = min(y),
                                             x.0 = mean(y),
                                             p = 1))
        }
        
        capture.output(fit, file = paste0(loc, "_dataFit.txt"))
        
}

miRNAData <- function(loc, cntl, filename = 'groupNames_miRNA.csv') {
        GetName()
        AggData(loc = loc, filename = filename)
        SubtractControl(ch = ch, cntl = cntl)
        SubtractControl(ch = 1, cntl = cntl)
        SubtractControl(ch = 2, cntl = cntl)
        SubtractControl(ch = "U", cntl = cntl)
        PlotRingData(cntl = "raw", ch = "U", splitPlot = TRUE)
        PlotRingData(cntl = cntl, ch = "U", splitPlot = TRUE)
        PlotRingData(cntl = cntl, ch = 1, splitPlot = FALSE)
        PlotRingData(cntl = "raw", ch = 1, splitPlot = FALSE)
        PlotRingData(cntl = cntl, ch = 2, splitPlot = FALSE)
        PlotRingData(cntl = "raw", ch = 2, splitPlot = FALSE)
        GetNetShifts(cntl = cntl, ch = 1, 
                     time1 = 22.5, time2 = 5, step = 20)
        GetNetShifts(cntl = cntl, ch = 1, 
                     time1 = 40, time2 = 22.5, step = 25)
        GetNetShifts(cntl = cntl, ch = 1, 
                     time1 = 57, time2 = 40, step = 30)
        GetNetShifts(cntl = cntl, ch = 1, 
                     time1 = 74, time2 = 57, step = 35)
        GetNetShifts(cntl = cntl, ch = 1, 
                     time1 = 92, time2 = 74, step = 40)
        GetNetShifts(cntl = cntl, ch = 1, 
                     time1 = 110, time2 = 92, step = 45)
        GetNetShifts(cntl = cntl, ch = 2, 
                     time1 = 22.5, time2 = 5, step = 20)
        GetNetShifts(cntl = cntl, ch = 2, 
                     time1 = 40, time2 = 22.5, step = 25)
        GetNetShifts(cntl = cntl, ch = 2, 
                     time1 = 57, time2 = 40, step = 30)
        GetNetShifts(cntl = cntl, ch = 2, 
                     time1 = 74, time2 = 57, step = 35)
        GetNetShifts(cntl = cntl, ch = 2, 
                     time1 = 92, time2 = 74, step = 40)
        GetNetShifts(cntl = cntl, ch = 2, 
                     time1 = 110, time2 = 92, step = 45)
        CombineNetShifts(ch = "U")
        PlotCombineNetShifts(ch = "U")
        Fit(loc = loc)
}

AnalyzeAllData <- function() {
        foldersList <- list.dirs(recursive = FALSE)
        for (i in foldersList){
                directory <- getwd()
                setwd(i)
                AnalyzeData()
                setwd(directory)
        }
}

AnalyzeData(loc = "plots", cntl = "Control", filename = "groupNames_miRNA.csv",
            ch = 1)
AnalyzeData(loc = "plots", cntl = "control", filename = "groupNames_miRNA.csv",
            ch = 2)
