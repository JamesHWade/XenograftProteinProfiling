GetName <- function(){
        # get the filename from the current working directory
        directory <- basename(getwd())
        
        # directory naming from MRR: "CHIPNAME_gaskGASK_DATE"
        # extracts and returns GASK from directory name
        name <- unlist(strsplit(directory, split = "_"))[2]
        
        # define name as global variable for use in other functions
        name <<- gsub('gask','',name) # removes "gask" from name
}

AggData <- function(loc, filename, fsr, fsrThresh = 3000) {
        # load relevant libraries
        library(tidyverse)
        
        # get working directory to reset at end of function
        directory <- getwd()
        
        # get information of chip layout from github repository
        if (!file.exists(filename)){
                git <- "https://raw.githubusercontent.com/"
                hub <- "JamesHWade/XenograftProteinProfiling/master/"
                github <- paste0(git, hub)
                url <- paste0(github, filename)
                filename <- basename(url)
                download.file(url, filename)
        }
        
        # read in recipe/chip layout
        recipe <- read_csv(filename)
        colnames(recipe)[1] <- "Target" # rename col & remove byte order mark
        targets <- recipe$Target
        
        # generate list of rings to analyze (gets all *.csv files)
        rings <- list.files(directory, pattern = "[[:digit:]].csv", recursive = FALSE)
        
        # add data to data frame corresponding for each ring in rings
        df <- lapply(rings, function(i){
                dat <- read_csv(i, col_names = FALSE)
                ringNum <- as.numeric(strsplit(i, "\\.")[[1]][1])
                recipeCol <- which(recipe$Ring == ringNum)
                tmp <- dat[,c(1,2)] # time and shift from raw data
                tmp$ring <- ringNum
                tmp$group <- recipe$Group[recipeCol]
                tmp$groupName <- as.character(recipe$Target[[recipeCol]])
                tmp$channel <- recipe$Channel[[recipeCol]]
                tmp$run <- name
                tmp$timePoint <- seq(1:nrow(dat))
                tmp
        })
        
        # correct for fsr
        if(fsr){
                for(i in seq_len(length(df))){
                        pointShift <- 0
                        for(j in seq_len(nrow(df[[i]]))){
                                shiftDiff <- pointShift - df[[i]][j, 2]
                                if(shiftDiff > fsrThresh){
                                        print("fuck, an fsr---kill it---got it")
                                        df[[i]][j, 2] <- df[[i]][j, 2] + 5980
                                }
                                pointShift <- df[[i]][j, 2]
                        }
                }
        }
        
        # combine data from list into single data frame
        df <- bind_rows(df)
        
        # renames columns in df
        names(df) <- c("Time", "Shift", "Ring", "Group", "Target", "Channel",
                       "Experiment", "Time Point")
        
        # creates "plots" directory
        dir.create(loc, showWarnings = FALSE)
        
        # saves aggregated data with name_allRings.csv
        write_csv(df, paste0(loc, '/', name, "_allRings.csv"))
}

SubtractControl <- function(loc, ch, cntl){
        #load relevant libraries
        library(tidyverse)
        
        # get ring data and filter by channel
        dat <- read_csv(paste0(loc, "/", name, "_", "allRings.csv"))
        if (ch != "U"){dat <- filter(dat, Channel == ch)}
        dat <- filter(dat, Target != "Ignore")
        
        # get thermal control averages
        controls <- filter(dat, Target == cntl) %>% group_by(`Time Point`) %>%
                summarise_at("Shift", mean) %>% select(Shift) %>% unlist()
        dat$Cntl <- rep(controls, length(unique(dat$Ring)))
        
        # subtracts thermal controls from each ring
        dat.cntl <- mutate(dat, Shift = Shift - Cntl)
        
        # remove control column and control rings
        dat.cntl <- filter(dat.cntl, Target != cntl & Target != "thermal")
        dat.cntl$Cntl <- NULL
        
        # save data to new file
        write_csv(dat.cntl, paste(loc,"/", name, "_", cntl, "Control", "_ch",
                                  ch, ".csv", sep = ''))   
}

PlotRingData <- function(cntl, ch, loc, splitPlot = FALSE){
        # loads relevant libraries and plot theme
        library(tidyverse)
        library(ggthemes)
        theme_set(theme_few(base_size = 16))
        
        # use thermally controlled data if desired
        if (cntl != "raw"){
                dat <- read_csv(paste(loc, "/", name, "_", cntl, "Control", 
                                      "_ch", ch,".csv", sep=''))
        } else if (cntl == "raw") {
                dat <- read_csv(paste(loc, "/", name, "_allRings.csv", sep=''))
                if (ch != "U") {dat <- filter(dat, Channel == ch)}
        }
        
        # configure plot and legend
        plots <- ggplot(dat, aes(x = Time, y = Shift,
                                 color = Target, group = Ring)) + 
                labs(x = "Time (min)", 
                     y = expression(paste("Relative Shift (",Delta,"pm)")),
                     color = "Target") +
                geom_line() + 
                ggtitle(paste(name, "Ch:", ch, "Control:", cntl, sep = " "))
        
        # alternative plots with averaged clusters
        
        dat.2 <- dat %>% group_by(Target, `Time Point`) %>% 
                summarise_at(vars(Time, Shift), funs(mean, sd))
        
        plot2 <- ggplot(dat.2, aes(x = Time_mean, y = Shift_mean, 
                                   color = Target)) +
                geom_line() +
                labs(x = "Time (min)", 
                     y = expression(paste("Relative Shift (",Delta,"pm)"))) +
                ggtitle(paste(name, "Ch:", ch, "Control:", cntl, sep = " "))
        
        plot3 <- plot2 + 
                geom_ribbon(aes(ymin = Shift_mean - Shift_sd,
                                ymax = Shift_mean + Shift_sd, 
                                linetype = NA),
                            fill = "slategrey", alpha = 1/8)
        
        if (splitPlot){
                plots <- plots + facet_grid(. ~ Channel)
        }
        
        # save plots
        ggsave(plot = plots, 
               file = paste0(loc, "/", name, "_", cntl, 
                             "Control_ch", ch, ".png"),
               width = 10, height = 6)
        ggsave(plot = plot2, 
               file = paste0(loc, "/", name, "_", cntl,
                             "Control", "_ch", ch, "_avg.png"),
               width = 10, height = 6)
        ggsave(plot = plot3, 
               file = paste0(loc, "/", name, "_", cntl,
                             "Control", "_ch", ch, "_avg.png"),
               width = 10, height = 6)
}

GetNetShifts <- function(cntl, ch, loc, time1, time2, step = 1){
        # load relevant libraries
        library(tidyverse)
        
        # use thermally controlled data if desired
        if (cntl != "raw"){
                dat <- read_csv(paste0(loc, "/", name, "_", cntl, "Control", 
                                       "_ch", ch, ".csv"))
        } else {
                dat <- read_csv(paste0(loc, "/", name, "_", "allRings.csv"))
        }
        
        # generate list of rings and empty dataframe to store net shift data
        ringList <- unique(dat$Ring)
        
        # locations for each time is determined using which, min, and abs func
        dat.rings <- lapply(ringList, function(i){
                dat.ring <- filter(dat, Ring == i)
                time1.loc <- which.min(abs(dat.ring$Time - time1))
                time1.val <- dat.ring$Shift[time1.loc]
                time2.loc <- which.min(abs(dat.ring$Time - time2))
                time2.val <- dat.ring$Shift[time2.loc]
                ring <- i
                group <- unique(dat.ring$Group)
                target <- unique(dat.ring$Target)
                experiment <- unique(dat.ring$Experiment)
                channel <- unique(dat.ring$Channel)
                data.frame(i, group, target, time1.val,
                           time2.val, experiment, channel, step)
        })
        
        # renames dat.rings columns
        dat.rings <- bind_rows(dat.rings)
        names(dat.rings) <- c("Ring", "Group", "Target", "Shift.1", "Shift.2", 
                              "Experiment", "Channel", "Step")
        
        # calculate nat shift and create new column in dataframe
        dat.rings <- dat.rings %>% 
                mutate(NetShift = Shift.1 - Shift.2)
        
        # save net shift data
        write_csv(dat.rings, paste0(loc, "/", name, "_netShifts_", cntl,
                                    "cntl_", "ch", ch, "_step", step, ".csv"))
}

PlotNetShifts <- function(cntl, ch, loc, step = 1){
        # load relevant libraries and plot theme
        library(tidyverse)
        library(ggthemes)
        theme_set(theme_few(base_size = 16))
        
        dat <- read_csv(paste0(loc, "/", name, "_netShifts_", cntl, "cntl_",
                               "ch", ch, "_step", step, ".csv"))
        
        # configure plot and legend
        dat.nothermal <- filter(dat, Target != "thermal")
        
        plots <- ggplot(dat.nothermal, 
                        aes(x = Target, y = NetShift, fill = Target)) +
                geom_boxplot() +
                theme(axis.text.x = element_text(angle = 45, hjust = 1),
                      legend.position="none") +
                ylab(expression(paste("Net Shift (",Delta,"pm)"))) +
                ggtitle(paste0(name, " Ch: ", ch, " Control: ", cntl))
        
        allRings <- ggplot(dat.nothermal, 
                           aes(x = factor(Ring), y= NetShift, 
                               fill = Target)) +
                geom_bar(stat = "identity") +
                theme(axis.text.x = 
                              element_text(angle = 90,
                                           hjust = 1, vjust = 0.5)) +
                ylab(expression(paste("Net Shift (",Delta,"pm)"))) + 
                xlab("Ring") +
                ggtitle(paste0(name, " Ch: ", ch, " Control: ", cntl))
        
        # save plot, uncomment to save
        ggsave(plot = plots, 
               file = paste0(loc, "/",  name, "_NetShift_", cntl, "cntl_",
                             "ch", ch, "_step", step, ".png"), 
               width = 10, height = 6)
        ggsave(plot = allRings, 
               file = paste0(loc, "/", name, "_IndyRings","_NetShift_", cntl, 
                             "cntl_", "ch", ch, "_step", step, ".png"), 
               width = 12, height = 6)
}

CheckRingQuality <- function(loc, time1, time2, nrings = 10) {
        # load relevant libraries
        library(tidyverse)
        library(ggthemes)
        
        # read in data and subset for a flat part of the run
        dat <- read_csv(paste0(loc,"/", name, "_allRings.csv"))
        dat <- subset(dat, Time > time1 & Time < time2)
        
        # fn to take absolute value of max signal
        absmax <- function(x) { x[which.max( abs(x) )]}
        
        # calculate variance and max signal for each ring
        dat.var <- dat %>% group_by(Ring) %>%
                summarise_at(vars(Shift), funs(Variance = var, Max = absmax))
        
        # plot Variance vs Max signal (on log axis)
        g1 <- ggplot(dat.var, aes(x = Variance, y = Max, 
                                  color = factor(Ring))) + 
                geom_point() + 
                scale_x_log10() + scale_y_log10()
        
        # create variables for rings with variance above/below given variance
        ringWinners <- arrange(dat.var, Variance) %>% select(Ring) %>% 
                head(nrings)
        ringLosers <- arrange(dat.var, Variance) %>% select(Ring) %>% 
                tail(nrings)
        
        # save files with list of good and bad rings base on given variance
        write_csv(ringWinners, paste0(loc, '/', name, "_ringWinners.csv"))
        write_csv(ringLosers, paste0(loc, '/', name, "_ringLosers.csv"))
        
        # save plot generated above
        ggsave(g1, filename = "Variance vs Max Signal.png",
               width = 8, height = 6)
}

AnalyzeData <- function(time1 = 51, time2 = 39, 
                        filename = "groupNames_LTBI.csv",
                        loc = "plots",
                        fsr = TRUE) {
        GetName()
        AggData(filename = filename, loc = loc, fsr = fsr)
        SubtractControl(ch = 1, cntl = "thermal", loc = loc)
        SubtractControl(ch = 2, cntl = "thermal", loc = loc)
        SubtractControl(ch = "U", cntl = "thermal", loc = loc)
        PlotRingData(cntl = "raw", ch = "U", splitPlot = TRUE, loc = loc)
        PlotRingData(cntl = "thermal", ch = "U", splitPlot = TRUE, loc = loc)
        PlotRingData(cntl = "thermal", ch = 1, splitPlot = FALSE, loc = loc)
        PlotRingData(cntl = "raw", ch = 1, splitPlot = FALSE, loc = loc)
        PlotRingData(cntl = "thermal", ch = 2, splitPlot = FALSE, loc = loc)
        PlotRingData(cntl = "raw", ch = 2, splitPlot = FALSE, loc = loc)
        GetNetShifts(cntl = "thermal", ch = 1, 
                     time1 = time1, time2 = time2, step = 1, loc = loc)
        GetNetShifts(cntl = "thermal", ch = 2,
                     time1 = time1, time2 = time2, step = 1, loc = loc)
        PlotNetShifts(cntl = "thermal", ch = 1, step = 1, loc = loc)
        PlotNetShifts(cntl = "thermal", ch = 2, step = 1, loc = loc)
        # CheckRingQuality(time1 = 20, time2 = 30)
        # shell.exec("https://youtu.be/dQw4w9WgXcQ")
}

AnalyzeAllData <- function() {
        foldersList <- list.dirs(recursive = FALSE)
        directory <- getwd()
        lapply(foldersList, function(i){
                setwd(i)
                AnalyzeData()
                setwd(directory)
        })
}
