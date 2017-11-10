### 1_ProcessRawRingData.R #####################################################
# 
# The purpose of this program is to process with the raw data from the 
# Maverick M1 detection system (Genalyte, Inc., San Diego, CA) and output simple
# line graphs, bar charts, and box plots. The functions also generate companion
# csv files containning processed the prcoessed data for subsequent analysis.
# The folder containing output from the M1 typically consists of:
#       1 a csv file for each ring and 
#       2 a comments file the describes the experimental run
# The comments file is not needed for this program. In addition to the csv files
# for each ring, a separate file containing the chip layout is required. An
# example of a chip layout file is provided in the "BaileyLabMRRs" repository
# located at https://github.com/BaileyLabUM/BaileyLabMRRs. See 
# "groupNames_allClusters.csv" for an example.
# 
# Note: This version of the software is optimized for the Bailey lab's HRP 
# assay. See dx.doi.org/10.1021/acscentsci.5b00250 for a description. However,
# input variables can be altered to accomodate many alternative experiments.
# 
# The following libraries are used in this program: tidyverse, ggthemes.
# To install these packages, run the following code in the console:
#       install.packages(c("tidyverse", "ggthemes"))
# 
# To use the program:
#       1. Ensure the you have the necessary libraries installed. See note above
#               for instructions on installing libraries.
#       2. Copy the chip layout file (e.g., "groupNames_XPP.csv") into the 
#               directory containing the raw ring data.
#               Note: This program has the highest chance of success if the
#               directoy only contains:
#                       1. raw ring data files (e.g., "03.csv") and
#                       2. the chip layout file (e.g., "groupNames_XPP.csv")
#       3. Source all of the code from this file. There are multiple ways, but
#               one method is to click the `Source` button on the window of the 
#               Source Console in RStudio after opening this file in RStudio.
#               
#       4. Set the working directory to the folder containing the raw data and
#              chip layout for example:
#              setwd("C:/Users/USERNAME/Documents/CHIPNAME_gaskGASK_DATE")
#              
#       5. Execute the code by running the AnalyzeData function. This function
#               requires 5 input variable:
#                       1. filename - the name of the chip layout file
#                       2. loc - the name of folder to store generated files
#                       3. fsr - a logical variable the corrects for free 
#                               spectral range shifts; most likely leave as
#                               FALSE if you don't know what this is
#                       4. time1 - the later time for net shift measurements
#                       5. time2 - the earlier time for net shift measurments
#               Note: to calculate net shift measurements, the relative shift
#               at time2 is subtracted from time1 (netshift = time1 - time2).
#               Here is an example of code to run:
#                      AnalyzeData(filename = "groupNames_XPP.csv", 
#                                  time1 = 51, time2 = 39,
#                                  loc = "plots",
#                                  fsr = FALSE)
#
################################################################################

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
                        fsr = FALSE) {
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
