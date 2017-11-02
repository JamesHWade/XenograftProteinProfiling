## Plot all data combined
PlotAllData <- function(){
        g <- ggplot(dat, aes(x = Target, y = NormLog))
        
        all.point <- g + 
                geom_point(aes(color = factor(Replicate)),
                           position = "jitter", alpha = 0.7) +
                facet_grid(Treatment ~ interaction(TimePoint, CellLine)) + 
                theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
                ggtitle("Full Dataset")
        
        
        ggsave(all.point, filename = "everything_point.png", 
               width = 20, height = 16)
        
        all.boxplot <- g + 
                geom_boxplot(aes(fill = Treatment)) +
                facet_grid(Treatment ~ interaction(TimePoint, CellLine)) +
                theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
                ggtitle("Full Dataset") +
                ylab("Normalized Response")
        
        ggsave(all.boxplot, filename = "everything_boxplot.png",
               width = 20, height = 16)
}

## Plot each target as Net Shift vs Treatment
PlotEachTarget <- function(){
        targetList <- unique(dat$Target)
        
        for(i in targetList) {
                
                dat.tar <- filter(dat, Target == i & Treatment != "(+)-Serum" & 
                                          Treatment != "(-)-Serum")
                
                g <- ggplot(dat.tar, aes(x = Treatment, y = NormLog))
                
                plotName <- unlist(strsplit(i, "/"))[1]
                
                target.point <- g + 
                        geom_point(aes(color = Treatment), 
                                   position = "jitter") + 
                        facet_grid(CellLine~TimePoint) +
                        ggtitle(paste0("Target: ", i)) + 
                        ylab("Normalized Response") +
                        theme(axis.text.x = element_text(angle = 45, hjust = 1))
                
                target.box <- g + 
                        geom_boxplot(aes(fill = Treatment)) + 
                        facet_grid(CellLine~TimePoint) +
                        ggtitle(paste0("Target: ", i)) + 
                        ylab("Normalized Response") +
                        theme(axis.text.x = element_text(angle = 45, hjust = 1))
                
                dir.create("Target Plots", showWarnings = FALSE)
                
                ggsave(target.point, 
                       filename = paste0("Target Plots/", 
                                         plotName, "_point.png"), 
                       width = 12, height = 8)
                
                ggsave(target.box, 
                       filename = paste0("Target Plots/", 
                                         plotName, "_boxplot.png"),
                       width = 12, height = 8)
        }
}

## Plot each treatment as Net Shift vs Target
PlotEachTreatment <- function(){
        treatmentList <- unique(dat$Treatment)
        
        for (i in treatmentList) {
                
                dat.rx <- filter(dat, Treatment == i)
                
                g <- ggplot(dat.rx, aes(x = Target, y = NormLog))
                
                fig4 <- g + 
                        geom_point(aes(color = Target), position = "jitter") + 
                        facet_grid(CellLine~TimePoint) +
                        ggtitle(paste0("Treatment: ", i)) +
                        ylab("Normalized Response") +
                        theme(axis.text.x = element_text(angle = 45, hjust = 1))
                
                fig5 <- g + 
                        geom_boxplot(aes(fill = Target)) + 
                        facet_grid(CellLine~TimePoint) +
                        ggtitle(paste0("Treatment: ", i)) +
                        ylab("Normalized Response") +
                        theme(axis.text.x = element_text(angle = 45, hjust = 1))
                
                dir.create("Treatment Plots", showWarnings = FALSE)
                
                ggsave(fig4, filename = paste0("Treatment Plots/",
                                               i, "_point.png"),
                       width = 12, height = 8)
                
                ggsave(fig5, filename = paste0("Treatment Plots/",
                                               i, "_boxplot.png"),
                       width = 12, height = 8)
        }
}

## Plot Treatments
PlotTreatment <- function(control, treatment, targets, cellLine){
        dat.cntl <- filter(dat, Treatment == control & TimePoint == "1h")
        dat.all <- rbind(filter(dat, Treatment == treatment), dat.cntl)
        dat.all$Treatment <- factor(dat.all$Treatment, 
                                    levels = c("DMSO", "(-)-Serum", 
                                               "Apitolisib", "Erlotinib",
                                               "Palbociclib", "GNE-317",
                                               "(+)-Serum"))
        
        # Plot all treatments for treatment
        g.all <- ggplot(dat.all, 
                        aes(x = interaction(TimePoint, Treatment, Target),
                            y = NormLog, 
                            fill = Target)) +
                geom_boxplot() + 
                facet_grid(CellLine~.) +
                labs(fill = "") + 
                scale_x_discrete(labels = rep(c("0 h", "1 h", "24 h"), 
                                              length(unique(dat.all$Target)))) +
                xlab("Treatment Time") + 
                ylab("Normalized Response") +
                ggtitle(paste0("Treatment: ", treatment)) +
                theme(axis.text.x = 
                              element_text(angle = 90, hjust = 1, vjust = 0.5))
        
        ggsave(g.all, 
               filename = paste0("Treatment Plots/", treatment, ".png"), 
               width = 12, height = 8)
        
        # Plot select treatments for treatment
        dat.rx <- filter(dat.all, Target %in% targets & 
                                 CellLine == cellLine)
        
        txt <- ggplot(dat.rx, 
                      aes(x = interaction(TimePoint, Treatment, 
                                          Target, CellLine), 
                          y = NormLog, 
                          fill = Target)) +
                geom_boxplot() + 
                labs(fill = "") + 
                scale_x_discrete(labels = rep(c("0 h", "1 h", "24 h"), 
                                              length(targets))) +
                xlab("Treatment Time") + 
                ylab("Normalized Response") +
                ggtitle(paste("Treatment:", treatment, cellLine)) +
                theme(axis.text.x = 
                              element_text(angle = 90, hjust = 1, vjust = 0.5))
        
        ggsave(txt, 
               filename = paste0("Treatment Plots/", treatment, "_", 
                                 cellLine, ".png"),
               width = 8, height = 6)
}

## Pairwise treatment comparisons
TreatmentComp <- function(treatments, targets, cellLine){
        dat.rx <- filter(dat, Treatment == treatments & Target %in% targets &
                                 CellLine == cellLine)
        print(head(dat.rx))
        g <- ggplot(dat.rx, aes(x = Target, 
                                group = interaction(Treatment, 
                                                    TimePoint,
                                                    Target),
                                y = NormLog))
        
        fig_comp <- g + 
                geom_boxplot(aes(fill = interaction(Treatment, TimePoint))) + 
                labs(fill = "") + 
                xlab("Target") + 
                ylab("Normalized Response") +
                ggtitle(paste0(treatments[1], " vs ", treatments[2])) +
                theme(axis.text.x = 
                              element_text(angle = 45, hjust = 1))
        
        fig_comp.2 <- g + 
                geom_boxplot(aes(fill = Treatment)) + 
                facet_wrap(~TimePoint) + 
                labs(fill = "") + 
                xlab("Target") + 
                ylab("Normalized Response") +
                ggtitle(paste0(treatments[1], " vs ", treatments[2])) +
                theme(axis.text.x = element_text(angle = 45, hjust = 1))
        
        dir.create(path = "Treatment Comparisons", showWarnings = FALSE)
        
        ggsave(fig_comp,
               filename = paste0("Treatment Comparisons/", treatments[1], "_", 
                                 treatments[2], "_", cellLine, ".png"), 
               width = 8, height = 6)
        ggsave(fig_comp.2, 
               filename = paste0("Treatment Comparisons/", treatments[1], "_", 
                                 treatments[2], "_", cellLine, "_2.png"), 
               width = 8, height = 6)
}

Cytometry <- function(){
        dir.create("Cytometry", showWarnings = FALSE)
        dat.rx <- filter(dat, !grepl("Serum", Treatment))
        pairList <- t(combn(unique(dat$Target), 2))
        for(i in seq_len(nrow(pairList))){
                targets <- pairList[i, ]
                cytDat <- filter(dat.rx, Target %in% pairList[i, ])
                cytCast <- dcast(data = cytDat, Treatment + CellLine + 
                                         TimePoint + Replicate + n ~ Target,
                                 value.var = "NormLog")
                
                plot <- ggplot(cytCast, aes(x = cytCast[, 6], 
                                    y = cytCast[, 7],
                                    color = interaction(CellLine, TimePoint,
                                                        Treatment))) +
                                    # color = TimePoint)) +
                        geom_point() +
                        labs(color = "", x = pairList[i,1], y = pairList[i,2])
                
                tar1 <- substr(pairList[i,1], 1, 9)
                tar2 <- substr(pairList[i,2], 1, 9)
                ggsave(plot, filename = paste0("Cytometry/", tar1, " vs ", 
                                               tar2, ".png"),
                       width = 8, height = 6)
        }
}

## Run all of the above functions to generate plots
PlotData <- function(){
        # Load libraries and set theme for all plots
        library(tidyverse)
        library(ggthemes)
        theme_set(theme_few(base_size = 16))
        
        # Load in data to make plots
        setwd("D:/Box Sync/Data/")
        dat <<- read_csv("compiledNormalized.csv")
        
        # Save current wd to return to later and setwd to plots folder
        directory <- getwd()
        setwd("../XPP_Plots/")
        
        PlotAllData()
        
        PlotEachTreatment()
        
        PlotEachTarget()
        
        control <- "DMSO"
        txtList <- unique(dat$Treatment)
        compTargets <- c("pAktSer473", "pS6Ser235/6", "pS6Ser240/4",
                         "pp70S6KThr389", "pRbSer780", "pRbSer807/11")
        
        lapply(txtList, function(i){
                PlotTreatment(control = control, treatment = i,
                              targets = compTargets, cellLine = "GBM6")
                PlotTreatment(control = control, treatment = i,
                              targets = compTargets, cellLine = "GBM26")
        })
        
        # Pairwise treatment list
        txtPairs <- combn(unique(dat$Treatment), 2, simplify = FALSE)
        
        # Run through pair-wise list to plot treatment comparisons
        lapply(txtPairs, function(i){
                TreatmentComp(treatments = as.vector(i), targets = compTargets,
                              cellLine = "GBM6")
                TreatmentComp(treatments = as.vector(i), targets = compTargets,
                              cellLine = "GBM26")
        })
        setwd(directory)
}
