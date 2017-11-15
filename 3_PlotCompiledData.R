## Plot all data combined
PlotAllData <- function(data){
        g <- ggplot(data, aes(x = Target, y = NormLog))
        
        allPoint <- g + 
                geom_point(aes(color = factor(TimePoint)),
                           position = "jitter", alpha = 0.7) +
                facet_grid(Treatment ~ CellLine) + 
                theme(axis.text.x = element_text(angle = 45, hjust = 1),
                      legend.position = "bottom") +
                ggtitle("Full Dataset") +
                labs(y = "Normalized Response", color = "Time Point")
        
        ggsave(allPoint, filename = "everything_point.png", 
               width = 20, height = 12)
        
        allBoxplot <- g + 
                geom_boxplot(aes(fill = Treatment)) +
                facet_grid(CellLine~TimePoint) +
                theme(axis.text.x = element_text(angle = 45, hjust = 1),
                      legend.position = "bottom") +
                ggtitle("Full Dataset") +
                ylab("Normalized Response")
        
        ggsave(allBoxplot, filename = "everything_boxplot.png",
               width = 20, height = 12)
        
        allTarget <- ggplot(data, aes(x = Treatment, y = NormLog)) +
                geom_boxplot(aes(fill = interaction(TimePoint, CellLine))) + 
                facet_wrap(~Target) +
                labs(y = "Normalized Response", fill = "") +
                theme(axis.text.x = element_text(angle = 45, hjust = 1),
                      legend.position = "bottom")
        
        ggsave(allTarget, filename = "everything_targetwrap.png",
               width = 20, height = 12)
}

## Plot each target as Net Shift vs Treatment
PlotEachTarget <- function(data){
        targetList <- unique(data$Target)
        
        for(i in targetList) {
                
                tarDat <- filter(data, Target == i & Treatment != "(+)-Serum" & 
                                          Treatment != "(-)-Serum")
                
                g <- ggplot(tarDat, aes(x = Treatment, y = NormLog))
                
                plotName <- unlist(strsplit(i, "/"))[1]
                
                targetPoint <- g + 
                        geom_point(aes(color = TimePoint), 
                                   position = "jitter") + 
                        facet_grid(CellLine~.) +
                        ggtitle(paste0("Target: ", i)) + 
                        labs(y = "Normalized Response", color = "") +
                        theme(axis.text.x = element_text(angle = 45, hjust = 1))
                
                targetBox <- g + 
                        geom_boxplot(aes(fill = TimePoint)) + 
                        facet_grid(CellLine~.) +
                        ggtitle(paste0("Target: ", i)) + 
                        labs(y = "Normalized Response", fill = "") +
                        theme(axis.text.x = element_text(angle = 45, hjust = 1))
                
                targetBox2 <- g + 
                        geom_boxplot(aes(fill = interaction(TimePoint, 
                                                             CellLine))) + 
                        ggtitle(paste0("Target: ", i)) + 
                        labs(y = "Normalized Response", fill = "") +
                        theme(axis.text.x = element_text(angle = 45, hjust = 1))
                
                dir.create("Target Plots", showWarnings = FALSE)
                
                ggsave(targetPoint, 
                       filename = paste0("Target Plots/", 
                                         plotName, "_point.png"), 
                       width = 8, height = 6)
                ggsave(targetBox, 
                       filename = paste0("Target Plots/", 
                                         plotName, "_boxplot.png"),
                       width = 8, height = 6)
                ggsave(targetBox2, 
                       filename = paste0("Target Plots/", 
                                         plotName, "_boxplot_2.png"),
                       width = 8, height = 6)
        }
}

## Plot each treatment as Net Shift vs Target
PlotEachTreatment <- function(data){
        treatmentList <- unique(data$Treatment)
        
        for (i in treatmentList) {
                
                rxDat <- filter(data, Treatment == i)
                
                g <- ggplot(rxDat, aes(x = Target, y = NormLog))
                
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
PlotTreatment <- function(data, control, treatment, targets, cellLine){
        cntlDat <- filter(data, Treatment == control & TimePoint == "1h")
        allDat <- rbind(filter(data, Treatment == treatment), cntlDat)
        allDat$Treatment <- factor(allDat$Treatment, 
                                    levels = c("DMSO", "(-)-Serum", 
                                               "Apitolisib", "Erlotinib",
                                               "Palbociclib", "GNE-317",
                                               "(+)-Serum"))
        
        # Plot all treatments for treatment
        gAll <- ggplot(allDat, 
                        aes(x = interaction(TimePoint, Treatment, Target),
                            y = NormLog, 
                            fill = Target)) +
                geom_boxplot() + 
                facet_grid(CellLine~.) +
                labs(fill = "") + 
                scale_x_discrete(labels = rep(c("0 h", "1 h", "24 h"), 
                                              length(unique(allDat$Target)))) +
                xlab("Treatment Time") + 
                ylab("Normalized Response") +
                ggtitle(paste0("Treatment: ", treatment)) +
                theme(axis.text.x = 
                              element_text(angle = 90, hjust = 1, vjust = 0.5))
        
        ggsave(gAll, 
               filename = paste0("Treatment Plots/", treatment, ".png"), 
               width = 12, height = 8)
        
        # Plot select treatments for treatment
        rxDat <- filter(allDat, Target %in% targets & 
                                 CellLine == cellLine)
        
        txt <- ggplot(rxDat, 
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
TreatmentComp <- function(data, treatments, targets, cellLine){
        rxDat <- filter(data, Treatment == treatments & Target %in% targets &
                                 CellLine == cellLine)
        print(head(rxDat))
        g <- ggplot(rxDat, aes(x = Target, 
                                group = interaction(Treatment, 
                                                    TimePoint,
                                                    Target),
                                y = NormLog))
        
        figComp <- g + 
                geom_boxplot(aes(fill = interaction(Treatment, TimePoint))) + 
                labs(fill = "") + 
                xlab("Target") + 
                ylab("Normalized Response") +
                ggtitle(paste0(treatments[1], " vs ", treatments[2])) +
                theme(axis.text.x = 
                              element_text(angle = 45, hjust = 1))
        
        figComp2 <- g + 
                geom_boxplot(aes(fill = Treatment)) + 
                facet_wrap(~TimePoint) + 
                labs(fill = "") + 
                xlab("Target") + 
                ylab("Normalized Response") +
                ggtitle(paste0(treatments[1], " vs ", treatments[2])) +
                theme(axis.text.x = element_text(angle = 45, hjust = 1))
        
        dir.create(path = "Treatment Comparisons", showWarnings = FALSE)
        
        ggsave(figComp,
               filename = paste0("Treatment Comparisons/", treatments[1], "_", 
                                 treatments[2], "_", cellLine, ".png"), 
               width = 8, height = 6)
        ggsave(figComp2, 
               filename = paste0("Treatment Comparisons/", treatments[1], "_", 
                                 treatments[2], "_", cellLine, "_2.png"), 
               width = 8, height = 6)
}

Cytometry <- function(data){
        dir.create("Cytometry", showWarnings = FALSE)
        rxDat <- filter(data, !grepl("Serum", Treatment))
        pairList <- t(combn(unique(data$Target), 2))
        for(i in seq_len(nrow(pairList))){
                targets <- pairList[i, ]
                cytDat <- filter(rxDat, Target %in% pairList[i, ])
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
        library(biosensor)
        theme_set(theme_few(base_size = 16))
        
        # Load in data to make plots
        x <- compLabelNorm()
        compDat <- x[[1]]
        
        # Save current wd to return to later and setwd to plots folder
        directory <- getwd()
        setwd("../XPP_Plots/")
        
        PlotAllData(data = compDat)
        
        PlotEachTreatment(data = compDat)
        
        PlotEachTarget(data = compDat)
        
        control <- "DMSO"
        txtList <- unique(compDat$Treatment)
        compTargets <- c("pAktSer473", "pS6Ser235/6", "pS6Ser240/4",
                         "pp70S6KThr389", "pRbSer780", "pRbSer807/11")
        
        lapply(txtList, function(i){
                PlotTreatment(data = compDat, control = control, treatment = i,
                              targets = compTargets, cellLine = "GBM6")
                PlotTreatment(data = compDat, control = control, treatment = i,
                              targets = compTargets, cellLine = "GBM26")
        })
        
        # Pairwise treatment list
        txtPairs <- combn(unique(compDat$Treatment), 2, simplify = FALSE)
        
        # Run through pair-wise list to plot treatment comparisons
        lapply(txtPairs, function(i){
                TreatmentComp(data = compDat, treatments = as.vector(i),
                              targets = compTargets, cellLine = "GBM6")
                TreatmentComp(data = compDat, treatments = as.vector(i),
                              targets = compTargets, cellLine = "GBM26")
        })
        setwd(directory)
}
