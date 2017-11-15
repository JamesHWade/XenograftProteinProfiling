combineNetShifts <- function(){
        netList <- grep('net', 
                        list.files(pattern = '.csv', recursive = TRUE), 
                        value = TRUE)
        
        netShifts <- lapply(netList, function(i){
                netShift <- read.csv(i)
                # convert Target and Experiment to character type to avoid
                # warnings when using bind_rows function below
                netShift$Target <- as.character(netShift$Target)
                netShift$Experiment <- as.character(netShift$Experiment)
                netShift
        })
        
        netComb <- dplyr::bind_rows(netShifts)
        netComb <- dplyr::filter(netComb, 
                                 !grepl("thermal|Ignore|Control", Target))
        
        netCombAvg <- dplyr::group_by(netComb, Target, Concentration)
        netCombAvg <- dplyr::summarise_at(netCombAvg, 
                                          dplyr::vars(NetShift),
                                          dplyr::funs(mean, stats::sd))

        readr::write_csv(netComb, path = "combinedNetShifts.csv")
        readr::write_csv(netCombAvg, path = "combinedNetShifts_Avg.csv")
        return(netComb)
}

plotCombineNetShifts <- function(data){
        ggplot2::theme_set(ggthemes::theme_few(base_size = 16))
        
        plot <- ggplot2::ggplot(data, 
                                ggplot2::aes(x = Concentration,
                                             y = NetShift,
                                             color = Target,
                                             group = Concentration)) +
                ggplot2::ggtitle(name) + 
                ggplot2::scale_x_log10() +
                ggplot2::labs(x = "Analyte Concentration (pg/mL)",
                              y = expression(paste("Relative Shift (",
                                                   Delta,"pm)")),
                              color = "Target")
        
        plot1 <- plot + ggplot2::geom_point()
        plot1a <- plot1 + ggplot2::facet_wrap(~Target)
        plot2 <- plot + ggplot2::geom_boxplot()
        plot2a <- plot2 + ggplot2::facet_wrap(~Target)

        ggplot2::ggsave(plot = plot1,
                        filename = paste0("netShiftsCombined_point.png"),
                        width = 16, height = 12)
        ggplot2::ggsave(plot = plot1a,
                        filename = paste0("netShiftsCombined_pointwrap.png"),
                        width = 16, height = 12)
        ggplot2::ggsave(plot = plot2,
                        filename = paste0("netShiftsCombined_box.png"),
                        width = 16, height = 12)
        ggplot2::ggsave(plot = plot2a,
                        filename = paste0("netShiftsCombined_boxwrap.png"),
                        width = 16, height = 12)
}

fitCalCurves <- function(data, loc = 'plots', tarList){
        ggplot2::theme_set(ggthemes::theme_few(base_size = 16))
        
        tarList <- unique(data$Target)
        fit <- list()
        for(i in seq_len(length(tarList))) {
                tar <- tarList[i]
                tarDat <- dplyr::filter(dat, Target == tar)
                y <- tarDat$NetShift
                x <- tarDat$Concentration
                print(tar)
                runNext <- TRUE
                tarFit <- tryCatch({fit.info <- nls(formula = y ~ A + (B - A) / 
                                                            (1 + (x / C) ^ D),
                                start = list(A = 10000,
                                             B = 100,
                                             C = 4000,
                                             D = 1))
                },
                         error = function(e) {"failed"},
                         finally = print(tar))
                if(tarFit[1] != "failed"){
                        fit[[i]] <- broom::tidy(fit.info)
                        fit[[i]]$Target <- unique(tarDat$Target)
                        A <- as.numeric(coef(fit.info)[1])
                        B <- as.numeric(coef(fit.info)[2])
                        C <- as.numeric(coef(fit.info)[3])
                        D <- as.numeric(coef(fit.info)[4])
                        
                        testFun <- function(x) {A + (B - A) / (1 + (x / C) ^ D)}
                        
                        plot <- ggplot2::ggplot(tarDat,
                                                  ggplot2::aes(x = Concentration,
                                                               y = NetShift,
                                                               group = Concentration)) +
                                ggplot2::geom_boxplot(fill = "red") +
                                ggplot2::stat_function(fun = testFun,
                                                       color = "blue", size = 1) +
                                ggplot2::scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                                                       labels = scales::trans_format("log10",
                                                                                     scales::math_format(10 ^ .x))) +
                                ggplot2::labs(x = "Analyte Concentration (pg/mL)",
                                              y = expression(paste("Relative Shift (",
                                                                   Delta,
                                                                   "pm)"))) +
                                ggplot2::annotation_logticks()
                        
                        ggplot2::ggsave(plot,
                                        filename = paste0(tar, "CalCurve.png"),
                                        width = 8, height = 6)
                }
        }
        
        fit <- dplyr::bind_rows(fit)
        
        capture.output(fit, file = "fitInfo.txt")
        readr::write_csv(fit, path = "fitInfo.csv")
}