# Import data ---------------------------------------------------

rm( list = ls() )
library(tidyverse)
library(reshape2)

setwd("D:/Box Sync/Data/")
datScaled <- read_csv("compiledNormalized.csv")
datScaled <- filter(datScaled, Target != "h-HIF1aPro564" &
                            Target != "pp53Ser15" &
                            Target != "pc-AblTyr204" &
                            Target != "pGSK3bSer9")
datSum <- read_csv("compiledSummed.csv")

# Assess data normality ---------------------------------------------------

shapiro.test(datScaled$LogTransformed)
shapiro.test(datScaled$NormLog)

par(mfrow=c(2,2))
hist(datScaled$Normalized, breaks = 25, main = "Normalized")
hist(datScaled$LogTransformed, breaks = 25, main = "Log Transformed")
hist(datScaled$NormLog, breaks = 25, main = "Normalized & Log Transformed")
hist(rnorm(4994), breaks = 25, main = "Normal Distribution")
dev.print(png,"Histograms of Data Transformations.png", width=800, height=600)
dev.off()

par(mfrow=c(2,2))
x = datScaled$Normalized
qqnorm(x, main = "Normalized"); qqline(x)
y = datScaled$LogTransformed
qqnorm(y, main = "Log Transformed"); qqline(y)
z = datScaled$NormLog
qqnorm(z, main = "Log Transformed & Normalized"); qqline(z)
qqnorm(rnorm(4994), main = "Normal Distribution"); qqline(z)
dev.print(png,"QQ Plots of Data Transformations.png", width=1200, height=800)
dev.off()

# Pair-wise data comparison -----------------------------------------------

library(GGally)
library(ggthemes)
theme_set(theme_few(base_size = 14))

cast <- dcast(data = datScaled, Treatment + CellLine + TimePoint +
                      Replicate + n ~ Target,
              value.var = "NormLog")

cast.gbm6 <- filter(cast, CellLine == "GBM6")
cast.gbm26 <- filter(cast, CellLine == "GBM26")

plot1 <- ggpairs(cast, aes(color = Treatment, alpha = 0.5), 
                 columns = c(6:ncol(cast)), 
                 upper = list(continuous = "density"))
plot2 <- ggpairs(cast, aes(color = TimePoint, alpha = 0.5), 
                 columns = c(6:ncol(cast)), 
                 upper = list(continuous = "density"))
plot3 <- ggpairs(cast, aes(color = CellLine, alpha = 0.5), 
                 columns = c(6:ncol(cast)),
                 upper = list(continuous = "density"))

ggsave(plot1, filename = "pairs_byTxt.png", width = 22, height = 20)
ggsave(plot2, filename = "pairs_byTimePoint.png", width = 22, height = 20)
ggsave(plot3, filename = "pairs_byCellLine.png", width = 22, height = 20)


cast.alt <- dcast(data = datScaled, Treatment + CellLine + TimePoint +
                          Replicate + n ~ Target,
                  value.var = "NetShift")

plot1.alt <- ggpairs(cast.alt, aes(color = Treatment, alpha = 0.5),
                     columns = c(6:ncol(cast)),
                     upper = list(continuous = "density"))

ggsave(plot1, filename = "pairs_byTxt_alt.png", width = 22, height = 20)


# MANOVA ------------------------------------------------------------------
library(broom)

targetList <- unique(datScaled$Target)
pairTarList <- t(combn(targetList, 2))

manList <- list()

for(i in seq_len(nrow(pairTarList))){
        pair <- paste(pairTarList[i, 1], pairTarList[i, 2])
        manResult <- 
                manova(cbind(eval(as.symbol(pairTarList[i, 1])),
                             eval(as.symbol(pairTarList[i, 2]))) ~ 
                               Treatment + 
                               CellLine +
                               TimePoint + 
                               Treatment * CellLine + 
                               Treatment * TimePoint + 
                               CellLine * Treatment + 
                               Treatment * CellLine * TimePoint,
                       data = cast) %>% tidy()
        manResult$Protein1 <- pairTarList[i, 1]
        manResult$Protein2 <- pairTarList[i, 2]
        manList[[i]] <- manResult
}

manDat <- bind_rows(manList)

# Hotelling T2 Test -------------------------------------------------------

library(Hotelling)

cellLineList <- unique(datScaled$CellLine)
timePointList <- unique(datScaled$TimePoint)
treatmentList <- unique(datScaled$Treatment)
treatmentList <- treatmentList[!treatmentList=="(-)-Serum"]

exp <- expand.grid(Treatment = treatmentList, 
                      CellLine = cellLineList,
                      TimePoint = timePointList)

hotel <- list()
for(i in seq_len(nrow(exp))){
        tmp <- exp[i,]
        cntl <- filter(datScaled, Treatment == "(-)-Serum" & 
                               TimePoint == exp[i, 3] &
                               CellLine == exp[i, 2])
        castCntl <- dcast(data = cntl, Replicate + n ~ Target,
                          value.var = "NetShift") %>% na.omit()
        castCntl <- castCntl[, -c(1,2)]
        smpl <- filter(datScaled, Treatment == exp[i, 1] &
                               CellLine == exp[i, 2] &
                               TimePoint == exp[i, 3])
        castSmpl <- dcast(data = smpl, Replicate + n ~ Target,
                          value.var = "NetShift") %>% na.omit()
        castSmpl <- castSmpl[, -c(1,2)]
        test <- hotelling.test(castSmpl, castCntl, TRUE)
        tmp$pval <- test$pval
        tmp$stat <- test$stats$statistic
        tmp$df1 <- test$stats$df[1]
        hotel[[i]] <- tmp
}

hotel <- bind_rows(hotel)

attach(datScaled)
interaction.plot(x.factor = interaction(Treatment, TimePoint, CellLine), 
                 trace.factor = Target, col = c(1:12), lty = 1,
                 response = NormLog)

ggplot(datSum, aes(x = interaction(Treatment, TimePoint, CellLine),
                   y = NormLog_mean, group = Target, color = Target)) +
        geom_line() + coord_flip()

# Protein Pairs -----------------------------------------------------------

pairsList <- t(combn(targetList, 2))

cellLine <- unique(dat$CellLine)
treatmentList <- unique(dat$Treatment)
timePoint <- unique(dat$TimePoint)
combo <- expand.grid(Treatment = treatmentList, 
                     CellLine = cellLine,
                     TimePoint = timePoint)

pairedDat <- list()

for(i in seq_len(nrow(combo))){
        txt <- combo[i, 1]
        cl <- combo[i, 2]
        tp <- combo[i, 3]
        pairDat <- filter(datSum, Treatment == txt & CellLine == cl &
                                  TimePoint == tp)
        print(paste0("i = ", i, "/", nrow(combo)))
        tmp <- data.frame()
        for(j in seq_len(nrow(pairsList))){
                print(paste0(j, "/", nrow(pairsList)))
                prot1 <- pairsList[j, 1]
                prot2 <- pairsList[j, 2]
                # print(paste(prot1, prot2))
                p1 <- filter(pairDat, Target == prot1)
                p2 <- filter(pairDat, Target == prot2)
                tmp[j, 1] <- prot1
                tmp[j, 2] <- prot2
                tmp[j, 3] <- txt
                tmp[j, 4] <- cl
                tmp[j, 5] <- tp
                tmp[j, 6] <- p1$NormLog_mean - p2$NormLog_mean
                tmp[j, 7] <- sqrt((p1$NormLog_sd)^2 + (p2$NormLog_sd)^2)
                tmp[j, 8] <- p1$NormLog_length
                tmp[j, 9] <- p2$NormLog_length
        }
        pairedDat[[i]] <- tmp
}

paired <- bind_rows(pairedDat)

names(paired) <- c("Target_1", "Target_2", "Treatment", "CellLine",
                    "TimePoint", "Mean", "SD", "N_1", "N_2")


# T-tests -----------------------------------------------------------

t.test2 <- function(m1,m2,s1,s2,n1,n2,m0=0,equal.variance=FALSE)
{
        if( equal.variance==FALSE ) 
        {
                se <- sqrt( (s1^2/n1) + (s2^2/n2) )
                # welch-satterthwaite df
                df <- ( (s1^2/n1 + s2^2/n2)^2 )/
                        ( (s1^2/n1)^2/(n1-1) + (s2^2/n2)^2/(n2-1) )
        } else
        {
                # pooled standard deviation, scaled by the sample sizes
                se <- sqrt( (1/n1 + 1/n2) * ((n1-1)*s1^2 + (n2-1)*s2^2)/
                                    (n1+n2-2) ) 
                df <- n1+n2-2
        }      
        t <- (m1-m2-m0)/se 
        dat <- c(m1-m2, se, t, 2*pt(-abs(t),df))    
        names(dat) <- c("Difference of means", "Std Error", "t", "p-value")
        return(dat) 
}

ttests <- list()
tmp <- data.frame()
for(i in seq_len(nrow(paired))){
        print(i)
        datRow <- paired[i,]
        cntlRow <- filter(paired, Treatment == "(-)-Serum", 
                          CellLine == datRow$CellLine,
                          Target_1 == datRow$Target_1,
                          Target_2 == datRow$Target_2,
                          TimePoint == datRow$TimePoint)
        comp <- t.test2(m1 = datRow$Mean, m2 = cntlRow$Mean,
                            s1 = datRow$SD, s2 = cntlRow$SD,
                            n1 = datRow$N_1, n2 = cntlRow$N_1)
        tmp[1, 1] <- comp[3]
        tmp[1, 2] <- comp[4]
        ttests[[i]] <- tmp
}

ttests <- bind_rows(ttests)

pairedTests <- cbind(paired, ttests)

pval <- 0.05
nprot <- length(unique(pairedTests$Target_1))
ntests <- factorial(nprot) / (2 * factorial(nprot-2))
alpha <- pval / ntests

