# Import data ---------------------------------------------------

rm( list = ls() )
library(tidyverse)
library(reshape2)

plotTheme <- theme_classic(base_size = 16) + 
        theme(panel.background = element_rect(fill = "transparent"),
              plot.background = element_rect(fill = "transparent"),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              legend.background = element_rect(fill = "transparent"),
              legend.box.background = element_rect(fill = "transparent"))

theme_set(plotTheme)

setwd("D:/Box Sync/XPP_Data/")
datScaled <- read_csv("compiledNormalized.csv")
datSum <- read_csv("compiledSummed.csv")

# Assess data normality ---------------------------------------------------

shapiro.test(datScaled$LogTransformed)
shapiro.test(datScaled$NormLog)
shapiro.test(datScaled$Normalized)

par(mfrow=c(2,2))
hist(datScaled$Normalized, breaks = 25, main = "Normalized")
hist(datScaled$LogTransformed, breaks = 25, main = "Log Transformed")
hist(datScaled$NormLog, breaks = 25, main = "Normalized & Log Transformed")
hist(rnorm(4994), breaks = 25, main = "Normal Distribution")
dev.print(png,"Histograms of Data Transformations.png", width=8, height=6,
          unit = "in", res = 200)
dev.off()

par(mfrow=c(2,2))
x = datScaled$Normalized
qqnorm(x, main = "Normalized"); qqline(x)
y = datScaled$LogTransformed
qqnorm(y, main = "Log Transformed"); qqline(y)
z = datScaled$NormLog
qqnorm(z, main = "Log Transformed & Normalized"); qqline(z)
qqnorm(rnorm(4994), main = "Normal Distribution"); qqline(z)
dev.print(png,"QQ Plots of Data Transformations.png", width=8, height=6,
          unit = "in", res = 200)
dev.off()

# Pair-wise data comparison -----------------------------------------------

cast <- dcast(data = datScaled, Treatment + CellLine + TimePoint +
                      Replicate + n ~ Target,
              value.var = "NormLog")

cast.gbm6 <- filter(cast, CellLine == "GBM6")
cast.gbm26 <- filter(cast, CellLine == "GBM26")

# plot1 <- ggpairs(cast, aes(color = Treatment, alpha = 0.5), 
#                  columns = c(6:ncol(cast)), 
#                  upper = list(continuous = "density"))
# plot2 <- ggpairs(cast, aes(color = TimePoint, alpha = 0.5), 
#                  columns = c(6:ncol(cast)), 
#                  upper = list(continuous = "density"))
# plot3 <- ggpairs(cast, aes(color = CellLine, alpha = 0.5), 
#                  columns = c(6:ncol(cast)),
#                  upper = list(continuous = "density"))
# 
# ggsave(plot1, filename = "pairs_byTxt.png", width = 22, height = 20)
# ggsave(plot2, filename = "pairs_byTimePoint.png", width = 22, height = 20)
# ggsave(plot3, filename = "pairs_byCellLine.png", width = 22, height = 20)
# 
# 
# cast.alt <- dcast(data = datScaled, Treatment + CellLine + TimePoint +
#                           Replicate + n ~ Target,
#                   value.var = "NetShift")
# 
# plot1.alt <- ggpairs(cast.alt, aes(color = Treatment, alpha = 0.5),
#                      columns = c(6:ncol(cast)),
#                      upper = list(continuous = "density"))
# 
# ggsave(plot1, filename = "pairs_byTxt_alt.png", width = 22, height = 20)


# MANOVA ------------------------------------------------------------------
library(broom)
library(MASS)

targetList <- unique(datScaled$Target)
pairTarList <- t(combn(targetList, 2))

manList <- list()
aovList <- list()
for(i in seq_len(nrow(pairTarList))){
        prot1 <- pairTarList[i, 1]
        prot2 <- pairTarList[i, 2]
        aovRes <- anova(lm(eval(as.symbol(prot1)) ~ 
                                  Treatment + 
                                  CellLine +
                                  TimePoint +
                                  Treatment * CellLine +
                                  Treatment * TimePoint +
                                  CellLine * Treatment +
                                  Treatment * CellLine * TimePoint,
                          data = cast))
        aovResult <- tidy(aovRes)
        aovResult$Protein <- pairTarList[i, 1]
        aovList[[i]] <- aovResult
        
        manRes <- 
                manova(cbind(eval(as.symbol(prot1)),
                             eval(as.symbol(prot2))) ~ 
                               Treatment + 
                               CellLine +
                               TimePoint + 
                               Treatment * CellLine + 
                               Treatment * TimePoint + 
                               CellLine * Treatment + 
                               Treatment * CellLine * TimePoint,
                       data = cast)
        manResult <- tidy(manRes)
        manResult$Protein1 <- prot1
        manResult$Protein2 <- prot2
        manList[[i]] <- manResult
}

# pairwise.t.test(cast$pAktSer473,cast$Treatment)#, p.adjust.method = "bonferroni")

aovDat <- bind_rows(aovList)
aovDat <- aovDat[complete.cases(aovDat), ]
manDat <- bind_rows(manList)

alpha = 0.001
bonferroni <- alpha / nrow(manDat)

manSig <- filter(manDat, p.value < bonferroni)
nrow(manSig) / nrow(manDat)


# Linear Discriminant Analysis --------------------------------------------

castLDA <- cast
castLDA$ID <- with(cast, paste(Treatment, TimePoint, CellLine))
castLDA <- castLDA[complete.cases(castLDA),]
castLDA2 <-castLDA[, -c(1:5)]


fit <- lda(ID ~ ., data=castLDA2)

ldaSVD <- fit$svd^2/sum(fit$svd^2)
ldaScale <- tibble::rownames_to_column(data.frame(fit$scaling))
meltLDAScale <- melt(data = ldaScale,
                measure.vars = c("LD1", "LD2"),
                id.vars = "rowname")

ldaScaling <- ggplot(meltLDAScale, aes(x = rowname, y = value, color = variable)) +
        geom_point(size = 2) + 
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        labs(x = "", y = "Loadings", color = "",
             legend.position = "bottom")
ggsave(plot = ldaScaling, filename = "LDA Scalings.png",
       width = 8, height = 6)


plot(ldaSVD, pch = 19,
     xlab = "Linear Discriminant",
     ylab = "Proportion Explained by LD")
dev.print(png,"SVD Plot from LDA.png", width=8, height=6, 
          unit = "in", res = 200)
dev.off()

fit.values <- predict(fit)
ldaFit <- as.data.frame(fit.values$x)
ldaFit$Treatment <- castLDA$Treatment
ldaFit$CellLine <- castLDA$CellLine
ldaFit$TimePoint <- castLDA$TimePoint

# what about classifying by type of target (e.g., these both hit PI3K)

ldaPlot1 <- ggplot(ldaFit, aes(x = LD1, y = LD2, 
                           color = Treatment,
                           shape = interaction(CellLine, TimePoint))) +
        geom_point(size = 3, alpha = 0.8) +
        labs(color = "", shape = "")
        
ggsave(plot = ldaPlot1, filename = "LDA_ColorByTreatment.png", 
       width = 10, height = 6)

ldaPlot2 <- ggplot(ldaFit, aes(x = LD1, y = LD2, 
                           color = interaction(CellLine, TimePoint))) +
        geom_point(size = 3, alpha = 0.8) +
        labs(color = "") + theme(legend.position = c(0.1, 0.9))

ggsave(plot = ldaPlot2, filename = "LDA_ColorByCLTP.png", 
       width = 10, height = 6)

ldaPlot3 <- ggplot(ldaFit, aes(x = LD1, y = LD2, 
                               color = CellLine)) +
        geom_point(size = 3, alpha = 0.8) +
        labs(color = "") + theme(legend.position = c(0.1, 0.9))

ggsave(plot = ldaPlot3, filename = "LDA_ColorByCL.png", 
       width = 10, height = 6)

ldaPlot4 <- ggplot(ldaFit, aes(x = LD1, y = LD2, 
                               color = TimePoint)) +
        geom_point(size = 3, alpha = 0.8) +
        labs(color = "") + theme(legend.position = c(0.1, 0.9))

ggsave(plot = ldaPlot4, filename = "LDA_ColorByTP.png", 
       width = 10, height = 6)

fit2 <- lda(ID ~ ., data=castLDA2, CV = TRUE)

# Assess the accuracy of the prediction
# percent correct for each category of G
ct <- table(castLDA$ID, fit2$class)
diag(prop.table(ct, 1))
# total percent correct
sum(diag(prop.table(ct)))


# PCA  --------------------------------------------------------------------
library(ggfortify)
library(reshape2)
castLDA3 <- castLDA2
castID <- castLDA2$ID
castLDA3$ID <- NULL
fitPCA <- prcomp(castLDA3)
loading <- tibble::rownames_to_column(data.frame(fitPCA$rotation))
svd <- svd(castLDA3)
pcaSVD <- svd$d^2/sum(svd$d^2)
plot(pcaSVD,
     ylab = "Proportion of Variance Explained",
     xlab = "Principle Component",
     pch = 19)
dev.print(png, "SVD Plot from PCA.png", width=8, height=6, 
          unit = "in", res = 200)
dev.off()

plot(log(pcaSVD)~log(ldaSVD), pch = 19, xlab = "LDA SVD", ylab = "PCA SVD")
qqline(y = x)
dev.print(png, "SVD Comparison.png", width=8, height=6, 
          unit = "in", res = 200)
dev.off()

loadDat <- melt(data = loading,
                measure.vars = c("PC1", "PC2"),
                id.vars = "rowname")

loadPlots <- ggplot(loadDat, aes(x = rowname, y = value, color = variable)) +
        geom_point(size = 2) + 
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        labs(x = "", y = "Loadings", color = "",
             legend.position = "bottom")

ggsave(plot = loadPlots, filename = "PCA Loadings.png",
       width = 8, height = 6)

pcaCL <- autoplot(prcomp(castLDA[, 6:20]), data = castLDA,
                  colour = "CellLine")
pcaTP <- autoplot(prcomp(castLDA[, 6:20]), data = castLDA,
                  colour = "TimePoint")
pcaRX <- autoplot(prcomp(castLDA[, 6:20]), data = castLDA,
                  colour = "Treatment")

ggsave(plot = pcaCL, filename = "PCA by Cell Line.png",
       width = 10, height = 6)
ggsave(plot = pcaTP, filename = "PCA by Time Point.png",
       width = 10, height = 6)
ggsave(plot = pcaRX, filename = "PCA by Treatment.png",
       width = 10, height = 6)

# Hotelling T2 Test -------------------------------------------------------

library(Hotelling)

cellLineList <- unique(datScaled$CellLine)
timePointList <- unique(datScaled$TimePoint)
treatmentList <- unique(datScaled$Treatment)
treatmentList <- treatmentList[!treatmentList=="(-)-Serum"]

combo <- expand.grid(Treatment = treatmentList, 
                   CellLine = cellLineList,
                   TimePoint = timePointList)

target1 <- rep(data.frame(pairTarList)$X1, nrow(combo)) %>% as.character()
target2 <- rep(data.frame(pairTarList)$X2, nrow(combo)) %>% as.character()
combo2 <- replicate(nrow(pairTarList), combo, simplify = FALSE) %>% bind_rows()
combo3 <- cbind(target1, target2, combo2) %>% mutate_all(as.character)
names(combo3) <- c("Target1", "Target2", "Treatment", "CellLine", "TimePoint")

hotel <- list()
for(i in seq_len(nrow(combo3))){
        tmp <- combo3[i,]
        cntl <- filter(datScaled, Treatment == "(-)-Serum" & 
                               TimePoint == combo3[i, 5] &
                               CellLine == combo3[i, 4] &
                                Target %in% combo3[i, c(1,2)])
        castCntl <- dcast(data = cntl, Replicate + n ~ Target,
                          value.var = "NormLog") %>% na.omit()
        castCntl <- castCntl[, -c(1,2)]
        smpl <- filter(datScaled, Treatment == combo3[i, 3] &
                               TimePoint == combo3[i, 5] &
                               CellLine == combo3[i, 4] &
                               Target %in% combo3[i, c(1,2)])
        castSmpl <- dcast(data = smpl, Replicate + n ~ Target,
                          value.var = "NormLog") %>% na.omit()
        castSmpl <- castSmpl[, -c(1,2)]
        test <- hotelling.test(castSmpl, castCntl, shrinkage = TRUE)
        tmp$pval <- test$pval
        tmp$stat <- test$stats$statistic
        tmp$df1 <- test$stats$df[1]
        hotel[[i]] <- tmp
}

hotel <- bind_rows(hotel)

# attach(datScaled)
# interaction.plot(x.factor = interaction(Treatment, TimePoint, CellLine), 
#                  trace.factor = Target, col = c(1:12), lty = 1,
#                  response = NormLog)

# ggplot(datSum, aes(x = interaction(Treatment, TimePoint,  CellLine),
#                   y = NormLog_mean, group = Target, color = Target)) +
#        geom_line() + coord_flip() + labs(y = "Z-score", x = "")

# Protein Pairs -----------------------------------------------------------

pairsList <- t(combn(targetList, 2))
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
        # names(dat) <- c("Difference of means", "Std Error", "t", "p-value")
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


# Calculate Score ---------------------------------------------------------

casting <- dcast(datScaled, TimePoint + CellLine + Treatment + Target ~
                         Replicate + n, 
                 value.var = "NormLog", fun.aggregate = mean)

rownames(casting) = paste(casting$Treatment,
                          casting$CellLine,
                          casting$TimePoint,
                          casting$Target)

mcast <-data.matrix(casting[, -c(1,2,3,4)])

mcastCor <- rcorr(t(mcast))

datCor <- mcastCor$r

correlation <- vector()

for (i in seq_len(nrow(paired))){
        txt <- paired[i, 3]
        cl <- paired[i, 4]
        tp <- paired[i, 5]
        tar1 <- paired[i, 1]
        tar2 <- paired[i, 2]
        whichRow <- paste(txt, cl, tp, tar1)
        whichCol <- paste(txt, cl, tp, tar2)
        
        correlation[i] <- datCor[whichRow, whichCol]
        # corList[[i]] <- correlation 
        # names(corList) <- "Correlation"
        
        # print(paste(txt, cl, tp))
}

updated <- cbind(paired, correlation)
