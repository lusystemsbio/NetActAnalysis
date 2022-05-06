library(DescTools)

library(reshape2)
library(ggplot2)
FANTOM_20 <- readRDS("ROC_FONTOM20.RDS")
colnames(FANTOM_20) <- c("Specificity", "Sensitivity")
FANTOM_50 <- readRDS("ROC_FONTOM50.RDS")
colnames(FANTOM_50) <- c("Specificity", "Sensitivity")
FANTOM_100 <- readRDS("ROC_FONTOM100.RDS")
colnames(FANTOM_100) <- c("Specificity", "Sensitivity")
Literature <- readRDS("ROC_hDB_v2.RDS")
colnames(Literature) <- c("Specificity", "Sensitivity")
Rcis10perMotif <- readRDS("ROC_hgs_Rcistarget_500bp-upstream_Target10_all.RDS")
colnames(Rcis10perMotif) <- c("Specificity", "Sensitivity")
Rcis50 <- readRDS("ROC_hgs_Rcistarget_50.RDS")
colnames(Rcis50) <- c("Specificity", "Sensitivity")
ChIP <- readRDS("ROC_ChIPDB.RDS")
colnames(ChIP) <- c("Specificity", "Sensitivity")
default <- data.frame(Specificity = seq(0,1,by=0.01),Sensitivity = seq(1,0,by=-0.01))
all <- reshape2::melt(default)
group <- rep("Random", 2*nrow(default))

all <-  melt(FANTOM_20)
group <- c(rep("FANTOM_20 (0.49)", 2*nrow(FANTOM_20)))

all <- rbind(all, melt(FANTOM_50))
group <- c(group,rep("FANTOM_50 (0.49)", 2*nrow(FANTOM_50)))

all <- rbind(all, melt(FANTOM_100))
group <- c(group,rep("FANTOM_100 (0.54)", 2*nrow(FANTOM_100)))

all <- rbind(all, melt(Rcis50))
group <-c(group,rep("Rcis50 (0.60)", 2*nrow(Rcis50)))

all <- rbind(all, melt(Rcis10perMotif))
group <- c(group,rep("Rcis10perMotif (0.65)", 2*nrow(Rcis10perMotif)))

all <- rbind(all, melt(ChIP))
group <- c(group,rep("ChIP (0.61)", 2*nrow(ChIP)))

all <- rbind(all, melt(Literature))
group <- c(group,rep("Literature (0.74)", 2*nrow(Literature)))

tmp <- cbind(group, all)
tmp2 <- tmp[tmp$variable=="Specificity",]
tmp3 <- tmp[tmp$variable=="Sensitivity",]
tmp4 <- cbind(tmp2,tmp3)
tmp4 <- tmp4[,c(-2,-4,-5)]
colnames(tmp4) <- c("Database", "Specificity", "Sensitivity")

ggplot(data = tmp4) +
  geom_line(data = default, aes(x=Specificity, y=Sensitivity), linetype = 2, show.legend = FALSE) +
  geom_line(aes(x=Specificity, y=Sensitivity, color = Database )) +
  labs(x="Specificity", y ="Sensitivity") +
  scale_x_reverse() +
  theme_bw() +
  guides(fill=guide_legend(title="New Legend Title")) +
  theme(text = element_text(size=20),
        legend.position = c(0.7, 0.2))



