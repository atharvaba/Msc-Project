## Clear variables in R directory
rm(list=ls())

## Set working directory
setwd("C:/Users/LENOVO/Downloads/final_project")

## Load libraries
library(dplyr)
library(ggplot2)
library(readr)
library(ggpubr)
theme_set(theme_bw())

####**** CASE ~ HUMAN ~ DD ****####
## Import dataset
case_human_dd <-  read_csv("slimprob_case_human_dd.csv")
head(case_human_dd)

## 1. Create N_Occ_UPC and E_Occ_UPC 
case_human_dd$N_Occ_UPC <- (case_human_dd$N_Occ * case_human_dd$N_UPC) / case_human_dd$N_Seq
case_human_dd$E_Occ_UPC <- (case_human_dd$E_Occ * case_human_dd$E_UPC) / case_human_dd$E_Seq
case_human_dd$Ratio <- case_human_dd$N_Occ_UPC/case_human_dd$E_Occ_UPC
head(case_human_dd)

## 2. Apply Poisson ppois
case_human_dd["pUnd_Occ_UPC"] <-    ppois(q = case_human_dd$N_Occ_UPC, lambda = case_human_dd$E_Occ_UPC, lower.tail = TRUE)
head(case_human_dd)

## Take threshold 0.000125 for dd, 0.00102040816 for md
case_human_dd <- filter(case_human_dd, pUnd_Occ_UPC < 0.05)
head(case_human_dd)

## 3. Apply log10(pUnd_Occ_UPC)
case_human_dd["log_pUnd_Occ_UPC"] <- log10(case_human_dd$pUnd_Occ_UPC)
head(case_human_dd)

## 4. Plot motif v/s log

plot_case_human_dd <- ggplot(case_human_dd, aes(x=Motif, y=log_pUnd_Occ_UPC)) + 
  geom_point(size=3) + 
  geom_segment(aes(x=Motif, 
                   xend=Motif, 
                   y=pUnd_Occ_UPC, 
                   yend=pUnd_Occ_UPC)) + 
  labs(title="Human proteins (cases)", 
       subtitle="DD motifs") + 
  theme(axis.text.x = element_text(angle=65, vjust=0.6))
plot_case_human_dd

####**** CONTROL ~ HUMAN ~ DD ****####
## Import dataset
control_human_dd <-  read_csv("slimprob_control_human_dd.csv")
head(control_human_dd)

## 1. Create N_Occ_UPC and E_Occ_UPC 
control_human_dd$N_Occ_UPC <- (control_human_dd$N_Occ * control_human_dd$N_UPC) / control_human_dd$N_Seq
control_human_dd$E_Occ_UPC <- (control_human_dd$E_Occ * control_human_dd$E_UPC) / control_human_dd$E_Seq
control_human_dd$Ratio <- control_human_dd$N_Occ_UPC/control_human_dd$E_Occ_UPC
head(control_human_dd)

## 2. Apply Poisson ppois
control_human_dd["pUnd_Occ_UPC"] <-    ppois(q = control_human_dd$N_Occ_UPC, lambda = control_human_dd$E_Occ_UPC, lower.tail = TRUE)
head(control_human_dd)

## Take threshold 0.000125 for dd, 0.00102040816 for md
control_human_dd <- filter(control_human_dd, pUnd_Occ_UPC < 0.05)
head(control_human_dd)

## 3. Apply -log10(pUnd_Occ_UPC)
control_human_dd["log_pUnd_Occ_UPC"] <- log10(control_human_dd$pUnd_Occ_UPC)
head(control_human_dd)

## 4. Plot motif v/s log

plot_control_human_dd <- ggplot(control_human_dd, aes(x=Motif, y=log_pUnd_Occ_UPC)) + 
  geom_point(size=3) + 
  geom_segment(aes(x=Motif, 
                   xend=Motif, 
                   y=pUnd_Occ_UPC, 
                   yend=pUnd_Occ_UPC)) + 
  labs(title="Human proteins (controls)", 
       subtitle="DD motifs") + 
  theme(axis.text.x = element_text(angle=65, vjust=0.6))
plot_control_human_dd

####**** CASE ~ HUMAN ~ DE ****####
## Import dataset
case_human_de <-  read_csv("slimprob_case_human_de.csv")
head(case_human_de)

## 1. Create N_Occ_UPC and E_Occ_UPC 
case_human_de$N_Occ_UPC <- (case_human_de$N_Occ * case_human_de$N_UPC) / case_human_de$N_Seq
case_human_de$E_Occ_UPC <- (case_human_de$E_Occ * case_human_de$E_UPC) / case_human_de$E_Seq
case_human_de$Ratio <- case_human_de$N_Occ_UPC/case_human_de$E_Occ_UPC
head(case_human_de)

## 2. Apply Poisson ppois
case_human_de["pUnd_Occ_UPC"] <-    ppois(q = case_human_de$N_Occ_UPC, lambda = case_human_de$E_Occ_UPC, lower.tail = TRUE)
head(case_human_de)

## Take threshold 0.000125 for dd, 0.00102040816 for md
case_human_de <- filter(case_human_de, pUnd_Occ_UPC < 0.05)
head(case_human_de)

## 3. Apply -log10(pUnd_Occ_UPC)
case_human_de["log_pUnd_Occ_UPC"] <- log10(case_human_de$pUnd_Occ_UPC)
head(case_human_de)

## 4. Plot motif v/s log

plot_case_human_de <- ggplot(case_human_de, aes(x=Motif, y=log_pUnd_Occ_UPC)) + 
  geom_point(size=3) + 
  geom_segment(aes(x=Motif, 
                   xend=Motif, 
                   y=pUnd_Occ_UPC, 
                   yend=pUnd_Occ_UPC)) + 
  labs(title="Human proteins (cases)", 
       subtitle="DE motifs") + 
  theme(axis.text.x = element_text(angle=65, vjust=0.6))
plot_case_human_de


####**** CONTROL ~ HUMAN ~ DE ****####
## Import dataset
control_human_de <-  read_csv("slimprob_control_human_de.csv")
head(control_human_de)

## 1. Create N_Occ_UPC and E_Occ_UPC 
control_human_de$N_Occ_UPC <- (control_human_de$N_Occ * control_human_de$N_UPC) / control_human_de$N_Seq
control_human_de$E_Occ_UPC <- (control_human_de$E_Occ * control_human_de$E_UPC) / control_human_de$E_Seq
control_human_de$Ratio <- control_human_de$N_Occ_UPC/control_human_de$E_Occ_UPC
head(control_human_de)

## 2. Apply Poisson ppois
control_human_de["pUnd_Occ_UPC"] <-    ppois(q = control_human_de$N_Occ_UPC, lambda = control_human_de$E_Occ_UPC, lower.tail = TRUE)
head(control_human_de)

## Take threshold 0.000125 for dd, 0.00102040816 for md
control_human_de <- filter(control_human_de, pUnd_Occ_UPC < 0.05)
head(control_human_de)

## 3. Apply -log10(pUnd_Occ_UPC)
control_human_de["log_pUnd_Occ_UPC"] <- log10(control_human_de$pUnd_Occ_UPC)
head(control_human_de)

## 4. Plot motif v/s log

plot_control_human_de <- ggplot(control_human_de, aes(x=Motif, y=log_pUnd_Occ_UPC)) + 
  geom_point(size=3) + 
  geom_segment(aes(x=Motif, 
                   xend=Motif, 
                   y=pUnd_Occ_UPC, 
                   yend=pUnd_Occ_UPC)) + 
  labs(title="Human proteins (controls)", 
       subtitle="DE motifs") + 
  theme(axis.text.x = element_text(angle=65, vjust=0.6))
plot_control_human_de

#### **** JOIN GRAPHS **** ####
final_figure_human <- ggarrange(plot_case_human_dd, plot_control_human_dd, plot_case_human_de, plot_control_human_de,
                                ncol = 2, nrow = 2)
final_figure_human

####**** CASE ~ BOVINE ~ DD ****####
## Import dataset
case_bovine_dd <-  read_csv("slimprob_case_bovine_dd.csv")
head(case_bovine_dd)

## 1. Create N_Occ_UPC and E_Occ_UPC 
case_bovine_dd$N_Occ_UPC <- (case_bovine_dd$N_Occ * case_bovine_dd$N_UPC) / case_bovine_dd$N_Seq
case_bovine_dd$E_Occ_UPC <- (case_bovine_dd$E_Occ * case_bovine_dd$E_UPC) / case_bovine_dd$E_Seq
case_bovine_dd$Ratio <- case_bovine_dd$N_Occ_UPC/case_bovine_dd$E_Occ_UPC
head(case_bovine_dd)

## 2. Apply Poisson ppois
case_bovine_dd["pUnd_Occ_UPC"] <-    ppois(q = case_bovine_dd$N_Occ_UPC, lambda = case_bovine_dd$E_Occ_UPC, lower.tail = TRUE)
head(case_bovine_dd)

## Take threshold 0.000125 for dd, 0.00102040816 for md
case_bovine_dd <- filter(case_bovine_dd, pUnd_Occ_UPC < 0.05)
head(case_bovine_dd)

## 3. Apply -log10(pUnd_Occ_UPC)
case_bovine_dd["log_pUnd_Occ_UPC"] <- log10(case_bovine_dd$pUnd_Occ_UPC)
head(case_bovine_dd)

## 4. Plot motif v/s log

plot_case_bovine_dd <- ggplot(case_bovine_dd, aes(x=Motif, y=log_pUnd_Occ_UPC)) + 
  geom_point(size=3) + 
  geom_segment(aes(x=Motif, 
                   xend=Motif, 
                   y=pUnd_Occ_UPC, 
                   yend=pUnd_Occ_UPC)) + 
  labs(title="Bovine proteins (cases)", 
       subtitle="DD motifs") + 
  theme(axis.text.x = element_text(angle=65, vjust=0.6))
plot_case_bovine_dd

####**** CONTROL ~ BOVINE ~ DD ****####
## Import dataset
control_bovine_dd <-  read_csv("slimprob_control_bovine_dd.csv")
head(control_bovine_dd)

## 1. Create N_Occ_UPC and E_Occ_UPC 
control_bovine_dd$N_Occ_UPC <- (control_bovine_dd$N_Occ * control_bovine_dd$N_UPC) / control_bovine_dd$N_Seq
control_bovine_dd$E_Occ_UPC <- (control_bovine_dd$E_Occ * control_bovine_dd$E_UPC) / control_bovine_dd$E_Seq
control_bovine_dd$Ratio <- control_bovine_dd$N_Occ_UPC/control_bovine_dd$E_Occ_UPC
head(control_bovine_dd)

## 2. Apply Poisson ppois
control_bovine_dd["pUnd_Occ_UPC"] <-    ppois(q = control_bovine_dd$N_Occ_UPC, lambda = control_bovine_dd$E_Occ_UPC, lower.tail = TRUE)
head(control_bovine_dd)

## Take threshold 0.000125 for dd, 0.00102040816 for md
control_bovine_dd <- filter(control_bovine_dd, pUnd_Occ_UPC < 0.05)
head(control_bovine_dd)

## 3. Apply -log10(pUnd_Occ_UPC)
control_bovine_dd["log_pUnd_Occ_UPC"] <- log10(control_bovine_dd$pUnd_Occ_UPC)
head(control_bovine_dd)

## 4. Plot motif v/s log

plot_control_bovine_dd <- ggplot(control_bovine_dd, aes(x=Motif, y=log_pUnd_Occ_UPC)) + 
  geom_point(size=3) + 
  geom_segment(aes(x=Motif, 
                   xend=Motif, 
                   y=pUnd_Occ_UPC, 
                   yend=pUnd_Occ_UPC)) + 
  labs(title="Bovine proteins (controls)", 
       subtitle="DD motifs") + 
  theme(axis.text.x = element_text(angle=65, vjust=0.6))
plot_control_bovine_dd

####**** CASE ~ BOVINE ~ DE ****####
## Import dataset
case_bovine_de <-  read_csv("slimprob_case_bovine_de.csv")
head(case_bovine_de)

## 1. Create N_Occ_UPC and E_Occ_UPC 
case_bovine_de$N_Occ_UPC <- (case_bovine_de$N_Occ * case_bovine_de$N_UPC) / case_bovine_de$N_Seq
case_bovine_de$E_Occ_UPC <- (case_bovine_de$E_Occ * case_bovine_de$E_UPC) / case_bovine_de$E_Seq
case_bovine_de$Ratio <- case_bovine_de$N_Occ_UPC/case_bovine_de$E_Occ_UPC
head(case_bovine_de)

## 2. Apply Poisson ppois
case_bovine_de["pUnd_Occ_UPC"] <-    ppois(q = case_bovine_de$N_Occ_UPC, lambda = case_bovine_de$E_Occ_UPC, lower.tail = TRUE)
head(case_bovine_de)

## Take threshold 0.000125 for dd, 0.00102040816 for md
case_bovine_de <- filter(case_bovine_de, pUnd_Occ_UPC < 0.05)
head(case_bovine_de)

## 3. Apply -log10(pUnd_Occ_UPC)
case_bovine_de["log_pUnd_Occ_UPC"] <- log10(case_bovine_de$pUnd_Occ_UPC)
head(case_bovine_de)

## 4. Plot motif v/s log

plot_case_bovine_de <- ggplot(case_bovine_de, aes(x=Motif, y=log_pUnd_Occ_UPC)) + 
  geom_point(size=3) + 
  geom_segment(aes(x=Motif, 
                   xend=Motif, 
                   y=pUnd_Occ_UPC, 
                   yend=pUnd_Occ_UPC)) + 
  labs(title="Bovine proteins (cases)", 
       subtitle="DE motifs") + 
  theme(axis.text.x = element_text(angle=65, vjust=0.6))
plot_case_bovine_de

####**** CONTROL ~ BOVINE ~ DE ****####
## Import dataset
control_bovine_de <-  read_csv("slimprob_control_bovine_de.csv")
head(control_bovine_de)

## 1. Create N_Occ_UPC and E_Occ_UPC 
control_bovine_de$N_Occ_UPC <- (control_bovine_de$N_Occ * control_bovine_de$N_UPC) / control_bovine_de$N_Seq
control_bovine_de$E_Occ_UPC <- (control_bovine_de$E_Occ * control_bovine_de$E_UPC) / control_bovine_de$E_Seq
control_bovine_de$Ratio <- control_bovine_de$N_Occ_UPC/control_bovine_de$E_Occ_UPC
head(control_bovine_de)

## 2. Apply Poisson ppois
control_bovine_de["pUnd_Occ_UPC"] <-    ppois(q = control_bovine_de$N_Occ_UPC, lambda = control_bovine_de$E_Occ_UPC, lower.tail = TRUE)
head(control_bovine_de)

## Take threshold 0.000125 for dd, 0.00102040816 for md
control_bovine_de <- filter(control_bovine_de, pUnd_Occ_UPC < 0.05)
head(control_bovine_de)

## 3. Apply -log10(pUnd_Occ_UPC)
control_bovine_de["log_pUnd_Occ_UPC"] <- log10(control_bovine_de$pUnd_Occ_UPC)
head(control_bovine_de)

## 4. Plot motif v/s log
theme_set(theme_bw())
plot_control_bovine_de <- ggplot(control_bovine_de, aes(x=Motif, y=log_pUnd_Occ_UPC)) + 
  geom_point(size=3) + 
  geom_segment(aes(x=Motif, 
                   xend=Motif, 
                   y=pUnd_Occ_UPC, 
                   yend=pUnd_Occ_UPC)) + 
  labs(title="Bovine proteins (controls)", 
       subtitle="DE motifs") + 
  theme(axis.text.x = element_text(angle=65, vjust=0.6))
plot_control_bovine_de

#### **** JOIN GRAPHS **** ####
final_figure_bovine <- ggarrange(plot_case_bovine_dd, plot_control_bovine_dd, plot_case_bovine_de, plot_control_bovine_de,
                                 ncol = 2, nrow = 2)
final_figure_bovine


## *** PCA

####**** HUMAN PCA****####
write.csv(case_human_dd, "1_case_human_dd.csv")
write.csv(control_human_dd, "1_control_human_dd.csv")
write.csv(case_human_de, "1_case_human_de.csv")
write.csv(control_human_de, "1_control_human_de.csv")

human_pca <- read_csv("human_pca.csv")

## Apply PCA
human.pca <- data.frame(human_pca)
library(Factoshiny)
result <- Factoshiny(human.pca)

## PCA code (alternative)
res.PCA<-PCA(human.pca[,-c(1)],quali.sup=c(1),quanti.sup=c(5,6,8),graph=FALSE)
plot.PCA(res.PCA,choix='var')
plot.PCA(res.PCA,invisible=c('ind','ind.sup'),habillage=1,title="PCA graph - Homo sapiens",cex=1.1,cex.main=1.1,cex.axis=1.1,label =c('quali'))



####**** BOVINE PCA****####

#Join datasets
write.csv(case_bovine_dd, "1_case_bovine_dd.csv")
write.csv(control_bovine_dd, "1_control_bovine_dd.csv")
write.csv(case_bovine_de, "1_case_bovine_de.csv")
write.csv(control_bovine_de, "1_control_bovine_de.csv")

bovine_pca <- read_csv("bovine_pca.csv")

## Apply PCA
bovine.pca <- data.frame(bovine_pca)
library(Factoshiny)
result <- Factoshiny(bovine.pca)

## PCA code (alternative)
res.PCA<-PCA(bovine.pca[,-c(1)],quali.sup=c(1),quanti.sup=c(5,6,8),graph=FALSE)
plot.PCA(res.PCA,choix='var')
plot.PCA(res.PCA,invisible=c('ind','ind.sup'),habillage=1,title="PCA graph - Bos taurus",cex=1.1,cex.main=1.1,cex.axis=1.1,label =c('quali'))



