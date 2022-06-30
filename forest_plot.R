######### TABLE #######
library(meta)
library(vegan)
library(randomForest)
library(MicrobeR)
library(tidyverse)


getwd()
setwd("gitignore")

################ IMPORT DATA ###########

metadata <- read.table("metadata.tsv", header = TRUE, sep = "\t")


#### + alpha diversity ####
alpha_div <- c('a_ACE', 'a_Chao1' , 'a_OTU' , 'a_Shannon' , 'a_Simpson' , 'a_Sobs')

a_ACE <- read.table("input_alpha_div/a_ACE.tsv", header = TRUE, sep = "\t")
a_Chao1 <- read.table("input_alpha_div/a_Chao1.tsv", header = TRUE, sep = "\t")
a_OTU <- read.table("input_alpha_div/a_OTU.tsv", header = TRUE, sep = "\t")
a_Shannon <- read.table("input_alpha_div/a_Shannon.tsv", header = TRUE, sep = "\t")
a_Simpson <- read.table("input_alpha_div/a_Simpson.tsv", header = TRUE, sep = "\t")
a_Sobs <- read.table("input_alpha_div/a_Sobs.tsv", header = TRUE, sep = "\t")

#### merge metadata to wanted data 



#### + bugs ######

P_Actinobacteria <- read.table("input_rela_abd/P_Actinobacteria.tsv", header = TRUE, sep = "\t")
P_Bacteroidetes <- read.table("input_rela_abd/P_Bacteroidetes.tsv", header = TRUE, sep = "\t")
P_Firmicutes <- read.table("input_rela_abd/P_Firmicutes.tsv", header = TRUE, sep = "\t")
P_Proteobacteria <- read.table("input_rela_abd/P_Proteobacteria.tsv", header = TRUE, sep = "\t")
P_Verrucomicrobia <- read.table("input_rela_abd/P_Verrucomicrobia.tsv", header = TRUE, sep = "\t")

G_Bacteroides <- read.table("input_rela_abd/G_Bacteroides.tsv", header = TRUE, sep = "\t")
G_Bifidobacterium <- read.table("input_rela_abd/G_Bifidobacterium.tsv", header = TRUE, sep = "\t")
G_Escherichia <- read.table("input_rela_abd/G_Escherichia.tsv", header = TRUE, sep = "\t")
G_Faecalibacterium <- read.table("input_rela_abd/G_Faecalibacterium.tsv", header = TRUE, sep = "\t")
G_Lactobacillus <- read.table("input_rela_abd/G_Lactobacillus.tsv", header = TRUE, sep = "\t")
G_Roseburia <- read.table("input_rela_abd/G_Roseburia.tsv", header = TRUE, sep = "\t")
G_Ruminococcus <- read.table("input_rela_abd/G_Ruminococcus.tsv", header = TRUE, sep = "\t")




################ META ANALYSIS ###########
###### + SMD #######
#run meta analysis
metacont_ACE <- metacont(i_num, i_mean, i_sd, c_num, c_mean, c_sd, 
         sm = "SMD", method.smd = "Hedges", id,
         data = alpha_ACE_meta)

meta_P_Firmicutes <- metacont(i_num, i_mean, i_sd, c_num, c_mean, c_sd, 
                         sm = "SMD", method.smd = "Hedges", id,
                         data = P_Firmicutes)


######### + ROM ########


##### ++ apha diversity ######
##grouped

gtROM_ACE <- gtROM(a_ACE)
gtROM_Chao1 <- gtROM(a_Chao1)
gtROM_OTU <- gtROM(a_OTU)
gtROM_Shannon <- gtROM(a_Shannon)
gtROM_Simpson <- gtROM(a_Simpson)
gtROM_Sobs <- gtROM(a_Sobs)

gtROM_Actinobacteria <- gtROM(P_Actinobacteria)
gtROM_Bacteroidetes <- gtROM(P_Bacteroidetes)
gtROM_Firmicutes <- gtROM(P_Firmicutes)
gtROM_Proteobacteria <- gtROM(P_Proteobacteria)
gtROM_Verrucomicrobia <- gtROM(P_Verrucomicrobia)

gtROM_Bacteroides <- gtROM(G_Bacteroides)
gtROM_Bifidobacterium <- gtROM(G_Bifidobacterium)
gtROM_Escherichia <- gtROM(G_Escherichia)
gtROM_Faecalibacterium <- gtROM(G_Faecalibacterium)
gtROM_Lactobacillus <- gtROM(G_Lactobacillus)
gtROM_Roseburia <- gtROM(G_Roseburia)
gtROM_Ruminococcus <- gtROM(G_Ruminococcus)

##when subgroup analysis is meaningless

ROM_Sobs <- metacont(i_num, i_mean, i_sd, c_num, c_mean, c_sd, 
                      sm = "ROM", id,  data = a_Sobs)
tROM_Sobs <- table_ROM(ROM_Sobs)



######## DRAW PLOT ##########3


#+ e skyblue3 / c pink3 / ns gray44
#+ na 1 1.5 / stat 16 2 / sum 18 3 / 

check_plot(gtROM_ACE)
ggsave("output_alpha/gtROM_ACE.png", device = "png", width = 4, height = 4, units = "in")
check_plot(gtROM_Chao1)
ggsave("output_alpha/gtROM_Chao1.png", device = "png", width = 4, height = 4, units = "in")
check_plot(gtROM_OTU)
ggsave("output_alpha/gtROM_OTU.png", device = "png", width = 4, height = 4, units = "in")
check_plot(gtROM_Shannon)
ggsave("output_alpha/gtROM_Shannon.png", device = "png", width = 4, height = 4, units = "in")
check_plot(gtROM_Simpson)
ggsave("output_alpha/gtROM_Simpson.png", device = "png", width = 4, height = 4, units = "in")
check_plot(gtROM_Sobs)
ggsave("output_alpha/gtROM_Sobs.png", device = "png", width = 4, height = 4, units = "in")


check_plot(gtROM_Actinobacteria)
ggsave("output_alpha/gtROM_Actinobacteria.png", device = "png", width = 4, height = 4, units = "in")
check_plot(gtROM_Bacteroidetes)
ggsave("output_alpha/gtROM_Bacteroidetes.png", device = "png", width = 4, height = 4, units = "in")

check_plot(gtROM_Firmicutes)
ggplot(data = gtROM_Firmicutes, aes(x=TE, y=studlab, color = signif, shape = stat)) +
  ggtitle(deparse(substitute(gtROM_Firmicutes)))+
  geom_vline(xintercept = 0, linetype="dashed", color="grey50") +
  geom_errorbarh(aes(xmin=lower, xmax=upper), height=0 ) +
  scale_color_manual(values=c("skyblue3", "pink3"))  +
  geom_point(aes(size = stat))+
  scale_shape_manual(values=c(18, 5, 16, 1))+
  scale_size_manual(values=c(2.5, 2, 1, 0.5))+
  scale_y_discrete(limits=rev)+
  theme(panel.border = element_blank(), axis.line = element_line()) +
  theme(axis.text.x=element_text(hjust=1))
ggsave("output_alpha/gtROM_Firmicutes.png", device = "png", width = 4, height = 4, units = "in")

check_plot(gtROM_Proteobacteria)
ggsave("output_alpha/gtROM_Proteobacteria.png", device = "png", width = 4, height = 4, units = "in")
check_plot(gtROM_Verrucomicrobia)
ggsave("output_alpha/gtROM_Verrucomicrobia.png", device = "png", width = 4, height = 4, units = "in")


check_plot(gtROM_Bacteroides)
ggsave("output_alpha/gtROM_Bacteroides.png", device = "png", width = 4, height = 4, units = "in")
check_plot(gtROM_Bifidobacterium)
ggsave("output_alpha/gtROM_Bifidobacterium.png", device = "png", width = 4, height = 4, units = "in")
check_plot(gtROM_Escherichia)
ggsave("output_alpha/gtROM_Escherichia.png", device = "png", width = 4, height = 4, units = "in")
check_plot(gtROM_Faecalibacterium)
ggsave("output_alpha/gtROM_Faecalibacterium.png", device = "png", width = 4, height = 4, units = "in")

check_plot(gtROM_Lactobacillus)
ggplot(data = gtROM_Lactobacillus, aes(x=TE, y=studlab, color = signif, shape = stat)) +
  ggtitle(deparse(substitute(gtROM_Lactobacillus)))+
  geom_vline(xintercept = 0, linetype="dashed", color="grey50") +
  geom_errorbarh(aes(xmin=lower, xmax=upper), height=0 ) +
  scale_color_manual(values=c("skyblue3", "pink3"))  +
  geom_point(aes(size = stat))+
  scale_shape_manual(values=c(18, 5, 16, 1))+
  scale_size_manual(values=c(2.5, 2, 1, 0.5))+
  scale_y_discrete(limits=rev)+
  theme(panel.border = element_blank(), axis.line = element_line()) +
  theme(axis.text.x=element_text(hjust=1))
ggsave("output_alpha/gtROM_Lactobacillus.png", device = "png", width = 4, height = 4, units = "in")

check_plot(gtROM_Roseburia)
ggsave("output_alpha/gtROM_Roseburia.png", device = "png", width = 4, height = 4, units = "in")
check_plot(gtROM_Ruminococcus)
ggsave("output_alpha/gtROM_Ruminococcus.png", device = "png", width = 4, height = 4, units = "in")



check_plot(tROM_ACE)
check_plot(gtROM_ACE)

ggsave("output_alpha/forest_ACE.png", device = "png", width = 4, height = 4, units = "in")


check_plot(tROM_Chao1)
check_plot(tROM_Chao1)

ggsave("output_alpha/forest_Chao1.png", device = "png", width = 3, height = 3, units = "in")



check_plot(tROM_OTU)

ggsave("output_alpha/forest_OTU.png", device = "png", width = 3, height = 3, #5로 나눈 정도
       units = "in")



check_plot(tROM_Shannon)

ggsave("output_alpha/forest_Shannon.png", device = "png", width = 3, height = 3, #5로 나눈 정도
       units = "in")

check_plot(tROM_Simpson)

ggsave("output_alpha/forest_Simpson.png", device = "png", width = 3, height = 3, #5로 나눈 정도
       units = "in")




check_plot(tROM_Sobs)

ggsave("output_alpha/forest_Sobs.png", device = "png", width = 3, height = 3, #5로 나눈 정도
       units = "in")








ggplot(data = tROM_ACE, aes(x=TE, y=studlab, color = signif, shape = stat)) +
  geom_vline(xintercept = 0, linetype="dashed", color="grey50") +
  geom_errorbarh(aes(xmin=lower, xmax=upper), height=0 ) +
  theme_MicrobeR() +
  scale_color_manual(values=c("skyblue3", "pink3", "gray44"))  +
  geom_point(aes(size = stat))+
  scale_shape_manual(values=c(16, 18))+
  scale_size_manual(values=c(1.5, 2))+
  scale_y_discrete(limits=rev)+
  theme(panel.border = element_blank(), axis.line = element_line()) +
  theme(axis.text.x=element_text(angle=45, hjust=1))





ggplot(data = tROM_Sobs, aes(x=TE, y=studlab, color = signif, shape = stat)) +
  geom_vline(xintercept = 0, linetype="dashed", color="grey50") +
  geom_errorbarh(aes(xmin=lower, xmax=upper), height=0 ) +
  theme_MicrobeR() +
  scale_color_manual(values=c( "gray44"))  +
  geom_point(aes(size = stat))+
  scale_shape_manual(values=c(16, 18))+
  scale_size_manual(values=c(1.5, 2))+
  scale_y_discrete(limits=rev)+
  theme(panel.border = element_blank(), axis.line = element_line()) +
  theme(axis.text.x=element_text(angle=45, hjust=1))














##