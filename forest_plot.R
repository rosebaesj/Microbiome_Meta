######### TABLE #######
library(meta)
library(vegan)
library(randomForest)
library(MicrobeR)
library(tidyverse)
library(forestplot)
library(ggpubr)
library(gridExtra)
library(mixmeta)
library(dplyr)
library(writexl)


getwd()
setwd("gitignore")

################ IMPORT DATA ###########

metadata <- read.table("metadata.tsv", header = TRUE, sep = "\t")

meta_div <- read.table("metadata_div.tsv", header = TRUE, sep = "\t")
meta_P <- read.table("metadata_P.tsv", header = TRUE, sep = "\t")
meta_G <- read.table("metadata_G.tsv", header = TRUE, sep = "\t")

list_div <- make_list(meta_div)
list_P <- make_list(meta_P)
list_G <- make_list(meta_G)


colnames(list) <- c("Study", "studlab")
subgroup_list <- c("Overall.(RE)", "Overall.(FE)" ,
                   "Disease", "--GI", "--Obesity", "--Others", 
                   "Treatment.type", "--EA", "--MA", "--MOX",
                   "Treatment.duration", "1", "2","3","4","5~",
                   "Species", "--Human", "--Rat", "--Mouse", "--Rabbit")

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
####### ++ total #####
####### +++ alpha diversity #####
table_total_Chao1 <-total_ROM(a_Chao1,  name = "Chao1")
table_total_Shannon <-total_ROM(a_Shannon,  name = "Shannon")
table_total_Simpson <-total_ROM(a_Simpson,  name = "Simpson")
table_total_ACE <-total_ROM(a_ACE,  name = "ACE")
table_total_OTU <-total_ROM(a_OTU,  name = "OTU")
table_total_Sobs <-total_ROM(a_Sobs,  name = "Sobs")

list_div <- make_list(meta_div)
p_list_div <- plot_list(list_div)


#### +++++in ROM form #######
p_Chao1 <-exp_check_plot_total(table_total_Chao1, list = list_div,name = "Chao1")
p_Shannon <- exp_check_plot_total(table_total_Shannon, list = list_div,name = "Shannon")
p_Simpson <- exp_check_plot_total(table_total_Simpson,list = list_div, name = "Simpson")
p_ACE <- exp_check_plot_total(table_total_ACE,list = list_div, name = "ACE")
p_OTU <- exp_check_plot_total(table_total_OTU, list = list_div,name = "OTU")
p_Sobs <- exp_check_plot_total(table_total_Sobs,list = list_div, name = "Sobs")

p_list <- plot_list(list_div) ###how to make p_list? see functions.R
png("output_forestplot/exp_total_alpha.png", width = 12, height = (29+3)/5, units = "in", res = 300)
grid.arrange(p_list, p_Chao1, p_Shannon, p_Simpson, p_ACE, p_OTU, p_Sobs, 
             ncol=7, 
             widths= c(30,11.667,11.667,11.667,11.667,11.667,11.667)#, heights = c(5,5,5,5,5,5,5)
             )
dev.off()


#### +++++in logROM form #######

p_Chao1 <- check_plot_total(table_total_Chao1, list = list_div,name = "Chao1")
p_Shannon <- check_plot_total(table_total_Shannon, list = list_div,name = "Shannon")
p_Simpson <- check_plot_total(table_total_Simpson,list = list_div, name = "Simpson")
p_ACE <- check_plot_total(table_total_ACE,list = list_div, name = "ACE")
p_OTU <- check_plot_total(table_total_OTU, list = list_div,name = "OTU")
p_Sobs <- check_plot_total(table_total_Sobs,list = list_div, name = "Sobs")

p_list <- plot_list(list_div) ###how to make p_list? see functions.R
png("output_forestplot/total_alpha.png", width = 12, height = (29+3)/5, units = "in", res = 300)
grid.arrange(p_list, p_Chao1, p_Shannon, p_Simpson, p_ACE, p_OTU, p_Sobs, 
             ncol=7, 
             widths= c(30,11.667,11.667,11.667,11.667,11.667,11.667)#, heights = c(5,5,5,5,5,5,5)
)
dev.off()



####### +++ Phylum #####

P_Actinobacteria <- read.table("input_rela_abd/P_Actinobacteria.tsv", header = TRUE, sep = "\t")
P_Bacteroidetes <- read.table("input_rela_abd/P_Bacteroidetes.tsv", header = TRUE, sep = "\t")
P_Firmicutes <- read.table("input_rela_abd/P_Firmicutes.tsv", header = TRUE, sep = "\t")
P_Proteobacteria <- read.table("input_rela_abd/P_Proteobacteria.tsv", header = TRUE, sep = "\t")
P_Verrucomicrobia <- read.table("input_rela_abd/P_Verrucomicrobia.tsv", header = TRUE, sep = "\t")


table_total_P_Actinobacteria <-total_ROM(P_Actinobacteria,  name = "P_Actinobacteria")
table_total_P_Bacteroidetes <-total_ROM(P_Bacteroidetes,  name = "P_Bacteroidetes")
table_total_P_Firmicutes <-total_ROM(P_Firmicutes,  name = "P_Firmicutes")
table_total_P_Proteobacteria <-total_ROM(P_Proteobacteria,  name = "P_Proteobacteria")
table_total_P_Verrucomicrobia <-total_ROM(P_Verrucomicrobia,  name = "P_Verrucomicrobia")

list_P <- make_list(meta_P)
p_list_P <- plot_list(list_P) ###how to make p_list? see functions.R



#### +++++in ROM form #######

pl_Actinobacteria <-exp_check_plot_total(table_total_P_Actinobacteria, list = list_P, name = "P_Actinobacteria")
pl_Bacteroidetes <- exp_check_plot_total(table_total_P_Bacteroidetes,list = list_P,  name = "P_Bacteroidetes")
pl_Firmicutes <- exp_check_plot_total(table_total_P_Firmicutes, list = list_P, name = "P_Firmicutes")
pl_Proteobacteria <- exp_check_plot_total(table_total_P_Proteobacteria, list = list_P, name = "P_Proteobacteria")
pl_Verrucomicrobia <- exp_check_plot_total(table_total_P_Verrucomicrobia,list = list_P,  name = "P_Verrucomicrobia")

png("output_forestplot/exp_total_phylum.png", width = 12, height = (17+3)/5, units = "in", res = 300)
grid.arrange(p_list_P, pl_Actinobacteria, pl_Bacteroidetes, pl_Firmicutes, pl_Proteobacteria, pl_Verrucomicrobia, 
                   ncol=6, 
                   widths= c(30,14,14,14,14,14)#, heights = c(5,5,5,5,5,5,5)
)
dev.off()


#### +++++in logROM form #######

pl_Actinobacteria <-check_plot_total(table_total_P_Actinobacteria, list = list_P, name = "P_Actinobacteria")
pl_Bacteroidetes <- check_plot_total(table_total_P_Bacteroidetes,list = list_P,  name = "P_Bacteroidetes")
pl_Firmicutes <- check_plot_total(table_total_P_Firmicutes, list = list_P, name = "P_Firmicutes")
pl_Proteobacteria <- check_plot_total(table_total_P_Proteobacteria, list = list_P, name = "P_Proteobacteria")
pl_Verrucomicrobia <- check_plot_total(table_total_P_Verrucomicrobia,list = list_P,  name = "P_Verrucomicrobia")

png("output_forestplot/total_phylum.png", width = 12, height = (17+3)/5, units = "in", res = 300)
grid.arrange(p_list_P, pl_Actinobacteria, pl_Bacteroidetes, pl_Firmicutes, pl_Proteobacteria, pl_Verrucomicrobia, 
             ncol=6, 
             widths= c(30,14,14,14,14,14)#, heights = c(5,5,5,5,5,5,5)
)
dev.off()




####### +++ Genus #####

G_Bacteroides <- read.table("input_rela_abd/G_Bacteroides.tsv", header = TRUE, sep = "\t")
G_Bifidobacterium <- read.table("input_rela_abd/G_Bifidobacterium.tsv", header = TRUE, sep = "\t")
G_Escherichia <- read.table("input_rela_abd/G_Escherichia.tsv", header = TRUE, sep = "\t")
G_Faecalibacterium <- read.table("input_rela_abd/G_Faecalibacterium.tsv", header = TRUE, sep = "\t")
G_Lactobacillus <- read.table("input_rela_abd/G_Lactobacillus.tsv", header = TRUE, sep = "\t")
G_Roseburia <- read.table("input_rela_abd/G_Roseburia.tsv", header = TRUE, sep = "\t")
G_Ruminococcus <- read.table("input_rela_abd/G_Ruminococcus.tsv", header = TRUE, sep = "\t")


table_total_G_Bacteroides <-total_ROM(G_Bacteroides,  name = "G_Bacteroides")
table_total_G_Bifidobacterium <-total_ROM(G_Bifidobacterium,  name = "G_Bifidobacterium")
table_total_G_Escherichia <-total_ROM(G_Escherichia,  name = "G_Escherichia")
table_total_G_Faecalibacterium <-total_ROM(G_Faecalibacterium,  name = "G_Faecalibacterium")
table_total_G_Lactobacillus <-total_ROM(G_Lactobacillus,  name = "G_Lactobacillus")
table_total_G_Roseburia <-total_ROM(G_Roseburia,  name = "G_Roseburia")
table_total_G_Ruminococcus <-total_ROM(G_Ruminococcus,  name = "G_Ruminococcus")

list_G <- make_list(meta_G)
p_list_G <- plot_list(list_G) ###how to make p_list? see functions.R


#### +++++in ROM form #######

pG_Bacteroides <-exp_check_plot_total(table_total_G_Bacteroides, list = list_G, name = "G_Bacteroides")
pG_Bifidobacterium <- exp_check_plot_total(table_total_G_Bifidobacterium,list = list_G,  name = "G_Bifidobacterium")
pG_Escherichia <- exp_check_plot_total(table_total_G_Escherichia, list = list_G, name = "G_Escherichia")
pG_Faecalibacterium <- exp_check_plot_total(table_total_G_Faecalibacterium, list = list_G, name = "G_Faecalibacterium")
pG_Lactobacillus <- exp_check_plot_total(table_total_G_Lactobacillus,list = list_G,  name = "G_Lactobacillus")
pG_Roseburia <- exp_check_plot_total(table_total_G_Roseburia,list = list_G,  name = "G_Roseburia")
pG_Ruminococcus <- exp_check_plot_total(table_total_G_Ruminococcus, list = list_G, name = "G_Ruminococcus")

png("output_forestplot/exp_total_genus.png", width = 12, height = (27+3)/5, units = "in", res = 300)
grid.arrange(p_list_G, pG_Bacteroides, pG_Bifidobacterium, pG_Escherichia, pG_Faecalibacterium, 
                   pG_Lactobacillus,  pG_Roseburia, pG_Ruminococcus,
                   ncol=8, 
                   widths= c(30,10,10,10,10,10,10,10)#, heights = c(5,5,5,5,5,5,5)
)
dev.off()

#### +++++in logROM form #######

pG_Bacteroides <-check_plot_total(table_total_G_Bacteroides, list = list_G, name = "G_Bacteroides")
pG_Bifidobacterium <- check_plot_total(table_total_G_Bifidobacterium,list = list_G,  name = "G_Bifidobacterium")
pG_Escherichia <- check_plot_total(table_total_G_Escherichia, list = list_G, name = "G_Escherichia")
pG_Faecalibacterium <- check_plot_total(table_total_G_Faecalibacterium, list = list_G, name = "G_Faecalibacterium")
pG_Lactobacillus <- check_plot_total(table_total_G_Lactobacillus,list = list_G,  name = "G_Lactobacillus")
pG_Roseburia <- check_plot_total(table_total_G_Roseburia,list = list_G,  name = "G_Roseburia")
pG_Ruminococcus <- check_plot_total(table_total_G_Ruminococcus, list = list_G, name = "G_Ruminococcus")

png("output_forestplot/total_genus.png", width = 12, height = (27+3)/5, units = "in", res = 300)
grid.arrange(p_list_G, pG_Bacteroides, pG_Bifidobacterium, pG_Escherichia, pG_Faecalibacterium, 
             pG_Lactobacillus,  pG_Roseburia, pG_Ruminococcus,
             ncol=8, 
             widths= c(30,10,10,10,10,10,10,10)#, heights = c(5,5,5,5,5,5,5)
)
dev.off()


####### ++ subgroup #####

##### +++ apha diversity ######


table_sub_Chao1 <- subgroup_ROM(a_Chao1, metadata)
table_sub_Shannon <- subgroup_ROM(a_Shannon, metadata)
table_sub_Simpson <- subgroup_ROM(a_Simpson, metadata)
table_sub_ACE <- subgroup_ROM(a_ACE, metadata)
table_sub_OTU <- subgroup_ROM(a_OTU, metadata)
table_sub_Sobs <- subgroup_ROM(a_Sobs, metadata)


exp_draw_sub(table_sub_Chao1, n=19,
         colorvector = c("white", "white", 
                                    "grey90","grey90","grey90","grey90",
                                    "white", "white", "white",
                                    "grey90","grey90","grey90","grey90",
                                    "white","white","white", "white", "white", "white"),
           filepath = "output_forestplot/exp_sub_Chao1.png")
exp_draw_sub(table_sub_Shannon, name = "Shannon",n=20,
         colorvector = c("white", "white", 
                                            "grey90","grey90","grey90","grey90","grey90",
                                            "white", "white", "white",
                                            "grey90","grey90","grey90","grey90",
                                            "white","white","white", "white", "white", "white"),
         filepath = "output_forestplot/exp_sub_Shannon.png")
exp_draw_sub(table_sub_Simpson, n=21,
         colorvector = c("white", "white", 
                         "grey90","grey90","grey90","grey90",
                         "white", "white", "white","white",
                         "grey90","grey90","grey90","grey90","grey90",
                         "white","white","white", "white", "white", "white"),
         filepath = "output_forestplot/exp_sub_Simpson.png")
exp_draw_sub(table_sub_ACE,n=21,
         colorvector = c("white", "white", 
                         "grey90","grey90","grey90","grey90",
                         "white", "white", "white","white",
                         "grey90","grey90","grey90","grey90","grey90",
                         "white","white","white", "white", "white", "white"),
         filepath = "output_forestplot/exp_sub_ACE.png")
exp_draw_sub(table_sub_OTU, n=16,
         colorvector = c("white", "white", 
                         "grey90","grey90","grey90",
                         "white", "white", "white",
                         "grey90","grey90","grey90",
                         "white","white","white", "white", "white"),
         filepath = "output_forestplot/exp_sub_OTU.png")
exp_draw_sub(table_sub_Sobs, n=14,
         colorvector = c("white", "white", 
                         "grey90","grey90","grey90",
                         "white", "white",  "white",
                         "grey90","grey90","grey90",
                         "white","white","white"),
         filepath = "output_forestplot/exp_sub_Sobs.png")


##### +++ phylum ######



table_sub_P_Act <- subgroup_ROM(P_Actinobacteria, metadata)
table_sub_P_Bac <- subgroup_ROM(P_Bacteroidetes, metadata)
table_sub_P_Fir <- subgroup_ROM(P_Firmicutes, metadata)
table_sub_P_Pro <- subgroup_ROM(P_Proteobacteria, metadata)
table_sub_P_Ver <- subgroup_ROM(P_Verrucomicrobia, metadata)


logexp_draw_sub(table_sub_P_Act, n=17,
         colorvector = c("grey90","grey90",
                         "white", "white", "white", "white", 
                         "grey90","grey90","grey90",
                         "white", "white", "white",
                         "grey90","grey90","grey90","grey90","grey90"),
         filepath = "output_forestplot/logexp_sub_P_Actinobacteria.png")

base19 <-  c("grey90","grey90",
           "white", "white", "white", "white", 
           "grey90","grey90","grey90","grey90",
           "white", "white", "white", "white",
           "grey90","grey90","grey90","grey90","grey90")

logexp_draw_sub(table_sub_P_Bac, name = "Shannon",n=19,
         colorvector = base19,
         filepath = "output_forestplot/logexp_sub_P_Bacteroidetes.png")
logexp_draw_sub(table_sub_P_Fir, n=19,
         colorvector = base19,
         filepath = "output_forestplot/logexp_sub_P_Firmicutes.png")
logexp_draw_sub(table_sub_P_Pro,n=19,
         colorvector = base19,
         filepath = "output_forestplot/logexp_sub_P_Proteobacteria.png")
logexp_draw_sub(table_sub_P_Ver, n=18,
         colorvector = c("grey90","grey90",
                         "white", "white", "white", "white", 
                         "grey90","grey90","grey90","grey90",
                         "white", "white", "white", "white",
                         "grey90","grey90","grey90","grey90"),
         filepath = "output_forestplot/logexp_sub_P_Verrucomicrobia.png")

##### +++ genus ######


#####ERROR HERE#####
table_sub_G_Bacteroides <- subgroup_ROM(G_Bacteroides, metadata)
table_sub_G_Bifidobacterium <- subgroup_ROM(G_Bifidobacterium, metadata)
table_sub_G_Escherichia <- subgroup_ROM(G_Escherichia, metadata)
table_sub_G_Faecalibacterium <- subgroup_ROM(G_Faecalibacterium, metadata)
table_sub_G_Lactobacillus <- subgroup_ROM(G_Lactobacillus, metadata)
table_sub_G_Roseburia <- subgroup_ROM(G_Roseburia, metadata)
table_sub_G_Ruminococcus <- subgroup_ROM(G_Ruminococcus, metadata)


exp_draw_sub(table_sub_G_Bacteroides, n=20,
         colorvector = c("grey90","grey90",
                         "white", "white", "white", "white", 
                         "grey90","grey90","grey90","grey90",
                         "white", "white", "white", "white",
                         "grey90","grey90","grey90","grey90","grey90","grey90"),
         filepath = "output_forestplot/exp_sub_G_Bacteroides.png")

exp_draw_sub(table_sub_G_Bifidobacterium, name = "Shannon",n=19,
         colorvector = base19,
         filepath = "output_forestplot/exp_sub_G_Bifidobacterium.png")

exp_draw_sub(table_sub_G_Escherichia, n=14,
         colorvector = c("grey90","grey90",
                         "white", "white", "white",
                         "grey90","grey90","grey90",
                         "white", "white", "white",
                         "grey90","grey90","grey90"),
         filepath = "output_forestplot/exp_sub_G_Escherichia.png")
exp_draw_sub(table_sub_G_Faecalibacterium,n=15,
         colorvector =  c("grey90","grey90",
                          "white", "white", "white","white",
                          "grey90","grey90","grey90",
                          "white", "white", "white",
                          "grey90","grey90"),
         filepath = "output_forestplot/exp_sub_G_Faecalibacterium.png")

exp_draw_sub(table_sub_G_Lactobacillus, n=19,
         colorvector = base19,
         filepath = "output_forestplot/exp_sub_G_Lactobacillus.png")

exp_draw_sub(table_sub_G_Roseburia, n=19,
         colorvector = base19,
         filepath = "output_forestplot/exp_sub_G_Roseburia.png")

exp_draw_sub(table_sub_G_Ruminococcus, n=16,
         colorvector = c("grey90","grey90",
                         "white", "white", "white","white",
                         "grey90","grey90","grey90","grey90",
                         "white", "white", "white",
                         "grey90","grey90","grey90"),
         filepath = "output_forestplot/exp_sub_G_Ruminococcus.png")





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



## subgroup analysis version 2 GI , obesity
ggtROM_ACE <- ggtROM(a_ACE)
ggtROM_Chao1 <- ggtROM(a_Chao1)
ggtROM_OTU <- ggtROM(a_OTU)
ggtROM_Shannon <- ggtROM(a_Shannon)
ggtROM_Simpson <- ggtROM(a_Simpson)
ggtROM_Sobs <- ggtROM(a_Sobs)

ggtROM_Actinobacteria <- ggtROM(P_Actinobacteria)
ggtROM_Bacteroidetes <- ggtROM(P_Bacteroidetes)
ggtROM_Firmicutes <- ggtROM(P_Firmicutes)
ggtROM_Proteobacteria <- ggtROM(P_Proteobacteria)
ggtROM_Verrucomicrobia <- ggtROM(P_Verrucomicrobia)

ggtROM_Bacteroides <- ggtROM(G_Bacteroides)
ggtROM_Bifidobacterium <- ggtROM(G_Bifidobacterium)
ggtROM_Escherichia <- ggtROM(G_Escherichia)
ggtROM_Faecalibacterium <- ggtROM(G_Faecalibacterium)
ggtROM_Lactobacillus <- ggtROM(G_Lactobacillus)
ggtROM_Roseburia <- ggtROM(G_Roseburia)
ggtROM_Ruminococcus <- ggtROM(G_Ruminococcus)









######## DRAW PLOT ##########3


#+ e skyblue3 / c pink3 / ns gray44
#+ na 1 1.5 / stat 16 2 / sum 18 3 / 

p_ACE<-check_plot(gtROM_ACE)
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


##ggR
check_plot(ggtROM_ACE)
ggsave("output_alpha/ggtROM_ACE.png", device = "png", width = 4, height = 4, units = "in")
check_plot(ggtROM_Chao1)
ggsave("output_alpha/ggtROM_Chao1.png", device = "png", width = 4, height = 4, units = "in")
check_plot(ggtROM_OTU)
ggsave("output_alpha/ggtROM_OTU.png", device = "png", width = 4, height = 4, units = "in")
check_plot(ggtROM_Shannon)
ggsave("output_alpha/ggtROM_Shannon.png", device = "png", width = 4, height = 4, units = "in")
check_plot(ggtROM_Simpson)
ggsave("output_alpha/ggtROM_Simpson.png", device = "png", width = 4, height = 4, units = "in")
check_plot(ggtROM_Sobs)
ggsave("output_alpha/ggtROM_Sobs.png", device = "png", width = 4, height = 4, units = "in")









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



################ RISK OF BIAS ###########
library(robvis)
library("dmetar")
#https://bookdown.org/MathiasHarrer/Doing_Meta_Analysis_in_R/risk-of-bias-plots.html
Bias_levels <- c("Sequence generation",
                 "Baseline characteristics",
                 "Allocation concealment",
                 "Random housing",
                 "Blinding (intervention)",
                 "Random outcome assessment",
                 "Blinding (outcome)",
                 "Incomplete outcome data",
                 "Selective outcome reporting",
                 "Other sources of bias")
Study_levels <- rownames(rob_all)
Risk_levels <- c("Low", "Unknown", "High", "Not applicable")


rob <- read.csv("input_ROB/ROB.csv", header = TRUE)
rownames(rob) <- rob[,1]
rob <- rob[,-1]

robsum_coord <- data.frame()
for (r in 1:nrow(rob)) {
  for(c in 1:ncol(rob)) {
    coord <- c(colnames(rob)[c], rownames(rob)[r], rob[r,c])
    robsum_coord <- rbind(robsum_coord, coord)
  }
}
colnames(robsum_coord) <- c("Bias", "Risk", "Number")

robsum_coord$Number <- as.numeric(robsum_coord$Number)
Bias_levels <- colnames(rob)
Risk_levels <- c("Low", "Unclear", "High", "Not.applicable")


robsum_coord$Bias <- factor(robsum_coord$Bias, levels = rev(Bias_levels), ordered = T)
robsum_coord$Risk <- factor(robsum_coord$Risk, levels = rev(Risk_levels), ordered = T)

n=23

ggplot(data=robsum_coord, aes(fill=Risk, x=Bias, y=Number))+
  geom_bar(position="fill", stat="identity", colour = 'black', width=0.7)+
  coord_flip()+
  scale_fill_manual(values = c("Low" = 'green3', "Unclear" = 'yellow2', 
                                       "High" = 'red2', "Not.applicable"='darkgrey'))+
  theme_classic()+
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), 
                     expand = c(0,0), labels = scales::percent)+
  # total 23
  theme(axis.ticks.x = element_line(colour = "black"),
        axis.title = element_blank(),
        legend.position = "none",
        plot.margin = unit(c(0.1, 0.3, 0.1, 0.1), "in")
        )
ggsave("output_ROB/ROB_summary.png", device = "png", width = 5, height = 2.2, units = "in")

 
# rob_table <- read.csv("input_ROB/ROB_table.csv", header = F)
# row.names(rob_table) <- rob_table[,1]
# rob_table <- rob_table[,-1]
# 
# rob_table[,1] <- as.numeric(rob_table[,1])
# rob_table[,2] <- as.numeric(rob_table[,2])
# rob_table[,3] <- as.numeric(rob_table[,3])
# rob_table[,4] <- as.numeric(rob_table[,4])
# rob_table[,5] <- as.numeric(rob_table[,5])
# rob_table[,6] <- as.numeric(rob_table[,6])
# rob_table[,7] <- as.numeric(rob_table[,7])
# rob_table[,8] <- as.numeric(rob_table[,8])
# rob_table[,9] <- as.numeric(rob_table[,9])
# rob_table[,10] <- as.numeric(rob_table[,10])
# 
# Bias<-rob$Bias
# tabl <- table(rob_table, Bias)
# 
# 
# barplot(rob_table,
#         horiz = TRUE)
##


rob_all <- data.frame(read.csv("input_ROB/ROB_all.csv", header = T))
rownames(rob_all) <- rob_all[,1]
rob_all <- rob_all[,-1]

rob_coord <- data.frame()
sign <- NULL
for (r in 1:nrow(rob_all)) {
  for(c in 1:ncol(rob_all)) {
    if (rob_all[r,c]=="Low") {sig <- c("+")
    } else if(rob_all[r,c]=="Unclear") {sig <- c("?")
    } else if(rob_all[r,c]=="High") {sig <- c("-")
    } else {sign <- c("NA")}
    
    coord <- c(colnames(rob_all)[c], rownames(rob_all)[r], rob_all[r,c], sig)
    rob_coord <- rbind(rob_coord, coord)
  }
}

colnames(rob_coord) <- c("Bias", "Study", "Risk", "Sign")

Bias_levels <- colnames(rob_all)
Study_levels <- rownames(rob_all)
Risk_levels <- c("Low", "Unclear", "High", "Not.applicable")


rob_coord$Bias <- factor(rob_coord$Bias, levels = Bias_levels, ordered = T)
rob_coord$Study <- factor(rob_coord$Study, levels = rev(Study_levels), ordered = T)
rob_coord$Risk <- factor(rob_coord$Risk, levels = rev(Risk_levels), ordered = T)



ggplot(data=rob_coord, aes(fill=Risk, x=Bias, y=Study, label = Sign))+
  geom_text()+
  geom_vline(xintercept = c(0.45, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5, 10.55), 
             linetype="solid", color="black") +
  geom_hline(yintercept = c(0.45, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 
                            10.5, 11.5, 12.5, 13.5, 14.5, 15.5, 16.5, 17.5, 18.5, 19.5,
                            20.5, 21.5, 22.5, 23.5, 23.55), 
             linetype="solid", color="black")+
  geom_point(shape = 21, size = 7)+
  scale_fill_manual(values = c("Low" = 'green3', "Unclear" = 'yellow2', 
                                           "High" = 'red2', "Not.applicable"='darkgrey'))+
  geom_text()+
  theme_classic()+
  #scale_y_continuous(expand = c(0,0))+
  scale_x_discrete(position = "top", expand = c(0.057,0.057))+
  scale_y_discrete( expand = c(0.0255,0.0255))+
  theme(axis.text.x = element_text(angle=45, hjust=0),
        axis.title = element_blank(),
        axis.line = element_blank(),
        legend.title = element_blank(),
        legend.position = c(1.3, 0.12),
        plot.margin = unit(c(0.1, 1.55, 0.1, 0.1), "in")
  )

ggsave("output_ROB/ROB.png", device = "png", width = 5, height = 7.2, units = "in")








####### FOREST PLOT #########
data <- ggtROM_Chao1[,c(1,2,4,5)]

colnames(data)
data$TE<-round(data$TE, 2)
data$lower<-round(data$lower, 2)
data$upper<-round(data$upper, 2)

png(filename = "output_forestplot/Forestplot.png", width=10, height=5, res = 300, units = "in")
data %>%
forestplot(
  lableltext = data,
  graph.pos=2,
  mean = c(data$TE),
  lower = c(data$lower),
  upper = c(data$upper),
  title = "ggtROM_Chao1",
  xlab = "                      <----control---- --intervention-->",
  hrxl_lines = list("3" = gpar(lwd=1, col="#99999922"), 
                    "7" = gpar(lwd=60, lineend="butt", columns=c(2:4), col="#99999922"),
                    "15" = gpar(lwd=60, lineend="butt", columns=c(2:4), col="#99999922")),
  txt_gp=fpTxtGp(label=gpar(cex=1.25),
                 ticks=gpar(cex=1.1),
                 xlab=gpar(cex = 1.2),
                 title=gpar(cex = 1.2)),
  col=fpColors(box="black", lines="black", zero = "gray50"),
  zero=0, cex=0.9, lineheight = "auto", boxsize=0.5, colgap=unit(6,"mm"),
  lwd.ci=2, ci.vertices=TRUE, ci.vertices.height = 0.2)
dev.off()
ggsave("output_ROB/ggforest_Chao1.png", device = "png", width = 6, height = 2, units = "in")


# data$colour <- rep(c("white", "gray95"), 10)
# ggplot(data, aes(x = rr, y = labels, xmin = rrlow, xmax = rrhigh)) +
#   geom_hline(aes(yintercept = labels, colour = colour), size = 7) + 
#   geom_pointrange(shape = 22, fill = "black") +
#   geom_vline(xintercept = 1, linetype = 3) +
#   xlab("Variable") +
#   ylab("Adjusted Relative Risk with 95% Confidence Interval") +
#   theme_classic() +
#   scale_colour_identity() +
#   scale_y_discrete(limits = rev(data$studlab)) +
#   scale_x_log10(limits = c(0.25, 4), 
#                 breaks = c(0.25, 0.5, 1, 2, 4), 
#                 labels = c("0.25", "0.5", "1", "2", "4"), expand = c(0,0)) +
#   theme(axis.text.y = element_blank(), axis.title.y = element_blank())




t_ACE <- data.frame()
t_ACE$TE <- round(gtROM_ACE$TE, 2)

ggplot(data=gtROM_ACE, aes(y=studlab))+
  #geom_hline(aes(yintercept = labels, colour = colour), size = 7) +
  geom_text(aes(x = 0, label = TE), hjust = 0) +
  geom_text(aes(x = 5, label = lower)) +
  geom_text(aes(x = 7, label = upper), hjust = 1) +
  #scale_colour_identity() +
  #theme_void() + 
  theme(plot.margin = margin(5, 0, 35, 0))

grid.arrange(t_ACE, p_ACE, ncol=2)


ACE <- left_join(list, gtROM_ACE, by = "studlab")



ACE$colour <- rep(c("white", "gray95"), 14)
ggplot(ACE, aes(x = rr, y = labels, xmin = rrlow, xmax = rrhigh)) +
  geom_hline(aes(yintercept = labels, colour = colour), size = 7) + 
  geom_pointrange(shape = 22, fill = "black") +
  geom_vline(xintercept = 1, linetype = 3) +
  xlab("Variable") +
  ylab("Adjusted Relative Risk with 95% Confidence Interval") +
  theme_classic() +
  scale_colour_identity() +
  scale_y_discrete(limits = rev(data$studlab)) +
  scale_x_log10(limits = c(0.25, 4), 
                breaks = c(0.25, 0.5, 1, 2, 4), 
                labels = c("0.25", "0.5", "1", "2", "4"), expand = c(0,0)) +
  theme(axis.text.y = element_blank(), axis.title.y = element_blank())



##




########## multivariate regression #######


Chao1 <- left_join(a_ACE, metadata)
mv_Chao1 <- merge(metadata, table_total_Chao1, by="id")
mv_Chao1 <- mv_Chao1 %>% filter(!is.na(seTE))

model <- mixmeta(TE ~ disease_group, S=seTE, data = mv_Chao1)
coef.mixmeta(model)
plot(fitted(model))






######## funnel plot #####
metait <-  metacont(i_num, i_mean, i_sd, c_num, c_mean, c_sd, 
                                sm = "ROM", id,  data = a_ACE)

metait$TE.random

tableit <-  data.frame(metacont(i_num, i_mean, i_sd, c_num, c_mean, c_sd, 
                       sm = "ROM", id,  data = a_ACE))

# logROM, SE of logROM

ggplot(data = tableit) + 
 geom_vline(xintercept = 0, linetype=1, color = "grey70")+
  geom_vline(xintercept = metait$TE.random, linetype="dashed", color = "grey70")+
  geom_abline(slope = 1/1.96, intercept = -(metait$TE.random)/1.96, linetype="dashed", color = "grey70")+
  geom_abline(slope = -1/1.96, intercept = (metait$TE.random)/1.96, linetype="dashed", color = "grey70")+
  geom_point( aes(x = TE, y = seTE))+
  scale_y_reverse("SE of log(ROM)", expand=c(0.1,0))	+
  scale_x_continuous("log(ROM)", expand=c(0.1,0))+
  theme_classic()
ggsave("output_funnel/ACE_SE.png", device = "png", width = 3, height = 3, units = "in")
summary(lm(data=tableit, formula = seTE~TE))
  #MD logROM,inverse n

ggplot(data = tableit)+
  geom_vline(xintercept = 0, linetype=1, color = "grey70")+
  geom_point(aes(x = TE, y = 1/sqrt(n.e)))+
  scale_y_reverse("inverse sqare root of n", expand=c(0.1,0))	+
  theme_classic()
ggsave("output_funnel/ACE_n.png", device = "png", width = 3, height = 3, units = "in")


metait <-  metacont(i_num, i_mean, i_sd, c_num, c_mean, c_sd, 
                    sm = "ROM", id,  data = G_Lactobacillus)
metait$TE.random
tableit <-  data.frame(metacont(i_num, i_mean, i_sd, c_num, c_mean, c_sd, 
                                sm = "ROM", id,  data = G_Lactobacillus))

# logROM, SE of logROM

ggplot(data = tableit) + 
  geom_vline(xintercept = 0, linetype=1, color = "grey70")+
  geom_vline(xintercept = metait$TE.random, linetype="dashed", color = "grey70")+
  geom_abline(slope = 1/1.96, intercept = -(metait$TE.random)/1.96, linetype="dashed", color = "grey70")+
  geom_abline(slope = -1/1.96, intercept = (metait$TE.random)/1.96, linetype="dashed", color = "grey70")+
  geom_point( aes(x = TE, y = seTE))+
  scale_color_manual(values=c('#80c269', '#eb6877','#f19149',  '#448aca', '#ae5da1')) +  
  scale_y_reverse("SE of log(ROM)", expand=c(0.05,0))	+
  scale_x_continuous("log(ROM)", expand=c(0.05,0))+
  theme_classic()
summary(lm(data=tableit, formula = seTE~TE))
ggsave("output_funnel/Lacto_SE.png", device = "png", width = 3, height = 3, units = "in")

ggplot(data = tableit)+
  geom_vline(xintercept = 0, linetype=1, color = "grey70")+
  geom_point(aes(x = TE, y = 1/sqrt(n.e)))+
  scale_y_reverse("inverse sqare root of n", expand=c(0.5,0))	+
  theme_classic()
summary(lm(data=tableit, formula = (1/sqrt(n.e))~TE))
ggsave("output_funnel/Lacto_n.png", device = "png", width = 3, height = 3, units = "in")


##### Leave one out

Wang_2012 <- metacont(i_num, i_mean, i_sd, c_num, c_mean, c_sd, 
         sm = "ROM", id,  data = G_Lactobacillus %>% filter(Study!="Wang(2012)"))
Xu_2013 <- metacont(i_num, i_mean, i_sd, c_num, c_mean, c_sd, 
                     sm = "ROM", id,  data = G_Lactobacillus %>% filter(Study!="Xu(2013)"))
Bao_2019 <- metacont(i_num, i_mean, i_sd, c_num, c_mean, c_sd, 
                    sm = "ROM", id,  data = G_Lactobacillus %>% filter(Study!="Bao(2019)"))
Huang_2019 <- metacont(i_num, i_mean, i_sd, c_num, c_mean, c_sd, 
                       sm = "ROM", id,  data = G_Lactobacillus %>% filter(Study!="Huang(2019)"))
Li_2020 <- metacont(i_num, i_mean, i_sd, c_num, c_mean, c_sd, 
                    sm = "ROM", id,  data = G_Lactobacillus %>% filter(Study!="Li(2020)"))
Xie_2020 <- metacont(i_num, i_mean, i_sd, c_num, c_mean, c_sd, 
                    sm = "ROM", id,  data = G_Lactobacillus %>% filter(Study!="Xie(2020)"))
Xu_2020 <- metacont(i_num, i_mean, i_sd, c_num, c_mean, c_sd, 
                    sm = "ROM", id,  data = G_Lactobacillus %>% filter(Study!="Xu(2020)"))


string <- c("Q","df.Q", "pval.Q","H", "I2", "tau2", "k",
            "TE.random", "lower.random", "upper.random", "pval.random")
loo <- data.frame()
loo <- data.frame(rbind(Wang_2012[string], Xu_2013[string], Bao_2019[string], Huang_2019[string],
             Li_2020[string],  Xie_2020[string], Xu_2020[string]))

loo <- loo %>% rowwise() %>%
  mutate(ROM = exp(TE.random),
         lower = exp(lower.random),
         upper = exp(upper.random), 
         I2 = I2)



write_csv (loo, "output_funnel/leave_one_out.csv")
