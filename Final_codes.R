##### FINAL CODES ######

library(meta)
library(vegan)
library(randomForest)
library(MicrobeR)
library(tidyverse)
library(forestplot)
library(ggpubr)
library(gridExtra)
library(mixmeta)

getwd()
setwd("gitignore")

metadata <- read.table("metadata.tsv", header = TRUE, sep = "\t")


#### 1. alpha diversity ####

#### + import data ####
a_ACE <- read.table("input_alpha_div/a_ACE.tsv", header = TRUE, sep = "\t")
a_Chao1 <- read.table("input_alpha_div/a_Chao1.tsv", header = TRUE, sep = "\t")
a_Shannon <- read.table("input_alpha_div/a_Shannon.tsv", header = TRUE, sep = "\t")
a_Simpson <- read.table("input_alpha_div/a_Simpson.tsv", header = TRUE, sep = "\t")
a_OTU <- read.table("input_alpha_div/a_OTU.tsv", header = TRUE, sep = "\t")

####### + overall analysis ######
####### ++ make ROM tables #####
table_total_ACE <-total_ROM(a_ACE,  name = "ACE") 
table_total_Chao1 <-total_ROM(a_Chao1,  name = "Chao1")
table_total_Shannon <-total_ROM(a_Shannon,  name = "Shannon")
table_total_Simpson <-total_ROM(a_Simpson,  name = "Simpson")
table_total_OTU <-total_ROM(a_OTU,  name = "OTU")

###### ++ import metadata #####
meta_div <- read.table("metadata_div.tsv", header = TRUE, sep = "\t")
list_div <- make_list(meta_div)
p_list_div <- plot_list(list_div)


#### ++ plot in logROM form #######
p_ACE <- check_plot_total(table_total_ACE,list = list_div, name = "ACE", overall = T)
p_Chao1 <- check_plot_total(table_total_Chao1, list = list_div,name = "Chao1", overall = T)
p_Shannon <- check_plot_total(table_total_Shannon, list = list_div,name = "Shannon", overall = T)
p_Simpson <- check_plot_total(table_total_Simpson,list = list_div, name = "Simpson", overall = T)
p_OTU <- check_plot_total(table_total_OTU, list = list_div,name = "OTU", overall = T)

p_list <- plot_list(list_div) ###how to make p_list? see functions.R
png("output_forestplot/total_alpha.png", width = 10, height = (nrow(list_div)+3)/5, units = "in", res = 300)
grid.arrange(p_list, p_ACE, p_Chao1, p_Shannon, p_Simpson, p_OTU, 
             ncol=6, 
             widths= c(25,14,14,14,14,14)#, heights = c(5,5,5,5,5,5,5)
)
dev.off()


####### + subgroup analysis ######
####### ++ make ROM tables #####
table_sub_ACE <- subgroup_ROM(a_ACE, metadata)
table_sub_Chao1 <- subgroup_ROM(a_Chao1, metadata)

table_sub_RE_ACE <- subgroup_RE_ROM(a_ACE, metadata)
table_sub_RE_Chao1 <- subgroup_RE_ROM(a_Chao1, metadata)


#### ++ plot in logROM form #######
logexp_draw_sub(table_sub_ACE,n=nrow(table_sub_ACE),name="ACE", 
             colorvector = c("grey90", "grey90", 
                             "white", "white", "white","white",
                             "grey90","grey90","grey90",
                             "white", "white", "white","white"
                             ),
             filepath = "output_forestplot/logexp_sub_ACE.png")

logexp_draw_sub(table_sub_Chao1, n=nrow(table_sub_Chao1),name = "Chao1",
             colorvector = c("grey90","grey90",
                             "white", "white", "white",
                             "grey90","grey90","grey90",
                             "white", "white", "white", "white"),
             filepath = "output_forestplot/logexp_sub_Chao1.png")

logexp_draw_sub(table_sub_RE_ACE,n=nrow(table_sub_RE_ACE),name="ACE", 
                colorvector = c("grey90",
                                "white", "white", "white","white",
                                "grey90","grey90","grey90",
                                "white", "white", "white","white"
                ),
                filepath = "output_forestplot/logexp_RE_sub_ACE.png")

logexp_draw_sub(table_sub_RE_Chao1, n=nrow(table_sub_RE_Chao1),name = "Chao1",
                colorvector = c("grey90",
                                "white", "white", "white",
                                "grey90","grey90","grey90",
                                "white", "white", "white", "white"),
                filepath = "output_forestplot/logexp_RE_sub_Chao1.png")

#### ++ plot Disease only  #######
table_sub_Disease_ACE <- subgroup_RE_Disease_ROM(a_ACE, metadata)
table_sub_Disease_Chao1 <- subgroup_RE_Disease_ROM(a_Chao1, metadata)

logexp_draw_sub(table_sub_Disease_ACE,n=nrow(table_sub_Disease_ACE),name="ACE", 
                colorvector = c("grey90",
                                "white", "white", "white","white"
                ),
                filepath = "output_forestplot/logexp_Disease_sub_ACE.png")

logexp_draw_sub(table_sub_Disease_Chao1, n=nrow(table_sub_Disease_Chao1),name = "Chao1",
                colorvector = c("grey90",
                                "white", "white", "white", "white"),
                filepath = "output_forestplot/logexp_Disease_sub_Chao1.png")



table_sub_FE_Disease_ACE <- subgroup_FE_Disease_ROM(a_ACE, metadata)

logexp_draw_sub(table_sub_FE_Disease_ACE,n=nrow(table_sub_FE_Disease_ACE),name="ACE", 
                colorvector = c("grey90",
                                "white", "white", "white","white"
                ),
                filepath = "output_forestplot/logexp_FE_Disease_sub_ACEs.png")

table_sub_FE_Disease_Chao1 <- subgroup_FE_Disease_ROM(a_Chao1, metadata)

logexp_draw_sub(table_sub_FE_Disease_Chao1,n=nrow(table_sub_FE_Disease_Chao1),name="Chao1", 
                colorvector = c("grey90",
                                "white", "white", "white","white"
                ),
                filepath = "output_forestplot/logexp_FE_Disease_sub_Chao1.png")



#### 2. phylum ####

#### + import data ####
P_Bacteroidetes <- read.table("input_rela_abd/P_Bacteroidetes.tsv", header = TRUE, sep = "\t")
P_Firmicutes <- read.table("input_rela_abd/P_Firmicutes.tsv", header = TRUE, sep = "\t")
P_Proteobacteria <- read.table("input_rela_abd/P_Proteobacteria.tsv", header = TRUE, sep = "\t")
P_Verrucomicrobia <- read.table("input_rela_abd/P_Verrucomicrobia.tsv", header = TRUE, sep = "\t")
P_Actinobacteria <- read.table("input_rela_abd/P_Actinobacteria.tsv", header = TRUE, sep = "\t")

####### + overall analysis ######
####### ++ make ROM tables #####
table_total_P_Bacteroidetes <-total_ROM(P_Bacteroidetes,  name = "P_Bacteroidetes")
table_total_P_Firmicutes <-total_ROM(P_Firmicutes,  name = "P_Firmicutes")
table_total_P_Proteobacteria <-total_ROM(P_Proteobacteria,  name = "P_Proteobacteria")
table_total_P_Verrucomicrobia <-total_ROM(P_Verrucomicrobia,  name = "P_Verrucomicrobia")
table_total_P_Actinobacteria <-total_ROM(P_Actinobacteria,  name = "P_Actinobacteria")




###### ++ import metadata #####
meta_P <- read.table("metadata_P.tsv", header = TRUE, sep = "\t")
list_P <- make_list(meta_P)
p_list_P <- plot_list(list_P)


####### +++ make big ROM table #####
table_total <- left_join(list_P[1], table_total_P_Bacteroidetes[,c(1,2,4,5)], by="id")
table_total <- left_join(table_total, table_total_P_Firmicutes[,c(1,2,4,5)], by="id") 
table_total <- left_join(table_total, table_total_P_Proteobacteria[,c(1,2,4,5)], by="id") 
table_total <- left_join(table_total, table_total_P_Verrucomicrobia[,c(1,2,4,5)], by="id") 
table_total <- left_join(table_total, table_total_P_Actinobacteria[,c(1,2,4,5)], by="id") 

rownames(table_total)<-table_total[,1]
table_total<-exp(table_total[,-1])


pvalue <- left_join(list_P[1], table_total_P_Bacteroidetes[,c(1,6)], by="id")
table_total <- left_join(pvalue, table_total_P_Firmicutes[,c(1,6)], by="id") 
table_total <- left_join(pvalue, table_total_P_Proteobacteria[,c(1,6)], by="id") 
table_total <- left_join(pvalue, table_total_P_Verrucomicrobia[,c(1,6)], by="id") 
table_total <- left_join(pvalue, table_total_P_Actinobacteria[,c(1,6)], by="id") 

rownames(pvalue)<-pvalue[,1]
pvalue<-pvalue[,-1]

  
#### ++ plot in logROM form #######
pl_Bacteroidetes <- check_plot_total(table_total_P_Bacteroidetes,list = list_P,  name = "P_Bacteroidetes",overall = T)
pl_Firmicutes <- check_plot_total(table_total_P_Firmicutes, list = list_P, name = "P_Firmicutes", overall = T)
pl_Proteobacteria <- check_plot_total(table_total_P_Proteobacteria, list = list_P, name = "P_Proteobacteria", overall = T)
pl_Verrucomicrobia <- check_plot_total(table_total_P_Verrucomicrobia,list = list_P,  name = "P_Verrucomicrobia", overall = F)
pl_Actinobacteria <-check_plot_total(table_total_P_Actinobacteria, list = list_P, name = "P_Actinobacteria", overall = F)


png("output_forestplot/total_phylum.png", width = 10, height = (nrow(list_P)+3)/5, units = "in", res = 300)
grid.arrange(p_list_P, pl_Bacteroidetes, pl_Firmicutes, pl_Proteobacteria, pl_Verrucomicrobia, pl_Actinobacteria, 
             ncol=6, 
             widths= c(25,14,14,14,14,14)#, heights = c(5,5,5,5,5,5,5)
)
dev.off()
library(lmerTest)

###### +++ F/B ratio graph #######
FBratio <- read.table("input_rela_abd/FBratio.tsv", header = TRUE, sep = "\t")
FBratio <- left_join(FBratio, metadata, by="id")

FBratio$disease_GIsep <- factor(FBratio$disease_GIsep, levels = c("GI", "GI-inflam", "Obesity", "Other"))
FBratio <- FBratio %>% arrange(disease_GIsep)
FBratio$id <- factor(FBratio$id, levels = rev(FBratio$id))

FBratio$disease_group <- factor(FBratio$disease_group, levels = c("GI", "Obesity", "Other"))
FBratio <- FBratio %>% arrange(disease_group)
FBratio$id <- factor(FBratio$id, levels = rev(FBratio$id))



colorvector = c("grey90","grey90","grey90","grey90",
                "white", "white", "white", "white", 
                "grey90","grey90","grey90","grey90","grey90","grey90","grey90","grey90")

FBratio_table <- ggplot(data=FBratio, aes(y=id))+
  ggtitle("    Study               Dis.     Treat.   Wks")+
  geom_rect( aes(ymin = as.numeric(id)-0.5, ymax= as.numeric(id)+0.5,NULL, NULL),
             xmin = -Inf, xmax = Inf,colour = NA,
             fill = rev(colorvector))+
  geom_text(aes(x = 0, label = Study), hjust = 0) +
  geom_text(aes(x = 1.1, label = disease_group)) +
  geom_text(aes(x = 1.6, label = i_type)) +
  geom_text(aes(x = 2, label = duration)) +
  theme_MicrobeR() +
  theme(panel.border = element_blank(), axis.line = element_line()) +
  theme(axis.text.x=element_text(color = "white"),
        axis.title.x = element_text(color = "white"),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.title.y = element_blank(),
        axis.line = element_blank(),
        legend.position = "none",
        plot.title = element_text(hjust=0, face = "bold"),
        plot.margin=unit(c(0.3,0,0.5,0), "cm"))


FBratio_forest <- ggplot(data = FBratio, aes(x=i_FB/c_FB, y=id))  +
  ggtitle("")+

  geom_rect( aes(ymin = as.numeric(id)-0.5, ymax= as.numeric(id)+0.5,NULL, NULL),
             xmin = -Inf, xmax = Inf,colour = NA,
             fill = rev(colorvector))+
  geom_vline(xintercept = 1, linetype="dashed", color="grey50") +
  # geom_hline(yintercept = 3, linetype="dashed", color="grey50")+
  # geom_errorbarh(aes(xmin=lower, xmax=upper), height=0 ) +
  # scale_color_manual(values=c(  "pink3","skyblue3","gray44"))  +
  geom_point()+#aes(size = stat))#+
  # scale_shape_manual(values=c(18, 5, 16, 1))+
  # scale_size_manual(values=c(2.5, 2, 1, 0.5))+
  # scale_y_discrete(limits=rev)+
  # scale_x_continuous(expand = c(0.05, 0.05))+
  # labs(x="log(ROM)")+
  theme_MicrobeR() +
  theme(panel.border = element_blank(), axis.line = element_line()) +
  theme(axis.text.x=element_text(hjust=1),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none",
        plot.title = element_text(hjust=0.5, face = "bold"),
        plot.margin=unit(c(0.3,0,0.5,0), "cm"))

png("output_forestplot/FBratio.png", width = 4, height = (nrow(FBratio)+3)/5, units = "in", res = 300)
grid.arrange(FBratio_table, FBratio_forest,
             ncol=2, 
             widths= c(25,14)#, heights = c(5,5,5,5,5,5,5)
)
dev.off()




#### 3. genus ####

#### + import data ####
G_Lactobacillus <- read.table("input_rela_abd/G_Lactobacillus.tsv", header = TRUE, sep = "\t")
G_Bifidobacterium <- read.table("input_rela_abd/G_Bifidobacterium.tsv", header = TRUE, sep = "\t")
G_Roseburia <- read.table("input_rela_abd/G_Roseburia.tsv", header = TRUE, sep = "\t")
G_Bacteroides <- read.table("input_rela_abd/G_Bacteroides.tsv", header = TRUE, sep = "\t")
G_Ruminococcus <- read.table("input_rela_abd/G_Ruminococcus.tsv", header = TRUE, sep = "\t")

####### + overall analysis ######
####### ++ make ROM tables #####
table_total_G_Lactobacillus <-total_ROM(G_Lactobacillus,  name = "G_Lactobacillus")
table_total_G_Bifidobacterium <-total_ROM(G_Bifidobacterium,  name = "G_Bifidobacterium")
table_total_G_Roseburia <-total_ROM(G_Roseburia,  name = "G_Roseburia")
table_total_G_Bacteroides <-total_ROM(G_Bacteroides,  name = "G_Bacteroides")
table_total_G_Ruminococcus <-total_ROM(G_Ruminococcus,  name = "G_Ruminococcus")

###### ++ import metadata #####
meta_G <- read.table("metadata_G.tsv", header = TRUE, sep = "\t")
list_G <- make_list(meta_G)
p_list_G <- plot_list(list_G)


####### +++ make big ROM table #####
table_total_G <- left_join(list_G[1], table_total_G_Lactobacillus[,c(1,2,4,5)], by="id")
table_total_G <- left_join(table_total_G, table_total_G_Bifidobacterium[,c(1,2,4,5)], by="id") 
table_total_G <- left_join(table_total_G, table_total_G_Roseburia[,c(1,2,4,5)], by="id") 
table_total_G <- left_join(table_total_G, table_total_G_Bacteroides[,c(1,2,4,5)], by="id") 
table_total_G <- left_join(table_total_G, table_total_G_Ruminococcus[,c(1,2,4,5)], by="id") 

rownames(table_total_G)<-table_total_G[,1]
table_total_G<-exp(table_total_G[,-1])






#### ++ plot in logROM form #######
pG_Lactobacillus <- check_plot_total(table_total_G_Lactobacillus,list = list_G,  name = "G_Lactobacillus", overall = T)
pG_Bifidobacterium <- check_plot_total(table_total_G_Bifidobacterium,list = list_G,  name = "G_Bifidobacterium", overall = T)
pG_Roseburia <- check_plot_total(table_total_G_Roseburia,list = list_G,  name = "G_Roseburia", overall = F)
pG_Bacteroides <-check_plot_total(table_total_G_Bacteroides, list = list_G, name = "G_Bacteroides", overall = F)
pG_Ruminococcus <- check_plot_total(table_total_G_Ruminococcus, list = list_G, name = "G_Ruminococcus", overall = F)


png("output_forestplot/total_genus.png", width = 10, height = (nrow(list_G)+3)/5, units = "in", res = 300)
grid.arrange(p_list_G, pG_Lactobacillus,  pG_Bifidobacterium,  pG_Roseburia, pG_Bacteroides, pG_Ruminococcus,
             ncol=6, 
             widths= c(25,14,14,14,14,14)#, heights = c(5,5,5,5,5,5,5)
)
dev.off()



####### + subgroup analysis ######
####### ++ make ROM tables #####
table_sub_G_Lactobacillus <- subgroup_ROM(G_Lactobacillus, metadata)


#### ++ plot in logROM form #######
logexp_draw_sub(table_sub_G_Lactobacillus, n=nrow(table_sub_G_Lactobacillus),name = "Lactobacillus",
             colorvector = c("grey90","grey90",
                             "white", "white", "white", "white", 
                             "grey90","grey90","grey90",
                             "white", "white", "white", "white"#,
                             #"grey90","grey90","grey90","grey90","grey90"
                             ),
             filepath = "output_forestplot/logexp_sub_G_Lactobacillus.png")


table_sub_RE_G_Lactobacillus <- subgroup_RE_ROM(G_Lactobacillus, metadata)
logexp_draw_sub(table_sub_RE_G_Lactobacillus, n=nrow(table_sub_G_Lactobacillus),name = "Lactobacillus",
                colorvector = c("grey90",
                                "white", "white", "white", "white", 
                                "grey90","grey90","grey90",
                                "white", "white", "white", "white"#,
                                #"grey90","grey90","grey90","grey90","grey90"
                ),
                filepath = "output_forestplot/logexp_sub_RE_G_Lactobacillus.png")


table_sub_Disease_Lactobacillus <- subgroup_RE_Disease_ROM(G_Lactobacillus, metadata)

logexp_draw_sub(table_sub_Disease_Lactobacillus,n=nrow(table_sub_Disease_Lactobacillus),name="Lactobacillus", 
                colorvector = c("grey90",
                                "white", "white", "white","white"
                ),
                filepath = "output_forestplot/logexp_Disease_sub_Lactobacillus.png")



table_sub_FE_Disease_Lactobacillus <- subgroup_FE_Disease_ROM(G_Lactobacillus, metadata)

logexp_draw_sub(table_sub_FE_Disease_Lactobacillus,n=nrow(table_sub_FE_Disease_Lactobacillus),name="Lactobacillus", 
                colorvector = c("grey90",
                                "white", "white", "white","white"
                ),
                filepath = "output_forestplot/logexp_FE_Disease_sub_Lactobacillus.png")



