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
alpha_ACE <- read.table("input_alpha_div/alpha_ACE.tsv", header = TRUE, sep = "\t")
P_Firmicutes <- read.table("input_rela_abd/P_Firmicutes.tsv", header = TRUE, sep = "\t")




#merge metadata to wanted data
alpha_ACE <- left_join(alpha_ACE, metadata, by="id")
P_Firmicutes <- left_join(P_Firmicutes, metadata, by="id")

#run meta analysis
metacont_ACE <- metacont(i_num, i_mean, i_sd, c_num, c_mean, c_sd, 
         sm = "SMD", method.smd = "Hedges", id,
         data = alpha_ACE_meta)

meta_P_Firmicutes <- metacont(i_num, i_mean, i_sd, c_num, c_mean, c_sd, 
                         sm = "SMD", method.smd = "Hedges", id,
                         data = P_Firmicutes)


######### ROM
ROM_P_Firmicutes <- metacont(i_num, i_mean, i_sd, c_num, c_mean, c_sd, 
                             sm = "ROM", id,
                             #backtransf # TRUE (default), results will be presented as ratio of means; otherwise log ratio of means will be shown.
                             data = P_Firmicutes)
tROM_P_Firmicutes <- table_ROM(ROM_P_Firmicutes)
check_plot(tROM_P_Firmicutes)


ggplot(data = tROM_P_Firmicutes, aes(x=TE, y=studlab, color = signif, shape = signif)) +
  geom_vline(xintercept = 0, linetype="dashed", color="grey50") +
  geom_errorbarh(aes(xmin=lower, xmax=upper, alpha=!is.na(pval)), height=0 ) +
  theme_MicrobeR() +
  scale_color_manual(values=c("skyblue3", "pink3", "skyblue3","pink3", "skyblue3","pink3"))  +
  geom_point(aes(size = signif))+
  scale_shape_manual(values=c(1, 1, 16, 16, 18, 18))+
  scale_size_manual(values=c(1.5, 1.5, 1.5, 1.5, 3, 3))+
  scale_y_discrete(limits=rev)+
  theme(panel.border = element_blank(), axis.line = element_line()) +
  theme(axis.text.x=element_text(angle=45, hjust=1))
ggsave("output/rela_abd/P_Firmicutes.png", device = "png", width = 3, height = 3, units = "in")






ROM_alpha_ACE <- metacont(i_num, i_mean, i_sd, c_num, c_mean, c_sd, 
                             sm = "ROM", id,
                             #backtransf # TRUE (default), results will be presented as ratio of means; otherwise log ratio of means will be shown.
                             data = alpha_ACE)
tROM_alpha_ACE <- table_ROM(ROM_alpha_ACE)

check_plot(tROM_alpha_ACE)
ggplot(data = tROM_alpha_ACE, aes(x=TE, y=studlab, color = signif, shape = signif)) +
  geom_vline(xintercept = 0, linetype="dashed", color="grey50") +
  geom_errorbarh(aes(xmin=lower, xmax=upper), height=0 ) +
  theme_MicrobeR() +
  scale_color_manual(values=c("skyblue3", "pink3", "gray44", "skyblue3","gray44"))  +
  geom_point()+
  scale_shape_manual(values=c(16, 16, 16, 18, 18))+
  scale_y_discrete(limits=rev)+
  theme(panel.border = element_blank(), axis.line = element_line()) +
  theme(axis.text.x=element_text(angle=45, hjust=1))
ggsave("output/alpha/forest_ACE.png", device = "png", width = 3, height = 3, units = "in")







##