###### functions #####
theme_set(
  theme_bw()+
    theme(axis.line = element_line(size=1),
          axis.ticks = element_line(size=1),
          panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.title = element_text(hjust = 0.5, face="bold"),
          plot.subtitle = element_text(hjust = 0.5)
    ))

## SMD results are $TE

table_SMD <- function(SMD){
  
  
  
  df <- as.data.frame(SMD)
  
  m <- df[,7:ncol(df)]
  
  fixed <- data.frame(c("Fixed_Effect", SMD[c("TE.fixed", "seTE.fixed", "lower.fixed", "upper.fixed", "statistic.fixed", "pval.fixed", "zval.fixed")], 100,NA))
  random <- data.frame(c("Random_Effect", SMD[c("TE.random", "seTE.random", "lower.random", "upper.random", "statistic.random", "pval.random", "zval.random")], NA,100))
  
  colnames(fixed) <- colnames(m)
  colnames(random) <- colnames(m)
  m <- rbind(m, fixed, random)
  
  signif <- NULL
  for (i in 1:nrow(m)) {
    if (is.na(m$pval[i])){
      signif <- rbind(signif, NA)
    } else if (m$pval[i]>=0.05) {
      signif <- rbind(signif, "ns")
    } else {
      if (m$TE[i]<1) {
        signif <- rbind(signif, "treatment")
      } else {
        signif <- rbind(signif, "control")
      }
    }
    
    if (is.na(m$TE[i])){
      m$TE[i] <- df$mean.e[i]-df$mean.c[i]
    }
  }
  
  signif <- factor(signif, levels = c("treatment", "control", "ns", NA), order=TRUE)
  m$signif <- signif
  m$studlab <- factor(m$studlab, levels = m$studlab, order=TRUE)
  return(m)
}




table_ROM <- function(ROM){

  df <- as.data.frame(ROM)
  m <- df[,7:ncol(df)]
  
  fixed <- data.frame(c("Fixed_Effect", ROM[c("TE.fixed", "seTE.fixed", "lower.fixed", "upper.fixed", "statistic.fixed", "pval.fixed", "zval.fixed")], NA,NA))
  random <- data.frame(c("Random_Effect", ROM[c("TE.random", "seTE.random", "lower.random", "upper.random", "statistic.random", "pval.random", "zval.random")], NA,NA))
  
  colnames(fixed) <- colnames(m)
  colnames(random) <- colnames(m)
  
  if (fixed$pval>=0.05) {
    fixed$stat <- "sum"    
    fixed$signif <- "ns"
  } else {
    fixed$stat <- "sum"
    fixed$signif <- ifelse(fixed$TE<0, "e", "c")
  }
  
  if (random$pval>=0.05) {
    random$stat <- "sum"    
    random$signif <- "ns"
  } else {
    random$stat <- "sum"
    random$signif <- ifelse(random$TE<0, "e", "c")
  }
  
  signif <- NULL
  stat <- NULL
  
  for (i in 1:nrow(m)) {
    if (is.na(m$pval[i])){
      stat <- rbind(stat,"na")
      signif <- rbind(signif, ifelse(m$TE[i]<0, "e", "c"))
    } else if (m$pval[i]>=0.05) {
      stat <- rbind(stat, "stat")
      signif <- rbind(signif, "ns")
    } else {
      stat <- rbind(stat, "stat")
      signif <- rbind(signif, ifelse(m$TE[i]<0, "e", "c"))
    }
  }

  m$stat <- stat
  m$signif <- signif
  
  m <- rbind(m, fixed, random)
  m$stat <- factor(m$stat, levels = c('na', 'stat', 'sum'), order = TRUE)
  m$signif <- factor(m$signif, levels = c('e', 'c', 'ns'), order = TRUE)
  
  # m$signif <- factor(m$signif, levels = c("na.e", "na.c", "na.ns",
  #                                     "stat.e", "stat.c", "stat.ns",
  #                                     "sum.e", "sum.c", "sum.ns"), order=TRUE)
  m$studlab <- factor(m$studlab, levels = m$studlab, order=TRUE)
  return(m)
}


### grouped table

gtROM <- function (aPG) {
  aPG <- left_join(aPG, metadata, by="id")
  gROM <- metacont(i_num, i_mean, i_sd, c_num, c_mean, c_sd, 
                   sm = "ROM", id, subgroup = disease,
                   #backtransf # TRUE (default), results will be presented as ratio of means; otherwise log ratio of means will be shown.
                   data = aPG)
  gtROM <- gtable_ROM(gROM)
  return(gtROM)
}


ggtROM <- function (aPG) {
  aPG <- left_join(aPG, metadata, by="id")
  gROM <- metacont(i_num, i_mean, i_sd, c_num, c_mean, c_sd, 
                   sm = "ROM", id, subgroup = disease_group,
                   #backtransf # TRUE (default), results will be presented as ratio of means; otherwise log ratio of means will be shown.
                   data = aPG)
  gtROM <- ggtable_ROM(gROM)
  return(gtROM)
}






gtable_ROM <- function(gROM){

  df <- as.data.frame(gROM)
  m <- df[,7:(ncol(df)-1)]
  
  #difference between fixed and random effect
  
  fixed <- data.frame(c("Fixed_Effect", gROM[c("TE.fixed", "seTE.fixed", "lower.fixed", "upper.fixed", "statistic.fixed", "pval.fixed", "zval.fixed")], 
                        NA,NA, "FE", "sum"))
  random <- data.frame(c("Random_Effect", gROM[c("TE.random", "seTE.random", "lower.random", "upper.random", "statistic.random", "pval.random", "zval.random")], 
                         NA,NA, "RE", "sum"))
  
  gf <- data.frame(c(gROM[c("bylevs", "TE.fixed.w", "seTE.fixed.w", 
                                "lower.fixed.w", "upper.fixed.w", "statistic.fixed.w", "pval.fixed.w", "zval.fixed.w")],  
                     gROM[c("w.fixed.w")],NA, gROM[c("bylevs")]))
  # gf$bylevs  <- str_c(gf$bylevs, "_f")
  # 
  # gr <- data.frame(c(gROM[c("bylevs", "TE.random.w", "seTE.random.w", 
  #                               "lower.random.w", "upper.random.w", "statistic.random.w", "pval.random.w", "zval.random.w")], 
  #                    NA,gROM[c("w.random.w")]))
  # gr$bylevs  <- str_c(gr$bylevs, "_r")
  ## assumed subgroup are fixed effect
  
  if (fixed$pval.fixed>=0.05) {
    fixed$signif <- "ns"
  } else {
    fixed$signif <- ifelse(fixed$TE.fixed<0, "e", "c")
  }
  
  if (random$pval.random>=0.05) {
    random$signif <- "ns"
  } else {
    random$signif <- ifelse(random$TE.random<0, "e", "c")
  }
  stat <- NULL
  signif <- NULL

  for (i in 1:nrow(gf)) {
    if (is.na(gf$pval.fixed.w[i])){
      stat <- rbind(stat,"sub")
      signif <- rbind(signif, ifelse(gf$TE.fixed.w[i]<0, "e", "c"))
    } else if (gf$pval.fixed.w[i]>=0.05) {
      stat <- rbind(stat, "sub")
      signif <- rbind(signif, "ns")
    } else {
      stat <- rbind(stat, "sub")
      signif <- rbind(signif, ifelse(gf$TE.fixed.w[i]<0, "e", "c"))
    }
  }
  
  gf$stat <- stat
  gf$signif <- signif
  
  signif <- NULL
  stat <- NULL
  
  for (i in 1:nrow(m)) {
    if (is.na(m$pval[i])){
      stat <- rbind(stat,"na")
      signif <- rbind(signif, ifelse(m$TE[i]<0, "e", "c"))
    } else if (m$pval[i]>=0.05) {
      stat <- rbind(stat, "stat")
      signif <- rbind(signif, "ns")
    } else {
      stat <- rbind(stat, "stat")
      signif <- rbind(signif, ifelse(m$TE[i]<0, "e", "c"))
    }
  }
  
  m$stat <- stat
  m$signif <- signif
  
  colnames(fixed) <- colnames(m)
  colnames(random) <- colnames(m)
  colnames(gf) <- colnames(m)
  #colnames(gr) <- colnames(m)
  
  m <- rbind(m, gf, fixed, random)
  m$stat <- factor(m$stat, levels = c('sum', 'sub', 'stat', 'na' ), order = TRUE)
  m$signif <- factor(m$signif, levels = c('ns', 'e', 'c' ), order = TRUE)
  m$subgroup <- factor(m$subgroup, levels = c('Obese', 'Colitis', 'IBS', 'Other_GI', 'Other', 'Healthy', 'FE', 'RE'), order = TRUE)
  # m$signif <- factor(m$signif, levels = c("na.e", "na.c", "na.ns",
  #                                     "stat.e", "stat.c", "stat.ns",
  #                                     "sum.e", "sum.c", "sum.ns"), order=TRUE)
  m <- m[order(m$subgroup),]
  m$studlab <- factor(m$studlab, levels = m$studlab, order=TRUE)
  
  return(m)
}

ggtable_ROM <- function(gROM){
  
  df <- as.data.frame(gROM)
  m <- df[,7:(ncol(df)-1)]
  
  #difference between fixed and random effect
  
  fixed <- data.frame(c("Fixed_Effect", gROM[c("TE.fixed", "seTE.fixed", "lower.fixed", "upper.fixed", "statistic.fixed", "pval.fixed", "zval.fixed")], 
                        NA,NA, "FE", "sum"))
  random <- data.frame(c("Random_Effect", gROM[c("TE.random", "seTE.random", "lower.random", "upper.random", "statistic.random", "pval.random", "zval.random")], 
                         NA,NA, "RE", "sum"))
  
  gf <- data.frame(c(gROM[c("bylevs", "TE.fixed.w", "seTE.fixed.w", 
                            "lower.fixed.w", "upper.fixed.w", "statistic.fixed.w", "pval.fixed.w", "zval.fixed.w")],  
                     gROM[c("w.fixed.w")],NA, gROM[c("bylevs")]))
  # gf$bylevs  <- str_c(gf$bylevs, "_f")
  # 
  # gr <- data.frame(c(gROM[c("bylevs", "TE.random.w", "seTE.random.w", 
  #                               "lower.random.w", "upper.random.w", "statistic.random.w", "pval.random.w", "zval.random.w")], 
  #                    NA,gROM[c("w.random.w")]))
  # gr$bylevs  <- str_c(gr$bylevs, "_r")
  ## assumed subgroup are fixed effect
  
  if (fixed$pval.fixed>=0.05) {
    fixed$signif <- "ns"
  } else {
    fixed$signif <- ifelse(fixed$TE.fixed<0, "e", "c")
  }
  
  if (random$pval.random>=0.05) {
    random$signif <- "ns"
  } else {
    random$signif <- ifelse(random$TE.random<0, "e", "c")
  }
  stat <- NULL
  signif <- NULL
  
  for (i in 1:nrow(gf)) {
    if (is.na(gf$pval.fixed.w[i])){
      stat <- rbind(stat,"sub")
      signif <- rbind(signif, ifelse(gf$TE.fixed.w[i]<0, "e", "c"))
    } else if (gf$pval.fixed.w[i]>=0.05) {
      stat <- rbind(stat, "sub")
      signif <- rbind(signif, "ns")
    } else {
      stat <- rbind(stat, "sub")
      signif <- rbind(signif, ifelse(gf$TE.fixed.w[i]<0, "e", "c"))
    }
  }
  
  gf$stat <- stat
  gf$signif <- signif
  
  signif <- NULL
  stat <- NULL
  
  for (i in 1:nrow(m)) {
    if (is.na(m$pval[i])){
      stat <- rbind(stat,"na")
      signif <- rbind(signif, ifelse(m$TE[i]<0, "e", "c"))
    } else if (m$pval[i]>=0.05) {
      stat <- rbind(stat, "stat")
      signif <- rbind(signif, "ns")
    } else {
      stat <- rbind(stat, "stat")
      signif <- rbind(signif, ifelse(m$TE[i]<0, "e", "c"))
    }
  }
  
  m$stat <- stat
  m$signif <- signif
  
  colnames(fixed) <- colnames(m)
  colnames(random) <- colnames(m)
  colnames(gf) <- colnames(m)
  #colnames(gr) <- colnames(m)
  
  m <- rbind(m, gf, fixed, random)
  m$stat <- factor(m$stat, levels = c('sum', 'sub', 'stat', 'na' ), order = TRUE)
  m$signif <- factor(m$signif, levels = c('ns', 'e', 'c' ), order = TRUE)
  ##########m$subgroup <- factor(m$subgroup, levels = c('Obesity', 'GI', 'Other', 'Healthy', 'FE', 'RE'), order = TRUE)
  # m$signif <- factor(m$signif, levels = c("na.e", "na.c", "na.ns",
  #                                     "stat.e", "stat.c", "stat.ns",
  #                                     "sum.e", "sum.c", "sum.ns"), order=TRUE)
  ###########m <- m[order(m$subgroup),]
  m$studlab <- factor(m$studlab, levels = m$studlab, order=TRUE)
  
  return(m)
}









# # # # # # # # # 

# plot_SMD <- function(SMD){
#   ggplot(data = SMD, aes(x=TE, y=studlab, color = signif)) +
#     geom_vline(xintercept = 0, linetype="dashed", color="grey50") +
#     geom_errorbarh(aes(xmin=lower, xmax=upper), height=0 ) +
#     geom_point() +
#     #facet_grid(~Metric, scales="free_x") +
#     theme_MicrobeR() +
#     scale_color_manual(values=c("skyblue", "pink", "black", "grey"))  +
#     scale_y_discrete(limits=rev)+
#     theme(panel.border = element_blank(), axis.line = element_line()) +
#     theme(axis.text.x=element_text(angle=45, hjust=1))
# }

plot_ROM <- function(tROM){
  ggplot(data = tROM, aes(x=TE, y=id, color = signif)) +
    geom_vline(xintercept = 0, linetype="dashed", color="grey50") +
    geom_errorbarh(aes(xmin=lower, xmax=upper), height=0 ) +
    #geom_point(aes(alpha=is.na(signif)), shape = 1) +
    #geom_point(aes(alpha=!is.na(signif)), shape = 16) +
    #facet_grid(~Metric, scales="free_x") +
    theme_MicrobeR() +
    scale_color_manual(values=c("blue", "red", "black",
                                "blue", "red", "black",
                                "blue", "red", "black"))  +
    geom_point(aes())+
    scale_shape_manual(values=c(1, 1, 1, 16, 16, 16, 18, 18, 18))+
    scale_y_discrete(limits=rev)+
    theme(panel.border = element_blank(), axis.line = element_line()) +
    theme(axis.text.x=element_text(angle=45, hjust=1))
}

color_factor <- data.frame()
color_factor <- cbind(c("na.e", "na.c", "na.ns",
                         "stat.e", "stat.c", "stat.ns",
                         "sum.e", "sum.c", "sum.ns"),
                         c("blue", "red", "black",
                        "blue", "red", "black",
                        "blue", "red", "black"),
                        c(1,1,1,16,16,16,18,18,18))
colnames(color_factor) <- c("signif", "color", "shape")



check_plot <- function(tROM, filepath, width = w, height = h){
  g<- ggplot(data = tROM, aes(x=TE, y=studlab, color = signif, shape = stat)) +
    ggtitle(deparse(substitute(tROM)))+
    geom_vline(xintercept = 0, linetype="dashed", color="grey50") +
    geom_errorbarh(aes(xmin=lower, xmax=upper), height=0 ) +
    scale_color_manual(values=c( "gray44", "skyblue3", "pink3"))  +
    geom_point(aes(size = stat))+
    scale_shape_manual(values=c(18, 5, 16, 1))+
    scale_size_manual(values=c(2.5, 2, 1, 0.5))+
    scale_y_discrete(limits=rev)+
    theme(panel.border = element_blank(), axis.line = element_line()) +
    theme(axis.text.x=element_text(hjust=1))

}
check_plot(gtROM_ACE)






########## RENEW FIGURES ################
########## + total RE, FE data table ################


total_ROM <- function(a_DATA, name = "DATA"){
   # name = "DATA"
   # a_DATA <- a_Chao1
  ROM_DATA <- metacont(i_num, i_mean, i_sd, c_num, c_mean, c_sd, 
                        sm = "ROM", id, #subgroup = disease,
                        #backtransf # TRUE (default), results will be presented as ratio of means; otherwise log ratio of means will be shown.
                        data = a_DATA)
  df_DATA <- as.data.frame(ROM_DATA)
  table_DATA <- data.frame(df_DATA$studlab)    
  table_DATA$TE <- df_DATA$TE
  table_DATA$seTE <- df_DATA$seTE
  table_DATA$lower <- df_DATA$lower
  table_DATA$upper <- df_DATA$upper
  table_DATA$pval <- df_DATA$pval
  table_DATA$w.common <- df_DATA$w.common
  table_DATA$w.random <- df_DATA$w.random
  
  fixed <- data.frame(c("Fixed.Effect", ROM_DATA[c("TE.fixed", "seTE.fixed", "lower.fixed", "upper.fixed", "pval.fixed")], NA,NA))
  random <- data.frame(c("Random.Effect", ROM_DATA[c("TE.random", "seTE.random","lower.random", "upper.random", "pval.random")], NA,NA))
  colnames(fixed) <-colnames(table_DATA)
  colnames(random) <- colnames(table_DATA)
  table_DATA <- rbind(table_DATA, fixed, random)
  
  signif <- NULL
  stat <- NULL
  
  for (i in 1:nrow(table_DATA)) {
    if (is.na(table_DATA$pval[i])){
      stat <- rbind(stat,"na")
      signif <- rbind(signif, ifelse(table_DATA$TE[i]<0, "e", "c"))
    } else if (table_DATA$pval[i]>=0.05) {
      stat <- rbind(stat, "stat")
      signif <- rbind(signif, "ns")
    } else {
      stat <- rbind(stat, "stat")
      signif <- rbind(signif, ifelse(table_DATA$TE[i]<0, "e", "c"))
    }
  }
  
  #signif <- factor(signif, levels = c("e", "c", "ns", NA), order=TRUE)  
  #stat <- factor(stat, levels = c("stat", "na", NA), order=TRUE)  
  table_DATA <- cbind(table_DATA, stat, signif)
  
  # coln <- paste(name, colnames(table_DATA), sep = ".")
  # coln[1]<- c("id")
  # colnames(table_DATA) <- c(coln)
  colnames(table_DATA)[1] <- c("id")
  
  # table_DATA$id <- factor(table_DATA$id, levels = table_DATA$id, order=TRUE)
  
  return(table_DATA)
}

check_plot_total <- function(total_ROM, list = list, name = "me", overall = T){
  # tROM = list_Chao1
  # name = "Chao1"
  tROM <- left_join(list, total_ROM, by = "id")

  tROM$signif <- factor(tROM$signif, levels = c("e", "c", "ns", NA), order=TRUE)  
  tROM$stat <- factor(tROM$stat, levels = c("stat", "na", NA), order=TRUE)
  tROM$id <- factor(tROM$id, levels = tROM$id, order=TRUE)
  if (overall == T) {
  ggplot(data = tROM, aes(x=TE, y=id, color = signif, shape = stat))  +
    ggtitle(name)+
    geom_vline(xintercept = 0, linetype="dashed", color="grey50") +
    geom_hline(yintercept = 3, linetype="dashed", color="grey50")+
    geom_errorbarh(aes(xmin=lower, xmax=upper), height=0 ) +
    scale_color_manual(values=c(  "pink3","skyblue3","gray44"))  +
    geom_point()+#aes(size = stat))#+
    scale_shape_manual(values=c(18, 5, 16, 1))+
    scale_size_manual(values=c(2.5, 2, 1, 0.5))+
    scale_y_discrete(limits=rev)+
    scale_x_continuous(expand = c(0.05, 0.05))+
      labs(x="log(ROM)")+
    theme_MicrobeR() +
    theme(panel.border = element_blank(), axis.line = element_line()) +
    theme(axis.text.x=element_text(hjust=1),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.y = element_blank(),
          legend.position = "none",
          plot.title = element_text(hjust=0.5, face = "bold"))
  } else {
    ggplot(data = tROM, aes(x=TE, y=id, color = signif, shape = stat, alpha = id))  +
      ggtitle(name)+
      geom_vline(xintercept = 0, linetype="dashed", color="grey50") +
      geom_hline(yintercept = 3, linetype="dashed", color="grey50")+
      geom_errorbarh(aes(xmin=lower, xmax=upper), height=0 ) +
      scale_color_manual(values=c(  "pink3","skyblue3","gray44"))  +
      geom_point()+#aes(size = stat))#+
      scale_shape_manual(values=c(18, 5, 16, 1))+
      scale_size_manual(values=c(2.5, 2, 1, 0.5))+
      scale_alpha_manual(values= c(Fixed.Effect = 0, Random.Effect = 0))+
      scale_y_discrete(limits=rev)+
      scale_x_continuous(expand = c(0.05, 0.05))+
      labs(x="log(ROM)")+
      theme_MicrobeR() +
      theme(panel.border = element_blank(), axis.line = element_line()) +
      theme(axis.text.x=element_text(hjust=1),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.title.y = element_blank(),
            legend.position = "none",
            plot.title = element_text(hjust=0.5, face = "bold"))
  }
  
}


check_plot_total(total_Chao1, name = "Chao1")

make_list <- function(metadata){
  
  list <- metadata[,c(1:7)]
  list <- rbind(list, 
                c("Blank", NA, NA, NA, NA, NA, NA),
                c( "Random.Effect","Overall.(RE)", NA, NA, NA, NA, NA), 
                c( "Fixed.Effect" , "Overall.(FE)",NA, NA, NA, NA, NA))
  list$id <- factor(list$id, levels = list$id, order=TRUE)
  return(list)
}

plot_list <- function(list){
#  list$id <- factor(list$id, levels = rev(list$id), order=TRUE)
  p_list <- 
    ggplot(data=list, aes(y=rev(id)))+
    ggtitle("    Study               Dis.     Treat.   Wks")+
    #geom_hline(aes(yintercept = labels, colour = colour), size = 7) +
    geom_text(aes(x = 0, label = Study), hjust = 0) +
    #geom_text(aes(x = 1, label = species)) +
    geom_text(aes(x = 1.1, label = disease_group)) +
    geom_text(aes(x = 1.6, label = i_type)) +
    geom_text(aes(x = 2, label = duration)) +
    #scale_colour_identity() +
    #theme_void() + 
    theme_MicrobeR() +
    theme(panel.border = element_blank(), axis.line = element_line()) +
    theme(axis.text.x=element_text(color = "white"),
          axis.title.x = element_text(color = "white"),
          axis.text.y = element_blank(),
          axis.ticks = element_blank(),
          axis.title.y = element_blank(),
          axis.line = element_blank(),
          legend.position = "none",
          plot.title = element_text(hjust=0, face = "bold"))
  
  return(p_list)
}

########## + subgroup RE, FE data table ################

subgroup_ROM <- function(a_DATA, metadata){
  
  #a_DATA <- G_Bifidobacterium
  g_DATA <- left_join(a_DATA, metadata)
  
  g_list <- c("Q","df.Q", "pval.Q","H", "I2", "tau2", 
              "TE.ROM", "lower.ROM", "upper.ROM", "pval.ROM")
  
  total_DATA <- metacont(i_num, i_mean, i_sd, c_num, c_mean, c_sd, 
                         sm = "ROM", id, #subgroup = disease_group,
                         #backtransf # TRUE (default), results will be presented as ratio of means; otherwise log ratio of means will be shown.
                         data = g_DATA)
  
  Fixed.Effect <- data.frame(c("Fixed.Effect",
                               total_DATA[c("Q","df.Q", "pval.Q","H", "I2", "tau2", "k",
                                            "TE.fixed", "lower.fixed", "upper.fixed", "pval.fixed")], 1))
  Random.Effect <- data.frame(c("Random.Effect",
                                total_DATA[c("Q","df.Q", "pval.Q","H", "I2", "tau2", "k",
                                             "TE.random", "lower.random", "upper.random", "pval.random")], 1))
  
  
  
  
  
  
  treatment_DATA <- metacont(i_num, i_mean, i_sd, c_num, c_mean, c_sd, 
                             sm = "ROM", id, subgroup = i_type,
                             #backtransf # TRUE (default), results will be presented as ratio of means; otherwise log ratio of means will be shown.
                             data = g_DATA)
  
  Treatment <- data.frame(c("Treatment", treatment_DATA[c("Q.b.random","df.Q.b", "pval.Q.b.random")],NA, NA, NA, NA,
                            NA, NA, NA, NA, 2))
  
  Treatment_subgroup <- data.frame(c(treatment_DATA[c("bylevs","Q.w")], NA, 
                                     treatment_DATA[c("pval.Q.w")], NA,
                                     treatment_DATA[c("I2.w", "tau2.w","k.w",
                                                      "TE.fixed.w", "lower.fixed.w", "upper.fixed.w", "pval.fixed.w")], 3))
  # species_DATA <- metacont(i_num, i_mean, i_sd, c_num, c_mean, c_sd, 
  #                          sm = "ROM", id, subgroup = species,
  #                          #backtransf # TRUE (default), results will be presented as ratio of means; otherwise log ratio of means will be shown.
  #                          data = g_DATA)
  # 
  # Species <- data.frame(c("Species",species_DATA[c("Q.b.random","df.Q.b", "pval.Q.b.random")],NA, NA, NA, NA,
  #                         NA, NA, NA, NA, 2))
  # 
  # Species_subgroup <- data.frame(c(species_DATA[c("bylevs","Q.w")], NA, 
  #                                  species_DATA[c("pval.Q.w")], NA,
  #                                  species_DATA[c("I2.w", "tau2.w","k.w",
  #                                                 "TE.fixed.w", "lower.fixed.w", "upper.fixed.w", "pval.fixed.w")], 3))
  # 
  duration_DATA <- metacont(i_num, i_mean, i_sd, c_num, c_mean, c_sd, 
                             sm = "ROM", id, subgroup = duration,
                             #backtransf # TRUE (default), results will be presented as ratio of means; otherwise log ratio of means will be shown.
                             data = g_DATA)
   
   Duration <- data.frame(c("Duration", duration_DATA[c("Q.b.random","df.Q.b", "pval.Q.b.random")],NA, NA, NA, 
                            NA, NA, NA, NA, NA, 2))
   
   Duration_subgroup <- data.frame(c(duration_DATA[c("bylevs","Q.w")], NA, 
                                     duration_DATA[c("pval.Q.w")], NA,
                                     duration_DATA[c("I2.w", "tau2.w","k.w",
                                                     "TE.fixed.w", "lower.fixed.w", "upper.fixed.w", "pval.fixed.w")], 3))
   
   disease_DATA <- metacont(i_num, i_mean, i_sd, c_num, c_mean, c_sd, 
                            sm = "ROM", id, subgroup = disease_group,
                            #backtransf # TRUE (default), results will be presented as ratio of means; otherwise log ratio of means will be shown.
                            data = g_DATA)
   
   Disease <- data.frame(c("Disease", 
                           disease_DATA[c("Q.b.random","df.Q.b", "pval.Q.b.random")],
                           NA, NA, NA, NA,NA, NA, NA, NA, 2))
   
   Disease_subgroup <- data.frame(c(disease_DATA[c("bylevs", "Q.w")], NA, 
                                    disease_DATA[c("pval.Q.w")], NA,
                                    disease_DATA[c("I2.w", "tau2.w","k.w",
                                                   "TE.fixed.w", "lower.fixed.w", "upper.fixed.w", "pval.fixed.w")], 3))
   
  table_sub_DATA <- t(cbind(t(Fixed.Effect), t(Random.Effect), 
                            
                            t(Treatment), t(Treatment_subgroup),
                            # t(Species), t(Species_subgroup),
                             t(Duration), t(Duration_subgroup),
                            t(Disease), t(Disease_subgroup)
                            
                            
  ))
  
  
  colnames(table_sub_DATA) <- c("Groups", "Q" , "df.Q" ,"pval.Q","H" , "I2" ,"tau2","k",
                                "TE", "lower", "upper", "pval" ,"level" )
  #colnames(lll) <- c("Groups")
  
  table_sub_DATA <- data.frame(table_sub_DATA)
  
  table_sub_DATA$Groups <- factor(table_sub_DATA$Groups, 
                                  levels =c("Fixed.Effect", "Random.Effect", 
                                            
                                            "Treatment", "EA", "MA", "MOX", 
                                            "Duration", "≤2", ">2",
                                            "Disease", "GI", "Obesity", "Other"), order=TRUE)
  table_sub_DATA<- table_sub_DATA %>% arrange(Groups)
  table_sub_DATA$Groups <- factor(table_sub_DATA$Groups, 
                                  levels = table_sub_DATA$Groups, order=TRUE)
  
  table_sub_DATA$Q <- as.numeric(table_sub_DATA$Q)
  table_sub_DATA$pval.Q <- as.numeric(table_sub_DATA$pval.Q)
  
  table_sub_DATA$TE <- as.numeric(table_sub_DATA$TE)
  table_sub_DATA$lower <- as.numeric(table_sub_DATA$lower)
  table_sub_DATA$upper <- as.numeric(table_sub_DATA$upper)
  table_sub_DATA$pval <- as.numeric(table_sub_DATA$pval)
  
  signif <- NULL
  stat <- NULL
  
  for (i in 1:nrow(table_sub_DATA)) {
    if (is.na(table_sub_DATA$pval[i])){
      stat <- rbind(stat,"na")
      signif <- rbind(signif, ifelse(table_sub_DATA$TE[i]<0, "e", "c"))
    } else if (table_sub_DATA$pval[i]>=0.05) {
      stat <- rbind(stat, "stat")
      signif <- rbind(signif, "ns")
    } else {
      stat <- rbind(stat, "stat")
      signif <- rbind(signif, ifelse(table_sub_DATA$TE[i]<0, "e", "c"))
    }
  }
  
  #signif <- factor(signif, levels = c("e", "c", "ns", NA), order=TRUE)  
  #stat <- factor(stat, levels = c("stat", "na", NA), order=TRUE)  
  table_sub_DATA <- cbind(table_sub_DATA, stat, signif)
  
  
  return(table_sub_DATA)
}

subgroup_RE_ROM <- function(a_DATA, metadata){
  
  #a_DATA <- G_Bifidobacterium
  g_DATA <- left_join(a_DATA, metadata)
  
  g_list <- c("Q","df.Q", "pval.Q","H", "I2", "tau2", 
              "TE.ROM", "lower.ROM", "upper.ROM", "pval.ROM")
  
  total_DATA <- metacont(i_num, i_mean, i_sd, c_num, c_mean, c_sd, 
                         sm = "ROM", id, #subgroup = disease_group,
                         #backtransf # TRUE (default), results will be presented as ratio of means; otherwise log ratio of means will be shown.
                         data = g_DATA)
  
    Random.Effect <- data.frame(c("Random.Effect",
                                total_DATA[c("Q","df.Q", "pval.Q","H", "I2", "tau2", "k",
                                             "TE.random", "lower.random", "upper.random", "pval.random")], 1))
  
  treatment_DATA <- metacont(i_num, i_mean, i_sd, c_num, c_mean, c_sd, 
                             sm = "ROM", id, subgroup = i_type,
                             #backtransf # TRUE (default), results will be presented as ratio of means; otherwise log ratio of means will be shown.
                             data = g_DATA)
  
  Treatment <- data.frame(c("Treatment", treatment_DATA[c("Q.b.random","df.Q.b", "pval.Q.b.random")],NA, NA, NA, NA,
                            NA, NA, NA, NA, 2))
  
  Treatment_subgroup <- data.frame(c(treatment_DATA[c("bylevs","Q.w")], NA, 
                                     treatment_DATA[c("pval.Q.w")], NA,
                                     treatment_DATA[c("I2.w", "tau2.w","k.w",
                                                      "TE.random.w", "lower.random.w", "upper.random.w", "pval.random.w")], 3))
  # species_DATA <- metacont(i_num, i_mean, i_sd, c_num, c_mean, c_sd, 
  #                          sm = "ROM", id, subgroup = species,
  #                          #backtransf # TRUE (default), results will be presented as ratio of means; otherwise log ratio of means will be shown.
  #                          data = g_DATA)
  # 
  # Species <- data.frame(c("Species",species_DATA[c("Q.b.random","df.Q.b", "pval.Q.b.random")],NA, NA, NA, NA,
  #                         NA, NA, NA, NA, 2))
  # 
  # Species_subgroup <- data.frame(c(species_DATA[c("bylevs","Q.w")], NA, 
  #                                  species_DATA[c("pval.Q.w")], NA,
  #                                  species_DATA[c("I2.w", "tau2.w","k.w",
  #                                                 "TE.fixed.w", "lower.fixed.w", "upper.fixed.w", "pval.fixed.w")], 3))
  # 
  duration_DATA <- metacont(i_num, i_mean, i_sd, c_num, c_mean, c_sd, 
                            sm = "ROM", id, subgroup = duration,
                            #backtransf # TRUE (default), results will be presented as ratio of means; otherwise log ratio of means will be shown.
                            data = g_DATA)
  
  Duration <- data.frame(c("Duration", duration_DATA[c("Q.b.random","df.Q.b", "pval.Q.b.random")],NA, NA, NA, 
                           NA, NA, NA, NA, NA, 2))
  
  Duration_subgroup <- data.frame(c(duration_DATA[c("bylevs","Q.w")], NA, 
                                    duration_DATA[c("pval.Q.w")], NA,
                                    duration_DATA[c("I2.w", "tau2.w","k.w",
                                                    "TE.random.w", "lower.random.w", "upper.random.w", "pval.random.w")], 3))
  
  disease_DATA <- metacont(i_num, i_mean, i_sd, c_num, c_mean, c_sd, 
                           sm = "ROM", id, subgroup = disease_group,
                           #backtransf # TRUE (default), results will be presented as ratio of means; otherwise log ratio of means will be shown.
                           data = g_DATA)
  
  Disease <- data.frame(c("Disease", 
                          disease_DATA[c("Q.b.random","df.Q.b", "pval.Q.b.random")],
                          NA, NA, NA, NA,NA, NA, NA, NA, 2))
  
  Disease_subgroup <- data.frame(c(disease_DATA[c("bylevs", "Q.w")], NA, 
                                   disease_DATA[c("pval.Q.w")], NA,
                                   disease_DATA[c("I2.w", "tau2.w","k.w",
                                                  "TE.random.w", "lower.random.w", "upper.random.w", "pval.random.w")], 3))
  
  table_sub_DATA <- t(cbind(t(Random.Effect), 
                            
                            t(Treatment), t(Treatment_subgroup),
                            # t(Species), t(Species_subgroup),
                            t(Duration), t(Duration_subgroup),
                            t(Disease), t(Disease_subgroup)
                            
                            
  ))
  
  
  colnames(table_sub_DATA) <- c("Groups", "Q" , "df.Q" ,"pval.Q","H" , "I2" ,"tau2","k",
                                "TE", "lower", "upper", "pval" ,"level" )
  #colnames(lll) <- c("Groups")
  
  table_sub_DATA <- data.frame(table_sub_DATA)
  
  table_sub_DATA$Groups <- factor(table_sub_DATA$Groups, 
                                  levels =c("Fixed.Effect", "Random.Effect", 
                                            
                                            "Treatment", "EA", "MA", "MOX", 
                                            "Duration", "≤2", ">2",
                                            "Disease", "GI", "Obesity", "Other"), order=TRUE)
  table_sub_DATA<- table_sub_DATA %>% arrange(Groups)
  table_sub_DATA$Groups <- factor(table_sub_DATA$Groups, 
                                  levels = table_sub_DATA$Groups, order=TRUE)
  
  table_sub_DATA$Q <- as.numeric(table_sub_DATA$Q)
  table_sub_DATA$pval.Q <- as.numeric(table_sub_DATA$pval.Q)
  
  table_sub_DATA$TE <- as.numeric(table_sub_DATA$TE)
  table_sub_DATA$lower <- as.numeric(table_sub_DATA$lower)
  table_sub_DATA$upper <- as.numeric(table_sub_DATA$upper)
  table_sub_DATA$pval <- as.numeric(table_sub_DATA$pval)
  
  signif <- NULL
  stat <- NULL
  
  for (i in 1:nrow(table_sub_DATA)) {
    if (is.na(table_sub_DATA$pval[i])){
      stat <- rbind(stat,"na")
      signif <- rbind(signif, ifelse(table_sub_DATA$TE[i]<0, "e", "c"))
    } else if (table_sub_DATA$pval[i]>=0.05) {
      stat <- rbind(stat, "stat")
      signif <- rbind(signif, "ns")
    } else {
      stat <- rbind(stat, "stat")
      signif <- rbind(signif, ifelse(table_sub_DATA$TE[i]<0, "e", "c"))
    }
  }
  
  #signif <- factor(signif, levels = c("e", "c", "ns", NA), order=TRUE)  
  #stat <- factor(stat, levels = c("stat", "na", NA), order=TRUE)  
  table_sub_DATA <- cbind(table_sub_DATA, stat, signif)
  
  
  return(table_sub_DATA)
}











p1_sub<- function(table_sub_DATA, colorvector){
  #  list$id <- factor(list$id, levels = rev(list$id), order=TRUE)
  #table_sub_DATA <- table_sub_Chao1
  p1_sub_DATA <- 
    ggplot(data=table_sub_DATA, aes(y=rev(Groups)))+
    ggtitle("Q       pval Q        I2      tau2     k  ")+
    #geom_hline(aes(yintercept = labels, colour = colour), size = 7) +
    geom_rect( aes(ymin = as.numeric(Groups)-0.5, ymax= as.numeric(Groups)+0.5,NULL, NULL),
               xmin = -Inf, xmax = Inf,colour = NA,
               fill = rev(colorvector))+
    geom_text(aes(x = 0, label = ifelse(level==3, paste("  ", as.character(Groups)), 
                                        as.character(Groups)), 
                  hjust = 0)) +
    
    geom_text(aes(x = 1, label = ifelse(as.numeric(Q)==0, NA, round(as.numeric(Q), digits = 2))))+                #round(as.numeric(Q), digits = 2))) +
      
    geom_text(aes(x = 1.5, label = ifelse(round(as.numeric(pval.Q), digits = 4)==0, 
                                          "<0.0001", round(as.numeric(pval.Q), digits = 4)))) +
    geom_text(aes(x = 2, label = ifelse(as.numeric(I2)==0, NA, round(as.numeric(I2), digits = 2)))) +
    geom_text(aes(x = 2.4, label = ifelse(as.numeric(tau2)==0, NA, round(as.numeric(tau2), digits = 2)))) +
    geom_text(aes(x = 2.7 , label = k)) +
    #scale_colour_identity() +
    #theme_void() + 
    theme_MicrobeR() +
    theme(panel.border = element_blank(), axis.line = element_line()) +
      # scale_x_continuous(position = "top", breaks = c(1,1.5,2,2.4,2.7),
      #                      labels = c("Q", "pval Q", "I2", "tau2", "k"))+
    theme(axis.text.x=element_text(color = "white"),
          axis.title.x = element_text(color = "white"),
          axis.text.y = element_blank(),
          axis.ticks = element_blank(),
          axis.title.y = element_blank(),
          axis.line = element_blank(),
          legend.position = "none",
          plot.title = element_text(hjust=1, size = 12, face="bold"),
          plot.margin=unit(c(0.5,0,0.5,0), "cm"))
  return(p1_sub_DATA)
}


p2_sub <- function(table_sub_DATA, colorvector, name = ""){#}, list = list, name = "me"){
  # table_sub_DATA <- table_sub_Chao1
  #tROM <- left_join(list, total_ROM, by = "id")
  # colorvector = c("white", "white", 
  #                 "grey90","grey90","grey90","grey90",
  #                 "white", "white", "white",
  #                 "grey90","grey90","grey90","grey90",
  #                 "white","white","white", "white", "white", "white")
  table_sub_DATA$signif <- factor(table_sub_DATA$signif, levels = c("e", "c", "ns", NA), order=TRUE)  
  table_sub_DATA$stat <- factor(table_sub_DATA$stat, levels = c("stat", "na", NA), order=TRUE)
  p2_sub_DATA <- 
  ggplot(data = table_sub_DATA, aes(x=TE, y=Groups, color = signif, shape = stat))  +
    ggtitle(name)+
    geom_rect( aes(ymin = as.numeric(Groups)-0.5, ymax= as.numeric(Groups)+0.5,NULL, NULL),
               xmin = -Inf, xmax = Inf,colour = NA,
               fill = rev(colorvector))+
    geom_vline(xintercept = 0, linetype="dashed", color="grey50") +
    geom_errorbarh(aes(xmin=as.numeric(lower), xmax=as.numeric(upper), color = signif), height=0) +
    scale_color_manual(values=c( "e"="pink3","c"="skyblue3","ns"="gray44","NA"= "gray"))  +
    geom_point(aes(shape = stat, size = stat))+#aes(size = stat))#+
    scale_shape_manual(values=c(18, 5, "stat" = 16, "ns" = 1))+
    scale_size_manual(values=c(2.5, 2, "stat" = 1, "ns" = 0.5))+
    scale_y_discrete(limits=rev)+
    labs(x="log(ROM)")+
    theme_MicrobeR() +
    #scale_x_continuous(position = "top")+
    theme(panel.border = element_blank(), axis.line = element_line()) +
    theme(axis.text.x=element_text(color = "black"),
          axis.title.x = element_text(color = "black"),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.y = element_blank(),
          axis.line.y = element_blank(),
          legend.position = "none",
          plot.title = element_text(hjust=0.5, size = 12, face="bold"),
          plot.margin=unit(c(0.5,0,0.5,0), "cm"))
  return(p2_sub_DATA)
  
}


p3_sub<- function(table_sub_DATA, colorvector){
  #  list$id <- factor(list$id, levels = rev(list$id), order=TRUE)
  # table_sub_DATA <- table_sub_Chao1
  p3_sub_DATA <- 
  ggplot(data=table_sub_DATA, aes(y=rev(Groups)))+
    ggtitle("ROM [95%CI]     pval")+
    #geom_hline(aes(yintercept = labels, colour = colour), size = 7) +
    geom_rect( aes(ymin = as.numeric(Groups)-0.5, ymax= as.numeric(Groups)+0.5,NULL, NULL),
               xmin = -Inf, xmax = Inf,colour = NA,
               fill = rev(colorvector))+
    geom_text(aes(x = 0.4, label = format(round(TE, digits =2)), nsmall=2),hjust = 1) +
    geom_text(aes(x = 0.4, label = ifelse(is.na(lower),"", 
                                        paste( " [",format(round((lower), digits =2), nsmall=2), ", ",
                                               format(round((upper), digits =2), nsmall=2), "]", sep=""))), hjust = 0) +                #round(as.numeric(Q), digits = 2))) +
    
    geom_text(aes(x = 1.7, label = ifelse(round(as.numeric(pval), digits = 3)==0, 
                                          "<0.0001", format(round(as.numeric(pval), digits = 3), nsmall=3)))) +
    theme_MicrobeR() +
    theme(panel.border = element_blank(), axis.line = element_line()) +
    # scale_x_continuous(position = "top", breaks = c(1,1.5,2,2.4,2.7),
    #                      labels = c("Q", "pval Q", "I2", "tau2", "k"))+
    xlim(c(0,2))+
    theme(axis.text.x=element_text(color = "white"),
          axis.title.x = element_text(color = "white"),
          axis.text.y = element_blank(),
          axis.ticks = element_blank(),
          axis.title.y = element_blank(),
          axis.line = element_blank(),
          legend.position = "none",
          plot.title = element_text(hjust=0.5, size = 12, face="bold"),
          plot.margin=unit(c(0.5,0,0.5,0), "cm"))
  return(p3_sub_DATA)
}

draw_sub <- function(table_sub_DATA, name="",  n=number, colorvector, filepath = ".png"){
p1_sub_DATA <- p1_sub(table_sub_DATA, colorvector)
p2_sub_DATA <- p2_sub(table_sub_DATA, colorvector, name = "")
p3_sub_DATA <- p3_sub(table_sub_DATA, colorvector)

png(filepath, width = 12, height = (n+3)/5, units = "in", res = 300)
grid.arrange(p1_sub_DATA, p2_sub_DATA, p3_sub_DATA, 
             ncol=3, 
             widths= c(30,10,20)#, heights = c(5,5,5,5,5,5,5)
)
dev.off()
}



exp_check_plot_total <- function(total_ROM, list = list, name = "me"){
  # tROM = list_Chao1
  # name = "Chao1"
  tROM <- left_join(list, total_ROM, by = "id")
  
  tROM$signif <- factor(tROM$signif, levels = c("e", "c", "ns", NA), order=TRUE)  
  tROM$stat <- factor(tROM$stat, levels = c("stat", "na", NA), order=TRUE)
  tROM$id <- factor(tROM$id, levels = tROM$id, order=TRUE)
  
  ggplot(data = tROM, aes(x=exp(TE), y=id, color = signif, shape = stat))  +
    ggtitle(name)+
    geom_vline(xintercept = 1, linetype="dashed", color="grey50") +
    geom_hline(yintercept = 3, linetype="dashed", color="grey50")+
    geom_errorbarh(aes(xmin=exp(lower), xmax=exp(upper) ), height=0 ) +
    scale_color_manual(values=c("e"="pink3","c" ="skyblue3","ns"= "gray44"))  +
    geom_point()+#aes(size = stat))#+
    scale_shape_manual(values=c(18, 5,"stat"=  16, "na" = 1))+
    scale_size_manual(values=c(2.5, 2,"stat"=  1, "na"=0.5))+
    scale_y_discrete(limits=rev)+
    scale_x_continuous(expand = c(0.05, 0.05))+
    theme_MicrobeR() +
    theme(panel.border = element_blank(), axis.line = element_line()) +
    theme(axis.text.x=element_text(hjust=1),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.y = element_blank(),
          legend.position = "none",
          plot.title = element_text(hjust=0.5))
}


exp_p1_sub<- function(table_sub_DATA, colorvector){
  #  list$id <- factor(list$id, levels = rev(list$id), order=TRUE)
  #table_sub_DATA <- table_sub_Chao1
  p1_sub_DATA <- 
    ggplot(data=table_sub_DATA, aes(y=rev(Groups)))+
    ggtitle("    Study                      n")+
    #geom_hline(aes(yintercept = labels, colour = colour), size = 7) +
    geom_rect( aes(ymin = as.numeric(Groups)-0.5, ymax= as.numeric(Groups)+0.5,NULL, NULL),
               xmin = -Inf, xmax = Inf,colour = NA,
               fill = rev(colorvector))+
    geom_text(aes(x = 0.1, label = ifelse(level==3, paste("  ", as.character(Groups)), 
                                        as.character(Groups)), 
                  hjust = 0)) +
    
    # geom_text(aes(x = 1, label = ifelse(as.numeric(Q)==0, NA, round(as.numeric(Q), digits = 2))))+                #round(as.numeric(Q), digits = 2))) +
    # 
    # geom_text(aes(x = 1.5, label = ifelse(round(as.numeric(pval.Q), digits = 3)==0, 
    #                                       "<0.001", round(as.numeric(pval.Q), digits = 3)))) +
    # geom_text(aes(x = 2, label = ifelse(as.numeric(I2)==0, NA, round(as.numeric(I2), digits = 2)))) +
    # geom_text(aes(x = 2.4, label = ifelse(as.numeric(tau2)==0, NA, round(as.numeric(tau2), digits = 2)))) +
    geom_text(aes(x = 1.35 , label = k)) +
    #scale_colour_identity() +
    #theme_void() + 
    xlim(c(0,1.5))+
    theme_MicrobeR() +
    theme(panel.border = element_blank(), axis.line = element_line()) +
    # scale_x_continuous(position = "top", breaks = c(1,1.5,2,2.4,2.7),
    #                      labels = c("Q", "pval Q", "I2", "tau2", "k"))+
    theme(axis.text.x=element_text(color = "white"),
          axis.title.x = element_text(color = "white"),
          axis.text.y = element_blank(),
          axis.ticks = element_blank(),
          axis.title.y = element_blank(),
          axis.line = element_blank(),
          legend.position = "none",
          plot.title = element_text(hjust=0, size = 12, face="bold"),
          plot.margin=unit(c(0.5,0,0.5,0), "cm"))
  return(p1_sub_DATA)
}


logexp_draw_sub(table_sub_ACE,n=21,name="ACE", 
                colorvector = c("white", "white", 
                                "grey90","grey90","grey90","grey90",
                                "white", "white", "white","white",
                                "grey90","grey90","grey90","grey90","grey90",
                                "white","white","white", "white", "white", "white"),
                filepath = "output_forestplot/logexp_sub_ACE.png")

exp_p2_sub <- function(table_sub_DATA, colorvector, name = ""){#}, list = list, name = "me"){
  # table_sub_DATA <- table_sub_Chao1
  #tROM <- left_join(list, total_ROM, by = "id")
  # colorvector = c("white", "white", 
  #                 "grey90","grey90","grey90","grey90",
  #                 "white", "white", "white",
  #                 "grey90","grey90","grey90","grey90",
  #                 "white","white","white", "white", "white", "white")
  table_sub_DATA$signif <- factor(table_sub_DATA$signif, levels = c("e", "c", "ns", NA), order=TRUE)  
  table_sub_DATA$stat <- factor(table_sub_DATA$stat, levels = c("stat", "na", NA), order=TRUE)
  p2_sub_DATA <- 
    ggplot(data = table_sub_DATA, aes(x=exp(TE), y=Groups, color = signif, shape = stat))  +
    ggtitle(name)+
    geom_rect( aes(ymin = as.numeric(Groups)-0.5, ymax= as.numeric(Groups)+0.5,NULL, NULL),
               xmin = -Inf, xmax = Inf,colour = NA,
               fill = rev(colorvector))+
    geom_vline(xintercept = 1, linetype="dashed", color="grey50") +
    geom_errorbarh(aes(xmin=exp(as.numeric(lower)), xmax=exp(as.numeric(upper)), color = signif), height=0) +
    scale_color_manual(values=c( "e"="pink3","c"="skyblue3","ns"="gray44","NA"= "gray"))  +
    geom_point(aes(shape = stat, size = stat))+#aes(size = stat))#+
    scale_shape_manual(values=c(18, 5, "stat" = 16, "ns" = 1))+
    scale_size_manual(values=c(2.5, 2, "stat" = 1, "ns" = 0.5))+
    scale_y_discrete(limits=rev)+
    theme_MicrobeR() +
    #scale_x_continuous(position = "top")+
    theme(panel.border = element_blank(), axis.line = element_line()) +
    theme(axis.text.x=element_text(color = "black"),
          axis.title.x = element_text(color = "black"),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.y = element_blank(),
          axis.line.y = element_blank(),
          legend.position = "none",
          plot.title = element_text(hjust=0.5, size = 12, face="bold"),
          plot.margin=unit(c(0.5,0,0.5,0), "cm"))
  return(p2_sub_DATA)
  
}


exp_p3_sub<- function(table_sub_DATA, colorvector){
  #  list$id <- factor(list$id, levels = rev(list$id), order=TRUE)
  # table_sub_DATA <- table_sub_Chao1
  p3_sub_DATA <- 
    ggplot(data=table_sub_DATA, aes(y=rev(Groups)))+
    ggtitle(bquote(bold("       ROM [95%CI]          Z-pval       I"^"2"*"(%)       Q-pval")))+
    #geom_hline(aes(yintercept = labels, colour = colour), size = 7) +
    geom_rect( aes(ymin = as.numeric(Groups)-0.5, ymax= as.numeric(Groups)+0.5,NULL, NULL),
               xmin = -Inf, xmax = Inf,colour = NA,
               fill = rev(colorvector))+
    geom_text(aes(x = 0.2, label = round(exp(TE), digits =2)),hjust = 1) +
    geom_text(aes(x = 0.2, label = ifelse(is.na(lower),"", 
                                          paste( " [",format(round(exp(lower), digits =2), nsmall=2), ", ",
                                                 round(exp(upper), digits =2), "]", sep=""))), hjust = 0) +                #round(as.numeric(Q), digits = 2))) +
    
    geom_text(aes(x = 1, label = ifelse(format(round(as.numeric(pval), digits = 3), nsmall=3)==0, 
                                          "<0.001", format(round(as.numeric(pval), digits = 3), nsmall=3)))) +
    geom_text(aes(x = 1.4, label = ifelse(as.numeric(I2)==0, NA, format((round(as.numeric(I2), digits = 4))*100, nsmall=2)))) +
    
    geom_text(aes(x = 1.8, label = ifelse(is.na(signif),
                                          ifelse(round(as.numeric(pval.Q), digits = 3)==0, 
                                                 "<0.001", format(round(as.numeric(pval.Q), digits = 3), nsmall=3)),
                                                 NA))) +
    
    theme_MicrobeR() +
    theme(panel.border = element_blank(), axis.line = element_line()) +
    # scale_x_continuous(position = "top", breaks = c(1,1.5,2,2.4,2.7),
    #                      labels = c("Q", "pval Q", "I2", "tau2", "k"))+
    xlim(c(0,1.9))+
    theme(axis.text.x=element_text(color = "white"),
          axis.title.x = element_text(color = "white"),
          axis.text.y = element_blank(),
          axis.ticks = element_blank(),
          axis.title.y = element_blank(),
          axis.line = element_blank(),
          legend.position = "none",
          plot.title = element_text(hjust=0, size = 12, face="bold"),
          plot.margin=unit(c(0.3,0,0.5,0), "cm"))
  return(p3_sub_DATA)
}


exp_draw_sub <- function(table_sub_DATA, name="",  n=number, colorvector, filepath = ".png"){
  p1_sub_DATA <- exp_p1_sub(table_sub_DATA, colorvector)
  p2_sub_DATA <- exp_p2_sub(table_sub_DATA, colorvector, name = "")
  p3_sub_DATA <- exp_p3_sub(table_sub_DATA, colorvector)
  
  png(filepath, width = 12, height = (n+3)/5, units = "in", res = 300)
  grid.arrange(p1_sub_DATA, p2_sub_DATA, p3_sub_DATA, 
               ncol=3, 
               widths= c(30,10,20)#, heights = c(5,5,5,5,5,5,5)
  )
  dev.off()
}

logexp_draw_sub <- function(table_sub_DATA, name="",  n=number, colorvector, filepath = ".png"){
  p1_sub_DATA <- exp_p1_sub(table_sub_DATA, colorvector)
  p2_sub_DATA <- p2_sub(table_sub_DATA, colorvector, name = name)
  p3_sub_DATA <- exp_p3_sub(table_sub_DATA, colorvector)
  
  png(filepath, width = 8, height = (n+3)/4.5, units = "in", res = 300)
  grid.arrange(p1_sub_DATA, p2_sub_DATA, p3_sub_DATA, 
               ncol=3, 
               widths= c(10,10,20)#, heights = c(5,5,5,5,5,5,5)
  )
  dev.off()
}



expROM_table_sub_Chao1 <- subgroup_expROM(a_Chao1, metadata)
expROM_table_sub_Shannon <- subgroup_expROM(a_Shannon, metadata)
expROM_table_sub_Simpson <- subgroup_expROM(a_Simpson, metadata)
expROM_table_sub_ACE <- subgroup_expROM(a_ACE, metadata)
expROM_table_sub_OTU <- subgroup_expROM(a_OTU, metadata)
expROM_table_sub_Sobs <- subgroup_expROM(a_Sobs, metadata)


exp_draw_sub(expROM_table_sub_Chao1, n=19,
             colorvector = c("white", "white", 
                             "grey90","grey90","grey90","grey90",
                             "white", "white", "white",
                             "grey90","grey90","grey90","grey90",
                             "white","white","white", "white", "white", "white"),
             filepath = "output_forestplot/expROM_sub_Chao1.png")
exp_draw_sub(expROM_table_sub_Shannon, name = "Shannon",n=20,
             colorvector = c("white", "white", 
                             "grey90","grey90","grey90","grey90","grey90",
                             "white", "white", "white",
                             "grey90","grey90","grey90","grey90",
                             "white","white","white", "white", "white", "white"),
             filepath = "output_forestplot/expROM_sub_Shannon.png")
exp_draw_sub(expROM_table_sub_Simpson, n=21,
             colorvector = c("white", "white", 
                             "grey90","grey90","grey90","grey90",
                             "white", "white", "white","white",
                             "grey90","grey90","grey90","grey90","grey90",
                             "white","white","white", "white", "white", "white"),
             filepath = "output_forestplot/expROM_sub_Simpson.png")
exp_draw_sub(expROM_table_sub_ACE,n=21,
             colorvector = c("white", "white", 
                             "grey90","grey90","grey90","grey90",
                             "white", "white", "white","white",
                             "grey90","grey90","grey90","grey90","grey90",
                             "white","white","white", "white", "white", "white"),
             filepath = "output_forestplot/expROM_sub_ACE.png")
exp_draw_sub(expROM_table_sub_OTU, n=16,
             colorvector = c("white", "white", 
                             "grey90","grey90","grey90",
                             "white", "white", "white",
                             "grey90","grey90","grey90",
                             "white","white","white", "white", "white"),
             filepath = "output_forestplot/expROM_sub_OTU.png")
exp_draw_sub(expROM_table_sub_Sobs, n=14,
             colorvector = c("white", "white", 
                             "grey90","grey90","grey90",
                             "white", "white",  "white",
                             "grey90","grey90","grey90",
                             "white","white","white"),
             filepath = "output_forestplot/expROM_sub_Sobs.png")






########### +dkdkdkd########
multi <- function(a_DATA, name = "DATA"){
  # name = "DATA"
  # a_DATA <- a_Chao1
  ROM_DATA <- metacont(i_num, i_mean, i_sd, c_num, c_mean, c_sd, 
                       sm = "ROM", id, #subgroup = disease,
                       #backtransf # TRUE (default), results will be presented as ratio of means; otherwise log ratio of means will be shown.
                       data = a_DATA)
  df_DATA <- as.data.frame(ROM_DATA)
  table_DATA <- data.frame(df_DATA$studlab)    
  table_DATA$TE <- df_DATA$TE
  table_DATA$lower <- df_DATA$lower
  table_DATA$upper <- df_DATA$upper
  table_DATA$pval <- df_DATA$pval
  table_DATA$w.common <- df_DATA$w.common
  table_DATA$w.random <- df_DATA$w.random
  
  fixed <- data.frame(c("Fixed.Effect", ROM_DATA[c("TE.fixed", "lower.fixed", "upper.fixed", "pval.fixed")], NA,NA))
  random <- data.frame(c("Random.Effect", ROM_DATA[c("TE.random", "lower.random", "upper.random", "pval.random")], NA,NA))
  colnames(fixed) <-colnames(table_DATA)
  colnames(random) <- colnames(table_DATA)
  table_DATA <- rbind(table_DATA, fixed, random)
  
  signif <- NULL
  stat <- NULL
  
  for (i in 1:nrow(table_DATA)) {
    if (is.na(table_DATA$pval[i])){
      stat <- rbind(stat,"na")
      signif <- rbind(signif, ifelse(table_DATA$TE[i]<0, "e", "c"))
    } else if (table_DATA$pval[i]>=0.05) {
      stat <- rbind(stat, "stat")
      signif <- rbind(signif, "ns")
    } else {
      stat <- rbind(stat, "stat")
      signif <- rbind(signif, ifelse(table_DATA$TE[i]<0, "e", "c"))
    }
  }
  
  #signif <- factor(signif, levels = c("e", "c", "ns", NA), order=TRUE)  
  #stat <- factor(stat, levels = c("stat", "na", NA), order=TRUE)  
  table_DATA <- cbind(table_DATA, stat, signif)
  
  # coln <- paste(name, colnames(table_DATA), sep = ".")
  # coln[1]<- c("id")
  # colnames(table_DATA) <- c(coln)
  colnames(table_DATA)[1] <- c("id")
  
  # table_DATA$id <- factor(table_DATA$id, levels = table_DATA$id, order=TRUE)
  
  return(table_DATA)
}

