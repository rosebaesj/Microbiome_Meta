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
  ggplot(data = tROM, aes(x=TE, y=studlab, color = signif)) +
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




