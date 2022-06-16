###### functions #####
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
    fixed$signif <- "sum.ns"
  } else {
    fixed$signif <- ifelse(fixed$TE<0, "sum.e", "sum.c")
  }
  
  if (random$pval>=0.05) {
    random$signif <- "sum.ns"
  } else {
    random$signif <- ifelse(random$TE<0, "sum.e", "sum.c")
  }
  
  signif <- NULL
  
  for (i in 1:nrow(m)) {
    if (is.na(m$pval[i])){
      signif <- rbind(signif, ifelse(m$TE[i]<0, "na.e", "na.c"))
    } else if (m$pval[i]>=0.05) {
      signif <- rbind(signif, "stat.ns")
    } else {
      signif <- rbind(signif, ifelse(m$TE[i]<0, "stat.e", "stat.c"))
    }
  }

 
  m$signif <- signif
  m <- rbind(m, fixed, random)
  m$signif <- factor(m$signif, levels = c("na.e", "na.c", "na.ns",
                                      "stat.e", "stat.c", "stat.ns",
                                      "sum.e", "sum.c", "sum.ns"), order=TRUE)
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

tROM_alpha_ACE<- left_join(tROM_alpha_ACE, color_factor, by = "signif", copy = TRUE)


check_plot <- function(tROM){
  ggplot(data = tROM, aes(x=TE, y=studlab, color = signif, shape = signif)) +
    geom_vline(xintercept = 0, linetype="dashed", color="grey50") +
    geom_errorbarh(aes(xmin=lower, xmax=upper), height=0 ) +
    theme_MicrobeR() +
    geom_point()+
    scale_y_discrete(limits=rev)+
    theme(panel.border = element_blank(), axis.line = element_line()) +
    theme(axis.text.x=element_text(angle=45, hjust=1))
}




