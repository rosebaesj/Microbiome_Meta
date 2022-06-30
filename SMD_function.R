
subgroup_SMD <- function(a_DATA, metadata){
  
  # a_DATA <- a_Chao1
  g_DATA <- left_join(a_DATA, metadata)
  
  g_list <- c("Q","df.Q", "pval.Q","H", "I2", "tau2", 
              "TE.SMD", "lower.SMD", "upper.SMD", "pval.SMD")
  
  total_DATA <- metacont(i_num, i_mean, i_sd, c_num, c_mean, c_sd, 
                         sm = "SMD", id, #subgroup = disease_group,
                         #backtransf # TRUE (default), results will be presented as ratio of means; otherwise log ratio of means will be shown.
                         data = g_DATA)
  
  Fixed.Effect <- data.frame(c("Fixed.Effect",
                               total_DATA[c("Q","df.Q", "pval.Q","H", "I2", "tau2", "k",
                                            "TE.fixed", "lower.fixed", "upper.fixed", "pval.fixed")], 1))
  Random.Effect <- data.frame(c("Random.Effect",
                                total_DATA[c("Q","df.Q", "pval.Q","H", "I2", "tau2", "k",
                                             "TE.random", "lower.random", "upper.random", "pval.random")], 1))
  
  disease_DATA <- metacont(i_num, i_mean, i_sd, c_num, c_mean, c_sd, 
                           sm = "SMD", id, subgroup = disease_group,
                           #backtransf # TRUE (default), results will be presented as ratio of means; otherwise log ratio of means will be shown.
                           data = g_DATA)
  
  Disease <- data.frame(c("Disease", 
                          disease_DATA[c("Q.b.random","df.Q.b", "pval.Q.b.random")],
                          NA, NA, NA, NA,NA, NA, NA, NA, 2))
  
  Disease_subgroup <- data.frame(c(disease_DATA[c("bylevs", "Q.w")], NA, 
                                   disease_DATA[c("pval.Q.w")], NA,
                                   disease_DATA[c("I2.w", "tau2.w","k.w",
                                                  "TE.fixed.w", "lower.fixed.w", "upper.fixed.w", "pval.fixed.w")], 3))
  
  treatment_DATA <- metacont(i_num, i_mean, i_sd, c_num, c_mean, c_sd, 
                             sm = "SMD", id, subgroup = i_type,
                             #backtransf # TRUE (default), results will be presented as ratio of means; otherwise log ratio of means will be shown.
                             data = g_DATA)
  
  Treatment <- data.frame(c("Treatment", treatment_DATA[c("Q.b.random","df.Q.b", "pval.Q.b.random")],NA, NA, NA, NA,
                            NA, NA, NA, NA, 2))
  
  Treatment_subgroup <- data.frame(c(treatment_DATA[c("bylevs","Q.w")], NA, 
                                     treatment_DATA[c("pval.Q.w")], NA,
                                     treatment_DATA[c("I2.w", "tau2.w","k.w",
                                                      "TE.fixed.w", "lower.fixed.w", "upper.fixed.w", "pval.fixed.w")], 3))
  species_DATA <- metacont(i_num, i_mean, i_sd, c_num, c_mean, c_sd, 
                           sm = "SMD", id, subgroup = species,
                           #backtransf # TRUE (default), results will be presented as ratio of means; otherwise log ratio of means will be shown.
                           data = g_DATA)
  
  Species <- data.frame(c("Species",species_DATA[c("Q.b.random","df.Q.b", "pval.Q.b.random")],NA, NA, NA, NA,
                          NA, NA, NA, NA, 2))
  
  Species_subgroup <- data.frame(c(species_DATA[c("bylevs","Q.w")], NA, 
                                   species_DATA[c("pval.Q.w")], NA,
                                   species_DATA[c("I2.w", "tau2.w","k.w",
                                                  "TE.fixed.w", "lower.fixed.w", "upper.fixed.w", "pval.fixed.w")], 3))
  
  duration_DATA <- metacont(i_num, i_mean, i_sd, c_num, c_mean, c_sd, 
                            sm = "SMD", id, subgroup = weeks,
                            #backtransf # TRUE (default), results will be presented as ratio of means; otherwise log ratio of means will be shown.
                            data = g_DATA)
  
  Duration <- data.frame(c("Duration", duration_DATA[c("Q.b.random","df.Q.b", "pval.Q.b.random")],NA, NA, NA, 
                           NA, NA, NA, NA, NA, 2))
  
  Duration_subgroup <- data.frame(c(duration_DATA[c("bylevs","Q.w")], NA, 
                                    duration_DATA[c("pval.Q.w")], NA,
                                    duration_DATA[c("I2.w", "tau2.w","k.w",
                                                    "TE.fixed.w", "lower.fixed.w", "upper.fixed.w", "pval.fixed.w")], 3))
  
  
  table_sub_DATA <- t(cbind(t(Fixed.Effect), t(Random.Effect), 
                            t(Disease), t(Disease_subgroup),
                            t(Treatment), t(Treatment_subgroup),
                            t(Species), t(Species_subgroup),
                            t(Duration), t(Duration_subgroup)
                            
  ))
  
  
  colnames(table_sub_DATA) <- c("Groups", "Q" , "df.Q" ,"pval.Q","H" , "I2" ,"tau2","k",
                                "TE", "lower", "upper", "pval" ,"level" )
  table_sub_DATA <- data.frame(table_sub_DATA)
  table_sub_DATA$Groups <- factor(table_sub_DATA$Groups, levels = table_sub_DATA$Groups, order=TRUE)
  
  
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



SMD_table_sub_Chao1 <- subgroup_SMD(a_Chao1, metadata)
SMD_table_sub_Shannon <- subgroup_SMD(a_Shannon, metadata)
SMD_table_sub_Simpson <- subgroup_SMD(a_Simpson, metadata)
SMD_table_sub_ACE <- subgroup_SMD(a_ACE, metadata)
SMD_table_sub_OTU <- subgroup_SMD(a_OTU, metadata)
SMD_table_sub_Sobs <- subgroup_SMD(a_Sobs, metadata)


draw_sub(SMD_table_sub_Chao1, n=19,
         colorvector = c("white", "white", 
                         "grey90","grey90","grey90","grey90",
                         "white", "white", "white",
                         "grey90","grey90","grey90","grey90",
                         "white","white","white", "white", "white", "white"),
         filepath = "output_forestplot/SMD_sub_Chao1.png")
draw_sub(SMD_table_sub_Shannon, name = "Shannon",n=20,
         colorvector = c("white", "white", 
                         "grey90","grey90","grey90","grey90","grey90",
                         "white", "white", "white",
                         "grey90","grey90","grey90","grey90",
                         "white","white","white", "white", "white", "white"),
         filepath = "output_forestplot/SMD_sub_Shannon.png")
draw_sub(SMD_table_sub_Simpson, n=21,
         colorvector = c("white", "white", 
                         "grey90","grey90","grey90","grey90",
                         "white", "white", "white","white",
                         "grey90","grey90","grey90","grey90","grey90",
                         "white","white","white", "white", "white", "white"),
         filepath = "output_forestplot/SMD_sub_Simpson.png")
draw_sub(SMD_table_sub_ACE,n=21,
         colorvector = c("white", "white", 
                         "grey90","grey90","grey90","grey90",
                         "white", "white", "white","white",
                         "grey90","grey90","grey90","grey90","grey90",
                         "white","white","white", "white", "white", "white"),
         filepath = "output_forestplot/SMD_sub_ACE.png")
draw_sub(SMD_table_sub_OTU, n=16,
         colorvector = c("white", "white", 
                         "grey90","grey90","grey90",
                         "white", "white", "white",
                         "grey90","grey90","grey90",
                         "white","white","white", "white", "white"),
         filepath = "output_forestplot/SMD_sub_OTU.png")
draw_sub(SMD_table_sub_Sobs, n=14,
         colorvector = c("white", "white", 
                         "grey90","grey90","grey90",
                         "white", "white",  "white",
                         "grey90","grey90","grey90",
                         "white","white","white"),
         filepath = "output_forestplot/SMD_sub_Sobs.png")

