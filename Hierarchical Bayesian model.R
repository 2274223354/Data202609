##===============================================================================
#####FD###################################################################
##===============================================================================
##### Load data ###################################################################

library(data.table)
library(dplyr)
library(R2jags)

C_datatotal_all <- fread("/data/Data.csv" )
colnames( C_datatotal_all  )

C_datatotal_all <- C_datatotal_all %>%
  filter(Carbon_Mg_ha > 0)
colnames(C_datatotal_all) <- gsub("_mean$", "", colnames(C_datatotal_all))
C_datatotal_all <- C_datatotal_all %>%
  filter(!is.na(CWD_stand_age)) %>%
  filter(!is.na(MAT_stand_age)) %>%
  filter(!is.na(MAP_stand_age)) 


C_bin <- read.csv("/data/Class_CWD_9bins.csv")
C_bin <- C_bin[, c("PLT_CN", "PDSI_ATA_bin")]
C_datatotal_all$PLT_CN <- as.character(C_datatotal_all$PLT_CN)
C_bin$PLT_CN            <- as.character(C_bin$PLT_CN)
C_datatotal_all <- merge(C_datatotal_all, C_bin, by = "PLT_CN")
C_datatotal_all$CWD1=C_datatotal_all$CWD_stand_age
# 标准化
C_datatotal_all <- C_datatotal_all %>%
  mutate(
    Carbon_Mg_ha = log(Carbon_Mg_ha+1),
    shannon_dbhR = scale(H_D_5),
           
           
    CWD          = scale(CWD_stand_age),  
    TEM1         = scale(MAT_stand_age),
    PRE1         = scale(MAP_stand_age),
    FDis         = scale(FDis),
    CWM          = scale(PC1.CWM),
    pH           = scale(nitrogen),
    bdod         = scale(bdod),
    STDAGE       = scale(STDAGE)
  )

data_sem2 <- C_datatotal_all

##### JAGS #########################################################


data_all <- data_sem2

y_CV      <- as.numeric(data_all$shannon_dbhR)
y_AGB     <- as.numeric(data_all$Carbon_Mg_ha)
trait     <- as.numeric(data_all$FDis)
trait_CWM <- as.numeric(data_all$CWM)
CWD       <- as.numeric(data_all$CWD)       
TEM1      <- as.numeric(data_all$TEM1)
PRE1      <- as.numeric(data_all$PRE1)
pH        <- as.numeric(data_all$pH)
bdod      <- as.numeric(data_all$bdod)
STDAGE    <- as.numeric(data_all$STDAGE)

bin   <- as.numeric(factor(data_all$PDSI_ATA_bin))
N     <- nrow(data_all)
N_bin <- max(bin)

cat("N =", N, "| N_bin =", N_bin, "\n")


FD_bin     <- tapply(trait,     bin, mean, na.rm = TRUE)
FI_bin     <- tapply(trait_CWM, bin, mean, na.rm = TRUE)
CWD_bin    <- tapply(CWD,       bin, mean, na.rm = TRUE)  
TEM1_bin   <- tapply(TEM1,      bin, mean, na.rm = TRUE)
PRE1_bin   <- tapply(PRE1,      bin, mean, na.rm = TRUE)
pH_bin     <- tapply(pH,        bin, mean, na.rm = TRUE)
bdod_bin   <- tapply(bdod,      bin, mean, na.rm = TRUE)
STDAGE_bin <- tapply(STDAGE,    bin, mean, na.rm = TRUE)



print(sapply(
  list(FD=FD_bin, FI=FI_bin, CWD=CWD_bin,
       TEM1=TEM1_bin, PRE1=PRE1_bin,
       pH=pH_bin, bdod=bdod_bin, STDAGE=STDAGE_bin),
  function(x) sum(is.na(x))
))


jags_data <- list(
  N            = N,
  N_bin        = N_bin,
  DBH_CV       = y_CV,
  Carbon_Mg_ha = y_AGB,
  trait        = trait,
  trait_CWM    = trait_CWM,
  CWD          = CWD,                        
  TEM1         = TEM1,
  PRE1         = PRE1,
  pH           = pH,
  bdod         = bdod,
  STDAGE       = STDAGE,
  bin          = bin,
  FD_bin       = as.numeric(FD_bin),
  FI_bin       = as.numeric(FI_bin),
  CWD_bin      = as.numeric(CWD_bin),        
  TEM1_bin     = as.numeric(TEM1_bin),
  PRE1_bin     = as.numeric(PRE1_bin),
  pH_bin       = as.numeric(pH_bin),
  bdod_bin     = as.numeric(bdod_bin),
  STDAGE_bin   = as.numeric(STDAGE_bin)
)



##### parameters##########################################


params <- c(
  "delta1", "delta2",
  "beta3",  "beta4",  "beta5",
  "gamma0_CV",
  "mu_alpha_AGB",
  "gamma_CV",
  "gamma_CV_bin",
  "gamma_AGB",
  "beta5_0",
  "beta5_1", "beta5_2",
  "beta5_3", "beta5_4",
  "beta5_5", "beta5_6",
  "beta5_7", "beta5_8",
  "sigma_CV",       "sigma_AGB",
  "sigma_alpha_CV", "sigma_alpha_AGB",
  "sigma_delta1",   "sigma_delta2",
  "sigma_beta3",    "sigma_beta4",    "sigma_beta5"
)



##### Run JAGS #################################################################


fit_sem_jags <- jags(
  model.file         = "/data/bayes_sem_model.txt",
  data               = jags_data,
  parameters.to.save = params,
  n.chains           = 3,
  n.iter             = 10000,
  n.burnin           = 2000,
  n.thin             = 10
)


##### Save #############################################################


post        <- as.data.frame(fit_sem_jags$BUGSoutput$sims.list)
summary_tab <- as.data.frame(fit_sem_jags$BUGSoutput$summary)

# 收敛诊断
rhats     <- fit_sem_jags$BUGSoutput$summary[, "Rhat"]
bad_rhats <- rhats[rhats > 1.1]
if (length(bad_rhats) > 0) {
  cat("（Rhat > 1.1）:\n"); print(bad_rhats)
} else {
  cat("Rhat < 1.1）\n")
}


saveRDS(fit_sem_jags,
        "/data/beiyesi_CWD9_DBH5_FD.rds")
write.csv(summary_tab,
          "/data/Beiyrsi_Summary_American_CWD9_DBH5_FD.csv",
          row.names = TRUE)
saveRDS(post,
        "/data/post_American_CWD9_DBH5_FD.rds")



##### Extract #####################################################



get_effect_by_bin <- function(post, param_prefix, path_label) {
  param_cols  <- grep(paste0("^", param_prefix, "\\."), names(post), value = TRUE)
  summary_mat <- t(apply(
    post[, param_cols, drop = FALSE], 2,
    function(x) c(mean  = mean(x),
                  lower = quantile(x, 0.025),
                  upper = quantile(x, 0.975))
  ))
  df           <- as.data.frame(summary_mat)
  colnames(df) <- c("mean", "lower", "upper")
  df$bin       <- as.numeric(sub(paste0("^", param_prefix, "\\."), "", param_cols))
  df$effect_type <- ifelse(
    df$lower > 0 & df$upper > 0, "positive significant",
    ifelse(df$lower < 0 & df$upper < 0, "negative significant",
           ifelse(df$mean > 0, "positive non-significant", "negative non-significant"))
  )
  df$path <- path_label
  return(df)
}


get_indirect_effect_by_bin <- function(post, delta_prefix, beta_prefix, path_label) {
  delta_cols  <- sort(grep(paste0("^", delta_prefix, "\\."), names(post), value = TRUE))
  beta_cols   <- sort(grep(paste0("^", beta_prefix,  "\\."), names(post), value = TRUE))
  sorted_bins <- as.numeric(sub(paste0("^", delta_prefix, "\\."), "", delta_cols))
  indirect    <- post[, delta_cols] * post[, beta_cols]
  summary_mat <- t(apply(indirect, 2,
                         function(x) c(mean  = mean(x),
                                       lower = quantile(x, 0.025),
                                       upper = quantile(x, 0.975))
  ))
  df           <- as.data.frame(summary_mat)
  colnames(df) <- c("mean", "lower", "upper")
  df$bin       <- sorted_bins
  df$effect_type <- ifelse(
    df$lower > 0 & df$upper > 0, "positive significant",
    ifelse(df$lower < 0 & df$upper < 0, "negative significant",
           ifelse(df$mean > 0, "positive non-significant", "negative non-significant"))
  )
  df$effect_type <- factor(df$effect_type, levels = c(
    "negative significant", "negative non-significant",
    "positive non-significant", "positive significant"))
  df$path <- path_label
  return(df)
}


all_effects <- bind_rows(
  get_effect_by_bin(post, "beta3",  "FD → AGB"),
  get_effect_by_bin(post, "beta4",  "CWM → AGB"),
  get_effect_by_bin(post, "beta5",  "DBHCV → AGB"),
  get_effect_by_bin(post, "delta1", "FD → DBHCV"),
  get_effect_by_bin(post, "delta2", "CWM → DBHCV"),
  get_indirect_effect_by_bin(post, "delta1", "beta5", "FD → DBHCV → AGB"),
  get_indirect_effect_by_bin(post, "delta2", "beta5", "CWM → DBHCV → AGB")
)


data_sem2_org     <- data_sem2
data_sem2_org$bin <- as.numeric(factor(data_sem2_org$PDSI_ATA_bin))

result <- data_sem2_org %>%
  group_by(bin) %>%
  dplyr::summarise(
    mean_MAT    = mean(mean_annual_temp,   na.rm = TRUE),
    mean_PRE    = mean(mean_annual_precip, na.rm = TRUE),
    mean_CWD    = mean(CWD1,           na.rm = TRUE),  
    mean_STDAGE = mean(STDAGE,             na.rm = TRUE),
    .groups = "drop"
  )

all_effects1 <- merge(all_effects, result, by = "bin")

write.csv(all_effects1,
          "/data/posterior_American_CWD9_DBH5_FD.csv",
          row.names = FALSE)


##### Proportional stacked bar chart #################################################################


library(scales)

get_effect_prop <- function(post, param_prefix, path_label) {
  param_cols <- grep(paste0("^", param_prefix, "\\."), names(post), value = TRUE)
  df <- t(apply(post[, param_cols, drop = FALSE], 2,
                function(x) c(mean  = mean(x),
                              lower = quantile(x, 0.025),
                              upper = quantile(x, 0.975))
  )) |> as.data.frame()
  colnames(df) <- c("mean", "lower", "upper")
  df$effect_type <- ifelse(
    df$lower > 0 & df$upper > 0, "positive significant",
    ifelse(df$lower < 0 & df$upper < 0, "negative significant",
           ifelse(df$mean > 0, "positive non-significant", "negative non-significant"))
  )
  dplyr::count(df, effect_type) %>%
    mutate(percent = n / sum(n) * 100, path = path_label)
}

get_indirect_prop <- function(post, delta_prefix, beta_prefix, path_label) {
  delta_cols <- sort(grep(paste0("^", delta_prefix, "\\."), names(post), value = TRUE))
  beta_cols  <- sort(grep(paste0("^", beta_prefix,  "\\."), names(post), value = TRUE))
  indirect   <- post[, delta_cols] * post[, beta_cols]
  df <- t(apply(indirect, 2,
                function(x) c(mean  = mean(x),
                              lower = quantile(x, 0.025),
                              upper = quantile(x, 0.975))
  )) |> as.data.frame()
  colnames(df) <- c("mean", "lower", "upper")
  df$effect_type <- ifelse(
    df$lower > 0 & df$upper > 0, "positive significant",
    ifelse(df$lower < 0 & df$upper < 0, "negative significant",
           ifelse(df$mean > 0, "positive non-significant", "negative non-significant"))
  )
  dplyr::count(df, effect_type) %>%
    mutate(percent = n / sum(n) * 100, path = path_label)
}

prop <- bind_rows(
  get_effect_prop(post, "beta3",  "FD → AGB"),
  get_effect_prop(post, "beta4",  "CWM → AGB"),
  get_effect_prop(post, "beta5",  "DBHCV → AGB"),
  get_effect_prop(post, "delta1", "FD → DBHCV"),
  get_effect_prop(post, "delta2", "CWM → DBHCV"),
  get_indirect_prop(post, "delta1", "beta5", "FD → DBHCV → AGB"),
  get_indirect_prop(post, "delta2", "beta5", "CWM → DBHCV → AGB")
)

write.csv(prop,
          "/data/bar_American_CWD9_DBH5_FD.csv",
          row.names = FALSE)


prop_plot <- prop %>%
  filter(path %in% c("FD → AGB", "FD → DBHCV → AGB",
                     "CWM → AGB", "CWM → DBHCV → AGB")) %>%
  mutate(
    path = factor(path, levels = c(
      "FD → AGB", "FD → DBHCV → AGB",
      "CWM → AGB", "CWM → DBHCV → AGB"
    )),
    effect_type = factor(effect_type, levels = c(
      "negative significant", "negative non-significant",
      "positive non-significant", "positive significant"
    ))
  )

fig_height <- ggplot(prop_plot,
                     aes(x = percent, y = path, fill = effect_type)) +
  geom_bar(stat = "identity", width = 0.65) +
  scale_x_continuous(labels = percent_format(scale = 1), expand = c(0, 0)) +
  scale_fill_manual(
    values = c(
      "negative significant"     = "#92C5DE",
      "negative non-significant" = "#dbebfa",
      "positive non-significant" = "#f9ebdf",
      "positive significant"     = "#F4A582"
    ),
    drop = FALSE
  ) +
  labs(x = "Proportion of bins (%)", y = NULL, fill = "Effect type") +
  ggtitle("Structural diversity → Carbon") +
  theme_classic(base_size = 13) +
  theme(axis.text.y  = element_text(size = 12),
        legend.position = "right")

fig_height
##===============================================================================
#####LMA###################################################################
##===============================================================================
##### Load data ###################################################################

library(data.table)
library(dplyr)
library(R2jags)

C_datatotal_all <- fread("/data/Data.csv" )
colnames( C_datatotal_all  )

C_datatotal_all <- C_datatotal_all %>%
  dplyr::select(-V1) %>%
  filter(Carbon_Mg_ha > 0)
colnames(C_datatotal_all) <- gsub("_mean$", "", colnames(C_datatotal_all))
C_datatotal_all <- C_datatotal_all %>%
  filter(!is.na(CWD_stand_age)) %>%
  filter(!is.na(MAT_stand_age)) %>%
  filter(!is.na(MAP_stand_age)) 

# 合并 bin 分组
C_bin <- read.csv("/data/Class_CWD_9bins")
C_bin <- C_bin[, c("PLT_CN", "PDSI_ATA_bin")]
C_datatotal_all$PLT_CN <- as.character(C_datatotal_all$PLT_CN)
C_bin$PLT_CN            <- as.character(C_bin$PLT_CN)
C_datatotal_all <- merge(C_datatotal_all, C_bin, by = "PLT_CN")
C_datatotal_all$CWD1=C_datatotal_all$CWD_stand_age
# 标准化
C_datatotal_all <- C_datatotal_all %>%
  mutate(
    Carbon_Mg_ha = log(Carbon_Mg_ha+1),
    shannon_dbhR = scale(H_D_5),
           
           
    CWD          = scale(CWD_stand_age),  
    TEM1         = scale(MAT_stand_age),
    PRE1         = scale(MAP_stand_age),
    FDis         = scale(LMA.FDis),
    CWM          = scale(LMA.CWM),
    pH           = scale(nitrogen),
    bdod         = scale(bdod),
    STDAGE       = scale(STDAGE)
  )

data_sem2 <- C_datatotal_all

##### JAGS #########################################################


data_all <- data_sem2

y_CV      <- as.numeric(data_all$shannon_dbhR)
y_AGB     <- as.numeric(data_all$Carbon_Mg_ha)
trait     <- as.numeric(data_all$FDis)
trait_CWM <- as.numeric(data_all$CWM)
CWD       <- as.numeric(data_all$CWD)       
TEM1      <- as.numeric(data_all$TEM1)
PRE1      <- as.numeric(data_all$PRE1)
pH        <- as.numeric(data_all$pH)
bdod      <- as.numeric(data_all$bdod)
STDAGE    <- as.numeric(data_all$STDAGE)

bin   <- as.numeric(factor(data_all$PDSI_ATA_bin))
N     <- nrow(data_all)
N_bin <- max(bin)

cat("N =", N, "| N_bin =", N_bin, "\n")


FD_bin     <- tapply(trait,     bin, mean, na.rm = TRUE)
FI_bin     <- tapply(trait_CWM, bin, mean, na.rm = TRUE)
CWD_bin    <- tapply(CWD,       bin, mean, na.rm = TRUE)  
TEM1_bin   <- tapply(TEM1,      bin, mean, na.rm = TRUE)
PRE1_bin   <- tapply(PRE1,      bin, mean, na.rm = TRUE)
pH_bin     <- tapply(pH,        bin, mean, na.rm = TRUE)
bdod_bin   <- tapply(bdod,      bin, mean, na.rm = TRUE)
STDAGE_bin <- tapply(STDAGE,    bin, mean, na.rm = TRUE)



print(sapply(
  list(FD=FD_bin, FI=FI_bin, CWD=CWD_bin,
       TEM1=TEM1_bin, PRE1=PRE1_bin,
       pH=pH_bin, bdod=bdod_bin, STDAGE=STDAGE_bin),
  function(x) sum(is.na(x))
))


jags_data <- list(
  N            = N,
  N_bin        = N_bin,
  DBH_CV       = y_CV,
  Carbon_Mg_ha = y_AGB,
  trait        = trait,
  trait_CWM    = trait_CWM,
  CWD          = CWD,                        
  TEM1         = TEM1,
  PRE1         = PRE1,
  pH           = pH,
  bdod         = bdod,
  STDAGE       = STDAGE,
  bin          = bin,
  FD_bin       = as.numeric(FD_bin),
  FI_bin       = as.numeric(FI_bin),
  CWD_bin      = as.numeric(CWD_bin),        
  TEM1_bin     = as.numeric(TEM1_bin),
  PRE1_bin     = as.numeric(PRE1_bin),
  pH_bin       = as.numeric(pH_bin),
  bdod_bin     = as.numeric(bdod_bin),
  STDAGE_bin   = as.numeric(STDAGE_bin)
)



##### parameters##########################################


params <- c(
  "delta1", "delta2",
  "beta3",  "beta4",  "beta5",
  "gamma0_CV",
  "mu_alpha_AGB",
  "gamma_CV",
  "gamma_CV_bin",
  "gamma_AGB",
  "beta5_0",
  "beta5_1", "beta5_2",
  "beta5_3", "beta5_4",
  "beta5_5", "beta5_6",
  "beta5_7", "beta5_8",
  "sigma_CV",       "sigma_AGB",
  "sigma_alpha_CV", "sigma_alpha_AGB",
  "sigma_delta1",   "sigma_delta2",
  "sigma_beta3",    "sigma_beta4",    "sigma_beta5"
)



##### Run JAGS #################################################################


fit_sem_jags <- jags(
  model.file         = "/data/bayes_sem_model.txt",
  data               = jags_data,
  parameters.to.save = params,
  n.chains           = 3,
  n.iter             = 10000,
  n.burnin           = 2000,
  n.thin             = 10
)



##### Save #############################################################


post        <- as.data.frame(fit_sem_jags$BUGSoutput$sims.list)
summary_tab <- as.data.frame(fit_sem_jags$BUGSoutput$summary)

# 收敛诊断
rhats     <- fit_sem_jags$BUGSoutput$summary[, "Rhat"]
bad_rhats <- rhats[rhats > 1.1]
if (length(bad_rhats) > 0) {
  cat("（Rhat > 1.1）:\n"); print(bad_rhats)
} else {
  cat("Rhat < 1.1）\n")
}


saveRDS(fit_sem_jags,
        "/data/beiyesi_CWD9_DBH5_LMA.rds")
write.csv(summary_tab,
          "/data/Beiyrsi_Summary_American_CWD9_DBH5_LMA.csv",
          row.names = TRUE)
saveRDS(post,
        "/data/post_American_CWD9_DBH5_LMA.rds")



##### Extract #####################################################



get_effect_by_bin <- function(post, param_prefix, path_label) {
  param_cols  <- grep(paste0("^", param_prefix, "\\."), names(post), value = TRUE)
  summary_mat <- t(apply(
    post[, param_cols, drop = FALSE], 2,
    function(x) c(mean  = mean(x),
                  lower = quantile(x, 0.025),
                  upper = quantile(x, 0.975))
  ))
  df           <- as.data.frame(summary_mat)
  colnames(df) <- c("mean", "lower", "upper")
  df$bin       <- as.numeric(sub(paste0("^", param_prefix, "\\."), "", param_cols))
  df$effect_type <- ifelse(
    df$lower > 0 & df$upper > 0, "positive significant",
    ifelse(df$lower < 0 & df$upper < 0, "negative significant",
           ifelse(df$mean > 0, "positive non-significant", "negative non-significant"))
  )
  df$path <- path_label
  return(df)
}


get_indirect_effect_by_bin <- function(post, delta_prefix, beta_prefix, path_label) {
  delta_cols  <- sort(grep(paste0("^", delta_prefix, "\\."), names(post), value = TRUE))
  beta_cols   <- sort(grep(paste0("^", beta_prefix,  "\\."), names(post), value = TRUE))
  sorted_bins <- as.numeric(sub(paste0("^", delta_prefix, "\\."), "", delta_cols))
  indirect    <- post[, delta_cols] * post[, beta_cols]
  summary_mat <- t(apply(indirect, 2,
                         function(x) c(mean  = mean(x),
                                       lower = quantile(x, 0.025),
                                       upper = quantile(x, 0.975))
  ))
  df           <- as.data.frame(summary_mat)
  colnames(df) <- c("mean", "lower", "upper")
  df$bin       <- sorted_bins
  df$effect_type <- ifelse(
    df$lower > 0 & df$upper > 0, "positive significant",
    ifelse(df$lower < 0 & df$upper < 0, "negative significant",
           ifelse(df$mean > 0, "positive non-significant", "negative non-significant"))
  )
  df$effect_type <- factor(df$effect_type, levels = c(
    "negative significant", "negative non-significant",
    "positive non-significant", "positive significant"))
  df$path <- path_label
  return(df)
}


all_effects <- bind_rows(
  get_effect_by_bin(post, "beta3",  "FD → AGB"),
  get_effect_by_bin(post, "beta4",  "CWM → AGB"),
  get_effect_by_bin(post, "beta5",  "DBHCV → AGB"),
  get_effect_by_bin(post, "delta1", "FD → DBHCV"),
  get_effect_by_bin(post, "delta2", "CWM → DBHCV"),
  get_indirect_effect_by_bin(post, "delta1", "beta5", "FD → DBHCV → AGB"),
  get_indirect_effect_by_bin(post, "delta2", "beta5", "CWM → DBHCV → AGB")
)


data_sem2_org     <- data_sem2
data_sem2_org$bin <- as.numeric(factor(data_sem2_org$PDSI_ATA_bin))

result <- data_sem2_org %>%
  group_by(bin) %>%
  dplyr::summarise(
    mean_MAT    = mean(mean_annual_temp,   na.rm = TRUE),
    mean_PRE    = mean(mean_annual_precip, na.rm = TRUE),
    mean_CWD    = mean(CWD1,           na.rm = TRUE),  
    mean_STDAGE = mean(STDAGE,             na.rm = TRUE),
    .groups = "drop"
  )

all_effects1 <- merge(all_effects, result, by = "bin")

write.csv(all_effects1,
          "/data/posterior_American_CWD9_DBH5_LMA.csv",
          row.names = FALSE)



##### Proportional stacked bar chart #################################################################


library(scales)

get_effect_prop <- function(post, param_prefix, path_label) {
  param_cols <- grep(paste0("^", param_prefix, "\\."), names(post), value = TRUE)
  df <- t(apply(post[, param_cols, drop = FALSE], 2,
                function(x) c(mean  = mean(x),
                              lower = quantile(x, 0.025),
                              upper = quantile(x, 0.975))
  )) |> as.data.frame()
  colnames(df) <- c("mean", "lower", "upper")
  df$effect_type <- ifelse(
    df$lower > 0 & df$upper > 0, "positive significant",
    ifelse(df$lower < 0 & df$upper < 0, "negative significant",
           ifelse(df$mean > 0, "positive non-significant", "negative non-significant"))
  )
  dplyr::count(df, effect_type) %>%
    mutate(percent = n / sum(n) * 100, path = path_label)
}

get_indirect_prop <- function(post, delta_prefix, beta_prefix, path_label) {
  delta_cols <- sort(grep(paste0("^", delta_prefix, "\\."), names(post), value = TRUE))
  beta_cols  <- sort(grep(paste0("^", beta_prefix,  "\\."), names(post), value = TRUE))
  indirect   <- post[, delta_cols] * post[, beta_cols]
  df <- t(apply(indirect, 2,
                function(x) c(mean  = mean(x),
                              lower = quantile(x, 0.025),
                              upper = quantile(x, 0.975))
  )) |> as.data.frame()
  colnames(df) <- c("mean", "lower", "upper")
  df$effect_type <- ifelse(
    df$lower > 0 & df$upper > 0, "positive significant",
    ifelse(df$lower < 0 & df$upper < 0, "negative significant",
           ifelse(df$mean > 0, "positive non-significant", "negative non-significant"))
  )
  dplyr::count(df, effect_type) %>%
    mutate(percent = n / sum(n) * 100, path = path_label)
}

prop <- bind_rows(
  get_effect_prop(post, "beta3",  "FD → AGB"),
  get_effect_prop(post, "beta4",  "CWM → AGB"),
  get_effect_prop(post, "beta5",  "DBHCV → AGB"),
  get_effect_prop(post, "delta1", "FD → DBHCV"),
  get_effect_prop(post, "delta2", "CWM → DBHCV"),
  get_indirect_prop(post, "delta1", "beta5", "FD → DBHCV → AGB"),
  get_indirect_prop(post, "delta2", "beta5", "CWM → DBHCV → AGB")
)

write.csv(prop,
          "/data/bar_American_CWD9_DBH5_LMA.csv",
          row.names = FALSE)


prop_plot <- prop %>%
  filter(path %in% c("FD → AGB", "FD → DBHCV → AGB",
                     "CWM → AGB", "CWM → DBHCV → AGB")) %>%
  mutate(
    path = factor(path, levels = c(
      "FD → AGB", "FD → DBHCV → AGB",
      "CWM → AGB", "CWM → DBHCV → AGB"
    )),
    effect_type = factor(effect_type, levels = c(
      "negative significant", "negative non-significant",
      "positive non-significant", "positive significant"
    ))
  )

fig_height <- ggplot(prop_plot,
                     aes(x = percent, y = path, fill = effect_type)) +
  geom_bar(stat = "identity", width = 0.65) +
  scale_x_continuous(labels = percent_format(scale = 1), expand = c(0, 0)) +
  scale_fill_manual(
    values = c(
      "negative significant"     = "#92C5DE",
      "negative non-significant" = "#dbebfa",
      "positive non-significant" = "#f9ebdf",
      "positive significant"     = "#F4A582"
    ),
    drop = FALSE
  ) +
  labs(x = "Proportion of bins (%)", y = NULL, fill = "Effect type") +
  ggtitle("Structural diversity → Carbon") +
  theme_classic(base_size = 13) +
  theme(axis.text.y  = element_text(size = 12),
        legend.position = "right")

fig_height
##===============================================================================
#####LDMC###################################################################
##===============================================================================
##### Load data ###################################################################

library(data.table)
library(dplyr)
library(R2jags)

C_datatotal_all <- fread("/data/Data.csv" )
colnames( C_datatotal_all  )

C_datatotal_all <- C_datatotal_all %>%
  dplyr::select(-V1) %>%
  filter(Carbon_Mg_ha > 0)
colnames(C_datatotal_all) <- gsub("_mean$", "", colnames(C_datatotal_all))
C_datatotal_all <- C_datatotal_all %>%
  filter(!is.na(CWD_stand_age)) %>%
  filter(!is.na(MAT_stand_age)) %>%
  filter(!is.na(MAP_stand_age)) 

# 合并 bin 分组
C_bin <- read.csv("/data/Class_CWD_9bins")
C_bin <- C_bin[, c("PLT_CN", "PDSI_ATA_bin")]
C_datatotal_all$PLT_CN <- as.character(C_datatotal_all$PLT_CN)
C_bin$PLT_CN            <- as.character(C_bin$PLT_CN)
C_datatotal_all <- merge(C_datatotal_all, C_bin, by = "PLT_CN")
C_datatotal_all$CWD1=C_datatotal_all$CWD_stand_age
# 标准化
C_datatotal_all <- C_datatotal_all %>%
  mutate(
    Carbon_Mg_ha = log(Carbon_Mg_ha+1),
    shannon_dbhR = scale(H_D_5),
           
           
    CWD          = scale(CWD_stand_age),  
    TEM1         = scale(MAT_stand_age),
    PRE1         = scale(MAP_stand_age),
    FDis         = scale(LDMC.FDis),
    CWM          = scale(LDMC.CWM),
    pH           = scale(nitrogen),
    bdod         = scale(bdod),
    STDAGE       = scale(STDAGE)
  )

data_sem2 <- C_datatotal_all

##### JAGS #########################################################


data_all <- data_sem2

y_CV      <- as.numeric(data_all$shannon_dbhR)
y_AGB     <- as.numeric(data_all$Carbon_Mg_ha)
trait     <- as.numeric(data_all$FDis)
trait_CWM <- as.numeric(data_all$CWM)
CWD       <- as.numeric(data_all$CWD)       
TEM1      <- as.numeric(data_all$TEM1)
PRE1      <- as.numeric(data_all$PRE1)
pH        <- as.numeric(data_all$pH)
bdod      <- as.numeric(data_all$bdod)
STDAGE    <- as.numeric(data_all$STDAGE)

bin   <- as.numeric(factor(data_all$PDSI_ATA_bin))
N     <- nrow(data_all)
N_bin <- max(bin)

cat("N =", N, "| N_bin =", N_bin, "\n")


FD_bin     <- tapply(trait,     bin, mean, na.rm = TRUE)
FI_bin     <- tapply(trait_CWM, bin, mean, na.rm = TRUE)
CWD_bin    <- tapply(CWD,       bin, mean, na.rm = TRUE)  
TEM1_bin   <- tapply(TEM1,      bin, mean, na.rm = TRUE)
PRE1_bin   <- tapply(PRE1,      bin, mean, na.rm = TRUE)
pH_bin     <- tapply(pH,        bin, mean, na.rm = TRUE)
bdod_bin   <- tapply(bdod,      bin, mean, na.rm = TRUE)
STDAGE_bin <- tapply(STDAGE,    bin, mean, na.rm = TRUE)



print(sapply(
  list(FD=FD_bin, FI=FI_bin, CWD=CWD_bin,
       TEM1=TEM1_bin, PRE1=PRE1_bin,
       pH=pH_bin, bdod=bdod_bin, STDAGE=STDAGE_bin),
  function(x) sum(is.na(x))
))


jags_data <- list(
  N            = N,
  N_bin        = N_bin,
  DBH_CV       = y_CV,
  Carbon_Mg_ha = y_AGB,
  trait        = trait,
  trait_CWM    = trait_CWM,
  CWD          = CWD,                        
  TEM1         = TEM1,
  PRE1         = PRE1,
  pH           = pH,
  bdod         = bdod,
  STDAGE       = STDAGE,
  bin          = bin,
  FD_bin       = as.numeric(FD_bin),
  FI_bin       = as.numeric(FI_bin),
  CWD_bin      = as.numeric(CWD_bin),        
  TEM1_bin     = as.numeric(TEM1_bin),
  PRE1_bin     = as.numeric(PRE1_bin),
  pH_bin       = as.numeric(pH_bin),
  bdod_bin     = as.numeric(bdod_bin),
  STDAGE_bin   = as.numeric(STDAGE_bin)
)



##### parameters##########################################


params <- c(
  "delta1", "delta2",
  "beta3",  "beta4",  "beta5",
  "gamma0_CV",
  "mu_alpha_AGB",
  "gamma_CV",
  "gamma_CV_bin",
  "gamma_AGB",
  "beta5_0",
  "beta5_1", "beta5_2",
  "beta5_3", "beta5_4",
  "beta5_5", "beta5_6",
  "beta5_7", "beta5_8",
  "sigma_CV",       "sigma_AGB",
  "sigma_alpha_CV", "sigma_alpha_AGB",
  "sigma_delta1",   "sigma_delta2",
  "sigma_beta3",    "sigma_beta4",    "sigma_beta5"
)



##### Run JAGS #################################################################


fit_sem_jags <- jags(
  model.file         = "/data/bayes_sem_model.txt",
  data               = jags_data,
  parameters.to.save = params,
  n.chains           = 3,
  n.iter             = 10000,
  n.burnin           = 2000,
  n.thin             = 10
)



##### Save #############################################################


post        <- as.data.frame(fit_sem_jags$BUGSoutput$sims.list)
summary_tab <- as.data.frame(fit_sem_jags$BUGSoutput$summary)

# 收敛诊断
rhats     <- fit_sem_jags$BUGSoutput$summary[, "Rhat"]
bad_rhats <- rhats[rhats > 1.1]
if (length(bad_rhats) > 0) {
  cat("（Rhat > 1.1）:\n"); print(bad_rhats)
} else {
  cat("Rhat < 1.1）\n")
}


saveRDS(fit_sem_jags,
        "/data/beiyesi_CWD9_DBH5_LDMC.rds")
write.csv(summary_tab,
          "/data/Beiyrsi_Summary_American_CWD9_DBH5_LDMC.csv",
          row.names = TRUE)
saveRDS(post,
        "/data/post_American_CWD9_DBH5_LDMC.rds")



##### Extract #####################################################



get_effect_by_bin <- function(post, param_prefix, path_label) {
  param_cols  <- grep(paste0("^", param_prefix, "\\."), names(post), value = TRUE)
  summary_mat <- t(apply(
    post[, param_cols, drop = FALSE], 2,
    function(x) c(mean  = mean(x),
                  lower = quantile(x, 0.025),
                  upper = quantile(x, 0.975))
  ))
  df           <- as.data.frame(summary_mat)
  colnames(df) <- c("mean", "lower", "upper")
  df$bin       <- as.numeric(sub(paste0("^", param_prefix, "\\."), "", param_cols))
  df$effect_type <- ifelse(
    df$lower > 0 & df$upper > 0, "positive significant",
    ifelse(df$lower < 0 & df$upper < 0, "negative significant",
           ifelse(df$mean > 0, "positive non-significant", "negative non-significant"))
  )
  df$path <- path_label
  return(df)
}


get_indirect_effect_by_bin <- function(post, delta_prefix, beta_prefix, path_label) {
  delta_cols  <- sort(grep(paste0("^", delta_prefix, "\\."), names(post), value = TRUE))
  beta_cols   <- sort(grep(paste0("^", beta_prefix,  "\\."), names(post), value = TRUE))
  sorted_bins <- as.numeric(sub(paste0("^", delta_prefix, "\\."), "", delta_cols))
  indirect    <- post[, delta_cols] * post[, beta_cols]
  summary_mat <- t(apply(indirect, 2,
                         function(x) c(mean  = mean(x),
                                       lower = quantile(x, 0.025),
                                       upper = quantile(x, 0.975))
  ))
  df           <- as.data.frame(summary_mat)
  colnames(df) <- c("mean", "lower", "upper")
  df$bin       <- sorted_bins
  df$effect_type <- ifelse(
    df$lower > 0 & df$upper > 0, "positive significant",
    ifelse(df$lower < 0 & df$upper < 0, "negative significant",
           ifelse(df$mean > 0, "positive non-significant", "negative non-significant"))
  )
  df$effect_type <- factor(df$effect_type, levels = c(
    "negative significant", "negative non-significant",
    "positive non-significant", "positive significant"))
  df$path <- path_label
  return(df)
}


all_effects <- bind_rows(
  get_effect_by_bin(post, "beta3",  "FD → AGB"),
  get_effect_by_bin(post, "beta4",  "CWM → AGB"),
  get_effect_by_bin(post, "beta5",  "DBHCV → AGB"),
  get_effect_by_bin(post, "delta1", "FD → DBHCV"),
  get_effect_by_bin(post, "delta2", "CWM → DBHCV"),
  get_indirect_effect_by_bin(post, "delta1", "beta5", "FD → DBHCV → AGB"),
  get_indirect_effect_by_bin(post, "delta2", "beta5", "CWM → DBHCV → AGB")
)


data_sem2_org     <- data_sem2
data_sem2_org$bin <- as.numeric(factor(data_sem2_org$PDSI_ATA_bin))

result <- data_sem2_org %>%
  group_by(bin) %>%
  dplyr::summarise(
    mean_MAT    = mean(mean_annual_temp,   na.rm = TRUE),
    mean_PRE    = mean(mean_annual_precip, na.rm = TRUE),
    mean_CWD    = mean(CWD1,           na.rm = TRUE),  
    mean_STDAGE = mean(STDAGE,             na.rm = TRUE),
    .groups = "drop"
  )

all_effects1 <- merge(all_effects, result, by = "bin")

write.csv(all_effects1,
          "/data/posterior_American_CWD9_DBH5_LDMC.csv",
          row.names = FALSE)



##### Proportional stacked bar chart #################################################################


library(scales)

get_effect_prop <- function(post, param_prefix, path_label) {
  param_cols <- grep(paste0("^", param_prefix, "\\."), names(post), value = TRUE)
  df <- t(apply(post[, param_cols, drop = FALSE], 2,
                function(x) c(mean  = mean(x),
                              lower = quantile(x, 0.025),
                              upper = quantile(x, 0.975))
  )) |> as.data.frame()
  colnames(df) <- c("mean", "lower", "upper")
  df$effect_type <- ifelse(
    df$lower > 0 & df$upper > 0, "positive significant",
    ifelse(df$lower < 0 & df$upper < 0, "negative significant",
           ifelse(df$mean > 0, "positive non-significant", "negative non-significant"))
  )
  dplyr::count(df, effect_type) %>%
    mutate(percent = n / sum(n) * 100, path = path_label)
}

get_indirect_prop <- function(post, delta_prefix, beta_prefix, path_label) {
  delta_cols <- sort(grep(paste0("^", delta_prefix, "\\."), names(post), value = TRUE))
  beta_cols  <- sort(grep(paste0("^", beta_prefix,  "\\."), names(post), value = TRUE))
  indirect   <- post[, delta_cols] * post[, beta_cols]
  df <- t(apply(indirect, 2,
                function(x) c(mean  = mean(x),
                              lower = quantile(x, 0.025),
                              upper = quantile(x, 0.975))
  )) |> as.data.frame()
  colnames(df) <- c("mean", "lower", "upper")
  df$effect_type <- ifelse(
    df$lower > 0 & df$upper > 0, "positive significant",
    ifelse(df$lower < 0 & df$upper < 0, "negative significant",
           ifelse(df$mean > 0, "positive non-significant", "negative non-significant"))
  )
  dplyr::count(df, effect_type) %>%
    mutate(percent = n / sum(n) * 100, path = path_label)
}

prop <- bind_rows(
  get_effect_prop(post, "beta3",  "FD → AGB"),
  get_effect_prop(post, "beta4",  "CWM → AGB"),
  get_effect_prop(post, "beta5",  "DBHCV → AGB"),
  get_effect_prop(post, "delta1", "FD → DBHCV"),
  get_effect_prop(post, "delta2", "CWM → DBHCV"),
  get_indirect_prop(post, "delta1", "beta5", "FD → DBHCV → AGB"),
  get_indirect_prop(post, "delta2", "beta5", "CWM → DBHCV → AGB")
)

write.csv(prop,
          "/data/bar_American_CWD9_DBH5_LDMC.csv",
          row.names = FALSE)


prop_plot <- prop %>%
  filter(path %in% c("FD → AGB", "FD → DBHCV → AGB",
                     "CWM → AGB", "CWM → DBHCV → AGB")) %>%
  mutate(
    path = factor(path, levels = c(
      "FD → AGB", "FD → DBHCV → AGB",
      "CWM → AGB", "CWM → DBHCV → AGB"
    )),
    effect_type = factor(effect_type, levels = c(
      "negative significant", "negative non-significant",
      "positive non-significant", "positive significant"
    ))
  )

fig_height <- ggplot(prop_plot,
                     aes(x = percent, y = path, fill = effect_type)) +
  geom_bar(stat = "identity", width = 0.65) +
  scale_x_continuous(labels = percent_format(scale = 1), expand = c(0, 0)) +
  scale_fill_manual(
    values = c(
      "negative significant"     = "#92C5DE",
      "negative non-significant" = "#dbebfa",
      "positive non-significant" = "#f9ebdf",
      "positive significant"     = "#F4A582"
    ),
    drop = FALSE
  ) +
  labs(x = "Proportion of bins (%)", y = NULL, fill = "Effect type") +
  ggtitle("Structural diversity → Carbon") +
  theme_classic(base_size = 13) +
  theme(axis.text.y  = element_text(size = 12),
        legend.position = "right")

fig_height
##===============================================================================
#####Drought.tolerance###################################################################
##===============================================================================
##### Load data ###################################################################

library(data.table)
library(dplyr)
library(R2jags)

C_datatotal_all <- fread("/data/Data.csv" )
colnames( C_datatotal_all  )

C_datatotal_all <- C_datatotal_all %>%
  dplyr::select(-V1) %>%
  filter(Carbon_Mg_ha > 0)
colnames(C_datatotal_all) <- gsub("_mean$", "", colnames(C_datatotal_all))
C_datatotal_all <- C_datatotal_all %>%
  filter(!is.na(CWD_stand_age)) %>%
  filter(!is.na(MAT_stand_age)) %>%
  filter(!is.na(MAP_stand_age)) 

# 合并 bin 分组
C_bin <- read.csv("/data/Class_CWD_9bins")
C_bin <- C_bin[, c("PLT_CN", "PDSI_ATA_bin")]
C_datatotal_all$PLT_CN <- as.character(C_datatotal_all$PLT_CN)
C_bin$PLT_CN            <- as.character(C_bin$PLT_CN)
C_datatotal_all <- merge(C_datatotal_all, C_bin, by = "PLT_CN")
C_datatotal_all$CWD1=C_datatotal_all$CWD_stand_age
# 标准化
C_datatotal_all <- C_datatotal_all %>%
  mutate(
    Carbon_Mg_ha = log(Carbon_Mg_ha+1),
    shannon_dbhR = scale(H_D_5),
           
           
    CWD          = scale(CWD_stand_age),  
    TEM1         = scale(MAT_stand_age),
    PRE1         = scale(MAP_stand_age),
    FDis         = scale(Drought.tolerance.FDis),
    CWM          = scale(Drought.tolerance.CWM),
    pH           = scale(nitrogen),
    bdod         = scale(bdod),
    STDAGE       = scale(STDAGE)
  )

data_sem2 <- C_datatotal_all

##### JAGS #########################################################


data_all <- data_sem2

y_CV      <- as.numeric(data_all$shannon_dbhR)
y_AGB     <- as.numeric(data_all$Carbon_Mg_ha)
trait     <- as.numeric(data_all$FDis)
trait_CWM <- as.numeric(data_all$CWM)
CWD       <- as.numeric(data_all$CWD)       
TEM1      <- as.numeric(data_all$TEM1)
PRE1      <- as.numeric(data_all$PRE1)
pH        <- as.numeric(data_all$pH)
bdod      <- as.numeric(data_all$bdod)
STDAGE    <- as.numeric(data_all$STDAGE)

bin   <- as.numeric(factor(data_all$PDSI_ATA_bin))
N     <- nrow(data_all)
N_bin <- max(bin)

cat("N =", N, "| N_bin =", N_bin, "\n")


FD_bin     <- tapply(trait,     bin, mean, na.rm = TRUE)
FI_bin     <- tapply(trait_CWM, bin, mean, na.rm = TRUE)
CWD_bin    <- tapply(CWD,       bin, mean, na.rm = TRUE)  
TEM1_bin   <- tapply(TEM1,      bin, mean, na.rm = TRUE)
PRE1_bin   <- tapply(PRE1,      bin, mean, na.rm = TRUE)
pH_bin     <- tapply(pH,        bin, mean, na.rm = TRUE)
bdod_bin   <- tapply(bdod,      bin, mean, na.rm = TRUE)
STDAGE_bin <- tapply(STDAGE,    bin, mean, na.rm = TRUE)



print(sapply(
  list(FD=FD_bin, FI=FI_bin, CWD=CWD_bin,
       TEM1=TEM1_bin, PRE1=PRE1_bin,
       pH=pH_bin, bdod=bdod_bin, STDAGE=STDAGE_bin),
  function(x) sum(is.na(x))
))


jags_data <- list(
  N            = N,
  N_bin        = N_bin,
  DBH_CV       = y_CV,
  Carbon_Mg_ha = y_AGB,
  trait        = trait,
  trait_CWM    = trait_CWM,
  CWD          = CWD,                        
  TEM1         = TEM1,
  PRE1         = PRE1,
  pH           = pH,
  bdod         = bdod,
  STDAGE       = STDAGE,
  bin          = bin,
  FD_bin       = as.numeric(FD_bin),
  FI_bin       = as.numeric(FI_bin),
  CWD_bin      = as.numeric(CWD_bin),        
  TEM1_bin     = as.numeric(TEM1_bin),
  PRE1_bin     = as.numeric(PRE1_bin),
  pH_bin       = as.numeric(pH_bin),
  bdod_bin     = as.numeric(bdod_bin),
  STDAGE_bin   = as.numeric(STDAGE_bin)
)



##### parameters##########################################


params <- c(
  "delta1", "delta2",
  "beta3",  "beta4",  "beta5",
  "gamma0_CV",
  "mu_alpha_AGB",
  "gamma_CV",
  "gamma_CV_bin",
  "gamma_AGB",
  "beta5_0",
  "beta5_1", "beta5_2",
  "beta5_3", "beta5_4",
  "beta5_5", "beta5_6",
  "beta5_7", "beta5_8",
  "sigma_CV",       "sigma_AGB",
  "sigma_alpha_CV", "sigma_alpha_AGB",
  "sigma_delta1",   "sigma_delta2",
  "sigma_beta3",    "sigma_beta4",    "sigma_beta5"
)



##### Run JAGS #################################################################


fit_sem_jags <- jags(
  model.file         = "/data/bayes_sem_model.txt",
  data               = jags_data,
  parameters.to.save = params,
  n.chains           = 3,
  n.iter             = 10000,
  n.burnin           = 2000,
  n.thin             = 10
)



##### Save #############################################################


post        <- as.data.frame(fit_sem_jags$BUGSoutput$sims.list)
summary_tab <- as.data.frame(fit_sem_jags$BUGSoutput$summary)

# 收敛诊断
rhats     <- fit_sem_jags$BUGSoutput$summary[, "Rhat"]
bad_rhats <- rhats[rhats > 1.1]
if (length(bad_rhats) > 0) {
  cat("（Rhat > 1.1）:\n"); print(bad_rhats)
} else {
  cat("Rhat < 1.1）\n")
}


saveRDS(fit_sem_jags,
        "/data/beiyesi_CWD9_DBH5_Drought.tolerance.rds")
write.csv(summary_tab,
          "/data/Beiyrsi_Summary_American_CWD9_DBH5_Drought.tolerance.csv",
          row.names = TRUE)
saveRDS(post,
        "/data/post_American_CWD9_DBH5_Drought.tolerance.rds")



##### Extract #####################################################



get_effect_by_bin <- function(post, param_prefix, path_label) {
  param_cols  <- grep(paste0("^", param_prefix, "\\."), names(post), value = TRUE)
  summary_mat <- t(apply(
    post[, param_cols, drop = FALSE], 2,
    function(x) c(mean  = mean(x),
                  lower = quantile(x, 0.025),
                  upper = quantile(x, 0.975))
  ))
  df           <- as.data.frame(summary_mat)
  colnames(df) <- c("mean", "lower", "upper")
  df$bin       <- as.numeric(sub(paste0("^", param_prefix, "\\."), "", param_cols))
  df$effect_type <- ifelse(
    df$lower > 0 & df$upper > 0, "positive significant",
    ifelse(df$lower < 0 & df$upper < 0, "negative significant",
           ifelse(df$mean > 0, "positive non-significant", "negative non-significant"))
  )
  df$path <- path_label
  return(df)
}


get_indirect_effect_by_bin <- function(post, delta_prefix, beta_prefix, path_label) {
  delta_cols  <- sort(grep(paste0("^", delta_prefix, "\\."), names(post), value = TRUE))
  beta_cols   <- sort(grep(paste0("^", beta_prefix,  "\\."), names(post), value = TRUE))
  sorted_bins <- as.numeric(sub(paste0("^", delta_prefix, "\\."), "", delta_cols))
  indirect    <- post[, delta_cols] * post[, beta_cols]
  summary_mat <- t(apply(indirect, 2,
                         function(x) c(mean  = mean(x),
                                       lower = quantile(x, 0.025),
                                       upper = quantile(x, 0.975))
  ))
  df           <- as.data.frame(summary_mat)
  colnames(df) <- c("mean", "lower", "upper")
  df$bin       <- sorted_bins
  df$effect_type <- ifelse(
    df$lower > 0 & df$upper > 0, "positive significant",
    ifelse(df$lower < 0 & df$upper < 0, "negative significant",
           ifelse(df$mean > 0, "positive non-significant", "negative non-significant"))
  )
  df$effect_type <- factor(df$effect_type, levels = c(
    "negative significant", "negative non-significant",
    "positive non-significant", "positive significant"))
  df$path <- path_label
  return(df)
}


all_effects <- bind_rows(
  get_effect_by_bin(post, "beta3",  "FD → AGB"),
  get_effect_by_bin(post, "beta4",  "CWM → AGB"),
  get_effect_by_bin(post, "beta5",  "DBHCV → AGB"),
  get_effect_by_bin(post, "delta1", "FD → DBHCV"),
  get_effect_by_bin(post, "delta2", "CWM → DBHCV"),
  get_indirect_effect_by_bin(post, "delta1", "beta5", "FD → DBHCV → AGB"),
  get_indirect_effect_by_bin(post, "delta2", "beta5", "CWM → DBHCV → AGB")
)


data_sem2_org     <- data_sem2
data_sem2_org$bin <- as.numeric(factor(data_sem2_org$PDSI_ATA_bin))

result <- data_sem2_org %>%
  group_by(bin) %>%
  dplyr::summarise(
    mean_MAT    = mean(mean_annual_temp,   na.rm = TRUE),
    mean_PRE    = mean(mean_annual_precip, na.rm = TRUE),
    mean_CWD    = mean(CWD1,           na.rm = TRUE),  
    mean_STDAGE = mean(STDAGE,             na.rm = TRUE),
    .groups = "drop"
  )

all_effects1 <- merge(all_effects, result, by = "bin")

write.csv(all_effects1,
          "/data/posterior_American_CWD9_DBH5_Drought.tolerance.csv",
          row.names = FALSE)



##### Proportional stacked bar chart #################################################################


library(scales)

get_effect_prop <- function(post, param_prefix, path_label) {
  param_cols <- grep(paste0("^", param_prefix, "\\."), names(post), value = TRUE)
  df <- t(apply(post[, param_cols, drop = FALSE], 2,
                function(x) c(mean  = mean(x),
                              lower = quantile(x, 0.025),
                              upper = quantile(x, 0.975))
  )) |> as.data.frame()
  colnames(df) <- c("mean", "lower", "upper")
  df$effect_type <- ifelse(
    df$lower > 0 & df$upper > 0, "positive significant",
    ifelse(df$lower < 0 & df$upper < 0, "negative significant",
           ifelse(df$mean > 0, "positive non-significant", "negative non-significant"))
  )
  dplyr::count(df, effect_type) %>%
    mutate(percent = n / sum(n) * 100, path = path_label)
}

get_indirect_prop <- function(post, delta_prefix, beta_prefix, path_label) {
  delta_cols <- sort(grep(paste0("^", delta_prefix, "\\."), names(post), value = TRUE))
  beta_cols  <- sort(grep(paste0("^", beta_prefix,  "\\."), names(post), value = TRUE))
  indirect   <- post[, delta_cols] * post[, beta_cols]
  df <- t(apply(indirect, 2,
                function(x) c(mean  = mean(x),
                              lower = quantile(x, 0.025),
                              upper = quantile(x, 0.975))
  )) |> as.data.frame()
  colnames(df) <- c("mean", "lower", "upper")
  df$effect_type <- ifelse(
    df$lower > 0 & df$upper > 0, "positive significant",
    ifelse(df$lower < 0 & df$upper < 0, "negative significant",
           ifelse(df$mean > 0, "positive non-significant", "negative non-significant"))
  )
  dplyr::count(df, effect_type) %>%
    mutate(percent = n / sum(n) * 100, path = path_label)
}

prop <- bind_rows(
  get_effect_prop(post, "beta3",  "FD → AGB"),
  get_effect_prop(post, "beta4",  "CWM → AGB"),
  get_effect_prop(post, "beta5",  "DBHCV → AGB"),
  get_effect_prop(post, "delta1", "FD → DBHCV"),
  get_effect_prop(post, "delta2", "CWM → DBHCV"),
  get_indirect_prop(post, "delta1", "beta5", "FD → DBHCV → AGB"),
  get_indirect_prop(post, "delta2", "beta5", "CWM → DBHCV → AGB")
)

write.csv(prop,
          "/data/bar_American_CWD9_DBH5_Drought.tolerance.csv",
          row.names = FALSE)


prop_plot <- prop %>%
  filter(path %in% c("FD → AGB", "FD → DBHCV → AGB",
                     "CWM → AGB", "CWM → DBHCV → AGB")) %>%
  mutate(
    path = factor(path, levels = c(
      "FD → AGB", "FD → DBHCV → AGB",
      "CWM → AGB", "CWM → DBHCV → AGB"
    )),
    effect_type = factor(effect_type, levels = c(
      "negative significant", "negative non-significant",
      "positive non-significant", "positive significant"
    ))
  )

fig_height <- ggplot(prop_plot,
                     aes(x = percent, y = path, fill = effect_type)) +
  geom_bar(stat = "identity", width = 0.65) +
  scale_x_continuous(labels = percent_format(scale = 1), expand = c(0, 0)) +
  scale_fill_manual(
    values = c(
      "negative significant"     = "#92C5DE",
      "negative non-significant" = "#dbebfa",
      "positive non-significant" = "#f9ebdf",
      "positive significant"     = "#F4A582"
    ),
    drop = FALSE
  ) +
  labs(x = "Proportion of bins (%)", y = NULL, fill = "Effect type") +
  ggtitle("Structural diversity → Carbon") +
  theme_classic(base_size = 13) +
  theme(axis.text.y  = element_text(size = 12),
        legend.position = "right")

fig_height
##===============================================================================
#####Nmass###################################################################
##===============================================================================
##### Load data ###################################################################

library(data.table)
library(dplyr)
library(R2jags)

C_datatotal_all <- fread("/data/Data.csv" )
colnames( C_datatotal_all  )

C_datatotal_all <- C_datatotal_all %>%
  dplyr::select(-V1) %>%
  filter(Carbon_Mg_ha > 0)
colnames(C_datatotal_all) <- gsub("_mean$", "", colnames(C_datatotal_all))
C_datatotal_all <- C_datatotal_all %>%
  filter(!is.na(CWD_stand_age)) %>%
  filter(!is.na(MAT_stand_age)) %>%
  filter(!is.na(MAP_stand_age)) 

# 合并 bin 分组
C_bin <- read.csv("/data/Class_CWD_9bins")
C_bin <- C_bin[, c("PLT_CN", "PDSI_ATA_bin")]
C_datatotal_all$PLT_CN <- as.character(C_datatotal_all$PLT_CN)
C_bin$PLT_CN            <- as.character(C_bin$PLT_CN)
C_datatotal_all <- merge(C_datatotal_all, C_bin, by = "PLT_CN")
C_datatotal_all$CWD1=C_datatotal_all$CWD_stand_age
# 标准化
C_datatotal_all <- C_datatotal_all %>%
  mutate(
    Carbon_Mg_ha = log(Carbon_Mg_ha+1),
    shannon_dbhR = scale(H_D_5),
           
           
    CWD          = scale(CWD_stand_age),  
    TEM1         = scale(MAT_stand_age),
    PRE1         = scale(MAP_stand_age),
    FDis         = scale(Nmass.FDis),
    CWM          = scale(Nmass.CWM),
    pH           = scale(nitrogen),
    bdod         = scale(bdod),
    STDAGE       = scale(STDAGE)
  )

data_sem2 <- C_datatotal_all

##### JAGS #########################################################


data_all <- data_sem2

y_CV      <- as.numeric(data_all$shannon_dbhR)
y_AGB     <- as.numeric(data_all$Carbon_Mg_ha)
trait     <- as.numeric(data_all$FDis)
trait_CWM <- as.numeric(data_all$CWM)
CWD       <- as.numeric(data_all$CWD)       
TEM1      <- as.numeric(data_all$TEM1)
PRE1      <- as.numeric(data_all$PRE1)
pH        <- as.numeric(data_all$pH)
bdod      <- as.numeric(data_all$bdod)
STDAGE    <- as.numeric(data_all$STDAGE)

bin   <- as.numeric(factor(data_all$PDSI_ATA_bin))
N     <- nrow(data_all)
N_bin <- max(bin)

cat("N =", N, "| N_bin =", N_bin, "\n")


FD_bin     <- tapply(trait,     bin, mean, na.rm = TRUE)
FI_bin     <- tapply(trait_CWM, bin, mean, na.rm = TRUE)
CWD_bin    <- tapply(CWD,       bin, mean, na.rm = TRUE)  
TEM1_bin   <- tapply(TEM1,      bin, mean, na.rm = TRUE)
PRE1_bin   <- tapply(PRE1,      bin, mean, na.rm = TRUE)
pH_bin     <- tapply(pH,        bin, mean, na.rm = TRUE)
bdod_bin   <- tapply(bdod,      bin, mean, na.rm = TRUE)
STDAGE_bin <- tapply(STDAGE,    bin, mean, na.rm = TRUE)



print(sapply(
  list(FD=FD_bin, FI=FI_bin, CWD=CWD_bin,
       TEM1=TEM1_bin, PRE1=PRE1_bin,
       pH=pH_bin, bdod=bdod_bin, STDAGE=STDAGE_bin),
  function(x) sum(is.na(x))
))


jags_data <- list(
  N            = N,
  N_bin        = N_bin,
  DBH_CV       = y_CV,
  Carbon_Mg_ha = y_AGB,
  trait        = trait,
  trait_CWM    = trait_CWM,
  CWD          = CWD,                        
  TEM1         = TEM1,
  PRE1         = PRE1,
  pH           = pH,
  bdod         = bdod,
  STDAGE       = STDAGE,
  bin          = bin,
  FD_bin       = as.numeric(FD_bin),
  FI_bin       = as.numeric(FI_bin),
  CWD_bin      = as.numeric(CWD_bin),        
  TEM1_bin     = as.numeric(TEM1_bin),
  PRE1_bin     = as.numeric(PRE1_bin),
  pH_bin       = as.numeric(pH_bin),
  bdod_bin     = as.numeric(bdod_bin),
  STDAGE_bin   = as.numeric(STDAGE_bin)
)



##### parameters##########################################


params <- c(
  "delta1", "delta2",
  "beta3",  "beta4",  "beta5",
  "gamma0_CV",
  "mu_alpha_AGB",
  "gamma_CV",
  "gamma_CV_bin",
  "gamma_AGB",
  "beta5_0",
  "beta5_1", "beta5_2",
  "beta5_3", "beta5_4",
  "beta5_5", "beta5_6",
  "beta5_7", "beta5_8",
  "sigma_CV",       "sigma_AGB",
  "sigma_alpha_CV", "sigma_alpha_AGB",
  "sigma_delta1",   "sigma_delta2",
  "sigma_beta3",    "sigma_beta4",    "sigma_beta5"
)



##### Run JAGS #################################################################


fit_sem_jags <- jags(
  model.file         = "/data/bayes_sem_model.txt",
  data               = jags_data,
  parameters.to.save = params,
  n.chains           = 3,
  n.iter             = 10000,
  n.burnin           = 2000,
  n.thin             = 10
)



##### Save #############################################################


post        <- as.data.frame(fit_sem_jags$BUGSoutput$sims.list)
summary_tab <- as.data.frame(fit_sem_jags$BUGSoutput$summary)

# 收敛诊断
rhats     <- fit_sem_jags$BUGSoutput$summary[, "Rhat"]
bad_rhats <- rhats[rhats > 1.1]
if (length(bad_rhats) > 0) {
  cat("（Rhat > 1.1）:\n"); print(bad_rhats)
} else {
  cat("Rhat < 1.1）\n")
}


saveRDS(fit_sem_jags,
        "/data/beiyesi_CWD9_DBH5_Nmass.rds")
write.csv(summary_tab,
          "/data/Beiyrsi_Summary_American_CWD9_DBH5_Nmass.csv",
          row.names = TRUE)
saveRDS(post,
        "/data/post_American_CWD9_DBH5_Nmass.rds")



##### Extract #####################################################



get_effect_by_bin <- function(post, param_prefix, path_label) {
  param_cols  <- grep(paste0("^", param_prefix, "\\."), names(post), value = TRUE)
  summary_mat <- t(apply(
    post[, param_cols, drop = FALSE], 2,
    function(x) c(mean  = mean(x),
                  lower = quantile(x, 0.025),
                  upper = quantile(x, 0.975))
  ))
  df           <- as.data.frame(summary_mat)
  colnames(df) <- c("mean", "lower", "upper")
  df$bin       <- as.numeric(sub(paste0("^", param_prefix, "\\."), "", param_cols))
  df$effect_type <- ifelse(
    df$lower > 0 & df$upper > 0, "positive significant",
    ifelse(df$lower < 0 & df$upper < 0, "negative significant",
           ifelse(df$mean > 0, "positive non-significant", "negative non-significant"))
  )
  df$path <- path_label
  return(df)
}


get_indirect_effect_by_bin <- function(post, delta_prefix, beta_prefix, path_label) {
  delta_cols  <- sort(grep(paste0("^", delta_prefix, "\\."), names(post), value = TRUE))
  beta_cols   <- sort(grep(paste0("^", beta_prefix,  "\\."), names(post), value = TRUE))
  sorted_bins <- as.numeric(sub(paste0("^", delta_prefix, "\\."), "", delta_cols))
  indirect    <- post[, delta_cols] * post[, beta_cols]
  summary_mat <- t(apply(indirect, 2,
                         function(x) c(mean  = mean(x),
                                       lower = quantile(x, 0.025),
                                       upper = quantile(x, 0.975))
  ))
  df           <- as.data.frame(summary_mat)
  colnames(df) <- c("mean", "lower", "upper")
  df$bin       <- sorted_bins
  df$effect_type <- ifelse(
    df$lower > 0 & df$upper > 0, "positive significant",
    ifelse(df$lower < 0 & df$upper < 0, "negative significant",
           ifelse(df$mean > 0, "positive non-significant", "negative non-significant"))
  )
  df$effect_type <- factor(df$effect_type, levels = c(
    "negative significant", "negative non-significant",
    "positive non-significant", "positive significant"))
  df$path <- path_label
  return(df)
}


all_effects <- bind_rows(
  get_effect_by_bin(post, "beta3",  "FD → AGB"),
  get_effect_by_bin(post, "beta4",  "CWM → AGB"),
  get_effect_by_bin(post, "beta5",  "DBHCV → AGB"),
  get_effect_by_bin(post, "delta1", "FD → DBHCV"),
  get_effect_by_bin(post, "delta2", "CWM → DBHCV"),
  get_indirect_effect_by_bin(post, "delta1", "beta5", "FD → DBHCV → AGB"),
  get_indirect_effect_by_bin(post, "delta2", "beta5", "CWM → DBHCV → AGB")
)


data_sem2_org     <- data_sem2
data_sem2_org$bin <- as.numeric(factor(data_sem2_org$PDSI_ATA_bin))

result <- data_sem2_org %>%
  group_by(bin) %>%
  dplyr::summarise(
    mean_MAT    = mean(mean_annual_temp,   na.rm = TRUE),
    mean_PRE    = mean(mean_annual_precip, na.rm = TRUE),
    mean_CWD    = mean(CWD1,           na.rm = TRUE),  
    mean_STDAGE = mean(STDAGE,             na.rm = TRUE),
    .groups = "drop"
  )

all_effects1 <- merge(all_effects, result, by = "bin")

write.csv(all_effects1,
          "/data/posterior_American_CWD9_DBH5_Nmass.csv",
          row.names = FALSE)



##### Proportional stacked bar chart #################################################################


library(scales)

get_effect_prop <- function(post, param_prefix, path_label) {
  param_cols <- grep(paste0("^", param_prefix, "\\."), names(post), value = TRUE)
  df <- t(apply(post[, param_cols, drop = FALSE], 2,
                function(x) c(mean  = mean(x),
                              lower = quantile(x, 0.025),
                              upper = quantile(x, 0.975))
  )) |> as.data.frame()
  colnames(df) <- c("mean", "lower", "upper")
  df$effect_type <- ifelse(
    df$lower > 0 & df$upper > 0, "positive significant",
    ifelse(df$lower < 0 & df$upper < 0, "negative significant",
           ifelse(df$mean > 0, "positive non-significant", "negative non-significant"))
  )
  dplyr::count(df, effect_type) %>%
    mutate(percent = n / sum(n) * 100, path = path_label)
}

get_indirect_prop <- function(post, delta_prefix, beta_prefix, path_label) {
  delta_cols <- sort(grep(paste0("^", delta_prefix, "\\."), names(post), value = TRUE))
  beta_cols  <- sort(grep(paste0("^", beta_prefix,  "\\."), names(post), value = TRUE))
  indirect   <- post[, delta_cols] * post[, beta_cols]
  df <- t(apply(indirect, 2,
                function(x) c(mean  = mean(x),
                              lower = quantile(x, 0.025),
                              upper = quantile(x, 0.975))
  )) |> as.data.frame()
  colnames(df) <- c("mean", "lower", "upper")
  df$effect_type <- ifelse(
    df$lower > 0 & df$upper > 0, "positive significant",
    ifelse(df$lower < 0 & df$upper < 0, "negative significant",
           ifelse(df$mean > 0, "positive non-significant", "negative non-significant"))
  )
  dplyr::count(df, effect_type) %>%
    mutate(percent = n / sum(n) * 100, path = path_label)
}

prop <- bind_rows(
  get_effect_prop(post, "beta3",  "FD → AGB"),
  get_effect_prop(post, "beta4",  "CWM → AGB"),
  get_effect_prop(post, "beta5",  "DBHCV → AGB"),
  get_effect_prop(post, "delta1", "FD → DBHCV"),
  get_effect_prop(post, "delta2", "CWM → DBHCV"),
  get_indirect_prop(post, "delta1", "beta5", "FD → DBHCV → AGB"),
  get_indirect_prop(post, "delta2", "beta5", "CWM → DBHCV → AGB")
)

write.csv(prop,
          "/data/bar_American_CWD9_DBH5_Nmass.csv",
          row.names = FALSE)


prop_plot <- prop %>%
  filter(path %in% c("FD → AGB", "FD → DBHCV → AGB",
                     "CWM → AGB", "CWM → DBHCV → AGB")) %>%
  mutate(
    path = factor(path, levels = c(
      "FD → AGB", "FD → DBHCV → AGB",
      "CWM → AGB", "CWM → DBHCV → AGB"
    )),
    effect_type = factor(effect_type, levels = c(
      "negative significant", "negative non-significant",
      "positive non-significant", "positive significant"
    ))
  )

fig_height <- ggplot(prop_plot,
                     aes(x = percent, y = path, fill = effect_type)) +
  geom_bar(stat = "identity", width = 0.65) +
  scale_x_continuous(labels = percent_format(scale = 1), expand = c(0, 0)) +
  scale_fill_manual(
    values = c(
      "negative significant"     = "#92C5DE",
      "negative non-significant" = "#dbebfa",
      "positive non-significant" = "#f9ebdf",
      "positive significant"     = "#F4A582"
    ),
    drop = FALSE
  ) +
  labs(x = "Proportion of bins (%)", y = NULL, fill = "Effect type") +
  ggtitle("Structural diversity → Carbon") +
  theme_classic(base_size = 13) +
  theme(axis.text.y  = element_text(size = 12),
        legend.position = "right")

fig_height
##===============================================================================
#####Pmass###################################################################
##===============================================================================
##### Load data ###################################################################

library(data.table)
library(dplyr)
library(R2jags)

C_datatotal_all <- fread("/data/Data.csv" )
colnames( C_datatotal_all  )

C_datatotal_all <- C_datatotal_all %>%
  dplyr::select(-V1) %>%
  filter(Carbon_Mg_ha > 0)
colnames(C_datatotal_all) <- gsub("_mean$", "", colnames(C_datatotal_all))
C_datatotal_all <- C_datatotal_all %>%
  filter(!is.na(CWD_stand_age)) %>%
  filter(!is.na(MAT_stand_age)) %>%
  filter(!is.na(MAP_stand_age)) 

# 合并 bin 分组
C_bin <- read.csv("/data/Class_CWD_9bins")
C_bin <- C_bin[, c("PLT_CN", "PDSI_ATA_bin")]
C_datatotal_all$PLT_CN <- as.character(C_datatotal_all$PLT_CN)
C_bin$PLT_CN            <- as.character(C_bin$PLT_CN)
C_datatotal_all <- merge(C_datatotal_all, C_bin, by = "PLT_CN")
C_datatotal_all$CWD1=C_datatotal_all$CWD_stand_age
# 标准化
C_datatotal_all <- C_datatotal_all %>%
  mutate(
    Carbon_Mg_ha = log(Carbon_Mg_ha+1),
    shannon_dbhR = scale(H_D_5),
           
           
    CWD          = scale(CWD_stand_age),  
    TEM1         = scale(MAT_stand_age),
    PRE1         = scale(MAP_stand_age),
    FDis         = scale(Pmass.FDis),
    CWM          = scale(Pmass.CWM),
    pH           = scale(nitrogen),
    bdod         = scale(bdod),
    STDAGE       = scale(STDAGE)
  )

data_sem2 <- C_datatotal_all

##### JAGS #########################################################


data_all <- data_sem2

y_CV      <- as.numeric(data_all$shannon_dbhR)
y_AGB     <- as.numeric(data_all$Carbon_Mg_ha)
trait     <- as.numeric(data_all$FDis)
trait_CWM <- as.numeric(data_all$CWM)
CWD       <- as.numeric(data_all$CWD)       
TEM1      <- as.numeric(data_all$TEM1)
PRE1      <- as.numeric(data_all$PRE1)
pH        <- as.numeric(data_all$pH)
bdod      <- as.numeric(data_all$bdod)
STDAGE    <- as.numeric(data_all$STDAGE)

bin   <- as.numeric(factor(data_all$PDSI_ATA_bin))
N     <- nrow(data_all)
N_bin <- max(bin)

cat("N =", N, "| N_bin =", N_bin, "\n")


FD_bin     <- tapply(trait,     bin, mean, na.rm = TRUE)
FI_bin     <- tapply(trait_CWM, bin, mean, na.rm = TRUE)
CWD_bin    <- tapply(CWD,       bin, mean, na.rm = TRUE)  
TEM1_bin   <- tapply(TEM1,      bin, mean, na.rm = TRUE)
PRE1_bin   <- tapply(PRE1,      bin, mean, na.rm = TRUE)
pH_bin     <- tapply(pH,        bin, mean, na.rm = TRUE)
bdod_bin   <- tapply(bdod,      bin, mean, na.rm = TRUE)
STDAGE_bin <- tapply(STDAGE,    bin, mean, na.rm = TRUE)



print(sapply(
  list(FD=FD_bin, FI=FI_bin, CWD=CWD_bin,
       TEM1=TEM1_bin, PRE1=PRE1_bin,
       pH=pH_bin, bdod=bdod_bin, STDAGE=STDAGE_bin),
  function(x) sum(is.na(x))
))


jags_data <- list(
  N            = N,
  N_bin        = N_bin,
  DBH_CV       = y_CV,
  Carbon_Mg_ha = y_AGB,
  trait        = trait,
  trait_CWM    = trait_CWM,
  CWD          = CWD,                        
  TEM1         = TEM1,
  PRE1         = PRE1,
  pH           = pH,
  bdod         = bdod,
  STDAGE       = STDAGE,
  bin          = bin,
  FD_bin       = as.numeric(FD_bin),
  FI_bin       = as.numeric(FI_bin),
  CWD_bin      = as.numeric(CWD_bin),        
  TEM1_bin     = as.numeric(TEM1_bin),
  PRE1_bin     = as.numeric(PRE1_bin),
  pH_bin       = as.numeric(pH_bin),
  bdod_bin     = as.numeric(bdod_bin),
  STDAGE_bin   = as.numeric(STDAGE_bin)
)



##### parameters##########################################


params <- c(
  "delta1", "delta2",
  "beta3",  "beta4",  "beta5",
  "gamma0_CV",
  "mu_alpha_AGB",
  "gamma_CV",
  "gamma_CV_bin",
  "gamma_AGB",
  "beta5_0",
  "beta5_1", "beta5_2",
  "beta5_3", "beta5_4",
  "beta5_5", "beta5_6",
  "beta5_7", "beta5_8",
  "sigma_CV",       "sigma_AGB",
  "sigma_alpha_CV", "sigma_alpha_AGB",
  "sigma_delta1",   "sigma_delta2",
  "sigma_beta3",    "sigma_beta4",    "sigma_beta5"
)



##### Run JAGS #################################################################


fit_sem_jags <- jags(
  model.file         = "/data/bayes_sem_model.txt",
  data               = jags_data,
  parameters.to.save = params,
  n.chains           = 3,
  n.iter             = 10000,
  n.burnin           = 2000,
  n.thin             = 10
)



##### Save #############################################################


post        <- as.data.frame(fit_sem_jags$BUGSoutput$sims.list)
summary_tab <- as.data.frame(fit_sem_jags$BUGSoutput$summary)

# 收敛诊断
rhats     <- fit_sem_jags$BUGSoutput$summary[, "Rhat"]
bad_rhats <- rhats[rhats > 1.1]
if (length(bad_rhats) > 0) {
  cat("（Rhat > 1.1）:\n"); print(bad_rhats)
} else {
  cat("Rhat < 1.1）\n")
}


saveRDS(fit_sem_jags,
        "/data/beiyesi_CWD9_DBH5_Pmass.rds")
write.csv(summary_tab,
          "/data/Beiyrsi_Summary_American_CWD9_DBH5_Pmass.csv",
          row.names = TRUE)
saveRDS(post,
        "/data/post_American_CWD9_DBH5_Pmass.rds")



##### Extract #####################################################



get_effect_by_bin <- function(post, param_prefix, path_label) {
  param_cols  <- grep(paste0("^", param_prefix, "\\."), names(post), value = TRUE)
  summary_mat <- t(apply(
    post[, param_cols, drop = FALSE], 2,
    function(x) c(mean  = mean(x),
                  lower = quantile(x, 0.025),
                  upper = quantile(x, 0.975))
  ))
  df           <- as.data.frame(summary_mat)
  colnames(df) <- c("mean", "lower", "upper")
  df$bin       <- as.numeric(sub(paste0("^", param_prefix, "\\."), "", param_cols))
  df$effect_type <- ifelse(
    df$lower > 0 & df$upper > 0, "positive significant",
    ifelse(df$lower < 0 & df$upper < 0, "negative significant",
           ifelse(df$mean > 0, "positive non-significant", "negative non-significant"))
  )
  df$path <- path_label
  return(df)
}


get_indirect_effect_by_bin <- function(post, delta_prefix, beta_prefix, path_label) {
  delta_cols  <- sort(grep(paste0("^", delta_prefix, "\\."), names(post), value = TRUE))
  beta_cols   <- sort(grep(paste0("^", beta_prefix,  "\\."), names(post), value = TRUE))
  sorted_bins <- as.numeric(sub(paste0("^", delta_prefix, "\\."), "", delta_cols))
  indirect    <- post[, delta_cols] * post[, beta_cols]
  summary_mat <- t(apply(indirect, 2,
                         function(x) c(mean  = mean(x),
                                       lower = quantile(x, 0.025),
                                       upper = quantile(x, 0.975))
  ))
  df           <- as.data.frame(summary_mat)
  colnames(df) <- c("mean", "lower", "upper")
  df$bin       <- sorted_bins
  df$effect_type <- ifelse(
    df$lower > 0 & df$upper > 0, "positive significant",
    ifelse(df$lower < 0 & df$upper < 0, "negative significant",
           ifelse(df$mean > 0, "positive non-significant", "negative non-significant"))
  )
  df$effect_type <- factor(df$effect_type, levels = c(
    "negative significant", "negative non-significant",
    "positive non-significant", "positive significant"))
  df$path <- path_label
  return(df)
}


all_effects <- bind_rows(
  get_effect_by_bin(post, "beta3",  "FD → AGB"),
  get_effect_by_bin(post, "beta4",  "CWM → AGB"),
  get_effect_by_bin(post, "beta5",  "DBHCV → AGB"),
  get_effect_by_bin(post, "delta1", "FD → DBHCV"),
  get_effect_by_bin(post, "delta2", "CWM → DBHCV"),
  get_indirect_effect_by_bin(post, "delta1", "beta5", "FD → DBHCV → AGB"),
  get_indirect_effect_by_bin(post, "delta2", "beta5", "CWM → DBHCV → AGB")
)


data_sem2_org     <- data_sem2
data_sem2_org$bin <- as.numeric(factor(data_sem2_org$PDSI_ATA_bin))

result <- data_sem2_org %>%
  group_by(bin) %>%
  dplyr::summarise(
    mean_MAT    = mean(mean_annual_temp,   na.rm = TRUE),
    mean_PRE    = mean(mean_annual_precip, na.rm = TRUE),
    mean_CWD    = mean(CWD1,           na.rm = TRUE),  
    mean_STDAGE = mean(STDAGE,             na.rm = TRUE),
    .groups = "drop"
  )

all_effects1 <- merge(all_effects, result, by = "bin")

write.csv(all_effects1,
          "/data/posterior_American_CWD9_DBH5_Pmass.csv",
          row.names = FALSE)



##### Proportional stacked bar chart #################################################################


library(scales)

get_effect_prop <- function(post, param_prefix, path_label) {
  param_cols <- grep(paste0("^", param_prefix, "\\."), names(post), value = TRUE)
  df <- t(apply(post[, param_cols, drop = FALSE], 2,
                function(x) c(mean  = mean(x),
                              lower = quantile(x, 0.025),
                              upper = quantile(x, 0.975))
  )) |> as.data.frame()
  colnames(df) <- c("mean", "lower", "upper")
  df$effect_type <- ifelse(
    df$lower > 0 & df$upper > 0, "positive significant",
    ifelse(df$lower < 0 & df$upper < 0, "negative significant",
           ifelse(df$mean > 0, "positive non-significant", "negative non-significant"))
  )
  dplyr::count(df, effect_type) %>%
    mutate(percent = n / sum(n) * 100, path = path_label)
}

get_indirect_prop <- function(post, delta_prefix, beta_prefix, path_label) {
  delta_cols <- sort(grep(paste0("^", delta_prefix, "\\."), names(post), value = TRUE))
  beta_cols  <- sort(grep(paste0("^", beta_prefix,  "\\."), names(post), value = TRUE))
  indirect   <- post[, delta_cols] * post[, beta_cols]
  df <- t(apply(indirect, 2,
                function(x) c(mean  = mean(x),
                              lower = quantile(x, 0.025),
                              upper = quantile(x, 0.975))
  )) |> as.data.frame()
  colnames(df) <- c("mean", "lower", "upper")
  df$effect_type <- ifelse(
    df$lower > 0 & df$upper > 0, "positive significant",
    ifelse(df$lower < 0 & df$upper < 0, "negative significant",
           ifelse(df$mean > 0, "positive non-significant", "negative non-significant"))
  )
  dplyr::count(df, effect_type) %>%
    mutate(percent = n / sum(n) * 100, path = path_label)
}

prop <- bind_rows(
  get_effect_prop(post, "beta3",  "FD → AGB"),
  get_effect_prop(post, "beta4",  "CWM → AGB"),
  get_effect_prop(post, "beta5",  "DBHCV → AGB"),
  get_effect_prop(post, "delta1", "FD → DBHCV"),
  get_effect_prop(post, "delta2", "CWM → DBHCV"),
  get_indirect_prop(post, "delta1", "beta5", "FD → DBHCV → AGB"),
  get_indirect_prop(post, "delta2", "beta5", "CWM → DBHCV → AGB")
)

write.csv(prop,
          "/data/bar_American_CWD9_DBH5_Pmass.csv",
          row.names = FALSE)


prop_plot <- prop %>%
  filter(path %in% c("FD → AGB", "FD → DBHCV → AGB",
                     "CWM → AGB", "CWM → DBHCV → AGB")) %>%
  mutate(
    path = factor(path, levels = c(
      "FD → AGB", "FD → DBHCV → AGB",
      "CWM → AGB", "CWM → DBHCV → AGB"
    )),
    effect_type = factor(effect_type, levels = c(
      "negative significant", "negative non-significant",
      "positive non-significant", "positive significant"
    ))
  )

fig_height <- ggplot(prop_plot,
                     aes(x = percent, y = path, fill = effect_type)) +
  geom_bar(stat = "identity", width = 0.65) +
  scale_x_continuous(labels = percent_format(scale = 1), expand = c(0, 0)) +
  scale_fill_manual(
    values = c(
      "negative significant"     = "#92C5DE",
      "negative non-significant" = "#dbebfa",
      "positive non-significant" = "#f9ebdf",
      "positive significant"     = "#F4A582"
    ),
    drop = FALSE
  ) +
  labs(x = "Proportion of bins (%)", y = NULL, fill = "Effect type") +
  ggtitle("Structural diversity → Carbon") +
  theme_classic(base_size = 13) +
  theme(axis.text.y  = element_text(size = 12),
        legend.position = "right")

fig_height
##===============================================================================
#####SLA###################################################################
##===============================================================================
##### Load data ###################################################################

library(data.table)
library(dplyr)
library(R2jags)

C_datatotal_all <- fread("/data/Data.csv" )
colnames( C_datatotal_all  )

C_datatotal_all <- C_datatotal_all %>%
  dplyr::select(-V1) %>%
  filter(Carbon_Mg_ha > 0)
colnames(C_datatotal_all) <- gsub("_mean$", "", colnames(C_datatotal_all))
C_datatotal_all <- C_datatotal_all %>%
  filter(!is.na(CWD_stand_age)) %>%
  filter(!is.na(MAT_stand_age)) %>%
  filter(!is.na(MAP_stand_age)) 

# 合并 bin 分组
C_bin <- read.csv("/data/Class_CWD_9bins")
C_bin <- C_bin[, c("PLT_CN", "PDSI_ATA_bin")]
C_datatotal_all$PLT_CN <- as.character(C_datatotal_all$PLT_CN)
C_bin$PLT_CN            <- as.character(C_bin$PLT_CN)
C_datatotal_all <- merge(C_datatotal_all, C_bin, by = "PLT_CN")
C_datatotal_all$CWD1=C_datatotal_all$CWD_stand_age
# 标准化
C_datatotal_all <- C_datatotal_all %>%
  mutate(
    Carbon_Mg_ha = log(Carbon_Mg_ha+1),
    shannon_dbhR = scale(H_D_5),
           
           
    CWD          = scale(CWD_stand_age),  
    TEM1         = scale(MAT_stand_age),
    PRE1         = scale(MAP_stand_age),
    FDis         = scale(SLA.FDis),
    CWM          = scale(SLA.CWM),
    pH           = scale(nitrogen),
    bdod         = scale(bdod),
    STDAGE       = scale(STDAGE)
  )

data_sem2 <- C_datatotal_all

##### JAGS #########################################################


data_all <- data_sem2

y_CV      <- as.numeric(data_all$shannon_dbhR)
y_AGB     <- as.numeric(data_all$Carbon_Mg_ha)
trait     <- as.numeric(data_all$FDis)
trait_CWM <- as.numeric(data_all$CWM)
CWD       <- as.numeric(data_all$CWD)       
TEM1      <- as.numeric(data_all$TEM1)
PRE1      <- as.numeric(data_all$PRE1)
pH        <- as.numeric(data_all$pH)
bdod      <- as.numeric(data_all$bdod)
STDAGE    <- as.numeric(data_all$STDAGE)

bin   <- as.numeric(factor(data_all$PDSI_ATA_bin))
N     <- nrow(data_all)
N_bin <- max(bin)

cat("N =", N, "| N_bin =", N_bin, "\n")


FD_bin     <- tapply(trait,     bin, mean, na.rm = TRUE)
FI_bin     <- tapply(trait_CWM, bin, mean, na.rm = TRUE)
CWD_bin    <- tapply(CWD,       bin, mean, na.rm = TRUE)  
TEM1_bin   <- tapply(TEM1,      bin, mean, na.rm = TRUE)
PRE1_bin   <- tapply(PRE1,      bin, mean, na.rm = TRUE)
pH_bin     <- tapply(pH,        bin, mean, na.rm = TRUE)
bdod_bin   <- tapply(bdod,      bin, mean, na.rm = TRUE)
STDAGE_bin <- tapply(STDAGE,    bin, mean, na.rm = TRUE)



print(sapply(
  list(FD=FD_bin, FI=FI_bin, CWD=CWD_bin,
       TEM1=TEM1_bin, PRE1=PRE1_bin,
       pH=pH_bin, bdod=bdod_bin, STDAGE=STDAGE_bin),
  function(x) sum(is.na(x))
))


jags_data <- list(
  N            = N,
  N_bin        = N_bin,
  DBH_CV       = y_CV,
  Carbon_Mg_ha = y_AGB,
  trait        = trait,
  trait_CWM    = trait_CWM,
  CWD          = CWD,                        
  TEM1         = TEM1,
  PRE1         = PRE1,
  pH           = pH,
  bdod         = bdod,
  STDAGE       = STDAGE,
  bin          = bin,
  FD_bin       = as.numeric(FD_bin),
  FI_bin       = as.numeric(FI_bin),
  CWD_bin      = as.numeric(CWD_bin),        
  TEM1_bin     = as.numeric(TEM1_bin),
  PRE1_bin     = as.numeric(PRE1_bin),
  pH_bin       = as.numeric(pH_bin),
  bdod_bin     = as.numeric(bdod_bin),
  STDAGE_bin   = as.numeric(STDAGE_bin)
)



##### parameters##########################################


params <- c(
  "delta1", "delta2",
  "beta3",  "beta4",  "beta5",
  "gamma0_CV",
  "mu_alpha_AGB",
  "gamma_CV",
  "gamma_CV_bin",
  "gamma_AGB",
  "beta5_0",
  "beta5_1", "beta5_2",
  "beta5_3", "beta5_4",
  "beta5_5", "beta5_6",
  "beta5_7", "beta5_8",
  "sigma_CV",       "sigma_AGB",
  "sigma_alpha_CV", "sigma_alpha_AGB",
  "sigma_delta1",   "sigma_delta2",
  "sigma_beta3",    "sigma_beta4",    "sigma_beta5"
)



##### Run JAGS #################################################################


fit_sem_jags <- jags(
  model.file         = "/data/bayes_sem_model.txt",
  data               = jags_data,
  parameters.to.save = params,
  n.chains           = 3,
  n.iter             = 10000,
  n.burnin           = 2000,
  n.thin             = 10
)



##### Save #############################################################


post        <- as.data.frame(fit_sem_jags$BUGSoutput$sims.list)
summary_tab <- as.data.frame(fit_sem_jags$BUGSoutput$summary)

# 收敛诊断
rhats     <- fit_sem_jags$BUGSoutput$summary[, "Rhat"]
bad_rhats <- rhats[rhats > 1.1]
if (length(bad_rhats) > 0) {
  cat("（Rhat > 1.1）:\n"); print(bad_rhats)
} else {
  cat("Rhat < 1.1）\n")
}


saveRDS(fit_sem_jags,
        "/data/beiyesi_CWD9_DBH5_SLA.rds")
write.csv(summary_tab,
          "/data/Beiyrsi_Summary_American_CWD9_DBH5_SLA.csv",
          row.names = TRUE)
saveRDS(post,
        "/data/post_American_CWD9_DBH5_SLA.rds")



##### Extract #####################################################



get_effect_by_bin <- function(post, param_prefix, path_label) {
  param_cols  <- grep(paste0("^", param_prefix, "\\."), names(post), value = TRUE)
  summary_mat <- t(apply(
    post[, param_cols, drop = FALSE], 2,
    function(x) c(mean  = mean(x),
                  lower = quantile(x, 0.025),
                  upper = quantile(x, 0.975))
  ))
  df           <- as.data.frame(summary_mat)
  colnames(df) <- c("mean", "lower", "upper")
  df$bin       <- as.numeric(sub(paste0("^", param_prefix, "\\."), "", param_cols))
  df$effect_type <- ifelse(
    df$lower > 0 & df$upper > 0, "positive significant",
    ifelse(df$lower < 0 & df$upper < 0, "negative significant",
           ifelse(df$mean > 0, "positive non-significant", "negative non-significant"))
  )
  df$path <- path_label
  return(df)
}


get_indirect_effect_by_bin <- function(post, delta_prefix, beta_prefix, path_label) {
  delta_cols  <- sort(grep(paste0("^", delta_prefix, "\\."), names(post), value = TRUE))
  beta_cols   <- sort(grep(paste0("^", beta_prefix,  "\\."), names(post), value = TRUE))
  sorted_bins <- as.numeric(sub(paste0("^", delta_prefix, "\\."), "", delta_cols))
  indirect    <- post[, delta_cols] * post[, beta_cols]
  summary_mat <- t(apply(indirect, 2,
                         function(x) c(mean  = mean(x),
                                       lower = quantile(x, 0.025),
                                       upper = quantile(x, 0.975))
  ))
  df           <- as.data.frame(summary_mat)
  colnames(df) <- c("mean", "lower", "upper")
  df$bin       <- sorted_bins
  df$effect_type <- ifelse(
    df$lower > 0 & df$upper > 0, "positive significant",
    ifelse(df$lower < 0 & df$upper < 0, "negative significant",
           ifelse(df$mean > 0, "positive non-significant", "negative non-significant"))
  )
  df$effect_type <- factor(df$effect_type, levels = c(
    "negative significant", "negative non-significant",
    "positive non-significant", "positive significant"))
  df$path <- path_label
  return(df)
}


all_effects <- bind_rows(
  get_effect_by_bin(post, "beta3",  "FD → AGB"),
  get_effect_by_bin(post, "beta4",  "CWM → AGB"),
  get_effect_by_bin(post, "beta5",  "DBHCV → AGB"),
  get_effect_by_bin(post, "delta1", "FD → DBHCV"),
  get_effect_by_bin(post, "delta2", "CWM → DBHCV"),
  get_indirect_effect_by_bin(post, "delta1", "beta5", "FD → DBHCV → AGB"),
  get_indirect_effect_by_bin(post, "delta2", "beta5", "CWM → DBHCV → AGB")
)


data_sem2_org     <- data_sem2
data_sem2_org$bin <- as.numeric(factor(data_sem2_org$PDSI_ATA_bin))

result <- data_sem2_org %>%
  group_by(bin) %>%
  dplyr::summarise(
    mean_MAT    = mean(mean_annual_temp,   na.rm = TRUE),
    mean_PRE    = mean(mean_annual_precip, na.rm = TRUE),
    mean_CWD    = mean(CWD1,           na.rm = TRUE),  
    mean_STDAGE = mean(STDAGE,             na.rm = TRUE),
    .groups = "drop"
  )

all_effects1 <- merge(all_effects, result, by = "bin")

write.csv(all_effects1,
          "/data/posterior_American_CWD9_DBH5_SLA.csv",
          row.names = FALSE)



##### Proportional stacked bar chart #################################################################


library(scales)

get_effect_prop <- function(post, param_prefix, path_label) {
  param_cols <- grep(paste0("^", param_prefix, "\\."), names(post), value = TRUE)
  df <- t(apply(post[, param_cols, drop = FALSE], 2,
                function(x) c(mean  = mean(x),
                              lower = quantile(x, 0.025),
                              upper = quantile(x, 0.975))
  )) |> as.data.frame()
  colnames(df) <- c("mean", "lower", "upper")
  df$effect_type <- ifelse(
    df$lower > 0 & df$upper > 0, "positive significant",
    ifelse(df$lower < 0 & df$upper < 0, "negative significant",
           ifelse(df$mean > 0, "positive non-significant", "negative non-significant"))
  )
  dplyr::count(df, effect_type) %>%
    mutate(percent = n / sum(n) * 100, path = path_label)
}

get_indirect_prop <- function(post, delta_prefix, beta_prefix, path_label) {
  delta_cols <- sort(grep(paste0("^", delta_prefix, "\\."), names(post), value = TRUE))
  beta_cols  <- sort(grep(paste0("^", beta_prefix,  "\\."), names(post), value = TRUE))
  indirect   <- post[, delta_cols] * post[, beta_cols]
  df <- t(apply(indirect, 2,
                function(x) c(mean  = mean(x),
                              lower = quantile(x, 0.025),
                              upper = quantile(x, 0.975))
  )) |> as.data.frame()
  colnames(df) <- c("mean", "lower", "upper")
  df$effect_type <- ifelse(
    df$lower > 0 & df$upper > 0, "positive significant",
    ifelse(df$lower < 0 & df$upper < 0, "negative significant",
           ifelse(df$mean > 0, "positive non-significant", "negative non-significant"))
  )
  dplyr::count(df, effect_type) %>%
    mutate(percent = n / sum(n) * 100, path = path_label)
}

prop <- bind_rows(
  get_effect_prop(post, "beta3",  "FD → AGB"),
  get_effect_prop(post, "beta4",  "CWM → AGB"),
  get_effect_prop(post, "beta5",  "DBHCV → AGB"),
  get_effect_prop(post, "delta1", "FD → DBHCV"),
  get_effect_prop(post, "delta2", "CWM → DBHCV"),
  get_indirect_prop(post, "delta1", "beta5", "FD → DBHCV → AGB"),
  get_indirect_prop(post, "delta2", "beta5", "CWM → DBHCV → AGB")
)

write.csv(prop,
          "/data/bar_American_CWD9_DBH5_SLA.csv",
          row.names = FALSE)


prop_plot <- prop %>%
  filter(path %in% c("FD → AGB", "FD → DBHCV → AGB",
                     "CWM → AGB", "CWM → DBHCV → AGB")) %>%
  mutate(
    path = factor(path, levels = c(
      "FD → AGB", "FD → DBHCV → AGB",
      "CWM → AGB", "CWM → DBHCV → AGB"
    )),
    effect_type = factor(effect_type, levels = c(
      "negative significant", "negative non-significant",
      "positive non-significant", "positive significant"
    ))
  )

fig_height <- ggplot(prop_plot,
                     aes(x = percent, y = path, fill = effect_type)) +
  geom_bar(stat = "identity", width = 0.65) +
  scale_x_continuous(labels = percent_format(scale = 1), expand = c(0, 0)) +
  scale_fill_manual(
    values = c(
      "negative significant"     = "#92C5DE",
      "negative non-significant" = "#dbebfa",
      "positive non-significant" = "#f9ebdf",
      "positive significant"     = "#F4A582"
    ),
    drop = FALSE
  ) +
  labs(x = "Proportion of bins (%)", y = NULL, fill = "Effect type") +
  ggtitle("Structural diversity → Carbon") +
  theme_classic(base_size = 13) +
  theme(axis.text.y  = element_text(size = 12),
        legend.position = "right")

fig_height
##===============================================================================
#####Leaf.longevity###################################################################
##===============================================================================
##### Load data ###################################################################

library(data.table)
library(dplyr)
library(R2jags)

C_datatotal_all <- fread("/data/Data.csv" )
colnames( C_datatotal_all  )

C_datatotal_all <- C_datatotal_all %>%
  dplyr::select(-V1) %>%
  filter(Carbon_Mg_ha > 0)
colnames(C_datatotal_all) <- gsub("_mean$", "", colnames(C_datatotal_all))
C_datatotal_all <- C_datatotal_all %>%
  filter(!is.na(CWD_stand_age)) %>%
  filter(!is.na(MAT_stand_age)) %>%
  filter(!is.na(MAP_stand_age)) 

# 合并 bin 分组
C_bin <- read.csv("/data/Class_CWD_9bins")
C_bin <- C_bin[, c("PLT_CN", "PDSI_ATA_bin")]
C_datatotal_all$PLT_CN <- as.character(C_datatotal_all$PLT_CN)
C_bin$PLT_CN            <- as.character(C_bin$PLT_CN)
C_datatotal_all <- merge(C_datatotal_all, C_bin, by = "PLT_CN")
C_datatotal_all$CWD1=C_datatotal_all$CWD_stand_age
# 标准化
C_datatotal_all <- C_datatotal_all %>%
  mutate(
    Carbon_Mg_ha = log(Carbon_Mg_ha+1),
    shannon_dbhR = scale(H_D_5),
           
           
    CWD          = scale(CWD_stand_age),  
    TEM1         = scale(MAT_stand_age),
    PRE1         = scale(MAP_stand_age),
    FDis         = scale(Leaf.longevity.FDis),
    CWM          = scale(Leaf.longevity.CWM),
    pH           = scale(nitrogen),
    bdod         = scale(bdod),
    STDAGE       = scale(STDAGE)
  )

data_sem2 <- C_datatotal_all

##### JAGS #########################################################


data_all <- data_sem2

y_CV      <- as.numeric(data_all$shannon_dbhR)
y_AGB     <- as.numeric(data_all$Carbon_Mg_ha)
trait     <- as.numeric(data_all$FDis)
trait_CWM <- as.numeric(data_all$CWM)
CWD       <- as.numeric(data_all$CWD)       
TEM1      <- as.numeric(data_all$TEM1)
PRE1      <- as.numeric(data_all$PRE1)
pH        <- as.numeric(data_all$pH)
bdod      <- as.numeric(data_all$bdod)
STDAGE    <- as.numeric(data_all$STDAGE)

bin   <- as.numeric(factor(data_all$PDSI_ATA_bin))
N     <- nrow(data_all)
N_bin <- max(bin)

cat("N =", N, "| N_bin =", N_bin, "\n")


FD_bin     <- tapply(trait,     bin, mean, na.rm = TRUE)
FI_bin     <- tapply(trait_CWM, bin, mean, na.rm = TRUE)
CWD_bin    <- tapply(CWD,       bin, mean, na.rm = TRUE)  
TEM1_bin   <- tapply(TEM1,      bin, mean, na.rm = TRUE)
PRE1_bin   <- tapply(PRE1,      bin, mean, na.rm = TRUE)
pH_bin     <- tapply(pH,        bin, mean, na.rm = TRUE)
bdod_bin   <- tapply(bdod,      bin, mean, na.rm = TRUE)
STDAGE_bin <- tapply(STDAGE,    bin, mean, na.rm = TRUE)



print(sapply(
  list(FD=FD_bin, FI=FI_bin, CWD=CWD_bin,
       TEM1=TEM1_bin, PRE1=PRE1_bin,
       pH=pH_bin, bdod=bdod_bin, STDAGE=STDAGE_bin),
  function(x) sum(is.na(x))
))


jags_data <- list(
  N            = N,
  N_bin        = N_bin,
  DBH_CV       = y_CV,
  Carbon_Mg_ha = y_AGB,
  trait        = trait,
  trait_CWM    = trait_CWM,
  CWD          = CWD,                        
  TEM1         = TEM1,
  PRE1         = PRE1,
  pH           = pH,
  bdod         = bdod,
  STDAGE       = STDAGE,
  bin          = bin,
  FD_bin       = as.numeric(FD_bin),
  FI_bin       = as.numeric(FI_bin),
  CWD_bin      = as.numeric(CWD_bin),        
  TEM1_bin     = as.numeric(TEM1_bin),
  PRE1_bin     = as.numeric(PRE1_bin),
  pH_bin       = as.numeric(pH_bin),
  bdod_bin     = as.numeric(bdod_bin),
  STDAGE_bin   = as.numeric(STDAGE_bin)
)



##### parameters##########################################


params <- c(
  "delta1", "delta2",
  "beta3",  "beta4",  "beta5",
  "gamma0_CV",
  "mu_alpha_AGB",
  "gamma_CV",
  "gamma_CV_bin",
  "gamma_AGB",
  "beta5_0",
  "beta5_1", "beta5_2",
  "beta5_3", "beta5_4",
  "beta5_5", "beta5_6",
  "beta5_7", "beta5_8",
  "sigma_CV",       "sigma_AGB",
  "sigma_alpha_CV", "sigma_alpha_AGB",
  "sigma_delta1",   "sigma_delta2",
  "sigma_beta3",    "sigma_beta4",    "sigma_beta5"
)



##### Run JAGS #################################################################


fit_sem_jags <- jags(
  model.file         = "/data/bayes_sem_model.txt",
  data               = jags_data,
  parameters.to.save = params,
  n.chains           = 3,
  n.iter             = 10000,
  n.burnin           = 2000,
  n.thin             = 10
)



##### Save #############################################################


post        <- as.data.frame(fit_sem_jags$BUGSoutput$sims.list)
summary_tab <- as.data.frame(fit_sem_jags$BUGSoutput$summary)

# 收敛诊断
rhats     <- fit_sem_jags$BUGSoutput$summary[, "Rhat"]
bad_rhats <- rhats[rhats > 1.1]
if (length(bad_rhats) > 0) {
  cat("（Rhat > 1.1）:\n"); print(bad_rhats)
} else {
  cat("Rhat < 1.1）\n")
}


saveRDS(fit_sem_jags,
        "/data/beiyesi_CWD9_DBH5_Leaf.longevity.rds")
write.csv(summary_tab,
          "/data/Beiyrsi_Summary_American_CWD9_DBH5_Leaf.longevity.csv",
          row.names = TRUE)
saveRDS(post,
        "/data/post_American_CWD9_DBH5_Leaf.longevity.rds")



##### Extract #####################################################



get_effect_by_bin <- function(post, param_prefix, path_label) {
  param_cols  <- grep(paste0("^", param_prefix, "\\."), names(post), value = TRUE)
  summary_mat <- t(apply(
    post[, param_cols, drop = FALSE], 2,
    function(x) c(mean  = mean(x),
                  lower = quantile(x, 0.025),
                  upper = quantile(x, 0.975))
  ))
  df           <- as.data.frame(summary_mat)
  colnames(df) <- c("mean", "lower", "upper")
  df$bin       <- as.numeric(sub(paste0("^", param_prefix, "\\."), "", param_cols))
  df$effect_type <- ifelse(
    df$lower > 0 & df$upper > 0, "positive significant",
    ifelse(df$lower < 0 & df$upper < 0, "negative significant",
           ifelse(df$mean > 0, "positive non-significant", "negative non-significant"))
  )
  df$path <- path_label
  return(df)
}


get_indirect_effect_by_bin <- function(post, delta_prefix, beta_prefix, path_label) {
  delta_cols  <- sort(grep(paste0("^", delta_prefix, "\\."), names(post), value = TRUE))
  beta_cols   <- sort(grep(paste0("^", beta_prefix,  "\\."), names(post), value = TRUE))
  sorted_bins <- as.numeric(sub(paste0("^", delta_prefix, "\\."), "", delta_cols))
  indirect    <- post[, delta_cols] * post[, beta_cols]
  summary_mat <- t(apply(indirect, 2,
                         function(x) c(mean  = mean(x),
                                       lower = quantile(x, 0.025),
                                       upper = quantile(x, 0.975))
  ))
  df           <- as.data.frame(summary_mat)
  colnames(df) <- c("mean", "lower", "upper")
  df$bin       <- sorted_bins
  df$effect_type <- ifelse(
    df$lower > 0 & df$upper > 0, "positive significant",
    ifelse(df$lower < 0 & df$upper < 0, "negative significant",
           ifelse(df$mean > 0, "positive non-significant", "negative non-significant"))
  )
  df$effect_type <- factor(df$effect_type, levels = c(
    "negative significant", "negative non-significant",
    "positive non-significant", "positive significant"))
  df$path <- path_label
  return(df)
}


all_effects <- bind_rows(
  get_effect_by_bin(post, "beta3",  "FD → AGB"),
  get_effect_by_bin(post, "beta4",  "CWM → AGB"),
  get_effect_by_bin(post, "beta5",  "DBHCV → AGB"),
  get_effect_by_bin(post, "delta1", "FD → DBHCV"),
  get_effect_by_bin(post, "delta2", "CWM → DBHCV"),
  get_indirect_effect_by_bin(post, "delta1", "beta5", "FD → DBHCV → AGB"),
  get_indirect_effect_by_bin(post, "delta2", "beta5", "CWM → DBHCV → AGB")
)


data_sem2_org     <- data_sem2
data_sem2_org$bin <- as.numeric(factor(data_sem2_org$PDSI_ATA_bin))

result <- data_sem2_org %>%
  group_by(bin) %>%
  dplyr::summarise(
    mean_MAT    = mean(mean_annual_temp,   na.rm = TRUE),
    mean_PRE    = mean(mean_annual_precip, na.rm = TRUE),
    mean_CWD    = mean(CWD1,           na.rm = TRUE),  
    mean_STDAGE = mean(STDAGE,             na.rm = TRUE),
    .groups = "drop"
  )

all_effects1 <- merge(all_effects, result, by = "bin")

write.csv(all_effects1,
          "/data/posterior_American_CWD9_DBH5_Leaf.longevity.csv",
          row.names = FALSE)



##### Proportional stacked bar chart #################################################################


library(scales)

get_effect_prop <- function(post, param_prefix, path_label) {
  param_cols <- grep(paste0("^", param_prefix, "\\."), names(post), value = TRUE)
  df <- t(apply(post[, param_cols, drop = FALSE], 2,
                function(x) c(mean  = mean(x),
                              lower = quantile(x, 0.025),
                              upper = quantile(x, 0.975))
  )) |> as.data.frame()
  colnames(df) <- c("mean", "lower", "upper")
  df$effect_type <- ifelse(
    df$lower > 0 & df$upper > 0, "positive significant",
    ifelse(df$lower < 0 & df$upper < 0, "negative significant",
           ifelse(df$mean > 0, "positive non-significant", "negative non-significant"))
  )
  dplyr::count(df, effect_type) %>%
    mutate(percent = n / sum(n) * 100, path = path_label)
}

get_indirect_prop <- function(post, delta_prefix, beta_prefix, path_label) {
  delta_cols <- sort(grep(paste0("^", delta_prefix, "\\."), names(post), value = TRUE))
  beta_cols  <- sort(grep(paste0("^", beta_prefix,  "\\."), names(post), value = TRUE))
  indirect   <- post[, delta_cols] * post[, beta_cols]
  df <- t(apply(indirect, 2,
                function(x) c(mean  = mean(x),
                              lower = quantile(x, 0.025),
                              upper = quantile(x, 0.975))
  )) |> as.data.frame()
  colnames(df) <- c("mean", "lower", "upper")
  df$effect_type <- ifelse(
    df$lower > 0 & df$upper > 0, "positive significant",
    ifelse(df$lower < 0 & df$upper < 0, "negative significant",
           ifelse(df$mean > 0, "positive non-significant", "negative non-significant"))
  )
  dplyr::count(df, effect_type) %>%
    mutate(percent = n / sum(n) * 100, path = path_label)
}

prop <- bind_rows(
  get_effect_prop(post, "beta3",  "FD → AGB"),
  get_effect_prop(post, "beta4",  "CWM → AGB"),
  get_effect_prop(post, "beta5",  "DBHCV → AGB"),
  get_effect_prop(post, "delta1", "FD → DBHCV"),
  get_effect_prop(post, "delta2", "CWM → DBHCV"),
  get_indirect_prop(post, "delta1", "beta5", "FD → DBHCV → AGB"),
  get_indirect_prop(post, "delta2", "beta5", "CWM → DBHCV → AGB")
)

write.csv(prop,
          "/data/bar_American_CWD9_DBH5_Leaf.longevity.csv",
          row.names = FALSE)


prop_plot <- prop %>%
  filter(path %in% c("FD → AGB", "FD → DBHCV → AGB",
                     "CWM → AGB", "CWM → DBHCV → AGB")) %>%
  mutate(
    path = factor(path, levels = c(
      "FD → AGB", "FD → DBHCV → AGB",
      "CWM → AGB", "CWM → DBHCV → AGB"
    )),
    effect_type = factor(effect_type, levels = c(
      "negative significant", "negative non-significant",
      "positive non-significant", "positive significant"
    ))
  )

fig_height <- ggplot(prop_plot,
                     aes(x = percent, y = path, fill = effect_type)) +
  geom_bar(stat = "identity", width = 0.65) +
  scale_x_continuous(labels = percent_format(scale = 1), expand = c(0, 0)) +
  scale_fill_manual(
    values = c(
      "negative significant"     = "#92C5DE",
      "negative non-significant" = "#dbebfa",
      "positive non-significant" = "#f9ebdf",
      "positive significant"     = "#F4A582"
    ),
    drop = FALSE
  ) +
  labs(x = "Proportion of bins (%)", y = NULL, fill = "Effect type") +
  ggtitle("Structural diversity → Carbon") +
  theme_classic(base_size = 13) +
  theme(axis.text.y  = element_text(size = 12),
        legend.position = "right")

fig_height
##===============================================================================
#####WD###################################################################
##===============================================================================
##### Load data ###################################################################

library(data.table)
library(dplyr)
library(R2jags)

C_datatotal_all <- fread("/data/Data.csv" )
colnames( C_datatotal_all  )

C_datatotal_all <- C_datatotal_all %>%
  dplyr::select(-V1) %>%
  filter(Carbon_Mg_ha > 0)
colnames(C_datatotal_all) <- gsub("_mean$", "", colnames(C_datatotal_all))
C_datatotal_all <- C_datatotal_all %>%
  filter(!is.na(CWD_stand_age)) %>%
  filter(!is.na(MAT_stand_age)) %>%
  filter(!is.na(MAP_stand_age)) 

# 合并 bin 分组
C_bin <- read.csv("/data/Class_CWD_9bins")
C_bin <- C_bin[, c("PLT_CN", "PDSI_ATA_bin")]
C_datatotal_all$PLT_CN <- as.character(C_datatotal_all$PLT_CN)
C_bin$PLT_CN            <- as.character(C_bin$PLT_CN)
C_datatotal_all <- merge(C_datatotal_all, C_bin, by = "PLT_CN")
C_datatotal_all$CWD1=C_datatotal_all$CWD_stand_age
# 标准化
C_datatotal_all <- C_datatotal_all %>%
  mutate(
    Carbon_Mg_ha = log(Carbon_Mg_ha+1),
    shannon_dbhR = scale(H_D_5),
           
           
    CWD          = scale(CWD_stand_age),  
    TEM1         = scale(MAT_stand_age),
    PRE1         = scale(MAP_stand_age),
    FDis         = scale(WD.FDis),
    CWM          = scale(WD.CWM),
    pH           = scale(nitrogen),
    bdod         = scale(bdod),
    STDAGE       = scale(STDAGE)
  )

data_sem2 <- C_datatotal_all

##### JAGS #########################################################


data_all <- data_sem2

y_CV      <- as.numeric(data_all$shannon_dbhR)
y_AGB     <- as.numeric(data_all$Carbon_Mg_ha)
trait     <- as.numeric(data_all$FDis)
trait_CWM <- as.numeric(data_all$CWM)
CWD       <- as.numeric(data_all$CWD)       
TEM1      <- as.numeric(data_all$TEM1)
PRE1      <- as.numeric(data_all$PRE1)
pH        <- as.numeric(data_all$pH)
bdod      <- as.numeric(data_all$bdod)
STDAGE    <- as.numeric(data_all$STDAGE)

bin   <- as.numeric(factor(data_all$PDSI_ATA_bin))
N     <- nrow(data_all)
N_bin <- max(bin)

cat("N =", N, "| N_bin =", N_bin, "\n")


FD_bin     <- tapply(trait,     bin, mean, na.rm = TRUE)
FI_bin     <- tapply(trait_CWM, bin, mean, na.rm = TRUE)
CWD_bin    <- tapply(CWD,       bin, mean, na.rm = TRUE)  
TEM1_bin   <- tapply(TEM1,      bin, mean, na.rm = TRUE)
PRE1_bin   <- tapply(PRE1,      bin, mean, na.rm = TRUE)
pH_bin     <- tapply(pH,        bin, mean, na.rm = TRUE)
bdod_bin   <- tapply(bdod,      bin, mean, na.rm = TRUE)
STDAGE_bin <- tapply(STDAGE,    bin, mean, na.rm = TRUE)



print(sapply(
  list(FD=FD_bin, FI=FI_bin, CWD=CWD_bin,
       TEM1=TEM1_bin, PRE1=PRE1_bin,
       pH=pH_bin, bdod=bdod_bin, STDAGE=STDAGE_bin),
  function(x) sum(is.na(x))
))


jags_data <- list(
  N            = N,
  N_bin        = N_bin,
  DBH_CV       = y_CV,
  Carbon_Mg_ha = y_AGB,
  trait        = trait,
  trait_CWM    = trait_CWM,
  CWD          = CWD,                        
  TEM1         = TEM1,
  PRE1         = PRE1,
  pH           = pH,
  bdod         = bdod,
  STDAGE       = STDAGE,
  bin          = bin,
  FD_bin       = as.numeric(FD_bin),
  FI_bin       = as.numeric(FI_bin),
  CWD_bin      = as.numeric(CWD_bin),        
  TEM1_bin     = as.numeric(TEM1_bin),
  PRE1_bin     = as.numeric(PRE1_bin),
  pH_bin       = as.numeric(pH_bin),
  bdod_bin     = as.numeric(bdod_bin),
  STDAGE_bin   = as.numeric(STDAGE_bin)
)



##### parameters##########################################


params <- c(
  "delta1", "delta2",
  "beta3",  "beta4",  "beta5",
  "gamma0_CV",
  "mu_alpha_AGB",
  "gamma_CV",
  "gamma_CV_bin",
  "gamma_AGB",
  "beta5_0",
  "beta5_1", "beta5_2",
  "beta5_3", "beta5_4",
  "beta5_5", "beta5_6",
  "beta5_7", "beta5_8",
  "sigma_CV",       "sigma_AGB",
  "sigma_alpha_CV", "sigma_alpha_AGB",
  "sigma_delta1",   "sigma_delta2",
  "sigma_beta3",    "sigma_beta4",    "sigma_beta5"
)



##### Run JAGS #################################################################


fit_sem_jags <- jags(
  model.file         = "/data/bayes_sem_model.txt",
  data               = jags_data,
  parameters.to.save = params,
  n.chains           = 3,
  n.iter             = 10000,
  n.burnin           = 2000,
  n.thin             = 10
)



##### Save #############################################################


post        <- as.data.frame(fit_sem_jags$BUGSoutput$sims.list)
summary_tab <- as.data.frame(fit_sem_jags$BUGSoutput$summary)

# 收敛诊断
rhats     <- fit_sem_jags$BUGSoutput$summary[, "Rhat"]
bad_rhats <- rhats[rhats > 1.1]
if (length(bad_rhats) > 0) {
  cat("（Rhat > 1.1）:\n"); print(bad_rhats)
} else {
  cat("Rhat < 1.1）\n")
}


saveRDS(fit_sem_jags,
        "/data/beiyesi_CWD9_DBH5_WD.rds")
write.csv(summary_tab,
          "/data/Beiyrsi_Summary_American_CWD9_DBH5_WD.csv",
          row.names = TRUE)
saveRDS(post,
        "/data/post_American_CWD9_DBH5_WD.rds")



##### Extract #####################################################



get_effect_by_bin <- function(post, param_prefix, path_label) {
  param_cols  <- grep(paste0("^", param_prefix, "\\."), names(post), value = TRUE)
  summary_mat <- t(apply(
    post[, param_cols, drop = FALSE], 2,
    function(x) c(mean  = mean(x),
                  lower = quantile(x, 0.025),
                  upper = quantile(x, 0.975))
  ))
  df           <- as.data.frame(summary_mat)
  colnames(df) <- c("mean", "lower", "upper")
  df$bin       <- as.numeric(sub(paste0("^", param_prefix, "\\."), "", param_cols))
  df$effect_type <- ifelse(
    df$lower > 0 & df$upper > 0, "positive significant",
    ifelse(df$lower < 0 & df$upper < 0, "negative significant",
           ifelse(df$mean > 0, "positive non-significant", "negative non-significant"))
  )
  df$path <- path_label
  return(df)
}


get_indirect_effect_by_bin <- function(post, delta_prefix, beta_prefix, path_label) {
  delta_cols  <- sort(grep(paste0("^", delta_prefix, "\\."), names(post), value = TRUE))
  beta_cols   <- sort(grep(paste0("^", beta_prefix,  "\\."), names(post), value = TRUE))
  sorted_bins <- as.numeric(sub(paste0("^", delta_prefix, "\\."), "", delta_cols))
  indirect    <- post[, delta_cols] * post[, beta_cols]
  summary_mat <- t(apply(indirect, 2,
                         function(x) c(mean  = mean(x),
                                       lower = quantile(x, 0.025),
                                       upper = quantile(x, 0.975))
  ))
  df           <- as.data.frame(summary_mat)
  colnames(df) <- c("mean", "lower", "upper")
  df$bin       <- sorted_bins
  df$effect_type <- ifelse(
    df$lower > 0 & df$upper > 0, "positive significant",
    ifelse(df$lower < 0 & df$upper < 0, "negative significant",
           ifelse(df$mean > 0, "positive non-significant", "negative non-significant"))
  )
  df$effect_type <- factor(df$effect_type, levels = c(
    "negative significant", "negative non-significant",
    "positive non-significant", "positive significant"))
  df$path <- path_label
  return(df)
}


all_effects <- bind_rows(
  get_effect_by_bin(post, "beta3",  "FD → AGB"),
  get_effect_by_bin(post, "beta4",  "CWM → AGB"),
  get_effect_by_bin(post, "beta5",  "DBHCV → AGB"),
  get_effect_by_bin(post, "delta1", "FD → DBHCV"),
  get_effect_by_bin(post, "delta2", "CWM → DBHCV"),
  get_indirect_effect_by_bin(post, "delta1", "beta5", "FD → DBHCV → AGB"),
  get_indirect_effect_by_bin(post, "delta2", "beta5", "CWM → DBHCV → AGB")
)


data_sem2_org     <- data_sem2
data_sem2_org$bin <- as.numeric(factor(data_sem2_org$PDSI_ATA_bin))

result <- data_sem2_org %>%
  group_by(bin) %>%
  dplyr::summarise(
    mean_MAT    = mean(mean_annual_temp,   na.rm = TRUE),
    mean_PRE    = mean(mean_annual_precip, na.rm = TRUE),
    mean_CWD    = mean(CWD1,           na.rm = TRUE),  
    mean_STDAGE = mean(STDAGE,             na.rm = TRUE),
    .groups = "drop"
  )

all_effects1 <- merge(all_effects, result, by = "bin")

write.csv(all_effects1,
          "/data/posterior_American_CWD9_DBH5_WD.csv",
          row.names = FALSE)



##### Proportional stacked bar chart #################################################################


library(scales)

get_effect_prop <- function(post, param_prefix, path_label) {
  param_cols <- grep(paste0("^", param_prefix, "\\."), names(post), value = TRUE)
  df <- t(apply(post[, param_cols, drop = FALSE], 2,
                function(x) c(mean  = mean(x),
                              lower = quantile(x, 0.025),
                              upper = quantile(x, 0.975))
  )) |> as.data.frame()
  colnames(df) <- c("mean", "lower", "upper")
  df$effect_type <- ifelse(
    df$lower > 0 & df$upper > 0, "positive significant",
    ifelse(df$lower < 0 & df$upper < 0, "negative significant",
           ifelse(df$mean > 0, "positive non-significant", "negative non-significant"))
  )
  dplyr::count(df, effect_type) %>%
    mutate(percent = n / sum(n) * 100, path = path_label)
}

get_indirect_prop <- function(post, delta_prefix, beta_prefix, path_label) {
  delta_cols <- sort(grep(paste0("^", delta_prefix, "\\."), names(post), value = TRUE))
  beta_cols  <- sort(grep(paste0("^", beta_prefix,  "\\."), names(post), value = TRUE))
  indirect   <- post[, delta_cols] * post[, beta_cols]
  df <- t(apply(indirect, 2,
                function(x) c(mean  = mean(x),
                              lower = quantile(x, 0.025),
                              upper = quantile(x, 0.975))
  )) |> as.data.frame()
  colnames(df) <- c("mean", "lower", "upper")
  df$effect_type <- ifelse(
    df$lower > 0 & df$upper > 0, "positive significant",
    ifelse(df$lower < 0 & df$upper < 0, "negative significant",
           ifelse(df$mean > 0, "positive non-significant", "negative non-significant"))
  )
  dplyr::count(df, effect_type) %>%
    mutate(percent = n / sum(n) * 100, path = path_label)
}

prop <- bind_rows(
  get_effect_prop(post, "beta3",  "FD → AGB"),
  get_effect_prop(post, "beta4",  "CWM → AGB"),
  get_effect_prop(post, "beta5",  "DBHCV → AGB"),
  get_effect_prop(post, "delta1", "FD → DBHCV"),
  get_effect_prop(post, "delta2", "CWM → DBHCV"),
  get_indirect_prop(post, "delta1", "beta5", "FD → DBHCV → AGB"),
  get_indirect_prop(post, "delta2", "beta5", "CWM → DBHCV → AGB")
)

write.csv(prop,
          "/data/bar_American_CWD9_DBH5_WD.csv",
          row.names = FALSE)


prop_plot <- prop %>%
  filter(path %in% c("FD → AGB", "FD → DBHCV → AGB",
                     "CWM → AGB", "CWM → DBHCV → AGB")) %>%
  mutate(
    path = factor(path, levels = c(
      "FD → AGB", "FD → DBHCV → AGB",
      "CWM → AGB", "CWM → DBHCV → AGB"
    )),
    effect_type = factor(effect_type, levels = c(
      "negative significant", "negative non-significant",
      "positive non-significant", "positive significant"
    ))
  )

fig_height <- ggplot(prop_plot,
                     aes(x = percent, y = path, fill = effect_type)) +
  geom_bar(stat = "identity", width = 0.65) +
  scale_x_continuous(labels = percent_format(scale = 1), expand = c(0, 0)) +
  scale_fill_manual(
    values = c(
      "negative significant"     = "#92C5DE",
      "negative non-significant" = "#dbebfa",
      "positive non-significant" = "#f9ebdf",
      "positive significant"     = "#F4A582"
    ),
    drop = FALSE
  ) +
  labs(x = "Proportion of bins (%)", y = NULL, fill = "Effect type") +
  ggtitle("Structural diversity → Carbon") +
  theme_classic(base_size = 13) +
  theme(axis.text.y  = element_text(size = 12),
        legend.position = "right")

fig_height
##===============================================================================
#####Plant_height###################################################################
##===============================================================================
##### Load data ###################################################################

library(data.table)
library(dplyr)
library(R2jags)

C_datatotal_all <- fread("/data/Data.csv" )
colnames( C_datatotal_all  )

C_datatotal_all <- C_datatotal_all %>%
  dplyr::select(-V1) %>%
  filter(Carbon_Mg_ha > 0)
colnames(C_datatotal_all) <- gsub("_mean$", "", colnames(C_datatotal_all))
C_datatotal_all <- C_datatotal_all %>%
  filter(!is.na(CWD_stand_age)) %>%
  filter(!is.na(MAT_stand_age)) %>%
  filter(!is.na(MAP_stand_age)) 

# 合并 bin 分组
C_bin <- read.csv("/data/Class_CWD_9bins")
C_bin <- C_bin[, c("PLT_CN", "PDSI_ATA_bin")]
C_datatotal_all$PLT_CN <- as.character(C_datatotal_all$PLT_CN)
C_bin$PLT_CN            <- as.character(C_bin$PLT_CN)
C_datatotal_all <- merge(C_datatotal_all, C_bin, by = "PLT_CN")
C_datatotal_all$CWD1=C_datatotal_all$CWD_stand_age
# 标准化
C_datatotal_all <- C_datatotal_all %>%
  mutate(
    Carbon_Mg_ha = log(Carbon_Mg_ha+1),
    shannon_dbhR = scale(H_D_5),
           
           
    CWD          = scale(CWD_stand_age),  
    TEM1         = scale(MAT_stand_age),
    PRE1         = scale(MAP_stand_age),
    FDis         = scale(Plant_height.FDis),
    CWM          = scale(Plant_height.CWM),
    pH           = scale(nitrogen),
    bdod         = scale(bdod),
    STDAGE       = scale(STDAGE)
  )

data_sem2 <- C_datatotal_all

##### JAGS #########################################################


data_all <- data_sem2

y_CV      <- as.numeric(data_all$shannon_dbhR)
y_AGB     <- as.numeric(data_all$Carbon_Mg_ha)
trait     <- as.numeric(data_all$FDis)
trait_CWM <- as.numeric(data_all$CWM)
CWD       <- as.numeric(data_all$CWD)       
TEM1      <- as.numeric(data_all$TEM1)
PRE1      <- as.numeric(data_all$PRE1)
pH        <- as.numeric(data_all$pH)
bdod      <- as.numeric(data_all$bdod)
STDAGE    <- as.numeric(data_all$STDAGE)

bin   <- as.numeric(factor(data_all$PDSI_ATA_bin))
N     <- nrow(data_all)
N_bin <- max(bin)

cat("N =", N, "| N_bin =", N_bin, "\n")


FD_bin     <- tapply(trait,     bin, mean, na.rm = TRUE)
FI_bin     <- tapply(trait_CWM, bin, mean, na.rm = TRUE)
CWD_bin    <- tapply(CWD,       bin, mean, na.rm = TRUE)  
TEM1_bin   <- tapply(TEM1,      bin, mean, na.rm = TRUE)
PRE1_bin   <- tapply(PRE1,      bin, mean, na.rm = TRUE)
pH_bin     <- tapply(pH,        bin, mean, na.rm = TRUE)
bdod_bin   <- tapply(bdod,      bin, mean, na.rm = TRUE)
STDAGE_bin <- tapply(STDAGE,    bin, mean, na.rm = TRUE)



print(sapply(
  list(FD=FD_bin, FI=FI_bin, CWD=CWD_bin,
       TEM1=TEM1_bin, PRE1=PRE1_bin,
       pH=pH_bin, bdod=bdod_bin, STDAGE=STDAGE_bin),
  function(x) sum(is.na(x))
))


jags_data <- list(
  N            = N,
  N_bin        = N_bin,
  DBH_CV       = y_CV,
  Carbon_Mg_ha = y_AGB,
  trait        = trait,
  trait_CWM    = trait_CWM,
  CWD          = CWD,                        
  TEM1         = TEM1,
  PRE1         = PRE1,
  pH           = pH,
  bdod         = bdod,
  STDAGE       = STDAGE,
  bin          = bin,
  FD_bin       = as.numeric(FD_bin),
  FI_bin       = as.numeric(FI_bin),
  CWD_bin      = as.numeric(CWD_bin),        
  TEM1_bin     = as.numeric(TEM1_bin),
  PRE1_bin     = as.numeric(PRE1_bin),
  pH_bin       = as.numeric(pH_bin),
  bdod_bin     = as.numeric(bdod_bin),
  STDAGE_bin   = as.numeric(STDAGE_bin)
)



##### parameters##########################################


params <- c(
  "delta1", "delta2",
  "beta3",  "beta4",  "beta5",
  "gamma0_CV",
  "mu_alpha_AGB",
  "gamma_CV",
  "gamma_CV_bin",
  "gamma_AGB",
  "beta5_0",
  "beta5_1", "beta5_2",
  "beta5_3", "beta5_4",
  "beta5_5", "beta5_6",
  "beta5_7", "beta5_8",
  "sigma_CV",       "sigma_AGB",
  "sigma_alpha_CV", "sigma_alpha_AGB",
  "sigma_delta1",   "sigma_delta2",
  "sigma_beta3",    "sigma_beta4",    "sigma_beta5"
)



##### Run JAGS #################################################################


fit_sem_jags <- jags(
  model.file         = "/data/bayes_sem_model.txt",
  data               = jags_data,
  parameters.to.save = params,
  n.chains           = 3,
  n.iter             = 10000,
  n.burnin           = 2000,
  n.thin             = 10
)



##### Save #############################################################


post        <- as.data.frame(fit_sem_jags$BUGSoutput$sims.list)
summary_tab <- as.data.frame(fit_sem_jags$BUGSoutput$summary)

# 收敛诊断
rhats     <- fit_sem_jags$BUGSoutput$summary[, "Rhat"]
bad_rhats <- rhats[rhats > 1.1]
if (length(bad_rhats) > 0) {
  cat("（Rhat > 1.1）:\n"); print(bad_rhats)
} else {
  cat("Rhat < 1.1）\n")
}


saveRDS(fit_sem_jags,
        "/data/beiyesi_CWD9_DBH5_Plant_height.rds")
write.csv(summary_tab,
          "/data/Beiyrsi_Summary_American_CWD9_DBH5_Plant_height.csv",
          row.names = TRUE)
saveRDS(post,
        "/data/post_American_CWD9_DBH5_Plant_height.rds")



##### Extract #####################################################



get_effect_by_bin <- function(post, param_prefix, path_label) {
  param_cols  <- grep(paste0("^", param_prefix, "\\."), names(post), value = TRUE)
  summary_mat <- t(apply(
    post[, param_cols, drop = FALSE], 2,
    function(x) c(mean  = mean(x),
                  lower = quantile(x, 0.025),
                  upper = quantile(x, 0.975))
  ))
  df           <- as.data.frame(summary_mat)
  colnames(df) <- c("mean", "lower", "upper")
  df$bin       <- as.numeric(sub(paste0("^", param_prefix, "\\."), "", param_cols))
  df$effect_type <- ifelse(
    df$lower > 0 & df$upper > 0, "positive significant",
    ifelse(df$lower < 0 & df$upper < 0, "negative significant",
           ifelse(df$mean > 0, "positive non-significant", "negative non-significant"))
  )
  df$path <- path_label
  return(df)
}


get_indirect_effect_by_bin <- function(post, delta_prefix, beta_prefix, path_label) {
  delta_cols  <- sort(grep(paste0("^", delta_prefix, "\\."), names(post), value = TRUE))
  beta_cols   <- sort(grep(paste0("^", beta_prefix,  "\\."), names(post), value = TRUE))
  sorted_bins <- as.numeric(sub(paste0("^", delta_prefix, "\\."), "", delta_cols))
  indirect    <- post[, delta_cols] * post[, beta_cols]
  summary_mat <- t(apply(indirect, 2,
                         function(x) c(mean  = mean(x),
                                       lower = quantile(x, 0.025),
                                       upper = quantile(x, 0.975))
  ))
  df           <- as.data.frame(summary_mat)
  colnames(df) <- c("mean", "lower", "upper")
  df$bin       <- sorted_bins
  df$effect_type <- ifelse(
    df$lower > 0 & df$upper > 0, "positive significant",
    ifelse(df$lower < 0 & df$upper < 0, "negative significant",
           ifelse(df$mean > 0, "positive non-significant", "negative non-significant"))
  )
  df$effect_type <- factor(df$effect_type, levels = c(
    "negative significant", "negative non-significant",
    "positive non-significant", "positive significant"))
  df$path <- path_label
  return(df)
}


all_effects <- bind_rows(
  get_effect_by_bin(post, "beta3",  "FD → AGB"),
  get_effect_by_bin(post, "beta4",  "CWM → AGB"),
  get_effect_by_bin(post, "beta5",  "DBHCV → AGB"),
  get_effect_by_bin(post, "delta1", "FD → DBHCV"),
  get_effect_by_bin(post, "delta2", "CWM → DBHCV"),
  get_indirect_effect_by_bin(post, "delta1", "beta5", "FD → DBHCV → AGB"),
  get_indirect_effect_by_bin(post, "delta2", "beta5", "CWM → DBHCV → AGB")
)


data_sem2_org     <- data_sem2
data_sem2_org$bin <- as.numeric(factor(data_sem2_org$PDSI_ATA_bin))

result <- data_sem2_org %>%
  group_by(bin) %>%
  dplyr::summarise(
    mean_MAT    = mean(mean_annual_temp,   na.rm = TRUE),
    mean_PRE    = mean(mean_annual_precip, na.rm = TRUE),
    mean_CWD    = mean(CWD1,           na.rm = TRUE),  
    mean_STDAGE = mean(STDAGE,             na.rm = TRUE),
    .groups = "drop"
  )

all_effects1 <- merge(all_effects, result, by = "bin")

write.csv(all_effects1,
          "/data/posterior_American_CWD9_DBH5_Plant_height.csv",
          row.names = FALSE)



##### Proportional stacked bar chart #################################################################


library(scales)

get_effect_prop <- function(post, param_prefix, path_label) {
  param_cols <- grep(paste0("^", param_prefix, "\\."), names(post), value = TRUE)
  df <- t(apply(post[, param_cols, drop = FALSE], 2,
                function(x) c(mean  = mean(x),
                              lower = quantile(x, 0.025),
                              upper = quantile(x, 0.975))
  )) |> as.data.frame()
  colnames(df) <- c("mean", "lower", "upper")
  df$effect_type <- ifelse(
    df$lower > 0 & df$upper > 0, "positive significant",
    ifelse(df$lower < 0 & df$upper < 0, "negative significant",
           ifelse(df$mean > 0, "positive non-significant", "negative non-significant"))
  )
  dplyr::count(df, effect_type) %>%
    mutate(percent = n / sum(n) * 100, path = path_label)
}

get_indirect_prop <- function(post, delta_prefix, beta_prefix, path_label) {
  delta_cols <- sort(grep(paste0("^", delta_prefix, "\\."), names(post), value = TRUE))
  beta_cols  <- sort(grep(paste0("^", beta_prefix,  "\\."), names(post), value = TRUE))
  indirect   <- post[, delta_cols] * post[, beta_cols]
  df <- t(apply(indirect, 2,
                function(x) c(mean  = mean(x),
                              lower = quantile(x, 0.025),
                              upper = quantile(x, 0.975))
  )) |> as.data.frame()
  colnames(df) <- c("mean", "lower", "upper")
  df$effect_type <- ifelse(
    df$lower > 0 & df$upper > 0, "positive significant",
    ifelse(df$lower < 0 & df$upper < 0, "negative significant",
           ifelse(df$mean > 0, "positive non-significant", "negative non-significant"))
  )
  dplyr::count(df, effect_type) %>%
    mutate(percent = n / sum(n) * 100, path = path_label)
}

prop <- bind_rows(
  get_effect_prop(post, "beta3",  "FD → AGB"),
  get_effect_prop(post, "beta4",  "CWM → AGB"),
  get_effect_prop(post, "beta5",  "DBHCV → AGB"),
  get_effect_prop(post, "delta1", "FD → DBHCV"),
  get_effect_prop(post, "delta2", "CWM → DBHCV"),
  get_indirect_prop(post, "delta1", "beta5", "FD → DBHCV → AGB"),
  get_indirect_prop(post, "delta2", "beta5", "CWM → DBHCV → AGB")
)

write.csv(prop,
          "/data/bar_American_CWD9_DBH5_Plant_height.csv",
          row.names = FALSE)


prop_plot <- prop %>%
  filter(path %in% c("FD → AGB", "FD → DBHCV → AGB",
                     "CWM → AGB", "CWM → DBHCV → AGB")) %>%
  mutate(
    path = factor(path, levels = c(
      "FD → AGB", "FD → DBHCV → AGB",
      "CWM → AGB", "CWM → DBHCV → AGB"
    )),
    effect_type = factor(effect_type, levels = c(
      "negative significant", "negative non-significant",
      "positive non-significant", "positive significant"
    ))
  )

fig_height <- ggplot(prop_plot,
                     aes(x = percent, y = path, fill = effect_type)) +
  geom_bar(stat = "identity", width = 0.65) +
  scale_x_continuous(labels = percent_format(scale = 1), expand = c(0, 0)) +
  scale_fill_manual(
    values = c(
      "negative significant"     = "#92C5DE",
      "negative non-significant" = "#dbebfa",
      "positive non-significant" = "#f9ebdf",
      "positive significant"     = "#F4A582"
    ),
    drop = FALSE
  ) +
  labs(x = "Proportion of bins (%)", y = NULL, fill = "Effect type") +
  ggtitle("Structural diversity → Carbon") +
  theme_classic(base_size = 13) +
  theme(axis.text.y  = element_text(size = 12),
        legend.position = "right")

fig_height
##===============================================================================
#####Shade.tolerance###################################################################
##===============================================================================
##### Load data ###################################################################

library(data.table)
library(dplyr)
library(R2jags)

C_datatotal_all <- fread("/data/Data.csv" )
colnames( C_datatotal_all  )

C_datatotal_all <- C_datatotal_all %>%
  dplyr::select(-V1) %>%
  filter(Carbon_Mg_ha > 0)
colnames(C_datatotal_all) <- gsub("_mean$", "", colnames(C_datatotal_all))
C_datatotal_all <- C_datatotal_all %>%
  filter(!is.na(CWD_stand_age)) %>%
  filter(!is.na(MAT_stand_age)) %>%
  filter(!is.na(MAP_stand_age)) 

# 合并 bin 分组
C_bin <- read.csv("/data/Class_CWD_9bins")
C_bin <- C_bin[, c("PLT_CN", "PDSI_ATA_bin")]
C_datatotal_all$PLT_CN <- as.character(C_datatotal_all$PLT_CN)
C_bin$PLT_CN            <- as.character(C_bin$PLT_CN)
C_datatotal_all <- merge(C_datatotal_all, C_bin, by = "PLT_CN")
C_datatotal_all$CWD1=C_datatotal_all$CWD_stand_age
# 标准化
C_datatotal_all <- C_datatotal_all %>%
  mutate(
    Carbon_Mg_ha = log(Carbon_Mg_ha+1),
    shannon_dbhR = scale(H_D_5),
           
           
    CWD          = scale(CWD_stand_age),  
    TEM1         = scale(MAT_stand_age),
    PRE1         = scale(MAP_stand_age),
    FDis         = scale(Shade.tolerance.FDis),
    CWM          = scale(Shade.tolerance.CWM),
    pH           = scale(nitrogen),
    bdod         = scale(bdod),
    STDAGE       = scale(STDAGE)
  )

data_sem2 <- C_datatotal_all

##### JAGS #########################################################


data_all <- data_sem2

y_CV      <- as.numeric(data_all$shannon_dbhR)
y_AGB     <- as.numeric(data_all$Carbon_Mg_ha)
trait     <- as.numeric(data_all$FDis)
trait_CWM <- as.numeric(data_all$CWM)
CWD       <- as.numeric(data_all$CWD)       
TEM1      <- as.numeric(data_all$TEM1)
PRE1      <- as.numeric(data_all$PRE1)
pH        <- as.numeric(data_all$pH)
bdod      <- as.numeric(data_all$bdod)
STDAGE    <- as.numeric(data_all$STDAGE)

bin   <- as.numeric(factor(data_all$PDSI_ATA_bin))
N     <- nrow(data_all)
N_bin <- max(bin)

cat("N =", N, "| N_bin =", N_bin, "\n")


FD_bin     <- tapply(trait,     bin, mean, na.rm = TRUE)
FI_bin     <- tapply(trait_CWM, bin, mean, na.rm = TRUE)
CWD_bin    <- tapply(CWD,       bin, mean, na.rm = TRUE)  
TEM1_bin   <- tapply(TEM1,      bin, mean, na.rm = TRUE)
PRE1_bin   <- tapply(PRE1,      bin, mean, na.rm = TRUE)
pH_bin     <- tapply(pH,        bin, mean, na.rm = TRUE)
bdod_bin   <- tapply(bdod,      bin, mean, na.rm = TRUE)
STDAGE_bin <- tapply(STDAGE,    bin, mean, na.rm = TRUE)



print(sapply(
  list(FD=FD_bin, FI=FI_bin, CWD=CWD_bin,
       TEM1=TEM1_bin, PRE1=PRE1_bin,
       pH=pH_bin, bdod=bdod_bin, STDAGE=STDAGE_bin),
  function(x) sum(is.na(x))
))


jags_data <- list(
  N            = N,
  N_bin        = N_bin,
  DBH_CV       = y_CV,
  Carbon_Mg_ha = y_AGB,
  trait        = trait,
  trait_CWM    = trait_CWM,
  CWD          = CWD,                        
  TEM1         = TEM1,
  PRE1         = PRE1,
  pH           = pH,
  bdod         = bdod,
  STDAGE       = STDAGE,
  bin          = bin,
  FD_bin       = as.numeric(FD_bin),
  FI_bin       = as.numeric(FI_bin),
  CWD_bin      = as.numeric(CWD_bin),        
  TEM1_bin     = as.numeric(TEM1_bin),
  PRE1_bin     = as.numeric(PRE1_bin),
  pH_bin       = as.numeric(pH_bin),
  bdod_bin     = as.numeric(bdod_bin),
  STDAGE_bin   = as.numeric(STDAGE_bin)
)



##### parameters##########################################


params <- c(
  "delta1", "delta2",
  "beta3",  "beta4",  "beta5",
  "gamma0_CV",
  "mu_alpha_AGB",
  "gamma_CV",
  "gamma_CV_bin",
  "gamma_AGB",
  "beta5_0",
  "beta5_1", "beta5_2",
  "beta5_3", "beta5_4",
  "beta5_5", "beta5_6",
  "beta5_7", "beta5_8",
  "sigma_CV",       "sigma_AGB",
  "sigma_alpha_CV", "sigma_alpha_AGB",
  "sigma_delta1",   "sigma_delta2",
  "sigma_beta3",    "sigma_beta4",    "sigma_beta5"
)



##### Run JAGS #################################################################


fit_sem_jags <- jags(
  model.file         = "/data/bayes_sem_model.txt",
  data               = jags_data,
  parameters.to.save = params,
  n.chains           = 3,
  n.iter             = 10000,
  n.burnin           = 2000,
  n.thin             = 10
)



##### Save #############################################################


post        <- as.data.frame(fit_sem_jags$BUGSoutput$sims.list)
summary_tab <- as.data.frame(fit_sem_jags$BUGSoutput$summary)

# 收敛诊断
rhats     <- fit_sem_jags$BUGSoutput$summary[, "Rhat"]
bad_rhats <- rhats[rhats > 1.1]
if (length(bad_rhats) > 0) {
  cat("（Rhat > 1.1）:\n"); print(bad_rhats)
} else {
  cat("Rhat < 1.1）\n")
}


saveRDS(fit_sem_jags,
        "/data/beiyesi_CWD9_DBH5_Shade.tolerance.rds")
write.csv(summary_tab,
          "/data/Beiyrsi_Summary_American_CWD9_DBH5_Shade.tolerance.csv",
          row.names = TRUE)
saveRDS(post,
        "/data/post_American_CWD9_DBH5_Shade.tolerance.rds")



##### Extract #####################################################



get_effect_by_bin <- function(post, param_prefix, path_label) {
  param_cols  <- grep(paste0("^", param_prefix, "\\."), names(post), value = TRUE)
  summary_mat <- t(apply(
    post[, param_cols, drop = FALSE], 2,
    function(x) c(mean  = mean(x),
                  lower = quantile(x, 0.025),
                  upper = quantile(x, 0.975))
  ))
  df           <- as.data.frame(summary_mat)
  colnames(df) <- c("mean", "lower", "upper")
  df$bin       <- as.numeric(sub(paste0("^", param_prefix, "\\."), "", param_cols))
  df$effect_type <- ifelse(
    df$lower > 0 & df$upper > 0, "positive significant",
    ifelse(df$lower < 0 & df$upper < 0, "negative significant",
           ifelse(df$mean > 0, "positive non-significant", "negative non-significant"))
  )
  df$path <- path_label
  return(df)
}


get_indirect_effect_by_bin <- function(post, delta_prefix, beta_prefix, path_label) {
  delta_cols  <- sort(grep(paste0("^", delta_prefix, "\\."), names(post), value = TRUE))
  beta_cols   <- sort(grep(paste0("^", beta_prefix,  "\\."), names(post), value = TRUE))
  sorted_bins <- as.numeric(sub(paste0("^", delta_prefix, "\\."), "", delta_cols))
  indirect    <- post[, delta_cols] * post[, beta_cols]
  summary_mat <- t(apply(indirect, 2,
                         function(x) c(mean  = mean(x),
                                       lower = quantile(x, 0.025),
                                       upper = quantile(x, 0.975))
  ))
  df           <- as.data.frame(summary_mat)
  colnames(df) <- c("mean", "lower", "upper")
  df$bin       <- sorted_bins
  df$effect_type <- ifelse(
    df$lower > 0 & df$upper > 0, "positive significant",
    ifelse(df$lower < 0 & df$upper < 0, "negative significant",
           ifelse(df$mean > 0, "positive non-significant", "negative non-significant"))
  )
  df$effect_type <- factor(df$effect_type, levels = c(
    "negative significant", "negative non-significant",
    "positive non-significant", "positive significant"))
  df$path <- path_label
  return(df)
}


all_effects <- bind_rows(
  get_effect_by_bin(post, "beta3",  "FD → AGB"),
  get_effect_by_bin(post, "beta4",  "CWM → AGB"),
  get_effect_by_bin(post, "beta5",  "DBHCV → AGB"),
  get_effect_by_bin(post, "delta1", "FD → DBHCV"),
  get_effect_by_bin(post, "delta2", "CWM → DBHCV"),
  get_indirect_effect_by_bin(post, "delta1", "beta5", "FD → DBHCV → AGB"),
  get_indirect_effect_by_bin(post, "delta2", "beta5", "CWM → DBHCV → AGB")
)


data_sem2_org     <- data_sem2
data_sem2_org$bin <- as.numeric(factor(data_sem2_org$PDSI_ATA_bin))

result <- data_sem2_org %>%
  group_by(bin) %>%
  dplyr::summarise(
    mean_MAT    = mean(mean_annual_temp,   na.rm = TRUE),
    mean_PRE    = mean(mean_annual_precip, na.rm = TRUE),
    mean_CWD    = mean(CWD1,           na.rm = TRUE),  
    mean_STDAGE = mean(STDAGE,             na.rm = TRUE),
    .groups = "drop"
  )

all_effects1 <- merge(all_effects, result, by = "bin")

write.csv(all_effects1,
          "/data/posterior_American_CWD9_DBH5_Shade.tolerance.csv",
          row.names = FALSE)



##### Proportional stacked bar chart #################################################################


library(scales)

get_effect_prop <- function(post, param_prefix, path_label) {
  param_cols <- grep(paste0("^", param_prefix, "\\."), names(post), value = TRUE)
  df <- t(apply(post[, param_cols, drop = FALSE], 2,
                function(x) c(mean  = mean(x),
                              lower = quantile(x, 0.025),
                              upper = quantile(x, 0.975))
  )) |> as.data.frame()
  colnames(df) <- c("mean", "lower", "upper")
  df$effect_type <- ifelse(
    df$lower > 0 & df$upper > 0, "positive significant",
    ifelse(df$lower < 0 & df$upper < 0, "negative significant",
           ifelse(df$mean > 0, "positive non-significant", "negative non-significant"))
  )
  dplyr::count(df, effect_type) %>%
    mutate(percent = n / sum(n) * 100, path = path_label)
}

get_indirect_prop <- function(post, delta_prefix, beta_prefix, path_label) {
  delta_cols <- sort(grep(paste0("^", delta_prefix, "\\."), names(post), value = TRUE))
  beta_cols  <- sort(grep(paste0("^", beta_prefix,  "\\."), names(post), value = TRUE))
  indirect   <- post[, delta_cols] * post[, beta_cols]
  df <- t(apply(indirect, 2,
                function(x) c(mean  = mean(x),
                              lower = quantile(x, 0.025),
                              upper = quantile(x, 0.975))
  )) |> as.data.frame()
  colnames(df) <- c("mean", "lower", "upper")
  df$effect_type <- ifelse(
    df$lower > 0 & df$upper > 0, "positive significant",
    ifelse(df$lower < 0 & df$upper < 0, "negative significant",
           ifelse(df$mean > 0, "positive non-significant", "negative non-significant"))
  )
  dplyr::count(df, effect_type) %>%
    mutate(percent = n / sum(n) * 100, path = path_label)
}

prop <- bind_rows(
  get_effect_prop(post, "beta3",  "FD → AGB"),
  get_effect_prop(post, "beta4",  "CWM → AGB"),
  get_effect_prop(post, "beta5",  "DBHCV → AGB"),
  get_effect_prop(post, "delta1", "FD → DBHCV"),
  get_effect_prop(post, "delta2", "CWM → DBHCV"),
  get_indirect_prop(post, "delta1", "beta5", "FD → DBHCV → AGB"),
  get_indirect_prop(post, "delta2", "beta5", "CWM → DBHCV → AGB")
)

write.csv(prop,
          "/data/bar_American_CWD9_DBH5_Shade.tolerance.csv",
          row.names = FALSE)


prop_plot <- prop %>%
  filter(path %in% c("FD → AGB", "FD → DBHCV → AGB",
                     "CWM → AGB", "CWM → DBHCV → AGB")) %>%
  mutate(
    path = factor(path, levels = c(
      "FD → AGB", "FD → DBHCV → AGB",
      "CWM → AGB", "CWM → DBHCV → AGB"
    )),
    effect_type = factor(effect_type, levels = c(
      "negative significant", "negative non-significant",
      "positive non-significant", "positive significant"
    ))
  )

fig_height <- ggplot(prop_plot,
                     aes(x = percent, y = path, fill = effect_type)) +
  geom_bar(stat = "identity", width = 0.65) +
  scale_x_continuous(labels = percent_format(scale = 1), expand = c(0, 0)) +
  scale_fill_manual(
    values = c(
      "negative significant"     = "#92C5DE",
      "negative non-significant" = "#dbebfa",
      "positive non-significant" = "#f9ebdf",
      "positive significant"     = "#F4A582"
    ),
    drop = FALSE
  ) +
  labs(x = "Proportion of bins (%)", y = NULL, fill = "Effect type") +
  ggtitle("Structural diversity → Carbon") +
  theme_classic(base_size = 13) +
  theme(axis.text.y  = element_text(size = 12),
        legend.position = "right")

fig_height
##===============================================================================
## Figure  ---------
##===============================================================================

library(dplyr)
library(tidyverse)
library(mgcv)

traits <- c("FD", "Nmass", "LDMC", "LMA",
            "WD", "Plant_height",
            "Shade.tolerance", "Drought.tolerance",
            "Leaf.longevity",
            
            "Pmass",
            "SLA")

base_dir <- "/data/"



get_slope_and_pred <- function(all_effects1, post, param_prefix, clim_col = "mean_CWD",
                               path_label, current_trait = "") {
  
  path_data <- all_effects1 %>%
    filter(path == path_label) %>%
    select(bin, !!sym(clim_col), mean, lower, upper) %>%
    arrange(bin) %>%
    filter(!is.na(!!sym(clim_col)))
  
  clim_vec <- path_data[[clim_col]]
  effect_vec <- path_data$mean
  clim_vec2 <- clim_vec^2
  sd_clim <- sd(clim_vec)
  
  
  fit_lin_mean <- lm(effect_vec ~ clim_vec)
  fit_quad_mean <- lm(effect_vec ~ clim_vec + clim_vec2)
  
  
  linear_exceptions <- list(
    list(path = "CWM → AGB", traits = c("FD", "WD", "LMA", "LDMC", "Leaf.longevity", 
                                        "Nmass", "Pmass", "SLA","Plant_height",  
                                        "Shade.tolerance", "Drought.tolerance")),
    list(path = "FD → AGB", traits = c("FD", "WD", "LMA", "LDMC", "Leaf.longevity", 
                                       "Nmass", "Pmass", "SLA","Plant_height",  
                                       "Shade.tolerance", "Drought.tolerance")),
    
    list(path = "FD → DBHCV → AGB", traits = c("FD","WD","Plant_height",  "Shade.tolerance",  "LMA", "SLA","LDMC", "Nmass",
                                               "Pmass", "Leaf.longevity","Drought.tolerance")),
    
    list(path = "CWM → DBHCV → AGB", traits = c("FD","Nmass", "Pmass", "WD", "Plant_height" ,
                                                "Drought.tolerance"))
  )
  
  is_linear <- FALSE
  for (exc in linear_exceptions) {
    if (path_label == exc$path && current_trait %in% exc$traits) {
      is_linear <- TRUE
      break
    }
  }
  
  is_nonlinear <- !is_linear
  
  
  param_cols <- sort(grep(paste0("^", param_prefix, "\\."), names(post), value = TRUE))
  
  if (length(param_cols) > 0) {
    posterior_mat <- as.matrix(post[, param_cols])
    n_samp <- nrow(posterior_mat)
    
    lin_slopes <- numeric(n_samp)
    quad_coefs <- numeric(n_samp)
    
    for (i in 1:n_samp) {
      coefs <- posterior_mat[i, ]
      
      if (is_nonlinear) {
        fit_i <- lm(coefs ~ clim_vec + clim_vec2)
        lin_slopes[i] <- coef(fit_i)[2] * sd_clim
        quad_coefs[i] <- coef(fit_i)[3] * sd_clim^2
      } else {
        fit_i <- lm(coefs ~ clim_vec)
        lin_slopes[i] <- coef(fit_i)[2] * sd_clim
        quad_coefs[i] <- NA
      }
    }
    
    lin_slopes_valid <- lin_slopes[!is.na(lin_slopes)]
    quad_coefs_valid <- quad_coefs[!is.na(quad_coefs)]
    
    lin_mean_s <- median(lin_slopes_valid, na.rm = TRUE)
    lin_lo_s <- quantile(lin_slopes_valid, 0.025, na.rm = TRUE)
    lin_hi_s <- quantile(lin_slopes_valid, 0.975, na.rm = TRUE)
    lin_sig <- (lin_lo_s > 0 | lin_hi_s < 0)
    
    if (is_nonlinear && length(quad_coefs_valid) > 0) {
      quad_mean_s <- median(quad_coefs_valid, na.rm = TRUE)
      quad_lo_s <- quantile(quad_coefs_valid, 0.025, na.rm = TRUE)
      quad_hi_s <- quantile(quad_coefs_valid, 0.975, na.rm = TRUE)
      quad_sig <- (quad_lo_s > 0 | quad_hi_s < 0)
      
      slope_sig <- lin_sig | quad_sig
      slope_lab <- sprintf(
        "Lin=%.2f(%.2f,%.2f)%s\nQuad=%.2f(%.2f,%.2f)%s",
        lin_mean_s, lin_lo_s, lin_hi_s, ifelse(lin_sig, "*", ""),
        quad_mean_s, quad_lo_s, quad_hi_s, ifelse(quad_sig, "*", "")
      )
    } else {
      slope_sig <- lin_sig
      slope_lab <- sprintf("Slope=%.2f(%.2f,%.2f)%s",
                           lin_mean_s, lin_lo_s, lin_hi_s,
                           ifelse(lin_sig, "*", ""))
    }
    
  } else {
    lin_coef <- coef(fit_lin_mean)[2]
    lin_se <- summary(fit_lin_mean)$coefficients[2, 2]
    lin_ci_lo <- lin_coef - 1.96 * lin_se
    lin_ci_hi <- lin_coef + 1.96 * lin_se
    lin_sig <- (lin_ci_lo > 0 | lin_ci_hi < 0)
    
    lin_mean_s <- lin_coef * sd_clim
    lin_lo_s   <- lin_ci_lo * sd_clim
    lin_hi_s   <- lin_ci_hi * sd_clim
    
    if (is_nonlinear) {
      quad_coef <- coef(fit_quad_mean)[3]
      quad_se <- summary(fit_quad_mean)$coefficients[3, 2]
      quad_ci_lo <- quad_coef - 1.96 * quad_se
      quad_ci_hi <- quad_coef + 1.96 * quad_se
      quad_sig <- (quad_ci_lo > 0 | quad_ci_hi < 0)
      
      slope_sig <- lin_sig | quad_sig
      quad_mean_s <- quad_coef * sd_clim^2
      quad_lo_s <- quad_ci_lo * sd_clim^2
      quad_hi_s <- quad_ci_hi * sd_clim^2
      
      slope_lab <- sprintf(
        "Lin=%.2f(%.2f,%.2f)%s\nQuad=%.2f(%.2f,%.2f)%s",
        lin_mean_s, lin_lo_s, lin_hi_s, ifelse(lin_sig, "*", ""),
        quad_mean_s, quad_lo_s, quad_hi_s, ifelse(quad_sig, "*", "")
      )
    } else {
      slope_sig <- lin_sig
      slope_lab <- sprintf("Slope=%.2f(%.2f,%.2f)%s",
                           lin_mean_s, lin_lo_s, lin_hi_s,
                           ifelse(lin_sig, "*", ""))
    }
  }
  
  
  clim_grid <- seq(min(clim_vec), max(clim_vec), length.out = 50)
  
  if (is_nonlinear) {
    
    fit_gam <- gam(effect_vec ~ s(clim_vec, k = 3), method = "REML")
    pred_data <- data.frame(clim_vec = clim_grid)
    pred_mean <- as.numeric(predict(fit_gam, newdata = pred_data))  # ★ 修复：as.numeric()
  } else {
    
    pred_mean <- coef(fit_lin_mean)[1] + coef(fit_lin_mean)[2] * clim_grid
  }
  
  
  lower_vals <- path_data$lower
  upper_vals <- path_data$upper
  mean_width <- mean((upper_vals - lower_vals) / 2, na.rm = TRUE)
  
  pred_lo <- pred_mean - mean_width
  pred_hi <- pred_mean + mean_width
  
  data.frame(
    clim_grid = clim_grid,
    pred_mean = pred_mean,
    pred_lo   = pred_lo,
    pred_hi   = pred_hi,
    path      = path_label,
    slope_lab = slope_lab,
    slope_sig = slope_sig,
    use_quad  = is_nonlinear,
    row.names = NULL  
  )
}


get_indirect_slope_and_pred <- function(all_effects1, post, delta_prefix, beta_prefix, 
                                        clim_col = "mean_CWD", path_label, current_trait = "") {
  
  path_data <- all_effects1 %>%
    filter(path == path_label) %>%
    select(bin, !!sym(clim_col), mean, lower, upper) %>%
    arrange(bin) %>%
    filter(!is.na(!!sym(clim_col)))
  
  clim_vec <- path_data[[clim_col]]
  effect_vec <- path_data$mean
  clim_vec2 <- clim_vec^2
  sd_clim <- sd(clim_vec)
  
  
  fit_lin_mean <- lm(effect_vec ~ clim_vec)
  fit_quad_mean <- lm(effect_vec ~ clim_vec + clim_vec2)
  
  
  linear_exceptions <- list(
    list(path = "CWM → AGB", traits = c("FD", "WD", "LMA", "LDMC", "Leaf.longevity", 
                                        "Nmass", "Pmass", "SLA","Plant_height",  
                                        "Shade.tolerance", "Drought.tolerance")),
    list(path = "FD → AGB", traits = c("FD", "WD", "LMA", "LDMC", "Leaf.longevity", 
                                       "Nmass", "Pmass", "SLA","Plant_height",  
                                       "Shade.tolerance", "Drought.tolerance")),
    list(path = "FD → DBHCV → AGB", traits = c("FD","WD","Plant_height",  "Shade.tolerance",  "LMA", "SLA","LDMC", "Nmass",
                                               "Pmass", "Leaf.longevity","Drought.tolerance")),
    
    list(path = "CWM → DBHCV → AGB", traits = c("FD","Nmass", "Pmass", "WD", "Plant_height" ,
                                                "Drought.tolerance"))
  )
  
  is_linear <- FALSE
  for (exc in linear_exceptions) {
    if (path_label == exc$path && current_trait %in% exc$traits) {
      is_linear <- TRUE
      break
    }
  }
  
  is_nonlinear <- !is_linear
  
  
  delta_cols <- sort(grep(paste0("^", delta_prefix, "\\."), names(post), value = TRUE))
  beta_cols  <- sort(grep(paste0("^", beta_prefix,  "\\."), names(post), value = TRUE))
  
  if (length(delta_cols) > 0 && length(beta_cols) > 0) {
    delta_mat <- as.matrix(post[, delta_cols])
    beta_mat  <- as.matrix(post[, beta_cols])
    indirect_mat <- delta_mat * beta_mat
    
    n_samp <- nrow(indirect_mat)
    
    lin_slopes <- numeric(n_samp)
    quad_coefs <- numeric(n_samp)
    
    for (i in 1:n_samp) {
      coefs <- indirect_mat[i, ]
      
      if (is_nonlinear) {
        fit_i <- lm(coefs ~ clim_vec + clim_vec2)
        lin_slopes[i] <- coef(fit_i)[2] * sd_clim
        quad_coefs[i] <- coef(fit_i)[3] * sd_clim^2
      } else {
        fit_i <- lm(coefs ~ clim_vec)
        lin_slopes[i] <- coef(fit_i)[2] * sd_clim
        quad_coefs[i] <- NA
      }
    }
    
    lin_slopes_valid <- lin_slopes[!is.na(lin_slopes)]
    quad_coefs_valid <- quad_coefs[!is.na(quad_coefs)]
    
    lin_mean_s <- median(lin_slopes_valid, na.rm = TRUE)
    lin_lo_s <- quantile(lin_slopes_valid, 0.025, na.rm = TRUE)
    lin_hi_s <- quantile(lin_slopes_valid, 0.975, na.rm = TRUE)
    lin_sig <- (lin_lo_s > 0 | lin_hi_s < 0)
    
    if (is_nonlinear && length(quad_coefs_valid) > 0) {
      quad_mean_s <- median(quad_coefs_valid, na.rm = TRUE)
      quad_lo_s <- quantile(quad_coefs_valid, 0.025, na.rm = TRUE)
      quad_hi_s <- quantile(quad_coefs_valid, 0.975, na.rm = TRUE)
      quad_sig <- (quad_lo_s > 0 | quad_hi_s < 0)
      
      slope_sig <- lin_sig | quad_sig
      slope_lab <- sprintf(
        "Lin=%.2f(%.2f,%.2f)%s\nQuad=%.2f(%.2f,%.2f)%s",
        lin_mean_s, lin_lo_s, lin_hi_s, ifelse(lin_sig, "*", ""),
        quad_mean_s, quad_lo_s, quad_hi_s, ifelse(quad_sig, "*", "")
      )
    } else {
      slope_sig <- lin_sig
      slope_lab <- sprintf("Slope=%.2f(%.2f,%.2f)%s",
                           lin_mean_s, lin_lo_s, lin_hi_s,
                           ifelse(lin_sig, "*", ""))
    }
    
  } else {
    lin_coef <- coef(fit_lin_mean)[2]
    lin_se <- summary(fit_lin_mean)$coefficients[2, 2]
    lin_ci_lo <- lin_coef - 1.96 * lin_se
    lin_ci_hi <- lin_coef + 1.96 * lin_se
    lin_sig <- (lin_ci_lo > 0 | lin_ci_hi < 0)
    
    lin_mean_s <- lin_coef * sd_clim
    lin_lo_s   <- lin_ci_lo * sd_clim
    lin_hi_s   <- lin_ci_hi * sd_clim
    
    if (is_nonlinear) {
      quad_coef <- coef(fit_quad_mean)[3]
      quad_se <- summary(fit_quad_mean)$coefficients[3, 2]
      quad_ci_lo <- quad_coef - 1.96 * quad_se
      quad_ci_hi <- quad_coef + 1.96 * quad_se
      quad_sig <- (quad_ci_lo > 0 | quad_ci_hi < 0)
      
      slope_sig <- lin_sig | quad_sig
      quad_mean_s <- quad_coef * sd_clim^2
      quad_lo_s <- quad_ci_lo * sd_clim^2
      quad_hi_s <- quad_ci_hi * sd_clim^2
      
      slope_lab <- sprintf(
        "Lin=%.2f(%.2f,%.2f)%s\nQuad=%.2f(%.2f,%.2f)%s",
        lin_mean_s, lin_lo_s, lin_hi_s, ifelse(lin_sig, "*", ""),
        quad_mean_s, quad_lo_s, quad_hi_s, ifelse(quad_sig, "*", "")
      )
    } else {
      slope_sig <- lin_sig
      slope_lab <- sprintf("Slope=%.2f(%.2f,%.2f)%s",
                           lin_mean_s, lin_lo_s, lin_hi_s,
                           ifelse(lin_sig, "*", ""))
    }
  }
  
  
  clim_grid <- seq(min(clim_vec), max(clim_vec), length.out = 50)
  
  if (is_nonlinear) {
    
    fit_gam <- gam(effect_vec ~ s(clim_vec, k = 3), method = "REML")
    pred_data <- data.frame(clim_vec = clim_grid)
    pred_mean <- as.numeric(predict(fit_gam, newdata = pred_data)) 
  } else {
    
    pred_mean <- coef(fit_lin_mean)[1] + coef(fit_lin_mean)[2] * clim_grid
  }
  
  lower_vals <- path_data$lower
  upper_vals <- path_data$upper
  mean_width <- mean((upper_vals - lower_vals) / 2, na.rm = TRUE)
  
  pred_lo <- pred_mean - mean_width
  pred_hi <- pred_mean + mean_width
  
  data.frame(
    clim_grid = clim_grid,
    pred_mean = pred_mean,
    pred_lo   = pred_lo,
    pred_hi   = pred_hi,
    path      = path_label,
    slope_lab = slope_lab,
    slope_sig = slope_sig,
    use_quad  = is_nonlinear,
    row.names = NULL  
  )
}


path_specs <- list(
  list(prefix = "beta3",  label = "FD → AGB"),
  list(prefix = "beta4",  label = "CWM → AGB"),
  list(prefix = "beta5",  label = "DBHCV → AGB"),
  list(prefix = "delta1", label = "FD → DBHCV"),
  list(prefix = "delta2", label = "CWM → DBHCV")
)

indirect_specs <- list(
  list(delta = "delta1", beta = "beta5", label = "FD → DBHCV → AGB"),
  list(delta = "delta2", beta = "beta5", label = "CWM → DBHCV → AGB")
)


all_pred <- list()
all_eff  <- list()

for (trait in traits) {
  
  rds_file <- file.path(base_dir, paste0("post_American_CWD9_DBH5_", trait, ".rds"))
  csv_file <- file.path(base_dir, paste0("posterior_American_CWD9_DBH5_", trait, ".csv"))
  
  if (!file.exists(rds_file) || !file.exists(csv_file)) {
    message("Skipping trait ", trait, ": file not found")
    next
  }
  
  post <- readRDS(rds_file)
  all_effects1 <- read.csv(csv_file)
  
  all_effects1 <- all_effects1 %>%
    filter(!is.na(mean_CWD)) %>%
    filter(!is.na(path), !is.na(mean), !is.na(lower), !is.na(upper))
  
  message(sprintf(" %s", trait))
  
  tryCatch({
    direct_pred <- bind_rows(lapply(path_specs, function(ps) {
      get_slope_and_pred(
        all_effects1 = all_effects1,
        post         = post,
        param_prefix = ps$prefix,
        clim_col     = "mean_CWD",
        path_label   = ps$label,
        current_trait = trait
      )
    }))
    
    indirect_pred <- bind_rows(lapply(indirect_specs, function(ps) {
      get_indirect_slope_and_pred(
        all_effects1  = all_effects1,
        post          = post,
        delta_prefix  = ps$delta,
        beta_prefix   = ps$beta,
        clim_col      = "mean_CWD",
        path_label    = ps$label,
        current_trait = trait
      )
    }))
    
    pred_trait <- bind_rows(direct_pred, indirect_pred) %>%
      mutate(trait = trait)
    
    all_effects1$trait <- trait
    
    all_pred[[trait]] <- pred_trait
    all_eff[[trait]]  <- all_effects1
    
    message(" ", trait)
    
  }, error = function(e) {
    message(" ", e$message)
  })
}

pred_combined <- bind_rows(all_pred)
eff_combined  <- bind_rows(all_eff)

write.csv(pred_combined, file.path(base_dir, "pred_combined_all_traits_CWD9.csv"), row.names = FALSE)
write.csv(eff_combined,  file.path(base_dir, "eff_combined_all_traits_CWD9.csv"),  row.names = FALSE)


##### Fig.4 #####################################################

library(ggplot2)
library(tidyverse)

pred_combined <- read.csv("E:/American_NFI/data_beiyesi/pred_combined_all_traits_CWD9.csv")
eff_combined  <- read.csv("E:/American_NFI/data_beiyesi/eff_combined_all_traits_CWD9.csv")


target_paths <- c(
  "FD → AGB",
  "FD → DBHCV → AGB",
  "CWM → AGB",
  "CWM → DBHCV → AGB"
)

pred_sub <- pred_combined %>% 
  filter(path %in% target_paths) %>%
  filter(trait %in% c("FD", "Drought.tolerance","LMA","LDMC","Nmass","Pmass"))

eff_sub <- eff_combined %>% 
  filter(path %in% target_paths) %>%
  filter(trait %in% c("FD", "Drought.tolerance","LMA","LDMC","Nmass","Pmass"))

# ── effect group ─────────────────────────────────────────
pred_sub <- pred_sub %>%
  mutate(
    effect_group = ifelse(path %in% c("FD → AGB", "CWM → AGB"),
                          "Direct effects",
                          "Indirect effects")
  )

eff_sub <- eff_sub %>%
  mutate(
    effect_group = ifelse(path %in% c("FD → AGB", "CWM → AGB"),
                          "Direct effects",
                          "Indirect effects"),
    sig = case_when(
      lower > 0 & upper > 0 ~ "significant",
      lower < 0 & upper < 0 ~ "significant",
      TRUE ~ "non-significant"
    )
  )

pred_sub$effect_group <- factor(pred_sub$effect_group,
                                levels = c("Direct effects", "Indirect effects"))
eff_sub$effect_group  <- factor(eff_sub$effect_group,
                                levels = c("Direct effects", "Indirect effects"))


path_levels <- c("FD → AGB", "FD → DBHCV → AGB",
                 "CWM → AGB", "CWM → DBHCV → AGB")

path_labels <- c("FD → AGC", "FD → SD → AGC",
                 "CWM → AGC", "CWM → SD → AGC")

pred_sub$path <- factor(pred_sub$path, levels = path_levels, labels = path_labels)
eff_sub$path  <- factor(eff_sub$path,  levels = path_levels, labels = path_labels)

# ── trait labels ─────────────────────────────────────────
trait_levels <- c("FD","Drought.tolerance","LMA","LDMC","Nmass","Pmass")

trait_labels <- c(
  "Multivariate FD\n/CWM PC1",
  "Drought\ntolerance",
  "Leaf mass\nper area",
  "Leaf dry\nmatter content",
  "Leaf nitrogen\nconcentration",
  "Leaf phosphorus\nconcentration"
)

pred_sub$trait <- factor(pred_sub$trait, levels = trait_levels, labels = trait_labels)
eff_sub$trait  <- factor(eff_sub$trait,  levels = trait_levels, labels = trait_labels)

# ── slope labels ─────────────────────────────────────────
slope_labels <- pred_sub %>%
  group_by(trait, path, effect_group) %>%
  slice(1) %>%
  select(trait, path, effect_group, slope_lab, slope_sig) %>%
  ungroup()

# ── facet range ─────────────────────────────
facet_range <- eff_sub %>%
  group_by(trait, effect_group) %>%
  summarise(
    y_max = max(upper, na.rm = TRUE),
    y_min = min(lower, na.rm = TRUE),
    y_range = y_max - y_min,
    .groups = "drop"
  )

col_xmin <- eff_sub %>%
  group_by(trait) %>%
  summarise(x_pos = min(mean_CWD, na.rm = TRUE), .groups = "drop")



slope_labels <- slope_labels %>%
  left_join(facet_range, by = c("trait", "effect_group")) %>%
  left_join(col_xmin, by = "trait") %>%
  mutate(
    y_pos = case_when(
      path == "FD → AGC"      ~ 0.65,
      path ==  "CWM → AGC"~ 0.55,
      
      path == "FD → SD → AGC"      ~ 0.8,
      path == "CWM → SD → AGC" ~ 0.65
      
    )
  ) %>%
  select(trait, path, effect_group, slope_lab, slope_sig, x_pos, y_pos)
# ── group（FD vs CWM）──────────────────────────────────
pred_sub <- pred_sub %>%
  mutate(group = ifelse(grepl("^FD", as.character(path)), "FD", "CWM"))

eff_sub <- eff_sub %>%
  mutate(group = ifelse(grepl("^FD", as.character(path)), "FD", "CWM"))

slope_labels <- slope_labels %>%
  mutate(group = ifelse(grepl("^FD", as.character(path)), "FD", "CWM"))

# ── plot ───────────────────────────────────────────────
fig_grid <- ggplot() +
  
  
  geom_ribbon(
    data = pred_sub,
    aes(
      x = clim_grid,
      ymin = pred_lo,
      ymax = pred_hi,
      fill = group
    ),
    alpha =0.6
  ) +
  scale_fill_manual(
    values = c(
      FD =  "#f5f5f5",
      CWM ="#EDEDE4"
    )
  )+
  scale_fill_manual(
    name = NULL,
    breaks = c("FD", "CWM"),
    labels = c("Functional diversity", "Community-weighted mean"),
    values = c(
      FD  = "#f5f5f5",
      CWM = "#EDEDE4"
    )
  ) +
  
  scale_color_manual(
    name = NULL,
    breaks = c("FD", "CWM"),
    labels = c("Functional diversity", "Community-weighted mean"),
    values = c(
      FD  = "#D37C5E",
      CWM = "#9DC7C6"
    )
  ) +
  
  
  
  
  geom_line(
    data = pred_sub,
    aes(x = clim_grid, y = pred_mean, color = group, group = path),
    linewidth = 0.6
  ) +
  
  geom_errorbar(
    data = eff_sub,
    aes(x = mean_CWD, ymin = lower, ymax = upper, color = group),
    width = 0, linewidth = 0.3
  ) +
  
  geom_point(
    data = eff_sub,
    aes(x = mean_CWD, y = mean, color = group),
    size = 1.5
  ) +
  
  geom_text(
    data = slope_labels,
    aes(x = x_pos, y = y_pos, label = slope_lab, color = group),
    hjust = 0, vjust = 0.5,
    size = 3, family = "serif", fontface = "bold"
  ) +
  
  geom_hline(yintercept = 0, linetype = "dashed",
             color = "#6a737b", linewidth = 0.5) +
  
  facet_grid(effect_group ~ trait, scales = "free_y") +
  
  
  
  scale_y_continuous(
    expand = expansion(mult = c(0.05, 0.30))
  ) +
  
  coord_cartesian(clip = "off") +
  
  geom_hline(yintercept = 0, color = "#6a737b", linetype = "dashed", linewidth = 0.5) +
  scale_y_continuous(
    expand = expansion(mult = c(0.05, 0.1))
  )+
  
  facet_grid(effect_group ~ trait, scales = "free_y") +
  theme_classic(base_size = 10) +
  labs(x = "Climatic water deficit (mm)", y = "Biodiversity effects on AGC") +
  theme(
    legend.position = "bottom",
    
    legend.text = element_text(size = 14, family = "serif"),
    legend.title = element_blank(),
    legend.text.align = 1,
    legend.key.size = unit(0.4, "cm"),
    panel.border = element_rect(
      colour = "#6a737b",
      fill = NA,
      linewidth = 0.5
    ),
    panel.grid = element_blank(),
    #panel.background = element_blank(),
    plot.background = element_rect(fill = "white", color = NA),
    strip.background = element_rect(
      fill = "white",
      colour = NA
    ),
    panel.background = element_rect(fill = "white", colour = NA),
    axis.line = element_blank( ),
    
    axis.ticks.x = element_line(size = 0.5),
    axis.ticks.y =  element_line(size = 0.5),
    
    axis.text.x = element_text(size = 12, family = "serif", color = "black"),
    axis.text.y = element_text(size = 12, family = "serif", color = "black"),
    
    axis.title = element_text(size = 14, family = "serif", color = "black"),
    
    plot.title = element_text(size = 14, family = "serif"),
    plot.subtitle = element_text(hjust = 0, size = 14, family = "serif"),
    
    strip.text = element_text(size = 14, family = "serif")#,
    #strip.background = element_blank()
  )

fig_grid

ggsave(
  "fig4.png",
  path = "/Figure/",
  width = 9.5, height = 6,
  dpi = 600
)

##### Fig.S other trait #####################################################

library(ggplot2)
library(tidyverse)

pred_combined <- read.csv("E:/American_NFI/data_beiyesi/pred_combined_all_traits_CWD9.csv")
eff_combined  <- read.csv("E:/American_NFI/data_beiyesi/eff_combined_all_traits_CWD9.csv")


target_paths <- c(
  "FD → AGB",
  "FD → DBHCV → AGB",
  "CWM → AGB",
  "CWM → DBHCV → AGB"
)

pred_sub <- pred_combined %>% 
  filter(path %in% target_paths) %>%
  filter(trait %in% c("WD", "Leaf.longevity",  "SLA", "Plant_height",  "Shade.tolerance" ))

eff_sub <- eff_combined %>% 
  filter(path %in% target_paths) %>%
  filter(trait %in% c("WD", "Leaf.longevity",  "SLA", "Plant_height",  "Shade.tolerance" ))

# ── effect group ─────────────────────────────────────────
pred_sub <- pred_sub %>%
  mutate(
    effect_group = ifelse(path %in% c("FD → AGB", "CWM → AGB"),
                          "Direct effects",
                          "Indirect effects")
  )

eff_sub <- eff_sub %>%
  mutate(
    effect_group = ifelse(path %in% c("FD → AGB", "CWM → AGB"),
                          "Direct effects",
                          "Indirect effects"),
    sig = case_when(
      lower > 0 & upper > 0 ~ "significant",
      lower < 0 & upper < 0 ~ "significant",
      TRUE ~ "non-significant"
    )
  )

pred_sub$effect_group <- factor(pred_sub$effect_group,
                                levels = c("Direct effects", "Indirect effects"))
eff_sub$effect_group  <- factor(eff_sub$effect_group,
                                levels = c("Direct effects", "Indirect effects"))


path_levels <- c("FD → AGB", "FD → DBHCV → AGB",
                 "CWM → AGB", "CWM → DBHCV → AGB")

path_labels <- c("FD → AGC", "FD → SD → AGC",
                 "CWM → AGC", "CWM → SD → AGC")

pred_sub$path <- factor(pred_sub$path, levels = path_levels, labels = path_labels)
eff_sub$path  <- factor(eff_sub$path,  levels = path_levels, labels = path_labels)

# ── trait labels ─────────────────────────────────────────
trait_levels <- c("WD", "Leaf.longevity",  "SLA", "Plant_height",  "Shade.tolerance" )

trait_labels <- c(
  "Wood density", 
  "Leaf longevity", "Specific leaf area" , 
  "Max height",
  "Shade tolerance"
)

pred_sub$trait <- factor(pred_sub$trait, levels = trait_levels, labels = trait_labels)
eff_sub$trait  <- factor(eff_sub$trait,  levels = trait_levels, labels = trait_labels)

# ── slope labels ─────────────────────────────────────────
slope_labels <- pred_sub %>%
  group_by(trait, path, effect_group) %>%
  slice(1) %>%
  select(trait, path, effect_group, slope_lab, slope_sig) %>%
  ungroup()


facet_range <- eff_sub %>%
  group_by(trait, effect_group) %>%
  summarise(
    y_max = max(upper, na.rm = TRUE),
    y_min = min(lower, na.rm = TRUE),
    y_range = y_max - y_min,
    .groups = "drop"
  )

col_xmin <- eff_sub %>%
  group_by(trait) %>%
  summarise(x_pos = min(mean_CWD, na.rm = TRUE), .groups = "drop")



slope_labels <- slope_labels %>%
  left_join(facet_range, by = c("trait", "effect_group")) %>%
  left_join(col_xmin, by = "trait") %>%
  mutate(
    y_pos = case_when(
      path == "FD → AGC"      ~ 0.75,
      path ==  "CWM → AGC"~ 0.65,
      
      path == "FD → SD → AGC"      ~ 0.8,
      path == "CWM → SD → AGC" ~ 0.65
      
    )
  ) %>%
  select(trait, path, effect_group, slope_lab, slope_sig, x_pos, y_pos)
# ── group（FD vs CWM）──────────────────────────────────
pred_sub <- pred_sub %>%
  mutate(group = ifelse(grepl("^FD", as.character(path)), "FD", "CWM"))

eff_sub <- eff_sub %>%
  mutate(group = ifelse(grepl("^FD", as.character(path)), "FD", "CWM"))

slope_labels <- slope_labels %>%
  mutate(group = ifelse(grepl("^FD", as.character(path)), "FD", "CWM"))

# ── plot ───────────────────────────────────────────────
fig_grid <- ggplot() +
  
  # geom_ribbon(
  #   data = pred_sub,
  #   aes(x = clim_grid, ymin = pred_lo, ymax = pred_hi, group = path),
  #   fill = "#F2F2F0", alpha = 0.6
  # ) +
  geom_ribbon(
    data = pred_sub,
    aes(
      x = clim_grid,
      ymin = pred_lo,
      ymax = pred_hi,
      fill = group
    ),
    alpha =0.6
  ) +
  scale_fill_manual(
    values = c(
      FD =  "#f5f5f5",
      CWM ="#EDEDE4"
    )
  )+
  scale_fill_manual(
    name = NULL,
    breaks = c("FD", "CWM"),
    labels = c("Functional diversity", "Community-weighted mean"),
    values = c(
      FD  = "#f5f5f5",
      CWM = "#EDEDE4"
    )
  ) +
  
  scale_color_manual(
    name = NULL,
    breaks = c("FD", "CWM"),
    labels = c("Functional diversity", "Community-weighted mean"),
    values = c(
      FD  = "#D37C5E",
      CWM = "#9DC7C6"
    )
  ) +
  
  
  
  
  geom_line(
    data = pred_sub,
    aes(x = clim_grid, y = pred_mean, color = group, group = path),
    linewidth = 0.6
  ) +
  
  geom_errorbar(
    data = eff_sub,
    aes(x = mean_CWD, ymin = lower, ymax = upper, color = group),
    width = 0, linewidth = 0.3
  ) +
  
  geom_point(
    data = eff_sub,
    aes(x = mean_CWD, y = mean, color = group),
    size = 1.5
  ) +
  
  geom_text(
    data = slope_labels,
    aes(x = x_pos, y = y_pos, label = slope_lab, color = group),
    hjust = 0, vjust = 0.5,
    size = 3, family = "serif", fontface = "bold"
  ) +
  
  geom_hline(yintercept = 0, linetype = "dashed",
             color = "#6a737b", linewidth = 0.5) +
  
  facet_grid(effect_group ~ trait, scales = "free_y") +
  
  
  # scale_color_manual(
  #   values = c(
  #     FD = "#D37C5E",
  #     CWM = "#9DC7C6"
  #   ),
  #   name = NULL,
  #   labels = c("Functional diversity", "Community-weighted mean")
  # ) +
  
  scale_y_continuous(
    expand = expansion(mult = c(0.05, 0.30))
  ) +
  
  coord_cartesian(clip = "off") +
  
  
  geom_hline(yintercept = 0, color = "#6a737b", linetype = "dashed", linewidth = 0.5) +
  scale_y_continuous(
    expand = expansion(mult = c(0.05, 0.1))
  )+
  
  facet_grid(effect_group ~ trait, scales = "free_y") +
  theme_classic(base_size = 10) +
  labs(x = "Climatic water deficit (mm)", y = "Biodiversity effects on AGC") +
  theme(
    legend.position = "bottom",
    
    legend.text = element_text(size = 14, family = "serif"),
    legend.title = element_blank(),
    legend.text.align = 1,
    legend.key.size = unit(0.4, "cm"),
    panel.border = element_rect(
      colour = "#6a737b",
      fill = NA,
      linewidth = 0.5
    ),
    panel.grid = element_blank(),
    #panel.background = element_blank(),
    plot.background = element_rect(fill = "white", color = NA),
    strip.background = element_rect(
      fill = "white",
      colour = NA
    ),
    panel.background = element_rect(fill = "white", colour = NA),
    axis.line = element_blank( ),
    
    axis.ticks.x = element_line(size = 0.5),
    axis.ticks.y =  element_line(size = 0.5),
    
    axis.text.x = element_text(size = 12, family = "serif", color = "black"),
    axis.text.y = element_text(size = 12, family = "serif", color = "black"),
    
    axis.title = element_text(size = 14, family = "serif", color = "black"),
    
    plot.title = element_text(size = 14, family = "serif"),
    plot.subtitle = element_text(hjust = 0, size = 14, family = "serif"),
    
    strip.text = element_text(size = 14, family = "serif")#,
    #strip.background = element_blank()
  )

fig_grid

ggsave(
  "figS_CWD9_othertrait.png",
  path = "/Figure/",
  width = 9.5, height = 6,
  dpi = 600
)
