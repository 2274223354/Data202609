##===============================================================================
## path analysis ---------
##===============================================================================
library(car)
library(lmerTest)
library(piecewiseSEM)
library(data.table)
library(pROC)
library(dplyr)
library(plyr)
library(doParallel)
####log and scale####
library(data.table)


C_datatotal_all <-fread("/data/Data.csv")
C_datatotal_all$PLT_CN <- as.character(C_datatotal_all$PLT_CN)


C_datatotal_all <- C_datatotal_all %>%
  filter(Carbon_Mg_ha > 0)

names(C_datatotal_all)[duplicated(names(C_datatotal_all))]


C_datatotal_all %>%
  filter(!is.na(NA_L2NAME)) %>%
  group_by(NA_L2NAME) %>%
  dplyr::summarise(n = n()) %>%
  arrange(desc(n))

range(C_datatotal_all$STDAGE)


A=as.data.frame(C_datatotal_all[,c("NA_L2NAME","STDAGE"  )])

colnames(C_datatotal_all) <- gsub("_mean$", "", colnames(C_datatotal_all))


C_datatotal_all <- C_datatotal_all %>%
  filter(!is.na(CWD_stand_age)) %>%
  filter(!is.na(MAT_stand_age)) %>%
  filter(!is.na(MAP_stand_age)) 



C_datatotal_all <- C_datatotal_all %>%
  mutate(
    Carbon_Mg_ha = log(Carbon_Mg_ha+1),
    shannon_dbhR = scale(H_D_5),
    
    CWD          = scale(CWD_stand_age), 
    TEM        = scale(MAT_stand_age),
    PRE         = scale(MAP_stand_age),
    FD         = scale(FDis),
    CWM          = scale(PC1.CWM),
    nitrogen_mean           = scale(nitrogen),
    bdod         = scale(bdod),
    STDAGE       = scale(STDAGE)
  )

C_datatotal_all <- C_datatotal_all[complete.cases(C_datatotal_all[, 
                                                                  c("shannon_dbhR","FD","CWM","Carbon_Mg_ha","CWD"  ,
                                                                    "nitrogen_mean","bdod","PRE","TEM","STDAGE",
                                                                    "NA_L2NAME", "LON" ,"LAT" ,  "COUNTYCD")]), ]

C_datatotal_all <- C_datatotal_all[, c("shannon_dbhR","FD","CWM","Carbon_Mg_ha","STDAGE",
                                       "CWD"  ,
                                       "nitrogen_mean","bdod","PRE","TEM","NA_L2NAME", "LON" ,"LAT",
                                       "COUNTYCD" )]

# SEM ------------------------------------------------------------------

data_sem2=C_datatotal_all
data_sem2$NA_L2NAME <- as.factor(data_sem2$NA_L2NAME)
data_sem2$shannon_dbhR <- as.numeric(data_sem2$shannon_dbhR)
data_sem2$FD <- as.numeric(data_sem2$FD)
data_sem2$CWM <- as.numeric(data_sem2$CWM)


library(lavaan)

model_formula <- '
  

  shannon_dbhR ~ FD + CWM+STDAGE + CWD + PRE + TEM+bdod+nitrogen_mean
  FD           ~ STDAGE + CWD + PRE + TEM+nitrogen_mean +STDAGE #+ bdod
  CWM          ~ STDAGE + CWD + PRE + TEM+bdod +nitrogen_mean 

  Carbon_Mg_ha ~ shannon_dbhR + FD + CWM+STDAGE + CWD + PRE + TEM+bdod+nitrogen_mean

  CWM ~~ FD
'

sem_mod <- lavaan::sem(model_formula, data = data_sem2)
#summary(sem_mod, fit.measures = TRUE, standardized = TRUE)
summary(sem_mod, fit.measures = TRUE, standardized = TRUE, rsquare = TRUE)
fitMeasures(sem_mod,c("chisq","df","pvalue","gfi","cfi","rmr","srmr","rmsea","AIC"))
standardizedSolution(sem_mod)



build_knn_weights <- function(coords, k = 100) {
  n <- nrow(coords)
  cat("n =", n, "\n")
  knn_result <- FNN::get.knn(coords, k = k)
  rows <- rep(1:n, each = k)
  cols <- as.vector(knn_result$nn.index)
  dists <- pmax(as.vector(knn_result$nn.dist), 1e-10)
  vals  <- 1 / dists
  W <- sparseMatrix(
    i    = c(rows, cols),
    j    = c(cols, rows),
    x    = c(vals, vals),
    dims = c(n, n)
  )
  W <- drop0(W)
  cat(":", nnzero(W), "\n")
  return(W)
}

coords <- cbind(data_sem2$LON, data_sem2$LAT)
W_knn  <- build_knn_weights(coords, k = 100)



lavResidualsY <- function(object,
                          ynames = lavNames(object, "ov.nox"),
                          xnames = lavNames(object, "ov.x")) {
  pred <- lavPredictY(object, ynames = ynames, xnames = xnames) |> as.data.frame()
  d    <- inspect(object, "data") |> as.data.frame()
  r    <- lapply(names(pred), function(x) pred[[x]] - d[[x]])
  res  <- do.call(cbind, r) |> as.data.frame()
  names(res) <- names(pred)
  res
}

moran_i_sparse <- function(x, W) {
  n    <- length(x)
  x    <- x - mean(x)
  row_sums <- Matrix::rowSums(W)
  row_sums[row_sums == 0] <- 1
  W_norm <- Diagonal(x = 1 / row_sums) %*% W
  Wx       <- as.vector(W_norm %*% x)
  observed <- (n / n) * sum(x * Wx) / sum(x^2)
  expected <- -1 / (n - 1)
  S1 <- sum(W_norm^2 + t(W_norm)^2) / 2
  S2 <- sum((Matrix::rowSums(W_norm) + Matrix::colSums(W_norm))^2)
  S0 <- n
  var_i <- (n^2 * S1 - n * S2 + 3 * S0^2) /
    (S0^2 * (n^2 - 1)) - expected^2
  sd_i  <- sqrt(var_i)
  z     <- (observed - expected) / sd_i
  p_val <- 2 * pnorm(abs(z), lower.tail = FALSE)
  list(observed = observed, expected = expected,
       sd = sd_i, p.value = p_val)
}

lavSpatialCorrect_fast <- function(obj, W, alpha = 0.05) {
  cat("...\n")
  resids <- lavResidualsY(obj)
  n      <- nrow(resids)
  cat(" Moran's I...\n")
  morans_i <- lapply(names(resids), function(vname) {
    cat("  →", vname, "\n")
    x  <- resids[[vname]]
    mi <- moran_i_sparse(x, W)
    mi$n.eff <- if (mi$p.value > alpha) {
      n
    } else {
      n * (1 - mi$observed) / (1 + mi$observed)
    }
    mi
  })
  names(morans_i) <- names(resids)
  cat("SE p ...\n")
  v <- diag(vcov(obj))
  params <- lapply(names(morans_i), function(acol) {
    idx   <- c(grep(paste0(acol, "~"), names(v)),
               grep(paste0("=~", acol), names(v)))
    v_idx <- v[idx] * n / morans_i[[acol]]$n.eff
    ret   <- data.frame(
      Parameter = names(v)[idx],
      Estimate  = coef(obj)[idx],
      n.eff     = round(morans_i[[acol]]$n.eff, 1),
      Std.err   = sqrt(v_idx)
    )
    ret[["Z-value"]] <- ret$Estimate / ret$Std.err
    ret[["P(>|z|)"]] <- 2 * pnorm(abs(ret[["Z-value"]]), lower.tail = FALSE)
    ret
  })
  names(params) <- names(morans_i)
  mi_out <- lapply(morans_i, function(m)
    data.frame(observed = m$observed, expected = m$expected,
               sd = m$sd, p.value = m$p.value, n.eff = m$n.eff))
  list(Morans_I = mi_out, parameters = params)
}



cat("\n...\n")
system.time({
  spatial_result <- lavSpatialCorrect_fast(obj = sem_mod, W = W_knn)
})



moran_df <- do.call(rbind, Map(
  function(df, nm) { df$variable <- nm; df },
  spatial_result$Morans_I,
  names(spatial_result$Morans_I)
))
cat("\n── Moran's I  ──────────────────────────────\n")
print(moran_df)


cat("\n── Carbon_Mg_ha  ──────────────────────\n")
print(spatial_result$parameters$Carbon_Mg_ha)

cat("\n── shannon_dbhR  ──────────────────────\n")
print(spatial_result$parameters$shannon_dbhR)

cat("\n── FD  ────────────────────────────────\n")
print(spatial_result$parameters$FD)

cat("\n── CWM  ───────────────────────────────\n")
print(spatial_result$parameters$CWM)


for (nm in names(spatial_result$parameters)) {
  write.csv(
    spatial_result$parameters[[nm]],
    file      = paste0("SpatialCorrected_", nm, ".csv"),
    row.names = FALSE
  )
}
#write.csv(moran_df, "MoransI_summary.csv", row.names = FALSE)

###############################################################################
## Extract path coefficients ##################################################


library(lavaan)
library(igraph)
library(data.table)

get_spatial_corrected_std <- function(sem_mod, spatial_result) {
  
  
  std_sol <- standardizedSolution(sem_mod)
  
  std_reg <- std_sol[std_sol$op == "~", ]
  
  
  neff_df <- do.call(rbind, lapply(names(spatial_result$Morans_I), function(v) {
    data.frame(
      lhs   = v,
      n.eff = spatial_result$Morans_I[[v]]$n.eff,
      stringsAsFactors = FALSE
    )
  }))
  
  
  std_reg <- merge(std_reg, neff_df, by = "lhs", all.x = TRUE)
  
  
  n <- nobs(sem_mod)
  std_reg$se_corrected <- std_reg$se * sqrt(n / std_reg$n.eff)
  std_reg$z_corrected  <- std_reg$est.std / std_reg$se_corrected
  std_reg$p_corrected  <- 2 * pnorm(abs(std_reg$z_corrected), lower.tail = FALSE)
  
  
  std_reg <- std_reg[, c("lhs", "op", "rhs",
                         "est.std", "se", "se_corrected",
                         "z_corrected", "p_corrected", "n.eff")]
  colnames(std_reg) <- c("lhs", "op", "rhs",
                         "std_coef", "se_orig", "se_corrected",
                         "z_corrected", "p_corrected", "n.eff")
  
  return(std_reg)
}


std_corrected <- get_spatial_corrected_std(sem_mod, spatial_result)


build_adj_matrix_from_std <- function(std_corrected) {
  
  all_vars <- unique(c(std_corrected$lhs, std_corrected$rhs))
  
  adj_matrix <- matrix(
    0,
    nrow = length(all_vars),
    ncol = length(all_vars),
    dimnames = list(all_vars, all_vars)
  )
  
  for (i in seq_len(nrow(std_corrected))) {
    from <- std_corrected$rhs[i]   
    to   <- std_corrected$lhs[i]   
    adj_matrix[from, to] <- std_corrected$std_coef[i]
  }
  
  return(adj_matrix)
}

adj_matrix <- build_adj_matrix_from_std(std_corrected)



find_all_paths <- function(adj_matrix, start_var,
                           end_var = "Carbon_Mg_ha", max_depth = 10) {
  g <- graph_from_adjacency_matrix(adj_matrix, weighted = TRUE, mode = "directed")
  all_paths <- all_simple_paths(g, from = start_var, to = end_var, mode = "out")
  
  all_paths <- all_paths[sapply(all_paths, length) > 2]
  all_paths <- all_paths[sapply(all_paths, length) <= max_depth]
  return(all_paths)
}

calculate_path_effect <- function(path, adj_matrix) {
  if (length(path) < 2) return(0)
  effect <- 1
  for (i in 1:(length(path) - 1)) {
    from_var <- names(path)[i]
    to_var   <- names(path)[i + 1]
    effect   <- effect * adj_matrix[from_var, to_var]
  }
  return(effect)
}

calc_indirect_path <- function(var_name, adj_matrix,
                               end_var = "Carbon_Mg_ha") {
  if (!var_name %in% rownames(adj_matrix)) return(NA)
  
  all_paths <- find_all_paths(adj_matrix, var_name, end_var)
  
  if (length(all_paths) == 0) return(0)
  
  total_effect  <- 0
  path_details  <- list()
  
  for (i in seq_along(all_paths)) {
    path   <- all_paths[[i]]
    effect <- calculate_path_effect(path, adj_matrix)
    if (!is.na(effect) && effect != 0) {
      total_effect  <- total_effect + effect
      path_details[[length(path_details) + 1]] <- list(
        path   = paste(names(path), collapse = " -> "),
        effect = effect
      )
    }
  }
  
  if (length(path_details) > 0) {
    cat("\n=== Indirect paths for", var_name, "===\n")
    for (d in path_details) {
      cat(sprintf("  %s : %.6f\n", d$path, d$effect))
    }
  }
  
  return(total_effect)
}

variables_to_test <- c(
  "shannon_dbhR", "STDAGE", "FD", "CWM",
  "CWD", "PRE", "bdod", "nitrogen_mean", "TEM"
)


indirect_effects <- list()
for (var in variables_to_test) {
  cat("\n>>> 变量:", var, "\n")
  indirect_effects[[var]] <- calc_indirect_path(var, adj_matrix)
  cat("Total indirect effect:", indirect_effects[[var]], "\n")
}


direct_effects_df <- std_corrected[std_corrected$lhs == "Carbon_Mg_ha",
                                   c("rhs", "std_coef", "p_corrected")]
colnames(direct_effects_df) <- c("Effect", "Value", "p_corrected")
direct_effects_df$mod <- "Direct effects"



indirect_effects_df <- data.frame(
  Effect      = variables_to_test,
  Value       = unlist(indirect_effects[variables_to_test]),
  p_corrected = NA_real_,    
  mod         = "Indirect effects",
  stringsAsFactors = FALSE
)

sens_gg_lavaan <- rbind(direct_effects_df, indirect_effects_df)
sens_gg_lavaan$DomSpecies <- "220_FDis"
sens_gg_lavaan$highlight  <- 0
sens_gg_lavaan$highlight  <- factor(sens_gg_lavaan$highlight)
sens_gg_lavaan$mod        <- factor(
  sens_gg_lavaan$mod,
  levels = c("Direct effects", "Indirect effects")
)



write.csv(sens_gg_lavaan,'/data/sens_gg2.csv', row.names = FALSE)




################################################################################
################################################################################

##===============================================================================
## Fig.1c,d  ---------
##===============================================================================
#### c：FD########################################################################
library(lavaan)

extract_sem_lines_lavaan <- function(sem_mod,
                                     spatial_result,
                                     data,
                                     focal     = "FD",
                                     mediator  = "shannon_dbhR",
                                     outcome   = "Carbon_Mg_ha",
                                     n_sim     = 10000,
                                     Climate_label = "CWD",   
                                     seed      = 123) {
  
  set.seed(seed)
  n <- lavaan::nobs(sem_mod)
  
  
  pe <- lavaan::parameterEstimates(sem_mod)
  get_est <- function(lhs, rhs) {
    r <- pe$est[pe$lhs == lhs & pe$op == "~" & pe$rhs == rhs]
    if (length(r) == 0) 0 else as.numeric(r[1])
  }
  a   <- get_est(mediator, focal)     # FD -> shannon_dbhR
  b   <- get_est(outcome,  mediator)  # shannon_dbhR -> Carbon
  c_p <- get_est(outcome,  focal)     # FD -> Carbon (direct)
  
  
  nm <- c(paste0(mediator, "~", focal),     # a
          paste0(outcome,  "~", mediator),  # b
          paste0(outcome,  "~", focal))     # c'
  V  <- lavaan::vcov(sem_mod)
  
  if (all(nm %in% rownames(V))) {
    Sigma <- V[nm, nm]
  } else {
    warning("")
    diag_v <- sapply(nm, function(x) if (x %in% rownames(V)) V[x, x] else 0)
    Sigma  <- diag(diag_v, 3, 3)
  }
  
 
  neff_of <- function(resp) {
    mi <- spatial_result$Morans_I[[resp]]
    if (is.null(mi)) return(n)
    ne <- mi$n.eff
    if (is.null(ne) || length(ne) == 0 || is.na(ne)) return(n)
    as.numeric(ne)[1]
  }
  k <- c(n / neff_of(mediator),   
         n / neff_of(outcome),    
         n / neff_of(outcome))   
  

  Sigma_corr <- Sigma * outer(sqrt(k), sqrt(k))
  

  sims   <- MASS::mvrnorm(n_sim, mu = c(a, b, c_p), Sigma = Sigma_corr)
  a_sim  <- sims[, 1]
  b_sim  <- sims[, 2]
  cp_sim <- sims[, 3]
  
  indirect_sim <- a_sim * b_sim         
  total_sim    <- cp_sim + indirect_sim  
  
 
  direct_slope <- c_p
  total_slope  <- mean(total_sim)
  
  
  indirect_mean <- mean(indirect_sim)
  indirect_ci   <- quantile(indirect_sim, c(0.025, 0.975))
  p_mc          <- 2 * min(mean(indirect_sim > 0), mean(indirect_sim < 0))
  total_ci      <- quantile(total_sim, c(0.025, 0.975))
  
  
  va  <- Sigma_corr[1, 1]; vb <- Sigma_corr[2, 2]; cab <- Sigma_corr[1, 2]
  se_indirect_delta <- sqrt(b^2 * va + a^2 * vb + 2 * a * b * cab)
  z_delta <- if (se_indirect_delta > 0) (a * b) / se_indirect_delta else NA
  p_delta <- if (!is.na(z_delta)) 2 * pnorm(abs(z_delta), lower.tail = FALSE) else NA
  
  cat(sprintf(" ",
              Climate_label, focal, mediator, outcome))
  cat(sprintf("  a (%s->%s)         = %.4f\n", focal, mediator, a))
  cat(sprintf("  b (%s->%s)   = %.4f\n", mediator, outcome, b))
  cat(sprintf("  c'(%s->%s)  = %.4f\n", focal, outcome, c_p))
  cat(sprintf("  ----------------------------------------------\n"))
  cat(sprintf("     = %.4f\n", direct_slope))
  cat(sprintf("    = %.4f  [95%% CI %.4f, %.4f]\n",
              indirect_mean, indirect_ci[1], indirect_ci[2]))
  cat(sprintf("        p(MC)=%.4g   p(delta)=%.4g\n", p_mc, p_delta))
  cat(sprintf("   (c'+a*b)         = %.4f  [95%% CI %.4f, %.4f]\n",
              total_slope, total_ci[1], total_ci[2]))
  
  
  x_vec   <- data[[focal]]
  x_span  <- max(x_vec, na.rm = TRUE) - min(x_vec, na.rm = TRUE)
  x_range <- seq(0, x_span, length.out = 100)
  

  
  direct_ci <- quantile(
    cp_sim,
    c(0.025,0.975)
  )
  
  total_ci <- quantile(
    total_sim,
    c(0.025,0.975)
  )
  
  
  direct_width <- max(
    abs(direct_ci - direct_slope)
  )
  
  
  total_width <- max(
    abs(total_ci - total_slope)
  )
  
  
  
  lines_df <- rbind(
    
    
    data.frame(
      
      FD_increase=x_range,
      
      Response=direct_slope*x_range,
      
      lower=direct_slope*x_range -
        max(abs(direct_ci-direct_slope)),
      
      upper=direct_slope*x_range +
        max(abs(direct_ci-direct_slope)),
      
      Group="Direct effect"
      
    ),
    
    
    
    data.frame(
      
      FD_increase=x_range,
      
      Response=total_slope*x_range,
      
      lower=total_slope*x_range -
        max(abs(total_ci-total_slope)),
      
      upper=total_slope*x_range +
        max(abs(total_ci-total_slope)),
      
      Group="Direct + indirect"
      
    )
    
  )
  
  summary_df <- data.frame(
    focal = focal, mediator = mediator, outcome = outcome,
    a = a, b = b, c_prime = c_p,
    direct        = direct_slope,
    indirect      = indirect_mean,
    indirect_lwr  = unname(indirect_ci[1]),
    indirect_upr  = unname(indirect_ci[2]),
    p_indirect_MC    = p_mc,
    p_indirect_delta = p_delta,
    total         = total_slope,
    total_lwr     = unname(total_ci[1]),
    total_upr     = unname(total_ci[2]),
    n = n, n_eff_mediator = neff_of(mediator), n_eff_outcome = neff_of(outcome),
    row.names = NULL
  )
  
  list(lines = lines_df, summary = summary_df)
}

res <- extract_sem_lines_lavaan(
  sem_mod        = sem_mod,
  spatial_result = spatial_result,
  data           = data_sem2,
  focal          = "FD",
  mediator       = "shannon_dbhR",
  outcome        = "Carbon_Mg_ha",
  n_sim          = 10000,
  Climate_label      = "CWD"      
)

final_df <- res$lines
print(head(final_df))
print(res$summary)





library(ggplot2)


plot_one_category_FD <- function(lines_df, summary_df,
                                 title    = "",
                                 subtitle = "") {
  df <- lines_df[lines_df$Response != 0, ]   
  
  direct_slope <- summary_df$direct     
  total_slope  <- summary_df$total    
  
  
  x_max     <- max(df$FD_increase, na.rm = TRUE)
  x_lab_pos <- x_max * 0.55
  y_direct  <- direct_slope * x_lab_pos
  y_total   <- total_slope  * x_lab_pos
  
  ggplot(df, aes(x = FD_increase,
                 y = Response,
                 color = Group,
                 linetype = Group)) +
    
    geom_ribbon(
      aes(
        ymin = lower,
        ymax = upper,
        fill = Group
      ),
      alpha = 0.20,
      color = NA
    ) +
    
    geom_line(
      size = 0.8
    ) +
    
    annotate("text",
             x = x_lab_pos, y = y_direct,
             label  = sprintf("%.2f", direct_slope),
             family = "serif", size = 4.5,
             hjust  = 1.15, vjust = -0.4) +
    
    annotate("text",
             x = x_lab_pos, y = y_total,
             label  = sprintf("%.2f", total_slope),
             family = "serif", size = 4.5,
             hjust  = -0.15, vjust = 1.2) +
    scale_color_manual(values = c("Direct effect"     = "black",
                                  "Direct + indirect" = "black")) +
    scale_linetype_manual(values = c("Direct effect"     = "dashed",
                                     "Direct + indirect" = "solid")) +
    scale_fill_manual(
      values = c(
        "Direct effect"     = "grey70",
        "Direct + indirect" = "grey70"
      )
    ) +
    guides(
      fill = "none",
      color = guide_legend(order = 1),
      linetype = guide_legend(order = 1)
    )+
    labs(x = "Functional diversity",
         y = expression(log("AGC ("*Mg~C~ha^{-1}*")")),
         title = title, subtitle = subtitle,
         color = NULL, linetype = NULL) +
    theme_bw() +
    theme(
      text             = element_text(family = "serif", size = 14),
      legend.position  = c(0.25, 0.85),
      panel.grid       = element_blank(),
      panel.background = element_blank(),
      plot.background  = element_rect(fill = "white", color = NA),
      axis.line        = element_blank(),#element_line(color = "black", size = 0.5),
      axis.ticks.x     = element_line(size = 0.5),
      axis.ticks.y     = element_line(size = 0.5),
      axis.text.x      = element_text(size = 14, family = "serif", color = "black"),
      axis.text.y      = element_text(size = 14, family = "serif", color = "black"),
      axis.title       = element_text(size = 14, family = "serif", color = "black"),
      strip.text       = element_text(size = 14, family = "serif"),
      strip.background = element_blank()
    )
}
p_single <- plot_one_category_FD(res$lines, res$summary)
print(p_single)


#################################################################################
#### d：CWM########################################################################
library(lavaan)


extract_sem_lines_lavaan <- function(sem_mod,
                                     spatial_result,
                                     data,
                                     focal     = "CWM",
                                     mediator  = "shannon_dbhR",
                                     outcome   = "Carbon_Mg_ha",
                                     n_sim     = 10000,
                                     Climate_label = "CWD",  
                                     seed      = 123) {
  
  set.seed(seed)
  n <- lavaan::nobs(sem_mod)
  
  
  pe <- lavaan::parameterEstimates(sem_mod)
  get_est <- function(lhs, rhs) {
    r <- pe$est[pe$lhs == lhs & pe$op == "~" & pe$rhs == rhs]
    if (length(r) == 0) 0 else as.numeric(r[1])
  }
  a   <- get_est(mediator, focal)     # CWM -> shannon_dbhR
  b   <- get_est(outcome,  mediator)  # shannon_dbhR -> Carbon
  c_p <- get_est(outcome,  focal)     # CWM -> Carbon (direct)
  
  
  nm <- c(paste0(mediator, "~", focal),     # a
          paste0(outcome,  "~", mediator),  # b
          paste0(outcome,  "~", focal))     # c'
  V  <- lavaan::vcov(sem_mod)
  
  if (all(nm %in% rownames(V))) {
    Sigma <- V[nm, nm]
  } else {
    warning(" ")
    diag_v <- sapply(nm, function(x) if (x %in% rownames(V)) V[x, x] else 0)
    Sigma  <- diag(diag_v, 3, 3)
  }
  
  
  neff_of <- function(resp) {
    mi <- spatial_result$Morans_I[[resp]]
    if (is.null(mi)) return(n)
    ne <- mi$n.eff
    if (is.null(ne) || length(ne) == 0 || is.na(ne)) return(n)
    as.numeric(ne)[1]
  }
  k <- c(n / neff_of(mediator),   
         n / neff_of(outcome),    
         n / neff_of(outcome))   
  
  
  Sigma_corr <- Sigma * outer(sqrt(k), sqrt(k))
  
  
  sims   <- MASS::mvrnorm(n_sim, mu = c(a, b, c_p), Sigma = Sigma_corr)
  a_sim  <- sims[, 1]
  b_sim  <- sims[, 2]
  cp_sim <- sims[, 3]
  
  indirect_sim <- a_sim * b_sim       
  total_sim    <- cp_sim + indirect_sim  
  
  
  direct_slope <- c_p
  total_slope  <- mean(total_sim)
  
  
  indirect_mean <- mean(indirect_sim)
  indirect_ci   <- quantile(indirect_sim, c(0.025, 0.975))
  p_mc          <- 2 * min(mean(indirect_sim > 0), mean(indirect_sim < 0))
  total_ci      <- quantile(total_sim, c(0.025, 0.975))
  
  
  va  <- Sigma_corr[1, 1]; vb <- Sigma_corr[2, 2]; cab <- Sigma_corr[1, 2]
  se_indirect_delta <- sqrt(b^2 * va + a^2 * vb + 2 * a * b * cab)
  z_delta <- if (se_indirect_delta > 0) (a * b) / se_indirect_delta else NA
  p_delta <- if (!is.na(z_delta)) 2 * pnorm(abs(z_delta), lower.tail = FALSE) else NA
  
  cat(sprintf("\n[%s] %s\n",
              Climate_label, focal, mediator, outcome))
  cat(sprintf("  a (%s->%s)         = %.4f\n", focal, mediator, a))
  cat(sprintf("  b (%s->%s)   = %.4f\n", mediator, outcome, b))
  cat(sprintf("  c'(%s->%s )  = %.4f\n", focal, outcome, c_p))
  cat(sprintf("  ----------------------------------------------\n"))
  cat(sprintf("              = %.4f\n", direct_slope))
  cat(sprintf("              = %.4f  [95%% CI %.4f, %.4f]\n",
              indirect_mean, indirect_ci[1], indirect_ci[2]))
  cat(sprintf("        p(MC)=%.4g   p(delta)=%.4g\n", p_mc, p_delta))
  cat(sprintf("   (c'+a*b)         = %.4f  [95%% CI %.4f, %.4f]\n",
              total_slope, total_ci[1], total_ci[2]))
  
  
  x_vec   <- data[[focal]]
  x_span  <- max(x_vec, na.rm = TRUE) - min(x_vec, na.rm = TRUE)
  x_range <- seq(0, x_span, length.out = 100)
  
  direct_ci <- quantile(
    cp_sim,
    c(0.025,0.975)
  )
  
  total_ci <- quantile(
    total_sim,
    c(0.025,0.975)
  )
  
  
  direct_width <- max(
    abs(direct_ci - direct_slope)
  )
  
  
  total_width <- max(
    abs(total_ci - total_slope)
  )
  
  
  
  lines_df <- rbind(
    
    
    data.frame(
      
      CWM_increase=x_range,
      
      Response=direct_slope*x_range,
      
      lower=direct_slope*x_range -
        max(abs(direct_ci-direct_slope)),
      
      upper=direct_slope*x_range +
        max(abs(direct_ci-direct_slope)),
      
      Group="Direct effect"
      
    ),
    
    
    
    data.frame(
      
      CWM_increase=x_range,
      
      Response=total_slope*x_range,
      
      lower=total_slope*x_range -
        max(abs(total_ci-total_slope)),
      
      upper=total_slope*x_range +
        max(abs(total_ci-total_slope)),
      
      Group="Direct + indirect"
      
    )
    
  )
  
  summary_df <- data.frame(
    focal = focal, mediator = mediator, outcome = outcome,
    a = a, b = b, c_prime = c_p,
    direct        = direct_slope,
    indirect      = indirect_mean,
    indirect_lwr  = unname(indirect_ci[1]),
    indirect_upr  = unname(indirect_ci[2]),
    p_indirect_MC    = p_mc,
    p_indirect_delta = p_delta,
    total         = total_slope,
    total_lwr     = unname(total_ci[1]),
    total_upr     = unname(total_ci[2]),
    n = n, n_eff_mediator = neff_of(mediator), n_eff_outcome = neff_of(outcome),
    row.names = NULL
  )
  
  list(lines = lines_df, summary = summary_df)
}

res <- extract_sem_lines_lavaan(
  sem_mod        = sem_mod,
  spatial_result = spatial_result,
  data           = data_sem2,
  focal          = "CWM",
  mediator       = "shannon_dbhR",
  outcome        = "Carbon_Mg_ha",
  n_sim          = 10000,
  Climate_label      = "CWD"       
)

final_df <- res$lines
print(head(final_df))
print(res$summary)


library(ggplot2)


plot_one_category_CWM <- function(lines_df, summary_df,
                                  title    = "",
                                  subtitle = "") {
  df <- lines_df[lines_df$Response != 0, ]   
  
  direct_slope <- summary_df$direct   
  total_slope  <- summary_df$total   
  
  
  x_max     <- max(df$CWM_increase, na.rm = TRUE)
  x_lab_pos <- x_max * 0.55
  y_direct  <- direct_slope * x_lab_pos
  y_total   <- total_slope  * x_lab_pos
  
  ggplot(df, aes(x = CWM_increase,
                 y = Response,
                 color = Group,
                 linetype = Group)) +
    
    geom_ribbon(
      aes(
        ymin = lower,
        ymax = upper,
        fill = Group
      ),
      alpha = 0.20,
      color = NA
    ) +
    
    geom_line(
      size = 0.8
    ) +
    
    annotate("text",
             x = x_lab_pos, y = y_direct,
             label  = sprintf("%.2f", direct_slope),
             family = "serif", size = 4.5,
             hjust  = 1.15, vjust = -0.4) +
    
    annotate("text",
             x = x_lab_pos, y = y_total,
             label  = sprintf("%.2f", total_slope),
             family = "serif", size = 4.5,
             hjust  = -0.15, vjust = 1.2) +
    scale_color_manual(values = c("Direct effect"     = "black",
                                  "Direct + indirect" = "black")) +
    scale_linetype_manual(values = c("Direct effect"     = "dashed",
                                     "Direct + indirect" = "solid")) +
    scale_fill_manual(
      values = c(
        "Direct effect"     = "grey70",
        "Direct + indirect" = "grey70"
      )
    ) +
    guides(
      fill = "none",
      color = guide_legend(order = 1),
      linetype = guide_legend(order = 1)
    )+
    
    labs(x = "CWM PC1",
         y = expression(log(AGC~(Mg~C~ha^{-1}))),
         title = title, subtitle = subtitle,
         color = NULL, linetype = NULL) +
    
    
    theme_bw() +
    theme(
      text             = element_text(family = "serif", size = 14),
      legend.position  = c(0.25, 0.85),
      panel.grid       = element_blank(),
      panel.background = element_blank(),
      plot.background  = element_rect(fill = "white", color = NA),
      axis.line        = element_blank(),#element_line(color = "black", size = 0.5),
      axis.ticks.x     = element_line(size = 0.5),
      axis.ticks.y     = element_line(size = 0.5),
      axis.text.x      = element_text(size = 14, family = "serif", color = "black"),
      axis.text.y      = element_text(size = 14, family = "serif", color = "black"),
      axis.title       = element_text(size = 14, family = "serif", color = "black"),
      strip.text       = element_text(size = 14, family = "serif"),
      strip.background = element_blank()
    )
}
p_single1 <- plot_one_category_CWM(res$lines, res$summary)
print(p_single1)
#####merge#######################################################################
library(ggpubr)
fig3 = ggarrange(p_single, p_single1,
                 ncol = 1, nrow = 2,
                 heights = c(1, 1), widths = c(1, 1),
                 common.legend = FALSE, align = "hv")
fig3

ggsave("fig.png",
       path = "/Figure/",
       width = 4, height = 8, units = "in",
       dpi = 600, plot = fig3)


#################################################################################
#################################################################################
##===============================================================================
## Fig.1b ---------
##===============================================================================
rm(list=ls())

sens_gg=read.csv('E:/FIA_DATA/森林演替/data/sens_gg2.csv' )
colnames(sens_gg)

sens_gg <- sens_gg[sens_gg$DomSpecies %in% c("220_FDis"), ]
library(dplyr)

net_effect <- sens_gg %>%
  group_by(Effect) %>% 
  dplyr::summarise(
    Value = sum(Value),    
    mod = "Total effects",
    .groups = "drop"
  )
sens_gg <- bind_rows(sens_gg, net_effect)


library(ggplot2)
####Define figure theme####
number_ticks <- function(n) {function(limits) pretty(limits, n)}
my.formula <- y ~ x      ####define linear fit in fig
Ftheme<-theme(              axis.line = element_line(colour = "black"),
                            panel.grid.major = element_blank(),
                            panel.grid.minor = element_blank(),
                            panel.border = element_rect(colour = NA,fill=NA),
                            panel.background = element_blank(),
                            plot.title = element_text(size = 12,family="serif"),
                            axis.text = element_text(size = 12,family="serif"),
                            axis.title = element_text(size = 12,family="serif"),
                            text = element_text(size = 12,family="serif"))+
  theme(strip.background = element_blank(),strip.placement = "in")#set a theme for Fig

#########
#sens_gg <- sens_gg[!sens_gg$Effect %in% c("CWM2", "INVYR_1"), ]
unique(net_effect$Effect)
net_effect$Effect<-factor(net_effect$Effect,
                          
                          
                          levels = c(
                            
                            "shannon_dbhR",
                            "FD","CWM", "STDAGE","CWD" ,
                            "PRE","TEM","bdod","nitrogen_mean"
                          ),
                          labels = c(
                            "Structural\ndiversity" ,  "Functional\ndiversity" ,"CWM PC1" , "Stand age",
                            "CWD",  "MAP", "MAT" ,"Soil bulk\ndensity","Soil\nnitrogen"
                          )) 
sens_gg$Effect<-factor(sens_gg$Effect,
                       
                       
                       levels = c(
                         
                         "shannon_dbhR",
                         "FD","CWM", "STDAGE" ,"CWD",
                         "PRE","TEM","bdod","nitrogen_mean"
                       ),
                       labels = c(
                         "Structural\ndiversity" ,  "Functional\ndiversity" ,
                         "CWM PC1", "Stand age", "CWD" , "MAP", "MAT" ,"Soil bulk\ndensity","Soil\nnitrogen"
                       ))







sens_plot=sens_gg

sens_plot$mod <- factor(
  sens_plot$mod,
  levels = c("Direct effects", "Indirect effects", "Total effects")
)
library(tidyr)
sens_plot <- sens_plot %>%
  complete(Effect, mod, fill = list(Value = 0))

P_A_220 = ggplot(sens_plot, 
                 aes(x = Effect,
                     y = Value,
                     fill = mod)) + 
  # coord_flip() +
  geom_col(
    position = position_dodge(width = 0.9), 
    width = 0.7,
    color = "black",
    linewidth = 0.3
  ) +
  
  scale_fill_manual(
    values = c(
      "Direct effects"   = "#ffc2e5",  
      "Indirect effects" = "#ebffac", 
      "Total effects"      =  "#76daff"   
    ),
    name = "Effect Type"
  ) +
  
  geom_hline(yintercept = 0, linetype = "longdash", color = "gray50", size = 0.5) +
  
  scale_y_continuous(limits = c(-0.5, 0.7),
                     breaks = c( -0.4,-0.2,0,0.2, 0.4,0.6)) +
  
  labs(x = "", y = expression(paste("Standardized effect size (", italic(beta), ")"))) +
  
  theme_bw() +
  
  theme(
    legend.position = c(0.8, 0.8),
    legend.text = element_text(size = 12, family = "serif", color = "black"),
    legend.title = element_blank(),
    panel.grid = element_blank(),
    legend.key.size = unit(0.4, "cm"),
    axis.text.x = element_text(
      angle = 90,
      hjust = 0.5,
      vjust = 0.5,
      size = 12,
      family = "serif", color = "black"
    ),
    
    axis.text.y = element_text(size = 12, family = "serif", color = "black"),
    axis.title = element_text(size = 12, family = "serif", color = "black"),
    # axis.line.y = element_line(size = 0.5),
    # axis.line.x = element_line(size = 0.5)#,
    axis.line = element_blank( )
    #axis.ticks.y= element_blank(),
    # panel.border = element_blank()
  )
P_A_220

ggsave(
  "fig.png",
  path = "/Figure/",
  width = 5,
  height = 3,
  units = "in",
  dpi = 600,
  plot = P_A_220
)


