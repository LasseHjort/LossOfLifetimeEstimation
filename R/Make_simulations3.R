##################### Simulations for loss of lifetime paper #####################

#Set up for linux computer
if(any(grepl("Linux|linux", Sys.info()))){
  project <- "~/LossOfLifetime/"
  setwd(project)
  
  #Load libraries
  library(rstpm2)
  library(cuRe)
  library(ggplot2)
  library(matrixStats)
  library(xtable)
  library(parallel)
  library(VGAM)
  
  #Figure and table directories
  fig.out <- "."
  tab.out <- "."
  data.out <- "."
  thesis <- FALSE
  
  #Table format
  tab.format <- "%.3f"
  
  #Colour palette for models
  cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  
  #Global project settings
  mai_par <- par("mai")
  mfrow_par <- par("mfrow")
  format <- "%.1f"
  ayear <- 365.24
}

get.knots <- function(data, k){
  quantile(data$FU_years[data$status ==1], seq(0, 1, length.out = k))
}

#Function for simulating survival data
sim_surv <- function(pars, age, n, type = "weibull"){
  #Generate general population survival. As implemented, age, gender, and year are fixed.
  sim_ages <- rnorm(n, mean = age, sd = 0)
  sim_gender <- factor(rbinom(n, size = 1, prob = 1), levels = 0:1, labels = c("male", "female"))
  dates <- as.Date("1980-01-01") #+ 1:as.numeric(as.Date("1980-01-01") - as.Date("1980-01-01"))
  diag_dates <- sample(dates, n, replace = T)
  D <- data.frame(age = sim_ages * ayear, sex = sim_gender, diag_date = diag_dates)
  S_exp <- survexp(~1, rmap = list(age = age, sex = sex, year = diag_date), 
                   data = D,
                   times = seq(0, (120 - age) * ayear, length.out = 1000),
                   ratetable = survexp.dk, 
                   scale = ayear)
  
  #Function to draw simulations from
  len <- length(S_exp$time)
  sim_fun_gen <- function(x){
    if(x > 1 - S_exp$surv[len]){
      S_exp$time[len]
    }else{
      S_exp$time[findInterval(x, vec = 1 - S_exp$surv) + 1] 
    }
  }
  
  
  #Set survival of the uncured and relative survival functions
  if(type == "weibull"){
    surv_can_fun <- function(pars, time) exp(-pars[3] * time ^ pars[2]) 
    dens_can_fun <- function(pars, time) exp(-pars[3] * time ^ pars[2]) * pars[3] * pars[2] * time ^ (pars[2] - 1)
    rel_surv <- function(time) pars[1] + (1 - pars[1]) * exp(-pars[3] * time ^ pars[2])
    qrel_surv <- function(q){
      res <- rep(Inf, length(q))
      wh <- q < 1 - pars[1]
      res[wh] <- qweibull(q[wh] / (1 - pars[1]), 
                          shape = pars[2], 
                          scale = 1 / pars[3] ^ (1 / pars[2]))
      res
    }
  }
  if(type == "gengamma"){
    surv_can_fun <- function(pars, time) pgengamma.stacy(time, scale = pars[2], 
                                                         d = pars[3], k = pars[4], lower.tail = F)
    dens_can_fun <- function(pars, time) dgengamma.stacy(time, scale = pars[2], 
                                                         d = pars[3], k = pars[4])
    rel_surv <- function(time) pars[1] + (1 - pars[1]) * pgengamma.stacy(time, scale = pars[2], 
                                                                         d = pars[3], k = pars[4], lower.tail = F)
    qrel_surv <- function(q){
      res <- rep(Inf, length(q))
      wh <- q < 1 - pars[1]
      res[wh] <- qgengamma.stacy(q[wh] / (1 - pars[1]), 
                                 scale = pars[2], 
                                 d = pars[3], 
                                 k = pars[4])
      res
    }
  }
  
  # unifind <- function(time, u) surv_can_fun(pars = pars, time) - u
  # 
  # sim_fun_rel <- function(x){
  #   res <- rep(Inf, length(x))
  #   for(i in 1:length(x)){
  #     #cat(i, "\n")
  #     if(x[i] >= pars[1]){
  #       u.new <- (x[i] - pars[1]) /(1 - pars[1])
  #       eval <- unifind(1e-15, u.new)
  #       if(eval < 0){
  #         res[i] <- .Machine$double.eps
  #       } else {
  #         res[i] <- uniroot(unifind, interval = c(1e-15, 300), u = u.new)$root 
  #       }
  #     }
  #   }
  #   res
  # }
  
  #Simulate uniform variable for general population and disease specific survival
  uni_sim1 <- runif(nrow(D))
  uni_sim2 <- runif(nrow(D))
  
  #Simulate from both distributions
  sim_gen <- sapply(uni_sim1, sim_fun_gen)
  sim_pop <- qrel_surv(uni_sim2)
  #sim_pop <- sim_fun_rel(uni_sim2)
  D$fu <- pmin(sim_gen, sim_pop)
  
  #Simulate from censoring distribution
  #Set parameters
  max <- 15
  sim_cens <- runif(n = nrow(D), min = 0, max = max)
  
  #Generated follow-up as the minimum of survival time and censoring time
  D$FU <- pmin(D$fu, sim_cens)
  D$status <- as.numeric(D$fu <= sim_cens)
  D$FU[D$FU < 1e-3] <- 1e-3
  # plot(survfit(Surv(FU, status) ~ 1, data = D))
  # curve(rel_surv, from = 0, to = 15, add = T, col = 2)
  #Follow-up in days
  D$FU <- D$FU *  ayear
  # rsfit <- rs.surv(Surv(FU, status) ~ 1 + ratetable(age = age, sex = sex, year = diag_date), 
  #                  data = D, ratetable = survexp.dk, method = "ederer2")
  # 
  # plot(rsfit)
  # abline(h = 0.5)
  D$FU_years <- D$FU / ayear
  #Get general population hazard
  D$exp_haz <- general.haz(time = "FU", age = "age", sex = "sex", year = "diag_date", 
                           data = D, ratetable = survexp.dk)
  list(D = D, pars = pars, rel_surv = rel_surv, S_exp = S_exp, 
       surv_can_fun = surv_can_fun, dens_can_fun = dens_can_fun)
}


gaussxw <- statmod::gauss.quad(100)

#Function for calculating the difference in loss of lifetime
calc_AUCs <- function(res_list, plot.haz = F, plot.surv = F){
  #Extract model parameters
  S_exp <- res_list$S_exp
  pars <- res_list$pars
  rel_surv <- res_list$rel_surv
  D <- res_list$D
  surv_can_fun <- res_list$surv_can_fun
  dens_can_fun <- res_list$dens_can_fun
  
  # rsfit <- rs.surv(Surv(FU, status) ~ 1 + ratetable(age = age, sex = sex, year = diag_date),
  #                  data = D, ratetable = survexp.dk, method = "ederer2")
  # rsfit$time <- rsfit$time / ayear
  # plot(rsfit)
  # f <- function(t) pars[1] + (1 - pars[1]) * exp(-pars[3] * t ^ pars[2])
  # curve(f, from = 0, to = 15, add = T, col = 2)
  
  exp_function2 <- function(t){
    s <- summary(S_exp, t)
    names(s$surv) <- s$time
    a <- s$surv[as.character(t)]
    names(a) <- NULL
    a
  }
  
  #Determine tau by as the time point when the general population survival is less than 1e-8
  wh <- min(which(S_exp$surv < 1e-8))
  tau <- S_exp$time[wh]
  
  #Function for calculating the true loss of lifetime
  calc_int <- function(time){
    scale <- (tau - time) / 2 
    scale2 <- (tau + time) / 2
    eval_gen_t <- exp_function2(time)
    eval_rel_t <- rel_surv(time)
    eval <- rep(NA, length(time))
    for(i in 1:length(time)){
      points <- scale[i] * gaussxw$nodes + scale2[i]
      eval_gen <- exp_function2(points)
      eval_rel <- rel_surv(points)
      eval[i] <- sum(gaussxw$weights * (eval_gen - eval_gen * eval_rel / eval_rel_t[i]))
    }
    scale * eval / eval_gen_t 
  }
  
  times <- 0:15
  LOL <- calc_int(times)
  
  # #Determine time points for AUC bias curve
  # times <- unique(sort(c(seq(0, 15, length.out = 1000), 1:14)))
  # #Determine time points for LOL calculation
  # t_new <- sort(unique(c(times, seq(0, tau, length.out = 1000))), decreasing = T)
  # df_time <- -diff(t_new)
  # 
  # #Calculate the loss of lifetime for the model we simulate from
  # exp_eval <- exp_function2(t_new)
  # surv_eval <- rel_surv(t_new) * exp_eval
  # surv_diff <- diff(surv_eval)
  # inner <- abs(surv_diff) / 2 + pmin(surv_eval[-length(surv_eval)], surv_eval[-1])
  # vals_pop <- cumsum(c(0, inner * df_time))
  # surv_diff <- diff(exp_eval)
  # inner <- abs(surv_diff) / 2 + pmin(exp_eval[-length(exp_eval)], exp_eval[-1])
  # vals_exp <- cumsum(c(0, inner * df_time))
  # these <- t_new %in% times
  # LOL <- rev(vals_exp[these]) / rev(exp_eval[these]) - rev(vals_pop[these]) / rev(surv_eval[these])
  
  # vals_pop <- cumsum(c(0, rel_surv(t_new[-length(t_new)]) * exp_function(t_new[-length(t_new)]) * df_time))
  # vals_pop <- rev(vals_pop[t_new %in% times])
  # vals_exp <- cumsum(c(0, exp_function(t_new[-length(t_new)]) * df_time))
  # vals_exp <- rev(vals_exp[t_new %in% times])
  # LOL <- vals_exp / exp_function(times) - vals_pop / (rel_surv(times) * exp_function(times))
  
  #Fit model by Nelson et al. 2007
  knots_nelson <- log(sort(quantile(D$FU_years[D$status == 1], c(0, 0.2, 0.4, 0.6, 0.8, 1))))
  fit_nelson <- stpm2(Surv(FU_years, status == 1) ~ -1, data = D, bhazard = D$exp_haz, 
                      smooth.formula = ~ cb(x = log(FU_years), knots = knots_nelson))
  
  #Fit models by Andersson et al. 2011
  add.knot <- 10
  knots_andersson1 <- log(sort(c(quantile(D$FU_years[D$status ==1], 
                                          c(0, 0.2, 0.4, 0.6, 0.8, 1)), add.knot)))
  fit_andersson1 <- stpm2(Surv(FU_years, status) ~ -1, data = D, bhazard = D$exp_haz, 
                          smooth.formula = ~cbc(x = log(FU_years), knots = knots_andersson1))
  
  add.knots <- c(10, 80)
  knots_andersson2 <- log(sort(c(quantile(D$FU_years[D$status ==1], 
                                          c(0, 0.2, 0.4, 0.6, 0.8)), add.knots)))
  fit_andersson2 <- stpm2(Surv(FU_years, status) ~ -1, data = D, bhazard = D$exp_haz, 
                          smooth.formula = ~cbc(x = log(FU_years), knots = knots_andersson2))
  
  #Fit flexible mixture cure models
  knots <- log(get.knots(D, 5))
  fit_flex_mix1 <- GenFlexCureModel(Surv(FU_years, status) ~ -1, data = D, bhazard = "exp_haz", 
                                    smooth.formula = ~ cb(x = log(FU_years), knots = knots), 
                                    covariance = F, verbose = F, ini.types = "cure")
  knots <- log(c(min(D$FU_years[D$status == 1]), 0.5, 1, 2, 5))
  fit_flex_mix2 <- GenFlexCureModel(Surv(FU_years, status) ~ -1, data = D, bhazard = "exp_haz", 
                                    smooth.formula = ~ cb(x = log(FU_years), knots = knots), 
                                    covariance = F, verbose = F, ini.types = "cure")
  # fit_flex_mix2 <- GenFlexCureModel(Surv(FU_years, status) ~ 1, data = D, 
  #                                   bhazard = "exp_haz", 
  #                                   knots = log(c(min(D$FU_years), 0.5, 1, 2, 5)), 
  #                                   covariance = F, verbose = F)
  
  #Assemble models
  all_models <- list(fit_nelson, fit_andersson1, fit_andersson2, fit_flex_mix1, fit_flex_mix2)

  if(plot.haz){
    t <- seq(0.1, tau, length.out = 100)
    #Plot hazards
    s_u <- surv_can_fun(pars, t)
    f_u <- dens_can_fun(pars, t)
    h_e1 <- (1 - pi) * f_u / (pi + (1 - pi) * s_u)
    
    b <- basis(fit_sim_cure3$knots, log(t)) %*% fit_sim_cure3$coefs.spline
    db <- dbasis(fit_sim_cure3$knots, log(t)) %*% fit_sim_cure3$coefs.spline
    h_u <- exp(b) / t * db
    s_u <- exp(-exp(b))
    pi <- logistic(fit_sim_cure3$coefs)
    h_e2 <- (1 - pi) * s_u * h_u / (pi + (1 - pi) * s_u)
    
    b <- basis(fit_sim_cure4$knots, log(t)) %*% fit_sim_cure4$coefs.spline
    db <- dbasis(fit_sim_cure4$knots, log(t)) %*% fit_sim_cure4$coefs.spline
    h_u <- exp(b) / t * db
    s_u <- exp(-exp(b))
    pi <- logistic(fit_sim_cure4$coefs)
    h_e3 <- (1 - pi) * s_u * h_u / (pi + (1 - pi) * s_u)
    
    b <- basis(fit_sim$knots, log(t)) %*% fit_sim$coefs
    db <- dbasis(fit_sim$knots, log(t)) %*% fit_sim$coefs
    h_e4 <- exp(b) / t * db
    
    plot(h_e1 ~ t, type = "l", ylim = range(c(h_e1, h_e2, h_e3, h_e4), na.rm = T), 
         xlab = "Time", ylab = "Excess mortality")
    lines(h_e2 ~ t, col = 2)
    lines(h_e3 ~ t, col = 3)
    lines(h_e4 ~ t, col = 4)
    legend("topright", fill = 1:4, legend = c("True", "CureModelknots", "CureModel", "FlexRel"))
  }
  
  if(plot.surv){
    t <- seq(0.1, tau, length.out = 100)
    s_u1 <- rel_surv(t)
    
    b <- basis(fit_sim_cure3$knots, log(t)) %*% fit_sim_cure3$coefs.spline
    pi <- logistic(fit_sim_cure3$coefs)
    s_u2 <- pi + (1 - pi) * exp(-exp(b))
    
    b <- basis(fit_sim_cure4$knots, log(t)) %*% fit_sim_cure4$coefs.spline
    pi <- logistic(fit_sim_cure4$coefs)
    s_u3 <- pi + (1 - pi) * exp(-exp(b))
    
    b <- basis(fit_sim$knots, log(t)) %*% fit_sim$coefs
    s_u4 <- exp(-exp(b))
    
    plot(s_u1 ~ t, type = "l", ylim = c(0, 1), ylab = "Relative survival", xlab = "Time")
    lines(s_u2 ~ t, col = 2)
    lines(s_u3 ~ t, col = 3)
    lines(s_u4 ~ t, col = 4)
    legend("topright", fill = 1:4, legend = c("True", "CureModelknots", "CureModel", "FlexRel"))
    
  }
  
  expected <- survexp(formula = ~ 1, rmap = list(year = diag_date),
                      data = D, ratetable = survexp.dk,
                      scale = ayear, times = seq(0, tau + 1, length.out = 1000) * ayear)
  
  # smooth.obj <- smooth.spline(x = expected$time, y = expected$surv, all.knots = T)
  # exp.fun <- function(time) predict(smooth.obj, x = time)$y
  
  exp.fun <- function(time){
    s <- summary(expected, time)
    names(s$surv) <- s$time
    survs <- s$surv[as.character(time)]
    names(survs) <- NULL
    survs
  }
  
  #Calculate the loss of lifetime
  LL_res <- lapply(all_models, function(fit){
    calc.LL(fit, time = times, tau = tau, var.type = "n", 
            rmap = list(year = diag_date), exp.fun = list(exp.fun))[[1]]
  })
  
  #Extract loss of lifetime estimates at certain time points
  LOL_combined <- do.call(cbind, LL_res)
  
  #Calculate loss of lifetime biases
  bias <- LOL_combined - LOL
  #bias <- lapply(LOL_combined, function(x) x - LOL)
  
  scale <- 15 / 2
  times_int <- scale * gaussxw$nodes + scale
  LL_res_t <- lapply(all_models, function(fit){
    calc.LL(fit, time = times_int, tau = tau, var.type = "n", 
            rmap = list(year = diag_date), exp.fun = list(exp.fun))[[1]]
  })
  TrueLOL_int <- calc_int(times_int)
  LL_res_t <- do.call(cbind, LL_res_t)
  abs_bias <- abs(LL_res_t - TrueLOL_int)
  AUCs <- scale * colSums(gaussxw$weights * abs_bias)
  
  #Calculate integrated loss of lifetime bias
  #dif_t <- diff(times)
  #abs_bias <- lapply(bias, function(x) abs(x))
  
  #AUCs <- sapply(abs_bias, function(x){
  #  inner <- abs(diff(x)) / 2 + pmin(x[-length(x)], x[-1])
  #  sum(c(0, inner * dif_t))
  #})
  
  #bias <- do.call(cbind, bias)
  #Output results
  list(AUC = AUCs, bias = bias, time = times)
}


########Weibull distribution############

create_pars.weibull <- function(age){
  L <- list(list(pars = c(0.4, 1, 1), age = age, n = n.obs2, type = "weibull"), 
            list(pars = c(0.4, 0.8, 0.5), age = age, n = n.obs2, type = "weibull"), 
            list(pars = c(0.75, 0.5, 0.5), age = age, n = n.obs2, type = "weibull"),
            list(pars = c(0, 1.2, 0.1), age = age, n = n.obs2, type = "weibull"),
            list(pars = c(0, 0.4, 0.1), age = age, n = n.obs2, type = "weibull"),
            list(pars = c(0, 1, 0.05), age = age, n = n.obs2, type = "weibull")
  )
}

#Plot excess hazards
if(FALSE){
  L <- create_pars.weibull(age = 50)
  n.sim2 <- 20
  for(j in 1:length(L)){
    pdf(file.path(fig.out, paste0("ehaz", j, ".pdf")), width = 8, height = 6)
    for(i in 1:n.sim2){
      cat("i = ", i, "\n")
      calc_AUCs(do.call(sim_wei, L[[j]]), plot.haz = T)
    }
    dev.off() 
  }
  
  for(j in 1:length(L)){
    pdf(file.path(fig.out, paste0("esurv", j, ".pdf")), width = 8, height = 6)
    for(i in 1:n.sim2){
      cat("i = ", i, "\n")
      calc_AUCs(do.call(sim_wei, L[[j]]), plot.surv = T)
    }
    dev.off() 
  } 
}

#Simulation parameters
n.obs2 <- 1000
ages <- c(50, 60, 70)
n.sim <- 500
n.cores <- 48


filename <- file.path(data.out, "Simulation_LOL_cases_weibull.RData")
if(file.exists(filename)){
  load(filename)
}else{
  set.seed(20180601, "L'Ecuyer")
  L_all_weibull <- vector("list", length(ages))
  names(L_all_weibull) <- ages
  for(age in ages){
    cat("Age = ", age, "\n")
    L <- create_pars.weibull(age = age)
    L_res <- L_res_time <- vector("list", length(L))
    for(i in 1:length(L)){
      cat(i, "\n")
      M <- mclapply(1:n.sim, function(j){
        #cat(j, "\n")
        df <- do.call(sim_surv, L[[i]])
        # rsfit <- rs.surv(Surv(FU, status) ~ 1 + ratetable(age = age, sex = sex, year = diag_date), 
        #                  ratetable = survexp.dk, data = df$D)
        # rsfit$time <- rsfit$time / ayear
        # plot(rsfit)
        # f <- function(time) L[[i]]$pars[1] + (1 - L[[i]]$pars[1]) * exp(-L[[i]]$pars[3] * time ^L[[i]]$pars[3])
        # curve(f, from = 0, to =15, add = T, col = 2)
        tmp <- calc_AUCs(df)
        if(is.na(tmp$AUC[5])){
          stop("Let's break out!")
        }
        tmp
      }, mc.cores = n.cores)
      L_res[[i]] <- t(sapply(M, function(x) x$AUC))
      L_res_time[[i]] <- lapply(M, function(x) x$bias)
    }
    L_all_weibull[[as.character(age)]] <- list(L_res, L_res_time, age = age)
  }
  save(L_all_weibull, file = filename)
}

n.models <- ncol(L_all_weibull[[1]][[1]][[1]])
n.cases <- length(L_all_weibull[[1]][[1]])
models <- LETTERS[1:n.models]

#Plot the loss of lifetime biases
plot_data <- lapply(1:n.cases, function(x){
  new_list <- lapply(L_all_weibull[[1]][[2]][[x]], function(y){
    D <- data.frame(biases = c(as.matrix(y)), Model = rep(models, each = length(0:15)), 
                    times = rep(0:15, length(models)))
    D
  })
  D <- do.call(rbind, new_list)
  D$Model <- factor(D$Model, levels = models)
  D$times <- factor(D$times)
  D
})
plot_data_all <- do.call(rbind, plot_data)
plot_data_all$case <- rep(paste0("Scenario ", 1:n.cases), each = nrow(plot_data[[1]]))
plot_data_all$case <- factor(plot_data_all$case)


p <- ggplot(data = plot_data_all[plot_data_all$times %in% c(0, 2, 5, 10),], 
            aes(x = times, y = biases, fill = Model)) + geom_boxplot() + 
  ylab("Loss of lifetime bias, D(t)") + xlab("Years survived") + 
  geom_vline(xintercept = 1:5 + 0.5, linetype = "dashed") + 
  theme_bw() + theme(legend.position = "bottom", 
                     legend.text=element_text(size=19), 
                     legend.title=element_text(size=19),
                     axis.title=element_text(size=20),
                     strip.text = element_text(size=15), 
                     axis.text = element_text(size = 16)) + 
  scale_fill_manual(values = cbPalette) + 
  facet_wrap(~case, ncol = 3) + coord_cartesian(ylim = c(-5, 5))


pdf(file.path(fig.out, "SimulationResults_weibull.pdf"), width = 11.5, height = 8.5)
print(p)
dev.off()

format_res <- function(medians, ranges) paste0(sprintf(format, medians), "(", ranges, ")" )

#Tabulate the loss of lifetime biases

L_tmp <- create_pars.weibull(age = 50)
cure_probs <- unlist(lapply(L_tmp, function(x) x$pars[1]))

res <- lapply(L_all_weibull, function(L_res){
  means <- t(sapply(L_res[[1]], colMeans))
  ranges <- t(sapply(L_res[[1]], function(x){
    tmp <- colRanges(x)
    paste0(sprintf(format, tmp[,1]), "-", sprintf(format, tmp[,2]))
  }))
  res_data <- data.frame(Scenario = rep(1:n.cases, each = 1),
                         pi = rep(cure_probs, each = 1),
                         n = sprintf("%.0f", rep(n.obs2, n.cases)))
  M <- matrix(nrow = nrow(means), ncol = ncol(means))
  for(i in 1:ncol(means)){
    M[,i] <- format_res(means[,i], ranges[,i])
  }
  M <- as.data.frame(M)
  names(M) <- paste0("Model ", LETTERS[1:ncol(M)])
  res_data <- cbind(res_data, M)
  res_data <- rbind(NA, res_data)
  age_df <- data.frame(Age = c(NA, as.character(L_res$age), rep(NA, nrow(res_data) - 2)))
  res_data <- cbind(age_df, res_data)
  #res_data <- cbind(c(paste0("Age=", L_res$age), rep(NA, nrow(res_data) - 1)), res_data)
  #names(res_data)[1] <- ""
  res_data
})
all_res <- do.call(rbind, res)
all_res <- subset(all_res[-1,], select = -n)
align <- paste0(rep("c", ncol(all_res) + 1), collapse = "")
if(thesis){
  print(xtable(all_res, caption = paste0("The integrated loss of lifetime bias in the Weibull scenario, 
                                       computed by integrating $|D(t)|$ from 0 to 15 years. 
                                       The loss of lifetime was computed for 50-, 
                                       60-, and 70-year-old patients. The mean and range from the ", n.sim, " 
                                       simulations are provided."), 
               label = "tab:sim_res_weibull", align = align), file = file.path(tab.out, "sim_res_weibull.tex"), 
        include.rownames = F, scalebox = 0.7, table.placement = "h!")#, floating.environment = "sidewaystable")
} else {
  print(xtable(all_res, caption = paste0("The integrated loss of lifetime bias in the Weibull scenario, 
                                       computed by integrating $|D(t)|$ from 0 to 15 years. 
                                       The loss of lifetime was computed for 50-, 
                                       60-, and 70-year-old patients. The mean and range from the ", n.sim, " 
                                       simulations are provided."), 
               label = "tab:sim_res_weibull", align = align), file = file.path(tab.out, "sim_res_weibull.tex"), 
        include.rownames = F)#, floating.environment = "sidewaystable")
}

#Plot the AUCs
# plot_data_age <- lapply(1:length(ages), function(y){
#   plot_data <- lapply(L_all_weibull[[y]][[1]], function(x){
#     D <- data.frame(bias = c(x), model = rep(models, each = n.sim))
#   })
#   plot_data <- do.call(rbind, plot_data)
#   plot_data$case <- rep(paste0("Case ", 1:length(L)), each = n.sim * length(models))
#   plot_data
# })
# 
# plot_data_age <- do.call(rbind, plot_data_age)
# plot_data_age$age <- rep(paste0(ages, " years of age"), each = length(models) * n.sim * length(L))
# 
# pdf(file.path(fig.out, "SimulationintWeibull.pdf"), width = 10, height = 10)
# ggplot(plot_data_age, aes(x = model, y = sqrt(bias), fill = model)) + geom_boxplot() + facet_grid(age~case) + 
#   theme_bw() + theme(legend.title=element_blank(), legend.position = "bottom",
#                      axis.title.x = element_blank(), 
#                      axis.text.x = element_blank(),
#                      axis.ticks.x = element_blank()) + 
#   ylab("Square rooted integrated loss of lifetime bias")
# dev.off()



############Gamma distributed###########
create_pars.gamma <- function(age){
  L <- list(list(pars = c(pi = 0.4, scale = 0.9, d = 0.9, k = 0.9), age = age, n = n.obs2, type = "gengamma"), 
            list(pars = c(pi = 0.4, scale = 3, d = 0.7, k = 0.5), age = age, n = n.obs2, type = "gengamma"),
            list(pars = c(pi = 0.75, scale = 7, d = 1.2, k = 0.9), age = age, n = n.obs2, type = "gengamma"),
            list(pars = c(pi = 0, scale = 7, d = 1.2, k = 0.9), age = age, n = n.obs2, type = "gengamma"),
            list(pars = c(pi = 0, scale = 25, d = 0.4, k = 1.8), age = age, n = n.obs2, type = "gengamma"),
            list(pars = c(pi = 0, scale = 20, d = 1.4, k = 0.7), age = age, n = n.obs2, type = "gengamma")
  )
}

filename <- file.path(data.out, "Simulation_LOL_cases_gamma.RData")
if(file.exists(filename)){
  load(filename)
}else{
  set.seed(20180601, "L'Ecuyer")
  L_all_gamma <- vector("list", length(ages))
  names(L_all_gamma) <- ages
  for(age in ages){
    cat("Age = ", age, "\n")
    L <- create_pars.gamma(age = age)
    L_res <- L_res_time <- vector("list", length(L))
    for(i in 1:length(L)){
      cat(i, "\n")
      M <- mclapply(1:n.sim, function(j){
        #cat(j, "\n")
        df <- do.call(sim_surv, L[[i]])
        # rsfit <- rs.surv(Surv(FU, status) ~ 1 + ratetable(age = age, sex = sex, year = diag_date), 
        #                  ratetable = survexp.dk, data = df$D)
        # rsfit$time <- rsfit$time / ayear
        # plot(rsfit)
        # f <- function(time) L[[i]]$pars[1] + (1 - L[[i]]$pars[1]) * exp(-L[[i]]$pars[3] * time ^L[[i]]$pars[3])
        # curve(f, from = 0, to =15, add = T, col = 2)
        tmp <- calc_AUCs(df)
        if(is.na(tmp$AUC[5])){
          stop("Let's break out!")
        }
        tmp
      }, mc.cores = n.cores)
      L_res[[i]] <- t(sapply(M, function(x) x$AUC))
      L_res_time[[i]] <- lapply(M, function(x) x$bias)
    }
    L_all_gamma[[as.character(age)]] <- list(L_res, L_res_time, age = age)
  }
  save(L_all_gamma, file = filename)
}


# filename <- "GeneratedData/Simulation_LOL_cases_gamma.RData"
# if(file.exists(filename)){
#   load(filename)
# }else{
#   set.seed(20170525)
#   L_all_gamma <- lapply(ages, function(age){
#     cat("Age = ", age, "\n")
#     L <- create_pars(age = age)
#     L_res <- L_res_time <- vector("list", length(L))
#     for(i in 1:length(L)){
#       cat(i, "\n")
#       M <- matrix(nrow = n.sim, ncol = 5)
#       M_time <- vector("list", n.sim)
#       for(j in 1:n.sim){
#         df <- do.call(sim_surv, L[[i]])
#         tmp <- calc_AUCs(df)
#         M[j,] <- tmp$AUC
#         M_time[[j]] <- tmp$bias[tmp$time %in% 0:15,]
#       }
#       L_res[[i]] <- M 
#       L_res_time[[i]] <- M_time
#     }
#     list(L_res, L_res_time, age = age)
#   })
#   save(L_all_gamma, file = filename)
# }

n.models <- ncol(L_all_gamma[[1]][[1]][[1]])
n.cases <- length(L_all_gamma[[1]][[1]])
models <- LETTERS[1:n.models]

#Plot the loss of lifetime biases
plot_data <- lapply(1:n.cases, function(x){
  new_list <- lapply(L_all_gamma[[1]][[2]][[x]], function(y){
    D <- data.frame(biases = c(as.matrix(y)), Model = rep(models, each = length(0:15)), 
                    times = rep(0:15, length(models)))
    D
  })
  D <- do.call(rbind, new_list)
  D$Model <- factor(D$Model, levels = models)
  D$times <- factor(D$times)
  D
})
plot_data_all <- do.call(rbind, plot_data)
plot_data_all$case <- rep(paste0("Scenario ", 1:n.cases), each = nrow(plot_data[[1]]))
plot_data_all$case <- factor(plot_data_all$case)

p <- ggplot(data = plot_data_all[plot_data_all$times %in% c(0, 2, 5, 10),], 
            aes(x = times, y = biases, fill = Model)) + geom_boxplot() + 
  ylab("Loss of lifetime bias, D(t)") + xlab("Years survived") + 
  geom_vline(xintercept = 1:5 + 0.5, linetype = "dashed") + 
  theme_bw() + theme(legend.position = "bottom", 
                     legend.text=element_text(size=19), 
                     legend.title=element_text(size=19),
                     axis.title=element_text(size=20),
                     strip.text = element_text(size=15), 
                     axis.text = element_text(size = 16)) +  
  scale_fill_manual(values = cbPalette) + 
  facet_wrap(~case, ncol = 3) + coord_cartesian(ylim = c(-5, 5))


pdf(file.path(fig.out, "SimulationResults_gamma.pdf"), width = 11.5, height = 8.5)
print(p)
dev.off()

format_res <- function(medians, ranges) paste0(sprintf(format, medians), "(", ranges, ")" )

#Tabulate the loss of lifetime biases
res <- lapply(L_all_gamma, function(L_res){
  means <- t(sapply(L_res[[1]], colMeans, na.rm = T))
  ranges <- t(sapply(L_res[[1]], function(x){
    tmp <- colRanges(x, na.rm = T)
    paste0(sprintf(format, tmp[,1]), "-", sprintf(format, tmp[,2]))
  }))
  res_data <- data.frame(Scenario = rep(1:n.cases, each = 1),
                         pi = rep(cure_probs, each = 1),
                         n = sprintf("%.0f", rep(n.obs2, n.cases)))
  M <- matrix(nrow = nrow(means), ncol = ncol(means))
  for(i in 1:ncol(means)){
    M[,i] <- format_res(means[,i], ranges[,i])
  }
  M <- as.data.frame(M)
  names(M) <- paste0("Model ", LETTERS[1:ncol(M)])
  res_data <- cbind(res_data, M)
  res_data <- rbind(NA, res_data)
  age_df <- data.frame(Age = c(NA, as.character(L_res$age), rep(NA, nrow(res_data) - 2)))
  res_data <- cbind(age_df, res_data)
  #res_data <- cbind(c(paste0("Age=", L_res$age), rep(NA, nrow(res_data) - 1)), res_data)
  #names(res_data)[1] <- ""
  res_data
})
all_res <- do.call(rbind, res)
all_res <- subset(all_res[-1,], select = -n)
align <- paste0(rep("c", ncol(all_res) + 1), collapse = "")
if(thesis){
  print(xtable(all_res, caption = paste0("The integrated loss of lifetime bias in the generalized gamma scenario, 
                                       computed by integrating $|D(t)|$ from 0 to 15 years.
                                       The loss of lifetime was simulated for 50-, 60-, and 70-year-old patients. 
                                       The mean and range from the ", n.sim, " simulations are provided."), 
               label = "tab:sim_res_gamma", align = align), file = file.path(tab.out, "sim_res_gamma.tex"), 
        include.rownames = F, scalebox = 0.7, table.placement = "h!")#, floating.environment = "sidewaystable")
} else {
  print(xtable(all_res, caption = paste0("The integrated loss of lifetime bias in the generalized gamma scenario, 
                                       computed by integrating $|D(t)|$ from 0 to 15 years.
                                         The loss of lifetime was simulated for 50-, 60-, and 70-year-old patients. 
                                         The mean and range from the ", n.sim, " simulations are provided."), 
               label = "tab:sim_res_gamma", align = align), file = file.path(tab.out, "sim_res_gamma.tex"), 
        include.rownames = F)#, floating.environment = "sidewaystable") 
}

#Plot the AUCs
# plot_data_age <- lapply(1:length(ages), function(y){
#   plot_data <- lapply(L_all_gamma[[y]][[1]], function(x){
#     D <- data.frame(bias = c(x), model = rep(models, each = n.sim))
#   })
#   plot_data <- do.call(rbind, plot_data)
#   plot_data$case <- rep(paste0("Case ", 1:length(L)), each = n.sim * length(models))
#   plot_data
# })
# 
# plot_data_age <- do.call(rbind, plot_data_age)
# plot_data_age$age <- rep(paste0(ages, " years of age"), each = length(models) * n.sim * length(L))
# 
# pdf(file.path(fig.out, "Simulationintgamma.pdf"), width = 10, height = 10)
# ggplot(plot_data_age, aes(x = model, y = sqrt(bias), fill = model)) + geom_boxplot() + facet_grid(age~case) + 
#   theme_bw() + theme(legend.title=element_blank(), legend.position = "bottom",
#                      axis.title.x = element_blank(), 
#                      axis.text.x = element_blank(),
#                      axis.ticks.x = element_blank()) + 
#   ylab("Square rooted integrated loss of lifetime bias")
# dev.off()

#Plot scenarios for Weibull and Gamma together
time <- seq(0, 30, length.out = 200)

L_wei <- list(c(0.4, 1, 1), 
              c(0.4, 0.8, 0.5),
              c(0.75, 0.5, 0.5),
              c(0, 1.2, 0.1),
              c(0, 0.4, 0.1),
              c(0, 1, 0.05))

L_gam <- list(c(pi = 0.4, scale = 0.9, d = 0.9, k = 0.9), 
              c(pi = 0.4, scale = 3, d = 0.7, k = 0.5),
              c(pi = 0.75, scale = 7, d = 1.2, k = 0.9),
              c(pi = 0, scale = 7, d = 1.2, k = 0.9),
              c(pi = 0, scale = 25, d = 0.4, k = 1.8),
              c(pi = 0, scale = 20, d = 1.4, k = 0.7))

rs_wei <- function(t, pars) pars[1] + (1 - pars[1]) * exp(-pars[3] * t ^ pars[2])
rs_gam <- function(t, pars) pars[1] + (1 - pars[1]) * pgengamma.stacy(time, scale = pars[2], 
                                                                      d = pars[3], k = pars[4], lower.tail = F)
res_wei <- lapply(L_wei, function(pars) data.frame(RS = rs_wei(time, pars), Time = time))
res_gam <- lapply(L_gam, function(pars) data.frame(RS = rs_gam(time, pars), Time = time))

res1 <- do.call(rbind, res_wei)
res1$scenario <- rep(rep(1:(length(L_wei) / 2), each = length(time)), length(L_wei) / 3)
res1$model <- rep(c("Mixture cure model", "Relative survival model"), each = length(time) * 3)

res2 <- do.call(rbind, res_gam)
res2$scenario <- rep(rep(1:(length(L_gam) / 2), each = length(time)), length(L_gam) / 3)
res2$model <- rep(c("Mixture cure model", "Relative survival model"), each = length(time) * 3)


res <- rbind(res1, res2)
dists <- c("Weibull", "Generalized gamma")
res$dist <- factor(rep(dists, c(nrow(res1), nrow(res2))), levels = dists)
res$Scenario <- factor(res$scenario)
res$model <- factor(res$model)



rs_wei <- function(t, pars) pars[1] + (1 - pars[1]) * exp(-pars[3] * t ^ pars[2])
rs_gam <- function(t, pars) pars[1] + (1 - pars[1]) * pgengamma.stacy(time, scale = pars[2], 
                                                                      d = pars[3], k = pars[4], lower.tail = F)
res_wei <- lapply(L_wei, function(pars) data.frame(RS = rs_wei(time, pars), Time = time))
res_gam <- lapply(L_gam, function(pars) data.frame(RS = rs_gam(time, pars), Time = time))

res1 <- do.call(rbind, res_wei)
res1$scenario <- rep(1:length(L_wei), each = length(time))

res2 <- do.call(rbind, res_gam)
res2$scenario <- rep(1:length(L_gam), each = length(time))


res <- rbind(res1, res2)
dists <- c("Weibull", "Generalized gamma")
res$dist <- factor(rep(dists, c(nrow(res1), nrow(res2))), levels = dists)
res$Scenario <- factor(res$scenario)

pdf(file.path(fig.out, "sim_models.pdf"), width = 10, height = 6)
ggplot(res, aes(x = Time, y = RS, colour = Scenario, linetype = Scenario)) + 
  geom_line() + facet_grid(~dist) +  
  scale_linetype_manual(values = c("solid", "solid", "dashed", 
                                   "dashed", "twodash", "twodash")) + 
  scale_color_manual(values = rep(c("black", "grey"), 3)) + 
  theme_bw() + ylab("Relative survival") + xlab("Time (years)") + 
  ylim(c(0, 1)) + geom_vline(xintercept = 15, linetype = "dotted") + 
  theme(legend.position = "bottom", 
        legend.text=element_text(size=15),
        legend.title=element_text(size=15),
        axis.title=element_text(size=17),
        axis.text = element_text(size = 13), 
        strip.text = element_text(size = 15), 
        legend.key.size = unit(2,"line")) + 
  guides(colour = guide_legend(ncol = 6))
dev.off()



ggplot(res, aes(x = Time, y = RS, colour = Scenario, linetype = Scenario)) + 
  geom_line() + facet_grid(~dist) +  
  scale_linetype_manual(values = c("solid", "solid", "dashed", 
                                   "dashed", "twodash", "twodash")) + 
  scale_color_manual(values = rep(c("black", "grey"), 3)) + 
  theme_bw() + ylab("Relative survival") + xlab("Time (years)") + 
  ylim(c(0, 1)) + geom_vline(xintercept = 15, linetype = "dotted") + 
  theme(legend.position = "bottom", 
        legend.text=element_text(size=15),
        legend.title=element_text(size=15),
        axis.title=element_text(size=17),
        axis.text = element_text(size = 13), 
        strip.text = element_text(size = 15), 
        legend.key.size = unit(2,"line")) + 
  guides(colour = guide_legend(ncol = 6))
