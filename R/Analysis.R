################# Lymphoma registry data ####################

get.knots <- function(data, k){
  quantile(data$FU_years[data$status ==1], seq(0, 1, length.out = k))
}

##Fit parametric relative survival models
#DLBCL
knots <- log(get.knots(DLBCL, 5))
fit_DLBCL <- GenFlexCureModel(Surv(FU_years, status) ~ -1, data = DLBCL, bhazard = "exp_haz", 
                              smooth.formula = ~cb(x = log(FU_years), knots = knots))

knots <- log(get.knots(DLBCL, 6))
fit_DLBCL2 <- stpm2(Surv(FU_years, status) ~ -1, data = DLBCL, bhazard = DLBCL$exp_haz, 
                    smooth.formula = ~cb(x = log(FU_years), knots = knots))

knots <- log(sort(c(get.knots(DLBCL, 6), 10)))
fit_DLBCL3 <- stpm2(Surv(FU_years, status) ~ -1, data = DLBCL, bhazard = DLBCL$exp_haz, 
                    smooth.formula = ~cbc(x = log(FU_years), knots = knots))

#FL
knots <- log(get.knots(FL, 5))
fit_FL <- GenFlexCureModel(Surv(FU_years, status) ~ -1, data = FL, bhazard = "exp_haz", 
                              smooth.formula = ~cb(x = log(FU_years), knots = knots))

knots <- log(get.knots(FL, 6))
fit_FL2 <- stpm2(Surv(FU_years, status) ~ -1, data = FL, bhazard = FL$exp_haz, 
                    smooth.formula = ~cb(x = log(FU_years), knots = knots))

knots <- log(sort(c(get.knots(FL, 6), 10)))
fit_FL3 <- stpm2(Surv(FU_years, status) ~ -1, data = FL, bhazard = FL$exp_haz, 
                    smooth.formula = ~cbc(x = log(FU_years), knots = knots))


#ML
knots <- log(get.knots(ML, 5))
fit_ML <- GenFlexCureModel(Surv(FU_years, status) ~ -1, data = ML, bhazard = "exp_haz", 
                              smooth.formula = ~cb(x = log(FU_years), knots = knots))

knots <- log(get.knots(ML, 6))
fit_ML2 <- stpm2(Surv(FU_years, status) ~ -1, data = ML, bhazard = ML$exp_haz, 
                    smooth.formula = ~cb(x = log(FU_years), knots = knots))

knots <- log(sort(c(get.knots(ML, 6), 10)))
fit_ML3 <- stpm2(Surv(FU_years, status) ~ -1, data = ML, bhazard = ML$exp_haz, 
                    smooth.formula = ~cbc(x = log(FU_years), knots = knots))


# Old implementation of the relative survival models
# fit_ML <- FlexCureModel(Surv(FU_years, status) ~ 1, data = ML, bhazard = "exp_haz",
#                                n.knots = 5)
# fit_ML2 <- stpm2(Surv(FU_years, status) ~ 1, data = ML, bhazard = ML$exp_haz, df = 5)
# knots_ML <- log(sort(c(quantile(ML$FU_years[ML$status ==1], c(0, 1, 0.2, 0.4, 0.6, 0.8)), 10)))
# fit_ML3 <- stpm2(Surv(FU_years, status) ~ -1, data = ML, bhazard = ML$exp_haz, 
#                     smooth.formula = ~basis_cure(knots = knots_ML, x = log(FU_years)))

#Assemble parametric models
fits <- list(DLBCL = list(fit_DLBCL, fit_DLBCL2, fit_DLBCL3, data = DLBCL),
             FL = list(fit_FL, fit_FL2, fit_FL3, data = FL), 
             ML = list(fit_ML, fit_ML2, fit_ML3, data = ML))

#Calculate non-parametric and parametric relative survival
plot_data <- lapply(fits, function(fit){
  rsfit <- rs.surv(Surv(FU, status) ~ 1 + ratetable(age = age, sex = sex, year = diag_date),
                   data = fit$data, ratetable = survexp.dk, method = "ederer2")
  rsfit$time <- rsfit$time / ayear
  D <- data.frame(RS = rsfit$surv, time = rsfit$time, ci.lower = rsfit$lower, ci.upper = rsfit$upper) 
  
  pred1 <- predict(fit[[1]], time = D$time, var.type = "n")[[1]]
  pred2 <- data.frame(Estimate = predict(fit[[2]], newdata = data.frame(FU_years = D$time)))
  pred3 <- data.frame(Estimate = predict(fit[[3]], newdata = data.frame(FU_years = D$time)))
  
  D_para <- rbind(pred1, pred2, pred3)
  
  D_para$time <- rep(D$time, 3)
  D_para$model <- rep(c("FMC", "NRS", "ARS"), each = nrow(D))
  list(D = D, D_para = D_para)
})

#Plot relative survival
para_plot_data <- do.call(rbind, lapply(plot_data, function(x) x$D_para))
para_plot_data$disease <- rep(c("DLBCL", "FL", "ML"), sapply(plot_data, function(x) nrow(x$D_para)))

npara_plot_data <- do.call(rbind, lapply(plot_data, function(x) x$D))
npara_plot_data$disease <- rep(c("DLBCL", "FL", "ML"), sapply(plot_data, function(x) nrow(x$D)))
colnames(npara_plot_data)[1] <- "Estimate"
npara_plot_data$model <- "Ederer II estimate"

p <- ggplot(data = npara_plot_data, aes(x = time, y = Estimate, group = model, colour = model)) + geom_step() + 
  facet_grid(.~disease) + geom_step(data = npara_plot_data, aes(x = time, y = ci.lower), linetype = "dashed") + 
  geom_step(data = npara_plot_data, aes(x = time, y = ci.upper), linetype = "dashed") + 
  geom_line(data = para_plot_data, aes(x = time, y = Estimate), size = 1) + ylim(c(0, 1.002)) + 
  scale_colour_manual(values = c("Ederer II estimate" = "black", "FMC" = "brown2", 
                                 "NRS" = "darkolivegreen3", "ARS" = "deepskyblue3"), 
                      breaks = c("Ederer II estimate", "NRS", 
                                 "ARS", "FMC")) + 
  theme_bw() + theme(legend.position = "bottom", 
                     legend.title = element_blank(), 
                     legend.text=element_text(size=15), 
                     axis.title=element_text(size=17),
                     strip.text = element_text(size=15), 
                     axis.text = element_text(size = 13)) + 
  xlab("Follow-up time (years)") + 
  ylab("Relative survival")

pdf(file.path(fig.out, "RSCombined.pdf"), width = 9.5, height = 5)
print(p)
dev.off()

png(file.path(fig.out, "RSCombined.png"), res = 200, width = 2300, height = 2300 * 5 / 11)
print(p)
dev.off()

#Function for calculating mean age
get_ages <- function(data){
  bdr <- floor(range(data$age_years))
  mean_age <- floor(median(data$age_years))
  paste0(mean_age, "(", bdr[1], "-", bdr[2], ")")
}

#Create table with age, relative survival, and loss of lifetime estimates
M <- matrix(nrow = 7, ncol = 5)
M[,1] <- c("Median age (range)", "5-year RS (95% CI)", NA, NA, "Loss of lifetime (95% CI)", NA, NA)
M[,2] <- c(NA, rep(c("NRS", "ARS", "FMC"), 2))

M[1, 2:5] <- c(NA, get_ages(DLBCL), get_ages(FL), get_ages(ML))

predict(fit_DLBCL, time = 5)
predict(fit_DLBCL2, newdata = data.frame(FU_years = 5))

get_rs <- function(fit, time = 5){
  if("cuRe" %in% class(fit)){
    pred <- round(predict(fit, time = 5)[[1]], 2) 
    paste0(pred$Estimate, "(", pred$lower, "-", pred$upper, ")")
  }else{
    pred <- round(predict(fit, newdata = data.frame(FU_years = 5), se.fit = T), 2)
    paste0(pred$Estimate, "(", pred$lower, "-", pred$upper, ")")
  }
}

M[2, 3:5] <- c(get_rs(fit_DLBCL2), get_rs(fit_FL2), get_rs(fit_ML2))
M[3, 3:5] <- c(get_rs(fit_DLBCL3), get_rs(fit_FL3), get_rs(fit_ML3))
M[4, 3:5] <- c(get_rs(fit_DLBCL), get_rs(fit_FL), get_rs(fit_ML))

#Function for calculating loss of lifetime estimates
get_LL <- function(fit){
  LL_res <- calc.LL(fit, time = 0, rmap = list(year = diag_date), smooth.exp = FALSE)
  LL_res <- sprintf("%.2f", LL_res[[1]])
  paste0(LL_res[1], "(", LL_res[2], "-", LL_res[3], ")")
}

M[5, 3:5] <- c(get_LL(fit_DLBCL2), get_LL(fit_FL2), get_LL(fit_ML2))
M[6, 3:5] <- c(get_LL(fit_DLBCL3), get_LL(fit_FL3), get_LL(fit_ML3))
M[7, 3:5] <- c(get_LL(fit_DLBCL), get_LL(fit_FL), get_LL(fit_ML))

M <- as.data.frame(M)
colnames(M) <- c("", "Model","DLBCL", "FL", "ML")

if(thesis){
  print(xtable(M, caption = "Median age, 5-year relative survival (RS), and loss of lifetime estimates at time zero 
             in Danish diffuse large B-cell lymphoma (DLBCL), follicular lymphoma (FL), 
               and mantle cell lymphoma (ML) patients.",
               label = "tab:sum", align = "cllccc"), include.rownames = F,
        file = file.path(tab.out, "SummaryMeasureTable.tex"), scalebox = 0.85)
}else{
  print(xtable(M, caption = "Median age, 5-year relative survival (RS), and loss of lifetime estimates at time zero 
             in Danish diffuse large B-cell lymphoma (DLBCL), follicular lymphoma (FL), 
             and mantle cell lymphoma (ML) patients.",
               label = "tab:sum", align = "cllccc"), include.rownames = F,
        file = file.path(tab.out, "SummaryMeasureTable.tex"))
}


#Choose time points for the loss of lifetime estimates
times <- seq(0, 10, length.out = 50)

#Calculate loss of lifetime for the three diseases using the three models
#DLBCL
res_DLBCL <- calc.LL(fit_DLBCL, time = times, rmap = list(year = diag_date), smooth.exp = F)[[1]]
res_DLBCL2 <- calc.LL(fit_DLBCL2, time = times, rmap = list(year = diag_date), smooth.exp = F)[[1]]
res_DLBCL3 <- calc.LL(fit_DLBCL3, time = times, rmap = list(year = diag_date), smooth.exp = F)[[1]]
#FL
res_FL <- calc.LL(fit_FL, time = times, rmap = list(year = diag_date), smooth.exp = F)[[1]]
res_FL2 <- calc.LL(fit_FL2, time = times, rmap = list(year = diag_date), smooth.exp = F)[[1]]
res_FL3 <- calc.LL(fit_FL3, time = times, rmap = list(year = diag_date), smooth.exp = F)[[1]]
#ML
res_ML <- calc.LL(fit_ML, time = times, rmap = list(year = diag_date), smooth.exp = F)[[1]]
res_ML2 <- calc.LL(fit_ML2, time = times, rmap = list(year = diag_date), smooth.exp = F)[[1]]
res_ML3 <- calc.LL(fit_ML3, time = times, rmap = list(year = diag_date), smooth.exp = F)[[1]]

#Combine results into a single data frame for plotting
res_all <- rbind(res_DLBCL, res_FL, res_ML, res_DLBCL2, res_FL2, res_ML2, res_DLBCL3, res_FL3, res_ML3)
res_all$disease <- rep(rep(c("DLBCL", "FL", "ML"), each = length(times)), 3)
res_all$Time <- rep(times, 9)
levs <- c("FMC", "NRS", "ARS")
res_all$model <- factor(rep(levs, each = length(times) * 3), levels = levs[c(2, 3, 1)])

#Plot the loss of lifetime curves
pdf(file.path(fig.out, "LOLLymphoma.pdf"), width = 10, height = 5.3)
ggplot(res_all, aes(x = Time, y = Estimate, group = model, linetype = model)) + geom_line(size = 1) +
  facet_grid(.~disease) + xlab("Follow-up time (years)") + ylab("Loss of lifetime (years)") +
  scale_x_continuous(breaks = seq(0, 12, by = 3)) +
  geom_hline(yintercept = 0, linetype = "dashed") + theme_bw() +
  theme(legend.position = "bottom", legend.title = element_blank(), 
        legend.key=element_rect(fill=NA), 
        legend.text=element_text(size=15),
        axis.title=element_text(size=17),
        strip.text = element_text(size=15), 
        axis.text = element_text(size = 14), 
        legend.key.size = unit(2,"line"))# +
  #scale_colour_manual(values = c("FMC" = "brown2", 
  #                               "NRS" = "darkolivegreen3", 
  #                               "ARS" = "deepskyblue3"))
dev.off()


################# Cancer registry data ####################
#The analyses are restricted to mean residual lifetime since the expected 
#survival is the same regardless of length of follow-up
#Compute the mean residual lifetime estimates using KM

gaussxw <- statmod::gauss.quad(100)
#The function allows for a flexible parametric version also and computes the integral
#using the rectangle rule
calc_LOLKM <- function(sfit, expected, time, tau, type = "KM"){
  if(type == "KM"){
    surv_fun <- function(t){
      s <- summary(sfit, t)
      names(s$surv) <- s$time
      a <- s$surv[as.character(t)]
      names(a) <- NULL
      a
    } 
  }else{
    surv_fun <- function(t) exp(-exp(basis(knots = sfit$knots, x = log(t)) %*% sfit$coefficients))
  }
  
  exp_fun <- function(t){
    s <- summary(expected, t)
    names(s$surv) <- s$time
    a <- s$surv[as.character(t)]
    names(a) <- NULL
    a
  }
  
  scale <- (tau - time) / 2
  scale2 <- (tau + time) / 2
  eval_gen_t <- exp_fun(time)
  eval_pop_t <- surv_fun(time)
  eval <- rep(NA, length(time))
  for(i in 1:length(time)){
    points <- scale[i] * gaussxw$nodes + scale2[i]
    eval_gen <- exp_fun(points)
    eval_pop <- surv_fun(points)
    inner_int <- eval_gen / eval_gen_t[i] - eval_pop / eval_pop_t[i]
    eval[i] <- sum(gaussxw$weights * inner_int)
  }
  scale * eval
  
  # t_new <- sort(unique(c(time, seq(0, tau, length.out = 5000))), decreasing = T)
  # df_time <- -diff(t_new)
  # mid_points <- t_new[-length(t_new)] + diff(t_new) / 2
  # vals_pop <- c(0, cumsum(surv_fun(mid_points) * df_time))
  # vals_pop <- rev(vals_pop[t_new %in% time])
  # vals_exp <- c(0, cumsum(exp_fun(mid_points) * df_time))
  # vals_exp <- rev(vals_exp[t_new %in% time])
  # vals_exp / exp_fun(time) - vals_pop / surv_fun(time)
}

#Set clinical subgroups and relevant time points
diseases <- levels(CR_tumor$disease)
ages <- levels(CR_tumor$age_group)
time <- seq(0, 15, length.out = 100)

#Compute the loss of lifetime function
L <- lapply(diseases, function(disease){
  LOLs <- lapply(ages, function(age_group){
    data_new <- CR_tumor[CR_tumor$disease == disease & CR_tumor$age_group == age_group, ]
    sfit <- survfit(Surv(FU, status) ~ 1, data = data_new)  
    #sfit <- flexsurvspline(Surv(FU, status) ~ 1, data = data_new, k = 6)
    expected <- survexp( ~ 1, rmap = list(age = age, sex = sex, year = diag_date),
                       ratetable = survexp.dk, scale = ayear, data = data_new,
                       time = seq(0, 70, length.out = 2000) * ayear)
    tau <- max(data_new$FU)
    if(tau > 40) tau <- 40
    #print(c(tau, summary(sfit, tau)$surv))
    calc_LOLKM(sfit, expected, time = time, tau = tau)
  })
  names(LOLs) <- ages
  LOLs
})
names(L) <- diseases

#Create new dataset with limited follow-up period
#Patients are censored at 16 years, i.e., the a new follow-up and status variable .
CR_tumor2 <- CR_tumor
CR_tumor2$last_followup <- CR_tumor2$D_STATDATO
CR_tumor2$last_followup[CR_tumor2$D_STATDATO > "1976-01-01"] <- "1976-01-01"
CR_tumor2$FU_days <- as.numeric(CR_tumor2$last_followup - CR_tumor2$diag_date)
CR_tumor2$FU_years <- CR_tumor2$FU_days / ayear
CR_tumor2$status[CR_tumor2$D_STATDATO > "1976-01-01"] <- 0
CR_tumor2$exp_haz <- general.haz(time = "FU_days", age = "age", sex = "sex", year = "diag_date", 
                                 data = CR_tumor2, ratetable = survexp.dk)

#Fit models for all diseases and age groups and compute mean residual lifetime
L_fit <- lapply(diseases, function(disease){
  cat(disease, "\n")
  LOLs <- lapply(ages, function(age_group){
    cat(age_group, "\n")
    data_new <- CR_tumor2[CR_tumor2$disease == disease & CR_tumor2$age_group == age_group, ]
    #Fit model by Nelson et al. 2007
    knots <- log(get.knots(data_new, 6))
    fit_nelson <- stpm2(Surv(FU_years, status) ~ -1, data = data_new, bhazard = data_new$exp_haz, 
                        smooth.formula = ~cb(x = log(FU_years), knots = knots))
        
    #Fit models by Andersson et al. 2011
    add.knot <- 10
    knots <- sort(c(knots, log(add.knot)))
    #knots_andersson1 <- log(sort(c(quantile(data_new$FU_years[data_new$status ==1], 
    #                                        c(0, 0.2, 0.4, 0.6, 0.8, 1)), add.knot)))
    fit_andersson1 <- stpm2(Surv(FU_years, status) ~ -1, data = data_new, bhazard = data_new$exp_haz, 
                            smooth.formula = ~cbc(x = log(FU_years), knots = knots))
    
    last.knot <- 80
    knots[length(knots)] <- log(last.knot)
    #knots_andersson2 <- log(sort(c(quantile(data_new$FU_years[data_new$status ==1], 
    #                                        c(0, 0.2, 0.4, 0.6, 0.8)), add.knots)))
    fit_andersson2 <- stpm2(Surv(FU_years, status) ~ -1, data = data_new, bhazard = data_new$exp_haz, 
                            smooth.formula = ~cbc(x = log(FU_years), knots = knots))
    
    #Fit flexible mixture cure models
    # fit_flex_mix1 <- FlexCureModel(Surv(FU_years, status) ~ 1, data = data_new, 
    #                                       bhazard = "exp_haz", n.knots = 5, 
    #                                       covariance = F, verbose = F)
    # fit_flex_mix2 <- FlexCureModel(Surv(FU_years, status) ~ 1, data = data_new, 
    #                                       bhazard = "exp_haz", 
    #                                       knots = log(c(min(data_new$FU_years), 0.5, 1, 2, 5)), 
    #                                       covariance = F, verbose = F)
    knots <- log(get.knots(data_new, 5))
    fit_flex_mix1 <- GenFlexCureModel(Surv(FU_years, status) ~ -1, data = data_new, bhazard = "exp_haz", 
                                      smooth.formula = ~cb(x = log(FU_years), knots = knots), 
                                      verbose = F, covariance = F)
    min.time <- min(data_new$FU_years[data_new$status == 1])
    knots <- log(c(min.time, 0.5, 1, 2, 5))
    fit_flex_mix2 <- GenFlexCureModel(Surv(FU_years, status) ~ -1, data = data_new, bhazard = "exp_haz", 
                                      smooth.formula = ~cb(x = log(FU_years), knots = knots), 
                                      verbose = F, covariance = F)
    
    #Plot the models
    # plot(fit_nelson, newdata = data.frame(age = 50), ylim = c(0, 1), ci = F, rug = F)
    # plot(fit_andersson1, newdata = data.frame(age = 50), ylim = c(0, 1), ci = F, rug = F, add = T, line.col = 2)
    # plot(fit_andersson2, newdata = data.frame(age = 50), ylim = c(0, 1), ci = F, rug = F, add = T, line.col = 3)
    # plot(fit_flex_mix1, time = seq(0, 15, length.out = 100), add = T, col = 4, ci = F)
    # plot(fit_flex_mix2, time = seq(0, 15, length.out = 100), add = T, col = 5, ci = F)
    
    tau <- 40
    #Calculate loss of lifetime esimates
    LOL1 <- calc.LL(fit_nelson, time = time, var.type = "n", tau = tau, rmap = list(year = diag_date), smooth.exp = F)[[1]]
    LOL2 <- calc.LL(fit_andersson1, time = time, var.type = "n", tau = tau, rmap = list(year = diag_date), smooth.exp = F)[[1]]
    LOL3 <- calc.LL(fit_andersson2, time = time, var.type = "n", tau = tau, rmap = list(year = diag_date), smooth.exp = F)[[1]]
    LOL4 <- calc.LL(fit_flex_mix1, time = time, var.type = "n", tau = tau, rmap = list(year = diag_date), smooth.exp = F)[[1]]
    LOL5 <- calc.LL(fit_flex_mix2, time = time, var.type = "n", tau = tau, rmap = list(year = diag_date), smooth.exp = F)[[1]]
    cbind(LOL1, LOL2, LOL3, LOL4, LOL5)
    })
  names(LOLs) <- ages
  LOLs
})
names(L_fit) <- diseases

#Assemble mean residual lifetime estimates and calculate biases
models <- LETTERS[1:5]
plot_data <- lapply(diseases, function(disease){
  biases <- lapply(ages, function(age_group){
    bias_matrix <- L_fit[[disease]][[age_group]] - L[[disease]][[age_group]]
    data.frame(bias = unlist(bias_matrix), time = rep(time, length(models)),
               Model = rep(models, each = nrow(bias_matrix)), age_group = age_group)
  })
  biases <- do.call(rbind, biases)
  biases$disease <- disease
  biases
})
#Assemble data in a data frame for plotting
plot_data <- do.call(rbind, plot_data)

#Plot LL biases
p <- ggplot(plot_data, aes(x = time, y = bias, colour = Model, group = Model)) + geom_line() + 
  facet_grid(age_group ~ disease) + geom_hline(yintercept = 0, linetype = "dashed") + 
  ylab("Loss of lifetime bias, D(t)") + xlab("Years survived") + theme_bw() + 
  coord_cartesian(ylim = c(-3, 2)) + 
  scale_color_manual(values = cbPalette) + 
  theme(legend.position = "bottom", 
        legend.text=element_text(size=15),
        legend.title=element_text(size=15),
        axis.title=element_text(size=17),
        strip.text = element_text(size=15), 
        axis.text = element_text(size = 13))

pdf(file.path(fig.out, "LOLBiasOldData.pdf"), width = 10.5, height = 12.5)
print(p)
dev.off()


#Create plots for PhD-defence

for(age in ages){
  p <- ggplot(plot_data[plot_data$age_group == age,], aes(x = time, y = bias, colour = Model, group = Model)) +
    geom_line() +
    facet_grid(age_group ~ disease) + geom_hline(yintercept = 0, linetype = "dashed") +
    ylab("Loss of lifetime bias, D(t)") + xlab("Years survived") + theme_bw() +
    coord_cartesian(ylim = c(-3, 2)) +
    scale_color_manual(values = c(cbPalette[1:4], "black")) +
    theme(legend.position = "bottom",
          legend.text=element_text(size=15),
          legend.title=element_text(size=15),
          axis.title=element_text(size=17),
          strip.text = element_text(size=15),
          axis.text = element_text(size = 13))

  file.out <- "C:/Users/sw1y/Dropbox/PhDDefence/"
  png(file.path(file.out, paste0("LOLBiasOldData_", age, ".png")), res = 200, width = 2800, height = 1000)
  print(p)
  dev.off()
}

old_pars <- par("mai")
png(file.path(file.out, "SurvOldData.png"), res = 200, width = 1600, height = 1300)
par(mfrow = c(2,2), mai = c(0.7, 0.7, 0.3, 0.1))
for(disease in diseases){
  sfit <- survfit(Surv(FU, status) ~ 1, data = CR_tumor[CR_tumor$disease == disease,])
  file.out <- "C:/Users/sw1y/Dropbox/PhDDefence/"
  plot(sfit, xlab = "Years since diagnosis", ylab = "Survival probability")
  title(main = disease, line = 0.5)
}
dev.off()


#Compute the LL estimates and assemble for plotting
models <- LETTERS[1:5]
plot_data <- lapply(diseases, function(disease){
  biases <- lapply(ages, function(age_group){
    all_data <- cbind(L_fit[[disease]][[age_group]], L[[disease]][[age_group]])
    data.frame(LOL = unlist(all_data), time = rep(time, length(models) + 1),
               Model = rep(c(models, "True LOL"), each = nrow(all_data)), 
               age_group = age_group)
  })
  biases <- do.call(rbind, biases)
  biases$disease <- disease
  biases
})
plot_data <- do.call(rbind, plot_data)

#Plot LL estimates
pdf(file.path(fig.out, "LOLOldData.pdf"), width = 12, height = 12)
ggplot(plot_data, aes(x = time, y = LOL, colour = Model, group = Model)) + geom_line() + 
  facet_grid(age_group ~ disease) + geom_hline(yintercept = 0, linetype = "dashed") + 
  ylab("Bias") + xlab("Follow-up time (years)") + theme_bw() + 
  scale_color_manual(values = cbPalette) + 
  theme(legend.position = "bottom", 
        legend.text=element_text(size=15), 
        axis.title=element_text(size=17),
        strip.text = element_text(size=15), 
        axis.text = element_text(size = 13), legend.title = element_blank())
dev.off()

#Calculate predicted survival curves all models in each age group of each disease
for(disease in diseases){
  cat(disease, "\n")
  plot_data <- lapply(ages, function(age_group){
    cat(age_group, "\n")
    #Full follow-up Kaplan-Meier
    data_new <- CR_tumor[CR_tumor$disease == disease & CR_tumor$age_group == age_group, ]
    sfit <- survfit(Surv(FU, status)~ 1, data = data_new)

    #Time points
    time <- seq(0, max(data_new$FU), length.out = 1000)
    
    #Restricted follow-up dataset
    data_new <- CR_tumor2[CR_tumor2$disease == disease & CR_tumor2$age_group == age_group, ]
    
    #Fit model by Nelson et al. 2007
    knots <- log(get.knots(data_new, 6))
    fit_nelson <- stpm2(Surv(FU_years, status) ~ -1, data = data_new, bhazard = data_new$exp_haz, 
                        smooth.formula = ~cb(x = log(FU_years), knots = knots))
    
    #Fit models by Andersson et al. 2011
    add.knot <- 10
    knots <- sort(c(knots, log(add.knot)))
    fit_andersson1 <- stpm2(Surv(FU_years, status) ~ -1, data = data_new, bhazard = data_new$exp_haz, 
                            smooth.formula = ~cbc(x = log(FU_years), knots = knots))
    
    last.knot <- 80
    knots[length(knots)] <- log(last.knot)
    fit_andersson2 <- stpm2(Surv(FU_years, status) ~ -1, data = data_new, bhazard = data_new$exp_haz, 
                            smooth.formula = ~cbc(x = log(FU_years), knots = knots))
    
    #Fit flexible mixture cure models
    knots <- log(get.knots(data_new, 5))
    fit_flex_mix1 <- GenFlexCureModel(Surv(FU_years, status) ~ -1, data = data_new, bhazard = "exp_haz", 
                                      smooth.formula = ~cb(x = log(FU_years), knots = knots), 
                                      verbose = F, covariance = F)
    min.time <- min(data_new$FU_years[data_new$status == 1])
    knots <- log(c(min.time, 0.5, 1, 2, 5))
    fit_flex_mix2 <- GenFlexCureModel(Surv(FU_years, status) ~ -1, data = data_new, bhazard = "exp_haz", 
                                      smooth.formula = ~cb(x = log(FU_years), knots = knots), 
                                      verbose = F, covariance = F)
    
    #Expected survival
    expected <- survexp( ~ 1, rmap = list(age = age, sex = sex, year = diag_date),
                         ratetable = survexp.dk, scale = ayear, data = data_new,
                         times = seq(0, 70, length.out = 2000) * ayear)
    
    #Predict survival probabilities
    model1 <- model2 <- model3 <- model4 <- model5 <- rep(1, length(time))
    model1[-1] <- predict(fit_nelson, newdata = data.frame(FU_years = time[time != 0]))
    model2[-1] <- predict(fit_andersson1, newdata = data.frame(FU_years = time[time != 0]))
    model3[-1] <- predict(fit_andersson2, newdata = data.frame(FU_years = time[time != 0]))
    model4[-1] <- predict(fit_flex_mix1, time = time[time != 0], var.type = "n")[[1]]$Estimate
    model5[-1] <- predict(fit_flex_mix2, time = time[time != 0], var.type = "n")[[1]]$Estimate

    #Merge to data frame
    models <- c(model1, model2, model3, model4, model5) * rep(summary(expected, time)$surv, 5)
    D <- data.frame(Est = models, time = rep(time, 5), Model = rep(LETTERS[1:5], each = length(time)))
    D <- rbind(D, data.frame(Est = sfit$surv, time = sfit$time, Model = "Full KM"))
    #D$Est[D$time == 0] <- 1
    D$age_group <- age_group
    D
  })
  
  #Assemble data and plot survival curves
  plot_data <- do.call(rbind, plot_data)
  p <- ggplot(plot_data, aes(x = time, y = Est, colour = Model, group = Model)) + geom_line() + 
    facet_wrap( ~ age_group, ncol = 2) + ylab("Survival probability") + 
    xlab("Follow-up time (years)") + ggtitle(disease) + theme_bw() + 
    scale_color_manual(values = cbPalette) + scale_x_continuous(breaks = seq(0, 50, by = 10)) + 
    theme(legend.position = "bottom", plot.title = element_text(hjust = 0.5, size = 18), 
          legend.text=element_text(size=15), 
          legend.title=element_text(size=15),
          axis.title=element_text(size=17),
          strip.text = element_text(size=15), 
          axis.text = element_text(size = 13)) 
  
  pdf(file.path(fig.out, paste0("SC_", gsub(" ", "_", disease), ".pdf")), width = 10, height = 7)
  print(p)
  dev.off()
  
  png(file.path(fig.out, paste0("SC_", gsub(" ", "_", disease), ".png")), 
      res = 200, width = 2000, height = 2000 * 6 / 9)
  print(p)
  dev.off()
}


################# Lymphoma registry revisited ####################
#Flexible mixture cure model

probs <- c(0, 1/3, 2/3, 1)
# 
# knots_DLBCL <- quantile(DLBCL$age_years, probs = probs)
# bs_DLBCL <- cuRe:::cb(x = DLBCL$age_years, knots = knots_DLBCL, intercept = FALSE, ortho = T)
# colnames(bs_DLBCL) <- paste0("bs", 1:(length(knots_DLBCL) - 1))
# DLBCL2 <- cbind(DLBCL, bs_DLBCL)
# 
# fit_DLBCL_time <- cuRe:::FlexCureModel(Surv(FU_years, status) ~ bs1 + bs2 + bs3, data = DLBCL2,
#                                        smooth.formula = ~1,
#                                        bhazard = "exp_haz",
#                                        n.knots = 5, n.knots.time = list(bs1 = 2, bs2 = 2, bs3 = 2))
knots.DLBCL <- log(get.knots(DLBCL, k = 5))
knots.time.DLBCL <- log(get.knots(DLBCL, k = 2))
knots.age.DLBCL <- quantile(DLBCL$age_years, probs = probs)

fit_DLBCL_time <- GenFlexCureModel(Surv(FU_years, status) ~ -1, 
                                   data = DLBCL, bhazard = "exp_haz",
                                   smooth.formula = ~ cb(x = log(FU_years), knots = knots.DLBCL),
                                   cr.formula = ~cb(x = age_years, knots = knots.age.DLBCL, intercept = F), 
                                   tvc.formula = ~cb(log(FU_years), knots = knots.time.DLBCL):
                                     cb(x = age_years, knots = knots.age.DLBCL, intercept = F))

fit_DLBCL_time$knots_age <- knots.age.DLBCL

# knots_FL <- quantile(FL$age_years, probs = probs)
# #knots_FL <- quantile(FL$age_years, probs = c(0, 0.333, 0.666, 1))
# bs_FL <- cuRe:::basis(x = FL$age_years, knots = knots_FL, intercept = FALSE, ortho = T)
# colnames(bs_FL) <- paste0("bs", 1:(length(knots_FL) - 1))
# FL2 <- cbind(FL, bs_FL) 
# 
# fit_FL_time <- FlexCureModel(Surv(FU_years, status) ~ bs1 + bs2 + bs3, data = FL2,
#                                     smooth.formula = ~1,
#                                     bhazard = "exp_haz",
#                                     n.knots = 5, n.knots.time = list(bs1 = 2, bs2 = 2, bs3 = 2))


knots.FL <- log(get.knots(FL, k = 5))
knots.time.FL <- log(get.knots(FL, k = 2))
knots.age.FL <- quantile(FL$age_years, probs = probs)

fit_FL_time <- GenFlexCureModel(Surv(FU_years, status) ~ -1, 
                                data = FL, bhazard = "exp_haz",
                                smooth.formula = ~ cb(x = log(FU_years), knots = knots.FL, ortho = F),
                                cr.formula = ~cb(x = age_years, knots = knots.age.FL, intercept = F, ortho = F), 
                                tvc.formula = ~cb(log(FU_years), knots = knots.time.FL, ortho = F):
                                  cb(x = age_years, knots = knots.age.FL, intercept = F, ortho = F))
fit_FL_time$covariance
fit_FL_time$knots_age <- knots.age.FL

# knots_ML <- quantile(ML$age_years, probs = probs)
# bs_ML <- cuRe:::basis(x = ML$age_years, knots = knots_ML, intercept = FALSE, ortho = T)
# colnames(bs_ML) <- paste0("bs", 1:(length(knots_ML) - 1))
# ML2 <- cbind(ML, bs_ML) 
# fit_ML_time <- FlexCureModel(Surv(FU_years, status) ~ bs1 + bs2 + bs3, data = ML2,
#                                     smooth.formula = ~1,
#                                     bhazard = "exp_haz",
#                                     n.knots = 5,
#                                     n.knots.time = list(bs1 = 2, bs2 = 2, bs3 = 2))


knots.ML <- log(get.knots(ML, k = 5))
knots.time.ML <- log(get.knots(ML, k = 2))
knots.age.ML <- quantile(ML$age_years, probs = probs)

fit_ML_time <- GenFlexCureModel(Surv(FU_years, status) ~ -1, 
                                data = ML, bhazard = "exp_haz",
                                smooth.formula = ~ cb(x = log(FU_years), knots = knots.ML),
                                cr.formula = ~cb(x = age_years, knots = knots.age.ML, intercept = F), 
                                tvc.formula = ~cb(log(FU_years), knots = knots.time.ML):
                                  cb(x = age_years, knots = knots.age.ML, intercept = F))
fit_ML_time$covariance
fit_ML_time$knots_age <- knots.age.ML


#Relative survival model
knots.DLBCL <- log(get.knots(DLBCL, k = 6))
knots.time.DLBCL <- log(get.knots(DLBCL, k = 3))
fit_DLBCL_time2 <- stpm2(Surv(FU_years, status) ~ -1, 
                         data = DLBCL, bhazard = DLBCL$exp_haz, 
                         smooth.formula = ~ cb(x = log(FU_years), knots = knots.DLBCL),
                         tvc.formula = ~cb(x = age_years, knots = knots.age.DLBCL, intercept = F):
                           cb(x = log(FU_years), knots = knots.time.DLBCL))

knots.FL <- log(get.knots(FL, k = 6))
knots.time.FL <- log(get.knots(FL, k = 3))
fit_FL_time2 <- stpm2(Surv(FU_years, status) ~ -1, 
                      data = FL, bhazard = FL$exp_haz, 
                      smooth.formula = ~ cb(x = log(FU_years), knots = knots.FL),
                      tvc.formula = ~ cb(x = age_years, knots = knots.age.FL, intercept = F) : 
                        cb(x = log(FU_years), knots = knots.time.FL))

knots.ML <- log(get.knots(ML, k = 6))
knots.time.ML <- log(get.knots(ML, k = 3))
fit_ML_time2 <- stpm2(Surv(FU_years, status) ~ -1, 
                      data = ML, bhazard = ML$exp_haz, 
                      smooth.formula = ~ cb(x = log(FU_years), knots = knots.ML),
                      tvc.formula = ~ cb(x = age_years, knots = knots.age.ML, intercept = F) : 
                        cb(x = log(FU_years), knots = knots.time.ML))

# fit_ML_time2 <- stpm2(Surv(FU_years, status) ~ 1 + ns(age_years, df = 3), data = ML, 
#                        bhazard = ML$exp_haz, df = 5, 
#                        tvc.formula = ~ns(age_years, df = 3):ns(log(FU_years), df = 2))


times <- c(0, 2, 5)
ages <- seq(50, 80, by = 2)


res_time <- lapply(list(fit_DLBCL_time, fit_FL_time, fit_ML_time), function(obj){
  lapply(ages, function(age){
    # bs <- cuRe:::basis(age, knots = obj[[1]]$knots_age, ortho = TRUE, 
    #                    intercept = FALSE, R.inv = attr(obj[[2]], "R.inv"))
    # colnames(bs) <- paste0("bs", 1:ncol(bs))
    # res <- calc.LL(obj[[1]], time = times, ci = F, 
    #         newdata = data.frame(bs, age = age * ayear, 
    #                              sex = "female", 
    #                              year = as.Date("2010-01-01")))$Ests[[1]]
    calc.LL(obj, time = times, var.type = "n", 
            newdata = data.frame(age = age * ayear, age_years = age, 
                                 sex = "female", year = as.Date("2010-01-01")), 
            smooth.exp = F)[[1]]
  })
})

res_time_nelson <- lapply(list(fit_DLBCL_time2, fit_FL_time2, fit_ML_time2), function(fit){
  lapply(ages, function(age){
    calc.LL(fit, time = times, var.type = "n", 
            newdata = data.frame(age_years = age, age = age * ayear, 
                                 sex = "female", 
                                 year = as.Date("2010-01-01")), 
            smooth.exp = F)[[1]]
  })
})

res_time2 <- lapply(res_time, function(x){
  D <- do.call(rbind, x)
  D$age <- rep(ages, each = length(times))
  D$Time <- rep(times, length(ages))
  D
})

res_time2 <- do.call(rbind, res_time2)
res_time2$disease <- rep(c("DLBCL", "FL", "ML"), each = length(ages) * length(times))
res_time2$Time <- factor(res_time2$Time)


res_time_nelson2 <- lapply(res_time_nelson, function(x){
  D <- do.call(rbind, x)
  D$age <- rep(ages, each = length(times))
  D$Time <- rep(times, length(ages))
  D
})

res_time_nelson2 <- do.call(rbind, res_time_nelson2)
res_time_nelson2$disease <- rep(c("DLBCL", "FL", "ML"), each = length(ages) * length(times))
res_time_nelson2$Time <- factor(res_time_nelson2$Time)


res_all <- rbind(res_time2, res_time_nelson2)
models <- c("FMC", "NRS")
res_all$Model <- factor(rep(models, c(nrow(res_time2), nrow(res_time_nelson2))), models)

p <- ggplot(res_all, aes(x = age, y = Estimate, linetype = Time, colour = Model)) + geom_line() + facet_grid(.~disease) +
  xlab("Age at diagnosis (years)") + ylab("Loss of lifetime (years)") +
  theme_bw() +
  theme(legend.position = "bottom", 
        legend.text=element_text(size=15),
        legend.title=element_text(size=15),
        axis.title=element_text(size=17),
        strip.text = element_text(size=15), 
        axis.text = element_text(size = 14), 
        legend.key.size = unit(2,"line")) + 
  scale_color_manual(name = "", values = c("black", "grey"))

pdf(file.path(fig.out, "Time_varyingLOL2.pdf"), width = 10, height = 5.3)
print(p)
dev.off()
