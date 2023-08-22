func_mag <- function(x, params, mb = 7, delta = 0.1) {
  # Function to caculate the magnitude scaling f the maximum
  # Input parameters:
  # x: value of predictor variable (magnitude)
  # params: vector of coefficients, length 3 (c1,c2,c3)
  # mb: value of break point (default 7)
  # delta: paameter controlling the smoothness of the transition (default 0.1)
  fm <- params[1] + params[2]*(x - mb) +
    (params[3] - params[2])*delta*log(1 + exp((x - mb)/delta))
  return(fm)
}

func_u <- function(x, mode, c, alpha, beta) {
  # function toc alculate slip dependent on U
  # Input parameters:
  # x: value of predictor variable (U_star or x/L)
  # mode: predicted maximum of average slip profile
  # c: coefficent controlling shape of profile
  # alpha: coefficient
  # beta: coefficient
  a <- mode - c * (alpha / (alpha + beta))^alpha * (beta / (alpha + beta))^beta
  mu <- a + c * x^(alpha) * (1 - x)^beta
  return(mu)
}

func_sd_mode_sigmoid <- function(x, params, mb = 7) {
  # function to calculate the standard deviation of the mode
  # which is constant at low M and constant at large M with a smooth transition
  # Input parameters:
  # x: value of predictor variable (magnitude)
  # params: vector of coefficients, length 3 (s_m1, s_m2, s_m3)
  sigma <- params[1] - params[2] / (1 + exp(-params[3] * (x - mb)))
  return(sigma)
}

func_sd_mode_bilinear <- function(x, params, delta = 0.1) {
  # function to calculate the standard deviation of the mode
  # with a bilinear function
  # Input parameters:
  # x: value of predictor variable (magnitude)
  # params: vector of coefficients, length 3 (s_m1, s_m2, s_m3)
  sigma <- params[1] + params[2] * (x - params[3]) -
    params[2] * delta * log(1 + exp((x - params[3]) / delta))
  return(sigma)
}

func_sd_u <- function(x, params, alpha, beta) {
  # function to calculate the standard deviation of the
  # observation likelihood dependent on U
  # Input parameters:
  # x: value of predictor variable
  # params: coefficients (s_1, s_2)
  # alpha, beta: parameters providing the minimum
  sd <- params[1] + params[2] * (x - alpha / (alpha + beta))^2
  return(sd)
}

func_ss <- function(mag, u, coeffs) {
  # function to calculate median prediction and standard deviation for strike-slip events
  # Inpu:
  # mag: magnitude
  # u: position along the rupture
  # coeffs: data frame with coefficients, named like in report
  mode <- func_mag(mag, c(coeffs$c1, coeffs$c2, coeffs$c3))
  
  # run for input u "peak on left"
  med_left <- func_u(u, mode, coeffs$c, coeffs$alpha, coeffs$beta)
  sd_mode <- func_sd_mode_bilinear(mag, c(coeffs$s_m.s1, coeffs$s_m.s2, coeffs$s_m.s3))
  sd_u_left <- func_sd_u(u, c(coeffs$s_s1, coeffs$s_s2), coeffs$alpha, coeffs$beta)
  sd_total_left <- sqrt(sd_mode^2 + sd_u_left^2)
  
  # run for complement u "peak on relative right"
  med_right <- func_u(1-u, mode, coeffs$c, coeffs$alpha, coeffs$beta)
  # sd_mode <- func_sd_mode_bilinear(mag, c(coeffs$s_m.s1, coeffs$s_m.s2, coeffs$s_m.s3))
  sd_u_right <- func_sd_u(1-u, c(coeffs$s_s1, coeffs$s_s2), coeffs$alpha, coeffs$beta)
  sd_total_right <- sqrt(sd_mode^2 + sd_u_right^2)
  
  # get log-space average median (note, do not do for sigmas)
  med_avg <- mean(c(med_left, med_right))
  
  # returns
  return(list(ln_median_left = med_left, 
              sigma_mode = sd_mode, 
              sigma_u_left = sd_u_left, 
              sigma_total_left = sd_total_left,
              ln_median_right = med_right, 
              sigma_u_right = sd_u_right, 
              sigma_total_right = sd_total_right,
              ln_median_average = med_avg)
         )
}

func_rev <- function(mag, u, coeffs, n_sample = 1000, seed = 1701) {
  # function to calculate median prediction and standard deviation for reverse events
  # Input:
  # mag: magnitude
  # u: position along the rupture
  # coeffs: data frame with coefficients, named like in report
  # seed: seed value for sampling of c
  # n_sample: how many c-values to sample
  set.seed(seed)
  mode <- func_mag(mag, c(coeffs$c1, coeffs$c2, coeffs$c3))
  
  # sample c values
  if(n_sample == 0) {
    c = coeffs$mu_c
  } else {
    c <- truncnorm::rtruncnorm(n_sample, a = 0, mean = coeffs$mu_c, sd = coeffs$sigma_c)
  }
  
  # run for input u "peak on left"
  med_c_left <- func_u(u, mode, c, coeffs$alpha, coeffs$beta)
  med_left <- median(med_c_left)
  sd_c_left <- sd(med_c_left)
  
  sd_mode <- coeffs$s_m.r
  sd_u_left <- func_sd_u(u, c(coeffs$s_r1, coeffs$s_r2), coeffs$alpha, coeffs$beta)
  sd_total_left <- sqrt(sd_mode^2 + sd_u_left^2 + sd_c_left^2)
  
  # run for complement u "peak on relative right"
  med_c_right <- func_u(1-u, mode, c, coeffs$alpha, coeffs$beta)
  med_right <- median(med_c_right)
  sd_c_right <- sd(med_c_right)
  
  # sd_mode <- coeffs$s_m.r
  sd_u_right <- func_sd_u(1-u, c(coeffs$s_r1, coeffs$s_r2), coeffs$alpha, coeffs$beta)
  sd_total_right <- sqrt(sd_mode^2 + sd_u_right^2 + sd_c_right^2)
  
  # get log-space average median (note, do not do for sigmas)
  med_avg <- mean(c(med_left, med_right))
  
  # returns
  return(list(ln_median_left = med_left, 
              sigma_mode = sd_mode, 
              sigma_c_left = sd_c_left,
              sigma_u_left = sd_u_left, 
              sigma_total_left = sd_total_left,
              ln_median_right = med_right,
              sigma_c_right = sd_c_right,
              sigma_u_right = sd_u_right, 
              sigma_total_right = sd_total_right,
              ln_median_average = med_avg)
         )
}

func_nm <- function(mag, u, coeffs, n_sample = 1000, seed = 1701) {
  # function to calculate median prediction and standard deviation for reverse events
  # Input:
  # mag: magnitude
  # u: position along the rupture
  # coeffs: data frame with coefficients, named like in report
  # seed: seed value for sampling of c
  # n_sample: how many c-values to sample
  set.seed(seed)
  mode <- func_mag(mag, c(coeffs$c1, coeffs$c2, coeffs$c3))
  
  # sample c values
  if(n_sample == 0) {
    c = coeffs$mu_c
  } else {
    c <- truncnorm::rtruncnorm(n_sample, a = 0, mean = coeffs$mu_c, sd = coeffs$sigma_c)
  }
  
  # run for input u "peak on left"
  med_c_left <- func_u(u, mode, c, coeffs$alpha, coeffs$beta)
  med_left <- median(med_c_left)
  sd_c_left <- sd(med_c_left)
  
  sd_mode <- func_sd_mode_sigmoid(mag, c(coeffs$s_m.n1, coeffs$s_m.n2, coeffs$s_m.n3))
  sd_u <- coeffs$sigma
  sd_total_left <- sqrt(sd_mode^2 + sd_u^2 + sd_c_left^2)
  
  # run for complement u "peak on relative right"
  med_c_right <- func_u(1-u, mode, c, coeffs$alpha, coeffs$beta)
  med_right <- median(med_c_right)
  sd_c_right <- sd(med_c_right)
  
  # sd_mode <- func_sd_mode_sigmoid(mag, c(coeffs$s_m.n1, coeffs$s_m.n2, coeffs$s_m.n3))
  # sd_u <- coeffs$sigma
  sd_total_right <- sqrt(sd_mode^2 + sd_u^2 + sd_c_right^2)
  
  # get log-space average median (note, do not do for sigmas)
  med_avg <- mean(c(med_left, med_right))
  
  # returns
  return(list(ln_median_left = med_left, 
              sigma_mode = sd_mode, 
              sigma_c_left = sd_c_left,
              sigma_u = sd_u, 
              sigma_total_left = sd_total_left,
              ln_median_right = med_right, 
              sigma_c_right = sd_c_right,
              sigma_total_right = sd_total_right,
              ln_median_average = med_avg)
         )
}

func_rev2 <- function(mag, u, coeffs, delta_med_c, sigma_c, u_c) {
  # function to calculate median prediction and standard deviation for reverse events
  # Input:
  # mag: magnitude
  # u: position along the rupture
  # coeffs: data frame with coefficients, named like in report
  # delta_med_c: median adjustment due to event terms c,
  # sigma_c: additional siga due to event terms c, same format as delta_med_c
  # u_c: values of U_star for which delta_med_c and sigma_c are calculated
  mode <- func_mag(mag, c(coeffs$c1, coeffs$c2, coeffs$c3))
  
  
  med <- func_u(u, mode, coeffs$mu_c, coeffs$alpha, coeffs$beta) + 
    approxfun(u_c, delta_med_c)(u)
  sd_c <- approxfun(u_c, sigma_c)(u)
  
  sd_mode <- coeffs$s_m.r
  sd_u <- func_sd_u(u, c(coeffs$s_r1, coeffs$s_r2), coeffs$alpha, coeffs$beta)
  sd_total <- sqrt(sd_mode^2 + sd_u^2 + sd_c^2)
  return(list(median = med, sigma_mode = sd_mode, sigma_c = sd_c,
              sigma_u = sd_u, sigma_total = sd_total))
}


func_nm2 <- function(mag, u, coeffs, delta_med_c, sigma_c, u_c) {
  # function to calculate median prediction and standard deviation for reverse events
  # Input:
  # mag: magnitude
  # u: position along the rupture
  # coeffs: data frame with coefficients, named like in report
  # delta_med_c: median adjustment due to event terms c,
  # sigma_c: additional siga due to event terms c, same format as delta_med_c
  # u_c: values of U_star for which delta_med_c and sigma_c are calculated
  
  mode <- func_mag(mag, c(coeffs$c1, coeffs$c2, coeffs$c3))
  med <- func_u(u, mode, coeffs$mu_c, coeffs$alpha, coeffs$beta) + 
    approxfun(u_c, delta_med_c)(u)
  sd_c <- approxfun(u_c, sigma_c)(u)
  
  sd_mode <- func_sd_mode_sigmoid(mag, c(coeffs$s_m.n1, coeffs$s_m.n2, coeffs$s_m.n3))
  sd_u <- coeffs$sigma
  sd_total <- sqrt(sd_mode^2 + sd_u^2 + sd_c^2)
  return(list(median = med, sigma_mode = sd_mode, sigma_c = sd_c,
              sigma_u = sd_u, sigma_total = sd_total))
}


func_ss2 <- function(mag, u, coeffs, delta_med_c, sigma_c, u_c) {
  # function to calculate median prediction and standard deviation for strike-slip events
  # Input:
  # mag: magnitude
  # u: position along the rupture
  # coeffs: data frame with coefficients, named like in report
  # delta_med_c: median adjustment due to event terms c, (not used)
  # sigma_c: additional siga due to event terms c, same format as delta_med_c (not used)
  # u_c: values of U_star for which delta_med_c and sigma_c are calculated (not used)
  mode <- func_mag(mag, c(coeffs$c1, coeffs$c2, coeffs$c3))
  
  # run for input u "peak on left"
  med <- func_u(u, mode, coeffs$c, coeffs$alpha, coeffs$beta)
  sd_mode <- func_sd_mode_bilinear(mag, c(coeffs$s_m.s1, coeffs$s_m.s2, coeffs$s_m.s3))
  sd_u <- func_sd_u(u, c(coeffs$s_s1, coeffs$s_s2), coeffs$alpha, coeffs$beta)
  sd_total <- sqrt(sd_mode^2 + sd_u^2)
  
  return(list(median = med, sigma_mode = sd_mode,
              sigma_u = sd_u, sigma_total = sd_total))
}
