dir_fun <- './R_code'
dir_results <- './RESULTS'
dir_coeff <- file.path(dir_results, 'COEFFS')

source(file.path(dir_fun, 'functions.R'))

posterior_ss <- read.csv(file.path(dir_coeff, "coefficients_posterior_SS.csv"))
posterior_rev <- read.csv(file.path(dir_coeff, "coefficients_posterior_REV.csv"))
posterior_nm <- read.csv(file.path(dir_coeff, "coefficients_posterior_NM.csv"))

coeffs_ss <- colMeans(posterior_ss)
coeffs_rev <- colMeans(posterior_rev)
coeffs_nm <- colMeans(posterior_nm)

u <- 0.4
mag <- 8
mag_list <- seq(5,8.5,0.1)
u_list <- seq(0,1,0.05)
seed <- 1701
n_sample <- 1000

# reverse faulting
post <- posterior_rev
coeffs <- as.data.frame(t(coeffs_rev))

func_rev(mag, u, coeffs, n_sample = 1000)

# normal faulting
post <- posterior_nm
coeffs <- as.data.frame(t(coeffs_nm))

func_nm(mag, u, coeffs, n_sample = 1000)

# strike-slip faulting
post <- posterior_ss
coeffs <- as.data.frame(t(coeffs_ss))

func_ss(mag, u, coeffs)
