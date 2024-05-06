
set.seed(6)
library(epiR)
library(tidyverse)
library(scales)

# Sensitivity Analysis for Collings et al. 2024

niche_diff <- function(Aaa, Aba, Abb, Aab) {
  ifelse(Aaa < 0 & Aba < 0 & Abb < 0 & Aab < 0, 
         1 - sqrt((Aba*Aab)/(Aaa*Abb)), NA)
}

fitness_ineq <- function(Aaa, Aba, Abb, Aab) {
  ifelse(Aaa < 0 & Aba < 0 & Abb < 0 & Aab < 0, 
         ifelse(sqrt((Aaa*Aab)/(Abb*Aba)) > sqrt((Abb*Aba)/(Aaa*Aab)), 
                sqrt((Aaa*Aab)/(Abb*Aba)), sqrt((Abb*Aba)/(Aaa*Aab))), NA)
}

fitness_ratio <- function(Aaa, Aba, Abb, Aab) {
  ifelse(Aaa < 0 & Aba < 0 & Abb < 0 & Aab < 0, 
         sqrt((Aaa*Aab)/(Abb*Aba)), NA)
}

cube_root <- function(x){
  ifelse(x < 0, -(abs(x)^(1/3)), x^(1/3))
}

std_params <- function(mymin, mymax, xvals) {
  m <- 2/(mymax - mymin)
  b <- 1- (2*mymax)/(mymax - mymin)
  y <- m*xvals+b
  return(y)
}

### Visual Sensitivity Analysis -----

# as a starting point, let's modify a bit the first live scenario in the 
# theory simulations

sce <- function(Sax, Sbx, Pxa, Pxb, 
                 Caa = -.06, Cba = -.05, Cbb = -.075, Cab = -.06) {
  Aaa <- Caa + Sax*Pxa
  Aba <- Cba + Sbx*Pxa
  Abb <- Cbb + Sbx*Pxb
  Aab <- Cab + Sax*Pxb
  
  dat1 <- cbind.data.frame(Aaa, Aba, Abb, Aab)
  dat1$niche <- niche_diff(dat1$Aaa, dat1$Aba, dat1$Abb, dat1$Aab)
  dat1$fit <- fitness_ratio(dat1$Aaa, dat1$Aba, dat1$Abb, dat1$Aab)
  
  f <- function(x) 1 / (-x + 1)
  ggplot(data = dat1, aes(x = niche, y = fit)) + 
    geom_point(size = 2, color = "#2F242C") + geom_abline(slope = -1, intercept = 1) + 
    theme_classic() + stat_function(fun=f) + xlim(-.5, .5) + ylim(.5, 1.5) + 
    xlab("Niche Difference") + ylab("Fitness Ratio")
}

sce(.002, .0025, 10.5, 10) # Getting to exclusion

# this is nearly identical... except now we have a genaralist microbe
# for the baseline conditions, let's assume species are pretty similar
# but not identical in their interactions with the microbe

Caa <- rep(-0.06, 1001)
Cba <- rep(-0.05, 1001)
Cbb <- rep(-0.075, 1001)
Cab <- rep(-0.055, 1001)
Sax <- rep(0.002, 1001)
Sbx <- rep(0.00225, 1001)
Pxa <- rep(10.5, 1001)
Pxb <- rep(10, 1001)

##### Niche Differences -----

Caa2 <- seq(-.1, 0, .0001)
Aaa <- Caa2 + Sax*Pxa
Aba <- Cba + Sbx*Pxa
Abb <- Cbb + Sbx*Pxb
Aab <- Cab + Sax*Pxb
ND_Caa <- niche_diff(Aaa, Aba, Abb, Aab)
FI_Caa <- fitness_ineq(Aaa, Aba, Abb, Aab)
FR_Caa <- fitness_ratio(Aaa, Aba, Abb, Aab)

Cba2 <- seq(-.1, 0, .0001)
Aaa <- Caa + Sax*Pxa
Aba <- Cba2 + Sbx*Pxa
Abb <- Cbb + Sbx*Pxb
Aab <- Cab + Sax*Pxb
ND_Cba <- niche_diff(Aaa, Aba, Abb, Aab)
FI_Cba <- fitness_ineq(Aaa, Aba, Abb, Aab)
FR_Cba <- fitness_ratio(Aaa, Aba, Abb, Aab)

Cbb2 <- seq(-.1, 0, .0001)
Aaa <- Caa + Sax*Pxa
Aba <- Cba + Sbx*Pxa
Abb <- Cbb2 + Sbx*Pxb
Aab <- Cab + Sax*Pxb
ND_Cbb <- niche_diff(Aaa, Aba, Abb, Aab)
FI_Cbb <- fitness_ineq(Aaa, Aba, Abb, Aab)
FR_Cbb <- fitness_ratio(Aaa, Aba, Abb, Aab)

Cab2 <- seq(-.1, 0, .0001)
Aaa <- Caa + Sax*Pxa
Aba <- Cba + Sbx*Pxa
Abb <- Cbb + Sbx*Pxb
Aab <- Cab2 + Sax*Pxb
ND_Cab <- niche_diff(Aaa, Aba, Abb, Aab)
FI_Cab <- fitness_ineq(Aaa, Aba, Abb, Aab)
FR_Cab <- fitness_ratio(Aaa, Aba, Abb, Aab)

Sax2 <- seq(-.01, 0.01, .00002)
Aaa <- Caa + Sax2*Pxa
Aba <- Cba + Sbx*Pxa
Abb <- Cbb + Sbx*Pxb
Aab <- Cab + Sax2*Pxb
ND_Sax <- niche_diff(Aaa, Aba, Abb, Aab)
FI_Sax <- fitness_ineq(Aaa, Aba, Abb, Aab)
FR_Sax <- fitness_ratio(Aaa, Aba, Abb, Aab)

Sbx2 <- seq(-.01, 0.01, .00002)
Aaa <- Caa + Sax*Pxa
Aba <- Cba + Sbx2*Pxa
Abb <- Cbb + Sbx2*Pxb
Aab <- Cab + Sax*Pxb
ND_Sbx <- niche_diff(Aaa, Aba, Abb, Aab)
FI_Sbx <- fitness_ineq(Aaa, Aba, Abb, Aab)
FR_Sbx <- fitness_ratio(Aaa, Aba, Abb, Aab)

Pxa2 <- seq(0, 20, .02)
Aaa <- Caa + Sax*Pxa2
Aba <- Cba + Sbx*Pxa2
Abb <- Cbb + Sbx*Pxb
Aab <- Cab + Sax*Pxb
ND_Pxa <- niche_diff(Aaa, Aba, Abb, Aab)
FI_Pxa <- fitness_ineq(Aaa, Aba, Abb, Aab)
FR_Pxa <- fitness_ratio(Aaa, Aba, Abb, Aab)

Pxb2 <- seq(0, 20, .02)
Aaa <- Caa + Sax*Pxa
Aba <- Cba + Sbx*Pxa
Abb <- Cbb + Sbx*Pxb2
Aab <- Cab + Sax*Pxb2
ND_Pxb <- niche_diff(Aaa, Aba, Abb, Aab)
FI_Pxb <- fitness_ineq(Aaa, Aba, Abb, Aab)
FR_Pxb <- fitness_ratio(Aaa, Aba, Abb, Aab)

sens_dat <- cbind.data.frame(NDs = c(ND_Caa, ND_Cba, 
                                   ND_Cbb, ND_Cab, 
                                   ND_Sax, ND_Sbx, 
                                   ND_Pxa, ND_Pxb),
                           FIs = c(FI_Caa, FI_Cba, 
                                   FI_Cbb, FI_Cab, 
                                   FI_Sax, FI_Sbx, 
                                   FI_Pxa, FI_Pxb),
                           FRs = c(FR_Caa, FR_Cba, 
                                   FR_Cbb, FR_Cab, 
                                   FR_Sax, FR_Sbx, 
                                   FR_Pxa, FR_Pxb), 
                           par = rep(c("Caa", "Cba", "Cbb", "Cab", 
                                         "Sax", "Sbx", "Pxa", "Pxb"), each = 1001), 
                           val = c(std_params(0, 0.1, abs(Caa2)), 
                                   std_params(0, 0.1, abs(Cba2)), 
                                   std_params(0, 0.1, abs(Cbb2)), 
                                   std_params(0, 0.1, abs(Cab2)), 
                                   std_params(-0.01, 0.01, Sax2), 
                                   std_params(-0.01, 0.01, Sbx2), 
                                   std_params(0, 20, Pxa2), 
                                   std_params(0, 20, Pxb2)), 
                           type = rep(c("plant", "microbe"), each = 1001*4))

ggplot(data = sens_dat, aes(x = val, y = NDs, color = par, linetype = type)) + 
  geom_line(size = 1.5) + ylim(-1, 1) + 
  theme_classic(base_size = 12) + ylab("Niche Differences") + 
  xlab("Parameter Value") + theme_classic(base_size = 15) + 
  theme(axis.text.x = element_blank()) +
  scale_color_manual(name = "Parameter", 
                     values = c("#844A9D", "#B990CB", 
                                "#462753", "#DDC7E5", 
                                "#F5D544", "#D5B00B", 
                                "#FF7F11", "#FFB370"))

ggplot(data = sens_dat, aes(x = val, y = FIs, color = par, linetype = type)) + 
  geom_line(size = 1.5) + ylim(1, 2) + 
  theme_classic(base_size = 12) + ylab("Fitness Inequalities") + 
  xlab("Parameter Value") + theme_classic(base_size = 15) + 
  theme(axis.text.x = element_blank()) + 
  scale_color_manual(name = "Parameter", 
                     values = c("#844A9D", "#B990CB", 
                                "#462753", "#DDC7E5", 
                                "#F5D544", "#D5B00B", 
                                "#FF7F11", "#FFB370"))

ggplot(data = sens_dat, aes(x = val, y = FRs, color = par, linetype = type)) + 
  geom_line(size = 1.5) + ylim(0, 2) + 
  theme_classic(base_size = 12) + ylab("Fitness Ratios") + 
  xlab("Parameter Value") + theme_classic(base_size = 15) + 
  theme(axis.text.x = element_blank()) + 
  scale_color_manual(name = "Parameter", 
                     values = c("#844A9D", "#B990CB", 
                                "#462753", "#DDC7E5", 
                                "#F5D544", "#D5B00B", 
                                "#FF7F11", "#FFB370"))

### Local Sensitivity Analysis -----

##### Niche Differences -----

ND_sens <- deriv(expression(1 - sqrt(((c_ij + sig_ix*phi_xj)*(c_ji + sig_jx*phi_xi))/
                            ((c_ii + sig_ix*phi_xi)*(c_jj + sig_jx*phi_xj)))), 
      c("c_ij", "c_ji", "c_jj", "c_ii", "sig_ix", "sig_jx", "phi_xi", "phi_xj"))

c_ii = -0.06; c_jj = -0.075; c_ji = -0.05; c_ij = -0.055
sig_ix = 0.002; sig_jx = 0.00225; phi_xi = 10.5; phi_xj = 10

ggplot(data = cbind.data.frame(vals = as.numeric(attr(eval(ND_sens), "gradient")), 
                               pars = c("c_ij", "c_ji", "c_jj", "c_ii", 
                                        "sig_ix", "sig_jx", "phi_xi", "phi_xj")), 
       aes(x = pars, y = vals)) + 
         geom_bar(stat = "identity") + 
  theme_classic(base_size = 15) + ylab("Sensitivity of Niche Differences") + xlab("Parameter")

##### Fitness Ratios -----

FR_sens <- deriv(expression(sqrt(((c_ij + sig_ix*phi_xj)*(c_ii + sig_ix*phi_xi))/
                                       ((c_ji + sig_jx*phi_xi)*(c_jj + sig_jx*phi_xj)))), 
                 c("c_ij", "c_ji", "c_jj", "c_ii", "sig_ix", "sig_jx", "phi_xi", "phi_xj"))

c_ii = -0.06; c_jj = -0.075; c_ji = -0.05; c_ij = -0.055
sig_ix = 0.002; sig_jx = 0.00225; phi_xi = 10.5; phi_xj = 10

FR_sens_values <- as.numeric(attr(eval(FR_sens), "gradient"))
ND_sens_values <- as.numeric(attr(eval(ND_sens), "gradient"))

ggplot(data = cbind.data.frame(vals = c(-FR_sens_values[1:4], FR_sens_values[5:8]), 
                               pars = c("c_ij", "c_ji", "c_jj", "c_ii", 
                                        "sig_ix", "sig_jx", "phi_xi", "phi_xj")), 
       aes(x = pars, y = vals)) + 
  geom_bar(stat = "identity") + 
  theme_classic(base_size = 15) + ylab("Sensitivity of Fitness Ratio") + xlab("Parameter")

ggplot(data = cbind.data.frame(vals = c(c(-FR_sens_values[1:4], FR_sens_values[5:8]), 
                                        c(-ND_sens_values[1:4], ND_sens_values[5:8])), 
                               pars = rep(c("c_ij", "c_ji", "c_jj", "c_ii", 
                                        "sig_ix", "sig_jx", "phi_xi", "phi_xj"), 2), 
                               metric = c(rep("FR", 8), rep("ND", 8))), 
       aes(x = pars, y = vals, fill = metric)) + 
  geom_bar(stat = "identity", position = "dodge") + 
  scale_y_continuous(trans= pseudo_log_trans(), 
                     breaks = c(-250, -100, -10, -1, 0, 1, 10, 100, 250, 500)) + 
  theme_classic(base_size = 15) + ylab("Pseudo-Log Sensitivity") + xlab("Parameter") + 
  scale_fill_manual(name = "Metric", labels = c("Fitness Ratio", "Niche Differences"), 
                    values = c("#F5D544", "#844A9D"))

ggplot(data = cbind.data.frame(vals = c(c(-FR_sens_values[1:4], FR_sens_values[5:8]), 
                                        c(-ND_sens_values[1:4], ND_sens_values[5:8])), 
                               pars = rep(c("c_ij", "c_ji", "c_jj", "c_ii", 
                                            "sig_ix", "sig_jx", "phi_xi", "phi_xj"), 2), 
                               metric = c(rep("FR", 8), rep("ND", 8))), 
       aes(x = pars, y = vals, fill = metric)) + 
  geom_bar(stat = "identity", position = "dodge") + 
  theme_classic(base_size = 15) + ylab("Sensitivity") + xlab("Parameter") + 
  scale_fill_manual(name = "Metric", labels = c("Fitness Ratio", "Niche Differences"), 
                    values = c("#F5D544", "#844A9D"))

### Regression Sensitivity Analysis -----

c_min <- -.1; c_max <- 0
sig_min <- -.01; sig_max <- .01
phi_min <- 0; phi_max <- 20

# some empty vectors to throw values into
cAAs <- cBAs <- cBBs <- cABs <- c()
sigAs <- sigBs <- phiAs <- phiBs <- c()
alphaAAs <- alphaBAs <- alphaBBs <- alphaABs <- c()
NDs <- FIs <- c()

cutoff <- 0
# generating data from 5000 runs
runs <- 100000
for(i in 1:runs){
  cAAs[i] <- runif(1, min = c_min, max = c_max)
  cBAs[i] <- runif(1, min = c_min, max = c_max)
  cBBs[i] <- runif(1, min = c_min, max = c_max)
  cABs[i] <- runif(1, min = c_min, max = c_max)
  
  sigAs[i] <- runif(1, min = sig_min, max = sig_max)
  sigBs[i] <- runif(1, min = sig_min, max = sig_max)
  
  phiAs[i] <- runif(1, min = phi_min, max = phi_max)
  phiBs[i] <- runif(1, min = phi_min, max = phi_max)
  
  alphaAAs[i] <- cAAs[i]+sigAs[i]*phiAs[i]
  alphaBAs[i] <- cBAs[i]+sigBs[i]*phiAs[i]
  alphaBBs[i] <- cBBs[i]+sigBs[i]*phiBs[i]
  alphaABs[i] <- cABs[i]+sigAs[i]*phiBs[i]
  
  cutoff <- -0.0075 #trying arbitrary cutoffs to tame FI
  NDs[i] <- ifelse(alphaABs[i] < cutoff & alphaBAs[i] < cutoff & alphaAAs[i] < cutoff & alphaBBs[i] < cutoff,
                   1 - sqrt((alphaABs[i]*alphaBAs[i])/(alphaAAs[i]*alphaBBs[i])), NA)
  FIs[i] <- ifelse(alphaABs[i] < cutoff & alphaBAs[i] < cutoff & alphaAAs[i] < cutoff & alphaBBs[i] < cutoff, sqrt((alphaAAs[i]*alphaABs[i])/(alphaBBs[i]*alphaBAs[i])), NA)
  
}

length(alphaAAs) - sum(alphaABs < cutoff & alphaBAs < cutoff & alphaAAs < cutoff & alphaBBs < cutoff)
(length(alphaAAs) - sum(alphaABs < cutoff & alphaBAs < cutoff & alphaAAs < cutoff & alphaBBs < cutoff))/length(alphaAAs)

PRCC_dat <- cbind.data.frame(cAA = scale(abs(cAAs)), cBA = scale(abs(cBAs)), 
                        cBB = scale(abs(cBBs)), cAB = scale(abs(cABs)),
                        sigA = scale(sigAs), sigB = scale(sigBs), 
                        phiA = scale(phiAs), phiB = scale(phiBs), 
                        alphaAA = alphaAAs, alphaBA = alphaBAs, 
                        alphaBB = alphaBBs, alphaAB = alphaABs, 
                        ND = NDs, FI = FIs)

PRCC_dat2 <- rbind.data.frame(epi.prcc(PRCC_dat[complete.cases(PRCC_dat),c(1:8,13)]), 
                         epi.prcc(PRCC_dat[complete.cases(PRCC_dat),c(1:8,14)]))

PRCC_dat2$metric <- c(rep("ND", 8), rep("FR", 8))

ggplot(data = PRCC_dat2, 
       aes(x = var, y = est, ymin = lower, ymax = upper, fill = metric)) + 
  geom_bar(stat = "identity", position = position_dodge(width = .9)) + 
  geom_errorbar(width = 0, position = position_dodge(width = .9)) + 
  theme_classic(base_size = 15) + ylab("PRCC Sensitivity") + xlab("Parameter") + 
  scale_fill_manual(name = "Metric", labels = c("Fitness Ratio", "Niche Differences"), 
                    values = c("#F5D544", "#844A9D"))

