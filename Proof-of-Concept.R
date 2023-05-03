##### PROOF-OF-CONCEPT #####

library(ggplot2)

niche_diff <- function(Aaa, Aba, Abb, Aab) {
  1 - sqrt((Aba*Aab)/(Aaa*Abb))
}

fitness_ineq <- function(Aaa, Aba, Abb, Aab) {
  sqrt((Aaa*Aab)/(Abb*Aba))
}

### Context-independent microbial scenarios -----

##### Baseline Scenario -----

sce1 <- function(Caa, Cba, Cbb, Cab) {
  Aaa <- Caa
  Aba <- Cba
  Abb <- Cbb
  Aab <- Cab
  
  dat1 <- cbind.data.frame(Aaa, Aba, Abb, Aab)
  dat1$niche <- niche_diff(dat1$Aaa, dat1$Aba, dat1$Abb, dat1$Aab)
  dat1$fit <- fitness_ineq(dat1$Aaa, dat1$Aba, dat1$Abb, dat1$Aab)
  
  f <- function(x) 1 / (-x + 1)
  ggplot(data = dat1, aes(x = niche, y = fit)) + 
    geom_point(size = 2, color = "#2F242C") + geom_abline(slope = -1, intercept = 1) + 
    theme_classic() + stat_function(fun=f) + xlim(-.5, .5) + ylim(.5, 1.5) + 
    xlab("Niche Difference") + ylab("Fitness Inequality")
}

# Trying some options
#sce1(-5, -6, -5.5, -5.25)
#sce1(-4, -6, -5.5, -5.25)
#sce1(-4, -6, -5, -4.5)
#sce1(-4, -5, -7.5, -6)
sce1(-6, -5, -7.5, -6) # Starting firmly in coexistence

# this will serve as the baseline parameterization

##### Live Scenario 1 -----

sce2 <- function(Sa1, Sb1, P1a, P1b, 
                 Caa = -6, Cba = -5, Cbb = -7.5, Cab = -6) {
  Aaa <- Caa + Sa1*P1a
  Aba <- Cba + Sb1*P1a
  Abb <- Cbb + Sb1*P1b
  Aab <- Cab + Sa1*P1b
  
  dat1 <- cbind.data.frame(Aaa, Aba, Abb, Aab)
  dat1$niche <- niche_diff(dat1$Aaa, dat1$Aba, dat1$Abb, dat1$Aab)
  dat1$fit <- fitness_ineq(dat1$Aaa, dat1$Aba, dat1$Abb, dat1$Aab)
  
  f <- function(x) 1 / (-x + 1)
  ggplot(data = dat1, aes(x = niche, y = fit)) + 
    geom_point(size = 2, color = "#2F242C") + geom_abline(slope = -1, intercept = 1) + 
    theme_classic() + stat_function(fun=f) + xlim(-.5, .5) + ylim(.5, 1.5) + 
    xlab("Niche Difference") + ylab("Fitness Inequality")
}

# Trying some options
#sce2(-1.1, 0, 1, 0)
#sce2(-2, 0, 1, 0) 
sce2(2, 0, 1, 0) # Getting to exclusion

##### Live Scenario 2 -----

sce3 <- function(Sa2, Sb2, P2a, P2b, 
                 Sa1 = 2, Sb1 = 0, P1a = 1, P1b = 0, 
                 Caa = -6, Cba = -5, Cbb = -7.5, Cab = -6) {
  Aaa <- Caa + Sa1*P1a + Sa2*P2a
  Aba <- Cba + Sb1*P1a + Sb2*P2a
  Abb <- Cbb + Sb1*P1b + Sb2*P2b
  Aab <- Cab + Sa1*P1b + Sa2*P2b
  
  dat1 <- cbind.data.frame(Aaa, Aba, Abb, Aab)
  dat1$niche <- niche_diff(dat1$Aaa, dat1$Aba, dat1$Abb, dat1$Aab)
  dat1$fit <- fitness_ineq(dat1$Aaa, dat1$Aba, dat1$Abb, dat1$Aab)
  
  f <- function(x) 1 / (-x + 1)
  ggplot(data = dat1, aes(x = niche, y = fit)) + 
    geom_point(size = 2, color = "#2F242C") + geom_abline(slope = -1, intercept = 1) + 
    theme_classic() + stat_function(fun=f) + xlim(-.5, .5) + ylim(.5, 1.5) + 
    xlab("Niche Difference") + ylab("Fitness Inequality")
}

sce3(0, 2.5, 0, 1) # Destabilized into priority effects

##### All together -----

Caa <- c(-6, -6, -6)
Cba <- c(-5, -5, -5)
Cbb <- c(-7.5, -7.5, -7.5)
Aab <- c(-6, -6, -6)

Sa1 <- c(0, 2, 2)
Sa2 <- c(0, 0, 0)
Sb1 <- c(0, 0, 0)
Sb2 <- c(0, 0, 2.5)

P1a <- c(0, 1, 1)
P2a <- c(0, 0, 0)
P1b <- c(0, 0, 0)
P2b <- c(0, 0, 1)

Aaa <- Caa + Sa1*P1a + Sa2*P2a
Aba <- Cba + Sb1*P1a + Sb2*P2a
Abb <- Cbb + Sb1*P1b + Sb2*P2b
Aab <- Cab + Sa1*P1b + Sa2*P2b

dat1 <- cbind.data.frame(Aaa, Aba, Abb, Aab)
dat1$niche <- niche_diff(dat1$Aaa, dat1$Aba, dat1$Abb, dat1$Aab)
dat1$fit <- fitness_ineq(dat1$Aaa, dat1$Aba, dat1$Abb, dat1$Aab)

f <- function(x) 1 / (-x + 1)
ggplot(data = dat1, aes(x = niche, y = fit)) + 
  geom_point(size = 2, color = "#2F242C") + geom_abline(slope = -1, intercept = 1) + 
  theme_classic() + stat_function(fun=f) + xlim(-.375, .375) + ylim(.75, 1.25) + 
  xlab("Niche Difference") + ylab("Fitness Inequality")

### Microbe-independent context-depency -----

eve <- function(mCaa, mCba, mCbb, mCab, 
                Caa = -6, Cba = -5, Cbb = -7.5, Cab = -6){
  v <- seq(0, 1, .01)
  
  # Adding linear scaling of C's with slopes mC's
  Aaa <- v*mCaa + rep(Caa, length(v))
  Aba <- v*mCba + rep(Cba, length(v))
  Abb <- v*mCbb + rep(Cbb, length(v))
  Aab <- v*mCab + rep(Cab, length(v))
  
  dat <- cbind.data.frame(Aaa, Aba, Abb, Aab, v)
  
  dat$niche <- niche_diff(dat$Aaa, dat$Aba, dat$Abb, dat$Aab)
  dat$fit <- fitness_ineq(dat$Aaa, dat$Aba, dat$Abb, dat$Aab)
  
  dat$niche[dat$Aaa > 0 | dat$Abb > 0 | dat$Aab > 0 | dat$Aba > 0] <- NaN
  dat$fit[dat$Aaa > 0 | dat$Abb > 0 | dat$Aab > 0 | dat$Aba > 0] <- NaN
  
  f <- function(x) 1 / (-x + 1)
  ggplot(data = dat, aes(x = niche, y = fit, label = v+2.5, color = v)) + 
    geom_point() + geom_abline(slope = -1, intercept = 1) + 
    theme_classic() + stat_function(fun=f) + xlim(-.375, .375) + ylim(.75, 1.25) + 
    xlab("Niche Difference") + ylab("Fitness Inequality") + 
    scale_color_gradient(low = "#E8CCF5", high = "#260A33") + 
    theme(legend.position = "none")
  
}

eve(2.5, 1, 1, 1.5) # Competitive effects on sp A are alleviated by 
# moisture much more than competitive effects on sp B

### Context-dependent microbes -----

mic_eve <- function(mSa1, mSa2, mSb1, mSb2,
                    Caa = -6, Cba = -5, Cbb = -7.5, Cab = -6, 
                    Sa1 = 2, Sb1 = 0, P1a = 1, P1b = 0,
                    Sa2 = 0, Sb2 = 2.5, P2a = 0, P2b = 1){
  v <- seq(0, 1, .01)
  
  Aaa <- rep(Caa, length(v)) + (rep(Sa1, length(v)) + v*mSa1)*(rep(P1a, length(v))) + (rep(Sa2, length(v)) + v*mSa2)*(rep(P2a, length(v)))
  Aba <- rep(Cba, length(v)) + (rep(Sb1, length(v)) + v*mSb1)*(rep(P1a, length(v))) + (rep(Sb2, length(v)) + v*mSb2)*(rep(P2a, length(v)))
  Abb <- rep(Cbb, length(v)) + (rep(Sb1, length(v)) + v*mSb1)*(rep(P1b, length(v))) + (rep(Sb2, length(v)) + v*mSb2)*(rep(P2b, length(v)))
  Aab <- rep(Cab, length(v)) + (rep(Sa1, length(v)) + v*mSa1)*(rep(P1b, length(v))) + (rep(Sa2, length(v)) + v*mSa2)*(rep(P2b, length(v)))
  
  dat <- cbind.data.frame(Aaa, Aba, Abb, Aab, v)
  
  dat$niche <- niche_diff(dat$Aaa, dat$Aba, dat$Abb, dat$Aab)
  dat$fit <- fitness_ineq(dat$Aaa, dat$Aba, dat$Abb, dat$Aab)
  
  dat$niche[dat$Aaa > 0 | dat$Abb > 0 | dat$Aab > 0 | dat$Aba > 0] <- NaN
  dat$fit[dat$Aaa > 0 | dat$Abb > 0 | dat$Aab > 0 | dat$Aba > 0] <- NaN
  
  f <- function(x) 1 / (-x + 1)
  ggplot(data = dat, aes(x = niche, y = fit, label = v+2.5, color = v)) + 
    geom_point() + geom_abline(slope = -1, intercept = 1) + 
    theme_classic() + stat_function(fun=f) + xlim(-.375, .375) + ylim(.75, 1.25) + 
    xlab("Niche Difference") + ylab("Fitness Inequality") + 
    scale_color_gradient(low = "#E8CCF5", high = "#260A33") + 
    theme(legend.position = "none")
  
}

mic_eve(-1.5, 0, 0, -3)

