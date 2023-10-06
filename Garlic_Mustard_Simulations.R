##### GARLIC MUSTARD SIMULATIONS #####

library(ggplot2)

coexist_plot <- function(d, sterile = FALSE){
  # d is a dataframe that includes the following:
  # Microbe independent competitive effects ("Caa", "Cab", "Cbb", "Cba")
  # Competition coefficients ("Aaa", "Aab", "Abb", "Aba")
  # Environmental gradient ("v")
  d$niche <- 1 - sqrt((d$Aba * d$Aab)/(d$Aaa * d$Abb))
  d$fit <- sqrt((d$Aaa * d$Aab) / (d$Abb * d$Aba))
  d$niche[d$Aaa > 0 | d$Abb > 0 | d$Aab > 0 | d$Aba > 0] <- NaN
  d$fit[d$Aaa > 0 | d$Abb > 0 | d$Aab > 0 | d$Aba > 0] <- NaN
  nicheS <- 1 - sqrt((d$Cba[1] * d$Cab[1])/(d$Caa[1] * d$Cbb[1]))
  fitS <- sqrt((d$Caa[1] * d$Cab[1]) / (d$Cbb[1] * d$Cba[1]))
  f <- function(x) 1 / (-x + 1)
  p <- ggplot(d, aes(x = niche, y = fit, color = v)) + 
    geom_point(size = 2) + geom_abline(slope = -1, intercept = 1) + 
    theme_classic(base_size = 12) + stat_function(fun=f) + xlim(-1, 1) + ylim(0, 2) + 
    xlab("Niche Difference") + ylab("Fitness Ratio") +
    scale_color_gradient(low = "#F5D544", high = "#844A9D")
  if(sterile == TRUE) p + geom_point(aes(x = nicheS, y = fitS), color = "#2F242C") else p
}

### Species Specific Pathogen -----

## Panel A: Conceptual Diagram for Species Specific Pathogen

## Panel B: Coexistence Plot for Species Specific Pathogen

Caa <- rep(-.5, 201)
Cbb <- rep(-.8, 201)
Cba <- rep(-.8, 201)
Cab <- rep(-.5, 201)

s <- rep(-.2, 201)
p <- rep(.2, 201)

v <- seq(-10, 10, .1)

dat <- cbind.data.frame(Caa, Cbb, Cba, Cab, s, p, v)

dat$p.scaled <- dat$p + dat$v

# Visualizing the scaling function

ggplot(data = dat, aes( x= v, y = p.scaled)) + geom_line() + theme_classic(base_size = 12) + 
  xlab("Rainfall") + ylab("Phi_A")

dat$Aaa <- dat$Caa + dat$s*dat$p.scaled
dat$Abb <- dat$Cbb
dat$Aba <- dat$Cba
dat$Aab <- dat$Cab

coexist_plot(dat, TRUE)

### Adding a Species Specific Mutualist -----

## Panel C: Conceptual Diagram for Pathogen + Species Specific Mutualist

## Panel D: Coexistence Plot for Pathogen + Species Specific Mutualist

Caa <- rep(-.5, 201)
Cbb <- rep(-.8, 201)
Cba <- rep(-.8, 201)
Cab <- rep(-.5, 201)

sa <- rep(-.2, 201)
pa <- rep(.2, 201)
sb <- rep(.2, 201)
pb <- rep(.2, 201)

v <- seq(-10, 10, .1)

dat <- cbind.data.frame(Caa, Cbb, Cba, Cab, sa, pa, sb, pb, v)

dat$pa.scaled <- dat$pa + dat$v
dat$pb.scaled <- dat$pb - dat$v

# Visualizing the scaling functions

ggplot(data = dat) + geom_line(aes(x = v, y = pa.scaled), color = "red") + 
  geom_line(aes(x = v, y = pb.scaled), color = "black") +
  theme_classic(base_size = 12) + xlab("Rainfall") + ylab("Phi")

dat$Aaa <- dat$Caa + dat$sa*dat$pa.scaled
dat$Abb <- dat$Cbb + dat$sb*dat$pb.scaled
dat$Aba <- dat$Cba
dat$Aab <- dat$Cab

coexist_plot(dat, TRUE)

### Adding a Generalist Decomposer

## Panel E: Conceptual Diagram for Pathogen + Mutualist + General Decomposer

## Panel F: Coexistence Plot for Pathogen + Mutualist + General Decomposer

Caa <- rep(-.5, 201)
Cbb <- rep(-.8, 201)
Cba <- rep(-.8, 201)
Cab <- rep(-.5, 201)

sP <- rep(-.2, 201)
pP <- rep(.2, 201)
sM <- rep(.2, 201)
pM <- rep(.2, 201)
sDa <- rep(.3, 201)
pDa <- rep(.2, 201)
sDb <- rep(.1, 201)
pDb <- rep(.2, 201)

v <- seq(-10, 10, .1)

dat <- cbind.data.frame(Caa, Cbb, Cba, Cab, sP, pP, sM, pM, sDa, pDa, sDb, pDb, v)

dat$pP.scaled <- dat$pP + dat$v
dat$pM.scaled <- dat$pM - .5*dat$v
dat$pDa.scaled <- dat$pDa + dat$v
dat$pDb.scaled <- dat$pDb + dat$v

# Visualizing the scaling functions

ggplot(data = dat) + geom_line(aes(x = v, y = pP.scaled+.1), color = "red") + 
  geom_line(aes(x = v, y = pM.scaled), color = "black") +
  geom_line(aes(x = v, y = pDa.scaled), color = "green") +
  geom_line(aes(x = v, y = pDb.scaled-.1), color = "blue") +
  theme_classic(base_size = 12) + xlab("Rainfall") + ylab("Phi")

dat$Aaa <- dat$Caa + dat$sP*dat$pP.scaled + dat$sDa*dat$pDa.scaled
dat$Abb <- dat$Cbb + dat$sM*dat$pM.scaled + dat$sDb*dat$pDb.scaled
dat$Aba <- dat$Cba + dat$sDb*dat$pDa.scaled
dat$Aab <- dat$Cab + dat$sDa*dat$pDb.scaled

coexist_plot(dat, TRUE)
