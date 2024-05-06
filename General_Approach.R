##### GENERAL APPROACH #####

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

# Panel A contains a simple, example visual model of the system

# Panel B shows the scaling function for sigma

v <- seq(-10, 10, .1)
sig <- .0002 + .001*v
dat <- cbind.data.frame(sig, v)
dat2 <- data.frame(v = c(-7.5, 0, 7.5), sig = c(.0002 + .001*-7.5, 
                                                .0002 + .001*0, 
                                                .0002 + .001*7.5))
ggplot(data = dat, aes(x = v, y = sig)) + geom_smooth(color = "black") + 
  ylab(expression(sigma)) + xlab("Environmental Variable") + theme_classic() +
  geom_point(data = dat2, aes(x = v, y = sig), shape = 1, color = "#844A9D", size = 7)

# Panel C is a triptych plot showing the downstream effect of scaling sigma on alphas

v <- seq(-10, 10, .1)

Caa <- c(rep(-.01, 201))
Cab <- c(rep(-.008, 201))
Cbb <- c(rep(-.0125, 201))
Cba <- c(rep(-.009, 201))

Pa <- c(rep(.2, 201))
Pb <- c(rep(.2, 201))

Sa <- c(rep(.0002, 201))
Sb <- c(rep(.0002, 201))

Sa2 <- Sa + .001*v

Aaa <- Caa + Pa*Sa2
Aba <- Cba + Pa*Sb
Abb <- Cbb + Pb*Sb
Aab <- Cab + Pa*Sa2

dat <- cbind.data.frame(v, Caa, Cab, Cbb, Cba, Aaa, Aba, Abb, Aab)

dat2 <- data.frame(v = c(-7.5, 0, 7.5), 
                   Aaa = c(-.01 + .2*(.0002+.001*-7.5), 
                           -.01 + .2*(.0002+.001*0), 
                           -.01 + .2*(.0002+.001*7.5)), 
                   Aba = rep(-.009 + .2*.0002, 3), 
                   Abb = rep(-.0125 + .2*.0002, 3), 
                   Aab = c(-.008 + .2*(.0002+.001*-7.5), 
                           -.008 + .2*(.0002+.001*0), 
                           -.008 + .2*(.0002+.001*7.5)))
dat2$niche <- 1 - sqrt((dat2$Aba * dat2$Aab)/(dat2$Aaa * dat2$Abb))
dat2$fit <- sqrt((dat2$Aaa * dat2$Aab) / (dat2$Abb * dat2$Aba))

strokes <- dat2[, 2:5]*2 # numerical value for stroke sizes in alpha plot; strokes defined in illustrator

# Panel D contains the coexistence plot with the environmental gradient

coexist_plot(dat, sterile = FALSE) + 
  geom_point(data = dat2, aes(x = niche, y = fit), color = "#844A9D", size = 5, shape = 1) + 
  xlim(-0.5, 0.5) + ylim(0.5, 1.5)

