##### NONLINEAR SCALING FUNCTIONS #####

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
    xlab("Niche Difference") + ylab("Fitness Inequality") +
    scale_color_gradient(low = "#E8CCF5", high = "#260A33")
  if(sterile == TRUE) p + geom_point(aes(x = nicheS, y = fitS), color = "#2F242C") else p
}

# nonlinear functions... 
# goal: a six-panel figure, a row for phi and a row for sigma
# columns are linear, sigmoidal, and parabolic scaling

# Sigmoid Function: Gompertz Equation
# y = ae^(-be^(-cx))
# a = upper limit (asymptope)
# b = x-axis displacement
# c = growth rate

# Parbolic Function
# y = ax^2 + bx + c
# vertex = (-b/(2a), (4ac - b^2)/(4a))
# focal length = 1/(4a)

# Panel A: Linear Scaling of Phi

# phi = .2 + v

Caa <- rep(-.5, 201)
Cbb <- rep(-.8, 201)
Cba <- rep(-.8, 201)
Cab <- rep(-.5, 201)

s <- rep(-.2, 201)
p <- rep(.2, 201)

v <- seq(-10, 10, .1)

dat <- cbind.data.frame(Caa, Cbb, Cba, Cab, s, p, v)

dat$p.scaled <- dat$p + dat$v

ggplot(data = dat, aes( x= v, y = p.scaled)) + geom_line() + theme_classic(base_size = 12) + 
  xlab("Rainfall") + ylab("Phi_A")

dat$Aaa <- dat$Caa + dat$s*dat$p.scaled
dat$Abb <- dat$Cbb
dat$Aba <- dat$Cba
dat$Aab <- dat$Cab

coexist_plot(dat, TRUE)

# try this constraining phi to be greater than or equal to 0

dat <- dat[which(dat$p.scaled >= 0), ]

ggplot(data = dat, aes( x= v, y = p.scaled)) + geom_line() + theme_classic(base_size = 12) + 
  xlab("Rainfall") + ylab("Phi_A")

dat$Aaa <- dat$Caa + dat$s*dat$p.scaled
dat$Abb <- dat$Cbb
dat$Aba <- dat$Cba
dat$Aab <- dat$Cab

coexist_plot(dat, TRUE)

# Panel B: Linear Scaling of Sigma

# sigma = -.2 - v

Caa <- rep(-.5, 201)
Cbb <- rep(-.8, 201)
Cba <- rep(-.8, 201)
Cab <- rep(-.5, 201)

s <- rep(-.2, 201)
p <- rep(.2, 201)

v <- seq(-10, 10, .1)

dat <- cbind.data.frame(Caa, Cbb, Cba, Cab, s, p, v)

dat$s.scaled <- dat$s - dat$v

ggplot(data = dat, aes( x= v, y = s.scaled)) + geom_line() + theme_classic(base_size = 12) + 
  xlab("Rainfall") + ylab("Phi_A")

dat$Aaa <- dat$Caa + dat$p*dat$s.scaled
dat$Abb <- dat$Cbb
dat$Aba <- dat$Cba
dat$Aab <- dat$Cab

coexist_plot(dat, TRUE)

# Panel C: Sigmoidal Scaling of Phi

# not constraining:
# phi = ((20)/(1-exp(-v*.5))) - 9.8 
# constraining:
# phi = 10.2/(1+exp(5.1-v))

Caa <- rep(-.5, 201)
Cbb <- rep(-.8, 201)
Cba <- rep(-.8, 201)
Cab <- rep(-.5, 201)

s <- rep(-.2, 201)
p <- rep(.2, 201)

v <- seq(-10, 10, .1)

dat <- cbind.data.frame(Caa, Cbb, Cba, Cab, s, p, v)

dat$p.scaled <- (20/(1+exp(-v*.5))) - 9.8

ggplot(data = dat, aes( x= v, y = p.scaled)) + geom_line() + theme_classic(base_size = 12) + 
  xlab("Rainfall") + ylab("Phi_A")

dat$Aaa <- dat$Caa + dat$s*dat$p.scaled
dat$Abb <- dat$Cbb
dat$Aba <- dat$Cba
dat$Aab <- dat$Cab

coexist_plot(dat, TRUE)

# trying without - phi values

dat <- dat[which(dat$p.scaled >= 0), ]

ggplot(data = dat, aes( x= v, y = p.scaled)) + geom_line() + theme_classic(base_size = 12) + 
  xlab("Rainfall") + ylab("Phi_A")

dat$Aaa <- dat$Caa + dat$s*dat$p.scaled
dat$Abb <- dat$Cbb
dat$Aba <- dat$Cba
dat$Aab <- dat$Cab

coexist_plot(dat, TRUE)

# Panel D: Sigmoidal Scaling of Sigma

# sigma = ((20)/(1+exp(v*.5))) - 10.2

Caa <- rep(-.5, 201)
Cbb <- rep(-.8, 201)
Cba <- rep(-.8, 201)
Cab <- rep(-.5, 201)

s <- rep(-.2, 201)
p <- rep(.2, 201)

v <- seq(-10, 10, .1)

dat <- cbind.data.frame(Caa, Cbb, Cba, Cab, s, p, v)

dat$s.scaled <- (20/(1+exp(v*.5))) - 10.8

ggplot(data = dat, aes( x= v, y = s.scaled)) + geom_line() + theme_classic(base_size = 12) + 
  xlab("Rainfall") + ylab("Phi_A")

dat$Aaa <- dat$Caa + dat$p*dat$s.scaled
dat$Abb <- dat$Cbb
dat$Aba <- dat$Cba
dat$Aab <- dat$Cab

coexist_plot(dat, TRUE)

# Panel E: Parabolic Scaling of Phi

# phi = -.1(v-10)^2+10.2

Caa <- rep(-.5, 201)
Cbb <- rep(-.8, 201)
Cba <- rep(-.8, 201)
Cab <- rep(-.5, 201)

s <- rep(-.2, 201)
p <- rep(.2, 201)

v <- seq(-10, 10, .1)

dat <- cbind.data.frame(Caa, Cbb, Cba, Cab, s, p, v)

dat$p.scaled <- -.1*(v-10)^2+10.2

ggplot(data = dat, aes( x= v, y = p.scaled)) + geom_line() + theme_classic(base_size = 12) + 
  xlab("Rainfall") + ylab("Phi_A")

dat$Aaa <- dat$Caa + dat$s*dat$p.scaled
dat$Abb <- dat$Cbb
dat$Aba <- dat$Cba
dat$Aab <- dat$Cab

coexist_plot(dat, TRUE)

# only + phi

dat <- dat[which(dat$p.scaled >= 0), ]

ggplot(data = dat, aes( x= v, y = p.scaled)) + geom_line() + theme_classic(base_size = 12) + 
  xlab("Rainfall") + ylab("Phi_A")

dat$Aaa <- dat$Caa + dat$s*dat$p.scaled
dat$Abb <- dat$Cbb
dat$Aba <- dat$Cba
dat$Aab <- dat$Cab

coexist_plot(dat, TRUE)

# Panel F: Parabolic Scaling of Sigma

# sigma = .1(v-10)^2-10.2

Caa <- rep(-.5, 201)
Cbb <- rep(-.8, 201)
Cba <- rep(-.8, 201)
Cab <- rep(-.5, 201)

s <- rep(-.2, 201)
p <- rep(.2, 201)

v <- seq(-10, 10, .1)

dat <- cbind.data.frame(Caa, Cbb, Cba, Cab, s, p, v)

dat$s.scaled <- .1*(v-10)^2-10.2

ggplot(data = dat, aes( x= v, y = s.scaled)) + geom_line() + theme_classic(base_size = 12) + 
  xlab("Rainfall") + ylab("Phi_A")

dat$Aaa <- dat$Caa + dat$p*dat$s.scaled
dat$Abb <- dat$Cbb
dat$Aba <- dat$Cba
dat$Aab <- dat$Cab

coexist_plot(dat, TRUE)



# I wonder if this would look different if I was able to show the entire range... 
# let's have different baselines

# Linear Phi
Caa <- rep(-.5, 401)
Cbb <- rep(-.8, 401)
Cba <- rep(-.8, 401)
Cab <- rep(-.5, 401)

s <- rep(-.2, 401)
p <- rep(.2, 401)

v <- seq(-2, 2, .01)

dat <- cbind.data.frame(Caa, Cbb, Cba, Cab, s, p, v)

dat$p.scaled <- dat$p + dat$v

ggplot(data = dat, aes( x= v, y = p.scaled)) + geom_line() + theme_classic(base_size = 12) + 
  xlab("Rainfall") + ylab("Phi_A")

dat$Aaa <- dat$Caa + dat$s*dat$p.scaled
dat$Abb <- dat$Cbb
dat$Aba <- dat$Cba
dat$Aab <- dat$Cab

coexist_plot(dat, TRUE)

# try this constraining phi to be greater than or equal to 0

dat <- dat[which(dat$p.scaled >= 0), ]

ggplot(data = dat, aes( x= v, y = p.scaled, color = v)) + geom_line(size = 5) + theme_classic(base_size = 12) + 
  xlab("Rainfall") + ylab("Phi_A") + scale_color_gradient(low = "#E8CCF5", high = "#260A33")

dat$Aaa <- dat$Caa + dat$s*dat$p.scaled
dat$Abb <- dat$Cbb
dat$Aba <- dat$Cba
dat$Aab <- dat$Cab

coexist_plot(dat, TRUE)

#Linear Sigma
Caa <- rep(-.5, 401)
Cbb <- rep(-.8, 401)
Cba <- rep(-.8, 401)
Cab <- rep(-.5, 401)

s <- rep(-.2, 401)
p <- rep(.2, 401)

v <- seq(-2, 2, .01)

dat <- cbind.data.frame(Caa, Cbb, Cba, Cab, s, p, v)

dat$s.scaled <- dat$s - dat$v

ggplot(data = dat, aes( x= v, y = s.scaled, color = v)) + geom_line(size = 5) + theme_classic(base_size = 12) + 
  xlab("Rainfall") + ylab("Sig_A") + scale_color_gradient(low = "#E8CCF5", high = "#260A33")

dat$Aaa <- dat$Caa + dat$p*dat$s.scaled
dat$Abb <- dat$Cbb
dat$Aba <- dat$Cba
dat$Aab <- dat$Cab

coexist_plot(dat, TRUE)

#Sigmoidal Phi

Caa <- rep(-.5, 401)
Cbb <- rep(-.8, 401)
Cba <- rep(-.8, 401)
Cab <- rep(-.5, 401)

s <- rep(-.2, 401)
p <- rep(.2, 401)

v <- seq(-2, 2, .01)

dat <- cbind.data.frame(Caa, Cbb, Cba, Cab, s, p, v)

dat$p.scaled <- (2/(1+exp(-2*v)))

ggplot(data = dat, aes( x= v, y = p.scaled)) + geom_line() + theme_classic(base_size = 12) + 
  xlab("Rainfall") + ylab("Phi_A")

dat$Aaa <- dat$Caa + dat$s*dat$p.scaled
dat$Abb <- dat$Cbb
dat$Aba <- dat$Cba
dat$Aab <- dat$Cab

coexist_plot(dat, TRUE)

# trying without - phi values

dat <- dat[which(dat$p.scaled >= 0), ]

ggplot(data = dat, aes( x= v, y = p.scaled, color = v)) + geom_line(size = 5) + theme_classic(base_size = 12) + 
  xlab("Rainfall") + ylab("Phi_A") + scale_color_gradient(low = "#E8CCF5", high = "#260A33")

dat$Aaa <- dat$Caa + dat$s*dat$p.scaled
dat$Abb <- dat$Cbb
dat$Aba <- dat$Cba
dat$Aab <- dat$Cab

coexist_plot(dat, TRUE)

#Sigmoidal Sigma

Caa <- rep(-.5, 401)
Cbb <- rep(-.8, 401)
Cba <- rep(-.8, 401)
Cab <- rep(-.5, 401)

s <- rep(-.2, 401)
p <- rep(.2, 401)

v <- seq(-2, 2, .01)


dat <- cbind.data.frame(Caa, Cbb, Cba, Cab, s, p, v)

dat$s.scaled <- (4/(1+exp(v*2))) - 2.2

ggplot(data = dat, aes( x= v, y = s.scaled, color = v)) + geom_line(size = 5) + theme_classic(base_size = 12) + 
  xlab("Rainfall") + ylab("Sig_A") + scale_color_gradient(low = "#E8CCF5", high = "#260A33")

dat$Aaa <- dat$Caa + dat$p*dat$s.scaled
dat$Abb <- dat$Cbb
dat$Aba <- dat$Cba
dat$Aab <- dat$Cab

coexist_plot(dat, TRUE)

#Parabolic Phi

Caa <- rep(-.5, 401)
Cbb <- rep(-.8, 401)
Cba <- rep(-.8, 401)
Cab <- rep(-.5, 401)

s <- rep(-.2, 401)
p <- rep(.2, 401)

v <- seq(-2, 2, .01)

dat <- cbind.data.frame(Caa, Cbb, Cba, Cab, s, p, v)

dat$p.scaled <- -.5*(v-2)^2+2.2

ggplot(data = dat, aes( x= v, y = p.scaled)) + geom_line() + theme_classic(base_size = 12) + 
  xlab("Rainfall") + ylab("Phi_A")

dat$Aaa <- dat$Caa + dat$s*dat$p.scaled
dat$Abb <- dat$Cbb
dat$Aba <- dat$Cba
dat$Aab <- dat$Cab

coexist_plot(dat, TRUE)

dat <- dat[which(dat$p.scaled >= 0), ]

ggplot(data = dat, aes( x= v, y = p.scaled, color = v)) + geom_line(size = 5) + theme_classic(base_size = 12) + 
  xlab("Rainfall") + ylab("Phi_A") + scale_color_gradient(low = "#E8CCF5", high = "#260A33")

dat$Aaa <- dat$Caa + dat$s*dat$p.scaled
dat$Abb <- dat$Cbb
dat$Aba <- dat$Cba
dat$Aab <- dat$Cab

coexist_plot(dat, TRUE)

#Parabolic Sigma

Caa <- rep(-.5, 401)
Cbb <- rep(-.8, 401)
Cba <- rep(-.8, 401)
Cab <- rep(-.5, 401)

s <- rep(-.2, 401)
p <- rep(.2, 401)

v <- seq(-2, 2, .01)


dat <- cbind.data.frame(Caa, Cbb, Cba, Cab, s, p, v)

dat$s.scaled <- (1/4)*(v-2)^2-2.2

ggplot(data = dat, aes( x= v, y = s.scaled, color = v)) + geom_line(size = 5) + theme_classic(base_size = 12) + 
  xlab("Rainfall") + ylab("Sig_A") + scale_color_gradient(low = "#E8CCF5", high = "#260A33")

dat$Aaa <- dat$Caa + dat$p*dat$s.scaled
dat$Abb <- dat$Cbb
dat$Aba <- dat$Cba
dat$Aab <- dat$Cab

coexist_plot(dat, TRUE)

