##### SENSITIVITY ANALYSIS #####

library(ggplot2)
library(brms)

set.seed(19)
#set.seed(20)
#set.seed(21)
# setting the range for the parameter distributions

c_min <- -2
c_max <- 0
sig_min <- -.1
sig_max <- .1
phi_min <- 0
phi_max <- 20

# some empty vectors to throw values into
cAAs <- c()
cBAs <- c()
cBBs <- c()
cABs <- c()
sigAs <- c()
sigBs <- c()
phiAs <- c()
phiBs <- c()
alphaAAs <- c()
alphaBAs <- c()
alphaBBs <- c()
alphaABs <- c()
NDs <- c()
FIs <- c()

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
  
  NDs[i] <- ifelse(alphaABs[i] < 0 & alphaBAs[i] < 0 & alphaAAs[i] < 0 & alphaBBs[i] < 0,
                   1 - sqrt((alphaABs[i]*alphaBAs[i])/(alphaAAs[i]*alphaBBs[i])), NA)
  FIs[i] <- ifelse(alphaABs[i] < 0 & alphaBAs[i] < 0 & alphaAAs[i] < 0 & alphaBBs[i] < 0, ifelse(sqrt((alphaAAs[i]*alphaABs[i])/(alphaBBs[i]*alphaBAs[i])) > sqrt((alphaBBs[i]*alphaBAs[i])/(alphaAAs[i]*alphaABs[i])),
                   sqrt((alphaAAs[i]*alphaABs[i])/(alphaBBs[i]*alphaBAs[i])), sqrt((alphaBBs[i]*alphaBAs[i])/(alphaAAs[i]*alphaABs[i]))), NA)
  
}

dat <- cbind.data.frame(cAA = scale(abs(cAAs)), cBA = scale(abs(cBAs)), 
                        cBB = scale(abs(cBBs)), cAB = scale(abs(cABs)),
                        sigA = scale(sigAs), sigB = scale(sigBs), 
                        phiA = scale(phiAs), phiB = scale(phiBs), 
                        alphaAA = alphaAAs, alphaBA = alphaBAs, 
                        alphaBB = alphaBBs, alphaAB = alphaABs, 
                        ND = NDs, FI = FIs)

# adding a binary coexistence variable for logistic regression
dat$coexist <- ifelse(dat$FI > -dat$ND + 1 & dat$FI < 1/(-dat$ND+1), 1, 0)

# running some models
fit.FI <- brm(abs(FI-1) ~  cAA + cBA + cBB + cAB + sigA + sigB + phiA + phiB, data = dat, cores = 4)
fit.ND <- brm(ND ~ cAA + cBA + cBB + cAB + sigA + sigB + phiA + phiB, data = dat, cores = 4)
fit.coexist <- brm(coexist ~ cAA + cBA + cBB + cAB + sigA + sigB + phiA + phiB, data = dat, 
                   family = bernoulli, cores = 4)

# extracting coefficients and CIs
dat2 <- cbind.data.frame(estimate = c(unlist(summary(fit.FI)[["fixed"]]["Estimate"]), 
                                      unlist(summary(fit.ND)[["fixed"]]["Estimate"]), 
                                      unlist(summary(fit.coexist)[["fixed"]]["Estimate"])), 
                         metric = c(rep("FI", length(unlist(summary(fit.FI)[["fixed"]]["Estimate"]))), 
                                    rep("ND", length(unlist(summary(fit.ND)[["fixed"]]["Estimate"]))), 
                                    rep("CX", length(unlist(summary(fit.coexist)[["fixed"]]["Estimate"])))), 
                         lower = c(unlist(summary(fit.FI)[["fixed"]]["l-95% CI"]), 
                                   unlist(summary(fit.ND)[["fixed"]]["l-95% CI"]), 
                                   unlist(summary(fit.coexist)[["fixed"]]["l-95% CI"])),
                         upper = c(unlist(summary(fit.FI)[["fixed"]]["u-95% CI"]), 
                                   unlist(summary(fit.ND)[["fixed"]]["u-95% CI"]), 
                                   unlist(summary(fit.coexist)[["fixed"]]["u-95% CI"])), 
                         par = c(rownames(summary(fit.FI)[["fixed"]]),
                                 rownames(summary(fit.ND)[["fixed"]]),
                                 rownames(summary(fit.coexist)[["fixed"]])))

ggplot(data = dat2[which(dat2$par != "Intercept"),], aes(x = par, y = estimate, fill = metric, 
                                                         ymin = lower, ymax = upper)) + 
  geom_col(position = position_dodge()) + theme_classic() + 
  scale_fill_manual(name = "Metric", values = c( "black", "#F5D544", "#583B64"), 
                    labels = c("Coexistence", "Fitness Inequalities", "Niche Differences")) + 
  theme(axis.text.x = element_text(angle = 290, vjust = 0.5, hjust=0, size = 15)) + 
  scale_x_discrete(labels = c(bquote(C[AA]), bquote(C[AB]), bquote(C[BA]), bquote(C[BB]), 
                              bquote(phi[A]), bquote(phi[B]), bquote(sigma[A]), bquote(sigma[B]))) + 
  xlab("Parameter") + ylab("Effect Size") + 
  geom_errorbar(position = position_dodge(width = .9), width = .25)

ggplot(dat[1:500,], aes(x = cAA, y = ND)) + geom_point(alpha = .7) +
  geom_smooth(method = "lm") + ylim(-5, 5)
ggplot(dat[1:500,], aes(x = cAB, y = ND)) + geom_point(alpha = .7) +
  geom_smooth(method = "lm") + ylim(-5, 5)

ggplot(dat[1:500,], aes(x = cAA, y = abs(FI-1))) + geom_point(alpha = .7) +
  geom_smooth(method = "lm") + ylim(-5, 5)
ggplot(dat[1:500,], aes(x = cAB, y = abs(FI-1))) + geom_point(alpha = .7) +
  geom_smooth(method = "lm") + ylim(-5, 5)

