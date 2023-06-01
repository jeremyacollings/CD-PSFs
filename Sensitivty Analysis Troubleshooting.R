set.seed(10)

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

ggplot(dat[1:500,], aes(x = cAA, y = ND)) + geom_point(alpha = .7) +
  geom_smooth(method = "lm") + ylim(-5, 5)

ggplot(dat[1:500,], aes(x = cAB, y = ND)) + geom_point(alpha = .7) +
  geom_smooth(method = "lm") + ylim(-5, 5)

ggplot(dat[1:500,], aes(x = cAA, y = abs(FI-1))) + geom_point(alpha = .7) +
  geom_smooth(method = "lm") + ylim(-.5, 5)
ggplot(dat[1:500,], aes(x = cAB, y = abs(FI-1))) + geom_point(alpha = .7) +
  geom_smooth(method = "lm") + ylim(-.5, 5)

i = 23
fi <- c()
nd <- c()
delta <- c()
for(j in seq(-2, 2, .01)){
alphaAA <- (cAAs[i]-j)#+sigAs[i]*phiAs[i]
alphaBA <- (cBAs[i]-0)#+sigBs[i]*phiAs[i]
alphaBB <- (cBBs[i]-0)#+sigBs[i]*phiBs[i]
alphaAB <- (cABs[i]-0)#+sigAs[i]*phiBs[i]

delta <- c(delta, j)
fi <- c(fi, ifelse(sqrt((alphaAA*alphaAB)/(alphaBB*alphaBA)) > sqrt((alphaBB*alphaBA)/(alphaAA*alphaAB)),
                 sqrt((alphaAA*alphaAB)/(alphaBB*alphaBA)), sqrt((alphaBB*alphaBA)/(alphaAA*alphaAB))))
nd <- c(nd, 1 - sqrt((alphaAB*alphaBA)/(alphaAA*alphaBB)))
}
plot(delta, fi)
plot(delta, nd)

i = 23
fi <- c()
nd <- c()
delta <- c()
for(j in seq(-2, 2, .01)){
  alphaAA <- (cAAs[i]-0)#+sigAs[i]*phiAs[i]
  alphaBA <- (cBAs[i]-j)#+sigBs[i]*phiAs[i]
  alphaBB <- (cBBs[i]-0)#+sigBs[i]*phiBs[i]
  alphaAB <- (cABs[i]-0)#+sigAs[i]*phiBs[i]
  
  delta <- c(delta, j)
  fi <- c(fi, ifelse(sqrt((alphaAA*alphaAB)/(alphaBB*alphaBA)) > sqrt((alphaBB*alphaBA)/(alphaAA*alphaAB)),
                     sqrt((alphaAA*alphaAB)/(alphaBB*alphaBA)), sqrt((alphaBB*alphaBA)/(alphaAA*alphaAB))))
  nd <- c(nd, 1 - sqrt((alphaAB*alphaBA)/(alphaAA*alphaBB)))
}
plot(delta, fi)
plot(delta, nd)
