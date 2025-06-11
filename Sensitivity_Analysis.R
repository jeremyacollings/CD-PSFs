
set.seed(6)
library(tidyverse)

# Visual Sensitivity Analysis for Collings et al. 2025

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

### Phi & Sigma Parameterization -----

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
Sax <- rep(-0.002, 1001)
Sbx <- rep(-0.00225, 1001)
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

### m Parameterization -----

# get min and max m values
min_Maa <- min(outer(Sax2, Pxa2)); max_Maa <- max(outer(Sax2, Pxa2))
min_Mab <- min(outer(Sax2, Pxb2)); max_Mab <- max(outer(Sax2, Pxb2))
min_Mba <- min(outer(Sbx2, Pxb2)); max_Mba <- max(outer(Sbx2, Pxb2))
min_Mbb <- min(outer(Sbx2, Pxb2)); max_Mbb <- max(outer(Sbx2, Pxb2))

Caa <- rep(-0.06, 1001)
Cba <- rep(-0.05, 1001)
Cbb <- rep(-0.075, 1001)
Cab <- rep(-0.055, 1001)
Maa <- rep(-0.002 * 10.5, 1001)
Mba <- rep(-0.00225 * 10.5, 1001)
Mbb <- rep(-0.00225 * 10, 1001)
Mab <- rep(-0.002 * 10, 1001)

##### Niche Differences -----

Caa2 <- seq(-.1, 0, .0001)
Aaa <- Caa2 + Maa
Aba <- Cba + Mba
Abb <- Cbb + Mbb
Aab <- Cab + Mab
ND_Caa <- niche_diff(Aaa, Aba, Abb, Aab)
FI_Caa <- fitness_ineq(Aaa, Aba, Abb, Aab)
FR_Caa <- fitness_ratio(Aaa, Aba, Abb, Aab)

Cba2 <- seq(-.1, 0, .0001)
Aaa <- Caa + Maa
Aba <- Cba2 + Mba
Abb <- Cbb + Mbb
Aab <- Cab + Mab
ND_Cba <- niche_diff(Aaa, Aba, Abb, Aab)
FI_Cba <- fitness_ineq(Aaa, Aba, Abb, Aab)
FR_Cba <- fitness_ratio(Aaa, Aba, Abb, Aab)

Cbb2 <- seq(-.1, 0, .0001)
Aaa <- Caa + Maa
Aba <- Cba + Mba
Abb <- Cbb2 + Mbb
Aab <- Cab + Mab
ND_Cbb <- niche_diff(Aaa, Aba, Abb, Aab)
FI_Cbb <- fitness_ineq(Aaa, Aba, Abb, Aab)
FR_Cbb <- fitness_ratio(Aaa, Aba, Abb, Aab)

Cab2 <- seq(-.1, 0, .0001)
Aaa <- Caa + Maa
Aba <- Cba + Mba
Abb <- Cbb + Mbb
Aab <- Cab2 + Mab
ND_Cab <- niche_diff(Aaa, Aba, Abb, Aab)
FI_Cab <- fitness_ineq(Aaa, Aba, Abb, Aab)
FR_Cab <- fitness_ratio(Aaa, Aba, Abb, Aab)

Maa2 <- seq(min_Maa, max_Maa, .0004)
Aaa <- Caa + Maa2
Aba <- Cba + Mba
Abb <- Cbb + Mbb
Aab <- Cab + Mab
ND_Maa <- niche_diff(Aaa, Aba, Abb, Aab)
FI_Maa <- fitness_ineq(Aaa, Aba, Abb, Aab)
FR_Maa <- fitness_ratio(Aaa, Aba, Abb, Aab)

Mba2 <- seq(min_Mba, max_Mba, .0004)
Aaa <- Caa + Maa
Aba <- Cba + Mba2
Abb <- Cbb + Mbb
Aab <- Cab + Mab
ND_Mba <- niche_diff(Aaa, Aba, Abb, Aab)
FI_Mba <- fitness_ineq(Aaa, Aba, Abb, Aab)
FR_Mba <- fitness_ratio(Aaa, Aba, Abb, Aab)

Mbb2 <- seq(min_Mbb, max_Mbb, .0004)
Aaa <- Caa + Maa
Aba <- Cba + Mba
Abb <- Cbb + Mbb2
Aab <- Cab + Mab
ND_Mbb <- niche_diff(Aaa, Aba, Abb, Aab)
FI_Mbb <- fitness_ineq(Aaa, Aba, Abb, Aab)
FR_Mbb <- fitness_ratio(Aaa, Aba, Abb, Aab)

Mab2 <- seq(min_Mab, max_Mab, .0004)
Aaa <- Caa + Maa
Aba <- Cba + Mba
Abb <- Cbb + Mbb
Aab <- Cab + Mab2
ND_Mab <- niche_diff(Aaa, Aba, Abb, Aab)
FI_Mab <- fitness_ineq(Aaa, Aba, Abb, Aab)
FR_Mab <- fitness_ratio(Aaa, Aba, Abb, Aab)

sens_dat2 <- cbind.data.frame(NDs = c(ND_Caa, ND_Cba, 
                                     ND_Cbb, ND_Cab, 
                                     ND_Maa, ND_Mba, 
                                     ND_Mbb, ND_Mab),
                             FIs = c(FI_Caa, FI_Cba, 
                                     FI_Cbb, FI_Cab, 
                                     FI_Maa, FI_Mba, 
                                     FI_Mbb, FI_Mab),
                             FRs = c(FR_Caa, FR_Cba, 
                                     FR_Cbb, FR_Cab, 
                                     FR_Maa, FR_Mba, 
                                     FR_Mbb, FR_Mab), 
                             par = rep(c("Caa", "Cba", "Cbb", "Cab", 
                                         "Maa", "Mba", "Mbb", "Mab"), each = 1001), 
                             val = c(std_params(0, 0.1, abs(Caa2)), 
                                     std_params(0, 0.1, abs(Cba2)), 
                                     std_params(0, 0.1, abs(Cbb2)), 
                                     std_params(0, 0.1, abs(Cab2)), 
                                     std_params(min_Maa, max_Maa, Maa2), 
                                     std_params(min_Mba, max_Mba, Mba2), 
                                     std_params(min_Mbb, max_Mbb, Mbb2), 
                                     std_params(min_Mab, max_Mab, Mab2)), 
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


# putting this into the right order for the figure

sens_dat3 <- rbind(sens_dat, sens_dat2[which(sens_dat2$par %in% 
                                  c("Maa", "Mba", "Mbb", "Mab")),])

sens_dat3$plot_index <- ifelse(sens_dat3$par %in% c("Caa", "Cba", 
                                                    "Cbb", "Cab"), "comp", 
                               ifelse(sens_dat3$par %in% c("Maa", "Mba", 
                                                           "Mbb", "Mab"), "m",
                                      ifelse(sens_dat3$par %in%
                                               c("Sax", "Sbx", 
                                                 "Pxa", "Pxb"), "micro", NA)))

# get little dataframe of baseline values
baseline_sens <- sens_dat3 %>%
  filter((par == "Caa" & val == std_params(-0.1, 0, Caa[1])) | 
           (par == "Cab" & val == std_params(-0.1, 0, Cab[1])) | 
           (par == "Cba" & val == std_params(-0.1, 0, Cba[1])) | 
           (par == "Cbb" & val == std_params(-0.1, 0, Cbb[1])) | 
           (par == "Sax" & val == std_params(-0.01, 0.01, Sax[1])) | 
           (par == "Sbx" & val == std_params(-0.01, 0.01, Sbx[1])) | 
           (par == "Pxa" & val == std_params(0, 20, Pxa[1])) | 
           (par == "Pax" & val == std_params(0, 20, Pxb[1])) | 
           (par == "Maa" & val == std_params(min_Maa, max_Maa, Maa[1])) | 
           (par == "Mab" & val == std_params(min_Mab, max_Mab, Mab[1])) | 
           (par == "Mba" & val == std_params(min_Mba, max_Mba, Mba[1])) | 
           (par == "Mbb" & val == std_params(min_Mbb, max_Mbb, Mbb[1])))

cbind.data.frame(par = c("Caa"), 
                 val = std_params(-0.1, 0, Caa[1]), 
                 NDs = niche_diff(Aaa[1], Aab[1], Aba[1], Abb[2]))
# ND

sens_dat3 %>%
  filter(plot_index == "comp") %>%
  ggplot(aes(x = val, y = NDs, color = par)) +
  geom_hline(yintercept = niche_diff(Aaa = Caa[1] + Maa[1], 
                                     Aab = Cab[1] + Mab[1], 
                                     Aba = Cba[1] + Mba[1], 
                                     Abb = Cbb[1] + Mbb[1]), 
             size = 1.5) + 
  geom_line(size = 1.5) + ylim(-1, 1) + 
  theme_classic(base_size = 12) + ylab("Niche Differences") + 
  xlab("Parameter Value") + theme_classic(base_size = 15) + 
  theme(axis.text.x = element_blank()) + 
  scale_color_manual(name = "Parameter", 
                     values = c("#844A9D", "#B990CB", 
                                "#DDC7E5", "#462753")) + 
  theme(legend.position = "none")

ggsave(file.path("Figures", "VSA_ND_comp.pdf"), 
       width = 5, height = 4, units = "in")

sens_dat3 %>%
  filter(plot_index == "micro") %>%
  ggplot(aes(x = val, y = NDs, color = par, linetype = type)) + 
  geom_hline(yintercept = niche_diff(Aaa = Caa[1] + Maa[1], 
                                     Aab = Cab[1] + Mab[1], 
                                     Aba = Cba[1] + Mba[1], 
                                     Abb = Cbb[1] + Mbb[1]), 
             size = 1.5) + 
  geom_line(size = 1.5) + ylim(-1, 1) + 
  theme_classic(base_size = 12) + ylab("Niche Differences") + 
  xlab("Parameter Value") + theme_classic(base_size = 15) + 
  theme(axis.text.x = element_blank()) + 
  scale_color_manual(name = "Parameter", 
                     values = c("#F5D544", "#D5B00B", 
                                "#FF7F11", "#FFB370")) + 
  theme(legend.position = "none")

ggsave(file.path("Figures", "VSA_ND_micro1.pdf"), 
       width = 5, height = 4, units = "in")

sens_dat3 %>%
  filter(plot_index == "m") %>%
  ggplot(aes(x = val, y = NDs, color = par, linetype = type)) + 
  geom_hline(yintercept = niche_diff(Aaa = Caa[1] + Maa[1], 
                                     Aab = Cab[1] + Mab[1], 
                                     Aba = Cba[1] + Mba[1], 
                                     Abb = Cbb[1] + Mbb[1]), 
             size = 1.5) + 
  geom_line(size = 1.5) + ylim(-1, 1) + 
  theme_classic(base_size = 12) + ylab("Niche Differences") + 
  xlab("Parameter Value") + theme_classic(base_size = 15) + 
  theme(axis.text.x = element_blank()) + 
  scale_color_manual(name = "Parameter", 
                     values = c("#3B747D","#58A4B0", 
                                "#9ECAD1", "#214045")) + 
  theme(legend.position = "none")

ggsave(file.path("Figures", "VSA_ND_micro2.pdf"), 
       width = 5, height = 4, units = "in")

# FR

sens_dat3 %>%
  filter(plot_index == "comp") %>%
  ggplot(aes(x = val, y = FRs, color = par, linetype = type)) + 
  geom_hline(yintercept = fitness_ratio(Aaa = Caa[1] + Maa[1], 
                                     Aab = Cab[1] + Mab[1], 
                                     Aba = Cba[1] + Mba[1], 
                                     Abb = Cbb[1] + Mbb[1]), 
             size = 1.5) + 
  geom_line(size = 1.5) + ylim(0, 2) + 
  theme_classic(base_size = 12) + ylab("Fitness Ratios") + 
  xlab("Parameter Value") + theme_classic(base_size = 15) + 
  theme(axis.text.x = element_blank()) + 
  scale_color_manual(name = "Parameter", 
                     values = c("#844A9D", "#B990CB", 
                                "#DDC7E5", "#462753")) + 
  theme(legend.position = "none")

ggsave(file.path("Figures", "VSA_FR_comp.pdf"), 
       width = 5, height = 4, units = "in")

sens_dat3 %>%
  filter(plot_index == "micro") %>%
  ggplot(aes(x = val, y = FRs, color = par, linetype = type)) + 
  geom_hline(yintercept = fitness_ratio(Aaa = Caa[1] + Maa[1], 
                                     Aab = Cab[1] + Mab[1], 
                                     Aba = Cba[1] + Mba[1], 
                                     Abb = Cbb[1] + Mbb[1]), 
             size = 1.5) + 
  geom_line(size = 1.5) + ylim(0, 2) + 
  theme_classic(base_size = 12) + ylab("Fitness Ratios") + 
  xlab("Parameter Value") + theme_classic(base_size = 15) + 
  theme(axis.text.x = element_blank()) + 
  scale_color_manual(name = "Parameter", 
                     values = c("#F5D544", "#D5B00B", 
                                "#FF7F11", "#FFB370")) + 
  theme(legend.position = "none")

ggsave(file.path("Figures", "VSA_FR_micro1.pdf"), 
       width = 5, height = 4, units = "in")

sens_dat3 %>%
  filter(plot_index == "m") %>%
  ggplot(aes(x = val, y = FRs, color = par, linetype = type)) + 
  geom_hline(yintercept = fitness_ratio(Aaa = Caa[1] + Maa[1], 
                                     Aab = Cab[1] + Mab[1], 
                                     Aba = Cba[1] + Mba[1], 
                                     Abb = Cbb[1] + Mbb[1]), 
             size = 1.5) + 
  geom_line(size = 1.5) + ylim(0, 2) + 
  theme_classic(base_size = 12) + ylab("Fitness Ratios") + 
  xlab("Parameter Value") + theme_classic(base_size = 15) + 
  theme(axis.text.x = element_blank()) + 
  scale_color_manual(name = "Parameter", 
                     values = c("#3B747D","#58A4B0", 
                                "#9ECAD1", "#214045")) + 
  theme(legend.position = "none")

ggsave(file.path("Figures", "VSA_FR_micro2.pdf"), 
       width = 5, height = 4, units = "in")


# FI (Appendix S3)

sens_dat3 %>%
  filter(plot_index == "comp") %>%
  ggplot(aes(x = val, y = FIs, color = par, linetype = type)) + 
  geom_line(size = 1.5) + ylim(1, 2) + 
  theme_classic(base_size = 12) + ylab("Fitness Inequalities") + 
  xlab("Parameter Value") + theme_classic(base_size = 15) + 
  theme(axis.text.x = element_blank()) + 
  scale_color_manual(name = "Parameter", 
                     values = c("#844A9D", "#B990CB", 
                                "#DDC7E5", "#462753")) + 
  theme(legend.position = "none")

ggsave(file.path("Figures", "VSA_FI_comp.pdf"), 
       width = 5, height = 4, units = "in")

sens_dat3 %>%
  filter(plot_index == "micro") %>%
  ggplot(aes(x = val, y = FIs, color = par, linetype = type)) + 
  geom_line(size = 1.5) + ylim(1, 2) + 
  theme_classic(base_size = 12) + ylab("Fitness Inequalities") + 
  xlab("Parameter Value") + theme_classic(base_size = 15) + 
  theme(axis.text.x = element_blank()) + 
  scale_color_manual(name = "Parameter", 
                     values = c("#F5D544", "#D5B00B", 
                                "#FF7F11", "#FFB370")) + 
  theme(legend.position = "none")

ggsave(file.path("Figures", "VSA_FI_micro1.pdf"), 
       width = 5, height = 4, units = "in")

sens_dat3 %>%
  filter(plot_index == "m") %>%
  ggplot(aes(x = val, y = FIs, color = par, linetype = type)) + 
  geom_line(size = 1.5) + ylim(1, 2) + 
  theme_classic(base_size = 12) + ylab("Fitness Inequalities") + 
  xlab("Parameter Value") + theme_classic(base_size = 15) + 
  theme(axis.text.x = element_blank()) + 
  scale_color_manual(name = "Parameter", 
                     values = c("#3B747D","#58A4B0", 
                                "#9ECAD1", "#214045")) + 
  theme(legend.position = "none")

ggsave(file.path("Figures", "VSA_FI_micro2.pdf"), 
       width = 5, height = 4, units = "in")

sens_dat3 %>%
  mutate(par = factor(par, levels = c("Caa", "Cab", "Cba", "Cbb", 
                                      "Pxa", "Pxb", "Sax", "Sbx", 
                                      "Maa", "Mab", "Mba", "Mbb"))) %>%
  ggplot(aes(x = val, y = FIs, color = par, linetype = type)) + 
  geom_line(size = 1.5) + ylim(1, 2) + 
  theme_classic(base_size = 12) + ylab("Fitness Inequalities") + 
  xlab("Parameter Value") + theme_classic(base_size = 15) + 
  theme(axis.text.x = element_blank()) + 
  scale_color_manual(name = "Parameter", 
                     values = c("#844A9D", "#B990CB", 
                                "#DDC7E5", "#462753",
                                "#F5D544", "#D5B00B", 
                                "#FF7F11", "#FFB370",
                                "#3B747D","#58A4B0", 
                                "#9ECAD1", "#214045"))

ggsave(file.path("Figures", "VSA_legend.pdf"), 
       width = 5, height = 4, units = "in")
