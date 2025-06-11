
##### Local Sensitivity Analysis on Multiple Structures #####

library(tidyverse)
library(scales)
library(patchwork)

# A function to make some plots

local_sensitivity <- function(ND_mod, FR_mod, 
                              params) {
  
  env <- new.env()
  list2env(params, env)
  
  # 1) Compute the gradient
  ND_gradient <- deriv(expr = ND_mod, namevec = names(params))
  FR_gradient <- deriv(expr = FR_mod, namevec = names(params))
  
  # 2) Evaluate the gradient
  ND_eval <- eval(ND_gradient, envir = env)
  FR_eval <- eval(FR_gradient, envir = env)
  
  ND_values <- attr(ND_eval, "gradient")
  FR_values <- attr(FR_eval, "gradient")
  
  # 3) Prepare the data for plotting
  sens_df <- data.frame(
    pars = c(colnames(ND_values), colnames(FR_values)),
    vals = c(as.numeric(ND_values), as.numeric(FR_values)),
    metric = rep(c("ND", "FR"), c(length(ND_values), length(FR_values)))
  )
  
  # 4) Plot
  sensitivity_plot <- ggplot(sens_df, aes(x = pars, y = vals, fill = metric)) +
    geom_bar(stat = "identity", position = "dodge") +
    theme_classic(base_size = 12) +
    scale_y_continuous(trans= pseudo_log_trans()) + 
    labs(x = "Parameters", y = "Sensitivity") + 
    theme(axis.text.x = element_text(angle = 270, hjust = 0, vjust = 0)) + 
    scale_fill_manual(name = "Metric", labels = c("Fitness Ratio", "Niche Differences"), 
                      values = c("#F5D544", "#844A9D"))
  
  return(sensitivity_plot)
}

### Define Parameters -----

c_ii = -0.06; c_ij = -0.055; c_ji = -0.05; c_jj = -0.075

phi_xi = 10.5; phi_yi = 10.6; phi_zi = 10.7
phi_xj = 10; phi_yj = 10.1; phi_zj = 10.2

sig_ix = -0.002; sig_iy = -0.0021; sig_iz = -0.0022
sig_jx = -0.00225; sig_jy = -0.00235; sig_jz = -0.00245

### Compare Microbial Species Richness -----

##### Single Species -----

####### Sigma & Phi -----

sp1_ND <- expression(1 - sqrt(((c_ij + phi_xj * sig_ix)*
                                 (c_ji + phi_xi * sig_jx))/
                                ((c_ii + phi_xi * sig_ix)*
                                   (c_jj + phi_xj * sig_jx))))
sp1_FR <- expression(sqrt(((c_ij + phi_xj * sig_ix)*
                             (c_ii + phi_xi * sig_ix))/
                                ((c_ji + phi_xi * sig_jx)*
                                   (c_jj + phi_xj * sig_jx))))

local_sensitivity(ND_mod = sp1_ND, FR_mod = sp1_FR,
                  params = list(c_ii = c_ii, c_ij = c_ij, 
                                c_ji = c_ji, c_jj = c_jj, 
                                phi_xi = phi_xi, phi_xj = phi_xj, 
                                sig_ix = sig_ix, sig_jx = sig_jx)) + 
  theme(legend.position = "none")

ggsave(file.path("Figures", "Richness1.pdf"), 
       width = 6, height = 5, units = "in")

local_sensitivity(ND_mod = sp1_ND, FR_mod = sp1_FR,
                  params = list(c_ii = c_ii, c_ij = c_ij, 
                                c_ji = c_ji, c_jj = c_jj, 
                                phi_xi = phi_xi, phi_xj = phi_xj, 
                                sig_ix = sig_ix, sig_jx = sig_jx))

ggsave(file.path("Figures", "Richness_legend.pdf"), 
       width = 6, height = 5, units = "in")

####### M -----

m_ND <- expression(1 - sqrt(((c_ij + m_ij)*
                                 (c_ji + m_ji))/
                                ((c_ii + m_ii)*
                                   (c_jj + m_jj))))

m_FR <- expression(sqrt(((c_ij + m_ij)*
                             (c_ii + m_ii))/
                            ((c_ji + m_ji)*
                               (c_jj + m_jj))))

sp1m_plot <- local_sensitivity(ND_mod = m_ND, FR_mod = m_FR,
                              params = list(c_ii = c_ii, c_ij = c_ij, 
                                            c_ji = c_ji, c_jj = c_jj,
                                            m_ii = phi_xi * sig_ix, 
                                            m_ij = phi_xj * sig_ix, 
                                            m_ji = phi_xi * sig_jx, 
                                            m_jj = phi_xj * sig_jx))

##### Two Species -----

####### Sigma & Phi -----

sp2_ND <- expression(1 - sqrt(((c_ij + phi_xj * sig_ix + 
                                  phi_yj * sig_iy)*
                                 (c_ji + phi_xi * sig_jx + 
                                    phi_yi * sig_jy))/
                                ((c_ii + phi_xi * sig_ix + 
                                    phi_yi * sig_iy)*
                                   (c_jj + phi_xj * sig_jx + 
                                      phi_yj * sig_jy))))
sp2_FR <- expression(sqrt(((c_ij + phi_xj * sig_ix + 
                              phi_yj * sig_iy)*
                             (c_ii + phi_xi * sig_ix + 
                                phi_yi * sig_iy))/
                            ((c_ji + phi_xi * sig_jx + 
                                phi_yi * sig_jy)*
                               (c_jj + phi_xj * sig_jx + 
                                  phi_yj * sig_jy))))

local_sensitivity(ND_mod = sp2_ND, FR_mod = sp2_FR,
                              params = list(c_ii = c_ii, c_ij = c_ij, 
                                            c_ji = c_ji, c_jj = c_jj, 
                                            phi_xi = phi_xi, phi_xj = phi_xj,
                                            phi_yi = phi_yi, phi_yj = phi_yj,
                                            sig_ix = sig_ix, sig_jx = sig_jx, 
                                            sig_iy = sig_iy, sig_jy = sig_jy)) + 
  theme(legend.position = "none")

ggsave(file.path("Figures", "Richness2.pdf"), 
       width = 6, height = 5, units = "in")

####### M -----

sp2m_plot_unscaled <- local_sensitivity(ND_mod = m_ND, FR_mod = m_FR,
                                       params = list(c_ii = c_ii, c_ij = c_ij, 
                                                     c_ji = c_ji, c_jj = c_jj, 
                                                     m_ii = (phi_xi * sig_ix) +
                                                       (phi_yi * sig_iy), 
                                                     m_ij = (phi_xj * sig_ix) +
                                                       (phi_yj * sig_iy), 
                                                     m_ji = (phi_xi * sig_jx) +
                                                       (phi_yi * sig_jy), 
                                                     m_jj = (phi_xj * sig_jx) +
                                                       (phi_yj * sig_jy)))

sp2m_plot <- local_sensitivity(ND_mod = m_ND, FR_mod = m_FR,
                               params = list(c_ii = c_ii, c_ij = c_ij, 
                                             c_ji = c_ji, c_jj = c_jj, 
                                             m_ii = (phi_xi/2 * sig_ix) +
                                               (phi_yi/2 * sig_iy), 
                                             m_ij = (phi_xj/2 * sig_ix) +
                                               (phi_yj/2 * sig_iy), 
                                             m_ji = (phi_xi/2 * sig_jx) +
                                               (phi_yi/2 * sig_jy), 
                                             m_jj = (phi_xj/2 * sig_jx) +
                                               (phi_yj/2 * sig_jy)))

##### Three Species -----

####### Sigma & Phi -----

sp3_ND <- expression(1 - sqrt(((c_ij + phi_xj * sig_ix + 
                                  phi_yj * sig_iy + 
                                  phi_zj * sig_iz)*
                                 (c_ji + phi_xi * sig_jx + 
                                    phi_yi * sig_jy + 
                                    phi_zi * sig_jz))/
                                ((c_ii + phi_xi * sig_ix + 
                                    phi_yi * sig_iy + 
                                    phi_zi * sig_iz)*
                                   (c_jj + phi_xj * sig_jx + 
                                      phi_yj * sig_jy + 
                                      phi_zj * sig_jz))))
sp3_FR <- expression(sqrt(((c_ij + phi_xj * sig_ix + 
                              phi_yj * sig_iy + 
                              phi_zj * sig_iz)*
                             (c_ii + phi_xi * sig_ix + 
                                phi_yi * sig_iy + 
                                phi_zi * sig_iz))/
                            ((c_ji + phi_xi * sig_jx + 
                                phi_yi * sig_jy + 
                                phi_zi * sig_iz)*
                               (c_jj + phi_xj * sig_jx + 
                                  phi_yj * sig_jy + 
                                  phi_zj * sig_jz))))

local_sensitivity(ND_mod = sp3_ND, FR_mod = sp3_FR,
                  params = list(c_ii = c_ii, c_ij = c_ij, 
                                c_ji = c_ji, c_jj = c_jj, 
                                phi_xi = phi_xi, phi_xj = phi_xj,
                                phi_yi = phi_yi, phi_yj = phi_yj,
                                phi_zi = phi_zi, phi_zj = phi_zj,
                                sig_ix = sig_ix, sig_jx = sig_jx, 
                                sig_iy = sig_iy, sig_jy = sig_jy,
                                sig_iz = sig_iz, sig_jz = sig_jz)) + 
  theme(legend.position = "none")

ggsave(file.path("Figures", "Richness3.pdf"), 
       width = 6, height = 5, units = "in")

sp3_plot <- local_sensitivity(ND_mod = sp3_ND, FR_mod = sp3_FR,
                              params = list(c_ii = c_ii, c_ij = c_ij, 
                                            c_ji = c_ji, c_jj = c_jj, 
                                            phi_xi = phi_xi/3, phi_xj = phi_xj/3,
                                            phi_yi = phi_yi/3, phi_yj = phi_yj/3,
                                            phi_zi = phi_zi/3, phi_zj = phi_zj/3,
                                            sig_ix = sig_ix, sig_jx = sig_jx, 
                                            sig_iy = sig_iy, sig_jy = sig_jy,
                                            sig_iz = sig_iz, sig_jz = sig_jz))

ggsave(file.path("Figures", "Richness1.pdf"), 
       width = 6, height = 5, units = "in")
sp1_plot / sp2_plot_unscaled + plot_layout(guides = "collect")
sp1_plot / sp2_plot / sp3_plot + plot_layout(guides = "collect")
