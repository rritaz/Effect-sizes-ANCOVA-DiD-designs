###-----------------------------------------------------------------###
###-----------------------------------------------------------------###
###-----------------------------------------------------------------###
### GAINS Effect size estimators 
###
### This file contains the functions one can use to compute the effect 
### size estimates and their variances 
###-----------------------------------------------------------------###
###-----------------------------------------------------------------###
###-----------------------------------------------------------------###

# Load libraries
library(tidyverse)
library(readxl)

# Read in GAINS data from Excel sheet 
# Note that this is a generic path. Make sure to change accordingly. 
GAINS_data <- read_excel("./Nutley.etal_data.xlsx")

ss_cor <- function(df){
    c <- 1-3/(4*df-1)
    return(c)
    
} # function to compute small sample correction factor with df degrees of freedom

###-----------------------------------------------------------------###
g <- function(mean_trt, mean_ctrl, n_trt, n_ctrl, post_sd_trt, post_sd_ctrl){
    diff_means = mean_trt - mean_ctrl
    s_y = sqrt(((n_trt-1)*post_sd_trt^2 + (n_ctrl-1)*post_sd_ctrl^2)/(n_trt+n_ctrl-2))
    g = ss_cor(n_trt + n_ctrl - 2)*(diff_means/s_y)
    return(g)
} # function to compute g

var_g <- function(g, n_trt, n_ctrl){
    n_tilda = (n_trt*n_ctrl)/(n_trt+n_ctrl)
    var_g = (1/n_tilda)+(g^2/(2*(n_trt+n_ctrl-2)))
    return(var_g)
} # function to compute var_g

g(mean_trt = GAINS_data$post_mean_trt, mean_ctrl = GAINS_data$post_mean_ctrl,
  n_trt = GAINS_data$n_trt, n_ctrl = GAINS_data$n_ctrl, 
  post_sd_trt = GAINS_data$post_sd_trt, post_sd_ctrl = GAINS_data$post_sd_ctrl)
# 1.560662

var_g(g = 1.560662, n_trt = GAINS_data$n_trt, n_ctrl = GAINS_data$n_ctrl)
# 0.107578

###-----------------------------------------------------------------###

gG1 <- function(gains_mean_trt, gains_mean_ctrl, n_trt, n_ctrl, post_sd_trt, post_sd_ctrl){
    diff_gain_means = gains_mean_trt - gains_mean_ctrl
    s_y = sqrt(((n_trt-1)*post_sd_trt^2 + (n_ctrl-1)*post_sd_ctrl^2)/(n_trt+n_ctrl-2))
    gG1 = ss_cor(n_trt + n_ctrl - 2)*(diff_gain_means/s_y)
    return(gG1)
} # function to compute g_G1

var_gG1 <- function(gG1, r, n_trt, n_ctrl){
    n_tilda = (n_trt*n_ctrl)/(n_trt+n_ctrl)
    var_gG1 = ((2*(1-r))/n_tilda)+((gG1^2)/(2*(n_trt+n_ctrl-2)))
    return(var_gG1)
} # function to compute var_gG1

gG1(gains_mean_trt = GAINS_data$gains_mean_trt, gains_mean_ctrl = GAINS_data$gains_mean_ctrl, 
    n_trt = GAINS_data$n_trt, n_ctrl = GAINS_data$n_ctrl, 
    post_sd_trt = GAINS_data$post_sd_trt, post_sd_ctrl = GAINS_data$post_sd_ctrl )
# 1.322228

var_gG1(gG1 = 1.322228, r = GAINS_data$r, n_trt = GAINS_data$n_trt , n_ctrl = GAINS_data$n_ctrl)
# 0.0734788

###-----------------------------------------------------------------###

gG2 <- function(gains_mean_trt, gains_mean_ctrl, n_trt, n_ctrl, gains_sd_trt, gains_sd_ctrl, r){
    diff_gain_means = gains_mean_trt - gains_mean_ctrl
    s_g = sqrt(((n_trt-1)*gains_sd_trt^2 + (n_ctrl-1)*gains_sd_ctrl^2)/(n_trt+n_ctrl-2))
    gG2 = ss_cor(n_trt + n_ctrl - 2)*(diff_gain_means/s_g)*sqrt(2*(1-r))
    return(gG2)
} # function to compute g_G2

var_gG2 <- function(gG2, r, n_trt, n_ctrl){
    n_tilda = (n_trt*n_ctrl)/(n_trt+n_ctrl)
    var_gG2 = ((2*(1-r))/n_tilda)+((gG2^2)/(2*(n_trt+n_ctrl-2)))
    return(var_gG2)
} # function to compute var_gG2

gG2(gains_mean_trt = GAINS_data$gains_mean_trt, gains_mean_ctrl = GAINS_data$gains_mean_ctrl, 
    n_trt = GAINS_data$n_trt, n_ctrl = GAINS_data$n_ctrl, 
    gains_sd_trt = GAINS_data$gains_mean_trt_sd, gains_sd_ctrl = GAINS_data$gains_mean_ctrl_sd, 
    r = GAINS_data$r)
# 1.443785

var_gG2(gG2 = 1.443785, r = GAINS_data$r, n_trt = GAINS_data$n_trt , n_ctrl = GAINS_data$n_ctrl)
# 0.07705569

###-----------------------------------------------------------------###

gPG <- function(gains_mean_trt, gains_mean_ctrl, n_trt, n_ctrl, pre_sd_trt, pre_sd_ctrl, 
                post_sd_trt, post_sd_ctrl){
    diff_gain_means = gains_mean_trt - gains_mean_ctrl
    s_p = sqrt(((n_trt-1)*post_sd_trt^2 + (n_ctrl-1)*post_sd_ctrl^2+
                    (n_trt-1)*pre_sd_trt^2 + (n_ctrl-1)*pre_sd_ctrl^2)/(2*(n_trt+n_ctrl-2)))
    gG2 = ss_cor(2*(n_trt + n_ctrl - 2))*(diff_gain_means/s_p)
    return(gG2)
} # function to compute g_PG

var_gPG <- function(gPG, r, n_trt, n_ctrl){
    n_tilda = (n_trt*n_ctrl)/(n_trt+n_ctrl)
    var_gPG = ((2*(1-r))/n_tilda)+((gPG^2)/(4*(n_trt+n_ctrl-2)))
    return(var_gPG)
} # function to compute var_gPG

gPG(gains_mean_trt = GAINS_data$gains_mean_trt, gains_mean_ctrl = GAINS_data$gains_mean_ctrl, 
    n_trt = GAINS_data$n_trt, n_ctrl = GAINS_data$n_ctrl, 
    pre_sd_trt = GAINS_data$pre_sd_trt, pre_sd_ctrl = GAINS_data$pre_sd_ctrl, 
    post_sd_trt = GAINS_data$post_sd_trt, post_sd_ctrl = GAINS_data$post_sd_ctrl)
# 1.437849

var_gPG(gPG = 1.437849, r = GAINS_data$r, n_trt = GAINS_data$n_trt , n_ctrl = GAINS_data$n_ctrl)
# 0.06587686

###-----------------------------------------------------------------###


