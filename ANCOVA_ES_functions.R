###-----------------------------------------------------------------###
###-----------------------------------------------------------------###
###-----------------------------------------------------------------###
### ANCOVA Effect size estimators 
###
### This file contains the functions one can use to compute the effect 
### size estimates and their variances 
###-----------------------------------------------------------------###
###-----------------------------------------------------------------###
###-----------------------------------------------------------------###
library(tidyverse)
library(readxl)

# Read in excel data in long format
# Note this contains a generic path. Make sure to change accordingly. 
ANCOVA_data <- read_excel("./Sprenger.etal_data.xlsx")

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

# DEMONSTRATION WITH NUTLEY DATA

g(mean_trt = ANCOVA_data$post_mean_trt, mean_ctrl = ANCOVA_data$post_mean_ctrl,
  n_trt = ANCOVA_data$n_trt, n_ctrl = ANCOVA_data$n_ctrl, 
  post_sd_trt = ANCOVA_data$post_sd_trt, post_sd_ctrl = ANCOVA_data$post_sd_ctrl)
# 0.608554

var_g(g = 0.608554, n_trt = ANCOVA_data$n_trt, n_ctrl = ANCOVA_data$n_ctrl)
# 0.03709139

###-----------------------------------------------------------------###
gA1 <- function(mean_trt_adj, mean_ctrl_adj, n_trt, n_ctrl, post_sd_trt, post_sd_ctrl){
    diff_adj_means = mean_trt_adj - mean_ctrl_adj
    s_y = sqrt(((n_trt-1)*post_sd_trt^2 + (n_ctrl-1)*post_sd_ctrl^2)/(n_trt+n_ctrl-2))
    gA1 = ss_cor(n_trt + n_ctrl - 2)*(diff_adj_means/s_y)
    return(gA1)
} # function to compute g_A1

var_gA1 <- function(gA1, r, n_trt, n_ctrl){
    n_tilda = (n_trt*n_ctrl)/(n_trt+n_ctrl)
    var_gA1 = ((1-r^2)/n_tilda)+((gA1^2)/(2*(n_trt+n_ctrl-2)))
    return(var_gA1)
} # function to compute var_gA1

gA1(mean_trt_adj = ANCOVA_data$adj_mean_trt, mean_ctrl_adj = ANCOVA_data$adj_mean_ctrl, 
    n_trt = ANCOVA_data$n_trt, n_ctrl = ANCOVA_data$n_ctrl, 
    post_sd_trt = ANCOVA_data$post_sd_trt, post_sd_ctrl = ANCOVA_data$post_sd_ctrl)
# 0.4505042

var_gA1(gA1 = 0.4505042, r = ANCOVA_data$r, n_trt = ANCOVA_data$n_trt, n_ctrl = ANCOVA_data$n_ctrl)
# 0.01641186

###-----------------------------------------------------------------###

gA2 <- function(mean_trt_adj, mean_ctrl_adj, n_trt, n_ctrl, t_adj, r){
    n_tilda = (n_trt*n_ctrl)/(n_trt+n_ctrl)
    diff_adj_means = mean_trt_adj - mean_ctrl_adj
    gA2 = ss_cor(n_trt + n_ctrl - 4)*t_adj*sqrt((1-r^2)/n_tilda)
    return(gA2)
} # function to compute g_A2

var_gA2 <- function(gA2, r, n_trt, n_ctrl){
    n_tilda = (n_trt*n_ctrl)/(n_trt+n_ctrl)
    var_gA2 = ((1-r^2)/n_tilda)+(gA2^2/(2*(n_trt+n_ctrl-4)))
    return(var_gA2)
} # function to compute var_gA2

gA2(mean_trt_adj = ANCOVA_data$adj_mean_trt, mean_ctrl_adj = ANCOVA_data$adj_mean_ctrl, 
    n_trt = ANCOVA_data$n_trt, n_ctrl = ANCOVA_data$n_ctrl, 
    t_adj = ANCOVA_data$t, r = ANCOVA_data$r)
# 0.4450716

var_gA2(gA2 = 0.4450716, r = ANCOVA_data$r, n_trt = ANCOVA_data$n_trt, n_ctrl = ANCOVA_data$n_ctrl)
# 0.01640631

###-----------------------------------------------------------------###

gPA <- function(mean_trt_adj, mean_ctrl_adj, n_trt, n_ctrl, pre_sd_trt, pre_sd_ctrl, post_sd_trt, post_sd_ctrl){
    diff_adj_means = mean_trt_adj - mean_ctrl_adj
    s_p = sqrt(((n_trt-1)*post_sd_trt^2 + (n_ctrl-1)*post_sd_ctrl^2+
                    (n_trt-1)*pre_sd_trt^2 + (n_ctrl-1)*pre_sd_ctrl^2)/(2*n_trt+2*n_ctrl-4))
    gPA = ss_cor(2*(n_trt + n_ctrl - 2))*(diff_adj_means/s_p)
    return(gPA)
} # function to compute g_PA

var_gPA <- function(gPA, r, n_trt, n_ctrl){
    n_tilda = (n_trt*n_ctrl)/(n_trt+n_ctrl)
    var_gPA = ((1-r^2)/n_tilda)+(((1+r^2)*gPA^2)/(4*(n_trt+n_ctrl-2)))
    return(var_gPA)
} # function to compute var_gPA

gPA(ANCOVA_data$adj_mean_trt, ANCOVA_data$adj_mean_ctrl, 
    ANCOVA_data$n_trt, ANCOVA_data$n_ctrl, 
    ANCOVA_data$pre_sd_trt, ANCOVA_data$pre_sd_ctrl, 
    ANCOVA_data$post_sd_trt, ANCOVA_data$post_sd_ctrl)
# 0.4464666

var_gPA(gPA = 0.4464666, r = ANCOVA_data$r, n_trt = ANCOVA_data$n_trt, n_ctrl = ANCOVA_data$n_ctrl)
# 0.01619913

###-----------------------------------------------------------------###

