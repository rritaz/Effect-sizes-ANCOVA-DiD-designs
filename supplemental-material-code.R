###-----------------------------------------------------------------###
###-----------------------------------------------------------------###
###-----------------------------------------------------------------###
### Effect size estimators and their properties
###
### This file contains the code used to conduct the simulation study in this paper
###-----------------------------------------------------------------###
###-----------------------------------------------------------------###
###-----------------------------------------------------------------###
library(MASS)
library(tidyverse)

# Create empty data frames for gains and ANCOVA
df_rho_cumulative_gains <- data.frame()
df_rho_cumulative_ancova <- data.frame()

# Vector of sample sizes considered - n
n_vector <- c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100)

for (j in 1:length(n_vector)){    
    
    df_rho_gains <- data.frame()
    df_rho_ancova <- data.frame()
    
    for (i in 1:10000){
        ###-----------------------------------------------------------------###
        ### Define values
        rho <- 0.4
        theta <- 0.2
        q <- 1 # number of covariates
        iterations <- 10000 
        ###-----------------------------------------------------------------###
        #X & Y Treatment
        n_trt <- n_vector[j]
        mu_trt <-  c(0,0.3)
        sigma <- matrix(c(1,rho,rho,1), ncol=2)
        A <- mvrnorm(n=n_trt, mu_trt, sigma)
        gain_trt<-(A[,2]- A[,1])
        within_trt_data <- cbind(A, gain_trt) 
        colnames(within_trt_data) <- c("x_trt", "y_trt", "gain_trt")
        Ancova_sims <- A %>% as_tibble() %>%  mutate(treatment_indicator = rep(1, n_trt))
        colnames(Ancova_sims) <- c("x", "y", "treatment_indicator")
        
        ###-----------------------------------------------------------------###
        #X & Y Control
        n_ctrl <- n_vector[j]
        mu_ctrl <- c(0,0.1)
        sigma <- matrix(c(1,rho,rho,1), ncol=2)
        B <- mvrnorm(n=n_ctrl, mu_ctrl, sigma)
        gain_ctrl <- (B[,2]- B[,1])
        within_ctrl_data <- cbind(B, gain_ctrl)
        colnames(within_ctrl_data) <- c("x_ctrl", "y_ctrl", "gain_ctrl")
        Bncova_sims <- B %>% as_tibble() %>%  mutate(treatment_indicator = rep(0, n_ctrl))
        colnames(Bncova_sims) <- c("x", "y", "treatment_indicator")
        
        ###-----------------------------------------------------------------###
        # Store data
        trtANDctrl_df<- cbind(within_trt_data, within_ctrl_data) %>% as.data.frame() 
        ancova_df <- rbind(Ancova_sims,Bncova_sims) %>% 
            mutate(treatment_indicator = as.character(treatment_indicator),
                   x = scale(x, center = FALSE, scale = FALSE),
                   y = scale(y, center = FALSE, scale = FALSE))
        
        ###-----------------------------------------------------------------###
        # Run ANCOVA
        ancova_mod <- lm(ancova_df$y ~ as.factor(ancova_df$treatment_indicator) + ancova_df$x)
        anc_mod_summary <- summary(ancova_mod)
        sd_adjusted <- anc_mod_summary$sigma 
        diffMeansAdjusted <- anc_mod_summary$coefficients[2,1]
        df_complete_ancova <- cbind(ancova_df, sd_adjusted) 
        colnames(df_complete_ancova) <- c("x", "y", "treatment_indicator", "sd_adjusted")
        
        ###-----------------------------------------------------------------###
        ###-----------------------------------------------------------------###
        
        # long dataframe for GAINS
        df_gains <- trtANDctrl_df %>% 
            mutate(
                diff_means_post = mean(y_trt-y_ctrl),
                sd_post = sqrt((var(y_trt)+var(y_ctrl))/2),
                
                # g
                t_1 = (diff_means_post/sd_post),
                unbiased_t1 = (1-(3/(4*(n_trt+n_ctrl-2)-1)))*(diff_means_post/sd_post),
                v_1 = ((1/((n_trt*n_ctrl)/(n_trt+n_ctrl)))+(t_1^2/(2*(n_trt+n_ctrl-2)))),
                
                diff_gains = gain_trt - gain_ctrl,
                mean_diff_gains = mean(diff_gains),
                var_gain_trt = var(gain_trt),
                var_gain_ctrl = var(gain_ctrl),
                sd_pooled = sqrt(((n_trt - 1)*var_gain_trt+(n_ctrl-1)*var_gain_ctrl)/(n_trt+n_ctrl-2)),
                cov_trt = cov(x_trt, y_trt),
                cov_ctrl = cov(x_ctrl, y_ctrl),
                
                combo_LB_t1 = unbiased_t1-1.96*sqrt(v_1),
                combo_UB_t1 = unbiased_t1+1.96*sqrt(v_1),
                combo_in_CI_t1 = ifelse((theta >= combo_LB_t1 & theta <= combo_UB_t1), 1, 0),
                
                ###-----------------------------------------------------------------###
                corr_pooled_trt = cov_trt/(sqrt(var(x_trt))*sqrt(var(y_trt))),
                corr_pooled_ctrl = cov_ctrl/(sqrt(var(x_ctrl))*sqrt(var(y_ctrl))),
                corr_pooled = (corr_pooled_trt+corr_pooled_ctrl)/2,
                part1_in_estT5 = mean_diff_gains/sd_pooled,
                part2_in_estT5 = sqrt(2*(1-corr_pooled)),
                ###-----------------------------------------------------------------###
                
                # gG1
                t_2 = (mean_diff_gains/sd_post),
                unbiased_t2 = (1-(3/(4*(n_trt+n_ctrl-2)-1)))*(mean_diff_gains/sd_post),
                v_2 =  (((2*(1-rho))/((n_trt*n_ctrl)/(n_trt+n_ctrl)))+((t_2^2)/(2*(n_trt+n_ctrl-2)))),
                
                combo_LB_t2 = unbiased_t2-1.96*sqrt(v_2),
                combo_UB_t2 = unbiased_t2+1.96*sqrt(v_2),
                combo_in_CI_t2 = ifelse((theta >= combo_LB_t2 & theta <= combo_UB_t2), 1, 0),
                
                # gG2
                t_5 = (part1_in_estT5*part2_in_estT5),
                unbiased_t5 = (1-(3/(4*(n_trt+n_ctrl-2)-1)))*(part1_in_estT5*part2_in_estT5),
                v_5 = ((2*(1-corr_pooled))/((n_trt*n_ctrl)/(n_trt+n_ctrl)))+((t_5^2)/(2*(n_trt+n_ctrl-2))),

                combo_LB_t5 = unbiased_t5-1.96*sqrt(v_5),
                combo_UB_t5 = unbiased_t5+1.96*sqrt(v_5),
                combo_in_CI_t5 = ifelse((theta >= combo_LB_t5 & theta <= combo_UB_t5), 1, 0),
                
                # gG3
                t_5_trueRho = (part1_in_estT5*sqrt(2*(1-rho))),
                unbiased_t5_trueRho = (1-(3/(4*(n_trt+n_ctrl-2)-1)))*(part1_in_estT5*sqrt(2*(1-rho))),
                v_5_trueRho = (((2*(1-rho))/((n_trt*n_ctrl)/(n_trt+n_ctrl)))+((t_5_trueRho^2)/(2*(n_trt+n_ctrl-2)))),
               
                combo_LB_t5_trueRho = unbiased_t5_trueRho-1.96*sqrt(v_5_trueRho),
                combo_UB_t5_trueRho = unbiased_t5_trueRho+1.96*sqrt(v_5_trueRho),
                combo_in_CI_t5_trueRho = ifelse((theta >= combo_LB_t5_trueRho & theta <= combo_UB_t5_trueRho), 1, 0),
                
                ###-----------------------------------------------------------------###
                
                # Pooled across trt and ctrl across pre and post
                sPOOLED = sqrt((var(y_trt)+var(y_ctrl)+var(x_trt)+var(x_ctrl))/4),
                
                # gP
                t_10 = diff_means_post/sPOOLED,
                unbiased_t10 = (1-(3/(4*(2*(n_trt+n_ctrl-2))-1)))*diff_means_post/sPOOLED,
                v_10 = (1/((n_trt*n_ctrl)/(n_trt+n_ctrl)))+((1+rho^2)*(t_10^2)/(4*(n_trt+n_ctrl-2))),
               
                combo_LB_t10 = unbiased_t10-1.96*sqrt(v_10),
                combo_UB_t10 = unbiased_t10+1.96*sqrt(v_10),
                combo_in_CI_t10 = ifelse((theta >= combo_LB_t10 & theta <= combo_UB_t10), 1, 0),
                
                # gPG
                t_11 = mean_diff_gains/sPOOLED,
                unbiased_t11 = (1-(3/(4*(2*(n_trt+n_ctrl-2))-1)))*(mean_diff_gains/sPOOLED),
                v_11 = ((2*(1-rho))/((n_trt*n_ctrl)/(n_trt+n_ctrl)))+((1+rho^2)*(t_11^2)/(4*(n_trt+n_ctrl-2))),
                
                combo_LB_t11 = unbiased_t11-1.96*sqrt(v_11),
                combo_UB_t11 = unbiased_t11+1.96*sqrt(v_11),
                combo_in_CI_t11 = ifelse((theta >= combo_LB_t11 & theta <= combo_UB_t11), 1, 0),
            )
        
        # long dataframe for ANCOVA
        df_adjusted <- df_complete_ancova %>% 
            group_by(treatment_indicator) %>% 
            mutate(corr_sample = cor(x,y),
                   var_post_group = var(y),
                   var_pre_group = var(x)) %>% 
            ungroup() %>% 
            mutate(
                sd_post_pooled = sqrt((var_post_group[1]+var_post_group[length(var_post_group)])/2),
                
                ###-----------------------------------------------------------------###
                
                # gA1
                t_3 = diffMeansAdjusted/sd_post_pooled,
                unbiased_t3 = (1-(3/(4*(n_trt+n_ctrl-2-q)-1)))*diffMeansAdjusted/sd_post_pooled,
                v_3 = ((1-rho^2)/((n_trt*n_ctrl)/(n_trt+n_ctrl)))+((t_3^2)/(2*(n_trt+n_ctrl-2))),
               
                combo_LB_t3 = unbiased_t3-1.96*sqrt(v_3),
                combo_UB_t3 = unbiased_t3+1.96*sqrt(v_3),
                combo_in_CI_t3 = ifelse((theta >= combo_LB_t3 & theta <= combo_UB_t3), 1, 0),

                # gA2
                part1_t9 = diffMeansAdjusted/sd_adjusted[1],
                corr_pooled = (corr_sample[1]+corr_sample[length(corr_sample)])/2,
                part2_t9 = sqrt(1-corr_pooled^2),
                
                t_9 = part1_t9*part2_t9,
                unbiased_t9 = (1-(3/(4*(n_trt+n_ctrl-2-q)-1)))*part1_t9*part2_t9,
                v_9 = ((1-corr_pooled^2)/((n_trt*n_ctrl)/(n_trt+n_ctrl)))+((t_9^2)/(2*(n_trt+n_ctrl-2-q))),
               
                combo_LB_t9 = unbiased_t9-1.96*sqrt(v_9),
                combo_UB_t9 = unbiased_t9+1.96*sqrt(v_9),
                combo_in_CI_t9 = ifelse((theta >= combo_LB_t9 & theta <= combo_UB_t9), 1, 0),
            
                # gA3
                part1_t9_trueRho = diffMeansAdjusted/sd_adjusted[1],
                part2_t9_trueRho = sqrt(1-rho^2),
                t_9_trueRho = part1_t9_trueRho*part2_t9_trueRho,
                unbiased_t9_trueRho = (1-(3/(4*(n_trt+n_ctrl-2-q)-1)))*part1_t9_trueRho*part2_t9_trueRho,
                v_9_trueRho = ((1-rho^2)/((n_trt*n_ctrl)/(n_trt+n_ctrl)))+((t_9_trueRho^2)/(2*(n_trt+n_ctrl-2-q))),
               
                combo_LB_t9_trueRho = unbiased_t9_trueRho-1.96*sqrt(v_9_trueRho),
                combo_UB_t9_trueRho = unbiased_t9_trueRho+1.96*sqrt(v_9_trueRho),
                combo_in_CI_t9_trueRho = ifelse((theta >= combo_LB_t9_trueRho & theta <= combo_UB_t9_trueRho), 1, 0),
                
                # gPA
                sPOOLED_anc = sqrt((var_post_group[1]+var_post_group[length(var_post_group)]+
                                        var_pre_group[1]+var_pre_group[length(var_pre_group)])/4),
                
                t_12 = diffMeansAdjusted/sPOOLED_anc,
                unbiased_t12 = (1-(3/(4*(2*(n_trt+n_ctrl-2))-1)))*diffMeansAdjusted/sPOOLED_anc,
                v_12 = (((1-rho^2))/((n_trt*n_ctrl)/(n_trt+n_ctrl)))+((1+rho^2)*(t_12^2)/(4*(n_trt+n_ctrl-2))),
                
                combo_LB_t12 = unbiased_t12-1.96*sqrt(v_12),
                combo_UB_t12 = unbiased_t12+1.96*sqrt(v_12),
                combo_in_CI_t12 = ifelse((theta >= combo_LB_t12 & theta <= combo_UB_t12), 1, 0)
            )
        
        ###-----------------------------------------------------------------###
        
        # Clean dataframe ANCOVA
        df_current_ancova <- df_adjusted %>% 
       #     select(-c(x, y, treatment_indicator, corr_sample)) %>% 
            distinct()

        # Clean dataframe GAINS
        df_current_gains <- df_gains %>% 
        #    select(-c(x_trt, y_trt, x_ctrl, y_ctrl, diff_gains, gain_trt, gain_ctrl)) %>% 
            distinct()

        # Final dataframes
        df_rho_gains <- rbind(df_rho_gains, df_current_gains)
        df_rho_ancova <- rbind(df_rho_ancova, df_current_ancova)
    }    
    
    ###-----------------------------------------------------------------###
    ###-----------------------------------------------------------------###
    ###-----------------------------------------------------------------###
    
    # Properties of estimators
    
    # GAINS
    df_rho_gains <- df_rho_gains %>% 
        mutate(
            # g
            mean_unbiased_t1 = mean(unbiased_t1),
            mean_v1 = mean(v_1),
            var_unbiased_t1 = var(unbiased_t1),
            bias_unbiased_t1 = mean_unbiased_t1 - theta,
            rel_bias_unbiased_t1 = bias_unbiased_t1/theta,
            combo_cov_prob_unbiased_t1 = sum(combo_in_CI_t1)/iterations,

            #gG1
            mean_unbiased_t2 = mean(unbiased_t2),
            mean_v2 = mean(v_2),
            var_unbiased_t2 = var(unbiased_t2),
            bias_unbiased_t2 = mean_unbiased_t2 - theta,
            rel_bias_unbiased_t2 = bias_unbiased_t2/theta,
            combo_cov_prob_unbiased_t2 = sum(combo_in_CI_t2)/iterations,
            
            # gG2
            mean_unbiased_t5 = mean(unbiased_t5),
            mean_v5 = mean(v_5),
            var_unbiased_t5 = var(unbiased_t5),
            bias_unbiased_t5 = mean_unbiased_t5 - theta,
            rel_bias_unbiased_t5 = bias_unbiased_t5/theta,
            combo_cov_prob_unbiased_t5 = sum(combo_in_CI_t5)/iterations,

            # gG3
            mean_unbiased_t5_trueRho = mean(unbiased_t5_trueRho),
            mean_v5_trueRho = mean(v_5_trueRho),
            var_unbiased_t5_trueRho = var(unbiased_t5_trueRho),
            bias_unbiased_t5_trueRho = mean_unbiased_t5_trueRho - theta,
            rel_bias_unbiased_t5_trueRho = bias_unbiased_t5_trueRho/theta,
            combo_cov_prob_unbiased_t5_trueRho = sum(combo_in_CI_t5_trueRho)/iterations,
           
            # gP
            mean_unbiased_t10 = mean(unbiased_t10),
            mean_v10 = mean(v_10),
            var_unbiased_t10 = var(unbiased_t10),
            bias_unbiased_t10 = mean_unbiased_t10 - theta,
            rel_bias_unbiased_t10 = bias_unbiased_t10/theta,
            combo_cov_prob_unbiased_t10 = sum(combo_in_CI_t10)/iterations,

            # gPG
            mean_unbiased_t11 = mean(unbiased_t11),
            mean_v11 = mean(v_11),
            var_unbiased_t11 = var(unbiased_t11),
            bias_unbiased_t11 = mean_unbiased_t11 - theta,
            rel_bias_unbiased_t11 = bias_unbiased_t11/theta,
            combo_cov_prob_unbiased_t11 = sum(combo_in_CI_t11)/iterations
        ) 
    
    df_rho_cumulative_gains <- rbind(df_rho_cumulative_gains, df_rho_gains)
    
    ###-----------------------------------------------------------------###
    
    # ANCOVA
    df_rho_ancova <- df_rho_ancova %>% 
        mutate(
            # gA1
            mean_unbiased_t3 = mean(unbiased_t3),
            mean_v3 = mean(v_3),
            var_unbiased_t3 = var(unbiased_t3),
            bias_unbiased_t3 = mean_unbiased_t3 - theta,
            rel_bias_unbiased_t3 = bias_unbiased_t3/theta,
            combo_cov_prob_unbiased_t3 = sum(combo_in_CI_t3)/iterations,

            # gA2
            mean_unbiased_t9 = mean(unbiased_t9),
            mean_v9 = mean(v_9),
            var_unbiased_t9 = var(unbiased_t9),
            bias_unbiased_t9 = mean_unbiased_t9 - theta,
            rel_bias_unbiased_t9 = bias_unbiased_t9/theta,
            combo_cov_prob_unbiased_t9 = sum(combo_in_CI_t9)/iterations,

            #gA3
            mean_unbiased_t9_trueRho = mean(unbiased_t9_trueRho),
            mean_v9_trueRho = mean(v_9_trueRho),
            var_unbiased_t9_trueRho = var(unbiased_t9_trueRho),
            bias_unbiased_t9_trueRho = mean_unbiased_t9_trueRho - theta,
            rel_bias_unbiased_t9_trueRho = bias_unbiased_t9_trueRho/theta,
            combo_cov_prob_unbiased_t9_trueRho = sum(combo_in_CI_t9_trueRho)/iterations,

            # gPA
            mean_unbiased_t12 = mean(unbiased_t12),
            mean_v12 = mean(v_12),
            var_unbiased_t12 = var(unbiased_t12),
            bias_unbiased_t12= mean_unbiased_t12 - theta,
            rel_bias_unbiased_t12 = bias_unbiased_t12/theta,
            combo_cov_prob_unbiased_t12 = sum(combo_in_CI_t12)/iterations
        ) 
    df_rho_cumulative_ancova <- rbind(df_rho_cumulative_ancova, df_rho_ancova)
}