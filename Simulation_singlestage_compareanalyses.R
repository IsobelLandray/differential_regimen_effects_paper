library(dplyr)
library(MASS)
library(ggplot2)
library(gridExtra)


#Algorithm 2 from Annette Kopp-Schneider paper (2023)
#(R1) select theta_e as parameter for external data
#(R2) select theta in theta_0 for T1E rate calculation
#(R3) select theta1 not in theta_0 for power calculation
#these are chosen at the end of this program as they are coded into a function...
trial1_simteststat_diffsig_singlestg_adj<-function(thetaE_chosen, nsim_noborrow, nsim_mc, delta=NA, empiricalbayes=NA, theta1, theta_c, sigma_c, theta0, bias, withbias,varfactor_e,varfactor_c,increase, c1,c2,alpha_nb){
  set.seed(1234)
  
  sigma_t<-sigma_c*1/sqrt(2)
  sigma_e_v<-sqrt(sigma_c^2*varfactor_e/increase)
  sigma_c_v<-sqrt(sigma_c^2*varfactor_c/increase)
  
  #(R4) repeat for a sufficient number of times nsim:
  #the following function simulates the prior for given weight parameter, mean of external and current data
  simulate_posterior_twoarm<-function(delta, mean_d_t, mean_d_c, mean_d_E){ #function to simulate the posterior mean and sd and then a posterior test statistic
    posterior_mean_control<-((delta*1*mean_d_E)/(sigma_e_v^2)+(1*mean_d_c)/(sigma_c_v^2))/(((delta*1)/(sigma_e_v^2)+1/(sigma_c_v^2)))
    posterior_var_control<-1/((delta*1)/(sigma_e_v^2)+(1/sigma_c_v^2))
    posterior_test<-(mean_d_t-posterior_mean_control-theta0)/sqrt(posterior_var_control+sigma_t^2/1) 
    return(list(posterior_test,mean_d_t-posterior_mean_control-theta0,sqrt(posterior_var_control+sigma_t^2/1)))
  } 
  
  #now just doing for one value of mean_d_E
  prob<-c()
  #(F3.b) determine T1E rate alpha_B(d_E)=max{E_theta[phi_B(D,d_E)]}
  t1e_stg1<-c()
  power_withborrowing_stg1<-c()
  power_withoutborrowing_calibrated_stg1<-c()
  diffinpowers_stg1<-c()
  posterior_test_h0_stg1<-c()
  posterior_mean_h0_stg1<-c()
  posterior_se_h0_stg1<-c()
  posterior_mean_h1_stg1<-c()
  posterior_se_h1_stg1<-c()
  posterior_test_h1_stg1<-c()
  withoutborrow_test_h1_stg1<-c()
  alpha_sa<-c()
  use_sa<-c()
  alpha<-c()
  use<-c()
  
  delta_hat_stg1<-c()
  
  for(k in 1:nsim_mc){
    #(R4.a) Generate one external data set d_E from DE~f_theta_e
    d_E_sim<-rnorm(n=1,mean=thetaE_chosen,sd=sigma_e_v)
    
    mean_dE_sim_stg1<-mean(d_E_sim)
    
    #(R4.b) Generate one current data set d_theta from D~f_theta, theta in theta_0 selected in (R2)
    #control data from current trial
    d_c_sim<-rnorm(n=1,mean=1*theta_c,sd=sigma_c_v)
    
    mean_dc_sim_stg1<-mean(d_c_sim)
    
    #treatment arm data from current trial under h0
    d_t_h0<-rnorm(n=1,mean=1*(theta0+bias*withbias),sd=sigma_t)
    
    mean_d_t_h0_stg1<-mean(d_t_h0)
    
    
    #(R4.c) Record test stat (algorithm says test decision but this is enough)
    #get the posterior test statistic using each external dataset with each of nsim current datasets  (first with current data simulated under h0)
    
    #IF EMPIRICAL BAYES
    if(empiricalbayes==1){
      delta_hat_stg1[k]<-(sigma_e_v^2/1)/(max(c((mean_dc_sim_stg1-mean_dE_sim_stg1)^2,(sigma_c_v^2/1+sigma_e_v^2/1)))-sigma_c_v^2/1)
    }else{
      delta_hat_stg1[k]<-delta #IF FIXED PRIOR, USE INPUTTED VALUE
    }
    
    posterior_test_h0_stg1[k]<-simulate_posterior_twoarm(delta=delta_hat_stg1[k],mean_d_E=mean_dE_sim_stg1, mean_d_t=mean_d_t_h0_stg1, mean_d_c=mean_dc_sim_stg1)[[1]]
    posterior_mean_h0_stg1[k]<-simulate_posterior_twoarm(delta=delta_hat_stg1[k],mean_d_E=mean_dE_sim_stg1, mean_d_t=mean_d_t_h0_stg1, mean_d_c=mean_dc_sim_stg1)[[2]]
    posterior_se_h0_stg1[k]<-simulate_posterior_twoarm(delta=delta_hat_stg1[k],mean_d_E=mean_dE_sim_stg1, mean_d_t=mean_d_t_h0_stg1, mean_d_c=mean_dc_sim_stg1)[[3]]    
    
    #(R4.d) Generate one current data set d_theta1 from D~f_theta1
    #treatment arm data from current trial under h1 & use same current control and external control data
    d_t_h1<-rnorm(n=1,mean=1*(theta1+bias*withbias),sd=sigma_t)
    
    mean_d_t_h1_stg1<-mean(d_t_h1)
    
    
    #(R4.e) Record test stat 
    #get the posterior test statistic using each external dataset with each of nsim current datasets  (second with current data simulated under h1)
    posterior_test_h1_stg1[k]<-simulate_posterior_twoarm(delta=delta_hat_stg1[k],mean_d_E=mean_dE_sim_stg1, mean_d_t=mean_d_t_h1_stg1, mean_d_c=mean_dc_sim_stg1)[[1]]
    posterior_mean_h1_stg1[k]<-simulate_posterior_twoarm(delta=delta_hat_stg1[k],mean_d_E=mean_dE_sim_stg1, mean_d_t=mean_d_t_h1_stg1, mean_d_c=mean_dc_sim_stg1)[[2]]
    posterior_se_h1_stg1[k]<-simulate_posterior_twoarm(delta=delta_hat_stg1[k],mean_d_E=mean_dE_sim_stg1, mean_d_t=mean_d_t_h1_stg1, mean_d_c=mean_dc_sim_stg1)[[3]]    
    
    
    #without borrowing:
    withoutborrow_test_h1_stg1[k]<-(mean_d_t_h1_stg1-mean_dc_sim_stg1-theta0)/sqrt(sigma_t^2/1+sigma_c_v^2/1)
    
    alpha[k]<-if_else(mean_dc_sim_stg1>mean_dE_sim_stg1,c1,c2)
    #for either comparison: if internal>external -> use c1,if not -> use c2
    if(alpha[k]==c1){use[k]<-"c1"}else{use[k]<-"c2"}
  }
  
  
  #(R5) Average test decisions from (R4.c) to obtain T1E
  t1e_stg1<-mean(posterior_test_h0_stg1>qnorm(1-alpha))
  #(R6) Average test decisions from (R4.e) to obtain power with borrowing
  power_withborrowing_stg1<-mean(posterior_test_h1_stg1>qnorm(1-alpha))
  power_withoutborrowing_calibrated_stg1<-mean(withoutborrow_test_h1_stg1>qnorm(1-t1e_stg1))
  diffinpowers_stg1<-power_withborrowing_stg1-power_withoutborrowing_calibrated_stg1
  
  xaxis<-(theta_c-thetaE_chosen)/sigma_c_v
  
  teststat_fixed_stg1<-c()
  teststat_fixed_stg1_h0<-c()
  effect_fixed_stg1_h1<-c()
  effect_fixed_stg1_h0<-c()
  se_fixed_stg1<-c()
  
  for(i in 1:nsim_noborrow){
    d_t_h1<-rnorm(n=1,mean=1*(theta1+bias*withbias),sd=sigma_t)
    
    mean_d_t_h1_stg1<-mean(d_t_h1)
    
    d_t_h0<-rnorm(n=1,mean=1*(theta0+bias*withbias),sd=sigma_t)
    
    mean_d_t_h0_stg1<-mean(d_t_h0)
    
    d_c_sim<-rnorm(n=1,mean=1*theta_c,sd=sigma_c_v)
    
    mean_d_c_sim_stg1<-mean(d_c_sim)
    
    teststat_fixed_stg1[i]<-(mean_d_t_h1_stg1-mean_d_c_sim_stg1-theta0)/sqrt(sigma_t^2/1+sigma_c_v^2/1)
    
    teststat_fixed_stg1_h0[i]<-(mean_d_t_h0_stg1-mean_d_c_sim_stg1-theta0)/sqrt(sigma_t^2/1+sigma_c_v^2/1)
    
    #added to be able to plot effect estimates and CIs
    effect_fixed_stg1_h1[i]<-(mean_d_t_h1_stg1-mean_d_c_sim_stg1-theta0)
    effect_fixed_stg1_h0[i]<-(mean_d_t_h0_stg1-mean_d_c_sim_stg1-theta0)
    se_fixed_stg1[i]<-sqrt(sigma_t^2/1+sigma_c_v^2/1)
    
  }
  
  power_without_borrowing_nocalib_stg1<-mean(teststat_fixed_stg1>qnorm(1-alpha_nb))
  
  t1e_without_borrowing_nocalib_stg1<-mean(teststat_fixed_stg1_h0>qnorm(1-alpha_nb))
  
  ##power and t1e if simple average of controls
  teststat_fixed_stg1_sah1<-c()
  teststat_fixed_stg1_sah0<-c()
  effect_sa_stg1_h1<-c()
  effect_sa_stg1_h0<-c()
  se_sa_stg1<-c()
  for(i in 1:nsim_noborrow){
    d_t_h1<-rnorm(n=1,mean=1*(theta1+bias*withbias),sd=sigma_t)
    d_t_h0<-rnorm(n=1,mean=1*(theta0+bias*withbias),sd=sigma_t)
    d_c_sim<-rnorm(n=1,mean=1*theta_c,sd=sigma_c_v)
    d_E_sim<-rnorm(n=1,mean=1*thetaE_chosen,sd=sigma_e_v)
    
    mean_d_t_h1_stg1<-mean(d_t_h1)
    mean_d_t_h0_stg1<-mean(d_t_h0)
    mean_d_c_sim_stg1<-mean(d_c_sim)
    mean_d_e_sim_stg1<-mean(d_E_sim)
    
    pooledvarcontrol_stg1<-1/(1/sigma_c_v^2+1/sigma_e_v^2)
    
    teststat_fixed_stg1_sah1[i]<-(mean_d_t_h1_stg1-(0.5*mean_d_c_sim_stg1+0.5*mean_d_e_sim_stg1)-theta0)/sqrt(sigma_t^2/1+pooledvarcontrol_stg1)
    teststat_fixed_stg1_sah0[i]<-(mean_d_t_h0_stg1-(0.5*mean_d_c_sim_stg1+0.5*mean_d_e_sim_stg1)-theta0)/sqrt(sigma_t^2/1+pooledvarcontrol_stg1)
    
    
    #added to be able to plot effect estimates and CIs
    effect_sa_stg1_h1[i]<-(mean_d_t_h1_stg1-(0.5*mean_d_c_sim_stg1+0.5*mean_d_e_sim_stg1)-theta0)
    effect_sa_stg1_h0[i]<-(mean_d_t_h0_stg1-(0.5*mean_d_c_sim_stg1+0.5*mean_d_e_sim_stg1)-theta0)
    se_sa_stg1[i]<-sqrt(sigma_t^2/1+pooledvarcontrol_stg1)
    
    alpha_sa[i]<-if_else(mean_d_c_sim_stg1>mean_d_e_sim_stg1,c1,c2)
    #for either comparison: if internal>external -> use c1,if not -> use c2
    if(alpha_sa[i]==c1){use_sa[i]<-"c1"}else{use_sa[i]<-"c2"}
  }
  
  
  power_sa_stg1<-mean(teststat_fixed_stg1_sah1>qnorm(1-alpha_sa))
  
  t1e_sa_stg1<-mean(teststat_fixed_stg1_sah0>qnorm(1-alpha_sa))
  
  #outputs all operating characteristics
  return(list(power_without_borrowing_nocalib_stg1=power_without_borrowing_nocalib_stg1,
              t1e_without_borrowing_nocalib_stg1=t1e_without_borrowing_nocalib_stg1,
              power_sa_stg1=power_sa_stg1,
              t1e_sa_stg1=t1e_sa_stg1,
              xaxis=xaxis, 
              t1e_stg1=t1e_stg1,
              power_withborrowing_stg1=power_withborrowing_stg1,
              power_withoutborrowing_calibrated_stg1=power_withoutborrowing_calibrated_stg1,
              posterior_mean_h0_stg1=posterior_mean_h0_stg1, 
              posterior_se_h0_stg1=posterior_se_h0_stg1, 
              posterior_mean_h1_stg1=posterior_mean_h1_stg1, 
              posterior_se_h1_stg1=posterior_se_h1_stg1,
              posterior_test_h0_stg1=posterior_test_h0_stg1, 
              posterior_test_h1_stg1=posterior_test_h1_stg1,
              effect_fixed_stg1_h1=effect_fixed_stg1_h1,
              effect_fixed_stg1_h0=effect_fixed_stg1_h0,
              se_fixed_stg1=se_fixed_stg1,
              effect_sa_stg1_h1=effect_sa_stg1_h1,
              effect_sa_stg1_h0=effect_sa_stg1_h0,
              se_sa_stg1=se_sa_stg1,
              delta_hat_stg1=delta_hat_stg1,
              alpha=alpha,
              use=use,
              alpha_sa=alpha_sa,
              use_sa=use_sa
  ))
}

#This function runs the function above for each comparison with differing levels of Differential regimen effect
datasetcreate_diffsig_adj<-function(sigma_c,c1,c2, bias_comp1,varfactor_pbo1,increase,nsim){
  
  df_bb<-data.frame()
  deltas <- vector("list", length = length(bias_comp1))
  uses <- vector("list", length = length(bias_comp1))
  
  for (i in 1:length(bias_comp1)) {
    
    eb_compar1 <- trial1_simteststat_diffsig_singlestg_adj(
      thetaE_chosen = 0,
      nsim_noborrow = nsim,
      nsim_mc = nsim,
      delta = NA,
      empiricalbayes=1,
      theta1 = 0.7506,
      theta_c = bias_comp1[i],
      sigma_c = sigma_c,
      c1 = c1,
      c2 = c2,
      theta0 = 0,
      bias = bias_comp1[i],
      withbias = 1,
      varfactor_c=varfactor_pbo1,
      varfactor_e=1,
      increase=increase,
      alpha_nb=0.025
      )
    
    eb_compar2 <- trial1_simteststat_diffsig_singlestg_adj(
      thetaE_chosen = bias_comp1[i],
      nsim_noborrow = nsim,
      nsim_mc = nsim,
      delta = NA,
      empiricalbayes=1,
      theta1 = 0.7506,
      theta_c = 0,
      sigma_c = sigma_c,
      c1 = c1,
      c2 = c2,
      theta0 = 0,
      bias = 0,
      withbias = 0,
      varfactor_c=1,
      varfactor_e=varfactor_pbo1,
      increase=increase,
      alpha_nb=0.025
      )
    
    #_eb means from empirical bayes prior
    #no _ means from fixed power prior with delta=0.5
    df_bb <- rbind(
      df_bb,
      data.frame(
        bias_compar1 = mean(bias_comp1[i]),
        xaxis=eb_compar1$xaxis, 
        power_without_borrowing_nocalib_stg1_c1_eb=eb_compar1$power_without_borrowing_nocalib_stg1,
        t1e_without_borrowing_nocalib_stg1_c1_eb=eb_compar1$t1e_without_borrowing_nocalib_stg1,
        power_sa_stg1_c1_eb=eb_compar1$power_sa_stg1,
        t1e_sa_stg1_c1_eb=eb_compar1$t1e_sa_stg1,
        t1e_stg1_c1_eb=eb_compar1$t1e_stg1,
        power_withborrowing_stg1_c1_eb=eb_compar1$power_withborrowing_stg1,
        power_withoutborrowing_calibrated_stg1_c1_eb=eb_compar1$power_withoutborrowing_calibrated_stg1,
        posterior_mean_h0_stg1_c1_eb=mean(eb_compar1$posterior_mean_h0_stg1), 
        posterior_se_h0_stg1_c1_eb=mean(eb_compar1$posterior_se_h0_stg1), 
        posterior_mean_h1_stg1_c1_eb=mean(eb_compar1$posterior_mean_h1_stg1), 
        posterior_se_h1_stg1_c1_eb=mean(eb_compar1$posterior_se_h1_stg1),
        posterior_test_h0_stg1_c1_eb=mean(eb_compar1$posterior_test_h0_stg1), 
        posterior_test_h1_stg1_c1_eb=mean(eb_compar1$posterior_test_h1_stg1),
        power_without_borrowing_nocalib_stg1_c2_eb=eb_compar2$power_without_borrowing_nocalib_stg1,
        t1e_without_borrowing_nocalib_stg1_c2_eb=eb_compar2$t1e_without_borrowing_nocalib_stg1,
        power_sa_stg1_c2_eb=eb_compar2$power_sa_stg1,
        t1e_sa_stg1_c2_eb=eb_compar2$t1e_sa_stg1,
        t1e_stg1_c2_eb=eb_compar2$t1e_stg1,
        power_withborrowing_stg1_c2_eb=eb_compar2$power_withborrowing_stg1,
        power_withoutborrowing_calibrated_stg1_c2_eb=eb_compar2$power_withoutborrowing_calibrated_stg1,
        posterior_mean_h0_stg1_c2_eb=mean(eb_compar2$posterior_mean_h0_stg1), 
        posterior_se_h0_stg1_c2_eb=mean(eb_compar2$posterior_se_h0_stg1), 
        posterior_mean_h1_stg1_c2_eb=mean(eb_compar2$posterior_mean_h1_stg1), 
        posterior_se_h1_stg1_c2_eb=mean(eb_compar2$posterior_se_h1_stg1),
        posterior_test_h0_stg1_c2_eb=mean(eb_compar2$posterior_test_h0_stg1), 
        posterior_test_h1_stg1_c2_eb=mean(eb_compar2$posterior_test_h1_stg1),
        varfactor_pbo1=varfactor_pbo1,
        effect_fixed_stg1_h1_c1=mean(eb_compar1$effect_fixed_stg1_h1),
        effect_fixed_stg1_h0_c1=mean(eb_compar1$effect_fixed_stg1_h0),
        se_fixed_stg1_c1=mean(eb_compar1$se_fixed_stg1),
        effect_sa_stg1_h1_c1=mean(eb_compar1$effect_sa_stg1_h1),
        effect_sa_stg1_h0_c1=mean(eb_compar1$effect_sa_stg1_h0),
        se_sa_stg1_c1=mean(eb_compar1$se_sa_stg1),
        effect_fixed_stg1_h1_c2=mean(eb_compar2$effect_fixed_stg1_h1),
        effect_fixed_stg1_h0_c2=mean(eb_compar2$effect_fixed_stg1_h0),
        se_fixed_stg1_c2=mean(eb_compar2$se_fixed_stg1),
        effect_sa_stg1_h1_c2=mean(eb_compar2$effect_sa_stg1_h1),
        effect_sa_stg1_h0_c2=mean(eb_compar2$effect_sa_stg1_h0),
        se_sa_stg1_c2=mean(eb_compar2$se_sa_stg1)
      ))
    
    deltas[[i]] <- list(bias_comp1[i], eb_compar1$delta_hat_stg1,
                        eb_compar2$delta_hat_stg1)
    
    uses[[i]] <- list(bias_comp1[i], mean(eb_compar1$use=="c1"),
                        mean(eb_compar2$use=="c1"))
    
    if (i %% 10 == 0 || i == length(bias_comp1)) { 
      cat("Completed", i, "out of", length(bias_comp1), "iterations\n")
    }
    
    
  }
  #EMPIRICAL BAYES
  df_bb$LL_h0c1_stg1_eb<-df_bb$posterior_mean_h0_stg1_c1_eb-1.96*df_bb$posterior_se_h0_stg1_c1_eb
  df_bb$UL_h0c1_stg1_eb<-df_bb$posterior_mean_h0_stg1_c1_eb+1.96*df_bb$posterior_se_h0_stg1_c1_eb
  df_bb$LL_h1c1_stg1_eb<-df_bb$posterior_mean_h1_stg1_c1_eb-1.96*df_bb$posterior_se_h1_stg1_c1_eb
  df_bb$UL_h1c1_stg1_eb<-df_bb$posterior_mean_h1_stg1_c1_eb+1.96*df_bb$posterior_se_h1_stg1_c1_eb
  df_bb$LL_h0c2_stg1_eb<-df_bb$posterior_mean_h0_stg1_c2_eb-1.96*df_bb$posterior_se_h0_stg1_c2_eb
  df_bb$UL_h0c2_stg1_eb<-df_bb$posterior_mean_h0_stg1_c2_eb+1.96*df_bb$posterior_se_h0_stg1_c2_eb
  df_bb$LL_h1c2_stg1_eb<-df_bb$posterior_mean_h1_stg1_c2_eb-1.96*df_bb$posterior_se_h1_stg1_c2_eb
  df_bb$UL_h1c2_stg1_eb<-df_bb$posterior_mean_h1_stg1_c2_eb+1.96*df_bb$posterior_se_h1_stg1_c2_eb
  
  
  #NO BORROWING - IE COMPARE TO N/2 PATIENTS
  df_bb$LL_h0c1_stg1_noborrow<-df_bb$effect_fixed_stg1_h0_c1-1.96*df_bb$se_fixed_stg1_c1
  df_bb$UL_h0c1_stg1_noborrow<-df_bb$effect_fixed_stg1_h0_c1+1.96*df_bb$se_fixed_stg1_c1
  df_bb$LL_h1c1_stg1_noborrow<-df_bb$effect_fixed_stg1_h1_c1-1.96*df_bb$se_fixed_stg1_c1
  df_bb$UL_h1c1_stg1_noborrow<-df_bb$effect_fixed_stg1_h1_c1+1.96*df_bb$se_fixed_stg1_c1
  df_bb$LL_h0c2_stg1_noborrow<-df_bb$effect_fixed_stg1_h0_c2-1.96*df_bb$se_fixed_stg1_c1
  df_bb$UL_h0c2_stg1_noborrow<-df_bb$effect_fixed_stg1_h0_c2+1.96*df_bb$se_fixed_stg1_c1
  df_bb$LL_h1c2_stg1_noborrow<-df_bb$effect_fixed_stg1_h1_c2-1.96*df_bb$se_fixed_stg1_c2
  df_bb$UL_h1c2_stg1_noborrow<-df_bb$effect_fixed_stg1_h1_c2+1.96*df_bb$se_fixed_stg1_c2
  
  #SCALED AVERAGE
  df_bb$LL_h0c1_stg1_sa<-df_bb$effect_sa_stg1_h0_c1-1.96*df_bb$se_sa_stg1_c1
  df_bb$UL_h0c1_stg1_sa<-df_bb$effect_sa_stg1_h0_c1+1.96*df_bb$se_sa_stg1_c1
  df_bb$LL_h1c1_stg1_sa<-df_bb$effect_sa_stg1_h1_c1-1.96*df_bb$se_sa_stg1_c1
  df_bb$UL_h1c1_stg1_sa<-df_bb$effect_sa_stg1_h1_c1+1.96*df_bb$se_sa_stg1_c1
  df_bb$LL_h0c2_stg1_sa<-df_bb$effect_sa_stg1_h0_c2-1.96*df_bb$se_sa_stg1_c1
  df_bb$UL_h0c2_stg1_sa<-df_bb$effect_sa_stg1_h0_c2+1.96*df_bb$se_sa_stg1_c1
  df_bb$LL_h1c2_stg1_sa<-df_bb$effect_sa_stg1_h1_c2-1.96*df_bb$se_sa_stg1_c2
  df_bb$UL_h1c2_stg1_sa<-df_bb$effect_sa_stg1_h1_c2+1.96*df_bb$se_sa_stg1_c2
  
  
  #Fully blinded has equivalent operating characteristics to the simple average method with Differential regimen effect=0
  
  return(list(df_bb,deltas,uses))
}


out_dta<-datasetcreate_diffsig_adj(sigma_c=0.2,c1=0.025,c2=0.025, bias_comp1=seq(0,0.25,0.01),varfactor_pbo1=1,increase=1,nsim=10000)
out_dta<-out_dta[[1]]

doubledummy_power_c1<-out_dta$power_sa_stg1_c1[out_dta$bias_compar1==0] #THIS VALUE IS THE "TIM LINE" THAT CAN BE ADDED TO THE PLOT AS REFERENCE FOR THE DUMMY control TRIAL - N/2 IN EACH control ARM BUT DUMMIED SO THE controls ARE THE SAME SO N IN OVERALL control GROUP REALLY !!
doubledummy_power_c2<-out_dta$power_sa_stg1_c2[out_dta$bias_compar1==0] #THIS VALUE IS THE "TIM LINE" THAT CAN BE ADDED TO THE PLOT AS REFERENCE FOR THE DUMMY control TRIAL - N/2 IN EACH control ARM BUT DUMMIED SO THE controls ARE THE SAME SO N IN OVERALL control GROUP REALLY !!
doubledummy_t1estg1_c1<-out_dta$t1e_sa_stg1_c1[out_dta$bias_compar1==0] #THIS VALUE IS THE "TIM LINE" THAT CAN BE ADDED TO THE PLOT AS REFERENCE FOR THE DUMMY control TRIAL - N/2 IN EACH control ARM BUT DUMMIED SO THE controls ARE THE SAME SO N IN OVERALL control GROUP REALLY !!
doubledummy_t1estg1_c2<-out_dta$t1e_sa_stg1_c2[out_dta$bias_compar1==0] #THIS VALUE IS THE "TIM LINE" THAT CAN BE ADDED TO THE PLOT AS REFERENCE FOR THE DUMMY control TRIAL - N/2 IN EACH control ARM BUT DUMMIED SO THE controls ARE THE SAME SO N IN OVERALL control GROUP REALLY !!


out_cvadj_sac1<-datasetcreate_diffsig_adj(sigma_c=0.20, c1=0.0032, c2=0.025,
                                          bias_comp1=seq(0,0.25,0.01),varfactor_pbo1=1,increase=1,nsim=10000)
out_cvadj_sac1_dta<-out_cvadj_sac1[[1]]
out_cvadj_sac1_dta$t1e_sa_stg1_c1_eb
out_cvadj_sac1_dta$power_sa_stg1_c1_eb
out_cvadj_sac1_dta$t1e_sa_stg1_c2_eb
out_cvadj_sac1_dta$power_sa_stg1_c2_eb

out_cvadj_sac1_dta_OBS<-out_cvadj_sac1_dta


out_cvadj_bbc1<-datasetcreate_diffsig_adj(sigma_c=0.20, c1=0.0065, c2=0.025,
                                          bias_comp1=seq(0,0.25,0.01),varfactor_pbo1=1,increase=1,nsim=10000)
out_cvadj_bbc1_dta<-out_cvadj_bbc1[[1]]
out_cvadj_bbc1_dta$t1e_stg1_c1_eb
out_cvadj_bbc1_dta$power_withborrowing_stg1_c1_eb
out_cvadj_bbc1_dta$t1e_stg1_c2_eb
out_cvadj_bbc1_dta$power_withborrowing_stg1_c2_eb

out_cvadj_bbc1_dta_OBS<-out_cvadj_bbc1_dta

power_c1<-ggplot(out_dta, aes(x=bias_compar1,y=power_withborrowing_stg1_c1_eb,color="Bayesian borrowing, Empirical Bayes"))+geom_smooth(se=F,linetype="dashed")+
  geom_smooth(aes(y=out_dta$power_without_borrowing_nocalib_stg1_c1, color="Comparison to only the N/2 matching controls"),se=F)+
  geom_smooth(aes(y=out_dta$power_sa_stg1_c1, color="Simple pooling of controls"),se=F,linetype="dashed")+
  geom_smooth(aes(y=out_cvadj_sac1_dta_OBS$power_sa_stg1_c1, color="Simple pooling of controls"),se=F)+
  geom_line(aes(y=doubledummy_power_c1,color="Fully blinded design with common control arm (size N)"),size=1.05)+
  xlab("Differential regimen effect (Regimen 1 - Regimen 2) (units per year)")+ylab("Power in comparison 1")+
  scale_color_manual(
    values = c(
      "Bayesian borrowing, Empirical Bayes"="orange",
      "Comparison to only the N/2 matching controls"="dodgerblue2",
      "Simple pooling of controls"="red",
      "Fully blinded design with common control arm (size N)"="magenta3"
    )
  )+
  theme_bw()+
  labs(title="Comparison 1 (active 1 vs control): Power",color = "Method")+
  guides(linetype=guide_legend(override.aes = list(color="black")))+ylim(c(0.8,1))+
  theme(legend.title=element_text(face="bold"),legend.direction="horizontal",legend.box.background=element_rect(color="black",size=1))+
  theme(plot.title=element_text(face="bold",size=10,hjust=0.5))+xlim(c(0,0.25))+
  geom_smooth(aes(y=out_cvadj_bbc1_dta_OBS$power_withborrowing_stg1_c1_eb, color="Bayesian borrowing, Empirical Bayes"),se=F)


power_c2<-ggplot(out_dta, aes(x=bias_compar1,y=power_withborrowing_stg1_c2_eb,color="Bayesian borrowing, Empirical Bayes"))+geom_smooth(se=F,linetype="dashed")+
  geom_smooth(aes(y=out_dta$power_without_borrowing_nocalib_stg1_c2, color="Comparison to only the N/2 matching controls"),se=F)+
  geom_smooth(aes(y=out_dta$power_sa_stg1_c2, color="Simple pooling of controls"),se=F,linetype="dashed")+
  geom_smooth(aes(y=out_cvadj_sac1_dta_OBS$power_sa_stg1_c2, color="Simple pooling of controls"),se=F)+
  geom_line(aes(y=doubledummy_power_c2,color="Fully blinded design with common control arm (size N)"),size=1.05)+
  xlab("Differential regimen effect (Regimen 1 - Regimen 2) (units per year)")+ylab("Power in comparison 2")+
  scale_color_manual(
    values = c(
      "Bayesian borrowing, Empirical Bayes"="orange",
      "Comparison to only the N/2 matching controls"="dodgerblue2",
      "Simple pooling of controls"="red",
      "Fully blinded design with common control arm (size N)"="magenta3"
    )
  )+
  theme_bw()+
  labs(title="Comparison 2 (active 2 vs control): Power",color = "Method")+
  guides(linetype=guide_legend(override.aes = list(color="black")))+ylim(c(0.8,1))+
  theme(plot.title=element_text(face="bold",size=10,hjust=0.5))+xlim(c(0,0.25))+
  geom_smooth(aes(y=out_cvadj_bbc1_dta_OBS$power_withborrowing_stg1_c2_eb, color="Bayesian borrowing, Empirical Bayes"),se=F)

t1e_stg1_c1<-ggplot(out_dta, aes(x=bias_compar1,y=t1e_stg1_c1_eb,color="Bayesian borrowing, Empirical Bayes"))+geom_smooth(se=F,linetype="dashed")+
  geom_smooth(aes(y=out_dta$t1e_without_borrowing_nocalib_stg1_c1, color="Comparison to only the N/2 matching controls"),se=F)+
  geom_smooth(aes(y=out_dta$t1e_sa_stg1_c1, color="Simple pooling of controls"),se=F,linetype="dashed")+
  geom_smooth(aes(y=out_cvadj_sac1_dta_OBS$t1e_sa_stg1_c1, color="Simple pooling of controls"),se=F)+
  geom_line(aes(y=doubledummy_t1estg1_c1,color="Fully blinded design with common control arm (size N)"),size=1.05)+
  xlab("Differential regimen effect (Regimen 1 - Regimen 2) (units per year)")+ylab("T1E in comparison 1")+
  scale_color_manual(
    values = c(
      "Bayesian borrowing, Empirical Bayes"="orange",
      "Comparison to only the N/2 matching controls"="dodgerblue2",
      "Simple pooling of controls"="red",
      "Fully blinded design with common control arm (size N)"="magenta3"
    )
  )+
  theme_bw()+
  labs(title="Comparison 1 (active 1 vs control): T1E",color = "Method")+
  guides(linetype=guide_legend(override.aes = list(color="black")))+geom_hline(yintercept=0.025,linetype="dashed")+ylim(c(0,0.1))+
  theme(plot.title=element_text(face="bold",size=10,hjust=0.5))+xlim(c(0,0.25))+
  geom_smooth(aes(y=out_cvadj_bbc1_dta_OBS$t1e_stg1_c1_eb, color="Bayesian borrowing, Empirical Bayes"),se=F)

t1e_stg1_c2<-ggplot(out_dta, aes(x=bias_compar1,y=t1e_stg1_c2_eb,color="Bayesian borrowing, Empirical Bayes"))+geom_smooth(se=F,linetype="dashed")+
  geom_smooth(aes(y=out_dta$t1e_without_borrowing_nocalib_stg1_c2, color="Comparison to only the N/2 matching controls"),se=F)+
  geom_smooth(aes(y=out_dta$t1e_sa_stg1_c2, color="Simple pooling of controls"),se=F,linetype="dashed")+
  geom_smooth(aes(y=out_cvadj_sac1_dta_OBS$t1e_sa_stg1_c2, color="Simple pooling of controls"),se=F)+
  geom_line(aes(y=doubledummy_t1estg1_c2,color="Fully blinded design with common control arm (size N)"),size=1.05)+
  xlab("Differential regimen effect (Regimen 1 - Regimen 2) (units per year)")+ylab("T1E in comparison 2")+
  scale_color_manual(
    values = c(
      "Bayesian borrowing, Empirical Bayes"="orange",
      "Comparison to only the N/2 matching controls"="dodgerblue2",
      "Simple pooling of controls"="red",
      "Fully blinded design with common control arm (size N)"="magenta3"
    )
  )+
  theme_bw()+
  labs(title="Comparison 2 (active 2 vs control): T1E",color = "Method")+
  guides(linetype=guide_legend(override.aes = list(color="black")))+geom_hline(yintercept=0.025,linetype="dashed")+
  ylim(c(0,0.1))+theme(plot.title=element_text(face="bold",size=10,hjust=0.5))+xlim(c(0,0.25))+
  geom_smooth(aes(y=out_cvadj_bbc1_dta_OBS$t1e_stg1_c2_eb, color="Bayesian borrowing, Empirical Bayes"),se=F)




get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}
legend<-get_legend(power_c1)

library("grid")
library("gridExtra")

maingrid_overallpower_stg1t1e<-grid.arrange(arrangeGrob(power_c1+theme(legend.position="none"),
                                                        power_c2+theme(legend.position="none"),
                                                        t1e_stg1_c1+theme(legend.position="none"),
                                                        t1e_stg1_c2+theme(legend.position="none")),legend,ncol=1,
                                            heights=c(2.7,0.3))


##TO GET 90% POWER:

#simple average 
out_aim90_sac1<-datasetcreate_diffsig_adj(sigma_c=0.20, c1=0.0028, c2=0.025,
                                          bias_comp1=seq(0,0.25,0.01),varfactor_pbo1=1,increase=1.6,nsim=1000)
out_aim90_sac1_dta<-out_aim90_sac1[[1]]
out_aim90_sac1_dta$t1e_sa_stg1_c1_eb
out_aim90_sac1_dta$power_sa_stg1_c1_eb
out_aim90_sac1_dta$t1e_sa_stg1_c2_eb
out_aim90_sac1_dta$power_sa_stg1_c2_eb

#bayesian borrowing
out_aim90_bbc1<-datasetcreate_diffsig_adj(sigma_c=0.20, c1=0.0044, c2=0.025,
                                          bias_comp1=seq(0,0.25,0.01),varfactor_pbo1=1,increase=1.45,nsim=1000)
out_aim90_bbc1_dta<-out_aim90_bbc1[[1]]
out_aim90_bbc1_dta$t1e_stg1_c1_eb
out_aim90_bbc1_dta$power_withborrowing_stg1_c1_eb
out_aim90_bbc1_dta$t1e_stg1_c2_eb
out_aim90_bbc1_dta$power_withborrowing_stg1_c2_eb

##TO GET 95% POWER:

#simple average 
out_aim95_sac1<-datasetcreate_diffsig_adj(sigma_c=0.20, c1=0.0025, c2=0.025,
                                          bias_comp1=seq(0,0.25,0.01),varfactor_pbo1=1,increase=2.3,nsim=1000)
out_aim95_sac1_dta<-out_aim95_sac1[[1]]
out_aim95_sac1_dta$t1e_sa_stg1_c1_eb
out_aim95_sac1_dta$power_sa_stg1_c1_eb
out_aim95_sac1_dta$t1e_sa_stg1_c2_eb
out_aim95_sac1_dta$power_sa_stg1_c2_eb

#bayesian borrowing
out_aim95_bbc1<-datasetcreate_diffsig_adj(sigma_c=0.20, c1=0.0055, c2=0.025,
                                          bias_comp1=seq(0,0.25,0.01),varfactor_pbo1=1,increase=2.3,nsim=1000)
out_aim95_bbc1_dta<-out_aim95_bbc1[[1]]
out_aim95_bbc1_dta$t1e_stg1_c1_eb
out_aim95_bbc1_dta$power_withborrowing_stg1_c1_eb
out_aim95_bbc1_dta$t1e_stg1_c2_eb
out_aim95_bbc1_dta$power_withborrowing_stg1_c2_eb
  
  
  
  