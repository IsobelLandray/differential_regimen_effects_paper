# differential_regimen_effects

Files contained in this repository can be used to reproduce the results of the simulation studies reported in the paper entitled
# I Landray, J Nicholas, M Burnell, J Carpenter (2026). *Patient-reported outcomes and different treatment regimens in Phase III multi-arm multi-stage trials: Can Bayesian Borrowing gain power?*

The file is **Simulation_singlestage_compareanalyses.R**. Details of its functions and their input/output are found in this README file.

## **Simulation_singlestage_compareanalyses.R**
The file "Simulation_singlestage_compareanalyses.R" contains the simulation functions: trial1_simteststat_diffsig_fourstgs_adj and datasetcreate_diffsig.The final line of the file has an example call of the functions.

#### Function: trial1_simteststat_diffsig_singlestg_adj

The function "trial1_simteststat_diffsig_singlestg_adj" takes the following inputs to simulate operating characteristics of a four-stage trial under the comparison to only N/2 matching placebos, simple average and Bayesian borrowing with Empirical Bayes methods:

- *thetaE_chosen* is the "external" ("other" regimen) placebo data true mean,

- *nsim_noborrow* is the number of simulations for the analysis without borrowing, 

- *nsim_mc* is the number of simulations for the other analyses, 

- *delta* can be specified for the fixed power prior, 

- *empiricalbayes* takes 1 (yes) or 0 (no) to specify if should use Empirical Bayes power prior rather than fixed power prior, 

- *theta1* is the mean under the alternative hypothesis, 

- *theta_c* is the "current" ("matching" regimen placebo arm) data true mean, 

- *sigma_c* is a vector of the standard deviation of the outcome in the "current" ("matching regimen) placebo data for each stage, 

- *alpha* is a vector of length 4 for the alpha levels at each stage, 

- *theta0* is the mean under the null hypothesis, bias is the size of the differential placebo effect, 

- *withbias* is an indicator (0 or 1) to specify if there is a differential placebo effect or not, 

- *rho_ij* is the correlation between the estimates at stage *i* and *j*, 

- *varfactor_e*, *varfactor_c* and *increase* control the relative variances between the placebo subgroups (but are all kept at 1 for the results presented in the manuscript),

- *c1* is the alpha-level to use for the comparison where the control is observed to have higher mean than the "other" control

- *c2* is the alpha-level to use for the comparison where the control is observed to have lower mean than the "other" control

- *alpha_nb* is the alpha-level to use for the analysis without borrowing

#### Function: datasetcreate_diffsig

The function "datasetcreate_diffsig" runs the above function across a range of differential placebo effect sizes, with a fewer number of values to specify:

- *sigma_c* is a vector of the standard deviation of the outcome in the "current" ("matching regimen) placebo data for each stage, 

- *c1* and *c2* are as described in the above (for non-adjusted analyses just enter the same alpha-level for both c1 and c2), 

- *bias_comp1* is the a vector of the differential placebo effect sizes to simulate,

- *varfactor_pbo1* controls how much larger the placebo subgroup 1's variance is compared to the other placebo subgroup (e.g. 1.1 means 10\% increase),

- *increase* is increase in placebo arm sample size (used when increase placebo arm to improve power).

**The dataset outputted will have MANY variables. The ones used in the figures in the paper are:**

**Type I error:**

- *t1e_without_borrowing_nocalib_stg1_c1* and *t1e_without_borrowing_nocalib_stg1_c2* for no borrowing method (comparison 1 and comparison 2, respectively)

- *t1e_sa_stg1_c1* and *t1e_sa_stg1_c2* for simple average method (comparison 1 and comparison 2, respectively)

- *t1e_stg1_c1_eb* and *t1e_stg1_c2_eb* for Bayesian borrowing with Empirical Bayes power prior method (comparison 1 and comparison 2, respectively) 

**Power:**

- *power_without_borrowing_nocalib_stg1_c1* and *power_without_borrowing_nocalib_stg1_c2* for no borrowing method (comparison 1 and comparison 2, respectively)

- *power_sa_stg1_c1* and *power_sa_stg1_c2* for simple average method (comparison 1 and comparison 2, respectively)

- *power_withborrowing_stg1_c1_eb* and *power_withborrowing_stg1_c2_eb* for Bayesian borrowing with Empirical Bayes power prior method (comparison 1 and comparison 2, respectively)

**Other:**

- *bias_compar1* is the true differential placebo effect comparing placebo 1 to placebo 2

- for double dummy design operating characteristics just take the simple average method's operating characteristics when the differential placebo effect size is 0 units per year

- the 95\% confidence interval limits for treatment effect estimates are calculated from lines 300-328. The effect estimates and variances for each comparison at each stage under the null/alternative hypotheses for each method are used in these calculations so I don't detail here each variable. The Root Mean Square Error can be calculated from the effect estimates and variances used in these lines.
