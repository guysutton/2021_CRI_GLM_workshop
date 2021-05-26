######################################################################
######################################################################
######################################################################
# - R workshop on linear models
# - Model #2: Full workthrough of modelling analysis 
#             (poisson/negative binomial GLM)
# - Script written: 20/05/2021
# - By: Guy F. Sutton
#   : Centre for Biological Control 
#   : Rhodes University
#   : email - g.sutton@ru.ac.za
######################################################################
######################################################################
######################################################################

######################################################################
# Aims:
# - 1. Understand when to use Poisson/Negative binomial distribution
# - 2. How to fit a Poisson GLM
# - 3. Understand how to evaluate Poisson GLM (model checking and diagnostics)
# - 4. Understand what an interaction term does
# - 5. How to fit a Negative Binomial GLM
# - 6. Understand how to evaluate Negative Binomial GLM (model checking and diagnostics)
######################################################################

######################################################################
# Step 1: Load packages ----------------------------------------------
######################################################################

# Load required packages
if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse, 
               tidyr,
               DHARMa,
               lmtest,
               glmmTMB,
               betareg,
               MuMIn,
               car,
               lmerTest,
               parameters,
               performance,
               lfe) 

# Change ggplot theme
# - This makes the graphs we will eventually make look pretty. 
# - Don't worry about this for now. 
theme_set(theme_classic() +
            theme(panel.border = element_rect(colour = "black", 
                                              fill = NA),
                  axis.text = element_text(colour = "black"),
                  axis.title.x = element_text(margin = unit(c(2, 0, 0, 0), 
                                                            "mm")),
                  axis.title.y = element_text(margin = unit(c(0, 4, 0, 0), 
                                                            "mm")),
                  legend.position = "none"))

###########################################################################
# Step 2: Import data -----------------------------------------------------
###########################################################################

# If your .csv uses commas 
raw_data <- readxl::read_xlsx("./data_raw/fly_trap_counts.xlsx")

# Check data 
head(raw_data)

###########################################################################
# Step 3: Clean data ------------------------------------------------------
###########################################################################

# Make the variable names easier to use and consistent
data <- raw_data %>%
  janitor::clean_names() %>%
  # Always convert <chr> columns to <fct> (factor)
  dplyr::mutate(across(where(is.character), as.factor)) %>%
  # Make fly densities into a count variable, just for illustration 
  dplyr::mutate(across(ftd_c_capitata_males_capilure:ftd_b_dorsalis_females_bl, ~ round(.x * 10), 0)) %>%
  # Let's just analyse the data for one fly group - ftd_c_capitata_males_capilure, 
  # and let's rename some of the columns to make our lives easier 
  dplyr::select(farm, 
                replicate, 
                bait = bait_density, 
                week = treatment_time, 
                fly_count = ftd_c_capitata_males_capilure)
head(data)

###########################################################################
# Step 4: Visualise (exploratory data analysis) ---------------------------
###########################################################################

# Visualise the distribution of our response variable (fly counts)
hist(data$fly_count)

# What does this tell us?:
# - Is our data likely to follow a normal (Gaussian distribution)?
# - ANSWER:

# - As such, which statistical distribution should we be fitting?
# - ANSWER: 

# Does the number of flies trapped differ between baits, farm and over time? 
data %>%
  dplyr::mutate(week = as.character(week)) %>%
  # Make week a numeric variable first 
  dplyr::mutate(week = dplyr::case_when(
    week == "Pre treatment" ~ "Week 0",
    TRUE ~ week)) %>%
  dplyr::mutate(week = readr::parse_number(week)) %>%
  dplyr::group_by(farm, bait, week) %>%
  dplyr::summarise(
    fly_count = mean(fly_count)) %>%
  dplyr::arrange(week) %>%
  ggplot(data = ., aes(x = week,
                       y = fly_count,
                       colour = bait,
                       group = bait)) +
  geom_line() +
  theme(legend.position = "right") +
  facet_wrap(~ farm, ncol = 2)

# You should be able to speculate as to what your statistical analyses 
# are going to say now:

# Is 'farm' going to be a significant variable? 
# - I.e. Do fly trap counts differ between farms, irrespective of time and baits?
# ANSWER: 

# Is 'bait' going to be a significant variable? 
# - I.e. Do fly trap counts differ between baits, irrespective of time and farm?
# ANSWER: 

# Does the effect of bait change over time? 
# ANSWER: 

# Does the effect of farm change over time? 
# ANSWER: 

# Does the effect of one or more baits depend on which farm it was deployed on? 
# ANSWER: 

###########################################################################
# Step 5: Model building --------------------------------------------------
###########################################################################

######################################
# Model #1: Standard Poisson GLM -----
######################################

# Let's now fit our first model:
# - Here, we are going to model the fly counts in traps as a function of 'farm', 
#   'bait' and 'week', using a Poisson distribution. 
# - We assume that the effect of bait, for example, is consistent across farms 
#   and weeks, and vice-versa. 

# Now fit simple additive model 
# - Poisson regression with a log link.
mod_poisson <- glm(fly_count ~ 
                     # What variables could explain the response variable? 
                     farm + bait + week,
                   data = data,
                   family = poisson)

# Was the model a good fit? 
mod1_res <- simulateResiduals(mod_poisson)
plot(mod1_res)

# - LHS: QQplot shows Poisson distribution was likely a bad choice.
# - RHS: Huge issues with homogeneity of variances. 

# Let's test for multicollinearity by calculating VIF's 
# - VIF > 2 is usually taken as evidence of multicollinearity.
# - If VIF > 2, remove that variable from the model, and refit. 
performance::check_collinearity(mod_poisson)

# - No issues with multicollinearity. 

# For a Poisson model, we can directly measure whether there was more variation 
# in the data than under the assumption of a Poisson distribution. 
# - This phenomen is called overdispersion. 
#    - Poisson models assume that mean and variance increase linearly
#    - I evalaute this by calculating an overdispersion statistic
#      from Gelman and Hill (2007)
#    - If the dispersion ratio is close to one = a Poisson model fits well to the data.
#      - Dispersion ratios larger than one indicate overdispersion, thus a negative binomial
#        model or another distribution that assumes larger variances might fit better to the data. 
#      - A p-value < .05 indicates overdispersion.
performance::check_overdispersion(mod_poisson)

# - Dispersion ratio is WAY larger than 1, indicating that there is much more
#   variation in our data than expected under the Poisson distribution. 

# Test for zero-inflation 
testZeroInflation(mod_poisson)

# - There are WAY too many zeroes in our data than expected under the Poisson distribution.

# - OVERALL: This model really sucks!!! 

######################################
# Model #2: Standard Poisson GLM, with interaction terms 
######################################

# Let's now fit our next model.
# - Here, we are going to model the fly counts in traps as a function of 'farm', 
#   'bait' and 'week', allowing the effect of farm and bait and week to vary within
#    each other, using a Poisson distribution. 
# - This model differs from the above, as we are allowing the effect of bait to vary
#   between farms and weeks, and vice-versa. 

# Fit model 
mod_poisson_int <- glm(fly_count ~ 
                         # What variables could explain the response variable? 
                         farm * bait * week,
                       data = data,
                       family = poisson)

# Was the model a good fit? 
mod1_res <- simulateResiduals(mod_poisson_int)
plot(mod1_res)

# - LHS: QQplot shows Poisson distribution was likely a bad choice.
# - RHS: Huge issues with homogeneity of variances. 

# Let's test for multicollinearity by calculating VIF's 
performance::check_collinearity(mod_poisson_int)

# - Issues. Eeeeek... 

# We test for overdispersion by calculating Gelman and Hill's overdispersion ratio. 
performance::check_overdispersion(mod_poisson_int)

# - Dispersion ratio is WAY larger than 1, indicating that there is much more
#   variation in our data than expected under the Poisson distribution. 

# Test for zero-inflation 
testZeroInflation(mod_poisson)

# - There are WAY too many zeroes in our data than expected under the Poisson distribution.

# - OVERALL: This model still really sucks!!! 

######################################
# Model #3: Standard negative binomial GLM
######################################

# Let's now fit our next model.
# - Here, we are going to model the fly counts in traps as a function of 'farm', 
#   'bait' and 'week', using a negative binomial distribution. 
#   - Remember from our earlier lecture that the negative binomial (NB) distribution
#     is an extension to the Poisson.
#     - The NB assumes alot more zeroes/smaller values, and much more variation.
#     - Remember, NB assumes variance > mean. 
# - Because there is no interaction term, the effect of bait, for example, 
#   is modelled as being consistent across farms and weeks, and vice-versa. 

# Fit model 
mod_nb <- glmmTMB::glmmTMB(fly_count ~ 
                             # What variables could explain the response variable? 
                             farm + bait + week,
                           data = data,
                           # Fit neg.bin (nbinom2) or quasi-poisson (nbinom1)
                           family = nbinom2,
                           # Use ML to fit model
                           REML = F)

# Was the model a good fit? 
mod1_res <- simulateResiduals(mod_nb)
plot(mod1_res)

# - LHS: QQplot shows NB distribution is a much better fit to the data.
# - RHS: Unfortunately, there is still an issue with homogeneity of variances. 

# Let's test for multicollinearity by calculating VIF's 
performance::check_collinearity(mod_nb)

# - No issues with multicollinearity. 

# Test for zero-inflation 
testZeroInflation(mod_nb)

# - The number of zeroes predicted by our model is in line with the 
#   expectation under the NB distribution. 

# - OVERALL: The NB model is definitely better, but still not a good fit.  

######################################
# Model #4: Standard negative binomial GLM, with interaction terms 
######################################

# Let's now fit our next model.
# - Here, we are going to model the fly counts in traps as a function of 'farm', 
#   'bait' and 'week', using a negative binomial distribution. 
# - This model differs from the above, as we are allowing the effect of bait to vary
#   between farms and weeks, and vice-versa. 

# Fit model 
mod_nb_int <- glmmTMB::glmmTMB(fly_count ~ 
                                 # What variables could explain the response variable? 
                                 farm * bait * week, 
                               data = data,
                               # Fit Negative Binomial distribution (nbinom2)
                               family = nbinom2, 
                               # Use ML to fit model
                               REML = F)

# What is model convergence? 
# - This warning tells us that R cannot find a numerical solution to the model
#   that we are trying to fit. 
#   - i.e. Simply, R cannot fit the model we are asking for.
#   - This usually happens when there is a mispecification in the model,
#     when the model is just too complex for the data, and/or 
#     the predictor variables are highly correlated (i.e. it cannot distinguish
#     between whether an effect should be attributed to farm/bait/week, in our example).

######################################
# Model #5: Standard negative binomial GLM, with interaction terms, different parameterisation
######################################

# Let's now fit our next model.
# - Here, we are going to model the fly counts in traps as a function of 'farm', 
#   'bait' and 'week', using a negative binomial distribution. 
# - Like the model above, we are allowing the effect of bait to vary
#   between farms and weeks, and vice-versa. 
# - However, we are adding the 'control' argument to tell R to estimate the parameters
#   of our model with slightly different settings. 

# Fit model 
mod_nb_int_control <- glmmTMB::glmmTMB(fly_count ~ 
                                         # What variables could explain the response variable? 
                                         farm * bait * week, 
                                       data = data,
                                       # Fit Negative Binomial distribution (nbinom2)
                                       family = nbinom2, 
                                       # Use ML to fit model
                                       REML = F,
                                       # Manually specify parameterisation
                                       control = 
                                         glmmTMBControl(optCtrl = list(iter.max = 1e3, 
                                                                       eval.max = 1e3)))


# Still no good.
# - The model we are trying to fix is way is too complex for the data we have. 
# - We are asking the model to identify whether fly counts differ for 
#   each week by farm by bait combination. 

# Also, we are not meeting our critical model assumption of data independence. 
# - Why are we failing this assumption? 
# - ANSWER: 

# What is the solution?
# - ANSWER: 






















######################################
# Model #6: Standard Poisson GLMM -----
######################################

# Let's now fit our next model.
# - Here, we are going to model the fly counts in traps as a function of 'farm' and 
#   'bait', while accounting for the non-independence between data points 
#   from repeated sampling in different weeks using a random intercept only, and
#   a Poisson distribution. 

# Fit model
mod_poisson_glmm <- glmer(fly_count ~ 
                       # What variables could explain the response variable? 
                       farm + bait +
                       # Random effects 
                       (1| week),             
                     # Family and link function
                     family = poisson(link = "log"),
                     # Data
                     data = data,
                     # Use an appropriate optimiser
                     control=glmerControl(optimizer="bobyqa",
                                          optCtrl=list(maxfun=2e5)),
                     verbose = FALSE,
                     nAGQ = 1)

# Was the model a good fit? 
mod1_res <- simulateResiduals(mod_poisson_glmm)
plot(mod1_res)

# Let's test for multicollinearity by calculating VIF's 
performance::check_collinearity(mod_poisson_glmm)

# We test for overdispersion by calculating Gelman and Hill's overdispersion ratio. 
performance::check_overdispersion(mod_poisson_glmm)

# Test for zero-inflation 
testZeroInflation(mod_poisson_glmm)

######################################
# Model #7: Standard Poisson GLMM, with interaction terms
######################################

# Let's now fit our next model.
# - Here, we are going to model the fly counts in traps as a function of 'farm' and 
#   'bait', while accounting for the non-independence between data points 
#   from repeated sampling in different weeks using a random intercept only, and
#   a Poisson distribution. 
#   - The only difference to the previous model is that we have now specified
#     an interaction term between 'farm' and 'bait', allowing the effect
#     of 'bait' to vary between farms. 
#   - NOTICE: This changes the hypothesis being tested!!! 

# Fit model
mod_poisson_glmm_int <- glmer(fly_count ~ 
                            # What variables could explain the response variable? 
                            farm * bait *
                            # Random effects 
                            (1| week),             
                          # Family and link function
                          family = poisson(link = "log"),
                          # Data
                          data = data,
                          # Use an appropriate optimiser
                          control=glmerControl(optimizer="bobyqa",
                                               optCtrl=list(maxfun=2e5)),
                          verbose = FALSE,
                          nAGQ = 1)

# Was the model a good fit? 
mod1_res <- simulateResiduals(mod_poisson_glmm_int)
plot(mod1_res)

# Let's test for multicollinearity by calculating VIF's 
performance::check_collinearity(mod_poisson_glmm_int)

# We test for overdispersion by calculating Gelman and Hill's overdispersion ratio. 
performance::check_overdispersion(mod_poisson_glmm_int)

# Test for zero-inflation 
testZeroInflation(mod_poisson_glmm_int)

######################################
# Model #8: Standard negative binomial GLMM
######################################

# Let's now fit our next model.
# - Here, we are going to model the fly counts in traps as a function of 'farm', 
#   'bait' and 'week', using a negative binomial distribution. 

# Fit model 
mod_nb_glmm <- glmmTMB::glmmTMB(fly_count ~ 
                             # What variables could explain the response variable? 
                             farm + bait +
                             # Random effects
                             (1| week),
                           data = data,
                           # Fit neg.bin (nbinom2) or quasi-poisson (nbinom1)
                           family = nbinom2,
                           # Use ML to fit model
                           REML = F)

# Was the model a good fit? 
mod1_res <- simulateResiduals(mod_nb_glmm)
plot(mod1_res)

# Let's test for multicollinearity by calculating VIF's 
performance::check_collinearity(mod_nb_glmm)

# Test for zero-inflation 
testZeroInflation(mod_nb_glmm)

######################################
# Model #9: Standard negative binomial GLMM, with interaction term
######################################

# Let's now fit our next model.
# - Here, we are going to model the fly counts in traps as a function of 'farm', 
#   'bait' and 'week', using a negative binomial distribution. 

# Fit model 
mod_nb_glmm_int <- glmmTMB::glmmTMB(fly_count ~ 
                                  # What variables could explain the response variable? 
                                  farm * bait +
                                  # Random effects
                                  (1| week),
                                data = data,
                                # Fit neg.bin (nbinom2) or quasi-poisson (nbinom1)
                                family = nbinom2,
                                # Use ML to fit model
                                REML = F)

# Was the model a good fit? 
mod1_res <- simulateResiduals(mod_nb_glmm_int)
plot(mod1_res)

# Let's test for multicollinearity by calculating VIF's 
performance::check_collinearity(mod_nb_glmm_int)

# Test for zero-inflation 
testZeroInflation(mod_nb_glmm_int)

######################################
# Model #9: Standard negative binomial GLMM, with interaction term and random slope (bait)
######################################

# Let's now fit our next model.
# - Here, we are going to model the fly counts in traps as a function of 'farm', 
#   'bait' and 'week', using a negative binomial distribution. 

# Fit model 
mod_nb_glmm_slope1 <- glmmTMB::glmmTMB(fly_count ~ 
                                      # What variables could explain the response variable? 
                                      farm * bait +
                                      # Random effects
                                      (1 + bait|week),
                                    data = data,
                                    # Fit neg.bin (nbinom2) or quasi-poisson (nbinom1)
                                    family = nbinom2,
                                    # Use ML to fit model
                                    REML = F)

######################################
# Model #9: Standard negative binomial GLMM, with interaction term and random slope (farm)
######################################

# Let's now fit our next model.
# - Here, we are going to model the fly counts in traps as a function of 'farm', 
#   'bait' and 'week', using a negative binomial distribution. 

# Fit model 
mod_nb_glmm_slope1 <- glmmTMB::glmmTMB(fly_count ~ 
                                         # What variables could explain the response variable? 
                                         farm * bait +
                                         # Random effects
                                         (1 + farm|week),
                                       data = data,
                                       # Fit neg.bin (nbinom2) or quasi-poisson (nbinom1)
                                       family = nbinom2,
                                       # Specify optimiser 
                                       control = 
                                         glmmTMBControl(optCtrl = list(iter.max = 1e3, 
                                                                       eval.max = 1e3)),
                                       # Use ML to fit model
                                       REML = F)

# Was the model a good fit? 
mod1_res <- simulateResiduals(mod_nb_glmm_slope1)
plot(mod1_res)

# Let's test for multicollinearity by calculating VIF's 
performance::check_collinearity(mod_nb_glmm_slope1)

# Test for zero-inflation 
testZeroInflation(mod_nb_glmm_slope1)

# Did the random slope for farm|week improve model fit 
# Perform likelihood ratio test 
anova(mod_nb_glmm_int, 
      mod_nb_glmm_slope1,
      test = "Chisq")

# - Yes. Random slope for farm|week provided a significantly better fit
#   than a random intercept for week alone. 
#   - Better fit indicated by lower AICc, lower logLik and P < 0.05. 
#   - This indicates that allowing the effect of time on fly counts 
#     to vary between farms was important. 

######################################
# Model #9: Standard negative binomial GLMM, with no interaction term and random slope (farm)
######################################

# Let's now fit our next model.
# - Here, we are going to model the fly counts in traps as a function of 'farm', 
#   'bait' and 'week', using a negative binomial distribution. 

# Fit model 
mod_nb_glmm_slope1_add <- glmmTMB::glmmTMB(fly_count ~ 
                                         # What variables could explain the response variable? 
                                         farm + bait +
                                         # Random effects
                                         (1 + farm|week),
                                       data = data,
                                       # Fit neg.bin (nbinom2) or quasi-poisson (nbinom1)
                                       family = nbinom2,
                                       # Specify optimiser 
                                       control = 
                                         glmmTMBControl(optCtrl = list(iter.max = 1e3, 
                                                                       eval.max = 1e3)),
                                       # Use ML to fit model
                                       REML = F)

# Was the model a good fit? 
mod1_res <- simulateResiduals(mod_nb_glmm_slope1_add)
plot(mod1_res)

# Let's test for multicollinearity by calculating VIF's 
performance::check_collinearity(mod_nb_glmm_slope1_add)

# Test for zero-inflation 
testZeroInflation(mod_nb_glmm_slope1_add)

# Did the removing the interaction between farm * bait improve model fit 
# Perform likelihood ratio test 
anova(mod_nb_glmm_slope1,
      mod_nb_glmm_slope1_add,
      test = "Chisq")

# - No. The model with the interaction term was supposedly a better fit, 
#   however, the variables were colinear, so we can't really trust that model. 

# - We can be quite happy with using the model without the interaction term,
#   it was a good fit to the data:
#   - QQplot showed the distribution used was okay - the significant KS test
#     is not too bad... GLMM's are able to compensate for small deviations 
#     from expectation... 
#   - Residual vs fitted plot showed no evidence for unequal variances. 
#   - VIF's showed no evidence for multicolinearity between predictors. 
#   - No evidence for too many zeroes in the data than assumed by the distribution
#     we used (negative binomial distribution). 

######################################################################
# Step 6: Interpret model output -------------------------------------
######################################################################

######################################
# Statistical significance -----------
######################################

# Let's see if any of our predictor variables were statistically significant
# predictors of fly trap counts, by calculating Wald's Chisq tests, using
# type II sum-of-squares
car::Anova(mod_nb_glmm_slope1_add, 
           type = "III",
           test.statistic = "Chisq")

# - Our results indicate that 'farm' and 'bait' were statistically significant.

#######################################
# Post-hoc tests ---------------------
######################################

# Perform post-hoc 
emmeans(mod_nb_glmm_slope1_add, 
        specs = pairwise ~ bait | farm,
        type = "response")

###########################################################################
# Alternative model inference for prediction ------------------------------
###########################################################################

# Below, I perform a model selection procedure to identify potentially
# important predictor variables from the models specified above. 
# - Model selection was performed using AICc and Akaike weights. 
aic_vals <- MuMIn::dredge(mod_nb_glmm_slope1_add)
aic_vals

# Save full AIC table
write.csv(aic_vals, "./results/full_aic_model_selection.csv")

# Subset only the plausible set of models (i.e. those models scoring
# within 2 AICc points of the top-performing model). 
plausible_mods <- subset(aic_vals, delta < 2)
plausible_mods

# Save plausible model AIC table
write.csv(plausible_mods, "./results/plausible_model_selection.csv")

# How many plausible models were retained? 
nrow(plausible_mods)

# Calculate cumulative support for each variable
# Sum of Akaike weights 
# Here, I was interested in which variables:
# (1) Retained in the most models (N containing models). 
#     - The more plausible models containing each model term, 
#       the higher the importance of that term.
# (2) Had the highest sum of weights. Weights closer to 1 indicate
#     higher importance, weights closer to 0 indicate lesser importance. 
importance_vals <- MuMIn::sw(plausible_mods)
importance_vals







