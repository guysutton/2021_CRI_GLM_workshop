######################################################################
######################################################################
######################################################################
# - R workshop on linear models
# - Model #1: Full workthrough of modelling analysis (logit/probit GLM)
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
# - 1. Import .csv and .xlsx files
# - 2. Make sure data imported correctly
# - 3. Tidy and visualise raw data
# - 4. Specify a logistic GLM
# - 5. Evaluate the fit of your model
# - 6. Interpret output from model
# - 7. Plot predictions from model
# - 8. Write-up your model results 
# - 9. Understand the difference between logit vs probit GLM
# - 10. Specify a probit GLM
# - 11. Calculate LD50/LD90
######################################################################

######################################################################
# Step 1: Load packages ----------------------------------------------
######################################################################

# Load required packages
if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse, 
               tidyr,
               DHARMa, 
               emmeans) 

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

######################################
# Import dataset (example using .xlsx)
######################################

# Import .xlsx file
raw_data <- readxl::read_excel("./data_raw/hypothetical_data_serfontein.xlsx")

# Check to make sure data imported properly
# - Check the first 6 rows of the data 
head(raw_data)

# Have a look at the structure of the data 
# - Check that the numbers are numeric (not characters), and vice-versa, 
#   and that there are as many columns and rows as expected from your .xlsx
#   or .csv file 
glimpse(raw_data)

######################################
# Import dataset (example using .csv)
######################################

# Import dataset (example using .csv)
# - Depending on how your PC works, some PC's encode .csv's with 
#   actual commas (as you would expect), but my PC uses a semi-colon 
#   (go figure!).

# Option #1: If your .csv uses semi-colons 
raw_data <- readr::read_csv2("./data_raw/hypothetical_data_serfontein.csv")

# Did that work?
head(raw_data)
dplyr::glimpse(raw_data)

# - Nope. Only one column in our dataset. Eeeeeek! 

# Option #2: If your .csv uses commas 
raw_data <- readr::read_csv("./data_raw/hypothetical_data_serfontein.csv") 

# Did that work?
head(raw_data)
dplyr::glimpse(raw_data)

# - Looks good, but all the column names are very strange. 

###########################################################################
# Step 3: Clean data-------------------------------------------------------
###########################################################################

# Make the variable names easier to use and consistent
data <- raw_data %>%
  janitor::clean_names()
head(data)

# But, we don't have a column containing proportion of flies that died. 

# Add a column with total no. of deaths after 72 hours (i.e. sum deaths 
# after 24, 48 and 72 hours)
data <- data %>%
  # Keep only the columns of interest
  # - Using the : keeps all the columns from test_rep to total_no_flies
  dplyr::select(test_rep:total_no_flies) %>%
  # Performs calculations for each row separately
  dplyr::rowwise() %>%
  # Now sum deaths in a new column
  dplyr::mutate(total_deaths_72 = sum(x24_h_deaths,
                                      x48_h_deaths,
                                      x72_h_deaths)) %>% 
  # Some replicates somehow had > 10 flies dead, maybe data transcription error?
  dplyr::mutate(total_deaths_72 = dplyr::case_when(
    # If total deaths > 10, then make it 10
    total_deaths_72 > 10 ~ 10,
    # If total deaths < 10, then keep the original value 
    TRUE ~ total_deaths_72))
head(data)

# Let's just select the columns of interest 
data <- data %>%
  dplyr::select(test_rep:sex, total_no_flies, total_deaths_72)
head(data)

# Usually a good idea to change columns that are character vectors (<chr>)
# into factor vectors (<fct>).
# - This will make your life much easier down the line. 
data <- data %>%
  dplyr::mutate(across(where(is.character), as.factor))
head(data)

###########################################################################
# Step 4: Visualise (exploratory data analysis) ---------------------------
###########################################################################

# Does the number of dead flies differ between different baits?
data %>%
  ggplot(data = ., aes(x = bait,
                       y = total_deaths_72)) +
  geom_boxplot()

# Does the number of dead flies differ between different trap ages?
data %>%
  ggplot(data = ., aes(x = trap_age,
                       y = total_deaths_72)) +
  geom_boxplot()

# Does the number of dead flies differ between different fly species?
data %>%
  ggplot(data = ., aes(x = fly_type,
                       y = total_deaths_72)) +
  geom_boxplot()

# Does the number of dead flies differ between males and female?
data %>%
  ggplot(data = ., aes(x = sex,
                       y = total_deaths_72)) +
  geom_boxplot()

# We can also look at two variables at once 
# e.g. Bait and trap age
data %>%
  ggplot(data = ., aes(x = bait,
                       y = total_deaths_72,
                       fill = trap_age)) +
  geom_boxplot() +
  # Add a legend
  theme(legend.position = "right")

# We can also look at three variables at once 
# e.g. Bait, fly species and sex
data %>%
  ggplot(data = ., aes(x = bait,
                       y = total_deaths_72,
                       fill = sex)) +
  geom_boxplot() +
  # Add a legend
  theme(legend.position = "right") +
  # Add a different panel for each fly species
  facet_wrap(~ fly_type, nrow = 2)

# For completeness, let's make the 'Control' boxes appear at the far left 
# of each graph. This improves readability. 
data %>%
  # Change the order of the x-axis bars 
  dplyr::mutate(bait = forcats::fct_relevel(bait, 
                                            "Control", 
                                            "CeraTrap",
                                            "GF120",
                                            "M3", 
                                            "MM")) %>%
  ggplot(data = ., aes(x = bait,
                       y = total_deaths_72,
                       fill = sex)) +
  geom_boxplot() +
  # Add a legend
  theme(legend.position = "right") +
  # Add a different panel for each fly species
  facet_wrap(~ fly_type, nrow = 2) +
  # Add nice axis and legend titles
  labs(x = "Bait type",
       y = "Total no. dead flies after 72h",
       fill = "Sex") +
  # Manually specify the y-axis breaks 
  scale_y_continuous(limits = c(0, 10),
                     breaks = seq(0, 10, 2))

# We can change the focus of our graphs easily.
# e.g. Instead of focusing on comparisons for each species (species per panel),
# let's rather focus on the different baits (bait per panel)
data %>%
  dplyr::mutate(bait = forcats::fct_relevel(bait, 
                                            "Control", 
                                            "CeraTrap",
                                            "GF120",
                                            "M3", 
                                            "MM")) %>%
  ggplot(data = ., aes(x = fly_type,
                       y = total_deaths_72,
                       fill = sex)) +
  geom_boxplot() +
  # Add a legend
  theme(legend.position = "right") +
  # Add a different panel for each fly species
  facet_wrap(~ bait, nrow = 1) +
  # Add nice axis and legend titles
  labs(x = "Fly species",
       y = "Total no. dead flies after 72h",
       fill = "Sex") +
  scale_y_continuous(limits = c(0, 10),
                     breaks = seq(0, 10, 2))

###########################################################################
# Step 5: Specify research question and/or hypotheses/expectations --------
###########################################################################

# What is our research question:
# ANSWER: DOES THE NUMBER OF FLIES THAT DIED DIFFER BETWEEN TRAP BAITS, 
#         TRAP AGES, FLY SPECIES AND/OR FLY GENDER? 

# You should be able to speculate as to what your statistical analyses 
# are going to say now:

# Is 'bait' going to be a significant variable? 
# - I.e. Does the proportion of flies that died differ between baits?
# ANSWER: 

# Is 'trap_age' going to be a significant variable? 
# - I.e. Does the proportion of flies that died differ between trap ages (fresh/aged)?
# ANSWER: 

# Is 'fly_type' going to be a significant variable? 
# - I.e. Does the proportion of flies that died differ between fly species?
# ANSWER: 

# Is 'sex' going to be a significant variable? 
# - I.e. Does the proportion of flies that died differ between sex (male/female)?
# ANSWER: 

###########################################################################
# Step 6: Model building --------------------------------------------------
###########################################################################

# Let's now fit our first model.
# - Here, we are going to model the proportion of flies that died (mortality) 
#   as a function of 'bait', 'trap_age', "fly_type' and 'sex, using a standard binomial GLM. 

# Performing a GLM in R is very simple. 
# - We will use the 'glm' function - which stands for general linear model.
#   - We have to specify a model formula (this style is consistent across 
#     most models in R, so get used to it).
#     - First, we specify our response variable (the thing we are trying to predict),
#       and then our explanatory variable(s). 
#     - Second, we have to tell R where to look for those variables. 
#     - Thirdly, we have to tell R which type of glm we want to fit.
#       - E.g. Gaussian (normal distiribution), Poisson, Binomial, ect... 
#     - Lastly, we store the linear regression in an object called 'mod1' 

# Now fit global model
# - Binomial regression with a logit link (logistic regression)
mod1 <- glm(cbind(total_deaths_72, total_no_flies - total_deaths_72) ~ 
              # What variables could explain the response variable? 
              bait + trap_age + fly_type + sex,
            data = data,
            family = binomial(link = "logit"))

# See how little time and code it took to actually run the model? 

# Wait. Before, we start looking at our results, it is always a good idea
# to see how good the model fit was to your data.
# - You don't want to look at the results and find a really cool result, and
#   then later check and find that the model was not a good fit. 

###########################################################################
# Step 7: Evaluate model fit ----------------------------------------------
###########################################################################

# Before looking at the outcome of our model,
# we should first check whether the model we specified is 
# a good fit to the data or not. 
# We do this by performing model diagnostics. 

# Basically, here we are going to check whether our 
# data and model meet the assumptions required to run a linear model.

# 1. Independence of data points
#    - This cannot be checked per se. 
#    - This is a property of your experimental design
#    - Make sure you design your experiment properly (pseudo-replication!!!). 

# 2. Was the statistical distribution (family) the correct choice?
#    - To evaluate, we calculate a QQplot of model residuals. 
#    - The QQplot is the plot on the left-hand side of this plot we produce below. 
#    - Residuals represent the difference between the observed data and the predicted 
#      value (obs - exp).
#    - This assumption is here to determine if the statistical distribution
#      we are using for our model is appropriate for the data. 
#    - If the model is a good fit, then the residuals should be fit approximately 
#      along the red 1:1 line and the KS test p-value should be > 0.05. 
#       - The KS test stands for a Kolmogorov-Smirnoff test which compares
#         the model residuals against the statistical distribution we assumed.
#       - E.g. If we fit family = binomial, the test will compare whether our 
#              model residuals approximate the binomial distribution, 
#              but if we fit family = poisson, the test will compare whether our 
#              model residuals approximate the poisson distribution.
#    - We don't want to see any systematic patterns in residuals, so
#      the points should fall approximately on the dashed line. 
mod1_res <- simulateResiduals(mod1)
plot(mod1_res)

# - Here, we can see that our points fall roughly along the red line (this is good),
#   and the KS test P > 0.05. 
#   - As such, the residuals from our model approximate the binomial distribution
#     relatively well.
#   - Therefore, the choice of statistical distribution (family = binomial)
#     appeared to be an appropriate choice. 

# 3. Homogeneity of variance
#   - Testing the assumption of equality of variances. 
#   - We must look at the residuals vs fitted plot on the right-hand side 
#     of our plots. 
#     - We want to see the three bold lines fall approximately on the y = 0.25,
#       0.50 and 0.75 lines.
#     - If there are issues with homogeneity of variances, there bold lines
#       will be coloured red automatically (bad), instead of black (good).
#     - Also, the title of the plot makes it very clear if there were problems or not. 
mod1_res <- simulateResiduals(mod1)
plot(mod1_res)

# - Here, we can see no significant departures from equal variances.
#   - We can be quite happy that our assumption of homogeneity of variances
#     was met. 

# - TAKE HOME: Our model was a good fit, so let's now look at the results 
#              from our model.  

######################################################################
# Step 8: Interpret model output -------------------------------------
######################################################################

######################################
# Statistical significance -----------
######################################

# Print summary of results 
# - Remember, our GLM is stored in the 'mod1' object.
# - Anyone who has run models in R before is probably familiar with 
#   extracting the model output by summary(<model_name>). 
# - I encourage you not to get too familiar with this output.
# - More times than not, these are NOT the results you want to look at. 
#   - This is due to the way in which R calculates the p-values. 
#   - The summary(...) output produces type I sum of squares, which calculates
#     p-values based on the order in which terms are added to the model.
#   - This is not great (but makes no difference when you only have
#     one predictor variable in the model).
#   - Also, all the p-values are comparisons to a baseline level
#     e.g. Here, R is comparing each row to the baseline level:
#     - Bait = Control, Trap age = Fresh, fly_type = B. dorsalis and sex = Female. 
#     - Not exactly a comparison that makes any sense, is it? 
summary(mod1)

# We need to specify how we want R to calculate p-values. 
# - When you have 1 term in the model, type I SS is fine. 
# - When you have 2 terms in the model (but no interaction),
#   type II SS is fine. 
# - When you have 2 terms in the model, and an interaction term, 
#   you must use type III SS. 
# - We won't go into more details now, but get into the habit of 
#   manually specifying the SS required. 

# Let's see if any of our predictor variables were statistically significant
# predictors of fly mortality, by performing likelihood ratio tests (LRT), using
# type II sum-of-squares
car::Anova(mod1, 
           type = "II",
           test.statistic = "LR")

# - Our results indicate that 'bait' and 'trap-age' were statistically significant,
#   but not fly_type and sex. 
#   - We will discuss how to report the statistics and results for a publication/thesis
#     later. 

######################################
# Post-hoc tests ---------------------
######################################

# Post-hoc analysis
# - Because fly-type and sex were not significant, we average our 
#   predictions over each level within fly_type and sex, and then
#   extract the post-hoc pairwise comparisons... 
post_hoc_results <- emmeans(
  # Give the variable name where your model is stored
  mod1, 
  # We want pairwise comparisons for each variable after ~ ... 
  specs = pairwise ~ bait + trap_age, type = "response")
print(post_hoc_results)

######################################################################
# Step 9: Make a figure ----------------------------------------------
######################################################################

# Let's plot the proportion of dead flies by bait and trap-age,
# averaged across 'fly-type' and 'sex'.
data %>%
  # First, we must convert the counts of dead flies into proportions 
  dplyr::group_by(bait, trap_age) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(
    prop_dead = total_deaths_72 / total_no_flies) %>%
  # Change the order of the x-axis bars 
  dplyr::mutate(bait = forcats::fct_relevel(bait, 
                                            "Control", 
                                            "CeraTrap",
                                            "GF120",
                                            "M3", 
                                            "MM")) %>%
  dplyr::ungroup() %>%
  # Now make the plot
  ggplot(data = ., aes(x = bait,
                       y = prop_dead,
                       fill = trap_age)) +
  # Make it a boxplot
  geom_boxplot() +
  # Add a legend
  theme(legend.position = "right") +
  # Make the bars grayscale
  scale_fill_manual(values = c("grey40", "grey80")) +
  # Add nice axis and legend titles
  labs(x = "Bait type",
       y = "Proportion of flies that died after 72h",
       fill = "Trap Age") +
  # Change the limits of the y-axis
  scale_y_continuous(limits = c(0, 1.1),
                     breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1.0)) +
  # Add significance stars (use the emmeans output/tukey)
  # Use x = "..." to indicate which group
  # x = 1 means first group,
  # x = 2 means second group... 
  # Manually play around with y-values
  # *** THE LABELS ARE NOT CORRECT - you can go make them correct, if you 
  #     want - this process takes a long-time to figure out the letters.
  annotate("text", x = 0.8, y = 0.5, label = "a") +
  annotate("text", x = 1.2, y = 0.5, label = "a") +
  annotate("text", x = 1.8, y = 1.05, label = "bd") + 
  annotate("text", x = 2.2, y = 1.05, label = "b") + 
  annotate("text", x = 2.8, y = 1.05, label = "b") + 
  annotate("text", x = 3.2, y = 1.05, label = "b") +
  annotate("text", x = 3.8, y = 0.85, label = "c") + 
  annotate("text", x = 4.2, y = 1.05, label = "bd") +
  annotate("text", x = 4.8, y = 1.05, label = "b") + 
  annotate("text", x = 5.2, y = 1.05, label = "e")

######################################################################
# Step 10: Reporting the statistics -----------------------------------
######################################################################

# Bait type (X2 = 294.42, d.f. = 4, P < 0.001) and trap age (X2 = 10.92, 
# d.f. = 1, P < 0.001) were statistically significant predictors of fly mortality
# (Table 1). There was no evidence for an effect of fly_type 
# (X2 = 0.76, d.f. = 3, P = 0.859) or sex (X2 = 294.42, d.f. = 1, P = 0.558)
# on fly mortality. 

# Then report results from your figure.
# - e.g. All baits produced higher fly mortality rates than the controls.
# - e.g. GF120 appeared to consistently result in the highest fly mortality rates.
# - e.g. In general, aged traps tended to produce slightly higher fly mortality rates
#        than fresh traps... 

######################################################################
# ADDITIONAL CREDITS: Probit analysis --------------------------------
######################################################################

# It seems that citrus research really likes the use of probit regression
# over the more common logistic regression (like we just performed above). 
# - The difference between the two is the shape of the curve fitted to the data. 

# How do we fit a probit GLM?
# - We simply change the link function 
mod_probit <- glm(cbind(total_deaths_72, total_no_flies - total_deaths_72) ~ 
              # What variables could explain the response variable? 
              bait + trap_age + fly_type + sex,
            data = data,
            # The only change is to make link = "probit" instead of link = "logit"
            family = binomial(link = "probit"))

# Check model meets assumptions for GLM
mod_probit_res <- simulateResiduals(mod_probit)
plot(mod_probit_res)

# Formal statistical test to evaluate whether logit or probit model was better? 
# - Likelihood ratio test 
anova(mod1, 
      mod_probit, 
      test = "Chisq")

# - ANSWER: No significant improvement of fit. Not better, but also not worse. 

####################################################
# Do the logit and probit produce different results? 
####################################################

# Recalculate results from logit model 
car::Anova(mod1, 
           type = "II",
           test.statistic = "LR")

# Now calculate results from the probit model
car::Anova(mod_probit, 
           type = "II",
           test.statistic = "LR")

# ANSWER: Not at all. The difference between the logit and probit is 
#         usually very little, if anything at all. 

######################################################################
# ADDITIONAL CREDITS: Extract LC50 and LC90 --------------------------
######################################################################

# Load R library
library(MASS)

# Look at the Lamprey toxicity data 
lamprey_tox <- lamprey_tox[lamprey_tox$nominal_dose != 0, ]
head(lamprey_tox)

# What do we have here?
# Lamprey mortality against different chemical concentrations (toxicity study)
# - Look at the doses applied
hist(lamprey_tox$dose)

# How do we fit this GLM?
# - If we want to calculate LC50 and LC90 of chemical dose on lamprey mortality. 
# - Fit Lamprey toxicity model 
mod_probit <- glm(cbind(response, total - response) ~ 
                    # What variables could explain the response variable? 
                    dose,
                  data = lamprey_tox,
                  # Let's fit a probit link function 
                  family = binomial(link = "probit"))

# Calculate LC50 and LC90
lc_data <- MASS::dose.p(mod_probit, 
                        p = c(0.5, 0.9))
lc_data

# Plot indicating LC50 and LC90
p1 <- ggplot(data = lamprey_tox,
            aes(x = dose, 
                y = (response / total))) +
 geom_point(size = 1) +
 geom_smooth(method = "glm",
             method.args = list(family = binomial(link = "probit")),
             aes(weight = total), 
             colour = "black", 
             fill = c("grey60"),
             fullrange=TRUE,
             se = TRUE) +
  # Add LC50 lines
  geom_segment(aes(x = 1.59 , 
                   y = 0, 
                   xend = 1.59, 
                   yend = 0.50),
               linetype = "dashed") +
  geom_segment(aes(x = 0.75, 
                   y = 0.50, 
                   xend = 1.59, 
                   yend = 0.50),
               linetype = "dashed") +  
  # Add LC90 lines
  geom_segment(aes(x = 5.81 , 
                   y = 0, 
                   xend = 5.81, 
                   yend = 0.90),
               linetype = "dotted") +
  geom_segment(aes(x = 0.75, 
                   y = 0.90, 
                   xend = 5.81, 
                   yend = 0.90),
               linetype = "dotted") +
  labs(x = "Pesticide (concentration)",
       y = "Proportion of population that died")
p1
