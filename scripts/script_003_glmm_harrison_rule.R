######################################################################
######################################################################
######################################################################
# - R workshop on linear models
# - Model #4: A gentle introduction to LMM/GLMM's 
# - Script written: 25/05/2021
# - By: Guy F. Sutton
#   : Centre for Biological Control 
#   : Rhodes University
#   : email - g.sutton@ru.ac.za
######################################################################
######################################################################
######################################################################

######################################################################
# Aims:
# - 1. Specify a LMM/GLMM model in R
# - 2. Understand how to specify random intercept only, 
#      random slope only and random intercept + random slope LMM. 
# - 3. How to statistically evaluate which model provides the best fit,
#      using AICc and Likelihood Ratio Tests. 
# - 4. Assess parameter/variable statistical significance once best model
#      is identified 
# - 5. Understand results obtained from a LMM (note: this only applies 
#      to Gaussian LMM - interpretation is different for GLMM [e.g. Poisson])
######################################################################

# Load required packages 
if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse, 
               lme4) 

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
                                                            "mm"))))

######################################################################
# Step 1: Generate simulated data (*NO NEED TO UNDERSTAND) -----------
######################################################################

# - Code and ideas from Michael Freeman (mfviz.com/hierarchical-models)

# Parameters for generating insect body mass data
site <- c('Grahamstown', 'Port Elizabeth', 'East London', 'Queenstown', 'Cradock')
base_mass <- c(40, 50, 60, 70, 80)
increase_mass <- c(9, 5, 5, 7, 5)
n_site <- 20
n_total <- n_site * length(site)

# Generate dataframe of body_mass and (random) sites
ids <- 1:n_total
site <- rep(site, n_site)
body_length <- floor(runif(n_total, 0, 10))
bases <- rep(base_mass, n_site) * runif(n_total, .9, 1.1) # noise
raises <- rep(increase_mass, n_site) * runif(n_total, .9, 1.1) # noise
df <- data.frame(ids, site, bases, body_length, raises)

# Generate body masses (base + experience * raise)
df <- df %>% 
  dplyr::mutate(body_mass = bases + body_length * raises,
                body_mass = (body_mass / 10) - 3) %>%
  dplyr::select(site, 
                plant_stem_diameter = body_length,
                insect_body_length = body_mass) %>%
  tibble::as_tibble()

###########################################################################
# Step 2: Visualise (exploratory data analysis) ---------------------------
###########################################################################

# Quick inspection of data
head(df)
dplyr::glimpse(df)

# Make plot 
df %>%
  ggplot(data = ., aes(x = plant_stem_diameter,
                       y = insect_body_length,
                       colour = site)) +
  geom_point() +
  labs(x = "Stem diameter (cm)",
       y = "Insect body length (mm)",
       colour = "Site")

###########################################################################
# Step 3: Model building --------------------------------------------------
###########################################################################

# Null model (standard Gaussian GLM)
m0 <- glm(insect_body_length ~ 
            plant_stem_diameter, 
          family = gaussian(link = "identity"),
          data = df)

# Model with random intercept
m1 <- lmer(insect_body_length ~ 
             plant_stem_diameter + (1|site), 
           data = df,
           REML = F)

# Model with random slope, no intercept 
m2 <- lmer(insect_body_length ~ 
             plant_stem_diameter + 
             (0 + plant_stem_diameter|site), 
           data = df,
           REML = F)

# Model with random slope and random intercept 
m3 <- lmer(insect_body_length ~ 
             plant_stem_diameter + 
             (1 + plant_stem_diameter|site), 
           data = df,
           REML = F)

###########################################################################
# Step 4: Model selection --------------------------------------------------
###########################################################################

# But, which model was best?
# - Use model diagnostics to check the fit of your data
#   - We did this yesterday, so not going to cover this again. 

# We can use Akaike's Information Criterion (AICc) to select see 
# which is the top-performing model (explains the most variation). 
# - AICc penalises models for being too complex.
# - It balances model complexity vs model performance. 
#   - The model that explains the model variation = lowest AICc
AICc(m0)
AICc(m1)
AICc(m2)
AICc(m3)

# - m3 is the best-performing model with the lowest AICc, by far. 

# We can also formally test for the performance of one model vs another using
# the likelihood ratio test. 
# - Does random slope & random intercept (m3) provide a better fit than 
#   a random slope only (m2)? 
anova(m2, 
      m3, 
      test = "Chisq")

# - Yes. The model with the random intercept & random slope (m3) provided
#   a significantly better fit (X2 = 126.33, P < 0.001) than the 
#   model with the random slope only (m2).

# Likelihood ratio test
# - Does random slope & random intercept (m3) provide a better fit than 
#   a random intercept only (m1)? 
anova(m1, 
      m3, 
      test = "Chisq")

# - Yes. The model with the random intercept & random slope (m3) provided
#   a significantly better fit (X2 = 55.83, P < 0.001) than the 
#   model with the random intercept only (m1).

###########################################################################
# Step 5: Model inference --------------------------------------------------
###########################################################################

# Now that you have your top model:
# Assessing parameter significance 
car::Anova(m3,
           test = "II",
           test.statistic = "Chisq")

# - After accounting for variation between sites, there is a statistically 
#   significant relationship between plant stem diameter and insect body length 
#   (X2 = 40.41, d.f. = 1, P < 0.001). 

# Extract fixed effect parameter estimate 
fixef(m3)

# - For every 1 unit increase in x (i.e. 1mm increase in plant stem diameter),
#   insect body length, on average, increases by 0.63 mm. 

# And calculate 95% confidence intervals
confint(m3, 
        level = 0.95,
        method = "boot")

# - The effect of increasing plant stem diameter by 1 mm will result in 
#   an increase of 0.43 to 0.82 cm in insect body length (95% CI). 

# - TAKEHOME: Really good support for Harrison's rule in our study system. 
#             (Rule = Parasite body size is positively correlated with host size).  
