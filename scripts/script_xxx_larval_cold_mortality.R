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
               betareg) 

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
raw_data <- readxl::read_xlsx("./data_raw/larval_cold_tolerance.xlsx")

# Check data 
head(raw_data)

###########################################################################
# Step 4: Clean data ------------------------------------------------------
###########################################################################

# Make the variable names easier to use and consistent
data <- raw_data %>%
  janitor::clean_names() %>%
  # Always convert <chr> columns to <fct> (factor)
  dplyr::mutate(across(where(is.character), as.factor))
head(data)

###########################################################################
# Step 4: Understand data/experimental design -----------------------------
###########################################################################

# How many different species were tested, how many reps per species, and how many cold exposure periods were there?
data %>%
  dplyr::count(species, replicate)

# Species - 3 species 
# Replicates - 4 per species
# Cold exposure periods - Usually 10 or 11 periods. 

###########################################################################
# Step 5: Visualise (exploratory data analysis) ---------------------------
###########################################################################

# Does the corrected mortality differ between species and over time (exposure periods)?
data %>%
  # First, calculate mean mortality across replicates for each species and 
  # exposure_period
  dplyr::group_by(species, exposure_period) %>%
  dplyr::summarise(
    mortality = mean(percent_corrected_mortality)) %>%
  # Now plot the graph
  ggplot(data = ., aes(x = exposure_period,
                       y = mortality,
                       group = species,
                       colour = species)) +
  geom_line(size = 1.2) +
  theme(legend.position = c(0.9, 0.2))

# You should be able to guess as to what your statistical analyses are going to say now:

# Is 'species' going to be a significant variable? 
# - I.e. Does corrected mortality differ between the different fly species, irrespective of exposure period?
# ANSWER: 

# Is 'exposure_period' going to be a significant variable? 
# - I.e. Does corrected mortality differ over time ('exposure_periods'), irrespective of which species is considered?
# ANSWER: 

# Does the effect of exposure period differ between species? 
# - I.e. Does the rate of mortality differ between fly species? 
# ANSWER: 

###########################################################################
# Step 3: Model building --------------------------------------------------
###########################################################################

# Before we fit the model, we need to decide on what an appropriate statistical distribution to use for this data. 
# - What do you think? 
# ANSWER:
# REASON: 

# Set global options to remove NA in model fitting otherwise 
# model selection dredge does not work later
options(na.action = "na.fail")

# Now fit global model
mod1 <- glm(percent_corrected_mortality ~ 
              # What variables could explain the response variable? 
              species + exposure_period,
            data = data,
            family = binomial(link = "logit"))

# When your Y variable is a proportion, R will give you a warning 
# 'non-integer #successes in a binomial glm!'
# - This is nothing to worry about. 

# Check model meets assumptions for GLM
mod1_res <- DHARMa::simulateResiduals(mod1)
plot(mod1_res)

# - LHS plot: Data conform to a bernoulli distribution. 
# - RHS plot: No significant departures from equal variances/linearity. 
# - Model diagnostics are good. 

car::Anova(mod1, type = "II")


