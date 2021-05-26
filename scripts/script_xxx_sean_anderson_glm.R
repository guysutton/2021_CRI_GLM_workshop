logit <- function(x) qlogis(x)
inverse_logit <- function(x) plogis(x)

xx <- seq(0, 5, length.out = 200)
plot(xx, exp(xx), type = "l")

A log link will make the data look linear:
  

plot(xx, log(exp(xx)), type = "l")

And if this is the curve we ultimately want to fit:
  
  ```{r}
xx <- seq(-5, 5, length.out = 200)
plot(xx, inverse_logit(xx), type = "l")
```

Then we can make it linear by applying the logit link:
  
```{r}
plot(xx, logit(inverse_logit(xx)), type = "l")






# Link functions 

Remember that there are two components to a GLM: a link function that describes what is happening with the mean of the data, and an error distribution that describes the variability around the mean. 

Put another way, the "link" is the transformation you need to apply to make the (mean of the) response data linear with respect to the predictors.

The error distribution describes the spread of the data around the raw untransformed mean.

The two most common link functions, and the only two we are going to work with in this workshop, are the log and logit links.

Let's look at those now. So if we want to fit a curve that looks like this:

```{r}
xx <- seq(0, 5, length.out = 200)
plot(xx, exp(xx), type = "l")
```

A log link will make the data look linear:

```{r}
plot(xx, log(exp(xx)), type = "l")
```

And if this is the curve we ultimately want to fit:

```{r}
xx <- seq(-5, 5, length.out = 200)
plot(xx, inverse_logit(xx), type = "l")
```

Then we can make it linear by applying the logit link:

```{r}
plot(xx, logit(inverse_logit(xx)), type = "l")
```

There are many ways you can specify a distribution family and link function in R. I'm going to try and be consistent and specify them like this `family(link = "link")`.

set.seed(111)
N <- 200
x <- runif(N, -1, 1)
a <- 0.5
b <- 1.3
d <- data_frame(x = x)
d

y_true <- exp(a + b * x)
shape <- 8
y <- rgamma(N, rate = shape / y_true, shape = shape)
d$y <- y
plot(x, y);lines(sort(x), y_true[order(x)])

(m_gamma <- glm(y ~ x, family = 
                  Gamma(link = "log"), data = d)) # exercise
ggeffects::ggpredict(m_gamma, "x", full.data = TRUE) %>%
  ggplot(aes(x, predicted)) +
  geom_line(colour = "red") + 
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .2) +
  geom_point(data = d, aes(y = y))

# Poisson GLM
y_true <- exp(a + b * x)
y <- rpois(N, lambda = y_true)
d$y <- y
plot(x, y);lines(sort(x), y_true[order(x)])

(m_poisson <- glm(y ~ x, family = 
                    poisson(link = "log"), data = d)) # exercise
ggeffects::ggpredict(m_poisson, "x", full.data = TRUE) %>%
  ggplot(aes(x, predicted)) +
  geom_line(colour = "red") + 
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .2) +
  geom_point(data = d, aes(y = y))

# Negative binomial, log link

#The negative binomial distribution with a log link can also model count data but allows the #variance to grow as a quadratic  of the mean. In real data sets, it's probably more common to see #the negative binomial than the Poisson.

y_true <- exp(a + b * x)
y <- MASS::rnegbin(N, mu = y_true, theta = 0.6)
d$y <- y
plot(x, y);lines(sort(x), y_true[order(x)])

Notice the much larger values on the right side of the graph.

(Also note that there is another common parameterization of the negative binomial which allows the variance to grow linearly with the mean.)

We have to use a special function to fit the negative binomial GLM in R

(m_nb <- MASS::glm.nb(y ~ x, data = d))
ggeffects::ggpredict(m_nb, "x", full.data = TRUE) %>%
  ggplot(aes(x, predicted)) +
  geom_line(colour = "red") + 
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .2) +
  geom_point(data = d, aes(y = y))

## Binomial, logit link

We can use a binomial response and logit link if we have response data represented by 0s and 1s. This is commonly referred to as logistic regression. 

```{r}
y_linear <- a + b * x
prob_true <- inverse_logit(y_linear)
y <- rbinom(N, 1, prob_true)
d$y <- y
plot(x, jitter(y, 0.1));lines(sort(x), prob_true[order(x)])
```

What does the true probabilities line indicate in the above plot and how does it correspond to the dots?
  
  In what scenario might you see this kind of data? 
  
  ```{r}
(m_bin <- glm(y ~ x, family = 
                binomial(link = "logit"), data = d)) # exercise
g <- ggeffects::ggpredict(m_bin, "x", full.data = TRUE) %>%
  ggplot(aes(x, predicted)) +
  geom_line(colour = "red") + 
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .2) +
  geom_point(data = d, aes(y = y))
g

