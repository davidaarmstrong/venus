---
title: "Adaptive MARS"
author: "Duncan Murdoch and David Armstrong"
date: "2019-07-09"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Adaptive MARS}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

## Introduction
Despite recent advances in semi- and non-parametric techniques, most quantitative social scientists are still using regression models to substantiate their theoretical claims.  In these models, there are generally two types of explanatory variables - those that are theoretically important and the subject of hypothesis testing (let's call these $X$) and those that are controls (let's call these $Z$).  As any regression/research design text (e.g., [cite]) tells us, the control variables are meant to capture all alternative explanations of the dependent variable.  One question that almost always remains unanswered is whether we are using the appropriate functional forms for these control variables.  We propose a method that allows users to be confident that the control variables have the optimal relationship with both $y$ and $X$, without making functional form assumptions and without interfering with inference on the variables of interest. 

The importance of these relationships cannot be overstated.  Often control variables are considered to be "less important" than other more theoretically interesting variables and as such are often not treated with similar care in their specification.  We argue that there are several problems that can arise from this.  Most prominently are: 
 
  1. If the relationship between the control variables ($Z$) and $y$ is wrong, then there is extra variance in $y$ left to be explained by $X$.  This means that we are likely overstating the effects of $X$ on $y$.  This is the classic omitted variable bias finding.  
  2. Less well-known is a similar problem where the relationship between $X$ and $Z$ cannot be appropriately estimated with the design matrix for the regression of $y$ on $(X,Z)$.  This would be the case if, for instance $z_{1}$ had a linear relationship with $y$ and a quadratic relationship with $x_{1}$.  A conventional model would not remove enough of the variability in $x_{1}$, leaving more of its variance to explain variance in $y$, thus likely inducing some bias. 

The first problem above can be easily handled by careful model specification of $y$ on $(X,Z)$.  The solution to the second problem is a bit less clear.  

We propose a solution based on Multivariate Adaptive Regression Splines (MARS) [cite].  The general idea is that we can use a flexible model to identify the appropriate design matrix for the regression of $y$ on $Z$, call this $m(Z|y)$.  We can also generate the same result for the relationship between $X$ and $Z$, generating $m(Z|X)$.  We can combine these two design matrices $(m(Z|y),m(Z|X))$ and use those to residualize both $y$ and $X$, leaving $e^{(y)}$ and $e^{(X)}$.  By regressing the former on the latter, we can get clean estimates of the effect of $X$ on $y$ controlling for $Z$.  We show through simulation that this estimator is unbiased and consistent while producing confidence intervals with approximately correct coverage in a wide range of circumstances. 

## The Model

The model we propose is based on the Multivariate Adaptive Regression Splines (MARS&trade;).  This is a model that uses piecewise linear "hinge functions" to model arbitrary non-linearity.  

$$(x-t)_{+} &= \begin{cases}
  x-t & \text{ if } x > t\\
  0   &   \text{ otherwise.}
\end{cases}$$

$$(t-x)_{+} &= \begin{cases}
  t-x & \text{ if } x < t\\
  0   &   \text{ otherwise.}$$


Our example in `help(venus)` has this setup, but with
larger variances:

<!-- ```{r}
  set.seed(12345)
  library(MASS)
  library(venus)
  sigma2.e.x <- 4.7*3*0.01
  sigma2.e <- 7.3/2*0.01
  ## Z are the control variables or nuisance variables
  Z <- mvrnorm(1000, c(0,0,0),
             matrix(c(1,.3,.3,.3,1,.3,.3,.3,1), ncol=3),
             empirical=TRUE)
  ## X is the variable of interest that has a non-linear relationship to
  ## one of the elements in Z 
  X <- Z[,1] + Z[,1]^2 + Z[,2] + rnorm(1000, 0, sqrt(sigma2.e.x))

  df <- data.frame(z1 = Z[,1], z2=Z[,2], z3 = Z[,3], x=scale(X))
  ## make y a linear function of x and the z variables,
  ## so the theoretical coefficient on x is 1
  df$y <- with(df, z1 + z2 + z3 + x) + rnorm(1000, 0, sqrt(sigma2.e))
  ```
  
Mathematically this can be expressed as:
\[
\begin{align}
Z & \sim & MVN\left(\left(\begin{array}{c}0 \\ 0 \\ 0\end{array}\right), 
              \left(\begin{array}{ccc}1 & 0.3 & 0.3 \\
                     0.3 & 1 & 0.3 \\
                     0.3 & 0.3 & 1 \end{array} \right) \right)\\
X & \sim & N(Z_1 + Z_1^2 + Z_2, \sigma^2_{ex}) \\
x & = & \mbox{Scaled, centred X} \\
y & \sim & N(Z_1 + Z_2 + Z_3 + x, \sigma^2_e)
\end{align}
\]

Pairwise plots:


```r
  pairs(df)
```

```
## Error in as.data.frame.default(x): cannot coerce class '"function"' to a data.frame
```

We partial out the nuisance variables using `venus(y ~ x, nuisance = y ~ .-x, data = df)`.  These are the steps
of that calculation:

First, fit the nuisance model:

```r
library(earth)
```

```
## Loading required package: Formula
```

```
## Loading required package: plotmo
```

```
## Loading required package: plotrix
```

```
## Loading required package: TeachingDemos
```

```r
nuisanceFit <- earth(y ~ .-x, data = df)
```

```
## Error in as.data.frame.default(data, optional = TRUE): cannot coerce class '"function"' to a data.frame
```

```r
yResids <- residuals(nuisanceFit)
```

```
## Error in residuals(nuisanceFit): object 'nuisanceFit' not found
```

Now extract the predictors that were chosen, and get residuals
on X:

```r
nuisanceModelmatrix <- model.matrix(nuisanceFit)
```

```
## Error in model.matrix(nuisanceFit): object 'nuisanceFit' not found
```

```r
xResids <- residuals(lm(x ~ nuisanceModelmatrix, data = df))
```

```
## Error in as.data.frame.default(data, optional = TRUE): cannot coerce class '"function"' to a data.frame
```

and plot again:


```r
plot(yResids ~ xResids)
```

```
## Error in eval(predvars, data, env): object 'yResids' not found
```

Here's the interesting linear model:


```r
library(xtable); options(xtable.type = "html")
```

```
## 
## Attaching package: 'xtable'
```

```
## The following object is masked from 'package:TeachingDemos':
## 
##     digits
```

```r
xtable(lm(yResids ~ xResids - 1))
```

```
## Error in eval(predvars, data, env): object 'yResids' not found
```

# What happens if $x = z_1^2$?

If $x$ is an error free function of one of the
nuisance variables, we would expect that it adds nothing to
the regression.  However, if it is a nonlinear function, we
might not recognize that.  We try it with our function:


```r
  X <- Z[,1]^2
```

```
## Error in eval(expr, envir, enclos): object 'Z' not found
```

```r
  df <- data.frame(z1 = Z[,1], z2=Z[,2], z3 = Z[,3], x=scale(X))
```

```
## Error in data.frame(z1 = Z[, 1], z2 = Z[, 2], z3 = Z[, 3], x = scale(X)): object 'Z' not found
```

```r
  ## make y a linear function of x and the z variables,
  ## so the theoretical coefficient on x is 1, but
  ## z1^2 gives the same information
  df$y <- with(df, z1 + z2 + z3 + x) + rnorm(1000, 0, sqrt(sigma2.e))
```

```
## Error in eval(substitute(expr), data, enclos = parent.frame()): invalid 'envir' argument of type 'closure'
```
 
  Pairwise plots:


```r
  pairs(df)
```

```
## Error in as.data.frame.default(x): cannot coerce class '"function"' to a data.frame
```

The `venus()` fit:


```r
  summary(venus(y ~ x, nuisance = y ~ .-x, data = df))

  # Compare to the pure linear fit

  summary(venus(y ~ x, nuisance = y ~ .-x, method = c("lm", "lm"), data=df))
``` -->
```

```
## Error: attempt to use zero-length variable name
```
