
---
title: "Introduction to Monte Carlo methods for estimation"
author: "Asael Alonzo Matamoros"
date: "2020-04-17"
layout: post
image: https://img.zhaohuabing.com/post-bg-2015.jpg
tags: ["R","Monte Carlo", "programming", "code"]
---

Introduction
============

Two major classes of numerical problems that arise in data analysis
procedures are optimization and integration problems. It is not always
possible to analytically compute the estimators associated with a given
model and we are often led to consider numerical solutions. One way to
afford that problematic is to use simulation. Monte Carlo estimation
refers to simulating hypothetical draws from a probability distribution,
in order to calculate important quantities of that distribution.

The basic idea of Monte Carlo consist of writing the integral as an
expected value with respect to some probability distribution, and then
approximated using the [method of moment
estimator](https://en.wikipedia.org/wiki/Method_of_moments_statistics.)  
$$E\[g(X)\] \\approx \\overline{g(X)} = \\dfrac{1}{n}\\sum g(X\_{i})$$

If we have a continuous function *g*(*θ*) and we want to integrated in
the interval (a,b), we can rewrite our integral as an expected value of
an uniform distribution *U* ∼ *U*\[*a*, *b*\], that is:

$$I = \\int\_{a}^{b}g(\\theta)d\\theta = \\int\_{a}^{b}\[g(\\theta)(b-a)\]\\dfrac{1}{(b-a)}d\\theta = E\_{U}\[(b-a)g(\\theta)\]$$
Using the method of moments estimator our integral approximation is:

$$I =\\int\_{a}^{b}g(\\theta)d\\theta \\approx \\dfrac{1}{n}\\sum\_{k=1}^{n}(b-a)g(\\theta\_{k})$$

Where the
*θ*<sub>1</sub>, *θ*<sub>2</sub>, *θ*<sub>3</sub>, *θ*<sub>4</sub>, …, *θ*<sub>*n*</sub>
are simulated values from an uniform distribution.

Example 1: Exponential integral approximation
=============================================

1.  Given a function *f*(*x*) = *e*<sup>*x*</sup> the integral in the
    interval \[3,5\] is:

*I* = ∫<sub>3</sub><sup>5</sup>*f*(*x*)*d**x* = ∫<sub>3</sub><sup>5</sup>*e*<sup>*x*</sup>*d**x* = *e*<sup>*x*</sup>\|<sub>3</sub><sup>5</sup> = *e*<sup>5</sup> − *e*<sup>3</sup> = 128.3276

The MonteCarlo approximation of the integral is:

``` r
#       Declaring the desired function
f = function(x){return(exp(x))}
#       Declaring the absolute error function
error = function(x,y){return(abs(x-y))} 
#       The actual integral answer
ans = exp(5)-exp(3) 

set.seed(6971)
#       number of iterations
n = 10^2 
#       simulated uniform data
x= runif(n,3,5) 
#       MonteCarlo approximation
MCa= (5-3)*mean(f(x)) 
#       Approximation error
e = error(ans,MCa)

rest =  data.frame(n = n,MCapprox = MCa,error = e)
set.seed(6971)
for(k in 3:6){
  n = 10^k
  x = runif(n,3,5)
  mca = (5-3)*mean(f(x))
  rest = rbind(rest,c(n,mca,error(ans,mca) ) )
}

kable(rest,digits = 5,align = 'c',caption = "Integral Monte Carlo approximation results",
      col.names =c("Number of simulations","Monte Carlo approximation","Error approximation"))
```

| Number of simulations | Monte Carlo approximation | Error approximation |
|:---------------------:|:-------------------------:|:-------------------:|
|         1e+02         |          127.2634         |       1.06424       |
|         1e+03         |          127.1676         |       1.16001       |
|         1e+04         |          128.9637         |       0.63608       |
|         1e+05         |          128.2422         |       0.08546       |
|         1e+06         |          128.3221         |       0.00552       |

Generalized Monte Carlo approximation
=====================================

In a general case, the integral approximation for a given distribution f
is:

$$I= \\int\_{a}^{b}g(\\theta)d\\theta =\\int\_{a}^{b}\\dfrac{g(\\theta)}{f(\\theta)} f(\\theta)d\\theta =E\_{f}\[g(\\theta)/f(\\theta)\] \\approx \\dfrac{1}{n}\\sum\_{k=1}^{n}\\dfrac{g(\\theta\_k)}{f(\\theta\_k)}$$

An algorithm for construction of *Î* can be described by the following
steps:

1.  Generate
    *θ*<sub>1</sub>, *θ*<sub>2</sub>, *θ*<sub>3</sub>, …, *θ*<sub>*n*</sub>
    from a f distribution

2.  Calculate:
    *g*(*θ*<sub>1</sub>)/*f*(*θ*<sub>1</sub>), *g*(*θ*<sub>2</sub>)/*f*(*θ*<sub>2</sub>), *g*(*θ*<sub>3</sub>)/*f*(*θ*<sub>3</sub>), …, *g*(*θ*<sub>*n*</sub>)/*f*(*θ*<sub>*n*</sub>)

3.  Obtain the sample mean:
    $$\\overline{I} = \\dfrac{1}{n}\\sum\_{k=1}^{n}\\dfrac{g(\\theta\_k)}{f(\\theta\_k)}$$

In the next chunk, the simple Monte Carlo approximation function is
presented to show how does the algorithm works, where a and b are the
uniform density parameters, n the number of desired simulations, and f
is the function the we want to integrate.

``` r
# The simple Monte Carlo function
MCaf = function(n,a,b,f){
x = runif(n,a,b)
MCa = (b-a)*mean(f(x))
return(MCa)
}
```

Monte Carlo methods in Bayesian data analysis
=============================================

The main idea of the bayesian data analysis is fitting a model (*such as
a regression or a time series model*) using a [bayesian
inference](https://en.wikipedia.org/wiki/Bayesian_inference) approach.
We assume that our parameters of interest have a theoretical
distribution, this distribution (*posterior*) is updated using the
distribution of the observed data (*likelihood*), and the previous or
external information about our parameters (*prior distribution*) by
using the [Bayes’
theorem](https://en.wikipedia.org/wiki/Bayes%27_theorem).

*P*(*θ*/*X*) *α* *P*(*X*/*θ*)*P*(*θ*)
Where:

*P*(*θ*/*X*) is the parameter [posterior
distribution](https://en.wikipedia.org/wiki/Posterior_probability)

*P*(*X*/*θ*) is the sampling distribution of the observed data (
[likelihood](https://en.wikipedia.org/wiki/Likelihood_function) ).

*P*(*θ*) is the parameter [prior
distribution](https://en.wikipedia.org/wiki/Prior_probability).

The main problematic in the bayesian approach is estimating the
posterior distribution. The markov chain monte carlo methods
([mcmc](https://en.wikipedia.org/wiki/Markov_chain_Monte_Carlo))
generates a sample of the posterior distribution, and approximate the
expected values, probabilities or quantiles using Monte Carlo methods.

In the next two sections, we provide two examples for approximating
probabilities and quantiles of a theoretical distribution. The procedure
presented above are the usual methodologies used in a bayesian approach.

Example 2: Probability approximation of a gamma distribution
------------------------------------------------------------

Lets suppose we want to calculate the probability that a random variable
*θ* is between zero and 5 *P*(0 \< *θ* \< 5), where *θ* has gamma
distribution with parameter a = 2 and b = 1/3
(*θ* ∼ *G**a**m**m**a*\[*a* = 2, *b* = 1/3\]), so the probability is:

$$P(0 \\leq \\theta \\leq 5) = E\[I\_{\[0,5\]}(\\theta)\] = \\int\_{0}^{\\infty} I\_{\[0,5\]}(\\theta) \\dfrac{b^{a}}{\\Gamma(a)}\\theta^{a -1}e^{-\\theta/ b}d\\theta  = \\int\_{0}^{5} \\dfrac{b^{a}}{\\Gamma(a)}\\theta^{a -1}e^{-\\theta/ b}d\\theta \\approx \\dfrac{1}{n}\\sum\_{k = 1}^{n} I\_{\[0,5\]}(\\theta\_{k})$$
Where *I*<sub>\[0, 5\]</sub>(*θ*) = 1 if *θ* belongs to the interval
\[0,5\].The idea of the Monte Carlo approximation, is count the number
of observations that belong to the interval \[0,5\], and divide it by
the total of the simulated data.

``` r
set.seed(6972)
#       number of iterations
n = 10^2 
#       simulated uniform data
x= rgamma(n,shape = 2,1/3) 
#       MonteCarlo approximation
MCa= mean(x <= 5) 
#       Approximation error
e = error(pgamma(5,2,1/3),MCa)

rest =  data.frame(n = n,MCapprox = MCa,error = e)

for(k in 3:6){
  n = 10^k
  x= rgamma(n,shape = 2,1/3) 
  mca= mean(x <= 5) 
 rest = rbind(rest,c(n,mca,error(pgamma(5,2,1/3),mca) ) )
}

kable(rest,digits = 5,align = 'c',caption = "Probability Monte Carlo approximation results",
      col.names =c("Number of simulations","Monte Carlo approximation","Error approximation"))
```

| Number of simulations | Monte Carlo approximation | Error approximation |
|:---------------------:|:-------------------------:|:-------------------:|
|         1e+02         |          0.51000          |       0.01367       |
|         1e+03         |          0.49300          |       0.00333       |
|         1e+04         |          0.48930          |       0.00703       |
|         1e+05         |          0.49753          |       0.00120       |
|         1e+06         |          0.49663          |       0.00030       |

Example 3: Quantile approximation of a normal distribution
----------------------------------------------------------

Lets suppose we want to calculate the 0.95 quantile of a random variable
*θ* that has normal distribution with parameters *μ* = 20 and *σ* = 3
(*θ* ∼ *n**o**r**m**a**l*(*μ* = 20, *σ*<sup>2</sup> = 9)), so the 0.95
quantile is:

$$\\int\_{-\\infty}^{q\_{95}}\\dfrac{1}{\\sqrt{2\\pi\\sigma^{2}}}  e^{\\dfrac{(\\theta - \\mu)^{2}}{2\\sigma^{2}} }d\\theta = 0.95$$
The main idea is found the largest sample value that gives a probabilty
equal or less than 0.95, the Monte Carlo quantile approximation is
estimate it using the *quantile()* function of the simulated data.

``` r
set.seed(6973)
#       number of iterations
n = 10^2 
#       simulated uniform data
  x= rnorm(n,20,3) 
#       MonteCarlo approximation
MCa= quantile(x,0.95) 
#       Approximation error
e = error(qnorm(0.95,20,3),MCa)

rest =  data.frame(n = n,MCapprox = MCa,error = e)

for(k in 3:6){
  n = 10^k
  x= rnorm(n,20,3)
  mca= quantile(x,0.95)
 rest = rbind(rest,c(n,mca,error(qnorm(0.95,20,3),mca) ) )
}

kable(rest,digits = 5,align = 'c',caption = "Quantile Monte Carlo approximation results",
      col.names =c("Number of simulations","Monte Carlo approximation","Error approximation"),row.names = FALSE)
```

| Number of simulations | Monte Carlo approximation | Error approximation |
|:---------------------:|:-------------------------:|:-------------------:|
|         1e+02         |          25.76068         |       0.82612       |
|         1e+03         |          24.83333         |       0.10123       |
|         1e+04         |          24.98454         |       0.04997       |
|         1e+05         |          24.94259         |       0.00803       |
|         1e+06         |          24.93760         |       0.00304       |

Discussions and conclusions
===========================

The Monte Carlo approximation methods offer an alternative tool for
integral approximation, and is a really important tool in the bayesian
inference approach, specially when we work with sophisticated and
complex models. As it seems in all our three examples, the Monte Carlo
methods offer a really good approximation, but it demands a huge number
of simulations for getting an approximation error close to zero.

References
==========

1.  *Introducing Monte Carlo methods with R*, Springer 2004,
    *Christian P. Robert and George Casella*.

2.  *Handbook of Markov Chain Monte Carlo*, Chapman and Hall, *Steve
    Brooks, Andrew Gelman, Galin L. Jones and Xiao-Li Meng*.

3.  *Introduction to mathematical Statistics*, Pearson, *Robert V. Hogg,
    Joseph W. Mckean, adn Allen T. Craig*.

4.  *Statistical Inference An Integred Approach*, Chapman and Hall,
    *Helio S. Migon, Dani Gamerman, Francisco Louzada*.
