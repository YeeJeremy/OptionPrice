# *OptionPrice* - R Package 

## Author
Jeremy Yee

## Description

This package prices financial options using various numerical methods.

Work in progress.

## Example

### Bermuda Put Option

Pricing a Bermuda put option using the method by Longstaff and Schwartz 2001.

~~~
library(OptionPrice)
set.seed(123)
## Parameters
step <- 0.02
mu <- 0.06 * step
vol <- 0.2 * sqrt(step)
n_dec <- 51
start <- 36
strike <- 40
btype <- "power"  ## power, laguerre
## Generate paths
n_path <- 100000
path <- GBM(start, mu, vol, n_dec, n_path, TRUE)  ## discouted price process
## Regression basis
basis <- matrix(c(1, 1, 1), nrow = 1)
## Least squares Monte Carlo
lsm <- BermudaPutLSM(path, strike, exp(-mu), basis, TRUE, btype)
~~~
The computed option value is
~~~
> mean(lsm$value)
[1] 4.478689
> sd(lsm$value)/sqrt(n_path)
[1] 0.009158486
~~~
The computed bounds using the duality approach (Rogers)
~~~
## Out of sample paths and 1-step subsimulation
n_path2 <- 1000
path2 <- GBM(start, mu, vol, n_dec, n_path2, TRUE)
n_subsim <- 1000
subsim <- NestedGBM(path2, mu, vol, n_subsim, TRUE)
## Primal and Duals
bounds <- BermudaPutBounds(path2, subsim, lsm$expected, strike, exp(-mu), basis, btype)
~~~
is given by
~~~
> mean(bounds$primal)
[1] 4.406204
> sd(bounds$primal)/sqrt(n_path2)
[1] 0.07428008
> mean(bounds$dual)
[1] 4.810256
> sd(bounds$dual)/sqrt(n_path2)
[1] 0.01830421
~~~
