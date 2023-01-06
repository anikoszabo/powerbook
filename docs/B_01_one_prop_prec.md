# One-sample proportion: precision-based sample size

## Example

An educator wants to conduct a survey of first-year medical residents. The survey contains yes/no questions for several specific medical procedures asking whether they would be comfortable performing it without the presence of an attending physician. How many residents should he recruit?

## Data generating mechanism

Independent identically distributed variables $X_1, \ldots, X_n \sim Bernoulli(p)$ with $p$, the event probability as the parameter of interest.

It will be estimated as the sample proportion $\hat{p} = \sum_{i=1}^n X_i / n$. Inference can be based on the exact distribution of $X = \sum_{i=1}^n X_i \sim Binom(n,p)$ or the asymptotic distribution of $\hat{p} \sim N(p, \frac{p(1-p)}{n})$.

## Precision measures

### Standard error

The standard error of $\hat{p}$ is

\begin{equation}
SE(\hat{p}) = \sqrt\frac{p(1-p)}{n} 
\end{equation}



```r
#' Calculate standard error for given sample size, or vice versa. Exactly one of 'n' and 'se' has to be NULL
#' p - expected event probability
#' n - sample size
#' se - standard error of estimated probability
prec.1prop <- function(p, n=NULL, se=NULL){
  if (is.null(n) & !is.null(se)){
    n <- p * (1-p) / se^2
  } else if (!is.null(n) & is.null(se)){
    se <- sqrt(p*(1-p)/n)
  } else stop("Exactly one of 'n' and 'se' has to be NULL")
  
  res <- list(p = p, n = n, se = se)
  res
}
```

### Confidence interval or margin of error

## Getting inputs, worst/best case scenarios

## Example revisited 

## Simulation study

## Functions in R packages