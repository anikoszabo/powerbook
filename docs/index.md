--- 
title: "Power and Sample Size Manual"
author: "Aniko Szabo"
date: "2023-06-28"
site: bookdown::bookdown_site
documentclass: book
bibliography: [book.bib, packages.bib]
url: https://anikoszabo.github.io/powerbook
# cover-image: path to the social sharing image like images/cover.jpg
description: |
  This is a reference book for power and sample size calculations using R
link-citations: yes
csl: ama.csl
github-repo: anikoszabo/powerbook
---

# Welcome {-}

This repository is for a [Bookdown](https://bookdown.org/) version of a reference book for power and sample size calculations using R.

It is very much a work in progress, and fixes/contributions in the form of pull requests are welcome.


## Guiding principles {-}

### Basic approach {-}

 * Each section describes the data assumptions, the statistical analysis methods (eg hypothesis test), and the sample size/power calculation formula or approach. 
 
 * R code or reference is provided for each analysis method and sample size calculation.
 
 * There is a worked example of a typical sample size calculation with this method.
 
 * The dependence of the power/sample size on the input parameters is visualized and discussed.
 
 * At least two simulation studies are included: one validating the analysis method (eg coverage for confidence intervals, type I error for tests), and one validating the sample size calculation.

### Coding style {-}

 * Use base R to the extent possible. The exceptions are:
 
    - packages implementing sample size / power calculations 
    - {ggplot2} for plotting
    
 * Use consistent naming conventions. Each section has a unique descriptive label (referred to as XX below) that is used in function names.
 
 * If multiple analysis approaches are possible, a wrapper function with options for the methods is provided. The function is called `XX_test` for hypothesis tests, and returns an object of class `htest`.
 
 * If multiple sample size approaches are possible, a wrapper function is provided. To the extent possible both calculating power for a given sample size, and the sample size for a given power are implemented. The function is called `power_XX_test` for hypothesis tests, and returns an object of class `power.htest`.
 
 * Simulation studies implement the following steps as separate functions to enable easy modification:
 
    - A function named `sim_XX_data` to generate multiple random data sets for a fixed set of parameters.
    - A function named `analyze_XX_sim` to analyze the simulated data that produces the summary statistic(s) of interest for each replicate. 
    - A function named `run_XX_sim` to loop through multiple settings for parameters that produces one dataset with all the results. 
    - The creation of simulation settings and the display of the simulation results can be in free code.

 * All new functions have Roxygen-style documentation
 
 * A simulation seed is set immediately before any code with randomness, so it can be reproduced as a stand-alone simulation


