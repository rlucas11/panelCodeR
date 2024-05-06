# Package Development Information

This page describes the structure of the package for future reference when adding features. The package is broken down into the following files, which include the following functions.

# Approach

The `panelcoder` command goes through the following steps to build and test a model:

* Get information about the data that will be analyzed
  - Determine whether to use a univariate or bivariate model
  - Identify the number of waves
  - Extract variable names
  - Identify phantom variables
  
* Determine which parts of the model to include  
* Build an initial parameter table (using functions in buildParTable.R)
* Add constraints based on model specifications
* Optionally add limits to help with convergence issues
* If lavaan: Run code from parameter table
* If Mplus: Build Mplus code from lavaan parameter table and run (with MplusAutomation)
* Create summary information for output
  - Key information in summary
  - Model code used to run model
  - Implied and observed correlations (only available for single-indicator models)
  - Full model output from lavaan or Mplus



# Files

## buildParTable.R

### .buildObserved

This function creates the code for the observed variables and the latent occasions. Specifically, for each variable at each observed occasion (phantom variables are created in `.buildPhantom`) it creates a new latent variable (labeled "l" plus the variable name and wave, separated with "_"). 

For measures with only one indicator, the observed variable is linked with this latent variable through a loading that is constrained to be 1. The variance for the observed variable is also set to 0, as the code creates a separate residual term. 

For measures with more than one indicator, each indicator loads on the latent trait. The loading of the first indicator is fixed to 1; each of the other loadings are labeled with "x" or "y", the wave number, and a letter representing the indicator (e.g., "b", or "c"). 

This function also sets the residual variance for the "l" variables at each wave to zero (as any residual variance can be captured in the "state" component if specified). 

### .buildPhantom

This function creates the code for waves where the observed variables are missing. This can only be used in cases where stationarity is imposed. 

For all missing waves, the code uses lavaan syntax to predict the "l" variable from itself, with a loading of 0 (this is how phantom variables are specified in parameter tables). 

### .buildResidVar

This function creates the residual terms when there is more than one indicator. 

For each indicator at each wave, a variance component is specified and labeled "x" or "y" plus the wave, the indicator, and a "v" (e.g., "x_1_1v" for the first indicator for the X variable at Wave 1). 

### .buildAr

This function creates the autoregressive components, along with the necessary loadings and variances. Eventually, I may constrain these variances to be 0 and create explicit "impulse" components, as per Zyphur et al. For now, however, this function sets the loadings and variances at the same time. 

#### Latent Variables and Loadings

First, the code specifies a new "a" variable for each variable at each wave. This is simply labeled with "a" and the wave, separated by "_". The "l" variables representing the latent occasion load on the wave-specific "a" variable with a loading that is constrained to 1. 

#### Variances

Next, the code creates variance terms for each "a" variable, which are labeled "xvar" with the wave (e.g., "xvar1" or "xvar2"). The labels are created, but constraints on these variables are added elsewhere. 

### .buildTrait

This function creates the "random intercept" or "stable trait" component. Note that when using the Dynamic Panel Model, the nature of this component differs, so a different function is used for the DPM. 

The function creates a latent trait for each variable in the model (up to two). Each latent occasion is specified to load on this latent trait with a loading of 1. In future versions, this constraint will be relaxed to allow for time-varying effects of the stable trait. 

The function also creates the variance term for this latent variable, setting the labels to be "x_tVar" (and "y_tVar" if it is a bivariate model). 

### .buildDpmTrait

This function builds the corresponding stable-trait term for the DPM. Instead of having the latent occasion variables load on this latent trait, each autoregressive variable loads on the latent trait. As above, the loadings are constraint to 1, but this can be relaxed in future versions. 

This function also creates the variance term, which is labeled in the same way as for the stable trait/random intercept component (e.g., "x_tVar" and "y_tVar"). 

### .buildStability

This function builds the stability paths for the autoregressive part of the model. Each "a" variable (beyond the first) is linked with the corresponding "a" variable from the prior wave (e.g., "a_x_2" is predicted from "a_x_1"). Each path has a loading that is labeled with either an "a" (for X variables) or a "b" (for Y variables) plus the wave of the outcome variable (e.g., "a2" when predicting "a_x_2" from "a_x_1"). 

No constraints on these stability paths are added at this point. 

### .buildCrossLag

This function builds the lag-1 cross-lagged paths for the autoregressive part of the model. Each "a" variable (beyond the first) is linked with the corresponding "a" variable from the other construct at the prior wave (e.g., "a_x_2" is predicted from "a_y_1"). Each path has a loading that is labeled with either a "c" (for Y predicted from X) or "d" (for X predicted from Y) plus the wave of the outcome variable (e.g., "c2" when predicting "a_y_2" from "a_x_1"). 

### .buildState

This function builds the state variance components for each variable at each wave. This can account for measurement error in single-indicator models and reliable occasion-specific variance in multiple-indicator models. 

The code creates state variables (named using "s", the variable name, and the wave number, separated by "_") and specifies that the latent occasion variables load on this new latent state variable with a loading constrained to be 1. 

The variances for these state variances are estimated (unless later constrained to be 0) and labeled with either "sx" or "sy" plus the wave number. These can be constrained in various ways later. 

### .buildCors

This function creates correlations across constructs for the various components in the models. These correlations are specified differently for the DPM, so there is a separate function for that. 

The following correlations are estimated:

* Correlations between stable traits/random intercepts (labeled "cov_txty").
* Correlations between wave-specific autoregressive disturbances (labeled with "cov_ar" and wave). This may change if I add explicit "impulse" terms. 
* Correlations between wave-specific state components (labeled with "cov_s" and wave). 

### .buildDpmCors

This function creates correlations across components when the DPM is selected. It is similar to the previous function, but in addition to specifying correlations between the two stable-trait components, it also allows these to correlate with the AR variables from the first waves. Both stable-trait variables correlate with each other and with both initial AR variables. 


### .buildResidCors

This function creates correlations between indicator-specific residuals across waves. Each correlation is labeled with "x" or "y", the higher wave number, and the length of the lag (e.g., "x2a_l1" for the correlation between Indicator A for the X variable across waves 1 and 2). These are not constrained in any way; equality constraints can be added later. 

### .buildTable

This function puts everything together to build the initial parameter table that lavaan can use to run a model (and that can be turned in to an mplus model). This function creates a list of all parts of the parameter table using the functions described above. Then, it selectively sets certain entries to be NULL, depending on the specific model being tested. The elements of this list are then combined using `rbind` to create the final parameter table. 

This function also identifies latent and observed variables. This is necessary to correctly specify intercepts and means when using FIML to deal with missing data. The function returns a list of tables, including the parameter table that is used to run the model, along with other tables that include subsets of these parameters (which are used elsewhere to specify other features of the model).

## buildConstraints

This file includes a set of functions used to build constraints that are either necessary to run the models or that can be changed depending on options specified by the user. 

### .buildConstraint

This is a generic function that takes an existing parameter table and adds a line reflecting the new constraint. The input to the function includes what is on the left-hand side of the equation, the operator, and what is on the right-hand side of the equation. 

### .constrainStability

This function constrains the stability coefficients to be equal. This is usually implemented as a part of a broader set of stationarity constraints. 

### .constrainCl

This function constrains the cross-lag paths from the same variable to be equal across waves. Again, this is currently only implemented when stationarity is specified to be TRUE. 

### .arStationarity

This function constrains the autoregressive variances to be equal, using appropriate non-linear constraints. The variances after the first wave are a function of the variance at the first wave, the stability coefficients, and the cross-lagged paths. 

This function can also be used to constrain these paths to be zero if cross-lagged paths are not wanted. 

### .constrainStateVar

This function constrains the state variances for each variable to be equal across waves. This is included when stationarity is set to TRUE. 

### .constrainResidVar

This function constrains the residual variance for a specific indicator to be equal across waves. 

### .constrainLoadings

This function constrains the loadings for the indicators to be equal across waves. 

### .constrainStateCors

This function constrains the correlations between state components to be the same across waves. This can also be used to constrain these to be zero. 

### .buildLimits

This function constrains certain variances to be greater than 0 and certain correlations to be less than 1. The following parameters are constrained:

* Autoregressive variance
* Trait variance
* Correlation between traits
* State variance
* Check: "x1_lv"

I will likely add additional parameters.
