# panelCodeR Package

This R package generates lavaan and mplus code for models for analyzing panel data. Currently, it can generate code for the bivariate STARTS model and a number of models nested under this more general model. The diagram below shows the STARTS model and its components, which include a stable-trait (ST) component, an autoregressive trait (ART) component, and a state (S) component for each variable. 

![Diagram of STARTS Model](images/startsTransparent.png)

Reduced form models can be specified by setting certain variance components equal to zero. For instance, omitting the state variance at each wave from the general model results in the RI-CLPM, and removing both the state and stable-trait variance results in the CLPM. 

![Diagram of RI-CLPM](images/riclpm.png =45%x)![Diagram of CLPM](clpmReduced.png =75%x)

Currently, it is possible to run bivariate and univariate versions of the following models: 

- STARTS (includes all components)
- RI-CLPM (excludes state component)
- ARTs (also known as a 'factor CLPM'; excludes stable-trait component)
- STS (excludes autoregressive component)
- CLPM (excludes stable-trait and state components)

A number of additional options are also available. By default, the model imposes stationarity constraints, which set the variances, stabilities, and cross-lagged paths to be equal across waves; but this constraint can often be removed. It is also possible to remove the lagged paths. Correlations between the components for the two variables in bivariate models can also be included or excluded. 

The package handles data with multiple indicators per wave (as long as the same indicators are available at each wave). In addition, most models can be run even if certain waves are missing. This is accomplished through the use of phantom variables for missing waves, but stationarity constraints must be included for these to work. 

Finally, there are a few helper functions. The package includes a function called `gen_starts()` that generates data from a STARTS model with user specified parameters (e.g., different amounts of variance for each component or different lagged paths). This can be useful for seeing how the models behave under different data-generating processes. There is a related function called `addIndicators()` that can take the data from `gen_starts()` and make multiple indicators for each wave. Finally, the function `parcels()` takes multiple-item scales and creates parcels based on the average loadings across waves. 

There are three functions that are useful after you have run a model. `panelPlot()` will plot the implied stabilities and actual stabilities for increasingly long waves (up to the length of the study). This can be useful for visualizing the ways that the model might not describe the underlying data well. The function `panelEstimates()` prints all estimates from the model. The function `modelCode()` prints a formatted version of lavaan or mplus code used to run the model. This can be useful for checking that everything was specified correctly or for modifying the basic code for variations that the package cannot handle. 


## Installation

You can install this package from from [GitHub](https://github.com/rlucas11/panelCodeR) with:
      
``` r
# install.packages("devtools")
devtools::install_github("rlucas11/panelCodeR")
```

## Data

One of the benefits of the package is that it flexibly creates model code based on the data you provide. In other words, the function determines how many variables (either one or two), waves, and indicators exist in your data, and then creates a model appropriate for those data. This means, however, that you must name your variables in a specific way for the code to work. In short, your variables must at least have a stem and a wave number, separated by an underscore (e.g. X_1, X_2, X_3). If you have multiple indicators, then you should add a number that indexes the indicator, again separated by an underscore (e.g., X_1_1, X_1_2, X_2_1, X_2_2). If there is a wave missing, just do not include those variables in the dataframe that is passed to the function (e.g., do not include X_3 if the variable was not assessed at Wave 3). The dataframe should not have any variables other than those to be used in the model (i.e., no ID variables).

**Note for Mplus**: Mplus restricts variable names to eight characters. Because of the requirement for wave and indicator indexes (along with separators), this means that the stem of your variable names must be very short. For instance, if you have more than nine waves, five of the eight characters are taken up by the wave and indicator indexes, which means that the stem cannot be more than 3 characters long. Keep this in mind when naming your variables. 


## Usage

The main function to build and run the code is `panelcoder()`. 

```R
panelcoder <- function(data,
                       title = "panelcoder",
                       panelModel = "starts",
                       program = "lavaan",
                       crossLag = TRUE,
                       traitCors = TRUE,
                       arCors = TRUE,
                       stateCors = FALSE,
                       residCors = FALSE,
                       residVar = FALSE,
                       limits = TRUE,
                       stationarity = TRUE,
                       invariance = TRUE,
                       mplusAnalysis = NULL,
                       mplusOutput = NULL,
                       mplusDirectory = "mplus",
                       constrainCors = TRUE,
                       run = TRUE,
                       ...
                       )
```

Hopefully, these options are self-explanatory. You will need to specify the (properly named) dataframe to use. Without specifying anything else, that will run a STARTS model using lavaan, with all the default settings. To change which model is run, use the `panelModel` option to specify whether to run the "starts", "riclpm", "clpm", "arts", or "sts" model. Eventually, there will be an option to run a dynamic panel model ("dpm"), but for now this requires using a different command (see below). To change which program to use, select either "lavaan" or "mplus" for the `program` option. 

The other options are described below (and in the R help functions).

- `crossLag` specifies whether the reciprocal lagged associations between the two variables in a bivariate model are included.
- `traitCors` specifies whether to include the correlation between the stable trait components.
- `arCors` specifies whether to include correlations between AR components in the same wave (in bivariate models).
- `stateCors` specifies whether to include correlations between state components in the same wave (in bivariate models).
- `residCors` specifies whether to include correlations between the residuals from a specific indicator at a specific wave with the residuals for the same indicator at other waves. 
- `residVar` specifies whether to constraint the residuals for specific indicators to be the same across waves
- `limits` specifies whether to restrict variances to be greater than zero and correlations to be between -1 and 1. This is sometimes needed to prevent inadmissible solutions when a variance component is very small or zero. 
- `stationarity` specifies whether to impose stationarity. If true, the following parameters are constrained
  - Stability coefficients
  - Lagged paths
  - Total autoregressive variance
  - State variance
- `mplusAnalysis` specifies the ANALYSIS command to use for mplus if different from the default.
- `mplusOutput` specifies the OUTPUT command to use for mplus if different from the default.
- `mplusDirectory` specifies the directory to store .inp and .out files from mplus. This needs to be created before running the model or the command will fail with a 'file not found' error. 
- `constrainCors` specifies whether to constrain correlations between residuals with equal lags to be equal (e.g., the correlation between X_1_1 and X_2_1 would be constrained to be equal to the correlation between X_2_1 and X_3_1. 
- `run` specifies whether to run the model or just to create and print the model code. 

When you call the function, it will present a summary of some of the most important results from the model, but it is best to save the output of the command to an object (e.g., `pcOutput <- panelcoder(data)`. This object includes some information used when constructing the model, along with the lavaan or mplus code, the lavaan or mplus output, and information to plot correlations. The following functions help examine this information (these are not yet implemented, but will be soon).

The function `panelplot()` will plot the implied and actual stabilities for increasingly long lags (up to the number of waves in the data). This can show whether the selected model accurately reproduces the actual stability coefficients from the data. Models like the CLPM often underestimate the long-term stability of the variables. 

The function `modelstatement()` will produce a formatted version of the mplus or lavaan code representing the model. This could be used to check that the model was specified correctly, and it can also be copied, pasted, and modified to run models that can't be specified with panelCodeR. 

The function `modelestimates()` will extract and print the lavaan or mplus estimates. The actual lavaan or mplus object that is created from running the model is stored in the panelCodeR output; you can access it by selecting the fourth element of the panelCodeR output (e.g., `pcOutput[[4]]`). 



# Old Version

Below is the documentation for the functions used in the early versions of this package. These work for now, but will eventually be deprecated. 

For both Lavaan and Mplus, there are functions just to build the model code or to build and run the code. It is often easiest to do the latter, but the former functions are especially useful if you want to build and then modify code. 


## Lavaan Commands

The basic (and most flexible) command to build the model is `buildLavaan()`. It has the following options and defaults:

```r
buildLavaan(waves,              # Number of total waves (e.g., 10)
           XVar = TRUE,         # Include X variables
           YVar = TRUE,         # Include Y variables
           xWaves = NULL,       # The actual waves for X (leave blank if no missing waves)
           yWaves = NULL,       # The actual waves for Y (leave blank if no missing waves)
           xIndicators = 1,     # How many indicators for X?
           yIndicators = 1,     # How many indicators for Y?
           stationarity = TRUE, # Stability, cross-lagged paths, and variances constrainted across waves)
           trait = TRUE,        # Include trait component
           AR = TRUE,           # Include AR component
           state = FALSE,       # Include state component
           crossLag = TRUE,     # Include cross-lagged paths
           stateCor = FALSE,    # Include correlation between wave-specific state components
           constrainCors = TRUE,# Constrain cross-wave indicators-specific cors to be equal
           limits = TRUE        # Limits variances > 0 and correlations < 1  
           )
```

Hopefully these options are pretty self-explanatory. If you only want a univariate model and your variables are named x1 through xw, you would set `YVar = FALSE`. If you wanted to exclude lagged directional paths, you could set `crossLag = FALSE` and maybe set `stateCor = TRUE` to allow for within wave correlations. If you set `AR = FALSE` you can specify a latent-trait model with no autocorrelations across waves. 

There are also some wrapper functions that pre-specify these options to create code for different models. These are: 
- `lavaanStartsX()` and `lavaanStartsY()`, which create code for univariate STARTS models (for variables named 'x' and 'y' respectively, in case you want to test univariate models for your two variables before testing bivariate models)
- `lavaanStarts2`, which creates code for bivariate STARTS models.
- `lavaanRiclpm`, which creates code for the RI-CLPM (i.e., it sets the state variance from the STARTS to 0).
- `lavaanClpm`, which creates code for the CLPM. 
- `lavaanArts`, which creates code for the bivariate ARTS model (i.e., CLPM with state component). 

For each of these, the default is to set `stationarity = TRUE` which constrains the parameter estimates for the stabilities, cross-lagged paths, variances, and covariances to be equal across waves. This can be changed by setting `stationarity = FALSE`. For some models, this may lead to identification problems. If you have missing waves, you will need to assume (and impose) stationarity.

### Running the models

The code generates the lavaan model code, but you still have to run the model. There are two options: Using the wrapper functions to run the code or generating the code using this script, saving it to a model object, checking the code, and then running it. Because the lavaan `sem()` command has some defaults that don't work with these models (e.g., allowing all latent variables to correlate), you should usually use the `lavaan()` command instead (or specify options when running the command). So you might use the following code:

```r
startsModel <- lavaanStarts2(10) ## Generate the lavaan model
cat(startsModel) ## Check the model code
startsFit <- lavaan(startsModel, data) ## Run lavaan on the model
summary(startsFit) ## Check your results!
```

#### Commands to run the models

The basic command to simultaneously build and run Lavaan code is `run_starts_lavaan()`. It has the same options as `buildLavaan()`, but with two additional options: One to specify the dataframe to use and one to pass any additional options on to the `lavaan()` command (such as different estimators or specifying how to handle missing data). So the syntax for `run_starts_lavaan()` is:

```r
run_starts_lavaan <- function(data,
                              waves,
                              XVar = TRUE,
                              YVar = TRUE,
                              xWaves = NULL,
                              yWaves = NULL,
                              xIndicators = 1,
                              yIndicators = 1,
                              trait = TRUE,
                              AR = TRUE,
                              state = TRUE,
                              crossLag = TRUE,
                              stateCor = FALSE,
                              stationarity = TRUE,
                              constrainCors = TRUE,
                              limits = TRUE,
                              ...)

```

The wrapper functions available are:

- `run_startsX_lavaan()`
- `run_startsY_lavaan()`
- `run_riclpm_lavaan()`
- `run_clpm_lavaan()`
- `run_arts_lavaan()`

I may add additional wrapper functions in the future.

### Missing Data

If you have missing data and want to use FIML, then you should add `missing="FIML", int.ov.free=TRUE, int.lv.free=FALSE` to the end of the command. For example, to test a STARTS model, you could use:

```R
startsFit <- run_starts_lavaan(data, 10, missing="FIML", int.ov.free=TRUE, int.lv.free=FALSE)
```

## Mplus commands

Because you can't call mplus directly from R, the structure of the mplus codes and the commands to call it differ somewhat from how this works for lavaan. There are two ways you can use the code with mplus. 

### run_starts_mplus()

I think the easiest requires that you have the `MplusAutomation` package installed. If you have this package, then the function `run_starts_mplus()` can build the model statement and create an mplus .inp file that is then run using mplus. It is best to direct this command to an object, as it will return all the info from the output file (e.g., parameter estimates, fit indices, etc.), which you can examine later (though the full output is still saved in a standard mplus .out file). You need to specify where the .inp and .out files should be written (the default is a directory called "mplus" in the directory in which you're currently working; if it doesn't exist, you should create it or specify an alternative location). The specific options are similar to those from the `buildLavaan()` command from above. Here are the options and defaults:

```R
run_starts_mplus <- function(data,
                             waves,               # Total number of waves (e.g., 10)
                             XVar = TRUE,         # Include X variable
                             YVar = TRUE,         # Include Y variable
                             xWaves = NULL,       # Which X waves exist (e.g., c(1:4, 6:10))
                             yWaves = NULL,       # Which Y waves exist
                             xIndicators = 1,     # How many indicators for X?
                             yIndicators = 1,     # How many indicators for Y?
                             trait = TRUE,        # Include Stable Trait
                             AR = TRUE,           # Include Autoregressive Trait
                             state = TRUE,        # Include State
                             crossLag = TRUE,     # Include cross-lagged paths
                             stateCor = FALSE,    # Include correlations between wave-specific states
                             stationarity = TRUE, # Impose stationarity
                             constrainCors = TRUE,# Constrain cross-wave indicators-specific cors to be equal
                             limits = TRUE,       # Limit variances and correlations to plausible values
                             dir="mplus",
                             title="test",
                             analysis = "MODEL=NOCOVARIANCES;\nCOVERAGE=.001;",
                             output = "stdyx; \n  cinterval; \n")
```

Like with `buildLavaan()`, `waves` is just the number of possible waves. The arguments `xWaves` and `yWaves` specify which of the possible waves you have data for. For example, in a 10-wave study where Wave 4 X variable is missing, you could specify `xWaves = c(1:3,5:10)`. If you don't enter anything for these arguments, the code generator assumes that all waves exist.

The 'constrainCor' argument specifies whether to constrain correlations between the same indicator at different waves to be equal for equal-length intervals. The default is to make this constraint, but it can be removed with this argument.

The 'limits' option specifies whether to constrain variances to be greater than 1 and correlations to be between 1 and -1. This is often necessary to get STARTS models to converge. 

The other different options are the directory where the mplus files should be stored and the title of the output (which is also used as the name of the files, so keep it short). I also have default "analysis" and "output" lines for mplus, but you can replace those with the option in the command. If you need to do things like increase the number of iterations, you can do so by changing the 'analysis' option. 

### buildMplus()

If you don't want to install `MplusAutomation` or just want to get the model statement that you can then paste into an Mplus input file, you can just call the `buildMplus()` command. This creates that part of the .inp file that would be specified on the `MODEL=` line. If you wanted it printed nicely from R, you can call: `cat(buildMplus())`. Because `run_starts_mplus()` calls `buildMplus()`, the options are mostly the same; they just omit anything that isn't relevant for the `Model=` line. Here are the options:

```R
buildMplus <- function(waves,
                       XVar = TRUE,
                       YVar = TRUE,
                       xWaves = NULL,
                       yWaves = NULL,
                       xIndicators = 1,
                       yIndicators = 1,
                       trait = TRUE,
                       AR = TRUE,
                       state = TRUE,
                       crossLag = TRUE,
                       stateCor = FALSE,
                       constrainCors = TRUE,
                       stationarity = TRUE)
```

If you want to do something that this code does not do, you can always call one of these two functions and then manually edit the resulting output/files. For instance, you can run `cat(buildMplus(waves=10))`, which prints the model. You can then use Mplus to run the .inp file, or you can use MplusAutomation to run the model. 

If you want to do the latter, you have to follow a couple other steps. First, you create an mplusObject. To do so, you can then write `modelStatement <- " {pasted model code} "`, where {pasted model code} is what you copied from the output of `buildMplus()`. 

If you have constraints in your model, including stationarity constraints or limits on the variances and covariances, you will also have to separately create a constraints statement using the `buildConstraints()` function. It's options are also quite similar to `run_starts_mplus()`:

```R
buildConstraints <- function(waves,
                             XVar = TRUE,
                             YVar = TRUE,
                             crossLag = TRUE,
                             trait = TRUE,
                             AR = TRUE,
                             state = TRUE,
                             stateCor = TRUE,
                             stationarity = TRUE,
                             constrainCor = stationarity,
                             limits = TRUE)
```

If you cat `buildConstraints()` you can edit the code and create a "constraints" object in the same way you created the model statement object.

You then use the code below to create the mplusObject:

```R
### Create the MplusObject
inp <- MplusAutomation::mplusObject(
    Title = "Title",
    rdata = data,
    usevariables = names(data),  ## Or whatever the names of the variables are
    ANALYSIS = "MODEL=NOCOVARIANCES;\nCOVERAGE=.001;", ## Add any options here
    MODEL = model, ## This is the edited model code that you get from buildMplus
	MODELCONSTRAINT = constraints, ## This is the edited constraints code that you get from buildConstraints
    OUTPUT = "stdyx;"
)

```

Finally, you can use the `mplusModeler()` function to actually run the model:

```R
MplusAutomation::mplusModeler(inp, modelout = "mplus/title.inp", run=1)
```

### Wrapper Functions

As with `buildLavaan()`, there are some wrapper functions that can be used to specify other models. Specifically, you can use `run_startsx_mplus()`, `run_startsy_mplus()`, `run_clpm_mplus()`, `run_riclpm_mplus()`, or `run_arts_mplus()` to get those variations on the STARTS model. I also created `run_sts_mplus()`, which runs a stable-trait + state model. These often have a limited set of options compared to `build_starts_mplus()`, given the constraints required to specify these models. For details, look at the code or help files for these models. 

### Comparing Models

> [!NOTE]
> This only works for Mplus output.

There is also one more function that can be used to compare a set of nested models for one variable (this is only available for Mplus right now; the function assumes that all variables are labelled "xw", where w is the wave number). The function `compareUnivariate` takes three arguments: `data` (the data file) `waves` (the number of waves), and `xWaves` (which waves actually exist). It then runs the univariate STARTS, ARTS, START, STS, and ART models. The univariate ART model corresponds to the bivariate CLPM; the univariate START model corresponds to the bivariate RI-CLPM model; the univariate ARTS model corresponds to the bivariate ARTS or factor CLPM model, and the univariate STARTS model corresponds to the bivariate STARTS model. The STS model is rarely used in a bivariate context, but it could be. 

### Dynamic Panel Models

You can also create code for and run dynamic panel models. Because of some of the differences between the DPM and the STARTS-based models, certain options are not available. For instance, it does not make sense to run a DPM without a stable trait (this would just be an ART model that can be specified using the original code). In addition, most implementations of the DPM do not impose strict stationarity, so I removed that option here (though stabilities and cross-lagged paths are constrained to be equal across waves; it is only the variance that is freed). This makes it difficult to use phantom variables because the model is typically not identified with missing waves and no stationarity. Most DPM models do not include state variance, so that component is not included. 

The options for the DPM are:
```R
run_dpm_mplus <- function(data,
                          waves,
                          XVar = TRUE,
                          YVar = TRUE,
                          xIndicators = 1,
                          yIndicators = 1,
                          crossLag = TRUE,
                          constrainCors = TRUE,
                          limits = TRUE,
                          dir = "mplus",
                          title = "dpm",
                          analysis = "MODEL=NOCOVARIANCES;\nCOVERAGE=.001;",
                          output = "stdyx; \n  cinterval; \n")
```

For lavaan, the commands are `buildLavaanDpm()` or `run_dpm_lavaan()` (with similar options to the above; see R documentation for details).


## Utilities

I've included a few utilities that can help you better understand these models. The function `gen_starts()` can generate data based on the STARTS model. As with the functions to test these models, the function to generate data can include or exclude different variance components. See the documentation for details. There is also a function called `addIndicators()` that (as the name suggests) takes the data generated above and adds indicators so you can test the latent-variable models. 

## Issues?

Feel free to create an "Issue" to document any problems and I'll try to get to it.
