# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this package does

`panelCodeR` is an R package that generates and runs lavaan and Mplus model code for panel data analysis. It supports a family of models rooted in the STARTS framework (Stable Trait, Autoregressive Trait, State): STARTS, RI-CLPM, ARTS, STS, CLPM, DPM (two variants), GCLM, LGCM, ALT, LCM-SR, and a measurement model. Models can be univariate or bivariate, with optional multiple indicators per wave and phantom variables for missing waves.

## Common commands

All development tasks are done interactively in R. There is no build script; use devtools:

```r
# Load package for interactive development
devtools::load_all()

# Build documentation from roxygen2 comments
devtools::document()

# Install locally
devtools::install()

# Run the informal test scripts
source("tests/tests.R")
```

The test file `tests/tests.R` is a script (not a formal testthat suite) that generates synthetic data with `gen_starts()` and runs models through `panelcoder()`.

## Architecture

### Data flow through the package

User data → `panelcoder()` → `.buildModel()` → `.buildTable()` (builds lavaan parameter table) → constraint functions → `lav2lavaan()` or `lav2mplus()` → lavaan/MplusAutomation → `.summarizeLavaan()` / `.summarizeMplus()` → returned `pcOutput` object.

The central intermediate representation is a **parameter table** (a data frame with columns `lhs`, `op`, `rhs`, `user`, `block`, `group`, `free`, `ustart`, `exo`, `label`) — the same format lavaan uses internally. All model construction operates on this table. Converting to lavaan syntax (`lav2lavaan()`) or Mplus syntax (`lav2mplus()`) happens at the end.

### Key files

- **`R/panelCodeR.R`** — the `panelcoder()` user-facing entry point and `.buildModel()` internal dispatcher. `.buildModel()` maps each `panelModel` string to a set of boolean flags (trait, ar, state, crossLag, dpm_c, etc.) and passes them to `.buildTable()`.
- **`R/buildParTable.R`** — `.buildTable()` constructs the initial parameter table from data info. Also contains `.buildObserved()` (measurement model rows) and related helpers.
- **`R/buildModel.R`** — string-based lavaan code builders (legacy path; still used by `buildLavaan()`, `lavaanStarts2()`, etc.). These build lavaan model strings directly rather than via parameter tables.
- **`R/buildDpm.R`** — string-based Mplus code builders for the DPM variants (`buildDpmLavaan.R` does the lavaan version of DPM).
- **`R/buildConstraints.R`** — post-processing constraint functions applied to the parameter table: `.constrainStability()`, `.constrainCl()`, `.constrainStateVar()`, `.constrainLoadings()`, `.constrainResidVar()`, `.constrainDpmLoadings()`, `.arStationarity()`, `.buildLimits()`, etc.
- **`R/fromLavaan.R`** — `lav2mplus()` and `lav2lavaan()` translate the parameter table to final model syntax strings.
- **`R/utils.R`** — `getInfo()` parses the data frame column names to extract variable names, wave numbers, and indicator counts. This is the first thing called in `panelcoder()` and its output (`info`) is threaded through most internal functions.
- **`R/compareModels.R`** — `compareModels()` loops `panelcoder()` over a list of model names and returns a comparison table.

### Data naming convention (critical)

`getInfo()` infers everything from column names split on `_`:
- Single indicator: `x_1`, `x_2`, `x_3` (stem `x`, wave index)
- Multiple indicators: `x_1_1`, `x_1_2`, `x_2_1`, `x_2_2` (stem, wave, indicator)
- Bivariate: columns with two distinct stems, e.g. `x_1`, `x_2`, `y_1`, `y_2`
- Mplus restricts names to 8 characters; keep stems very short when targeting Mplus
- The data frame must contain **only** the model variables (no id columns)

### Return object structure

`panelcoder()` returns a list of class `pcOutput` with 8 elements:
1. `pcSum` — summary object (class `pcmObject` or `pclObject` depending on program)
2. `info` — list from `getInfo()`: `$gen`, `$x`, `$y` sub-lists with wave/indicator metadata
3. `model` — the parameter table used to build the model
4. `fit` — the raw lavaan or MplusAutomation fit object
5. `corSummary` — lag correlations for plotting
6. `modelCode` — the model code string
7. `warningM` — captured warning message (or NULL)
8. `errorM` — captured error message (or NULL)

Access the raw lavaan object via `pcOutput[[4]]`; use `modelCode(pcOutput)` to print model syntax.

### Model variants and the panelModel flag

The `panelModel` argument in `panelcoder()` maps to these strings:
`"starts"`, `"riclpm"` / `"start"`, `"arts"`, `"clpm"` / `"art"`, `"sts"`, `"dpm_c"`, `"dpm_p"`, `"alt"`, `"gclm"`, `"lgcm"`, `"lcmsr"`, `"measurement"`

Key constraints: lags > 1 are not allowed with random-intercept models; DPM/GCLM cannot use phantom variables; LGCM requires `slope != "none"`.
