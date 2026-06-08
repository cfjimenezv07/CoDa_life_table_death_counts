<br />
<div align="center">
  <a href="https://github.com/cfjimenezv07/Density_Functional_Panel">
    <img src="UoY.png" alt="York Logo" height="150">
    <img src="KAUST.png" alt="KAUST Logo" height="150">
    <img src="MQ.png" alt="MQ Logo" height="150">
  </a>

<h3 align="center">Forecasting Density-Valued Functional Panel Data</h3>
</div>

## Abstract
<p align="justify">
We introduce a statistical method for modelling and forecasting functional panel data represented by multiple densities. Density functions are non-negative and have a constrained integral, and thus do not constitute a linear vector space. We implement a centre log-ratio transformation to transform densities into unconstrained functions. These functions exhibit cross-sectional correlation and temporal dependence. Via a functional analysis-of-variance decomposition, we decompose the unconstrained functional panel data into a deterministic trend component and a time-varying residual component. To produce forecasts for the time-varying component, a functional time series forecasting method, based on the estimation of the long-run covariance, is implemented. By combining the forecasts of the time-varying residual component with the deterministic trend component, we obtain $h$-step-ahead forecast curves for multiple populations. Illustrated by age- and sex-specific life-table death counts in the United States, we apply our proposed method to generate forecasts of the life-table death counts for 51 states.
</p>

### Main Results
The R script files in the `R Code` folder should be used in the following order:

#### 1. Setup, Package Initialization, and Data Infrastructure
* **MITS_class.R**: Defines custom object structures and classes to organize and manage high-dimensional or multivariate/panel functional time-series data streams.
* **CoDa_transformations.R**: Handles the mathematical center log-ratio (clr) mappings and their corresponding inverse transformations to map density curves into an unconstrained space.
* **method.FPE.R**: Selects the optimal number of functional principal components ($K$) to retain using the Functional Prediction Error (FPE) criterion.

#### 2. Covariance Estimation and Core Forecasting Methods
* **New_cov_all.R**: Calculates panel covariance surfaces across the different cross-sections to rigorously account for cross-sectional dependence.
* **forecast_Arima.R**: Fits univariate ARIMA models to time-varying functional principal component scores to project them into future horizons.
* **nonparametric_fof_regression.R**: Executes function-on-function non-parametric regressions.
* **New_point_forecast2.R**: Serves as the primary master workflow that orchestrates the execution and extracts the out-of-sample point-forecasted density curves.

#### 3. Uncertainty Quantification and Joint Prediction Bands
* **CoDa_nonparametric_boot2.R**: Implements the non-parametric bootstrap routine.
* **sieve_code.R**: Implements the functional sieve bootstrap algorithm designed for generating robust non-parametric bootstrap prediction bands.
* **Aux_Uniform_prediction_band.R**: Provides supporting backend mathematical functions needed to construct joint uniform prediction bands for density curves.
* **Uniform_Prediction_bands.R**: Constructs simultaneous, continuous confidence envelopes around the forecasted density paths.

#### 4. Baseline Benchmarking and Alternative Models
* **Alternative_methods_2.R**: Evaluates baseline alternative forecasting models to compare against the paper's proposed framework.
* **Alternative_methods_wo_clr.R**: Runs alternative benchmarking forecasting models directly on the data without applying the Center Log-Ratio transformation.
* **Aux_MEM.R**: Contains auxiliary functions for computing forecast errors using the Maximum Entropy Mortality (MEM) model framework of Pascariu, Lenart, and Canudas-Romo (2019).
* **MEM19.R**: Runs specific configurations or legacy versions of the underlying Maximum Entropy Mortality framework.
* **model.MEM.R**: Contains the core statistical framework and algorithmic logic for fitting the primary Maximum Entropy Mortality factor model.
* **aux_TNH.R**: Serves as the auxiliary helper suite for the High-Dimensional Functional Time Series (HDFTS) method established by Tavakoli, et al. (2022) in the *Journal of Time Series Analysis* (JTSA).

#### 5. Master Empirical Execution, Error Evaluation, and Visualization
* **CoDa_case_USA.R**: Functions as the master execution script that runs the entire empirical analysis on the US mortality dataset.
* **forecast_errors.R**: Computes structural error metrics like integrated squared error to quantify the accuracy of the density forecasts.
* **Comparisons_actual_cruves.R**: Evaluates the differences and visual match between the true observed density curves and the out-of-sample forecasts.
* **Gini_plots.R**: Generates visual plots analyzing lifespan inequality trends via Gini coefficients calculated from the mortality distributions.

### Shiny Application
The directory contains a **`Shiny_app`** folder housing all the R code, UI layout, and server logic required to run the interactive dashboard presented in the paper, allowing users to dynamically visualize density forecasts, prediction bands, and comparative models.

## How to Cite
If you use this code or methodology in your research, please cite our paper:

**APA Style:**
> Jiménez-Varón, C. F., Sun, Y., & Shang, H. L. (2025). Forecasting Density-Valued Functional Panel Data. *Australian & New Zealand Journal of Statistics*, 67(3), 401–415. https://doi.org/10.1111/anzs.70013

**BibTeX:**
```bibtex
@article{jimenezvaron2025forecasting,
  author    = {Jim{\'e}nez-Var{\'o}n, Cristian F. and Sun, Yanrong and Shang, Han Lin},
  title     = {Forecasting Density-Valued Functional Panel Data},
  journal   = {Australian \& New Zealand Journal of Statistics},
  volume    = {67},
  number    = {3},
  pages     = {401--415},
  year      = {2025},
  doi       = {10.1111/anzs.70013},
  url       = {[https://doi.org/10.1111/anzs.70013](https://doi.org/10.1111/anzs.70013)}
}
