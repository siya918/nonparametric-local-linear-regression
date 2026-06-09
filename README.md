# nonparametric-local-linear-regression
Local linear regression smoother implemented from scratch in SAS IML  and Python with Leave-One-Out Cross Validation for optimal bandwidth  selection. Includes both implementations for direct methodology comparison.


# Non-Parametric Local Linear Regression

## Overview
Implementation of a non-parametric local linear regression smoother 
built from scratch in both SAS IML and Python. At each focal point in 
the dataset, the algorithm selects the m nearest neighbours, fits a 
local ordinary least squares regression, and extracts the fitted value 
at the focal point. Repeating this across all observations produces a 
smooth non-parametric estimate of the underlying regression function.

Leave-One-Out Cross Validation (LOOCV) is used to select the optimal 
bandwidth m — the number of nearest neighbours that minimises 
out-of-sample prediction error.

This approach makes no assumptions about the global functional form of 
the relationship between x and y — the fit is driven entirely by local 
data structure.

---

## What is Local Linear Regression?
Standard parametric regression assumes a fixed global functional form 
such as y = β₀ + β₁x and estimates a single set of coefficients across 
the entire dataset. This works well when the true relationship is 
globally linear but fails when the relationship changes across the range 
of x.

Local linear regression relaxes this assumption entirely. Instead of 
one global model, it fits a separate linear model in the neighbourhood 
of each focal point using only the nearest observations. The fitted 
value at each point comes from its own local model, producing a smoother 
that adapts to whatever functional form the data exhibits locally.

**Algorithm at each focal point xᵢ:**
1. Compute the absolute distance from xᵢ to all other observations
2. Select the m nearest neighbours
3. Fit OLS regression using only those m neighbours
4. Extract the intercept — this is the fitted value at xᵢ since the 
   local distance x0 = 0 at the focal point
5. Repeat for all n observations

**The local OLS solution:**
Where X is the design matrix of the m nearest neighbours with an 
intercept column. The fitted value at xᵢ is bh[0].

---

## Why This Model Matters

### 1. No Functional Form Assumptions
Parametric models require you to specify the relationship before seeing 
the data — linear, quadratic, exponential, and so on. If you guess 
wrong, your model is misspecified and all inference is unreliable. Local 
linear regression lets the data speak for itself. It will fit a linear 
relationship where one exists, a curve where the relationship bends, and 
a flat line where there is no signal — all without being told what to 
look for.

### 2. Handles Complex Non-Linear Relationships
Many real-world relationships are non-linear and non-monotonic — they 
increase in some regions and decrease in others. A global polynomial 
can approximate this but is sensitive to the degree chosen and behaves 
poorly at the boundaries. Local linear regression handles boundary 
regions naturally and adapts to any shape.

### 3. Interpretability
Unlike black-box models such as neural networks or random forests, local 
linear regression is fully transparent. At any focal point you can 
inspect the local slope and intercept, understand which observations 
influenced the fit, and explain the prediction in terms of a simple 
linear model. This makes it particularly valuable in regulated 
industries where model explainability is required.

### 4. Bandwidth Selection is Principled
The choice of m (bandwidth) directly controls the bias-variance 
tradeoff. Too small and the fit is noisy. Too large and the fit is 
oversmoothed. LOOCV provides a statistically principled, data-driven 
method for selecting m rather than relying on arbitrary choices.

### 5. Foundation for Advanced Methods
Local linear regression is the conceptual foundation for kernel 
regression, LOESS/LOWESS smoothing, Gaussian process regression, and 
spline-based methods. Understanding it from scratch builds the 
intuition needed for these more advanced non-parametric techniques.

---

## Where This Model is Used

### Economics and Econometrics
- Estimating non-linear wage-experience profiles
- Modelling consumption functions without assuming linearity
- Regression discontinuity designs where local linearity is assumed 
  at the cutoff
- Density estimation for income and wealth distributions

### Finance
- Estimating non-linear relationships between risk factors and returns
- Volatility surface smoothing in options pricing
- Credit risk modelling where default probability curves are non-linear
- Yield curve smoothing in fixed income markets

### Public Health and Epidemiology
- Dose-response curve estimation without parametric constraints
- Age-specific disease incidence smoothing
- Environmental exposure modelling where effects are non-linear
- Growth curve analysis in child development studies

### Engineering and Signal Processing
- Sensor data smoothing where the underlying signal is unknown
- Calibration curve fitting for measurement instruments
- Fault detection where normal operating curves are non-parametric
- Time series trend extraction

### Climate and Environmental Science
- Temperature trend estimation over time
- Species distribution modelling
- Pollution concentration surface estimation
- Rainfall pattern smoothing across geographic regions

### Social Sciences
- Survey response smoothing
- Educational outcome modelling across demographic groups
- Crime rate surface estimation
- Polling data trend extraction

---

## Bandwidth Selection — LOOCV
The bandwidth m controls how many nearest neighbours are used at each 
focal point. Selecting m is the most important methodological decision 
in local linear regression.

**Leave-One-Out Cross Validation procedure:**
1. Remove one observation from the dataset
2. Fit the local linear regression on the remaining n-1 observations
3. Predict the removed observation using the fitted model
4. Record the squared prediction error
5. Repeat for all n observations
6. Compute mean squared error across all leave-one-out predictions
7. Repeat for all candidate values of m
8. Select the m that minimises LOOCV MSE

**Bias-variance tradeoff:**

| m | Bias | Variance | Fit |
|---|---|---|---|
| Too small | Low | High | Noisy, overfit |
| Optimal | Balanced | Balanced | Best generalisation |
| Too large | High | Low | Oversmoothed, underfit |

The LOOCV error curve is U-shaped — the bottom of the U identifies 
the optimal bandwidth.

---

## Dataset
Simulated dataset of 100 observations with a non-linear data generating 
process:

The sinusoidal signal with additive linear trend and noise provides a 
meaningful test case — a global linear model fits this data poorly while 
local linear regression recovers the underlying non-linear trend.

---

## Outputs

### Graph 1 — Raw Scatter Plot
Visualises the raw observed data points before fitting, showing the 
noisy non-linear relationship between x and y.

### Graph 2 — Observed vs Fitted Overlay
Overlays the local linear regression fitted values on the observed data 
points, showing how the smoother adapts to local data structure without 
imposing a global functional form.

### Graph 3 — LOOCV MSE vs Bandwidth
U-shaped curve showing cross validation error across candidate bandwidth 
values. The minimum identifies the optimal m.

### Graph 4 — Initial vs Optimal Bandwidth Comparison
Overlays the fit at the initial bandwidth (m=30) against the 
CV-optimal fit, showing the practical difference bandwidth choice makes.

---

## Implementation Details

### SAS vs Python Comparison

| Step | SAS | Python |
|---|---|---|
| Random seed | `call streaminit(42)` | `np.random.seed(42)` |
| Uniform draws | `rand('uniform') * 10` | `np.random.uniform(0, 10, n)` |
| Normal noise | `rand('normal', 0, 1.2)` | `np.random.normal(0, 1.2, n)` |
| Read data | `read all into xy` | `df[['x','y']].values` |
| Sort by distance | `call sort(axy, {1})` | `axy[axy[:,0].argsort()]` |
| OLS solution | `inv(x'*x)*x'*y` | `np.linalg.inv(X.T @ X) @ X.T @ y` |
| Scatter plot | `proc gplot / interpol=none` | `plt.scatter` |
| Line plot | `proc gplot / interpol=line` | `plt.plot` |

Both implementations produce identical results from the same random 
seed and algorithm logic.

---
