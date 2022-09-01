import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
import statsmodels.api as sm
from statsmodels.tools.eval_measures import rmse
from sklearn.model_selection import cross_validate
from sklearn.linear_model import LinearRegression
from sklearn.base import BaseEstimator, RegressorMixin


# A wrapper for statsmodel, for use within sklearn
# Not really needed since LinearRegression produces the same results as sm.OLS
class statsmodel_wrapper(BaseEstimator, RegressorMixin):
    """A universal sklearn-style wrapper for statsmodels regressors"""

    def __init__(self, model_class, fit_intercept=True):
        self.model_class = model_class
        self.fit_intercept = fit_intercept

    def fit(self, X, y):
        if self.fit_intercept:
            X = sm.add_constant(X)
        self.model_ = self.model_class(y, X)
        self.results_ = self.model_.fit()
        return self

    def predict(self, X):
        if self.fit_intercept:
            X = sm.add_constant(X)
        return self.results_.predict(X)
        return self.result.predict(X)


# ols fitting
def ols_fit(x, y):
    """
    Multiple linear ordinary least squares regresson
    """
    x_train = sm.add_constant(x)
    reg_ols = sm.OLS(y, x_train).fit()
    return reg_ols


# F test for models comparison
def f_resid(df, reg1, reg2):
    """
    F test for models comparison
    """
    return (
        abs((reg1.resid**2).sum() - (reg2.resid**2).sum())
        / (len(reg2.params) - len(reg1.params))
        / ((reg2.resid**2).sum() / (len(df) - len(reg2.params)))
    )


def f_to_p(df, reg1, reg2):
    """
    convert F-value to associated P-value
    """
    f = f_resid(df, reg1, reg2)
    return stats.f.cdf(
        f, abs(len(reg1.params) - len(reg2.params)), len(df) - len(reg2.params)
    )


def compare_models(x, y):
    """
    Regress y on x-1 and compare with regressions of y on x
    """
    # Regress y on x
    upper_reg = ols_fit(x, y)

    results = pd.DataFrame(
        columns=["f_ratio", "p_value", "r_ratio"],
        index=np.arange(0, x.shape[1], 1),
        dtype=float,
    )
    parameters = {i: None for i in results.index}
    models = {i: None for i in results.index}

    # Regress y on x-1 for all permutations
    for col_drop in results.index:
        x_reduced = x.drop(x.columns[col_drop], axis=1)
        models[col_drop] = ols_fit(x_reduced, y)
        # Compare with y regressed on x
        results.loc[col_drop, "f_ratio"] = f_resid(x, models[col_drop], upper_reg)
        results.loc[col_drop, "p_value"] = f_to_p(x, models[col_drop], upper_reg)
        results.loc[col_drop, "r2"] = models[col_drop].rsquared

        parameters[col_drop] = pd.concat(
            [models[col_drop].params, models[col_drop].pvalues], axis=1
        )
        parameters[col_drop].columns = ["coefficients", "p_value"]

    return results, parameters, models


def get_model(x, y):
    """
    Compare regression results and select the statistically best model
    """
    # Get models
    results, parameters, models = compare_models(x, y)

    # Select the model with the lowest p vaue
    max_index = results["p_value"].idxmin()
    x_data = x.drop(x.columns[max_index], axis=1)
    return (
        x_data,
        parameters[max_index],
        models[max_index],
    )


def optimise_regression(data: pd.DataFrame, y_predict: str, exclude=None, split=0.15):
    if exclude:
        remove = [y_predict] + exclude
    else:
        remove = [y_predict]
    x = data.drop(columns=remove)
    y = data[y_predict].copy()

    parameters_total = pd.DataFrame(dtype=float, columns=["const"] + list(x.columns))
    errors = pd.DataFrame(dtype=float, columns=["calibration", "validation", "delta"])

    for _ in x.columns[1:]:
        x_data, parameters, model = get_model(x, y)

        # Copy coefficients
        n_parameters = len(x_data.columns)
        parameters_total.loc[n_parameters, parameters.index] = parameters[
            "coefficients"
        ]

        # Calculate errors
        errors.loc[n_parameters, "calibration"] = rmse(
            y, model.predict(sm.add_constant(x_data))
        )
        cross_validation = cross_validate(
            statsmodel_wrapper(sm.OLS),  # LinearRegression(),
            x_data,
            y,
            cv=int(1 / split),
            scoring=("neg_mean_squared_error"),
        )
        errors.loc[n_parameters, "validation"] = np.sqrt(
            abs(cross_validation["test_score"])
        ).mean()

        x = x_data.copy()

        # if (parameters["p_value"] > 0.0001).any() or (r2 < 0.96):
        #     break
    errors["delta"] = abs(errors["validation"] - errors["calibration"])

    return parameters_total, errors
