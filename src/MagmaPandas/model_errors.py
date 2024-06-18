import numpy as np
import scipy.optimize as opt
from scipy import interpolate


def _error_func(data, a, b, c, d):
    """
    second order polynomial + exponential function.
    """
    return a * data + b * data**2 + c ** (data + d)


def _reject_outliers(data: np.ndarray, sigma: int = 2):
    """
    Filter data outside ``sigma`` standard deviations
    """
    return data[abs(data - np.median(data)) < (sigma * np.std(data))]


def _running_stddev(x, y, boxsize: int = 30, sigma: int = 2):
    """
    Calculate the standard deviations in a moving window with size ``boxsize``.
    Data outside ``sigma`` standard deviations are filtered.
    """
    # remove nan's and inf's
    idx = ~(np.isinf(y) | np.isnan(y))
    ydata = y[idx].copy()
    xvals = x[idx].copy()

    fr = boxsize // 2
    to = len(xvals) - fr

    xvals = xvals[fr:-fr]

    stddev = np.array(
        [
            np.std(_reject_outliers(ydata[j - fr : j + fr - 1], sigma=sigma))
            for j in np.arange(fr, to, 1)
        ]
    )

    return xvals, stddev


def calculate_error_parameters(y_true, y_model, boxsize=30, sigma=2, func=_error_func):
    """
    Fit function ``func`` to the running standard deviation of a dataset.
    """
    xvals, stddev = _running_stddev(y_true, y_model, boxsize=boxsize, sigma=sigma)
    f_results = opt.curve_fit(f=func, xdata=xvals, ydata=stddev)
    error_params = f_results[0]

    xlimits = xvals.min(), xvals.max()

    return error_params, xlimits


def calculate_spline_parameters(y_true, y_model, boxsize=30, sigma=2, s=3, fix0=True):

    xvals, stddev = _running_stddev(y_true, y_model, boxsize=boxsize, sigma=sigma)
    weights = np.onex(len(xvals))
    if fix0:
        weights[0] = 10
    params = interpolate.splrep(xvals, stddev, w=weights, s=s)

    xlimits = xvals.min(), xvals.max()

    return params, xlimits
