import numpy as np

"""
 For each Fe3Fe2 model, parameters a, b, c, d for functions fitted to the running standard deviation in the 1 bar validation dataset. The window for the running standard deviation was 30 and the fitted function is a second order polynomial + exponential function:

a+x + b*x + x**2 + c**(x+d)
"""

error_params_1bar = {
    "borisov": np.array([0.07868941, 0.02633911, 0.44011372, 3.51168686]),
    "kressCarmichael": np.array(
        [7.12393565e-02, 2.41635659e-02, 9.87971302e-01, 2.24839355e02]
    ),
    "jayasuriya": np.array(
        [1.38974880e-01, 2.05773269e-02, 9.82551705e-01, 1.50581835e02]
    ),
    "putirka2016_6b": np.array(
        [5.55903859e-02, 5.17449983e-02, 9.85793596e-01, 1.81777331e02]
    ),
    "putirka2016_6c": np.array(
        [-3.41628312e-02, 7.55501213e-02, 9.84711223e-01, 1.53008554e02]
    ),
    "Deng2020": np.array(
        [2.07930660e-01, 3.68177087e-03, 9.85367141e-01, 2.45645691e02]
    ),
    "Oneill2006": np.array(
        [2.40986501e-01, 1.13929610e-02, 9.86720388e-01, 2.13904257e02]
    ),
    "Oneill2018": np.array(
        [1.09493618e-02, 8.16800011e-02, 7.83183748e-01, 1.10172901e01]
    ),
    "Armstrong2019": np.array(
        [1.86137099e-01, 2.98435869e-02, 9.80371940e-01, 1.35011645e02]
    ),
    "Zhang2017": np.array(
        [1.71521144e-01, 5.24650894e-03, 9.85646185e-01, 2.42106215e02]
    ),
    "Hirschmann2022": np.array([0.06473992, 0.02796151, 0.57637993, 4.86026864]),
    "Sun2024": np.array([0.05594865, 0.04124761, 0.48202404, 3.70602711]),
}

"""
For each Fe3Fe2 model, parameters for the spline representatition of the running standard deviation in the > 1 bar validation dataset. The window for the running standard deviation was 30 datapoints and splines were fitted with scipy.iterpolate.splrep with s=3. Resulting splines can be evaluated with scipy.interpolate.splev(xvals, parameters)
"""

error_params_high_pressure = {
    "borisov": (
        np.array(
            [
                0.05263158,
                0.05263158,
                0.05263158,
                0.05263158,
                0.36602451,
                2.16064117,
                2.16064117,
                2.16064117,
                2.16064117,
            ]
        ),
        np.array(
            [
                -0.26851486,
                0.17700856,
                2.82276119,
                0.29727467,
                0.998269,
                0.0,
                0.0,
                0.0,
                0.0,
            ]
        ),
        3,
    ),
    "kressCarmichael": (
        np.array(
            [
                0.05263158,
                0.05263158,
                0.05263158,
                0.05263158,
                2.16064117,
                2.16064117,
                2.16064117,
                2.16064117,
            ]
        ),
        np.array([-0.15803373, 2.4029921, 0.24315807, 0.55458778, 0.0, 0.0, 0.0, 0.0]),
        3,
    ),
    "jayasuriya": (
        np.array(
            [
                0.05263158,
                0.05263158,
                0.05263158,
                0.05263158,
                2.16064117,
                2.16064117,
                2.16064117,
                2.16064117,
            ]
        ),
        np.array([-0.1593344, 2.59665428, -0.05032735, 0.73183981, 0.0, 0.0, 0.0, 0.0]),
        3,
    ),
    "putirka2016_6b": (
        np.array(
            [
                0.05263158,
                0.05263158,
                0.05263158,
                0.05263158,
                0.36602451,
                2.16064117,
                2.16064117,
                2.16064117,
                2.16064117,
            ]
        ),
        np.array(
            [
                -0.19827421,
                0.15232086,
                2.23039298,
                0.09093898,
                0.71149013,
                0.0,
                0.0,
                0.0,
                0.0,
            ]
        ),
        3,
    ),
    "putirka2016_6c": (
        np.array(
            [
                0.05263158,
                0.05263158,
                0.05263158,
                0.05263158,
                2.16064117,
                2.16064117,
                2.16064117,
                2.16064117,
            ]
        ),
        np.array([-0.19486923, 2.25964404, 0.62214883, 0.51085643, 0.0, 0.0, 0.0, 0.0]),
        3,
    ),
    "Deng2020": (
        np.array(
            [
                0.05263158,
                0.05263158,
                0.05263158,
                0.05263158,
                2.16064117,
                2.16064117,
                2.16064117,
                2.16064117,
            ]
        ),
        np.array(
            [-0.05797588, 1.26548832, -0.22312187, 0.41829024, 0.0, 0.0, 0.0, 0.0]
        ),
        3,
    ),
    "Oneill2006": (
        np.array(
            [
                0.05263158,
                0.05263158,
                0.05263158,
                0.05263158,
                2.16064117,
                2.16064117,
                2.16064117,
                2.16064117,
            ]
        ),
        np.array([0.09243536, 0.81805066, -0.00252054, 0.75367289, 0.0, 0.0, 0.0, 0.0]),
        3,
    ),
    "Oneill2018": (
        np.array(
            [
                0.05263158,
                0.05263158,
                0.05263158,
                0.05263158,
                0.36602451,
                2.16064117,
                2.16064117,
                2.16064117,
                2.16064117,
            ]
        ),
        np.array(
            [
                -0.17439819,
                -0.01519068,
                3.92851877,
                0.09896842,
                1.23105731,
                0.0,
                0.0,
                0.0,
                0.0,
            ]
        ),
        3,
    ),
    "Armstrong2019": (
        np.array(
            [
                0.05263158,
                0.05263158,
                0.05263158,
                0.05263158,
                2.16064117,
                2.16064117,
                2.16064117,
                2.16064117,
            ]
        ),
        np.array([0.0936293, 0.81094331, 0.01593797, 0.79183567, 0.0, 0.0, 0.0, 0.0]),
        3,
    ),
    "Zhang2017": (
        np.array(
            [
                0.05263158,
                0.05263158,
                0.05263158,
                0.05263158,
                2.16064117,
                2.16064117,
                2.16064117,
                2.16064117,
            ]
        ),
        np.array([-0.01421484, 0.76319884, 0.00203087, 0.60638093, 0.0, 0.0, 0.0, 0.0]),
        3,
    ),
    "Hirschmann2022": (
        np.array(
            [
                0.05263158,
                0.05263158,
                0.05263158,
                0.05263158,
                2.16064117,
                2.16064117,
                2.16064117,
                2.16064117,
            ]
        ),
        np.array([0.02126475, 1.02796781, -0.27739497, 0.56921383, 0.0, 0.0, 0.0, 0.0]),
        3,
    ),
    "Sun2024": (
        np.array(
            [
                0.05263158,
                0.05263158,
                0.05263158,
                0.05263158,
                2.16064117,
                2.16064117,
                2.16064117,
                2.16064117,
            ]
        ),
        np.array([0.01892744, 1.14826304, -0.40097838, 0.41416629, 0.0, 0.0, 0.0, 0.0]),
        3,
    ),
}
