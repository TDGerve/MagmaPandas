import numpy as np

class iterate():

    def __init__(self, initial_guess, converge, iterate_var):
        self.initial = initial_guess
        self.converge = converge
        self.iterate_var = iterate_var

    def __call__(self, func):
        def wrapper(*args):
            results_guess = np.repeat(self.initial, len(args[0]))
            results_equilibrium = func(*args)
            results_delta = abs(results_guess - results_equilibrium) / results_guess
            iterate = results_delta > self.converge

            while sum(iterate) > 1:
                args = [arg[iterate] for arg in args]
