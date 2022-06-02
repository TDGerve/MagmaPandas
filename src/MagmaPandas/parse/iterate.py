import numpy as np

class iterate():

    def __init__(self, initial_guess, converge, iterate_var):
        self.initial = initial_guess
        self.converge = converge
        self.iterate_var = iterate_var

    def __call__(self, func):
        def wrapper(**kwargs):
            results_guess = kwargs[self.iterate_var]
            results_equilibrium = func(**kwargs)
            results_delta = abs(results_guess - results_equilibrium) / results_guess
            iterate = results_delta > self.converge

            while sum(iterate) > 1:
                kwargs[self.iterate_var] = results_equilibrium
                kwargs = {key: val[iterate] for key, val in kwargs.items()}
                results_equilibrium = func(**kwargs)
                results_delta = abs(results_guess - results_equilibrium) / results_guess
                iterate = results_delta > self.converge

            return results_equilibrium

        return wrapper

