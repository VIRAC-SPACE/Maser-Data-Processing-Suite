from dataphile.statistics.regression.modeling import Parameter, Model, CompositeModel, AutoGUI
from dataphile.datasets import SyntheticDataset
from dataphile.statistics.distributions import polynomial1D
import numpy as np

m = Parameter(value=2, bounds=(1,3), label='slope')
b = Parameter(value=1, bounds=(0,2), label='intercept')

xdata, ydata = SyntheticDataset(polynomial1D, [100, -0.01, -1e-5], (0, 2400), linspace=True, noise=0, samples=2400).generate()

def linear(x, *p):
    return p[0] + p[1] * x

model = Model(linear, b, m)
model.fit(xdata, ydata)


A     = Parameter(value=1,   bounds=(1/4, 2), label='amplitude')
x0    = Parameter(value=1,   bounds=(0, 2),   label='centroid')
sigma = Parameter(value=1/2, bounds=(1/4, 1), label='width')

def gaussian(x, *p):
    return p[0] * np.exp(-0.5 * (x - p[1])**2 / p[2]**2)

model = CompositeModel(Model(linear, b, m),
                       Model(gaussian, A, x0, sigma),
                       label='gaussian_with_background')
model.fit(xdata, ydata)