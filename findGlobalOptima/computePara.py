import numpy as np


E = 1.8e5
h = 1.6e-3

kb = E * np.pi * h**4/4
ks = E * np.pi * h**2

rho = 1180
g = 10

Lgb = E * h**2/(8*rho * g)
Lgb = Lgb**(1/3)

normKs = ks * Lgb**2/kb

ls = 1

print(ls/Lgb)

print(normKs, Lgb)
