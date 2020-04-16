import matplotlib.pyplot as plt
from numpy import sqrt
import numpy as np
from qutip import *

w  = 1.0
w0 = 1.0

g = 1.0
gc = sqrt(w * w0)/2 # critical coupling strength

kappa = 0.05
gamma = 0.15
M = 2
N = 4
j = N/2.0
n = 2*j + 1

a  = tensor(destroy(M), qeye(n))
Jp = tensor(qeye(M), jmat(j, '+'))
Jm = tensor(qeye(M), jmat(j, '-'))
Jz = tensor(qeye(M), jmat(j, 'z'))

H0 = w * a.dag() * a + w0 * Jz
H1 = 1.0 / sqrt(N) * (a + a.dag()) * (Jp + Jm)
H = H0 + g * H1

# fig, ax = plt.subplots(1, 1, figsize=(10,10))
# hinton(H, ax=ax)
# plt.show()