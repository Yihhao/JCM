import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
from qutip import *
from numpy import pi


wc = 2.0 * 2 * pi  # cavity frequency
delta = 0          # detuning
wa = wc + delta    # atom frequency
g = wa             # coupling strength
N = 5
j = 2

sm = tensor(qeye(N), destroy(j))
sz = sm.dag() * sm - sm * sm.dag()  # [S^+, S^-]= S^z
a = tensor(destroy(N), qeye(j))

H0 = wc * (a.dag() * a + 0.5 * sz)
Hint = 0.5 * delta * sz + g * (sm.dag() * a + sm * a.dag())

# plot_energy_levels([H0, Hint], figsize=(8,4))
# plt.show()

# fig, ax = plt.subplots(1, 1, figsize=(8,8))
# xvec = np.linspace(-5,5,200)
# W = wigner(rho, xvec, xvec)
# wlim = abs(W).max()
# ax.contourf(xvec, xvec, W, 100, norm=mpl.colors.Normalize(-wlim, wlim), cmap=plt.cm.get_cmap('RdBu'))
# ax.set_xlabel(r'$x_1$', fontsize=16)
# ax.set_ylabel(r'$x_2$', fontsize=16)
# plt.show()
