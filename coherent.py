import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
from qutip import *

N = 20

rho_coherent = coherent_dm(N, np.sqrt(2))
rho_thermal = thermal_dm(N, 2)
rho_fock = fock_dm(N, 2)

fig, axes = plt.subplots(1, 3, figsize=(12,3))
bar0 = axes[0].bar(np.arange(0, N)-.5, rho_coherent.diag())
lbl0 = axes[0].set_title("Coherent state")
lim0 = axes[0].set_xlim([-.5, N])
bar1 = axes[1].bar(np.arange(0, N)-.5, rho_thermal.diag())
lbl1 = axes[1].set_title("Thermal state")
lim1 = axes[1].set_xlim([-.5, N])
bar2 = axes[2].bar(np.arange(0, N)-.5, rho_fock.diag())
lbl2 = axes[2].set_title("Fock state")
lim2 = axes[2].set_xlim([-.5, N])
plt.show()

fig, axes = plt.subplots(1, 3, figsize=(12,3))
plot_fock_distribution(rho_coherent, fig=fig, ax=axes[0], title="Coherent state")
plot_fock_distribution(rho_thermal, fig=fig, ax=axes[1], title="Thermal state")
plot_fock_distribution(rho_fock, fig=fig, ax=axes[2], title="Fock state")
fig.tight_layout()
plt.show()

rho1 = coherent(15, -3)
rho2 = coherent(15, 2)
rho3 = coherent(15, 2j)

# fig, axes = plt.subplots(3, 2, figsize=(8, 12))
# plot_wigner_fock_distribution(rho1, fig=fig, axes=axes[0])
# plot_wigner_fock_distribution(rho2, fig=fig, axes=axes[1])
# plot_wigner_fock_distribution(rho3, fig=fig, axes=axes[2])
# fig.tight_layout()
# plt.show()

fig, axes = plt.subplots(4, 2, figsize=(8, 16))
plot_wigner_fock_distribution(rho2, alpha_max=3, fig=fig, axes=axes[0])
plot_wigner_fock_distribution(rho2, fig=fig, axes=axes[1])
plot_wigner_fock_distribution(rho3, alpha_max=3, fig=fig, axes=axes[2])
plot_wigner_fock_distribution(rho3, fig=fig, axes=axes[3])
plt.savefig('./test')
plt.show()
# fig, ax = plt.subplots(1, 1, figsize=(8,8))
# xvec = np.linspace(-5,5,200)
# W = wigner(rho1, xvec, xvec)
# wlim = abs(W).max()
# ax.contourf(xvec, xvec, W, 100, norm=mpl.colors.Normalize(-wlim, wlim), cmap=plt.cm.get_cmap('RdBu'))
# ax.set_xlabel(r'$x_1$', fontsize=16)
# ax.set_ylabel(r'$x_2$', fontsize=16)
# plt.show()
