import matplotlib.pyplot as plt
import numpy as np
from numpy import sqrt, pi, real
from qutip import *
from package.JCM import JCM_Hamiltonian, initial_coherent_state, operator, initial_fock_state


N = 30                 # number of cavity fock states
z = sqrt(4)            # fock state occupy number of cavity
z = 5
wc = 2.0 * 2 * pi      # cavity frequency
wa = 3.0 * 2 * pi      # atom frequency
chi = 0.025 * 2 * pi   # parameter in the dispersive hamiltonian >> g**2/delta
delta = abs(wc - wa)   # detuning
g = sqrt(delta * chi)  # coupling strength that is consistent with chi
# g = 0.5 * 2 * pi
use_rwa = True        # rwa: rotating wave approximation

H = JCM_Hamiltonian(N, wc, wa, g, use_rwa)
# H = JCM_Hamiltonian(N, wc, wa, g, eff=True)
sm, sz, a, I = operator(N)
# H = wc * (a.dag() * a + I/2.0) + (wa / 2.0) * sz + chi * (a.dag() * a + I/2) * sz

psi0 = initial_coherent_state(N, z, wav=(1, 0))
# psi0 = initial_fock_state(N, n=0)

taulist = np.linspace(0, 1000, 10000)
# collapse operator
kappa = 0.75
n_th = 2.00  # bath temperature in terms of excitation number
c_ops = [sqrt(kappa*(1+n_th)) * a, sqrt(kappa*n_th) * a.dag()]

# start with a coherent state
plot_fock_distribution(psi0)
plt.show()

# first calculate the occupation number as a function of time
# n, xc = mesolve(H, rho0, taulist, c_ops, [a.dag() * a, a.dag() + a]).expect
n, xc = mesolve(H, psi0, taulist, [], [a.dag() * a, a.dag() + a]).expect

# calculate the correlation function G1 and normalize with n to obtain g1
# G1 = correlation(H, rho0, None, taulist, c_ops, a.dag(), a)
G1 = correlation(H, psi0, None, taulist, [], a.dag(), a)
g1 = G1 / sqrt(n[0] * n)

fig, axes = plt.subplots(4, 1, sharex=True, figsize=(12,6))
plt.subplots_adjust(hspace=0.2)

axes[0].plot(taulist, real(g1), 'b', label=r'First-order coherence function $g^{(1)}(\tau)$')
axes[0].set_ylabel("correlation", fontsize=16)
axes[0].legend()

axes[1].plot(taulist, real(xc),  'g', label=r'position $x(\tau)$')
axes[1].set_ylabel("position", fontsize=16)
axes[1].legend()

axes[2].plot(taulist, real(G1),  'm', label=r'position $x(\tau)$')
axes[2].set_ylabel("G1", fontsize=16)
axes[2].legend()

axes[3].plot(taulist, n,  'r', label=r'occupation number $n(\tau)$')
axes[3].set_ylim(0, 5)
axes[3].set_ylabel("n", fontsize=16)
axes[3].set_xlabel(r'$\tau$', fontsize=16)
axes[3].legend()

plt.xlim(0, 50)
plt.tight_layout()
plt.show()

# w1, S1 = spectrum_correlation_fft(taulist, g1)
# fig, ax = plt.subplots(1, 1, figsize=(9, 3))
# ax.plot(w1 / (2 * pi), abs(S1))
# ax.set_xlabel(r'$\omega$', fontsize=18)
# ax.set_xlim(wc/(2*pi)-.5, wc/(2*pi)+.5)
# plt.show()

fig, axes = plt.subplots(1, 1, sharex=True, figsize=(12,6))

axes.plot(taulist, n,  'r', label=r'occupation number $n(\tau)$')
axes.set_ylim(0, 5)
axes.set_ylabel("n", fontsize=16)
axes.set_xlabel(r'$\tau$', fontsize=16)
axes.legend()

plt.xlim(0, 50)
plt.tight_layout()
plt.show()