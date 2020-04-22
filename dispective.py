from qutip import *
from numpy import pi, sqrt, real
from package.JCM import JCM_Hamiltonian, initial_coherent_state, operator, initial_fock_state
import numpy as np
import matplotlib.pyplot as plt


"""Parameters"""
N = 16

wr = 2.0 * 2 * pi      # resonator frequency
wq = 3.0 * 2 * pi      # qubit frequency
chi = 0.025 * 2 * pi   # parameter in the dispersive hamiltonian >> g**2/delta

delta = abs(wr - wq)        # detuning
g = sqrt(delta * chi)  # coupling strength that is consistent with chi
# compare detuning and g, the first should be much larger than the second

#operators
sm, sz, a, I = operator(N)

nc = a.dag() * a
xc = a + a.dag()

nq = sm.dag() * sm
xq = sm + sm.dag()
sx = (sm + sm.dag())
# dispersive hamiltonian
H = JCM_Hamiltonian(N, wr, wq, g, eff=True)
# H = JCM_Hamiltonian(N, wr, wq, g)

# """
# Try different initial state of the resonator,
# and see how the spectrum further down in the notebook reflects the photon distribution chosen here.
# """
# psi0 = tensor(coherent(N, sqrt(6)), (basis(2,0)+basis(2,1)).unit())
# psi0 = tensor(thermal_dm(N, 3), ket2dm(basis(2,0)+basis(2,1))).unit()
# psi0 = tensor(coherent(N, sqrt(4)), (basis(2,0)+basis(2,1)).unit())
psi0 = initial_coherent_state(N, z=sqrt(4), wav=(1, 1))
# psi0 = initial_fock_state(N, n=0)
plot_fock_distribution(psi0)
plt.show()

tlist = np.linspace(0, 250, 1000)
res = mesolve(H, psi0, tlist, [], [])
# """
# Excitation numbers
# We can see that the systems do not exchange any energy,
# because of they are off resonance with each other.
# """
nc_list = expect(nc, res.states)
nq_list = expect(nq, res.states)
fig, ax = plt.subplots(1, 1, sharex=True, figsize=(12,4))

ax.plot(tlist, nc_list, 'r', linewidth=2, label="cavity")
# ax.plot(tlist, nq_list, 'b--', linewidth=2, label="qubit")
ax.set_ylim(0, 7)
ax.set_xlim(0, 20)
ax.set_ylabel("n", fontsize=16)
ax.set_xlabel("Time (ns)", fontsize=16)
ax.legend()

fig.tight_layout()
plt.show()

xc_list = expect(xc, res.states)
fig, ax = plt.subplots(1, 1, sharex=True, figsize=(12,4))

ax.plot(tlist, xc_list, 'r', linewidth=2, label="cavity")
ax.set_ylabel("x", fontsize=16)
ax.set_xlabel("Time (ns)", fontsize=16)
ax.set_xlim(0,50)
ax.legend()
fig.tight_layout()
plt.show()

tlist = np.linspace(0, 1000, 10000)
corr_vec = correlation(H, psi0, None, tlist, [], a.dag(), a)

fig, ax = plt.subplots(1, 1, sharex=True, figsize=(12,4))

ax.plot(tlist, real(corr_vec), 'r', linewidth=2, label="resonator")
ax.set_ylabel("correlation", fontsize=16)
ax.set_xlabel("Time (ns)", fontsize=16)
ax.legend()
ax.set_xlim(0,50)
fig.tight_layout()
plt.show()

w, S = spectrum_correlation_fft(tlist, corr_vec)
fig, ax = plt.subplots(figsize=(9,3))
ax.plot(w / (2 * pi), abs(S))
ax.set_xlabel(r'$\omega$', fontsize=18)
ax.set_xlim(wr/(2*pi)-.5, wr/(2*pi)+.5)
plt.show()

fig, ax = plt.subplots(figsize=(9,3))
ax.plot((w-wr)/chi, abs(S))
ax.set_xlabel(r'$(\omega-\omega_r)/\chi$', fontsize=18)
ax.set_xlim(-2,2)
plt.show()
