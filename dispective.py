from qutip import *
from numpy import pi, sqrt, real
from package.JCM import JCM_Hamiltonian
import numpy as np
import matplotlib.pyplot as plt


"""Parameters"""
N = 2

wr = 2.0 * 2 * pi      # resonator frequency
wq = 3.0 * 2 * pi      # qubit frequency
chi = 0.025 * 2 * pi   # parameter in the dispersive hamiltonian >> g**2/delta

delta = abs(wr - wq)        # detuning
g = sqrt(delta * chi)  # coupling strength that is consistent with chi
# compare detuning and g, the first should be much larger than the second

# cavity operators
a = tensor(destroy(N), qeye(2))
nc = a.dag() * a
xc = a + a.dag()

# atomic operators
sm = tensor(qeye(N), destroy(2))
sz = tensor(qeye(N), sigmaz())
sx = tensor(qeye(N), sigmax())
nq = sm.dag() * sm
xq = sm + sm.dag()

I = tensor(qeye(N), qeye(2))
# dispersive hamiltonian
H0, H1 = JCM_Hamiltonian(N, wr, wq, g, eff=True)
H = H0 + H1
# H = wr * (a.dag() * a + I/2.0) + (wq / 2.0) * sz + chi * (a.dag() * a + I/2) * sz

"""
Try different initial state of the resonator, 
and see how the spectrum further down in the notebook reflects the photon distribution chosen here.
"""
#psi0 = tensor(coherent(N, sqrt(6)), (basis(2,0)+basis(2,1)).unit())
#psi0 = tensor(thermal_dm(N, 3), ket2dm(basis(2,0)+basis(2,1))).unit()
psi0 = tensor(coherent(N, sqrt(4)), (basis(2,0)+basis(2,1)).unit())

tlist = np.linspace(0, 250, 1000)
res = mesolve(H, psi0, tlist, [], [], options=Odeoptions(nsteps=5000))

"""
Excitation numbers
We can see that the systems do not exchange any energy, 
because of they are off resonance with each other.
"""
nc_list = expect(nc, res.states)
nq_list = expect(nq, res.states)
fig, ax = plt.subplots(1, 1, sharex=True, figsize=(12,4))

ax.plot(tlist, nc_list, 'r', linewidth=2, label="cavity")
ax.plot(tlist, nq_list, 'b--', linewidth=2, label="qubit")
ax.set_ylim(0, 7)
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

