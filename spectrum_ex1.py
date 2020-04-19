import numpy as np
from qutip import *
import pylab as plt
from package.JCM import JCM_Hamiltonian, inital_fock_state

wc = 1.0  * 2 * np.pi  # cavity frequency
wa = 1.0  * 2 * np.pi  # atom frequency
g  = 0.01 * 2 * np.pi #0.1 * 2 * np.pi  # coupling strength
kappa = 0.005          # cavity dissipation rate
gamma = 0.05           # atom dissipation rate
N = 15                 # number of cavity fock states
n_th_a = 0.0           # temperature in frequency units
use_rwa = True


# Jaynes-Cummings Hamiltonian
a  = tensor(destroy(N), qeye(2))
sm = tensor(qeye(N), destroy(2))
H = JCM_Hamiltonian(N, wc, wa, g, use_rwa)

# collapse operators
n_th = 0.25
c_ops = [np.sqrt(kappa * (1 + n_th)) * a, np.sqrt(kappa * n_th) * a.dag(), np.sqrt(gamma) * sm]

# calculate the correlation function using the mesolve solver, and then fft to
# obtain the spectrum. Here we need to make sure to evaluate the correlation
# function for a sufficient long time and sufficiently high sampling rate so 
# that the discrete Fourier transform (FFT) captures all the features in the
# resulting spectrum.
tlist = np.linspace(0, 100, 5000)
corr = correlation_ss(H, tlist, c_ops, a.dag(), a)
wlist1, spec1 = spectrum_correlation_fft(tlist, corr)


# calculate the power spectrum using spectrum, which internally uses essolve
# to solve for the dynamics (by default)
# wlist2 = np.linspace(0.25, 1.75, 200) * 2 * np.pi
# spec2 = spectrum(H, wlist2, c_ops, a.dag(), a)

# plot the spectra
fig, ax = plt.subplots(1, 1)
ax.plot(wlist1 / (2 * np.pi), spec1, 'b', lw=2, label='eseries method')
# ax.plot(wlist2 / (2 * np.pi), spec2, 'r--', lw=2, label='me+fft method')
ax.legend()
ax.set_xlabel('Frequency')
ax.set_ylabel('Power spectrum')
ax.set_title('Vacuum Rabi splitting g = %s* 2$\pi$' % (g/(2*np.pi)))
ax.set_xlim(0.25, 1.75)
plt.savefig('./fig/Vacuum_Rabi_splitting_g=%s.png' % (g/(2*np.pi)))
plt.show()

psi0 = inital_fock_state(N, n=0)
output = mesolve(H, psi0, tlist, c_ops, [a.dag() * a, sm.dag() * sm])
fig, ax = plt.subplots(figsize=(8,5))
ax.plot(tlist, output.expect[0], label="Cavity")
ax.plot(tlist, output.expect[1], label="Atom excited state")
ax.legend()
ax.set_xlabel('Time')
ax.set_ylabel('Occupation probability')
ax.set_title('Vacuum Rabi oscillations g = %s * 2$\pi$' % (g/(2*np.pi)))
plt.savefig('./fig/Vacuum_Rabi_oscillations_g=%s.png' % (g/(2*np.pi)))
plt.show()