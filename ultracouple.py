# setup the matplotlib graphics library and configure it to show figures inline in the notebook
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
# make qutip available in the rest of the notebook
from qutip import *
from package.JCM import JCM_Hamiltonian
from numpy import pi, sqrt ,real, zeros, array, linspace, where, arange

wc = 1.0 * 2 * pi  # cavity frequency
wa = 1.0 * 2 * pi  # atom frequency

N = 15  # number of cavity fock states
use_rwa = False

# operators
a = tensor(destroy(N), qeye(2))
sm = tensor(qeye(N), destroy(2))

na = sm.dag() * sm  # atom
nc = a.dag() * a  # cavity

# decoupled Hamiltonian
g_vec = linspace(0, 2.0, 101) * 2 * pi  # coupling strength vector

psi_list = []
for index, g in enumerate(g_vec):
    H = JCM_Hamiltonian(N, wc, wa, g, use_rwa=use_rwa)
    # find the groundstate and its energy
    gnd_energy, gnd_state = H.groundstate()

    # store the ground state
    psi_list.append(gnd_state)

na_expt = expect(na, psi_list)  # qubit  occupation probability
nc_expt = expect(nc, psi_list)  # cavity occupation probability

fig, axes = plt.subplots(1, 1, sharex=True, figsize=(8, 4))

axes.plot(g_vec / (2 * pi), nc_expt, 'r', linewidth=2, label="cavity")
axes.plot(g_vec / (2 * pi), na_expt, 'b', linewidth=2, label="atom")
axes.set_ylabel("Occupation probability", fontsize=16)
axes.set_xlabel("coupling strenght", fontsize=16)
axes.legend(loc=0)

fig.tight_layout()
plt.show()

g_idx = where([g_vec == 2 * pi * g for g in [0.0, 0.5, 1.0, 1.5, 2.0]])[1]
psi_sublist = array(psi_list, dtype=object)[g_idx]

xvec = linspace(-5, 5, 200)

fig_grid = (2, len(psi_sublist) * 2)
fig = plt.figure(figsize=(3 * len(psi_sublist), 6))

for idx, psi in enumerate(psi_sublist):
    rho_cavity = ptrace(psi, 0)
    W = wigner(rho_cavity, xvec, xvec)
    ax = plt.subplot2grid(fig_grid, (0, 2 * idx), colspan=2)
    ax.contourf(xvec, xvec, W, 100, norm=mpl.colors.Normalize(-.125, .125), cmap=plt.get_cmap('RdBu'))
    ax.set_title(r"$g = %.1f$" % (g_vec[g_idx][idx] / (2 * pi)), fontsize=16)

# plot the cavity occupation probability in the ground state
ax = plt.subplot2grid(fig_grid, (1, 1), colspan=(fig_grid[1] - 2))
ax.plot(g_vec / (2 * pi), nc_expt, label="Cavity")
ax.plot(g_vec / (2 * pi), na_expt, label="Atom excited state")
ax.legend(loc=0)
ax.set_xlabel('coupling strength')
ax.set_ylabel('Occupation probability')
plt.show()

entropy_cavity = zeros(shape(g_vec))
entropy_atom = zeros(shape(g_vec))

for idx, psi in enumerate(psi_list):
    rho_cavity = ptrace(psi, 0)
    entropy_cavity[idx] = entropy_vn(rho_cavity, 2)

    rho_atom = ptrace(psi, 1)
    entropy_atom[idx] = entropy_vn(rho_atom, 2)
fig, axes = plt.subplots(1, 1, figsize=(12,6))
axes.plot(g_vec/(2*pi), entropy_cavity, 'b', label="cavity", linewidth=2)
axes.plot(g_vec/(2*pi), entropy_atom, 'r--', label="atom", linewidth=2)
axes.set_ylim(0,1)
axes.set_ylabel("entropy", fontsize=16)
axes.set_xlabel("coupling strength", fontsize=16)
axes.legend(loc=0)
plt.show()

H = JCM_Hamiltonian(N, wc, wa, 1.0 * 2 * pi, use_rwa=use_rwa)

psi0 = tensor(basis(N,1), basis(2,0))
tlist = linspace(0, 20, 1000)
output = mesolve(H, psi0, tlist, [], [a.dag() * a, sm.dag() * sm])
fig, axes = plt.subplots(1, 1, sharex=True, figsize=(8,4))

axes.plot(tlist, real(output.expect[0]), 'r', linewidth=2, label="cavity")
axes.plot(tlist, real(output.expect[1]), 'b', linewidth=2, label="atom")
axes.legend(loc=0)

fig.tight_layout()
plt.show()

tlist = linspace(0, 0.35, 8)
output = mesolve(H, psi0, tlist, [], [])
rho_ss_sublist = output.states  # [::4]

xvec = linspace(-5, 5, 200)

fig, axes = plt.subplots(2, len(rho_ss_sublist), figsize=(2 * len(rho_ss_sublist), 4))

for idx, rho_ss in enumerate(rho_ss_sublist):
    # trace out the cavity density matrix
    rho_ss_cavity = ptrace(rho_ss, 0)

    # calculate its wigner function
    W = wigner(rho_ss_cavity, xvec, xvec)

    # plot its wigner function
    axes[0, idx].contourf(xvec, xvec, W, 100, norm=mpl.colors.Normalize(-.25, .25),
                          cmap=plt.get_cmap('RdBu'))

    # plot its fock-state distribution
    axes[1, idx].bar(arange(0, N), real(rho_ss_cavity.diag()), color="blue", alpha=0.6)
    axes[1, idx].set_ylim(0, 1)
    axes[1, idx].set_xlim(0, N)
plt.show()

kappa = 0.25
tlist = linspace(0, 20, 1000)
output = mesolve(H, psi0, tlist, [sqrt(kappa) * a], [a.dag() * a, sm.dag() * sm])
fig, axes = plt.subplots(1, 1, sharex=True, figsize=(8,4))
axes.plot(tlist, output.expect[0], 'r', linewidth=2, label="cavity")
axes.plot(tlist, output.expect[1], 'b', linewidth=2, label="atom")
axes.legend(loc=0)
plt.show()

tlist = linspace(0, 10, 8)
output = mesolve(H, psi0, tlist, [sqrt(kappa) * a], [])
xvec = linspace(-5, 5, 200)

fig, axes = plt.subplots(2, len(output.states), figsize=(2 * len(output.states), 4))

for idx, rho_ss in enumerate(output.states):
    # trace out the cavity density matrix
    rho_ss_cavity = ptrace(rho_ss, 0)

    # calculate its wigner function
    W = wigner(rho_ss_cavity, xvec, xvec)

    # plot its wigner function
    axes[0, idx].contourf(xvec, xvec, W, 100,
                          norm=mpl.colors.Normalize(-.25, .25), cmap=plt.get_cmap('RdBu'))

    # plot its fock-state distribution
    axes[1, idx].bar(arange(0, N), real(rho_ss_cavity.diag()), color="blue", alpha=0.6)
    axes[1, idx].set_ylim(0, 1)
    axes[1, idx].set_xlim(0, N)
plt.show()

tlist = linspace(0, 30, 50)

psi0 = H.groundstate()[1]

output = mesolve(H, psi0, tlist, [sqrt(kappa) * a], [])
entropy_tot = zeros(shape(tlist))
entropy_cavity = zeros(shape(tlist))
entropy_atom = zeros(shape(tlist))

for idx, rho in enumerate(output.states):
    entropy_tot[idx] = entropy_vn(rho, 2)

    rho_cavity = ptrace(rho, 0)
    entropy_cavity[idx] = entropy_vn(rho_cavity, 2)

    rho_atom = ptrace(rho, 1)
    entropy_atom[idx] = entropy_vn(rho_atom, 2)
fig, axes = plt.subplots(1, 1, figsize=(12, 6))
axes.plot(tlist, entropy_tot, 'k', label="total", linewidth=2)
axes.plot(tlist, entropy_cavity, 'b', label="cavity", linewidth=2)
axes.plot(tlist, entropy_atom, 'r--', label="atom", linewidth=2)
axes.set_ylabel("entropy", fontsize=16)
axes.set_xlabel("coupling strength", fontsize=16)
axes.set_ylim(0, 1.5)
axes.legend(loc=0)
plt.show()
