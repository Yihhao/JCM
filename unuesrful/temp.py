import matplotlib.pyplot as plt
import numpy as np
from numpy import pi
from qutip import *


def Xmi_ab(na, eigen_states, eigen_energies, level_num, iState, fState, wr, g, ):
    shift_iState = 0
    shift_fState = 0

    # iState chi
    for idx in range(level_num):
        if (idx == iState):
            continue
        trans_energy = eigen_energies[idx] - eigen_energies[iState]
        element = na.matrix_element(eigen_states[iState], eigen_states[idx])
        shift_iState += abs(element) ** 2 * 2.0 * trans_energy / (trans_energy ** 2 - wr ** 2)

    # fState chi
    for idx in range(level_num):
        if (idx == fState):
            continue
        trans_energy = eigen_energies[idx] - eigen_energies[fState]
        element = na.matrix_element(eigen_states[fState], eigen_states[idx])
        shift_fState += abs(element) ** 2 * 2.0 * trans_energy / (trans_energy ** 2 - wr ** 2)

    x = (g ** 2) * (shift_iState - shift_fState)
    return x


def charge_dispersive_shift(N, level_num, E_l, E_c, E_j, phi_extlist, iState, fState, wr, g, ):
    a = tensor(destroy(N))
    phi = (a + a.dag()) * (8.0 * E_c / E_l) ** (0.25) / np.sqrt(2.0)
    na = 1.0j * (a.dag() - a) * (E_l / (8 * E_c)) ** (0.25) / np.sqrt(2.0)

    xi = np.zeros(len(phi_extlist))
    for idx, phi_ext in enumerate(phi_extlist):
        ope = 1.0j * (phi + phi_ext)
        H = 4.0 * E_c * na ** 2.0 + 0.5 * E_l * phi ** 2.0 - 0.5 * E_j * (ope.expm() + (-ope).expm())
        # Eigenstates, eigenvectors
        eigen_energies, eigen_states = H.eigenstates()
        xi[idx] = Xmi_ab(na, eigen_states, eigen_energies, level_num, iState, fState, wr, g)
    return xi


# charge matrix element
def charge_matrix_element(N, El, Ec, Ej, iState, fState):
    a = tensor(destroy(N))
    phi = (a + a.dag()) * (8.0 * Ec / El) ** (0.25) / np.sqrt(2.0)
    na = 1.0j * (a.dag() - a) * (El / (8 * Ec)) ** (0.25) / np.sqrt(2.0)
    element = np.zeros((len(phi_extlist), N))

    for idx, phi_ext in enumerate(phi_extlist):
        ope = 1.0j * (phi + phi_ext)
        H = 4.0 * Ec * na ** 2.0 + 0.5 * El * phi ** 2.0 - 0.5 * Ej * (ope.expm() + (-ope).expm())
        eigen_energies, eigen_states = H.eigenstates()
        element[idx] = (abs(na.matrix_element(eigen_states[iState], eigen_states[fState])))

    return element


def specturm(N, E_l, E_c, E_j, phi_extlist):
    a = tensor(destroy(N))
    phi = (a + a.dag()) * (8.0 * E_c / E_l) ** (0.25) / np.sqrt(2.0)
    na = 1.0j * (a.dag() - a) * (E_l / (8 * E_c)) ** (0.25) / np.sqrt(2.0)

    evals_mat = np.zeros((len(phi_extlist), N))
    for index, phi_ext in enumerate(phi_extlist):
        ope = 1.0j * (phi + phi_ext)
        # evaluate the Hamiltonian
        H = 4.0 * E_c * na ** 2.0 + 0.5 * E_l * phi ** 2.0 - 0.5 * E_j * (ope.expm() + (-ope).expm())
        # find the energy eigenvalues of the composite system
        evals, ekets = H.eigenstates()
        evals_mat[index, :] = np.real(evals)

    return evals_mat


def plot_energy_b(delta, evals_mat, plot_line,
                  fig=None, ax=None, figsize=(8, 12)):
    color = ['b', 'r', 'k']
    if not fig and not ax:
        fig, ax = plt.subplots(1, 1, figsize=figsize)
    for n in np.arange(plot_line):
        if n == 0:
            c = color[2]
        else:
            c = color[n % 2]
        ax.plot(delta/(2*pi), (evals_mat[:, n] - evals_mat[:, 0]), c)

    return fig, ax


N = 80
level_num = 2
E_l = 0.5
E_c = 2.5
E_j = 10
iState = 0
fState = 1
wr = 4
g = 0.01
phi_ext = -1

phi_extlist = np.linspace(-2 * pi, 0, 100)
# evals_mat = specturm(N, E_l, E_c, E_j, phi_extlist)
#
# plot_energy_b(phi_extlist, evals_mat, plot_line=5)
# plt.show()

# element = charge_matrix_element(N, E_l, E_c, E_j, iState, fState=2)
# xtest = -1 * g * element ** 2
# plt.plot(phi_extlist, element)
# plt.show()

fig, ax = plt.subplots(figsize=(8, 8))
Xmi = charge_dispersive_shift(N, level_num, E_l, E_c, E_j, phi_extlist, iState, fState, wr, g)
ax.plot(phi_extlist / (2 * pi), Xmi, 'r', label='xi')
ax.set_xlabel('Phi external')
ax.set_ylabel('Dispersive shift')
ax.set_title('charge dispersive shift')
plt.legend()
plt.show()
