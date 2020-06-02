import matplotlib.pyplot as plt
import numpy as np
from numpy import pi
from qutip import *


def Xtest(wr, g, na, eigen_energies, eigen_states):
    trans_energy = eigen_energies[0] - eigen_energies[2]
    element = na.matrix_element(eigen_states[0], eigen_states[2])
    xmi = -1 * g**2 * abs(element) ** 2 * 2.0 * trans_energy / (trans_energy ** 2 - wr ** 2)
    print(trans_energy)
    return xmi


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
    xtest = np.zeros(len(phi_extlist))
    for idx, phi_ext in enumerate(phi_extlist):
        ope = 1.0j * (phi + phi_ext)
        H = 4.0 * E_c * na ** 2.0 + 0.5 * E_l * phi ** 2.0 - 0.5 * E_j * (ope.expm() + (-ope).expm())
        # Eigenstates, eigenvectors
        eigen_energies, eigen_states = H.eigenstates()
        xi[idx] = Xmi_ab(na, eigen_states, eigen_energies, level_num, iState, fState, wr, g)
        xtest[idx] = Xtest(wr, g, na, eigen_energies, eigen_states)
    return xi, xtest


N = 80
level_num = 2
E_l = 0.5
E_c = 2.5
E_j = 10
iState = 1
fState = 0
wr = 4
g = 0.01
phi_ext = -1

phi_extlist = np.linspace(-2 * pi, 0, 100)

fig, ax = plt.subplots(figsize=(8, 8))
Xmi, xtest = charge_dispersive_shift(N, level_num, E_l, E_c, E_j, phi_extlist, iState, fState, wr, g)
ax.plot(phi_extlist / (2 * pi), Xmi, 'r', label='xi')
ax.plot(phi_extlist / (2 * pi), xtest, 'b', label='test')
ax.set_xlabel('Phi external')
ax.set_ylabel('Dispersive shift')
ax.set_title('charge dispersive shift')
plt.legend()
plt.show()
