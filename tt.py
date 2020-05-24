# % matplotlib
# inline
import matplotlib.pyplot as plt
import numpy as np
from numpy import pi
from qutip import *


def test(N, level_num, E_l, E_c, E_j, phi_extlist, iState, fState, wr, g, ):
    a = tensor(destroy(N))
    phi = (a + a.dag()) * (8.0 * E_c / E_l) ** (0.25) / np.sqrt(2.0)
    na = 1.0j * (a.dag() - a) * (E_l / (8 * E_c)) ** (0.25) / np.sqrt(2.0)

    xi = np.zeros(len(phi_extlist))
    x_test = np.zeros(len(phi_extlist))
    for idx, phi_ext in enumerate(phi_extlist):
        ope = 1.0j * (phi + phi_ext)
        H = 4.0 * E_c * na ** 2.0 + 0.5 * E_l * phi ** 2.0 - 0.5 * E_j * (ope.expm() + (-ope).expm())
        # Eigenstates, eigenvectors
        eigen_energies, eigen_states = H.eigenstates()
        x_test[idx] = compute(na, eigen_states, eigen_energies, level_num, iState, fState, wr, g, )
    return x_test


def compute(na, eigen_states, eigen_energies, level_num, iState, fState, wr, g, ):
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
    # print(x)
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


        shift_iState = 0
        shift_fState = 0
        # iState chi
        for idx in range(level_num):
            if (idx == iState):
                continue
            trans_energy = eigen_energies[idx] - eigen_energies[iState]
            element = na.matrix_element(eigen_states[iState], eigen_states[idx])
            shift_iState += abs(element) ** 2 * 2.0 * trans_energy / (trans_energy ** 2 - wr ** 2)
        xtest[idx] = compute(na, eigen_states, eigen_energies, level_num, iState, fState, wr, g, )
        # fState chi
        for idx in range(level_num):
            if (idx == fState):
                continue
            trans_energy = eigen_energies[idx] - eigen_energies[fState]
            element = na.matrix_element(eigen_states[fState], eigen_states[idx])
            shift_fState += abs(element) ** 2 * 2.0 * trans_energy / (trans_energy ** 2 - wr ** 2)

        xi[idx] = g ** 2 * (shift_iState - shift_fState)

    return xi, xtest

N = 80
level_num = 4
E_l = 10
E_c = 2.5
E_j = 0.5
iState = 0
fState = 1
wr = 1
# g = 0.001
g = 1000
phi_ext = -1

phi_extlist = np.linspace(-2 * pi, 0, 100)
fig, ax = plt.subplots(figsize=(8, 8))
# Xmi, x = charge_dispersive_shift(N, level_num, E_l, E_c, E_j, phi_extlist, iState, fState, wr, g)
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

    xtest[idx] = compute(na, eigen_states, eigen_energies, level_num, iState, fState, wr, g, )

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

    xi[idx] = g ** 2 * (shift_iState - shift_fState)
    x = g ** 2 * (shift_iState - shift_fState)

ax.plot(phi_extlist / (2 * pi), xtest, 'b', label='test')
ax.plot(phi_extlist / (2 * pi), xi, 'r', label='xi')
ax.scatter(phi_extlist[-1], x, c='k')
ax.set_xlabel('Phi external')
ax.set_ylabel('Dispersive shift')
ax.set_title('charge dispersive shift')
plt.legend()
plt.show()
