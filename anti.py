import numpy as np
import matplotlib.pyplot as plt
from numpy import sqrt, array
from package.JCM import JCM_Hamiltonian


# initial parameters
wc = 5.4  # cavity frequency
g = 1e-2  # coupling strength
N = 2  # number of cavity fock states
use_rwa = False  # rwa: rotating wave approximation

energy = []
wa_list = np.linspace(4.5, 6.5, 100)
energies = array([JCM_Hamiltonian(N, wc, wa, g).eigenenergies() for wa in wa_list])

fig, axes = plt.subplots(1,1, figsize=(8,6))

for n in [1, 2]:
    plt.plot(wa_list, energies[:, n])
axes.set_xlim(wa_list[0], wa_list[-1])
axes.set_xlabel(r'$w_a$', fontsize=18)
axes.set_ylabel(r'$E_n$', fontsize=18)
plt.show()
