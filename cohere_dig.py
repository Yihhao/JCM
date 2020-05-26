from qutip import *
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np

N = 20
z = 1.2
psi0 = coherent(N, z)
fig, ax = plot_fock_distribution(psi0)
ax.set_xticks(np.arange(N+1), minor=False)
plt.savefig('./fig/coherent.png')
plt.show()

plot_wigner(psi0)
plt.xlabel('X', fontsize=16)
plt.ylabel('P', fontsize=16)
plt.savefig('./fig/wigner_coh.png')
plt.show()
