from qutip import *
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np

N = 20
z = 1.2
psi0 = coherent(N, z)
fig, ax = plot_fock_distribution(psi0)
# ax.set_xticks([1.5, 2.5], minor=True)
ax.set_xticks(np.arange(N+1), minor=False)
# plot_wigner_fock_distribution(psi0)
plt.savefig('./fig/coherent.png')
plt.show()


xvec = np.linspace(-5, 5, 200)
fig, ax = plt.subplots(figsize=(6, 6))
W = wigner(psi0, xvec, xvec)
ax.contourf(xvec, xvec, W, 100, norm=mpl.colors.Normalize(-.125, .125), cmap=plt.get_cmap('RdBu'))
plt.xlabel('X', fontsize=16)
plt.ylabel('P', fontsize=16)
plt.title('Phase space \n z = %s' % z)
plt.tight_layout()
plt.savefig('./fig/wigner_coh.png')
plt.show()

parfor()