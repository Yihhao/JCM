from package.plot import plot_fock_number
from package.state import fock_state, coherent_state
from package import *

N = 20
z = 0
f = arange(N)
# psi = fock_state(20, 0.2)
psi = coherent_state(N, z)
fig, ax = plt.subplots(figsize=(12, 8))
plot_fock_number(psi, fig=fig, ax=ax)
ax.set_xticks(f)
ax.set_xlim(-.5, 9)
ax.set_xticklabels(f)
plt.tight_layout()
plt.savefig('./fig/pure.png')
plt.show()
