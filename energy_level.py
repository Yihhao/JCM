import matplotlib.pyplot as plt
from package.JCM import JCM_Hamiltonian
from qutip import *
from numpy import pi


wc = 2.0 * 2 * pi  # cavity frequency
delta = 0          # detuning
wa = wc + delta    # atom frequency
g = 0.1 * 2 * pi         # coupling strength
N = 3
j = 2
usr_rwa = False

sm = tensor(qeye(N), destroy(j))
sz = sm.dag() * sm - sm * sm.dag()  # [S^+, S^-]= S^z
a = tensor(destroy(N), qeye(j))


H0, H1 = JCM_Hamiltonian(N, wc, wa, g, usr_rwa)

H = H0 + H1
eval = H.eigenenergies()
print(eval)

plot_energy_levels([H0, H1], N=5, labels=('|g>', '|g,0>'), figsize=(8,4))
plt.title("Energy level")
plt.savefig('./fig/energy_level')
plt.show()
