import matplotlib.pyplot as plt
from qutip import *
from numpy import pi


wc = 2.0 * 2 * pi  # cavity frequency
delta = 0          # detuning
wa = wc + delta    # atom frequency
g = 1 * 2 * pi         # coupling strength
N = 3
j = 2

sm = tensor(qeye(N), destroy(j))
sz = sm.dag() * sm - sm * sm.dag()  # [S^+, S^-]= S^z
a = tensor(destroy(N), qeye(j))

# H0 = wc * (a.dag() * a + 0.5 * sz)
# Hint = 0.5 * delta * sz + g * (sm.dag() * a + sm * a.dag())

H0 = wc * a.dag() * a + 0.5 * wa * sz
Hint = g * (sm.dag() * a + sm * a.dag())

H_p = H0 + Hint
eval = H_p.eigenenergies()

print(eval)
plot_energy_levels([H0, Hint], N=4, labels=('|g>', '|g,0>'), figsize=(8,4))
plt.title("Energy level")
plt.savefig('./fig/energy_level')
plt.show()
