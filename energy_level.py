import matplotlib.pyplot as plt
from qutip import *
from numpy import pi


wc = 2.0 * 2 * pi  # cavity frequency
delta = 0          # detuning
wa = wc + delta    # atom frequency
g = 2 * pi         # coupling strength
N = 5
j = 2

sm = tensor(qeye(N), destroy(j))
sz = sm.dag() * sm - sm * sm.dag()  # [S^+, S^-]= S^z
a = tensor(destroy(N), qeye(j))

# H0 = wc * (a.dag() * a + 0.5 * sz)
# Hint = 0.5 * delta * sz + g * (sm.dag() * a + sm * a.dag())

H0 = wc * a.dag() * a + 0.5 * wa * sz
Hint = g * (sm.dag() * a + sm * a.dag())

plot_energy_levels([H0, Hint], figsize=(8,4))
plt.show()
