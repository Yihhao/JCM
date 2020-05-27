import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
from qutip import *
from qutip.ipynbtools import version_table
from package.JCM import JCM_Hamiltonian
from numpy import pi, cos, sin


# def H1_coeff(t, args):
#     # tp_G = 2.000  # Gaussian pulse parameter
#     # Om_G = 1.0  # driving strength
#     # t_offset_G = 5
#     # pulse_shape_G = Om_G / 2 * np.exp(-(t - t_offset_G) ** 2 /
#     #                                   (2 * tp_G ** 2))
#     w = 5 * pi / 6
#     return t + cos(2 * w * t)

def H1_coeff(t, args):
    w = 5 * pi / 6
    return -2 * cos(w * t)

# def H2_coeff(t, args):
#     w = 5 * pi / 6
#     return -2 * sin(w * t)

# delta = 0.5 * 2 * np.pi
# v = 2.0 * 2 * np.pi  # sweep rate
# H0 = delta / 2.0 * sigmax()
# H1 = v / 2.0 * sigmaz()
# H = [H0, [H1, H1_coeff]]
# expt_list = mesolve(H, psi0, tlist, [], expt_ops).expect


# parameters
N = 3
wc = 1.0 * 2 * pi
wa = 1.0 * 2 * pi
g = 0.05 * 2 * pi
use_rwa = False


sm = destroy(2)
sx = sigmax()
sy = sigmay()
sz = sigmaz()

psi0 = basis(2, 0)

H0 = wa / 2 * sz
# H1 = g / 4 * sx
# H = [[H1, H1_coeff]]
H1 = g / 2 * (sm + sm.dag())
# H2 = sm.dag()
H = [H0, [H1, H1_coeff]]

tlist = np.linspace(0, 10.0, 150)
res = mesolve(H, psi0, tlist, [], [])

expt_ops = [sm.dag() * sm, sx, sy, sz]

# result = []
# for idx, _ in enumerate(tlist):
#     result.append(ptrace(res.states[idx], 1))

expt_list = expect(expt_ops, res.states)

fig = plt.figure(figsize=(12, 12))
ax = fig.gca(projection='3d')
b = Bloch(fig=fig, axes=ax)
## normalize colors to times in tlist ##
# nrm = mpl.colors.Normalize(-2, 10)
# colors = cm.cool(nrm(tlist))

## add data points from expectation values ##
b.add_points([expt_list[2],
              expt_list[1],
              -expt_list[3]], 'l')

## customize sphere properties ##
# b.point_color = list(colors)
# b.point_marker = ['o']
# b.point_size = [20]

# b.zlpos = [1.1, -1.2]
b.show()

# left, bottom, width, height = [0.98, 0.05, 0.05, 0.9]
# ax2 = b.fig.add_axes([left, bottom, width, height])
#
# mpl.colorbar.ColorbarBase(ax2, cmap=cm.cool,
#                           norm=nrm,
#                           orientation='vertical')
plt.tight_layout()
plt.show()

# plt.plot(tlist, expt_list[1], 'r', label='x')
# plt.plot(tlist, expt_list[2], 'b.', label='y')
# plt.plot(tlist, expt_list[3], 'k:', label='z')
# plt.legend()
# plt.show()
