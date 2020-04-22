import matplotlib.pyplot as plt
import numpy as np
from numpy import pi
from qutip import *
from qutip.piqs import *
from numpy import sqrt
from mpl_toolkits.mplot3d import Axes3D


N = 10
psi = coherent(N, 1j)
n = np.arange(N)
psi0 = psi.full()


# rho0 = psi * psi.dag()
# rho = abs(rho0.full())**2
#
# xv = 7.5
# xvec = np.linspace(-7.5, 7.5, 100)
# M = np.arange(N)
# X, Y = np.meshgrid(xvec, xvec)
# A = sqrt(2) * (X + 1j * Y)
# B = abs(A) ** 2
#
# plt.imshow(B, origin=' lower ', extent=[-7.5, 7.5, -7.5, 7.5])
# plt.show()
#
# w0 = abs((2*rho.data[0,-1])*np.ones_like(A))
# plt.imshow(w0, origin=' lower ', extent=[-7.5, 7.5, -7.5, 7.5])
# plt.show()

# psi = fock(N, 5)
psi = coherent(N, 2j)
plot_wigner(psi)
plt.show()


fig = plt.figure()
axis = fig.gca(projection='3d')
xvec = np.linspace(-7.5, 7.5, 100)
X, Y = np.meshgrid(xvec, xvec)
W = wigner(psi, xvec, xvec)

# 畫圖
surface = axis.plot_surface(X, Y, W, rstride=1, cstride=1, cmap='coolwarm_r')
fig.colorbar(surface, shrink=1.0, aspect=20)
plt.tight_layout()
plt.show()
