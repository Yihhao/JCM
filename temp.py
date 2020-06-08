from qutip import *
from package import *
from package.state import coherent_state, fock_state
from package.operator import destroy_op, sigmax_op, sigmaz_op, sigmam_op
from package.plot import plot_fock_number
import matplotlib as mpl
from package.plot_energy import plot_energy_b
from package.JCM import JCM_Hamiltonian, operator
from package.time_evolution import expect_value, H_evolution
from scipy import integrate
from numpy import pi, inf, exp, linspace, meshgrid

def wigner_function(rho, range=5):
    xvec = linspace(-range, range, 100)
    x,y = np.meshgrid(xvec, xvec)
    if rho.shape[0] != rho.shape[1]:
        rho0 = density_matrix(rho)

    z = -np.hypot(x, y)

    plt.contourf(x, y, z, 100)
    # plt.colorbar()

    # plt.axhline(0, color='white')
    # plt.axvline(0, color='white')

    plt.show()
    # invexp = lambda y: exp(-2*1.0j*p*y)
    #
    # print(integrate.quad(invexp, 0, inf))
    return 0


