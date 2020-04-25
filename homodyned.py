import matplotlib.pyplot as plt
import numpy as np
from qutip import *

# define system operators
gamma = 1                           # decay rate
sm_TLS = destroy(2)                 # dipole operator
c_op_TLS = [np.sqrt(gamma)*sm_TLS]  # represents spontaneous emission

# choose range of driving strengths to simulate
Om_list_TLS = gamma*np.logspace(-2, 1, 300)

# calculate steady-state density matricies for the driving strengths
rho_ss_TLS = []
for Om in Om_list_TLS:
    H_TLS = Om * (sm_TLS + sm_TLS.dag())
    rho_ss_TLS.append(steadystate(H_TLS, c_op_TLS))

# decompose the emitted light into the coherent and incoherent
# portions
I_c_TLS = expect(sm_TLS.dag(), rho_ss_TLS)*expect(sm_TLS, rho_ss_TLS)
I_inc_TLS = expect(sm_TLS.dag()*sm_TLS, rho_ss_TLS) - I_c_TLS

plt.semilogx(Om_list_TLS, abs(I_c_TLS),
             label='TLS $I_\mathrm{c}$')
plt.semilogx(Om_list_TLS, abs(I_inc_TLS),
             'r', label='TLS $I_\mathrm{inc}$')
plt.xlabel('Driving strength [$\Gamma$]')
plt.ylabel('Normalized flux [$\Gamma$]')
plt.legend(loc=2)
plt.ylim(0, )
plt.show()