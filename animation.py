import numpy as np
from matplotlib import cm
import matplotlib.pyplot as plt
from qutip.ipynbtools import plot_animation
from qutip import *


def jc_integrate(N, wc, wa, g, kappa, gamma, psi0, use_rwa, tlist):
    # Hamiltonian
    idc = qeye(N)
    ida = qeye(2)

    a = tensor(destroy(N), ida)
    sm = tensor(idc, destroy(2))

    if use_rwa:
        # use the rotating wave approxiation
        H = wc * a.dag() * a + wa * sm.dag() * sm + g * (a.dag() * sm + a * sm.dag())
    else:
        H = wc * a.dag() * a + wa * sm.dag() * sm + g * (a.dag() + a) * (sm + sm.dag())

    # collapse operators
    c_op_list = []

    n_th_a = 0.0  # zero temperature

    rate = kappa * (1 + n_th_a)
    if rate > 0.0:
        c_op_list.append(np.sqrt(rate) * a)

    rate = kappa * n_th_a
    if rate > 0.0:
        c_op_list.append(np.sqrt(rate) * a.dag())

    rate = gamma
    if rate > 0.0:
        c_op_list.append(np.sqrt(rate) * sm)

    # evolve and calculate return state vectors
    result = mesolve(H, psi0, tlist, c_op_list, [])

    return result

def plot_setup(result):
    fig, axes = plt.subplots(figsize=(12, 6))
    return fig, axes


cb = None
def plot_result(result, n, fig=None, axes=None):
    global cb

    if fig is None or axes is None:
        fig, ax = plot_setup(result)

    axes.cla()

    # trace out the atom
    rho_cavity = ptrace(result.states[n], 0)
    fig, axes = plot_wigner(rho_cavity)


    return fig, axes

# parameters
wc = 1.0 * 2 * np.pi   # cavity frequency
wa = 1.0 * 2 * np.pi   # atom frequency
g  = 0.05 * 2 * np.pi  # coupling strength
kappa = 0.05        # cavity dissipation rate
gamma = 0.15        # atom dissipation rate
N = 10              # number of cavity fock states

use_rwa = True

# initial state
psi0 = tensor(basis(N,0), basis(2,1))  # start with an excited atom
tlist = np.linspace(0, 30, 150)
result = jc_integrate(N, wc, wa, g, kappa, gamma, psi0, use_rwa, tlist)
plot_animation(plot_setup, plot_result, result)