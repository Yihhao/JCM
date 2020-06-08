from package import *


def plot_fock_number(rho0, fig=None, ax=None, figsize=(8, 6)):
    if fig is None and ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    if len(rho0.shape) == 1:
        rho0.reshape([len(rho0), 1])
    if rho0.shape[0] != rho0.shape[1]:
        rho0 = density_matrix(rho0)
    N = rho0.shape[0]
    for i in arange(N):
        ax.bar(i, rho0.real[i, i],
               color="green", alpha=0.6, width=0.8)
    ax.set_ylim(0, 1)
    ax.set_xlim(-.5, N)
    ax.set_xlabel('Fock number', fontsize=12)
    ax.set_ylabel('Occupation probability', fontsize=12)
    return fig, ax