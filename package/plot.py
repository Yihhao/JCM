import numpy as np
import matplotlib.pyplot as plt
from qutip import *
from package.JCM import operator

def plot_fock(N, psi):
    rho = []
    n = np.arange(N)
    psi0 = psi.full()
    for i in n:
        prob = abs(psi0[i, 0])**2
        rho.append(prob)
    plt.bar(n, rho)
    plt.ylabel('Occupation probability')
    plt.xlabel('Fock number')
    plt.xlim(0, )
    plt.ylim(0., 1.0)
    plt.show()


def expect_fig(res=None, N=None, taulist=None, e_ops=None, **kwargs):
    if not isinstance(N, int):
        raise ValueError("N must be integral")
    if res is None or taulist is None:
        raise ValueError("res or taulist must input data")

    if e_ops is None:
        #         e_ops = {'Cavity_n': nc, 'Atom_excited': excited,
        #                  'Cavity_position': xc, 'Atom_ground': ground}
        eval_dict = expect_operator(res, N, e_ops)

    kwargs.setdefault('column', 1)
    kwargs.setdefault('row', 1)
    kwargs.setdefault('sharex', True)
    kwargs.setdefault('figsize', (12, 6))
    kwargs.setdefault('plotnum', [2, 1, 1])

    fig, axes = plt.subplots(kwargs['row'], kwargs['column'],
                             sharex=kwargs['sharex'], figsize=kwargs['figsize'])

    key = list(eval_dict.keys())
    color = ['r', 'b']
    if kwargs['row'] == 1 and kwargs['column'] == 1:
        axes.plot(taulist, eval_dict[key[0]], 'b', label=key[0])
        axes.legend(loc=1)
    else:
        i = 0
        while kwargs['row'] > i:
            j = 1
            label, ylabel = key[i].split('_')
            axes[i].plot(taulist, eval_dict[key[i]], color[i % 2], label=label)
            while kwargs['plotnum'][i] > j:
                label, ylabel = key[i+j].split('_')
                axes[i].plot(taulist, eval_dict[key[i+j]], color[(i+j) % 2], label=label)
                j += 1
                i += 1
            axes[i].set_ylabel(ylabel, fontsize=16)
            axes[i].legend(loc=1)
            i += 1
        axes[-1].set_xlabel('Time', fontsize=16)
    # fig.suptitle(text_tilte + detuning, fontsize=16)
    # axes[0].plot(taulist, nc_list, 'b', label="Cavity")
    # axes[0].plot(taulist, na_list, 'r', label="Atom")
    # axes[0].set_ylabel("n", fontsize=16)
    # axes[0].legend(loc=1)
    #
    # axes[1].plot(taulist, na_list, 'r', label="Atom")
    # axes[1].set_ylabel("n", fontsize=16)
    # axes[1].legend(loc=1)
    #
    # axes[2].plot(taulist, xc_list, 'b', label="Cavity")
    # axes[2].set_ylabel("position", fontsize=16)
    # axes[2].set_xlabel('Time', fontsize=16)
    # axes[2].legend(loc=1)
    #
    plt.xlim(0, taulist[-1])
    plt.tight_layout()
    plt.subplots_adjust(top=0.9)
    # # plt.savefig(filename, dpi=720)
    # plt.show()
    return 0


def expect_operator(res=None, N=None, e_ops=None):
    sm, sz, a, I = operator(N)

    nc = a.dag() * a
    xc = a.dag() + a
    excited = sm.dag() * sm
    ground = sm * sm.dag()

    if e_ops is None:
        e_ops = {'Cavity_n': nc, 'Atom_excited': excited,
                 'Cavity_position': xc, 'Atom_ground': ground}

    eval = expect(list(e_ops.values()), res.states)
    eval_dict = e_ops.copy()

    d = 0
    for i in e_ops:
        eval_dict.update({i: eval[d]})
        d += 1
    return eval_dict


if __name__ == '__main__':
    n = 5
    H = tensor(qeye(n), qeye(2))
    psi = tensor(fock(n, 0), fock(2, 0))
    tlist = np.linspace(0, 5, 100)
    out = mesolve(H, psi, tlist, [], [])
    # e = expect_operator(out, n)
    expect_fig(out, n, tlist, row=3)
    plt.show()
