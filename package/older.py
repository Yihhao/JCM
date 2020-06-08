from package import *
from qutip import *


def plot_energy_specturm(H_list, N=0, labels=None, show_ylabels=False,
                         figsize=(8, 12), fig=None, ax=None):
    """
    Plot the energy level diagrams for a list of Hamiltonians. Include
    up to N energy levels. For each element in H_list, the energy
    levels diagram for the cummulative Hamiltonian sum(H_list[0:n]) is plotted,
    where n is the index of an element in H_list.

    Parameters
    ----------

        H_list : List of Qobj
            A list of Hamiltonians.

        labels : List of string
            A list of labels for each Hamiltonian

        show_ylabels : Bool (default False)
            Show y labels to the left of energy levels of the initial
            Hamiltonian.

        N : int
            The number of energy levels to plot

        figsize : tuple (int,int)
            The size of the figure (width, height).

        fig : a matplotlib Figure instance
            The Figure canvas in which the plot will be drawn.

        ax : a matplotlib axes instance
            The axes context in which the plot will be drawn.

    Returns
    -------
    fig, ax : tuple
        A tuple of the matplotlib figure and axes instances used to produce
        the figure.

    Raises
    ------

        ValueError
            Input argument is not valid.

    """

    if not isinstance(H_list, list):
        raise ValueError("H_list must be a list of Qobj instances")

    if not fig and not ax:
        fig, ax = plt.subplots(1, 1, figsize=figsize)

    H = H_list[0]
    N = H.shape[0] if N == 0 else min(H.shape[0], N)

    xticks = []
    yticks = []

    x = 0
    evals0 = H.eigenenergies(eigvals=N)
    gnd_energy0 = evals0[0]
    for e_idx, e in enumerate(evals0[:N]):
        energy = e - gnd_energy0
        ax.plot([x, x + 2], array([1, 1]) * energy, 'b', linewidth=2)
        yticks.append(energy)
    xticks.append(x + 1)
    x += 2

    for H1 in H_list[1:]:

        H = H + H1
        evals1 = H.eigenenergies()
        gnd_energy1 = evals1[0]
        for e_idx, e in enumerate(evals1[:N]):
            energy = e - gnd_energy1
            ax.plot([x, x + 1], np.array([evals0[e_idx] - gnd_energy0, energy]), 'k:')
        x += 1

        for e_idx, e in enumerate(evals1[:N]):
            energy = e - gnd_energy1
            ax.plot([x, x + 2], np.array([1, 1]) * energy, 'b', linewidth=2)
        xticks.append(x + 1)
        x += 2

        evals0 = evals1

    ax.set_frame_on(False)

    if show_ylabels:
        yticks = np.unique(np.around(yticks, 1))
        ax.set_yticks(yticks)
    else:
        ax.axes.get_yaxis().set_visible(False)

    if labels:
        ax.get_xaxis().tick_bottom()
        ax.set_xticks(xticks)
        ax.set_xticklabels(labels, fontsize=16)
    else:
        ax.axes.get_xaxis().set_visible(False)

    return fig, ax


def plot_dispective(H_list, plot_line=5, wa=0, wc=0, g=0,
                    labels=None, show_ylabels=True,
                    figsize=(8, 12), fig=None, ax=None):
    """
    Plot the energy level diagrams for a list of Hamiltonians. Include
    up to N energy levels. For each element in H_list, the energy
    levels diagram for the cummulative Hamiltonian sum(H_list[0:n]) is plotted,
    where n is the index of an element in H_list.let ground energy is equal to zero.

    Parameters
    ----------

        H_list : List of Qobj
            A list of Hamiltonians.
        wa : float
            The freq of JCM atom
        wc : float
            The freq of JCM cavity
        g : float
            The freq of JCM couple strength
        labels : List of string
            A list of labels for each Hamiltonian

        show_ylabels : Bool (default False)
            Show y labels to the left of energy levels of the initial
            Hamiltonian.

        plot_line : int
            The number of energy levels to plot

        figsize : tuple (int,int)
            The size of the figure (width, height).

        fig : a matplotlib Figure instance
            The Figure canvas in which the plot will be drawn.

        ax : a matplotlib axes instance
            The axes context in which the plot will be drawn.

    Returns
    -------
    fig, ax : tuple
        A tuple of the matplotlib figure and axes instances used to produce
        the figure.

    Raises
    ------

        ValueError
            Input argument is not valid.

    """
    if labels is None:
        labels = [r'$g=0$', r'$g\neq 0$']
    fig, ax = plot_energy_specturm(H_list, N=plot_line,
                                   labels=labels, show_ylabels=show_ylabels,
                                   figsize=figsize, fig=fig, ax=ax)
    delta = wa - wc

    y_ticks = []
    y_position = []

    if delta == 0:
        for i in np.arange(0, plot_line // 2 + 1):
            if i == 0:
                y_position.append(i * wc)
                y_ticks.append('|%s, g>' % i)
                continue
            dispective = np.sqrt(i) * g
            plt.annotate(s='', xy=(4, i * wc - dispective), xytext=(4, i * wc + dispective),
                         arrowprops=dict(arrowstyle='<->'))
            y = 0.05 * i - 0.2  # plot text position
            plt.text(4.2, i * wc + y, r'$\sqrt{%s}2g$' % i,
                     fontdict={'size': 16, 'color': 'k'})
            y_position.append(i * wc)
            y_ticks.append('|%s, e>, |%s g>' % (i - 1, i))
        plt.yticks(y_position, y_ticks, fontsize=16)
    else:
        i, j, k = 0, 0, 0
        while k < plot_line:
            if i * wc < j * wc + wa:
                y_position.append(i * wc)
                y_ticks.append('|%s, g>' % i)
                i += 1
            else:
                y_position.append(wa + j * wc)
                y_ticks.append('|%s, e>' % j)
                j += 1
            k += 1
        plt.yticks(y_position, y_ticks, fontsize=16)
        plt.annotate(s='', xy=(0.5, wa), xytext=(0.5, wc),
                     arrowprops=dict(arrowstyle='<->'))
        plt.text(0.7, delta / 2 + wc - 0.1, r'$\Delta$',
                 fontdict={'size': 16, 'color': 'k'})

    plt.annotate(s='', xy=(1, wc), xytext=(1, 0),
                 arrowprops=dict(arrowstyle='<->'))
    plt.text(1.2, wc / 2, r'$\omega_c$',
             fontdict={'size': 16, 'color': 'k'})
    plt.text(1.7, wa / 2, r'$\omega_a$',
             fontdict={'size': 16, 'color': 'k'})
    plt.annotate(s='', xy=(1.5, wa), xytext=(1.5, 0),
                 arrowprops=dict(arrowstyle='<->'))
    return fig, ax


def initial_fock_state(N, n=0, wav=(0, 1)):
    """
    default initial state is |0,1> it mean cavity |0> state and atom |1> state.

    Parameters
    ----------
    N:
        Cavity total state
    n:
        occupy of state cavity
    wav:
        a list of two level system. e.g. (0,1) means |1>
    Returns
    -------
    psi:
        a tensor of state |N> * |2>
    """
    # TLS: two level system
    TLS = (wav[0] * basis(2, 0) + wav[1] * basis(2, 1)).unit()  # initial state  0*|0> + 1*|1>
    fock_st = basis(N, n)  # fock state for cavity
    psi = tensor(fock_st, TLS).unit()
    return psi


def initial_coherent_state(N, z=0., wav=(1, 0)):
    """
    default initial state is |0,0> it mean cavity coherent state and atom |0> state.

    Parameters
    ----------
    N:
        Cavity total state, it will influence a and a.dag dimension.
    z:
        coherent state of exp(za+za*)
    wav:
        a list of two level system. e.g. (1,0) means |0>

    Returns
    -------
    psi:
        a tensor of state |N> * |2>
    """

    TLS = wav[0] * basis(2, 0) + wav[1] * basis(2, 1)
    cavity = coherent(N, z)
    psi = tensor(cavity, TLS).unit()
    return psi


def operator(N, number=1):
    """
    define operator.
    return tensor of sm, sz, a and I, the dimension is  |N> * |2>

    sm :
        sigmaminus operator
    sz :
        sigmaz operator
    a  :
        destroy operator
    I  :
        unit operator
    """
    j = number + 1  # mean number of 1/2-spin energy level. e.g. 1 paticle of 1/2-spin have 2 energy level
    sm = tensor(qeye(N), destroy(j))
    sz = sm.dag() * sm - sm * sm.dag()  # [S^+, S^-]= S^z
    a = tensor(destroy(N), qeye(j))
    I = tensor(qeye(N), qeye(j))

    return sm, sz, a, I


def JCM_Hamiltonian(N, wc, wa, g=None, use_rwa=True, eff=False, tuple=False):
    """
    Jaynes–Cummings model Hamiltonian

    Parameters
    ----------
    N:
        Cavity total state
    wc:
        frequence of Cavity
    wa:
        frequence of atom
    g:
        couple strength
    use_rwa:
        use rotation wave?
    eff:
        Hamiltonian is effective ?
    tuple:
        return a Hamiltonian or [H0, Hint]?
    """

    sm, sz, a, I = operator(N, number=1)
    # define Hamiltonian

    if eff is True:
        delta = wa - wc
        xmi = g ** 2 / delta
        H0 = wc * (a.dag() * a + 0.5 * I) + 0.5 * wa * sz
        H1 = xmi * (a.dag() * a + 0.5 * I) * sz
    elif use_rwa is True:
        H0 = wc * a.dag() * a + 0.5 * wa * sz
        H1 = (sm.dag() * a + sm * a.dag())
    else:
        H0 = wc * a.dag() * a + 0.5 * wa * sz
        H1 = (sm.dag() + sm) * (a + a.dag())

    if tuple is True or g is None:
        return H0, H1
    else:
        H = H0 + g * H1
        return H

def initial_state(text, N=None, **kwargs):
    if N is None:
        raise ValueError("N must input data")
    kwargs.setdefault('wav', (0, 1))
    kwargs.setdefault('z', 0)

    wav = kwargs['wav']
    z = kwargs['z']
    if text == 'coherent':
        fig_tilte = 'z'
        psi0 = initial_coherent_state(N, z, wav)
    elif text == 'fock':
        fig_tilte = 'n'
        psi0 = initial_fock_state(N, z, wav)
    else:
        psi0 = 0
        fig_tilte = ''
    return psi0, fig_tilte
