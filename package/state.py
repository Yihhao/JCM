from package.JCM import initial_coherent_state, initial_fock_state


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
    return psi0, fig_tilte
