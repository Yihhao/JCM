from package import *


def H_evolution(H, rho0, tlist):
    """
    H is independent time
    :param H: Hamilton
    :param rho0: initial state
    :param tlist: a list of time
    :return:
    """
    shape = [len(tlist)]
    two_d = False
    for i in rho0.shape:
        shape.append(i)
    if len(rho0.shape) == 1:
        rho0.reshape([len(rho0), 1])
    if rho0.shape[0] == rho0.shape[1]:
        two_d = True
    result = zeros(shape, dtype=complex)
    for index, dt in enumerate(tlist):
        U = expm(-1.0j * H * dt)
        if two_d:
            temp = dot(dagger(U), rho0)
            result[index] = dot(temp, U)
        else:
            result[index] = dot(U, rho0)
    return result