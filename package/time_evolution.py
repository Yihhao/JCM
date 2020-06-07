from package import *


def H_evolution(H, rho0, tlist):
    """
    H is time-independent.
    :param H: Hamilton
    :param rho0: initial state
    :param tlist: a list of time
    :return: density state of time evolution
    """

    shape = [len(tlist)]
    for i in rho0.shape:
        shape.append(i)
    if len(rho0.shape) == 1:
        rho0.reshape([len(rho0), 1])
    if rho0.shape[0] != rho0.shape[1]:
        rho0 = density_matrix(rho0)
    result = []
    dt_list = np.diff(tlist)
    result.append(rho0)

    for index, dt in enumerate(dt_list):
        U = expm(-1.0j * H * dt)
        temp = dot(dagger(U), rho0)
        temp = dot(temp, U)
        rho0 = temp
        result.append(temp)
    return result


def expect_value(operator, rho):
    if not isinstance(rho, list) and not isinstance(operator, list):
        if operator.shape[0] != operator.shape[0]:
            raise TypeError("input operator dimension must be same [N,N]")
        value = trace(tensordot(rho, operator, axes=[[1], [0]]))
    elif not isinstance(operator, list):
        value = zeros([len(rho)], dtype=complex)
        if operator.shape[0] != operator.shape[0]:
            raise TypeError("input operator dimension must be same [N,N]")
        for idx, r in enumerate(rho):
            value[idx] = trace(tensordot(r, operator, axes=[[1], [0]]))
    else:
        value = zeros([len(operator), len(rho)], dtype=complex)
        for iop, op in enumerate(operator):
            if op.shape[0] != op.shape[0]:
                raise TypeError("input operator dimension must be same [N,N]")
            for idx, r in enumerate(rho):
                value[iop, idx] = trace(tensordot(r, op, axes=[[1], [0]]))
    return value