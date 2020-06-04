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
    result = []
    for index, dt in enumerate(tlist):
        U = expm(-1.0j * H * dt)
        if two_d:
            temp = dot(dagger(U), rho0)
            result.appned(dot(temp, U))
        else:
            result.append(dot(U, rho0))
    return result


def expect_value(operator, state):
    if not isinstance(state, list) and not isinstance(operator, list):
        if operator.shape[0] != operator.shape[0]:
            raise TypeError("input operator dimension must be same [N,N]")
        temp = tensordot(operator, state, axes=[[1], [0]])
        value = tensordot(state, temp, axes=[[0], [0]])
    elif not isinstance(operator, list):
        value = zeros([len(state)], dtype=complex)
        if operator.shape[0] != operator.shape[0]:
            raise TypeError("input operator dimension must be same [N,N]")
        for idx, s in enumerate(state):
            temp = tensordot(operator, s, axes=[[1], [0]])
            value[idx] = tensordot(s, temp, axes=[[0], [0]])
    else:
        value = zeros([len(operator), len(state)], dtype=complex)
        for iop, op in enumerate(operator):
            if op.shape[0] != op.shape[0]:
                raise TypeError("input operator dimension must be same [N,N]")
            for idx, s in enumerate(state):
                temp = tensordot(op, s, axes=[[1], [0]])
                value[iop, idx] = tensordot(s, temp, axes=[[0], [0]])
    return value