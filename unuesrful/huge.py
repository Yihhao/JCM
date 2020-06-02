def Xtest(na, eigen_energies, eigen_states):
    trans_energy = eigen_energies[2] - eigen_energies[0]
    element = na.matrix_element(eigen_states[0], eigen_states[2])
    xmi = -1 * g**2 * abs(element) ** 2 * 2.0 * trans_energy / (trans_energy ** 2 - wr ** 2)
    return xmi