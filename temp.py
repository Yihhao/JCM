# kappa = 0.05  # cavity dissipation rate
# gamma = 0.05  # atom dissipation rate
# n_th = 1  # temperature in frequency unit

# sm, sz, a = operator(N)
# c_ops = [np.sqrt(kappa * (1 + n_th)) * a, np.sqrt(kappa * n_th) * a.dag(), np.sqrt(gamma) * sm]
# e_ops = [sz, sm.dag() * sm, sm * sm.dag(), sm.dag() * sm + sm * sm.dag()]
# output = mesolve(H, psi, t, c_ops, e_ops)
#
# n_a = output.expect[0]
# n_b = output.expect[1]
# n_c = output.expect[2]
# n_d = output.expect[3]
#
# fig, axes = plt.subplots(1, 1, figsize=(10, 6))
# axes.plot(g * t, n_a, label="pop inversion")
# axes.legend(loc=0)
# axes.set_xlabel('Time')
# axes.set_ylabel('expectation value')
# axes.set_title(text_tilte)
# axes.text(0, 0, 'delta=%.2f' % delta)
# axes.axis([0, 100, -1, 1])
# plt.savefig(filename+'_decay_z')
# plt.show()
#
# fig, axes = plt.subplots(1, 1, figsize=(10, 6))
# axes.plot(g * t, n_b, label="excited")
# axes.plot(g * t, n_c, label="ground")
# axes.plot(g * t, n_d, label="norm")
# axes.legend(loc=0)
# axes.set_xlabel('Time')
# axes.set_ylabel('expectation value')
# axes.set_title(text_tilte)
# axes.axis([0, 100, 0, 1.5])
# plt.savefig(filename+'_expect')
# plt.show()