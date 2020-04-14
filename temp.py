from qutip import *
import numpy as np
# import matplotlib.pyplot as plt
from  JCM import inital_fock_state
from matplotlib.pyplot import figure, show, legend, xlabel, ylabel, plot


N = 10
a = destroy(N)
x = a.dag() + a

H = a.dag() * a

c_ops = lambda kappa: [np.sqrt(kappa) * a]
times = np.linspace(0, 10.0, 200)

corr1 = correlation_ss(H, times, c_ops(0.5), x, x)
corr2 = correlation_ss(H, times, c_ops(1.0), x, x)
corr3 = correlation_ss(H, times, c_ops(2.0), x, x)

figure()
plot(times, np.real(corr1), times, np.real(corr2), times, np.real(corr3))
legend(['0.5', '1.0', '2.0'])
xlabel(r'Time $t$')
ylabel(r'Correlation $\left<x(t)x(0)\right>$')
show()

psi0 = fock(N, 1)
_, eket = H.eigenstates()
rho = ket2dm(eket[0])
output = essolve(H, rho, times, c_ops(2.0), [x, x])
plot(times, output.expect[0])
plot(times, output.expect[1])
xlabel('Time')
ylabel('Occupation probability')
# title('Vacuum Rabi oscillations')
show()
