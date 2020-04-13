import matplotlib.pyplot as plt
import numpy as np
from numpy import pi
from qutip import *
from qutip.piqs import *


N = 6
ntls = N
nds = num_dicke_states(ntls)
[jx, jy, jz] = jspin(N)
jp = jspin(N, "+")
jm = jp.dag()
w0 = 1
gE = 0.1
gD = 0.01
gP = 0.1
gCP = 0.1
gCE = 0.1
gCD = 0.1
h = w0 * jz
system = Dicke(N, emission = 1, pumping = 2)