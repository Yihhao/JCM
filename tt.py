import numpy as np
import numpy.linalg as la

A = np.array([[0, 2.12],
              [0.59, 0]])

eval, eket = la.eig(A)