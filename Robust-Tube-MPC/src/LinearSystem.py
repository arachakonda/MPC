import control as ct
import numpy as np
from utils import *
class LinearSystem:
    def __init__(self, A, B, Q, R):
        self.A = A
        self.B = B
        self.Q = Q
        self.R = R
        self.nx = A.shape[0]
        self.nu = B.shape[1]
        #calculate the LQR controller
        k,p,e = ct.lqr(A, B, Q, R)
        self.K = -k
        self.P = p
        self.Ak = self.A + self.B @ self.K


    def propagate(self, x, u):
        x_new = self.A @ x + self.B @ u
        return x_new

    def compute_MPIset(self, Xc, Uc):
        F, G, nc = convert_Poly2Mat(Xc, Uc)
        Fpi = lambda i: np.dot((F + np.dot(G, self.K)), np.linalg.matrix_power(self.Ak, i))
        
        def Xpi(i):
            return pc.Polyhedron(Fpi(i), np.ones((Fpi(i).shape[0], 1)))
        
        Xmpi = Xpi(0)
        i = 0
        while True:
            i += 1
            Xmpi_tmp = Xmpi.intersect(Xpi(i))
            if Xmpi_tmp == Xmpi:
                break
            else:
                Xmpi = Xmpi_tmp
        return Xmpi
    
