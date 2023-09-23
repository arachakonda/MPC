import control as ct

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


