import numpy as np

def convert_Poly2Mat(X, U):
    def poly2ineq(poly):
        return poly.A / np.tile(poly.b[:, np.newaxis], (1, poly.A.shape[1]))
    
    F_tmp = poly2ineq(X) if X.A.size != 0 else np.zeros((0, X.Dim))
    G_tmp = poly2ineq(U) if U.A.size != 0 else np.zeros((0, U.Dim))
    
    F = np.vstack([F_tmp, np.zeros((G_tmp.shape[0], X.Dim))])
    G = np.vstack([np.zeros((F_tmp.shape[0], U.Dim)), G_tmp])
    nc = F.shape[0]
    
    return F, G, nc

