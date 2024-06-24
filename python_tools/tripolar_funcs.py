import numpy as np

def a(y,a_coeff):
    '''
    '''
    a0, a1, a2, a3, a4 = a_coeff
    return a0+y*(a1+y*(a2+y*(a3+y*a4)))

def b(y,b_coeff):
    '''    
    '''
    b0, b1, b2, b3, b4, b5 = b_coeff
    return b0+y*(b1+y*(b2+y*(b3+y*(b4+y*b5))))


def dady(y,a_coeff):
    '''
    '''
    a0, a1, a2, a3, a4 = a_coeff
    return a1+y*(2*a2+y*(3*a3+4*y*a4))

def dbdy(y,b_coeff):
    '''
    '''
    b0, b1, b2, b3, b4, b5 = b_coeff
    return b1+y*(2*b2+y*(3*b3+y*(4*b4+5*y*b5)))


def dydn(n,y,dtheta,fe,le):
    '''
    Equation 36 in Bentsen 20?? 
    '''
    return dtheta*(1-(1-fe)*np.exp(-(y/le)**4))


def dpdn(n,p,dtheta,fe,le,a_coeff,b_coeff):
    '''
    
    '''
    r        = np.zeros(2)
    psi_t    = p[0]
    y_t      = p[1]
    a_t      = a(y_t,a_coeff)
    dady_t   = dady(y_t,a_coeff)
    b_t      = b(y_t,b_coeff)
    dbdy_t   = dbdy(y_t,b_coeff)
    dydn_t   = dydn(n,y_t,dtheta,fe,le)
    cospsi_t = np.cos(psi_t)
    sinpsi_t = np.sin(psi_t)
    #
    r[0] = cospsi_t*sinpsi_t*(a_t*dady_t-b_t*dbdy_t)*dydn_t/((a_t*sinpsi_t)**2+(b_t*cospsi_t)**2)
    r[1] = dydn_t
    #
    return r[0], r[1]
