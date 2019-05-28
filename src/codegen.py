
import numpy as np
import sympy as sp
from sympy import *
from sympy.physics.quantum.spin import Rotation

def carttospher(x,y,z):
    # INPUTS: (x,y,z): This represents a position vector which is given using Cartesian coordinates.
    # RETURNS: The corresponding polar coordinates - radial, azimuthal, polar.

    r = np.sqrt(x**2+y**2+z**2)
    beta = np.arccos(z/r)
    if x != 0:
        alpha = np.arctan2(y,x)
    else:
        if y > 0:
            alpha = np.pi/2
        elif y < 0:
            alpha = -np.pi/2
        else:
            alpha = 0
    return r, alpha, beta

def Gsym(l_1,l_2,m_1,m_2,M):
    alpha, beta = symbols('alpha, beta')
    # s = signmatrix(l_1, l_2)[m_1+l_1, m_2+l_2]
    return str(G2(l_1, l_2, m_1, m_2, M, alpha, beta))


def G2(l_1,l_2,m_1,m_2,M,alpha,beta):

   # INPUTS: l_1: Angular quantum number of atom 1.
   #         l_2: Angular quantum number of atom 2.
   #         m_1: Magnetic quantum number of atom 1.
   #         m_2: Magnetic quantum number of atom 2.
   #         M: A float which represents the type of symmetry bond we are considering, e.g. M = 0
   #            for sigma, M = 1 for pi, M = 2 for delta, etc.
   #         alpha: Azithmuthal coordinate.
   #         beta: Polar coordinate.
   # RETURNS: This function returns the relevant coefficient that should be are used in writing the
   #          Slater-Koster transformations.

   def tau(m):
       if m >= 0:
           return 1.0
       else:
           return 0.0

   def A(m,alpha):
       if m == 0:
           return np.sqrt(2)/2
       else:
           return ((-1)**m)*(tau(m)*sp.cos(np.absolute(m)*alpha)+tau(-m)*sp.sin(np.absolute(m)*alpha))

   def B(m,alpha):
       if m == 0:
           return 0.0
       else:
           return ((-1)**m)*(tau(-m)*sp.cos(np.absolute(m)*alpha)-tau(m)*sp.sin(np.absolute(m)*alpha))

   def S(l,m,M,alpha,beta):
       W1 = Rotation.d(l,np.absolute(m),M,beta).doit()
       W2 = Rotation.d(l,np.absolute(m),-M,beta).doit()
       return A(m,alpha)*(((-1)**M)*W1+W2)

   def T(l,m,M,alpha,beta):
       W1 = Rotation.d(l,np.absolute(m),M,beta).doit()
       W2 = Rotation.d(l,np.absolute(m),-M,beta).doit()
       return B(m,alpha)*((-1)**M*W1-W2)

   if M == 0:
       return 2*A(m_1,alpha)*A(m_2,alpha)*Rotation.d(l_1,np.absolute(m_1),0,beta).doit()*Rotation.d(l_2,np.absolute(m_2),0,beta).doit()
   else:
       return (S(l_1,m_1,M,alpha,beta)*S(l_2,m_2,M,alpha,beta) + T(l_1,m_1,M,alpha,beta)*T(l_2,m_2,M,alpha,beta))

def G(l_1,l_2,m_1,m_2,M,alpha,beta):

    # INPUTS: 'l_1': Angular quantum number of atom 1.
    # 'l_2': Angular quantum number of atom 2.
    # 'm_1': Magnetic quantum number of atom 1.
    # 'm_2': Magnetic quantum number of atom 2.
    # 'M': A float which represents the type of symmetry bond we are considering, e.g. M = 0
    # for sigma, M = 1 for pi, M = 2 for delta, etc.
    # 'alpha': Azithmuthal coordinate.
    # 'beta': Polar coordinate.
    # RETURNS: This function returns the relevant coefficient that should be are used in writing the
    # Slater-Koster transformations.

    def tau(m):
        if m >= 0:
            return 1.0
        else:
            return 0.0

    def A(m,alpha):
        if m == 0:
            return np.sqrt(2)/2
        else:
            return ((-1)**m)*(tau(m)*sp.cos(np.absolute(m)*alpha)+tau(-m)*sin(np.absolute(m)*alpha))

    def B(m,alpha):
        if m == 0:
            return 0.0
        else:
            return ((-1)**m)*(tau(-m)*sp.cos(np.absolute(m)*alpha)-tau(m)*sin(np.absolute(m)*alpha))

    def S(l,m,M,alpha,beta):
        W1 = Rotation.d(l,np.absolute(m),M,beta).doit()
        W2 = Rotation.d(l,np.absolute(m),-M,beta).doit()
        return A(m,alpha)*(((-1)**M)*W1+W2)

    def T(l,m,M,alpha,beta):
        W1 = Rotation.d(l,np.absolute(m),M,beta).doit()
        W2 = Rotation.d(l,np.absolute(m),-M,beta).doit()
        return B(m,alpha)*((-1)**M*W1-W2)

    if M == 0:
        return 2*A(m_1,alpha)*A(m_2,alpha)*Rotation.d(l_1,np.absolute(m_1),0,beta).doit()*Rotation.d(l_2,np.absolute(m_2),0,beta).doit()
    else:
        return (S(l_1,m_1,M,alpha,beta)*S(l_2,m_2,M,alpha,beta) + T(l_1,m_1,M,alpha,beta)*T(l_2,m_2,M,alpha,beta))

def signmatrix(l_1,l_2):

    # INPUTS: 'l_1': Angular quantum number of atom 1.
    # 'l_2': Angular quantum number of atom 2.
    # RETURNS: A numpy array that will be used to correct for sign errors present when implementing the code into FHI-aims.

    S = np.ones((int(2*l_1+1),int(2*l_2+1)))

    for i in range(int(2*l_1+1)):
        for j in range(int(2*l_2+1)):
            if i - l_1 > 0 and (i - l_1) % 2 == 1:
                S[i,j] *= -1
            if j - l_2 > 0 and (j - l_2) % 2 == 1:
                S[i,j] *= -1

    return S

def gmatrix(l_1,l_2,sym,labels,Rs):

    # INPUTS: 'l_1': Angular quantum number of atom 1.
    # 'l_2': Angular quantum number of atom 2.
    # 'sym': A float which represents the type of symmetry bond we are considering, e.g. sym = 0
    # for sigma, sym = 1 for pi, sym = 2 for delta, etc.
    # 'labels': Represents which two atoms we are working with.
    # 'Rs': The position vectors of all the atoms within the molecule.
    # RETURNS: This function returns a matrix of Slater-Koster coefficients that is used to transform
    # the submatrices of a matrix between a Cartesian basis to a rotationally invariant
    # basis.

    v1 = Rs[labels[0]-1]
    v2 = Rs[labels[1]-1]
    vecdiff = np.subtract(v2,v1)

    x = vecdiff[0]
    y = vecdiff[1]
    z = vecdiff[2]

    g = np.empty((int(2*l_1+1),int(2*l_2+1)))

    for m in range(-int(l_1),int(l_1+1)):
        for M in range(-int(l_2),int(l_2+1)):
            alpha = carttospher(x,y,z)[1]
            beta = carttospher(x,y,z)[2]
            g[int(l_1+m),int(l_2+M)] = re(simplify(G(int(l_1),int(l_2),m,M,sym,alpha,beta)))

    return np.multiply(g,signmatrix(l_1,l_2))
