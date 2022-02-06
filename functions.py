"""
This file contains all of the functions needed to run the form factor calculations.
"""

import numpy as np
import cmath as cm
import math

PI = math.pi
EPS = 1e-6

def kallen(a,b,c):
  return a**2 + b**2 + c**2 - 2*(a*b + a*c + b*c)

def sig(m,s):
  return cm.sqrt(1 - 4*m**2/s)

def A_pion(k_squared,a):
    m=.6
    M=.3
    alpha = a / 4 / PI**2
    k1 = cm.sqrt( kallen(k_squared,m,M) )
    n1 = np.log(m**2 / M**2)
    n2 = (k_squared - m**2 - M**2) / k1 * cm.log( (k_squared - m**2 - M**2 - k1) / (k_squared - m**2 - M**2 + k1) )
    return   alpha * 4. / 3. * ( n1 + n2 )

def B_pion(k_squared,a):
    m=.6
    M=.3
    alpha = a / 4 / PI**2
    k1 = cm.sqrt( kallen(k_squared,m,M) )
    n1 = np.log(M**2 / m**2)
    n2 = (k_squared - m**2 + M**2) / k1 * cm.log( (k_squared - m**2 - M**2 + k1) / (k_squared - m**2 - M**2 - k1) )
    return  M * alpha * 4. / 3. / k_squared * ( n1 + n2 )

def pion_disc(x,s,a):
    mq=.1
    mpi=.14
    mu=0.5
    M2 = mq**2 - mpi**2 - mu**2
    z = x[0]
    c = s / m**2
    k_squared = mq**2 + mpi**2 - s/2*( 1 - z * sig(mpi,s) * sig(mq,s) )
    k = cm.sqrt(k_squared)
    a_a = A_pion(k_squared,a)
    b_b = B_pion(k_squared,a)
    # num1 = np.cosh( b_b * k * c) + mu * np.sinh( b_b * k * c) / k
    num1 = 1.
    num2 = c**a_a
    num3 = (s + 2.*(mu-mq)**2 - 2.*mpi**2)*M2/(s + 4.*mpi**2) + (mu - mq)**2 - mpi**2
    den = k_squared - mu**2
    return np.real( num1 * num2 * num3 / den  * sig(mq,s) )

def pion_ff(x,S,a):
    z = [x[0]]
    y = x[1]
    s = math.tan(y)
    J = 1 / math.cos(y)**2
    delta = pion_disc(z,s,a)
    return np.real( delta * J / ( s - S - 1j*EPS) )
