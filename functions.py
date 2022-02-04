"""
This file contains all of the functions needed to run the form factor calculations.
"""

import numpy as np
import cubature
import cmath as cm
import math

PI = math.pi
EPS = 1e-6

def kallen(a,b,c):
  return a**2 + b**2 + c**2 - 2*(a*b + a*c + b*c)

def sig(m,s):
  return cm.sqrt(1 - 4*m**2/s)

def pion_ff(x,S,a):
  mq=.1, mpi=.14, mu=0.5, m=.6, M=.3
  M2 = mq**2 - mpi**2 - mu**2
  alpha = a / 4 / PI**2
  z = x[0]
  y = x[1]
  s = math.tan(y)
  c = s / m**2
  J = 1 / math.cos(y)**2
  k_squared = mq**2 + mpi**2 - s/2*( 1 - z * sig(mpi,s) * sig(mq,s) )
  k = cm.sqrt(k_squared)
  def A(kk):
    k1 = cm.sqrt( kallen(kk,m,M) )
    n1 = np.log(m**2 / M**2)
    n2 = (kk - m**2 - M**2) / k1 * cm.log( (kk - m**2 - M**2 - k1) / (kk - m**2 - M**2 + k1) )
    return   alpha * 4. / 3. * ( n1 + n2 )
  def B(kk):
    k1 = cm.sqrt( kallen(kk,m,M) )
    n1 = np.log(M**2 / m**2)
    n2 = (kk - m**2 + M**2) / k1 * cm.log( (kk - m**2 - M**2 + k1) / (kk - m**2 - M**2 - k1) )
    return  M * alpha * 4. / 3. / kk * ( n1 + n2 )
  a_a = A(k_squared)
  b_b = B(k_squared)
  num1 = np.cosh( b_b * k * c) + mu * np.sinh( b_b * k * c) / k
  num2 = c**a_a
  num3 = (s + 2.*(mu-mq)**2 - 2.*mpi**2)*M2/(s + 4.*mpi**2) + (mu - mq)**2 - mpi**2
  den = k_squared - mu**2
  return np.real( num1 * num2 * num3 / den * J / ( s - S - 1j*EPS) * sig(mq,s) )
  
