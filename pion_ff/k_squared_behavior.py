import numpy as np
import matplotlib.pyplot as ply
import functions

mpi = 0.14
mq = 0.1
ss = [1e3, 1e4, 1e5, 1e6]
num_z = 1000
zs = np.linspace(-1,1,num_z)

def k_squared(z,s):
  return -s/2*( 1 - z * functions.sig(mpi,s) * functions.sig(mq,s) ) + mpi**2 + mq**2

plt.figure(figsize=(13,7))
plt.tight_layout()
for i in range(len(ss)):
  plt.subplots(2,2,i)
  plt.plot( zs, 
            [k_squared(z,ss[i]) for z in zs],
            label=r'$k^2(z,s=${}$)$'.format( ss[i] ) 
          )
  plt.xscale('log')
  plt.legend()
plt.show()
