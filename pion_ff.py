import functions
import matplotlib.pyplot as plt
import numpy as np
from cubature import cubature
import time 

smin = 1e3
smax = 1e6
ss = np.logspace( np.log10(smin), np.log10(smax) )
alphas = {.3,.4,.5,.6}

def pion_ff_integrator(SS,i):
  fs = []
  mq = 0.1
  mpi = 0.14
  err = 1e-2
  xmin = [ -1., np.arctan(4*mq**2) ]
  xmax = [ 1., np.pi/2. ]
  for j,s in enumerate(ss):
    f = cubature( functions.pion_ff, 2, 1, xmin, xmax, args=[s,alphas[i]], abserr=err  )
    print("Done with iteration {} out of {}".format(j+1,len(ss)))
  return fs

Fs1 = pion_ff_integrator(ss,0)

plt.figure(figsize=(12,6))
plt.tight_layout()
plt.subplot(1,1,1)
plt.plot(ss,Fs,label=r'Reggiezed $F(s)$ with $\alpha=${}'.format(alphas[0]))
plt.xscale('log')
plt.yscale('log')
plt.savefig('images/ReggiezedPionFormFactor.pdf')
plt.show()
  
