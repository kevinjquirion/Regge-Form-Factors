"""
Functin that plots the Reggiezed pion form factor. The main function called is
pion_ff_integrator which takes as an argument an array of s values, s being the
mandelstam variable. An important parameter for this function is err which
determines the relative error of the cubature integrator, and should no be set
to greater than 1e-5. Increasing err wi
"""

import functions
import matplotlib.pyplot as plt
import numpy as np
from cubature import cubature
import time

smin = 1e3
smax = 1e6
num_points = 20
ss = np.logspace( np.log10(smin), np.log10(smax), num_points )
alphas = [.3,.4,.5,.6]

TIC = time.time()

def pion_ff_integrator(SS,i):
  fs = []
  mq = 0.1
  mpi = 0.14
  err = 1e-6
  xmin = [ -1., np.arctan(4*mq**2) ]
  xmax = [ 1., np.pi/2. ]
  f0 = cubature( functions.pion_ff, 2, 1, xmin, xmax, args=[SS[0],alphas[i]], relerr=err  )[0]
  for j,s in enumerate(ss):
    tic = time.time()
    f = cubature( functions.pion_ff, 2, 1, xmin, xmax, args=[s,alphas[i]], relerr=err  )[0]
    fs.append(  f / f0  )
    toc = time.time()
    print("Done with iteration {} out of {}. Time taken: {} s".format(j+1,len(ss), round(toc-tic,2)))
  return fs

f1 = [ ss[0] / s for s in ss ]
f2 = [ ss[0] / np.log(ss[0]) * np.log(s) / s for s in ss]
f3 = [ ss[0] / np.log(np.log(ss[0])) * np.log(np.log(s)) / s for s in ss]

Fs1 = pion_ff_integrator(ss,3)
# print(Fs1)

TOC = time.time()
print("Total time taken: {} minutes".format( round( (TOC-TIC)/60 ,2) ) )

plt.figure(figsize=(12,6))
plt.tight_layout()
plt.subplot(1,1,1)
plt.plot(ss, f1, label=r'$ f(s) = \frac{1}{s}$', c='b')
plt.plot(ss, f2, label=r'$ f(s) = \frac{\ln s}{s}$', c='g')
plt.plot(ss, f3, label=r'$ f(s) = \frac{\ln(\ln s)}{s}$', c='r')
plt.scatter(ss,Fs1,label=r'Reggiezed $F(s)$ with $\alpha=${}'.format(alphas[3]), c='k', s=5)
plt.xscale('log')
plt.yscale('log')
plt.legend()
plt.savefig('images/ReggiezedPionFormFactor.pdf')
plt.show()
