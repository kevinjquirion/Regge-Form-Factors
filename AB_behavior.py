""" 
This file plots the behavior of the A and B functions for the pion form factor discontinuity
"""

import numpy as np
import functions
import matplotlib.pyplot as plt

kmin = 1e-2
kmax = 1e2
num_k = 1000
alphas = [.3,.4,.5,.6]
k2s = np.logspace(np.log10(kmin), np.log10(kmax), num_k)

As = [ np.real( functions.A_pion(-k2,alphas[0]) ) for k2 in k2s ]
Bs = [ np.real( functions.B_pion(-k2,alphas[0]) ) for k2 in k2s ]

plt.figure(figsize=(13,7))
plt.tight_layout()

for i in range(len(alphas)):
    plt.subplot(2,2,i+1)
    plt.plot( k2s,
                [ np.real( functions.A_pion(-k2,alphas[i]) ) for k2 in k2s ],
                label=r"$A(k^2,\alpha=${}$)$".format( alphas[i] ) )
    plt.plot( k2s,
                [ np.real( functions.B_pion(-k2,alphas[i]) ) for k2 in k2s ],
                label=r"$B(k^2,\alpha=${}$)$".format( alphas[i] ) )
    plt.xlabel(r'$-k^2$')
    plt.legend(loc='best')
    plt.xscale('log')
plt.savefig('images/AB_behavior.pdf')
plt.show()
