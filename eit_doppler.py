# -*- coding: utf-8 -*-
"""
Created on Mon Dec 14 09:11:18 2020

@author: naomy
"""


import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integrate

c = 3*10**8 
deltac = 2*np.pi*0.1*10**6
omegapp = 2*np.pi*0.1*10**6#para ser usado na parte analítica
omegap = 0.1*10**6
omegac = 50*10**6
gammap = 2*np.pi*6.06*10**6#valor tirado do data sheet do Rb87
g = 0
N = 5*10**14
kp = 2*np.pi/(780*10**-9)
kc = 2*np.pi/(480*10**-9)
z = 0.1 #7cm
epsilon0 = 8.854*10**(-12)
hbar = 1.054*10**-34
mu12 = 3.584*10**-29
u = 290.5 #root mean square velocity
x = np.arange(-70*10**6, 70*10**6, 1*10**6) #\Delta_{p}
w0p = 140*10**-6 #integral é feita na região da cintura do feixe de prova
w0c = 180*10**-6

result_array = np.empty((0))
for i in range (-70*10**6, 70*10**6, 1*10**6):
    #feixe Gaussiano com alargamento doppler. A integral é feita somente na distribuição de velocidade dos átomos
    A = integrate.quad(lambda v: ((np.exp(-v**2/u**2)*8*N*kp*gammap*mu12**2)/(u*np.sqrt(np.pi)*epsilon0*hbar))*((i+deltac-kp*v+kc*v)**2/(4*gammap**2*(i+deltac-kp*v+kc*v)**2 + (omegac**2-4*(i-kp*v)*(i+deltac-kp*v+kc*v))**2)), -np.inf, np.inf)
  
    result = np.exp(-A[0]*z)
    result_array = np.append(result_array, [result], axis=0)
    #result = A[0]
    #--------------------------------------------------LG---------------------------------------
    #feixe LG com alargamento doppler. A integral é feita na distribuição de velocidades e no raio do feixe rosquinha
result_arrayLG = np.empty((0))
for n in range (-70*10**6, 70*10**6, 1*10**6):
    f = lambda v, r: r*((np.exp(-v**2/u**2)*8*N*kp*gammap*mu12**2)/(u*np.sqrt(np.pi)*epsilon0*hbar))*((n+deltac-kp*v+kc*v)**2/(4*gammap**2*(n+deltac-kp*v+kc*v)**2 + ((omegac*np.sqrt(2)*(r/w0c)*np.exp(-(r/w0c)**2))**2-4*(n-kp*v)*(n+deltac-kp*v+kc*v))**2))
    A_LG = integrate.dblquad(f, 0, w0p, lambda v: -np.inf, lambda v: np.inf)#integro primeiro em v e por último em r
    result2 = np.exp(-A_LG[0]*100000000*z)
    result_arrayLG = np.append(result_arrayLG, [result2], axis=0)
    
    


plt.plot(x, result_array, label='Gaussian', color='r')
plt.plot(x, result_arrayLG, label='LG', color='g')
plt.xlabel('$\Delta_{p}$(MHz)')
plt.ylabel('Transmission')
plt.legend()
#plt.savefig('EIT.pdf', format='pdf', dpi=1000)
plt.show()
