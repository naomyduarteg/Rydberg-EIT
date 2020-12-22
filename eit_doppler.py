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
omegac = 30*10**6
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


result_array = np.empty((0))
for i in range (-70*10**6, 70*10**6, 1*10**6):
    #B = (8*N*kp*gammap*mu12**2)/(epsilon0*hbar) #constantes que entram no cálculo da absorção
    A = integrate.quad(lambda v: ((np.exp(-v**2/u**2)*8*N*kp*gammap*mu12**2)/(u*np.sqrt(np.pi)*epsilon0*hbar))*((i+deltac-kp*v+kc*v)**2/(4*gammap**2*(i+deltac-kp*v+kc*v)**2 + (omegac**2-4*(i-kp*v)*(i+deltac-kp*v+kc*v))**2)), -500, 500)
    #T_LG = integrate.quad(lambda t: np.exp(-B*(gammap*g**2 - g*i*(i + deltac) + g*(omegac*np.sqrt(2)*(t/w0)*np.exp(-(t/w0)**2))**2 + (i + deltac)*(i*g + gammap*(i + deltac)))/(2*((gammap*g - i*(i + deltac) + (omegac*np.sqrt(2)*(t/w0)*np.exp(-(t/w0)**2))**2)**2 + (i*g + gammap*(i + deltac))**2)))*t, 0, w0)
    result = np.exp(-A[0]*z)
 
    #print(A)
    result_array = np.append(result_array, [result], axis=0)
#------------------------------------------------------------FEIXE LG--------------------------------------------------------------------------
#feixe LG com alargamento doppler. A integral é feita na distribuição de velocidades e no raio do feixe rosquinha
for n in range (-70*10**6, 70*10**6, 1*10**6):
    f = lambda v, r: r*((np.exp(-v**2/u**2)*8*N*kp*gammap*mu12**2)/(u*np.sqrt(np.pi)*epsilon0*hbar))*((n+deltac-kp*v+kc*v)**2/(4*gammap**2*(n+deltac-kp*v+kc*v)**2 + ((omegac*np.sqrt(2)*(r/w0)*np.exp(-(r/w0)**2))**2-4*(n-kp*v)*(n+deltac-kp*v+kc*v))**2))
    A_LG = integrate.dblquad(f, 0, w0, lambda v: -800, lambda v: 800)#integro primeiro em v e por último em r
    result2 = np.exp(-A_LG[0]*z)
    result_arrayLG = np.append(result_arrayLG, [result2], axis=0)
    

plt.plot(x, result_array, label='Gaussian', color='r')
plt.plot(x, result_arrayLG, label='LG', color='r')
plt.xlabel('$\Delta_{p}$(MHz)')
plt.ylabel('Transmission')
plt.legend()
#plt.savefig('EIT.pdf', format='pdf', dpi=1000)
plt.show()
