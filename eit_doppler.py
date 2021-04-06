# -*- coding: utf-8 -*-
"""
Created on Mon Dec 14 09:11:18 2020

@author: naomy
"""


import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integrate
from operator import sub

#NESTE PROGRAMA, USAMOS OS VALORES MULTIPLICADOS POR 2PI PARA CALCULAR 
c = 3*10**8 
deltap = 0.1*10**6
omegac = 2*np.pi*6.15*10**6
gammap = 2*np.pi*6.06*10**6#valor tirado do data sheet do Rb87
g = 0
N = 5*10**14
kp = 2*np.pi/(780*10**-9)
kc = 2*np.pi/(480*10**-9)
z = 0.8 #7cm
epsilon0 = 8.854*10**(-12)
hbar = 1.054*10**-34
mu12 = (3.584*10**-29)*np.sqrt(5/28)##momento de dipolo alterado pelo campo magnetico
u = 290.5 #root mean square velocity
x = np.arange(-100*10**6, 100*10**6, 1*10**6) #\Delta_{p}
w0p = 170*10**-6 #integral é feita na região da cintura do feixe de prova
w0c = 190*10**-6
diff = 0#10*2.8*0.6182*10**6


result_array = np.empty((0))
for i in range (-100*10**6, 100*10**6, 1*10**6):
    #feixe Gaussiano com alargamento doppler. A integral é feita somente na distribuição de velocidade dos átomos
    A1 = integrate.quad(lambda v: ((np.exp(-v**2/u**2)*8*N*kp*gammap*mu12**2)/(u*np.sqrt(np.pi)*epsilon0*hbar))*((i+deltap-kp*v+kc*v)**2/(4*gammap**2*(i+deltap-kp*v+kc*v)**2 + (omegac**2-4*(deltap+diff-kp*v)*(i+deltap-kp*v+kc*v))**2)), -np.inf, np.inf)
    #absorçãodo feixe de prova sigma+
    A2 = integrate.quad(lambda v: ((np.exp(-v**2/u**2)*8*N*kp*gammap*mu12**2)/(u*np.sqrt(np.pi)*epsilon0*hbar))*((i+deltap-kp*v+kc*v)**2/(4*gammap**2*(i+deltap-kp*v+kc*v)**2 + (omegac**2-4*(deltap-diff-kp*v)*(i+deltap-kp*v+kc*v))**2)), -np.inf, np.inf)
    #absorção do feixe de prova sigma-
    A = tuple(map(sub, A1, A2))
    result = A1[0]
    result_array = np.append(result_array, [result], axis=0)
    #result = A[0]
    
    
#--------------------------------------------------LG---------------------------------------
    #feixe LG com alargamento doppler. A integral é feita na distribuição de velocidades e no raio do feixe rosquinha
result_arrayLG = np.empty((0))
for n in range (-100*10**6, 100*10**6, 1*10**6):
    f1 = lambda v, r: r*((np.exp(-v**2/u**2)*8*N*kp*gammap*mu12**2)/(u*np.sqrt(np.pi)*epsilon0*hbar))*((n+deltap-kp*v+kc*v)**2/(4*gammap**2*(n+deltap-kp*v+kc*v)**2 + ((omegac*np.sqrt(2)*(r/w0c)*np.exp(-(r/w0c)**2))**2-4*(deltap+diff-kp*v)*(n+deltap-kp*v+kc*v))**2))
    f2 = lambda v, r: r*((np.exp(-v**2/u**2)*8*N*kp*gammap*mu12**2)/(u*np.sqrt(np.pi)*epsilon0*hbar))*((n+deltap-kp*v+kc*v)**2/(4*gammap**2*(n+deltap-kp*v+kc*v)**2 + ((omegac*np.sqrt(2)*(r/w0c)*np.exp(-(r/w0c)**2))**2-4*(deltap-diff-kp*v)*(n+deltap-kp*v+kc*v))**2))
    A_LG1 = integrate.dblquad(f1, 0, w0p, lambda v: -np.inf, lambda v: np.inf)#integro primeiro em v e por último em r
    A_LG2 = integrate.dblquad(f2, 0, w0p, lambda v: -np.inf, lambda v: np.inf)
    A_LG = tuple(map(sub, A_LG1, A_LG2))
    result2 = A_LG1[0]*(2/w0p**2)#fator 2/w0p**2 é normalização da curva
    result_arrayLG = np.append(result_arrayLG, [result2], axis=0)
    
    


plt.plot(x, result_array, label='Gaussian', color='b')
plt.plot(x, result_arrayLG, label='LG', color='black')
plt.xlabel('$\Delta_{p}$(MHz)')
plt.ylabel('Absorption')
plt.legend()
#plt.savefig('EIT_doppler.pdf', format='pdf', dpi=1000)
plt.show()
