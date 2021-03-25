# -*- coding: utf-8 -*-
"""
Created on Thu Mar 25 10:23:54 2021

@author: naomy
"""


import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integrate
from operator import sub

c = 3*10**8 
deltac = 0.1*10**6
#omegapp = 2*np.pi*0.1*10**6#para ser usado na parte analítica
omegap = 0.1*10**6
omegac = 80*10**6
gammap = 2*np.pi*6.06*10**6#valor tirado do data sheet do Rb87
g = 12*10**6
N = 5*10**14
kp = 2*np.pi/(780*10**-9)
kc = 2*np.pi/(776*10**-9)
z = 0.8 #7cm
epsilon0 = 8.854*10**(-12)
hbar = 1.054*10**-34
mu12 = (3.584*10**-29)*np.sqrt(5/28)##momento de dipolo alterado pelo campo magnetico
u = 290.5 #root mean square velocity
x = np.arange(-100*10**6, 100*10**6, 1*10**6) #\Delta_{p}
w0p = 1*10**-3 #integral é feita na região da cintura do feixe de prova
w0c = 1.5*10**-3
diff = 10*2.8*0.6182*10**6


result_array = np.empty((0))
for i in range (-100*10**6, 100*10**6, 1*10**6):
    #feixe Gaussiano com alargamento doppler. A integral é feita somente na distribuição de velocidade dos átomos
    A2 = integrate.quad(lambda v: ((np.exp(-v**2/u**2)*2*N*kp*mu12**2)/(u*np.sqrt(np.pi)*epsilon0*hbar))*(-4*(i+deltac-kp*v+kc*v)*((i-kp*v+diff)*g + gammap*(i+deltac-kp*v+kc*v)) + g*(omegac**2+g*gammap-4*(i-kp*v+diff)*(i+deltac-kp*v+kc*v)))/((omegac**2 + g*gammap -4*(i-kp*v+diff)*(i+deltac-kp*v+kc*v))**2 + 4*((i-kp*v+diff)*g + gammap*(i+deltac-kp*v+kc*v))**2), -np.inf, np.inf)
    #absorçãodo feixe de prova sigma+
    A1 = integrate.quad(lambda v: ((np.exp(-v**2/u**2)*2*N*kp*mu12**2)/(u*np.sqrt(np.pi)*epsilon0*hbar))*(-4*(i+deltac-kp*v+kc*v)*((i-kp*v-diff)*g + gammap*(i+deltac-kp*v+kc*v)) + g*(omegac**2+g*gammap-4*(i-kp*v-diff)*(i+deltac-kp*v+kc*v)))/((omegac**2 + g*gammap -4*(i-kp*v-diff)*(i+deltac-kp*v+kc*v))**2 + 4*((i-kp*v-diff)*g + gammap*(i+deltac-kp*v+kc*v))**2), -np.inf, np.inf)
    #absorção do feixe de prova sigma-
    A = tuple(map(sub, A1, A2))
    result = A[0]
    result_array = np.append(result_array, [result], axis=0)
    
plt.plot(x, result_array, label='Gaussian', color='b')
