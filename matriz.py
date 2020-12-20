# -*- coding: utf-8 -*-
"""
Created on Fri Aug 21 15:30:53 2020
@author: naomy
"""
#NESTE CÓDIGO, CALCULO, ATRAVÉS DE MATRIZES, O VETOR DENSIDADE DE ESTADOS 
# n = [22 33 12 21 13 31 23 32], com \rho_{11} = 1 - \rho_{22} - \rho{33}
import numpy as np
from numpy.linalg import inv
import matplotlib.pyplot as plt

deltac = 2*np.pi*0.1*10**6
omegapp = 2*np.pi*0.1*10**6 #para ser usado na parte analítica
omegap = (1j*2*np.pi*0.1*10**6)/2
omegac = (1j*20.7*10**6)/2
omegacc = 20.7*10**6
g = 0#2*np.pi*0.01*10**6 #\gamma_{deph}
gammap = 2*np.pi*6.06*10**6#valor tirado do data sheet do Rb87
N = 5*10**14
kp = 2*np.pi/(780*10**-9)
z = 0.005 #7cm
epsilon0 = 8.854*10**(-12)
hbar = 1.054*10**-34
mu12 = 3.584*10**-29
x = np.arange(-70*10**6, 70*10**6, 1*10**5) #\Delta_{p}


#------------------------------------- ANALITICO -------------------------------
Nom = 4*gammap*(x + deltac)**2
Den = 4*((x + deltac)**2)*gammap**2 + (omegacc**2 - 4*x*(x + deltac))**2
A = (2*N*kp*Nom*(mu12**2))/(epsilon0*hbar*Den)#absorção para o feixe gaussiano

rho22_num = (omegapp**2)*Nom/(Den*gammap)
#--------------------------------------------------------------------------------

#------------------------------------ NUMERICO -------------------------------------
rho_array = np.empty((0))
for i in range (-70*10**6, 70*10**6, 1*10**5):#loop para calcular todos os valores das funções
    #para os diferentes valores de deltap
    g12 = (gammap + 2j*i)/2
    g13 = 1j*(i+deltac)
    g23 = (omegap + 2j*deltac)/2
    #agora os complexos conjugados:
    g12c = (gammap - 2j*i)/2
    g13c = -1j*(i+deltac)
    g23c = (omegap - 2j*deltac)/2 
    

    M = np.array([[-gammap, 0, -omegap, omegap, 0, 0, omegac, -omegac], [0, 0, 0, 0, 0, 0, -omegac, omegac], [-2*omegap, -omegap, -g12, 0, omegac, 0, 0, 0], [2*omegap, omegap, 0, -g12c, 0, -omegac, 0, 0], [0, 0, omegac, 0, -g13, 0, -omegap, 0], [0, 0, 0, -omegac, 0, -g13c, 0, omegap], [omegac, -omegac, 0, 0, -omegap, 0, -g23, 0], [-omegac, omegac, 0, 0, 0, omegap, 0, -g23c]])
    Z = inv(M) #inversa de uma matriz
    b = np.array([0, 0, omegap, -omegap, 0, 0, 0, 0])
    n = np.matmul(-Z, b) #multiplicação de matriz
    rho = (n[2].imag)*(2*N*kp*(mu12**2))/(epsilon0*hbar*omegapp) #\rho_{12}
    #rho = n[0].real #\rho_{22}
    #rho = n[1].real #\rho_{33}
    #rho = 1 - n[0].real - n[1].real #\rho_{11}
    rho_array = np.append(rho_array, [rho], axis=0)

plt.plot(x, rho_array, color='g', label='Numerical')
plt.plot(x, rho22_num, color='r', label='Analytical', linestyle=':')
plt.xlabel('$\Delta_{p}$(MHz)')
plt.ylabel('Absorption')
plt.legend(loc='upper right')
#plt.savefig('EIT_comparacao_matriz_e_aproximacao.pdf', format='pdf', dpi=1000)
plt.show()
#O RESULTADO DESSE GRÁFICO DEVE SER O MESMO DAQUELE COM APROXIMAÇÃO, PORQUE 
#OMEGAC >> OMEGAP 