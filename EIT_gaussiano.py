# -*- coding: utf-8 -*-
"""
Created on Thu Jul 30 18:49:55 2020

@author: naomy
"""



import matplotlib.pyplot as plt
import numpy as np
import scipy.integrate as integrate

#AQUI EU VOU CALCULAR A INTEGRAL DA TRANSMISSÃO COMO FUNÇÃO DE R, FAZENDO UM LOOP NOS VALORES DE X
deltac = 2*np.pi*0.1*10**6
omegap = 2*np.pi*1*10**6
omegac = 3*10.7*10**6
g = 0#2*np.pi*0.01*10**6 #\gamma_{deph}
gammap = 2*np.pi*6*10**6#valor tirado do data sheet do Rb87
x = np.arange(-70*10**6, 70*10**6, 1*10**5) #\Delta_{p}
N = 1*10**10
kp = 2*np.pi/(780*10**-9)
z = 0.001 #7cm
epsilon0 = 8.854*10**(-12)
hbar = 1.054*10**-34
mu12 = 3.584*10**-29

teste1 = omegap/gammap
teste2 = omegac/gammap


Nom = gammap*g**2 - g*x*(x + deltac) + g*omegac**2 + (x + deltac)*(x*g + gammap*(x + deltac))
#Nom = gamma13*(4*x**2+x*deltac-2*gammap*gamma13-omegac**2) + 4*(x+deltac)*(x*gamma13+2*x*gammap-2*deltac*gammap)
Den = 2*((gammap*g - x*(x + deltac) + omegac**2)**2 + (x*g + gammap*(x + deltac))**2)
A = (N*kp*Nom*(mu12**2))/(epsilon0*hbar*Den)#absorção para o feixe gaussiano
T = np.exp(-A*z) #transmissão do feixe gaussiano



#--------------PARTE QUE CALCULA O EIT PARA O FEIXE LG, INTEGRANDO A TRANSMISSÃO EM R-----------------

w0 = 270*10**(-6)#cintura do feixe LG

result_array = np.empty((0))
for i in range (-70*10**6, 70*10**6, 1*10**5): #LOOP PARA INTEGRAR EM TODOS OS VALORES DE X
#Nom = gammap*g**2 - g*x*(x + deltac) + g*(omegac*np.sqrt(2)*(t/w0)**2) + (x + deltac)*(x*g + gammap*(x + deltac))
#Nom = gamma13*(4*x**2+x*deltac-2*gammap*gamma13-omegac**2) + 4*(x+deltac)*(x*gamma13+2*x*gammap-2*deltac*gammap)
#Den = 2*((gammap*g - x*(x + deltac) + (omegac*np.sqrt(2)*(t/w0)*np.exp((-t/w0)**2))**2)**2 + (x*g + gammap*(x + deltac))**2)
    B = (N*kp*z*mu12**2)/(epsilon0*hbar) #constantes que entram no cálculo da absorção

    T_LG = integrate.quad(lambda t: np.exp(-B*(gammap*g**2 - g*i*(i + deltac) + g*(omegac*np.sqrt(2)*(t/w0)*np.exp(-(t/w0)**2))**2 + (i + deltac)*(i*g + gammap*(i + deltac)))/(2*((gammap*g - i*(i + deltac) + (omegac*np.sqrt(2)*(t/w0)*np.exp(-(t/w0)**2))**2)**2 + (i*g + gammap*(i + deltac))**2)))*t, 0, w0)
    #t é o r, sei lá por que usei t 
    result = 2*T_LG[0]/w0**2

    print(result)
    result_array = np.append(result_array, [result], axis=0)
    

#NESSA PARTE EU VOU PLOTAR O DELTAP E A TRANSMISSÃO CALCULADA PARA CADA VALOR E PARA CADA FEIXE
print(teste1, teste2)
plt.plot(x, result_array, label='LG', color='r')
plt.plot(x, T, label='Gaussiano', color='g')
plt.xlabel('$\Delta_{p}$(MHz)')
plt.ylabel('Transmissão')
plt.legend()
plt.savefig('EIT.pdf', format='pdf', dpi=1000)
plt.show()