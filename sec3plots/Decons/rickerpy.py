"""
script para generar una ondicula de ricker.

input:
f0: frecuencia central
leng: longitud en tiempo
dt: intervalo de muestreo en tiempo

output:
rick: ondícula de Ricker
"""

import numpy as np

def rickerpy(leng,dt,f0):

    M = round(leng/(2*dt)) 
    
    N = int(2*M+1)
    
    rick = np.zeros((N,1))
    for n in range(0,N):
        t = dt*(n-M-1)/1000 # para que este en ms
        aux = np.pi*f0*t
        rick[n]= (1-2*aux**2)*np.exp(-aux**2)

    return rick


