# -*- coding: utf-8 -*-
"""
Created on Thu Oct 31 09:20:41 2019

@author: Gabriel
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.signal as ss

def rickerpy(leng,dt,f0,fase=0):

    M = round(leng/(2*dt)) 
    
    N = int(2*M+1)
    
    rick = np.zeros((N,1))
    for n in range(0,N):
        t = dt*(n-M-1)/1000 # para que este en ms
        aux = np.pi*f0*t
        rick[n]= (1-2*aux**2)*np.exp(-aux**2)

    angles = np.radians(fase)
    w_rot = np.zeros((len(rick),1))
    w_rot = np.cos(angles)*rick[:,0]-np.sin(angles)*np.imag(ss.hilbert(rick[:,0]))
    
    return w_rot

def ondi_ormsby(leng, dt, f1,f2,f3,f4,fase=0,taper=True):
    
    M = round(leng/(2*dt))
    N = int(2*M+1)
    
    w = np.zeros(N)
    
    def numerator(f, t):
        return (np.sinc(f * t)**2) * ((np.pi * f) ** 2)

    #t = np.arange(-length/2, length/2,dt)
    for n in range(0,N):        
        t = dt*(n-M-1)/1000
        pf43 = (np.pi * f4) - (np.pi * f3)
        pf21 = (np.pi * f2) - (np.pi * f1)
        w[n] = ((numerator(f4, t)/pf43) - (numerator(f3, t)/pf43) -
             (numerator(f2, t)/pf21) + (numerator(f1, t)/pf21))
    
    if taper == True:
        w = w * np.hamming(N)
        
    angles = np.radians(fase)
    w_rot = np.zeros(N)
    w_rot = np.cos(angles)*w-np.sin(angles)*np.imag(ss.hilbert(w))
    
    
    return w_rot

if __name__ == '__main__':
    fase = 30
    f1 = 5
    f2 = 10
    f3 = 65
    f4 = 70
    dt = 2
    l = 208
    #ondicula = ondi_ormsby(l,dt,f1,f2,f3,f4,fase,taper=True)
    ondicula = rickerpy(l,dt,f3,fase)
    plt.plot(ondicula)