# -*- coding: utf-8 -*-
"""
Created on Thu May  7 21:10:46 2020

@author: Gabriel R. Gelpi
"""

import numpy as np



def espc2d(dato,normalizado=1):
    ns = dato.shape[0]
    nt = dato.shape[1]
    dt = 0.0005

    espe=np.abs(np.fft.fft(dato,axis=0))

    espe=(np.sum(espe,axis=1)/nt)   
    f = np.arange(0,ns,1)/(ns*dt)
    
    if normalizado == 1:
        espe = espe/np.max(np.abs(espe))
    
    return f,espe

def espc1d(dato,normalizado=1):
    ns = len(dato)
    dt = 0.002

    espe=np.abs(np.fft.fft(dato))
    
    if normalizado == 1:
        espe = espe/np.max(np.abs(espe))
    
    f = np.arange(0,ns,1)/(ns*dt)
    
    return f,espe