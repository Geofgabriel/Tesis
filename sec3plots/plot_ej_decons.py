# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import numpy as np
import matplotlib.pyplot as plt
from espe_fun import espc1d
import ondiculas as on

from matplotlib import rc, font_manager

#rc('text', usetex=True)
#plt.rc('font', family='serif')
#plt.rcParams.update({'font.size': 15})

path_d = r'/home/gaby/Documents/Tesis_v6/Graficos_suplementarios/Deconvolucion/v6/sec3plots/data.txt'
path_rori = r'/home/gaby/Documents/Tesis_v6/Graficos_suplementarios/Deconvolucion/v6/sec3plots/refle_.txt'
path_rq = r'/home/gaby/Documents/Tesis_v6/Graficos_suplementarios/Deconvolucion/v6/sec3plots/refle_l2.txt'
path_rnq = r'/home/gaby/Documents/Tesis_v6/Graficos_suplementarios/Deconvolucion/v6/sec3plots/refle_cauchy.txt'


dato = np.loadtxt(path_d)[27:-27]
r = np.loadtxt(path_rori)
if len(r)!=len(dato):
    print('r',len(r))
    print('data',len(dato))
    raise ValueError

r_l2 = np.loadtxt(path_rq)
r_l1 = np.loadtxt(path_rnq)

dt = 0.002
t_axis = np.arange(0,len(r))*dt


f,d_f = espc1d(dato)
f,r_f = espc1d(r)
f,rl1_f = espc1d(r_l1)
f,rl2_f = espc1d(r_l2)


f1_1 = 5;f2_1 = 25;f3_1 = 100;f4_1 = 120 # 2,30,120,150
w1 = on.ondi_ormsby(101,2,f1_1,f2_1,f3_1,f4_1,fase=0)
w1 = w1/np.max(np.abs(w1))
r_l1_w = np.convolve(r_l1,w1,mode='full')[26:-24]
f,r_l1_w_f = espc1d(r_l1_w)

#=======================================================
# ESPECTROS

limf = 250
fig1,ax1 = plt.subplots(figsize=(7,6))
ax1.plot(f,d_f,'k',label='dato')
ax1.plot(f,rl2_f,'k',alpha=0.7,label='q')
ax1.plot(f,r_l1_w_f,'k',alpha=0.3,label='nq')
ax1.set_xlim(0,limf)
ax1.set_xlabel('Frecuencia(Hz)',fontsize=15)
ax1.set_ylabel('Amplitud',fontsize=15)
ax1.legend()

fig1.savefig('ejemplo_deconvoluciones_espectros_tesisv6_v2.pdf',bbox_inches='tight')
plt.show()
#========================================================
# DECONVOLUCIONES

fig,ax = plt.subplots(1,5,figsize=(10,8))
labels_list=[r'a)',r'b)',r'c)',r'd)',r'e)']

for i in range(0,5): 
    ax[i].set_xlabel('Amplitud')
    ax[i].set_ylim(0,.6)
    #ax[i].set_xticks([])
    
    if i > 0 :
        ax[i].set_yticks([])


ax[0].plot(r,t_axis,'k')
ax[0].invert_yaxis()
ax[0].set_ylabel('Tiempo(s)')
ax[0].text(-0.1,1.03, labels_list[0], 
           transform=ax[0].transAxes) 
ax[0].set_yticks([])

ax[1].plot(dato,t_axis,'k')
ax[1].invert_yaxis()
ax[1].text(-0.1,1.03, labels_list[1], 
           transform=ax[1].transAxes) 

ax[2].plot(r_l2,t_axis,'k')
ax[2].invert_yaxis()
ax[2].text(-0.1,1.03, labels_list[2], 
           transform=ax[2].transAxes) 

ax[3].plot(r_l1,t_axis,'k')
ax[3].invert_yaxis()
ax[3].text(-0.1,1.03, labels_list[3], 
           transform=ax[3].transAxes) 

ax[4].plot(r_l1_w,t_axis,'k')
ax[4].invert_yaxis()
ax[4].text(-0.1,1.03, labels_list[4], 
           transform=ax[4].transAxes) 

fig.savefig('ejemplo_deconvoluciones_tesisv6_v2.pdf',bbox_inches='tight')
plt.show()
