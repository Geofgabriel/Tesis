# -*- coding: utf-8 -*-
"""
Created on Sun Jun 14 17:37:13 2020

@author: Gabriel R. Gelpi
"""

import numpy as np
import matplotlib.pyplot as plt
import espe_fun as ef
import ondiculas as on
import espe_fun as ef
from matplotlib import rc, font_manager
rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rcParams.update({'font.size': 12})

# CARGAMOS LOS DATOS
dato = np.loadtxt('sismica_input_medio_ms.txt')[108:-108,:]
r_hat = np.loadtxt('refle_out_lj1002.txt')

name = 'inline_jl1002' #guarda el dato con el realce en freq.

rs1 = r_hat.shape[0] # muestas
rt1 = r_hat.shape[1] # trazas

dt = 0.0005
f1_1 = 5
f2_1 = 15
f3_1 = 70
f4_1 = 90

f1_2 = 5
f2_2 = 15
f3_2 = 100
f4_2 = 120

f1_3 = 5
f2_3 = 15
f3_3 = 130
f4_3 = 150



w1 = on.ondi_ormsby(108,.5,f1_1,f2_1,f3_1,f4_1,fase=90)
w2 = on.ondi_ormsby(108,.5,f1_2,f2_2,f3_2,f4_2,fase=90)
w3 = on.ondi_ormsby(108,.5,f1_3,f2_3,f3_3,f4_3,fase=90)

#w = on.rickerpy(108,.5,40,90)
w1 = w1/np.max(np.abs(w1))
w2 = w2/np.max(np.abs(w2))
w3 = w3/np.max(np.abs(w3))


new_section1 = np.zeros((rs1,rt1))
new_section2 = np.zeros((rs1,rt1))
new_section3 = np.zeros((rs1,rt1))

#new_section_in = np.zeros((rs,rt))

for i in range(0,rt1): # Genero la seccion de alta resolucion.
    new_section1[:,i] = np.convolve(r_hat[:,i],w1,mode='same')
    new_section2[:,i] = np.convolve(r_hat[:,i],w2,mode='same')
    new_section3[:,i] = np.convolve(r_hat[:,i],w3,mode='same')

new_section1 = new_section1/np.amax(np.abs(new_section1))    
new_section2 = new_section2/np.amax(np.abs(new_section2))    
new_section3 = new_section3/np.amax(np.abs(new_section3))   

#np.savetxt(f'realce_inline_{f1_1}_{f2_1}_{f3_1}_{f4_1}.txt',new_section1) 
#np.savetxt(f'realce_inline_{f1_2}_{f2_2}_{f3_2}_{f4_2}.txt',new_section2) 
#np.savetxt(f'realce_inline_{f1_3}_{f2_3}_{f3_3}_{f4_3}.txt',new_section3) 



fd,espe_dato = ef.espc2d(dato)
fo1,espe_ondi1 = ef.espc1d(w1)
fo2,espe_ondi2 = ef.espc1d(w2)
fo3,espe_ondi3 = ef.espc1d(w3)

fr1,espe_realce1 = ef.espc2d(new_section1)
fr2,espe_realce2 = ef.espc2d(new_section2)
fr3,espe_realce3 = ef.espc2d(new_section3)

fig, ax =plt.subplots(3,1)


ax[0].plot(fd,espe_dato,'k',label='bw')
ax[0].plot(fr1,espe_realce1,'r')
#ax[0].legend(f'{f1_1}{f1_1}{f1_1}{f1_1}')
ax[0].set_ylabel(f'{f1_1}-{f2_1}-{f3_1}-{f4_1}')

ax[0].set_xlim(0,150)
ax[0].set_xticks([])

ax[1].plot(fd,espe_dato,'k')
ax[1].plot(fr2,espe_realce2,'r')

ax[1].set_xlim(0,150)
ax[1].set_xticks([])
ax[1].set_ylabel(f'{f1_2}-{f2_2}-{f3_2}-{f4_2}')

ax[2].plot(fd,espe_dato,'k')
ax[2].plot(fr3,espe_realce3,'r')

ax[2].set_xlim(0,150)
ax[2].set_xlabel('F(Hz)')

ax[2].set_ylabel(f'{f1_3}-{f2_3}-{f3_3}-{f4_3}')
#ax[2].set_xticks([])

fig.savefig('realces_espectros_inline.pdf',bbox_inches='tight')


fig2,ax1 = plt.subplots(1,4)

for i in range(0,4): ax1[i].set_xlabel('Traza')
ax1[0].imshow(dato,aspect='auto',cmap='seismic',vmin=-1,vmax=1)
ax1[0].set_title('Dato')
ax1[0].set_yticks(np.arange(0,rs1,200))
ax1[0].set_yticklabels(np.arange(0,rs1,200)*dt)

ax1[1].imshow(new_section1,aspect='auto',cmap='seismic',vmin=-1,vmax=1)
ax1[1].set_title(f'{f1_1}-{f2_1}-{f3_1}-{f4_1}')
ax1[1].set_yticks([])


ax1[2].imshow(new_section2,aspect='auto',cmap='seismic',vmin=-1,vmax=1)
ax1[2].set_title(f'{f1_2}-{f2_2}-{f3_2}-{f4_2}')
ax1[2].set_yticks([])

ax1[3].imshow(new_section3,aspect='auto',cmap='seismic',vmin=-1,vmax=1)
ax1[3].set_title(f'{f1_3}-{f2_3}-{f3_3}-{f4_3}')
ax1[3].set_yticks([])

fig2.savefig('realces_inline.pdf',dpi=500,bbox_inches='tight')

