"""
Script para generar dato sismico sintetico
by Gabriel R. Gelpi

El dato se genera realizando la convolución entre la reflectividad
y ondícula mas algún ruido aditivo.

s = r*w+n
"""

import numpy as np
import matplotlib.pyplot as plt
from rickerpy import rickerpy
from conv_mat import convmtx
import scipy.signal as ss
import scipy.stats as se
from matplotlib import rcParams
rcParams['text.usetex'] = False
plt.rc('font', family='serif')#serif,monospace,sans
plt.rcParams.update({'font.size': 20})
plt.rcParams.update({'lines.linewidth':1.})


#===============================================================
# ONDICULA
#===============================================================

# cargo la ondicula
#wavelet = np.float32(np.loadtxt(argv[1],ndmin=2))    #cargo la ondicula de un archivo .txt

# o genero la ondicula --> Rickerpy.py
longi = 108
f0 = 30
dt = 2

wavelet = rickerpy(longi,dt,f0)                      # genero la ondicula
wl = wavelet.shape[0]                                # largo de la ondicula
   
fase = -30                                           # rotacion a la fase si                                                                                     # es que necesito. 
angles = np.radians(fase)
wavelet_rot = np.zeros((wl,1))
wavelet_rot = np.cos(angles)*wavelet[:,0]-np.sin(angles)*np.imag(ss.hilbert(wavelet[:,0])) #THG
#np.savetxt('ondicula_sintetico.txt',wavelet_rot)
#=================================================================
# MODELO DEL SUBSUELO - REFLECTIVIDAD 
#=================================================================
nr = 300                                             # Nro de muestras de la reflectividad
ntr = 20                                             # numero de reflectividades

R = np.zeros((nr,ntr)) 
refle = r'/home/gaby/Documents/Tesis_v6/Graficos_suplementarios/Deconvolucion/v6/sec3plots/Decons/R_marmousi.txt'
R = np.loadtxt(refle)# uso una porción del modelo Marmousi2.
nr = R.shape[0]
ntr = R.shape[1]

#==============================================================
# DATO SINTETICO
#==============================================================
dl = nr+wl-1 # nro de elementos que va a tener las trazas generadas.

d = np.zeros((dl,ntr))
for j in range(0,ntr):
    d[:,j] = np.convolve(R[:,j],wavelet_rot) # convolucion traza x traza.

#si queremos usar la forma matricial:
#X = convmtx(wavelet_rot,nr) # Matriz de convolucion
#d = np.dot(X,R) # Genero el dato: dato=Wx.

SNR = 5
ruido = np.random.randn(dl,ntr)
ruido = ruido/np.amax(np.abs(ruido))
alfa = np.max(np.abs(d))/(np.max(np.abs(ruido))*SNR)            
d = d + alfa*ruido

nro = 'marmo'
#np.savetxt('dato_%s.txt'%nro,d)

#====================================================
# Graficos
#====================================================
dt=0.002
m=np.int16((wl-1)/2)
colormap = 'gray'# defino el mapa de colores.

xfin = d.shape[0]

plot = 1
if plot==1:
    fig, ax = plt.subplots(1,4,figsize=(15,8))  # base x altura
    plt.subplots_adjust(wspace =0.3,hspace=0.9) # espacio entre graficos
    
    ax[0].set_title('Dato')
    ax[1].set_title('Ruido')
    ax[2].set_title('Reflectividad')
    ax[3].set_title( r'Ond\'icula')
    [ax[i].set_xlabel('Trazas') for i in range(0,3)]
    [ax[i].set_yticklabels([]) for i in range(0,4)]

    ax[0].imshow(d[m:-1-m,:],aspect='auto',cmap=colormap)
    ax[1].imshow(ruido[m:dl-m,:] ,aspect='auto', cmap=colormap)
    ax[1].set_ylabel(r'$\parallel$')
    ax[2].imshow(R, aspect='auto', cmap=colormap)
    ax[2].set_ylabel('+')
    ax[3].plot(wavelet,np.arange(-m,m+1)*dt,'k')
    ax[3].set_ylabel('*')
else:
    trace = 1 # traza a graficar si hay mas de una
    
    fig, ax = plt.subplots(1,4,figsize=(15,8))# base x altura
    plt.subplots_adjust(wspace =0.3,hspace=0.9) # espacio entre graficos
    
    ax[0].set_title('Reflectividad')
    ax[1].set_title(r'Ond\'icula')
    ax[2].set_title('Ruido')
    ax[3].set_title('Dato')
   
    ax[0].plot(R[:,trace],np.arange(0,nr)*dt,'k')
    ax[0].set_yticks([])
    ax[0].set_xticks([])
    ax[0].spines['top'].set_visible(False)
    ax[0].spines['right'].set_visible(False)
    ax[0].spines['bottom'].set_visible(False)
    ax[0].spines['left'].set_visible(False)
    
    ax[1].plot(wavelet,np.arange(-m,m+1)*dt,'k')
    ax[1].set_ylabel('*')
    ax[1].set_yticks([])
    ax[1].set_xticks([])
    ax[1].spines['top'].set_visible(False)
    ax[1].spines['right'].set_visible(False)
    ax[1].spines['bottom'].set_visible(False)
    ax[1].spines['left'].set_visible(False)
                                     
    
    ax[2].plot(ruido[m:dl-m,trace],np.arange(m,dl-m)*dt,'k')
    ax[2].set_ylabel('+')
    ax[2].set_yticks([])
    ax[2].set_xticks([])
    ax[2].spines['top'].set_visible(False)
    ax[2].spines['right'].set_visible(False)
    ax[2].spines['bottom'].set_visible(False)
    ax[2].spines['left'].set_visible(False)
   
    d1 = d[m:dl-m,trace]
    x1 = np.arange(m,dl-m)*dt
    ax[3].plot(d1,x1,'k')
    ax[3].fill_betweenx(x1,d1,0, where=d1 >= 0, facecolor='black')
    ax[3].set_ylabel(r'$\parallel$')
    ax[3].set_yticks([])
    ax[3].set_xticks([])
    
    ax[3].spines['top'].set_visible(False)
    ax[3].spines['right'].set_visible(False)
    ax[3].spines['bottom'].set_visible(False)
    ax[3].spines['left'].set_visible(False)
    #ax[3].set_ylim([-1, 1])
    #ax[3].set_xlim([0, xfin*dt])

#fig.savefig('sintetico.pdf', bbox_inches='tight')
plt.show()
