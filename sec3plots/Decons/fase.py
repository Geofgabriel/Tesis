"""
Script para estimar la ondícula sísmica, utilizando el método de la 
curtosis para la estimación de la fase, y realizar la deconvolución
de tipo sparse via IRLS.
by Gabriel R. Gelpi
"""

import numpy as np
import scipy.signal as ss
import scipy.stats as se
import matplotlib.pyplot as plt
from numpy.linalg import inv
from conv_mat import convmtx
from matplotlib import rcParams
rcParams['text.usetex'] = False
plt.rc('font', family='serif')#serif,monospace,sans
plt.rcParams.update({'font.size': 18})

#=============================================================================
# Datos de entrada
#=============================================================================
datopath = r'/home/gaby/Documents/Tesis_v6/Graficos_suplementarios/Deconvolucion/v6/sec3plots/Decons/dato_marmo.txt'
data = np.loadtxt(datopath)
ns,nt=data.shape # dimensiones del dato(muestras y offset).

dt = 0.002 # intervalo de muestreo

#=============================================================================
# Estimamos una ondicula de fase cero
#=============================================================================
# Espectro de amplitud de todas las trazas disponibles.
espc=np.abs(np.fft.fft(data,axis=0))
# Promedio del espectro de amplitud sobre todas las trazas disponibles.
espc=(np.sum(espc,axis=1)/nt)
# Taper en frecuencia si es necesario--> ruido de banda limitada.
LT= 80
espc[LT:ns-LT]=0.0 # observar que hay maneras mas elegantes que esta.
f = np.arange(0,ns,1)/(ns*dt)

# Vuelvo al dominio del tiempo quedandome con la parte real unicamente.
espc = np.real(np.fft.ifft(espc))
# reordeno la senal
split= np.split(espc,2)
espc = np.append(split[1],split[0])

# Aplico una ventana de Hamming y se obtiene la ondicula de fase cero.
nw=55 # Longitud de la ventana de Hamming...o long de la ondicula.
wav = np.float32(np.hamming(nw)*espc[np.int_(ns/2)-np.int16(nw/2):np.int_(ns/2)+np.int_(nw/2)+1])
wav = wav/np.max(np.abs(wav))

#np.savetxt('zero_phase_wavelet.txt',wav) # guardo la ondicula en un archivo .txt

#==================================================================================
# Estimacion de la fase mediante el método basado en la máximización de la curtosis
#==================================================================================

nang=181    # tamano de la grilla de busqueda.
minang=-90  # minimo valor de la grilla.
maxang=90   # maximo valor de la grilla.
angles = np.radians(np.linspace(minang,maxang,nang)) # genero la grilla.
#obs: usando radians y degree de numpy no hace falta escribir conversiones de grad a rad.

kt = np.zeros((nang,nt))
kt_prom = np.zeros(nang)

for i in range(0,nang):    
    for j in range(0,nt):
        
        #Rotacion del dato
        trace_rot = np.cos(angles[i])*data[:,j]-np.sin(angles[i])*np.imag(ss.hilbert(data[:,j]))
 
        #Kurtosis del dato rotado
        kt[i][j]=se.kurtosis(trace_rot,fisher=True,bias=True)            

    kt_prom[i]=np.sum(kt[i,:])/nt

kt_fase = np.degrees(-angles[np.argmax(kt_prom)])

print("Fase kt:",kt_fase) 
print('Observar que la fase con la que generamos el dato es -30')

#np.savetxt('curva_kt.txt',kt_prom) # descomentar si se quiere guardar la fase

#=====================================================================
# La ondicula estimada es:
#=====================================================================
# Combino la informacion de fase y amplitud:
#kt_fase = 45 # fase utilizada para mostrar que pasa con la decon cuando la ondícula es incorrecta
wavelet_kt = np.float32(np.cos(np.radians(kt_fase))*wav-np.sin(np.radians(kt_fase))*np.imag(ss.hilbert(wav)))


#np.savetxt('full_wavelet.txt',wavelet_kt)
#========================================================
# Deconvolucion con regularizacion L1
#========================================================
hlambda = 0.05          # Parametro de Trade-off de la decon
epsilon = 0.00001      # parametro del IRLS
it = 50                # iteracion maxima
tol = 0.001            # Tolerancia
nr=ns-nw+1             # longitud de la reflectividad estimada

X = convmtx(wavelet_kt,nr) # Matriz de convolucion

rh = np.zeros((nr,nt))
i = 0
for i in range(0,nt):
   
    I = np.eye(nr)
    h = np.dot(X.T, X)
   
    delta = 0
    rnorm_old = 0
    rnorm = 0
    j = 0
    info = False

    while j<it and info == False:
        rh[:,i]=np.dot(inv(h+hlambda*I),np.dot(X.T,data[:,i]))
        rnorm = np.sqrt(np.dot(rh[:,i],rh[:,i]))          # calculo la norma l2
        delta = np.abs(rnorm-rnorm_old)/rnorm
        info = delta<tol
        rnorm_old = rnorm
        
        I = np.diag(1.0/(np.abs(rh[:,i])+epsilon))
        j = j+1

#np.savetxt('refle_estimada.txt',rh) # descomentar para guardar la reflectividad estimada
#====================================================
# Graficos
#====================================================
dt=0.002
m=np.int16((nw-1)/2)
print('arranco a graficar')
colormap = 'gray'# seismic color option de matplotlib
refle = r'/home/gaby/Documents/Tesis_v6/Graficos_suplementarios/Deconvolucion/v6/sec3plots/Decons/R_marmousi.txt'
R = np.loadtxt(refle) # cargo la refle si la conozco

fig, ax = plt.subplots(1,3,figsize=(12,8))# base x alturad
ax[0].set_title('Dato')
ax[1].set_title(r'Deconvoluci\'on')
ax[2].set_title('Reflectividad')
ax[0].set_yticks(np.arange(0,nr+1,50))
ax[1].set_yticks(np.arange(0,nr+1,50))
ax[2].set_yticks(np.arange(0,nr+1,50))
#ax[0].set_yticklabels(np.arange(0,ns+1,50)*dt)
ax[1].set_yticklabels([])
ax[2].set_yticklabels([])
ax[0].set_ylabel('Tiempo (s)')
[ax[i].set_xlabel('Traza') for i in range(0,3)]

data = np.column_stack((data,data[:,-1])) # repito la ultima columna por cuestiones esteticas.
rh = np.column_stack((rh,rh[:,-1])) 
R = np.column_stack((R,R[:,-1]))
ax[0].imshow(data[m:-1-m,:],aspect='auto',cmap=colormap)
ax[1].imshow(rh,aspect='auto',cmap=colormap)
# flechas para mostrar diferencias entre la deconvolución "correcta" e incorrecta.
#ax[1].annotate('', xy=(60,140 ), xytext=(68,133),
#            arrowprops=dict(facecolor='red', shrink=0.05, width=10, headwidth=18,frac=0.3),
#            )
#ax[1].annotate('', xy=(60,95 ), xytext=(68,100),
#            arrowprops=dict(facecolor='red', shrink=0.05, width=10, headwidth=18,frac=0.3),
#            )
#ax[1].annotate('', xy=(60,60 ), xytext=(68,55),
#            arrowprops=dict(facecolor='red', shrink=0.05, width=10, headwidth=18,frac=0.3),
#            )
ax[2].imshow(R,aspect='auto',cmap=colormap)

#fig.savefig('fase_cor.pdf', bbox_inches='tight')
plt.show()




