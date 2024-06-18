import numpy as np
import scipy.signal as ss
import scipy.stats as se
import matplotlib.pyplot as plt
from numpy.linalg import inv
from conv_mat import convmtx
from matplotlib import rcParams
from rickerpy import rickerpy

nr = 300                                             # Nro de muestras de la reflectividad
ntr = 1                                             # numero de reflectividades

R = np.zeros((nr,ntr)) 

R[47,:] = -0.02
R[65,:] = 0.2
R[82,:] = 0.22
R[121,:] = 0.05
R[150,:] = -0.09
R[172,:] = -0.1
R[205,:] = 0.08
R[247,:] = 0.18
R[265,:] = 0.05

print('long refle',R.shape)

dt = 0.002 # intervalo de muestreo
l = 108
f0 = 30

o = rickerpy(l,2,f0)
#o = np.reshape(-1,1)
print('long ondicula',len(o))
print('tipo ondi',o.shape)

dl = nr+o.shape[0]-1 # nro de elementos que va a tener las trazas generadas.
d = np.zeros((dl,ntr))
d = np.convolve(R[:,0],o[:,0],mode='full')

SNR = 20
ruido = np.random.randn(dl)
ruido = ruido/np.amax(np.abs(ruido))
alfa = np.max(np.abs(d))/(np.max(np.abs(ruido))*SNR)            
d = d + alfa*ruido

# DECON l2
hmu = 0.0005
X = convmtx(o[:,0],nr)    # Matriz de convolucion
rh = np.zeros((nr,dl))
I = np.eye(nr)
h = np.dot(X.T, X)

r_hat = np.dot(inv(h+hmu*I),np.dot(X.T,d))

# DECON l1
hlambda = 0.0001          # Parametro de Trade-off de la decon
epsilon = 0.00001         # parametro del IRLS con la norma l1
sc = 0.0001               # parametro del IRLS con la norma de cauchy                        
it = 60                   # iteracion maxima
tol = 0.0001              # Tolerancia

X = convmtx(o[:,0],nr)    # Matriz de convolucion
rh = np.zeros((nr,dl))
I = np.eye(nr)
h = np.dot(X.T, X)
   
delta = 0
rnorm_old = 0
rnorm = 0
j = 0
info = False

while j<it and info == False:
        rh = np.dot(inv(h+hlambda*I),np.dot(X.T,d))
        rnorm = np.sqrt(np.dot(rh,rh))          # calculo la norma l2
        delta = np.abs(rnorm-rnorm_old)/rnorm
        info = delta<tol
        rnorm_old = rnorm
        
        #I = np.diag(1.0/(np.abs(rh)+epsilon)) #1/(abs(r)+eps)
        I = np.diag(1.0/(rh**2+sc**2))         # 1/(sigma**2+r**2)
        j = j+1


# guardamos los archivos:
plots = 2
if plots == 1:
    #np.savetxt('refle_.txt',R)
    #np.savetxt('refle_cauchy.txt',rh)
    #np.savetxt('refle_l2.txt',r_hat)
    #np.savetxt('data.txt',d)
    print('hay graficos para guardar')

#plt.plot(R)
#plt.show()
#plt.plot(o)
#plt.show()
plt.plot(d[26:-25],'r')
plt.plot(r_hat,'b')
plt.show()
#plt.plot(rh)
#plt.show()
#plt.plot(R,'r')
#plt.plot(rh,'b')

plt.show()