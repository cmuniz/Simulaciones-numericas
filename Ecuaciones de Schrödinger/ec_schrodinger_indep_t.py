# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
import numpy as np
import numpy.linalg as la
from scipy.linalg import eig,eigh

# Datos

Na = 5              # numero de atomos
Nx = 100           # numero de subintervalos
Nf = Na              # numero de funciones
xi = 0              # posicion inicial
xf = 10             # posicion final
Lx = xf-xi
dx = (xf - xi)/Nx   # incremento de x
z = 1               # numero de protones
delta = 0.1         #0.01
beta = 1.63
gamma = 1.0
x = np.zeros(Nx)

ibF = 2             # Selector base: 0-Slayer, 1-Gauss, 2-Ondas planas

for i in range(Nx):
    x[i] = xi + i*dx
#x = np.linspace(xi, xf, Nx)
xa = np.linspace(xi,xf, Na + 2)[1:-1]

# Funciones

Va = lambda x, xa : -z/(np.abs(x-xa)**2 + delta)

def V(x,xa):
    """Potencial"""
    v=np.zeros_like(x)
    for ia in xa:
        v += Va(x,ia)
    return v

def STO(x,xa):
    """Slayer"""
    phi = np.zeros((Nf,Nx),complex)
    for i in range(Na):
        phi[i] = np.exp(-beta*abs(x-xa[i]))
    return phi
    
def GTO(x,xa):
    """Gauss"""
    phi = np.zeros((Nf,Nx),complex)
    for i in range(Na):
        phi[i] = np.exp(-beta*abs(x-xa[i])**2)
    return phi

def ondasPlanas(x,xa):
    """ondas planas"""
    phi = np.zeros((Nf,Nx),complex)
    for iF in range(1,Nf):
        G = 2 * np.pi / Lx * (int(-Nf/2) +iF-1) #redondear nF/2
        for ix in range(Nx):
            phi[iF,ix]=np.exp(1j*G*(x[ix]))
    return phi

            
def integral(X,fx):    
    #f es un vector formado por f(x) para cada x en X, no una funcion
    I = 0
    for i in range(len(X)-1): #-2 porque queremos que i+1 sea el ultimo
        I = I + dx*(fx[i+1]+fx[i])/2
    return I

def derivada(X,fx):
    #fx y X son vectores, deben tener la misma dimension
    
    dfx = []
    #el primer termino del vector dfx lo hacemos con los dos siguientes
    dfx.append((fx[1]-fx[0])/dx) 
    for i in range(1,len(fx)-1):     
        dfxanterior = (fx[i]-fx[i-1])/dx
        dfxsiguiente = (fx[i+1]-fx[i])/dx
        dfx.append((dfxanterior + dfxsiguiente)/2)
    #el ultimo termino del vector dfx lo hacemos con los dos anteriores
    dfx.append((fx[-1]-fx[-2])/dx)
    return dfx

def derivadasegunda(X,fx):
    #fx y X son vectores, deben tener la misma dimension
    
    ddfx = []
    #el primer termino del vector dfx lo hacemos con los dos siguientes
    ddfx.append(0) 
    for i in range(1,len(fx)-1):     
        ddfx.append((fx[i+1]+fx[i-1]-2*fx[i])/(dx*dx))
    #el ultimo termino del vector dfx lo hacemos con los dos anteriores
    ddfx.append(0)
    return ddfx


# Calcular y pintar potencial
v = V(x,xa)
plt.plot(x,v)
plt.plot(xa,np.zeros_like(xa),"o")
plt.plot(x,np.zeros_like(x),"--")
plt.show()

# obtener y Pintar base

if ibF == 0:
    phi = STO(x,xa)
elif ibF == 1:
    phi = GTO(x,xa)
elif ibF == 2:
    phi = ondasPlanas(x,xa)
else:
    phi = np.zeros((Nf,Nx),complex)
    
#NORMALIZACION
for j in phi: 
    fx = np.array(j)
    norma = integral(x,fx.conjugate()*fx)
    for i in range(len(j)):
        j[i] = j[i] / np.sqrt(norma)


#REPRESENTAMOS LA PARTE REAL
for i in range(Na):
    plt.plot(x,np.real(phi[i]))
plt.title('Parte real de las funciones de la base')
plt.show()

#REPRESENTAMOS LA PARTE IMAGINARIA
for i in range(Nf):    
    plt.plot(x,np.imag(phi[i]))
plt.title('Parte imaginaria de las funciones de la base') 
plt.show()


"""
solape y hamiltoniano
"""
S = np.zeros((Nf,Nf))
H = np.zeros((Nf,Nf))

for iF in range(Nf):
    for jF in range(Nf):
        for ix in range(1,Nf-1):
            S[iF,jF] = S[iF,jF] + np.conj(phi[iF,ix])*phi[iF,ix]*dx
            # falta el hamiltoniano


#H=np.zeros((nfunciones,nfunciones))
#for j in range(len(BASE)):
#    for i in range(len(BASE)):    
#        ddFj = np.array(derivadasegunda(X,BASE[j]))   
#        Vx = np.array(V(X)) 
#        Fj = np.array(BASE[j])
#        HFj = -0.5*ddFj+Vx*Fj   
#        Fi = np.array(BASE[i])
#        Fi_HFj = Fi.conjugate()*HFj
#        H[i][j] = integral(X,Fi_HFj)
#print('\n'+'H=')
#print(H)
#
#
#S = np.zeros((nfunciones,nfunciones))
#for i in range(len(BASE)):
#    for j in range(len(BASE)):
#        Fj=np.array(BASE[j])
#        Fi=np.array(BASE[i])
#        Fj_Fi= Fj.conjugate()*Fi
#        S[j][i]=integral(X,Fj_Fi)
#
#print('\n'+'S=')
#print(S)

#Eval, Evec = eig(H,S) #Autovalores = E y autovectoes
#
#"""
#normalizar autovector
#"""
#ind = np.argsort(np.real(Eval)) #ordena autovalores
#Eval = np.real(Eval[ind])
#Evec = Evec[:,ind] # ordena autvectores igual que autovalores
#for iV in range(nF):
#    norm2 = 0.0 #calcular norma
#    for jF in range(nF):
#        for kF in range(nF):
#            norm2=norm2+np.conj(Evec[kF,iV])*np.conj(S[jF,kF])*Evec[jF,iV]
#        for jF in range(nF):
#            Evec[jF,iV]=Evec[jF,iV]/np.sqrt(norm2) #normaliza el vector dividiendo por la norma
