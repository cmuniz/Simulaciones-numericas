# -*- coding: utf-8 -*-
import numpy as np
import numpy.linalg as la
from scipy.linalg import eig,eigh
from matplotlib import pyplot as plt


#ondas planas
#Data

Xmin=0.0
Xmax=10.0
nX=1000
Xoff=0.0
Lx=Xmax-Xmin
Xstep=Lx/(nX-1)
xx=np.zeros(nX)
for ix in range(nX):
    xx[ix]=ix*Xstep
nA=5
nEi=nA
nF=10
inF=2
zeta=-1.0
delta=1
e=1
beta=1.64
gamma=1.0
m=np.linspace(-nF/2.,nF/2.,nF+1)
print(m)

#potencial
XstepA=(Lx-2.0*Xoff)/(nA+1)
xa=np.zeros(nA)
for ia in range(nA):
    xa[ia]=Xmin+Xoff+Xstep*(ia*1.0)
vt=np.zeros(nX)
for ia in range(nA):
    for ix in range(nX):
        vt[ix]=vt[ix]+zeta/((xx[ix]-xa[ia]))**2+delta

plt.plot(xx,vt)
plt.show()

def potencial(x,at):  
    cte = 1
    z = 1
    V=np.zeros(len(x))
    for k in range(len(x)):
        for j in range(len(at)):
            Vj=-(e**2/cte)*(z/(abs(x[k]-at[j])+delta))
            V[k]=V[k]+Vj
    return V

#V=potencial(mallado,atomos)
V=potencial(xx,xa)
fig=plt.figure()
plt.plot(xx,V,'-',)
plt.show()
#funcion de onda

def fiF(x,Lx,Nf):
    nif=0.5*np.linspace(-Nf/2.,Nf/2.,Nf+1)
    real=np.zeros((len(x),len(nif)))
    imaginaria=np.zeros((len(x),len(nif)))
    for j in range(len(x)):
        for k in range(len(nif)):
            real[j][k]=np.cos(2*np.pi*nif[k]*x[j]/Lx)
            imaginaria[j][k]=np.sin(2*np.pi*nif[k]*x[j]/Lx)
    return nif,real,imaginaria

natomos,R,I=fiF(xx,Lx,nF)
print(len(R[1,:]))
print(len(xx))
for i in range(len(m)):
    print(-nF+i)
    plt.plot(xx,R[:,i],'r',xx,I[:,i],'b')
    #poner titulo a cada n
    plt.show()
    
#dividir el hamiltoniano entre la parte real y la imaginaria


def derivacion(x,F,delta):      #delta es Xstep
    #obtenemos la derivada segunda
    d1=np.zeros((len(F[:,1])-1,len(F[1,:])))        #la dimension sera en una menos que la F en la componente de las x 
    d2=np.zeros((len(d1[:,1])-1,len(d1[1,:])))          #analogo, menos dimension que d1 en la misma componente
    for k in range(1,len(x)):
        d1[k-1,:]=(F[k,:]-F[k-1,:])/delta
    for j in range(1,len(x)-1):
        d2[j-1,:]=(d1[j,:]-d1[j-1,:])/delta
    return d2
derivreal=derivacion(xx,R,Xstep) 
derivim=derivacion(xx,I,Xstep)       

#hamiltoniano
def hamiltoniano(V,re,im):
    Hr=np.zeros((len(re[:,1]),len(re[1,:])))
    Hi=np.zeros((len(im[:,1]),len(im[1,:])))
    for i in range(0,len(re[:,1])):
        for j in range(0,len(re[1,:])):
            Hr[i,j]=-re[i,j]+V[i]       #MULTIPLICAR POR LAS FUNCIONES DE ONDA
            Hi[i,j]=-im[i,j]+V[i]
    return Hr,Hi
print('HAMILTONIANOS')
Hamilre,Hamilim=hamiltoniano(V,derivreal,derivim)
for l in range(nF+1):
    print(-nF/2.+l)
    plt.plot(xx[0:998],Hamilre[:,l],'r-',xx[0:998],Hamilim[:,l],'b-')
    plt.show()
#funcion de onda total, la suma de las componentes por su constante, aplicacion del hamiltoniano 8derivar)

def funcionda(x,fireal,fimaginaria,nif):        #nif natomos
    fondareal=np.zeros(len(x))
    fondaimaginaria=np.zeros(len(x))
    sumafireal=0
    sumafimaginaria=0

    for k in range(len(x)):
        for j in range(len(nif)):
            sumafireal=sumafireal+fireal[k][j]
            sumafimaginaria=sumafimaginaria+fimaginaria[k][j]
        fondareal[k]=sumafireal
        fondaimaginaria[k]=sumafimaginaria
        sumafireal=0
        sumafimaginaria=0
        
    return fondareal,fondaimaginaria

FHIR,FHII=funcionda(xx,R,I,natomos)

print(len(FHIR))
plt.plot(xx,FHIR,'r-',xx,FHII,'b-')
plt.title('Funcion de onda')
plt.show()

#fi segunda va desde 1 hasta Nx-2, fi desde 2 hasta Nx-1, y V desde 2 hasta Nx-1    
#normalizar, integrar con trapecio. dividimos la funcion de onda entre la raiz de la integral
#usar  Eval,Evec=eig(Hmat,Smat), da los autovalores desordenados, hemos de ordenarlos y normalizarlos.
#ind=np.argsort(np.real(Eval))
#Eval=np.real(Eval[ind])
#Evec=Evec[:,ind]
#for iV in range(nF):
    #norm2=0.0
    #for jF in range(nF):
    #   for kF in range(nF):
    #       norm2=norm2+np.conj(Evec[kF,iV])*Smat[jF,kF]*Evec[jF,iV]
    #for jF in range(nF):
    # Evec[jF,iV]=Evec[jF,iV]/np.sqrt(norm2)
