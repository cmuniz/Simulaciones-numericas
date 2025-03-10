# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import eig

'''DEFINICIONES DE FUNCIONES QUE NOS SERAN UTILES MAS ADELANTE'''
#SON FUNCIONES QUE CALCULAN INTEGRALES Y DERIVADAS NUMERICAS, ES DECIR
  #LOS DATOS QUE USAN SON LISTAS DE VALORES, NO CALCULAN DERIVs E INTEGs SIMBOLICAS
def integral(X,fx):    
    #f es un vector formado por f(x) para cada x en X, no una funcion
    dx = X[1]-X[0]
    I = 0
    for i in range(len(X)-1): #-2 porque queremos que i+1 sea el ultimo
        I = I + dx*(fx[i+1]+fx[i])/2
    return I

def derivada(X,fx):
    #fx y X son vectores, deben tener la misma dimension
    dx = X[1]-X[0]
    
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
    dx = X[1]-X[0]
    
    ddfx = []
    #el primer termino del vector dfx lo hacemos con los dos siguientes
    ddfx.append(0) 
    for i in range(1,len(fx)-1):     
        ddfx.append((fx[i+1]+fx[i-1]-2*fx[i])/(dx*dx))
    #el ultimo termino del vector dfx lo hacemos con los dos anteriores
    ddfx.append(0)
    return ddfx

'''DATOS, variables'''
###############################################################################
#Todas las distancias vienen dadas en 1 Bohr = 0,529 Amstons                  #
me = 1                                                                        #
hb = 1 #h barra                                                               #
e = 1                                                                         #
#Al tomar estos valores, el operador T=-1/2(hb/e*me)*F" (Energia Cinetica)    #
  #vendra dado por -1/2 F"                                                    #
mallado = 1000                                                                #
natomos = 5 #numero de atomos                                                 #
Z = 1                                                                         #
                                                                              #
Xmax = 10  #longitud de la caja                                               # 
Xatomos = np.linspace(1+2/3,8+1/3,natomos)                                    #
                                                                              #
'''LA BASE QUE NOSOTROS VAYAMOS A UTILIZAR LA DEFINIMOS AQUI'''               #
tipobase = 'base ondas planas'                                                            #
###############################################################################

#Definimos el potencial que tiene un atomo debido a los demas
dd = 0.1 #arreglo para que el potencial no diverja
def V(x):
    Suma = 0
    for i in range(natomos):
        Suma = Suma-Z/((x-Xatomos[i])**2+dd)**0.5
    return Suma


nfunciones = natomos #numero de funciones de la base = numero de atomos

X = np.linspace(0,Xmax,mallado)        
representacion = plt.figure()
eje = representacion.add_subplot(111)

eje.plot(X,V(X),'r-')
plt.title('Potencial atC3mico')
plt.show()



'''DEFINICION DE DISTINTAS BASES'''
B=1.0

#Slater
Si = np.zeros((nfunciones,mallado),complex)
for i in range(len(Xatomos)):
    Si[i] = np.exp(-B*abs(X-Xatomos[i]))
    
#Gauss
Gi = np.zeros((nfunciones,mallado),complex)
for i in range(len(Xatomos)):
    Gi[i] = np.exp(-B*(X-Xatomos[i])**2)

#Base ondas planas
BOP = np.zeros((nfunciones,mallado),complex)
mi = range(int(-(nfunciones/2-0.5)),int(nfunciones/2-0.5+1),1)
for i in range(len(Xatomos)):
    #L = X[-1]-X[0]
    L=10
    ki = 2*3.1416*mi[i]/L
    BOP[i] = np.exp(1j*ki*X)  

baseencontrada = 0
if tipobase=='slater':
    BASE = Si
    baseencontrada = 1
if tipobase=='gauss':
    BASE = Gi
    baseencontrada = 1
if tipobase=='base ondas planas':
    BASE = BOP
    baseencontrada = 1
if baseencontrada ==0:
    print('')
    print('No se ha encontrado la base, por favor, seleccione una base distinta')
    print('')

#NORMALIZACION: normalizamos la base
for j in BASE: 
    #PRINTEAMOS LA NORMA DE LA FUNCION SIN NORMALIZAR
    fx = np.array(j)
    norma = integral(X,fx.conjugate()*fx)
    print('Sin normalizar:'+str(norma))
    for i in range(len(j)):
        j[i] = j[i]/((norma)**0.5)

    #comprobamos que este normalizada:
    #PRINTEAMOS LA NORMA DE LA FUNCION UNA VEZ NORMALIZADA
    J = np.array(j) 
    print('Normalizada:'+str(integral(X,J.conjugate()*J)))

#REPRESENTAMOS LA PARTE REAL
representacion2 = plt.figure()    
partereal = representacion2.add_subplot(111)
for i in range(natomos):
    partereal.plot(X,BASE[i])
plt.title('Parte real de las funciones de la base')
plt.xlabel(tipobase)
plt.show()

#REPRESENTAMOS LA PARTE IMAGINARIA
representacion3 = plt.figure()
parteimaginaria = representacion3.add_subplot(111)    
for i in range(natomos):    
    parteimaginaria.plot(X,np.imag(BASE[i]))

plt.title('Parte imaginaria de las funciones de la base') 
plt.xlabel(tipobase)
plt.show()



'''
Ahora vamos a ver como actua el hamiltonmiano H sobre los vectores de la base
H = T + V  donde T = hbarra/2m * d^2/dx^2 que en unidades atomicas sera:
    T = 1/2 * d^2/dx^2
para hacer la segunda derivada hay que hacerlo con uno de los metodos numericos
vistos en el primer cuatrimestre
'''

#DERIVADA NUMERICA PARA CUANDO QUEREMOS DERIVAR LOS VALORES EN UNA CISTA



#REPRESENTAMOS LA segunda DERIVADA
representacion4 = plt.figure()
TRep = representacion4.add_subplot(111)    
for i in range(natomos):  
    ddfx = derivadasegunda(X,BASE[i])
    ddfx = np.array(ddfx)
    T = ddfx/2 # T = 1/2 * d^2/dx^2
    X = np.array(X)
    T = np.array(T)
    TRep.plot(X[1:len(X)-1],T[1:len(X)-1]) 

plt.title('Segunda derivada de las funciones (T)')
plt.xlabel(tipobase)    
plt.show()




'''  TEORIA SOBRE EL HAMILTONIANO COMO OPERADOR CUANTICO:
    
Hij = <i|H|j>=integral(dx Fi* H Fj)
donde por teoria sabemos que H*Fi funciona como el siguiente operador:
                     
                   H*|Fi> = -0.5Fi"+V*Fi

Queremos hallar la representacion matricial de H conociendo cC3mo actua sobre los Fis

Ademas, calcularemos tambien la matriz S que opera del siguiente modo
representacion matricial:   Sji = <j|i>
'''

H=np.zeros((nfunciones,nfunciones))
for j in range(len(BASE)):
    for i in range(len(BASE)):    
        ddFj = np.array(derivadasegunda(X,BASE[j]))   
        Vx = np.array(V(X)) 
        Fj = np.array(BASE[j])
        HFj = -0.5*ddFj+Vx*Fj   
        Fi = np.array(BASE[i])
        Fi_HFj = Fi.conjugate()*HFj
        H[i][j] = integral(X,Fi_HFj)
print('\n'+'H=')
print(H)


S = np.zeros((nfunciones,nfunciones))
for i in range(len(BASE)):
    for j in range(len(BASE)):
        Fj=np.array(BASE[j])
        Fi=np.array(BASE[i])
        Fj_Fi= Fj.conjugate()*Fi
        S[j][i]=integral(X,Fj_Fi)

print('\n'+'S=')
print(S)


evalor, evector = eig(H,S)  #Nos da los autovalores y los vectores de la matriz H

'''ESTO ESTA BIEN PERO ANTES TENEMOS QUE ORDENAR LOS VECTORES Y AUTOVALORES'''
Hdiag = np.zeros((len(evalor),len(evalor))) #diagnoalizamos H 
for i in range(len(Hdiag)):
    Hdiag[i][i] = evalor[i]
print('\n'+'Hdiag = ')
print(Hdiag)























