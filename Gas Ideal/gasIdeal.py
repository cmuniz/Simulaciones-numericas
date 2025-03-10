# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt

# Definicion de las constates del sistema
boxsize = 15    # Lado de la caja 2D
m = 1           # Masa de la particula
Nparticulas = np.arange(10, 500, 50)     # Numero de particulas
#v0 = 1         # velocidad inicial
R = 0.25 /2     # Radio de las particulas
KB = 0.01       # Constante de Boltzman
T = 100         # Temperatura
beans = 10
# Funciones
def choque(r12,v1,v2):
    mr12 = np.sqrt(r12[0]**2 + r12[1]**2)                       # Modulo de r12
    uParalelo = r12 / mr12                                      # vector unitario paralelo
    uPerpendicular = np.array([uParalelo[1],-uParalelo[0]])     # vector unitario perpendicular
    
    v1Paralelo = np.dot(v1,uParalelo)
    v1Perpendicular = np.dot(v1,uPerpendicular)
    
    v2Paralelo = np.dot(v2,uParalelo)
    v2Perpendicular = np.dot(v2,uPerpendicular)
    
    if v1Paralelo - v2Paralelo < 0:
        return np.array([v1,v2])
    
    v1_ = v2Paralelo * uParalelo + v1Perpendicular * uPerpendicular
    v2_ = v1Paralelo * uParalelo + v2Perpendicular * uPerpendicular
    return np.array([v1_, v2_])

# Simulacion temporal

deltat = 0.01
nsteps = 10000
npause = 10
tpause = 0.01



# Simulación para distintos numero de partículas
gIdeal = []
gvdW = []
for Npart in Nparticulas:
    
    print("Número de partículas: "+ str(Npart))
    
    # Definicion de las condiciones iniciales del sistema
    
    E = np.random.exponential(KB*T, Npart)
    r = np.random.rand(2,Npart)*boxsize
    tetha0 = 2 * np.pi * np.random.rand(Npart)
    v0 = np.sqrt(2/m*E)
    v = np.array([v0 * np.cos(tetha0), v0 *np.sin(tetha0)])
    
    fig = plt.figure(figsize=(12,6))
    ax = fig.add_subplot(121)
    ax1 = fig.add_subplot(122)
    ax.plot(r[0],r[1],"ro")
    ax.set_xlim([0,boxsize])
    ax.set_ylim([0,boxsize])
    E = 1/2 * m * (v[0]**2 + v[1]**2)
    ax1.cla()
    ax1.hist(E,beans,(0,2*np.mean(E)))
    T = np.mean(E)/KB
    print("Temperatura inicial: " + str(T))
    presion = []
    for istep in range(nsteps):
        #Actualizar posicion
        r = r + v * deltat
        fuerza = 0
        mDeltaV = 0
        for i in range(Npart):
            mDeltaV = 0
            # Comprobar colisiones con las paredes
            if (r[0,i] > boxsize) or r[0,i] < 0:
                v[0,i] = -v[0,i]
                
                mDeltaV = np.abs(2*v[0,i])
            if (r[1,i] > boxsize) or r[1,i] < 0:
                v[1,i] = -v[1,i]
                
                mDeltaV = np.abs(2*v[1,i])
            
            fuerza += mDeltaV/deltat
            
            for j in range(Npart-i):
                if i == j+i: 
                    continue
                r12 = r[:,j+i]-r[:,i]
                mr12 = np.sqrt(r12[0]**2 + r12[1]**2)  
                if mr12 < 2*R:
                    v1, v2 = choque(r12,v[:,i],v[:,j+i])
                    v[:,i] = v1
                    v[:,j+i] = v2
        presion.append(fuerza/(4*boxsize))        
        # Actualizar la representacion si procede
        if(istep % npause == 0):
            ax.cla()
            ax.plot(r[0],r[1], "ro")
            ax.set_xlim([0,boxsize])
            ax.set_ylim([0,boxsize])
            E = 1/2 * m * (v[0]**2 + v[1]**2)
            T = np.mean(E)/KB
    
            ax1.cla()
            ax1.hist(E,beans,(0,2*np.mean(E)))
            plt.pause(tpause)
            
    E = 1/2 * m * (v[0]**2 + v[1]**2)
    ax1.cla()
    ax1.hist(E,beans,(0,2*np.mean(E)))
    T = np.mean(E)/KB
    print("Temperatura final: " + str(T))
   
    
    A = boxsize**2
    pMedia = np.mean(presion)
    
    gI = pMedia*A/(Npart*KB*T) # gas ideal
    gIdeal.append(gI)
    print("Gas ideal: " + str(gI))
    # PA = NKbT
    #PA /NKBT = 1 # gas ideal
    #P(A-N·a)/NKBT = 1 gas de van der Waals
    a = np.pi*R**2
    gW = pMedia*(A-a*Npart)/(Npart*KB*T)# gas de van der Waals
    gvdW.append(gW)
    print("Gas Van der Waals: " + str(gW))
    print()
    
# Tabla g. Ideal vs g. van der Waals
print("{:^10}  {:^10}  {:^10}".format("N. part.", "gIdeal", "gvdW"))
print("-"*35)
for i in range(len(Nparticulas)):
    print("{:^10.0f}  {:^10.2f}  {:^10.2f}".format(Nparticulas[i], gIdeal[i], gvdW[i]))


#Representacion Gas ideal y de van der Waals vs Numero particulas
plt.figure(2)
plt.title("gIdeal y gvdW vs N. partículas")
plt.grid()
plt.plot(Nparticulas,gIdeal, "o-", label="gIdeal")
plt.plot(Nparticulas,gvdW, "o-",  label="gvdW")
plt.legend()

