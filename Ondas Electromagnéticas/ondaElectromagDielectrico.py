# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
# Definicion de los parametros

c = 3.0e8                       #Velocidad de la luz
wl = 632.9e-9                   #longitud de onda
tp = wl/c                       # Periodo
Nz = 200                        # Puntos espaciales
Nt = 700                        # Puntos temporales
k0 = 50                         # Posicion de la fuente
kd = 100                        # Posicion del dilectrico
dz = wl/10                      # Particion espacial
dt = dz/(2*c)                   # Particion temporal
t = 0.0                         # Valor del tiempo inicial
z = dz*np.linspace(1,Nz,Nz)     # Valores z del espacio real
zi=dz/1e-6                      # Valor de z inicial (micrometros)(microns)
zf=Nz*dz/1e-6                   # Valor final de z (microns)
t0 = 2.5*tp                     #Time lag Tiempo de desfase
sigma = 0.5*tp                  #spread #Dispersion
ta = 4.0                        # Tiempo de aceleracion
pldim = "2d"                    # 2D or 3D plot
Er = 4                          # Permitividad relativa


def cd(k):
    if k < kd:
        return 0.5
    else:
        return 0.5 / Er

# Inicializacion de los campos
ex=np.zeros(Nz)                 # Valores iniciales del campo electrico
hy=np.zeros(Nz)                 # Valores iniciales del campo magnetico
ze=np.zeros(Nz)                 # Array with zeros tu use in 3D plot

# Plot definitions
fig = plt.figure(figsize=(8,6))
if pldim == "2d":
    ax = fig.add_subplot(111)
else:
    ax = fig.add_subplot(111,projection="3d")

# Variables auxiliares conciones de contorno
ei1 = ei2 = 0
ef1 = ef2 = ef3 = ef4 = 0

for n in range(1,Nt+1):

    # Campo electrico
    for k in range(1,Nz-1):

        ex[k]=ex[k]+cd(k)*(hy[k-1]-hy[k])

    # Fuente
    ex[k0] = ex[k0] + np.exp(-0.5*((t-t0)/sigma)**2.0)
    
    # Conciciones de contorno absorbentes
    ex[0], ei1, ei2 = ei1, ei2, ex[1]
    ex[-1], ef1, ef2, ef3, ef4 = ef1, ef2, ef3, ef4, ex[-2]

    
    #Campo magnetico
    for k in range(0,Nz-1):
        hy[k]=hy[k]+0.5*(ex[k]-ex[k+1])
        
    # Representación
    if np.ceil(n/ta)==n/ta:
        ax.cla()
        ax.set_xlim([zi,zf])
        ax.set_ylim([-2.1,2.1])
        ax.set_xlabel("z (microns)")
        ax.axvline(x=dz*kd/1e-6)
        if pldim != "3d":
            ax.set_ylabel("Ex(z,t), Hy(z,t)")
            ax.set_title("n = %i time steps" %n)
            ax.plot(z/1e-6,ex)
            ax.plot(z/1e-6,hy)
        else: # 3D
            ax.set_zlim([-2.1,2.1])
            ax.set_ylabel("Ex(z,t)")
            ax.set_zlabel("Hy(z,t)")
            ax.set_title("n = %i time steps" %n, loc= "left")
            ax.view_init(45,225)
            ax.plot(z/1e-6,ex,0)
            ax.plot(z/1e-6,ze,hy)
        plt.pause(0.01)
    
    # Time
    t = n*dt
