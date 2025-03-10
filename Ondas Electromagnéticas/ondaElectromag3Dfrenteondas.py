# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
# Definicion de los parametros

c = 3.0e8                       #Velocidad de la luz
wl = 632.9e-9                   #longitud de onda
tp = wl/c                       # Periodo
Nx = 60
Nz = 60                        # Puntos espaciales
Nt = 500                        # Puntos temporales
k0x = 30                        # Posicion de la fuente
k0z = 30
dz = wl/10                      # Particion espacial
dt = dz/(2*c)                   # Particion temporal
t = 0.0                         # Valor del tiempo inicial
z = dz*np.linspace(1,Nz,Nz)     # Valores z del espacio real
x = dz*np.linspace(1,Nx,Nx)
X, Z = np.meshgrid(x,z)
zi=dz/1e-6                      # Valor de z inicial (micrometros)(microns)
zf=Nz*dz/1e-6                   # Valor final de z (microns)
t0 = 2.5*tp                     #Time lag Tiempo de desfase
sigma = 0.5*tp                  #spread #Dispersion
ta = 4.0                        # Tiempo de aceleracion
pldim="3d"                       # 2D or 3D plot

# Inicializacion de los campos
ey=np.zeros((Nx,Nz))                 # Valores iniciales del campo electrico
hx=np.zeros((Nx,Nz))
hz=np.zeros((Nx,Nz))                 # Valores iniciales del campo magnetico
ze=np.zeros((Nx,Nz))                 # Array with zeros tu use in 3D plot

# Plot definitions
fig = plt.figure(figsize=(12,6))
ax = fig.add_subplot(111,projection="3d")


ei1x = ei2x = np.zeros(Nx)
ef1x = ef2x = np.zeros(Nx)
ei1z = ei2z = np.zeros(Nz)
ef1z = ef2z = np.zeros(Nz)
for n in range(1,Nt+1):
    
    ei = ey[1]
    ef = ey[Nz-2]
    # Campo electrico
    for k in range(1,Nz-1):
        for l in range(1,Nx-1):
            # Obstaculo
            if k > 45 and k < 50 and l > 45 and l < 50:
                ey[l,k] = 0
            else:
                ey[l,k]=ey[l,k]+0.5*(hx[l,k]-hx[l,k-1]-hz[l,k] + hz[l-1,k])
                
            
                   
    
    # Fuente
    ey[k0x-3:k0x+3,k0z] = ey[k0x,k0z] + np.exp(-0.5*((t-t0)/sigma)**2.0) # blanda
#    ey[k0x-3:k0x+3,k0z] = np.exp(-0.5*((t-t0)/sigma)**2.0) # dura
    
    # Conciciones de contorno absorbentes
    ey[0:Nx,0], ei1x, ei2x = ei1x, ei2x, ey[0:Nx,1]
    ey[0:Nx,-1], ef1x, ef2x = ef1x, ef2x, ey[0:Nx,-2]
    
    ey[0,0:Nz], ei1z, ei2z = ei1z, ei2z, ey[1,0:Nz]
    ey[-1,0:Nz], ei1z, ei2z = ei1z, ei2z, ey[-2,0:Nz]
    
#    ey[:,0], ei1x, ei2x = ei1x, ei2x, ey[:,1]
#    ey[:,-1], ef1x, ef2x = ef1x, ef2x, ey[:,-2]
#    
#    ey[0,:], ei1z, ei2z = ei1z, ei2z, ey[1,:]
#    ey[-1,:], ei1z, ei2z = ei1z, ei2z, ey[-2,:]

 
    #Campo magnetico
    for k in range(0,Nz-1):
        for l in range(0,Nx-1):
            hx[l,k]=hx[l,k]+0.5*(ey[l,k+1]-ey[l,k])
            hz[l,k]=hz[l,k]+0.5*(ey[l,k]-ey[l+1,k])
        
    if np.ceil(n/ta)==n/ta: 
        ax.cla()
        ax.set_zlim([-2.1,2.1])
#        ax.set_xlim()
        ax.set_xlabel("x (microns)")
        ax.set_ylabel("z (microns)")
        ax.set_zlabel("Ey")
        ax.set_title("n = %i time steps" %n, loc= "left")


#        ax1.cla()
#        ax2.cla()
        ax.plot_surface(X/1e-6,Z/1e-6,ey)
#        ax1.plot_surface(X/1e-6,Z/1e-6,hx)
#        ax2.plot_surface(X/1e-6,Z/1e-6,hz)
        plt.pause(0.01)
    
    # Time
    t = n*dt
