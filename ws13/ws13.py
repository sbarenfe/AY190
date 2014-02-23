#!/usr/bin/env python

import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d as mpl3d

# global constants
G = 6.67e-8
msun  = 1.99e33
seconds_per_year = 24.*3600*365 # roughly
cm_per_pc = 3.1e18
distance_to_sgrAstar = 8e3 * cm_per_pc
cm_per_AU=1.49597871e13

# system parameters
initial_data_file = "sgrAstar.asc"

#units for Sun-Earth:
#distance_unit_to_cm = 1.
#time_unit_to_s = 1.
#mass_unit_to_g = 1.

#units for the Galactic center
distance_unit_to_cm = 8e3*cm_per_AU
time_unit_to_s = seconds_per_year
mass_unit_to_g = msun

Nsteps = 10000
t0 = 0
t1 = 100 * seconds_per_year
dt = (t1-t0)/Nsteps

final_data_file = "final_positions2.asc"

def NbodyRHS(u,mass,time):
    RHS=np.zeros(np.shape(u))
#Positions updated based on velocity
    RHS[:,0]=u[:,3]
    RHS[:,1]=u[:,4]
    RHS[:,2]=u[:,5]
    
    RHS[:,3]=0
    RHS[:,4]=0
    RHS[:,5]=0
#Velocities updated base on sum of forces
    for i in range(len(u[:,3])):
	for j in range(len(u[:,3])):
	    if j==i:
		continue
	    RHS[i][3]+=-G*mass[j]/((u[i][0]-u[j][0])**2+(u[i][1]-u[j][1])**2+(u[i][2]-u[j][2])**2)*(u[i][0]-u[j][0])/((u[i][0]-u[j][0])**2+(u[i][1]-u[j][1])**2+(u[i][2]-u[j][2])**2)**0.5	    
	    RHS[i][4]+=-G*mass[j]/((u[i][0]-u[j][0])**2+(u[i][1]-u[j][1])**2+(u[i][2]-u[j][2])**2)*(u[i][1]-u[j][1])/((u[i][0]-u[j][0])**2+(u[i][1]-u[j][1])**2+(u[i][2]-u[j][2])**2)**0.5	    
	    RHS[i][5]+=-G*mass[j]/((u[i][0]-u[j][0])**2+(u[i][1]-u[j][1])**2+(u[i][2]-u[j][2])**2)*(u[i][2]-u[j][2])/((u[i][0]-u[j][0])**2+(u[i][1]-u[j][1])**2+(u[i][2]-u[j][2])**2)**0.5	    
   
    return RHS

#RK2 Integrator
def NbodyRK2(u,mass,time,dt):
    k1=dt*NbodyRHS(u,mass,time)
    unew = u + dt*NbodyRHS(u+0.5*k1,mass,time+0.5*dt)
    return unew

#Total energy calculation
def TotalEnergy(u,mass,time):
    Ekin=0
    Epot=0
    for i in range(len(u[:,0])):
	Ekin+=0.5*mass[i]*(u[i][3]**2+u[i][4]**2+u[i][5]**2)
	for j in range(len(u[:,0])):
	    if i==j:
		continue
	    Epot+=-G*mass[i]*mass[j]/((u[i][0]-u[j][0])**2+(u[i][1]-u[j][1])**2+(u[i][2]-u[j][2])**2)**0.5
    Etot=Ekin+Epot
    return Etot

# main program
plt.ion()

(x,y,z,vx,vy,vz,mass) = np.loadtxt(initial_data_file, unpack = True)


# convert from unitis in initial data file to cgs
x *= distance_unit_to_cm
y *= distance_unit_to_cm
z *= distance_unit_to_cm
vx *= distance_unit_to_cm / time_unit_to_s
vy *= distance_unit_to_cm / time_unit_to_s
vz *= distance_unit_to_cm / time_unit_to_s
mass *= mass_unit_to_g

xmin = np.amin(x)
xmax = np.amax(x)
ymin = np.amin(y)
ymax = np.amax(y)
zmin = np.amin(z)
zmax = np.amax(z)
rmax = 2.5*max(abs(xmin),abs(xmax),abs(ymin),abs(ymax),abs(zmin),abs(zmax))

# use a single state vector to simplify the ODE code
# indices:
# u[:,0] = x
# u[:,1] = y
# u[:,2] = z
# u[:,3] = vx
# u[:,4] = vy
# u[:,5] = vz
u = np.array((x,y,z,vx,vy,vz)).transpose()
r_orb=np.zeros(Nsteps)
E=np.zeros(Nsteps)

t=np.zeros(Nsteps)
usave=np.zeros((np.shape(u)[0],np.shape(u)[1],Nsteps))
for it in range(0, Nsteps):
#integrate orbits, calculate orbital radius (for Earth), and total energy
    time = t0 + it * dt
    t[it]=time
#    r_orb[it]=(u[1][0]**2+u[1][1]**2+u[1][2]**2)**0.5
    E[it]=TotalEnergy(u,mass,time)
    u = NbodyRK2(u,mass,time,dt)
    if it % max(1,Nsteps/100) == 0:
      print "it = %d, time = %g years, energy = %g" % \
            (it, time / seconds_per_year,
             TotalEnergy(u,mass,time))
#plot current postions at each iteration
      plt.clf()
      fig = plt.gcf()
      ax = mpl3d.Axes3D(fig)
      ax.scatter(u[:,0],u[:,1],u[:,2])
      ax.set_xlim((-rmax,rmax))
      ax.set_ylim((-rmax,rmax))
      ax.set_zlim((-rmax,rmax))
      plt.draw()
#save each position
    usave[:,:,it]=u
#to plot the orbital path
    for i in range(np.shape(usave)[0]):
      ax.plot(xs=usave[i,0],ys=usave[i,1],zs=usave[i,2])

# output result, and plot the entire orbit
plt.clf()
fig = plt.figure()
ax = fig.gca(projection='3d')
ax.set_xlim((-rmax,rmax))
ax.set_ylim((-rmax,rmax))
ax.set_zlim((-rmax,rmax))
for i in range(np.shape(usave)[0]):
    ax.plot(xs=usave[i,0],ys=usave[i,1],zs=usave[i,2])
plt.draw()
plt.savefig('final_orbits_GC.pdf')
file_header = "1:x 2:y 3:z 4:vx 5:vy 6:vz 7:mass"
np.savetxt(final_data_file, u, header=file_header)

#plot the orbital radius (for Earth)
#plt.clf()
#t=np.linspace(t0,t1,Nsteps)
#plt.plot(t/seconds_per_year,r_orb/cm_per_AU,'k')
#plt.xlabel('Time (years)',fontsize=20)
#plt.ylabel('Earth Orbital Radius (AU)',fontsize=20)
#plt.savefig('R_E.pdf')

#plot the relative energy
plt.clf()
plt.plot(t/seconds_per_year,E/E[0],'k')
plt.xlabel('Time (years)',fontsize=20)
plt.ylabel('Energy Relative to E$_0$ (ergs)',fontsize=20)
plt.savefig('E_GC.pdf')
