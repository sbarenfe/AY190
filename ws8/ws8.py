import numpy as np
import matplotlib.pyplot as plt

#!/usr/bin/env python

import numpy as np
import scipy as sp


# global constants
ggrav = 6.67e-8
msun  = 1.99e33

# EOS parameters
# for white dwarfs:
polyG = 4.0/3.0
polyK = 1.244e15*0.5**polyG


#######################################
# function definitions
def tov_RHS(rad,p,rho,m):
    
    # RHS function
    
    rhs = np.zeros(2)
    if(rad > 1.0e-10):
        rhs[0] = -ggrav*m*rho/rad**2
        rhs[1] = 4*np.pi*rho*rad**2
    else:
        rhs[0] = 0.0
        rhs[1] = 0.0

    return rhs

def tov_integrate_FE(rad,dr,p,rho,m,method):

    # Forward-Euler Integrator

    new = np.zeros(2)
    old = np.zeros(2)
    old[0] = p
    old[1] = m

    # forward Euler integrator
    if method=='Euler':
        new = old + dr*tov_RHS(rad,p,rho,m)
    if method=='RK2':
	k1=dr*tov_RHS(rad,p,rho,m)
        rho = ((p+0.5*k1[0])/polyK)**(1./polyG)
        new = old + dr*tov_RHS(rad+0.5*dr,p+0.5*k1[0],rho,m+0.5*k1[1])
    # assign outputs
    pnew = new[0]
    mnew = new[1]
    
    return (pnew,mnew)

#######################################

# set up grid
npoints = 2000
radmax = 2.0e8 # 2000 km
radius = np.linspace(0,radmax,npoints)
dr = radius[1]-radius[0]

# set up variables
press = np.zeros(npoints)
rho   = np.zeros(npoints)
mass  = np.zeros(npoints)

# set up central values
rho[0]   = 1.0e10
press[0] = polyK * rho[0]**polyG
mass[0]  = 0.0

# set up termination criterion
press_min = 1.0e-10 * press[0]

nsurf = 0
for n in range(npoints-1):
    
    (press[n+1],mass[n+1]) = tov_integrate_FE(radius[n],
                                              dr,
                                              press[n],
                                              rho[n],mass[n],'RK2')
    # check for termination criterion
    if(press[n+1] < press_min and nsurf==0):
        nsurf = n

    if(n+1 > nsurf and nsurf > 0):
        press[n+1] = press[nsurf]
        rho[n+1]   = rho[nsurf]
        mass[n+1]  = mass[nsurf]

    # invert the EOS to get density
    rho[n+1] = (press[n+1]/polyK)**(1./polyG)


print radius[nsurf]/1.0e5
print mass[nsurf]/msun
print rho[0]
print press[0]

plt.clf()
fig=plt.figure()
ax=fig.add_subplot(111)
ax.set_xlabel(r'$r$ (km)',fontsize=20)
ax.set_ylabel(r'$\frac{\rho(r)}{\rho_{core}}$, $\frac{P(r)}{P_{core}}$',fontsize=20)
ax.set_xlim(0,radius[nsurf]/1.0e5)
ax.set_ylim(0,1)
a,=ax.plot(radius/1.0e5,rho/rho[0],'k-')
b,=ax.plot(radius/1.0e5,press/press[0],'k--')

ax2=ax.twinx()
ax2.set_ylabel(r'$\frac{M(r)}{M_{\odot}}$',fontsize=20)
ax2.set_xlim(0,radius[nsurf]/1.0e5)
ax2.set_ylim(0,mass[nsurf]/msun)
c,=ax2.plot(radius/1.0e5,mass/msun,'k-.')
plt.legend([a,b,c],[r'$\frac{\rho(r)}{\rho_{core}}$',r'$\frac{P(r)}{P_{core}}$',r'$\frac{M(r)}{M_{\odot}}$'],loc=7)
plt.subplots_adjust(right=0.85)
plt.savefig('star.pdf')


