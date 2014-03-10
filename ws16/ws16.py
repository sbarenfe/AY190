#!/usr/bin/env python
import sys,math
import matplotlib.pyplot as mpl
import numpy as np

class mydata:
    def __init__(self,nzones):
        """ Primary 'container' of hydrodynamic data 
            rho, vel, eps, press are the PRIMITIVE vars
            
            The CONSERVED vars (see Euler eqns.) are
            rho, rho*v, rho*eps + 0.5*rho*vel**2 and they 
            are stored in self.q[0,:], q[1,:], and q[2,:]
            respectively.
           
            The variables ending with p and m are the
            reconstructed states at the outer cell edge of 
            cell i (plus) and at the inner cell edge of zone
            i (minus).
        """

        self.x     = np.zeros(nzones) # cell centers
        self.xi    = np.zeros(nzones) # cell LEFT interfaces
        self.rho   = np.zeros(nzones)
        self.rhop   = np.zeros(nzones)
        self.rhom   = np.zeros(nzones)
        self.vel    = np.zeros(nzones)
        self.velp   = np.zeros(nzones)
        self.velm   = np.zeros(nzones)
        self.eps    = np.zeros(nzones)
        self.epsp   = np.zeros(nzones)
        self.epsm   = np.zeros(nzones)
        self.press  = np.zeros(nzones)
        self.pressp = np.zeros(nzones)
        self.pressm = np.zeros(nzones)
        self.q     = np.zeros((3,nzones)) # conserved quantities
        self.qp    = np.zeros((3,nzones))  
        self.qm    = np.zeros((3,nzones))  
        self.n     = nzones
        self.g     = 3



    def setup_grid(self,xmin,xmax):
        dx = (xmax - xmin) / (self.n - self.g*2 - 1)
        xmin = xmin - self.g*dx
        xmax = xmax + self.g*dx
        for i in range(self.n):
            # cell centers
            self.x[i] = xmin + (i)*dx
        # cell LEFT interfaces
        for i in range(self.n):
            self.xi[i] = self.x[i] - 0.5*dx

    
    def setup_ID(self,rchange,gamma):
        # Shocktube initial data
        rho1 = 10.0
        rho2 = 1.0
        eps1 = 2.5
        eps2 = 1.795
        press1 = (gamma-1.0)*rho1*eps1
        press2 = (gamma-1.0)*rho2*eps2
        for i in range(self.n):
            if self.x[i] < rchange:
                self.rho[i] = rho1
                self.press[i] = press1
                self.eps[i] = eps1
                self.vel[i] = 0.0
            else:
                self.rho[i] = rho2
                self.press[i] = press2
                self.eps[i] = eps2
                self.vel[i] = 0.0


def prim2con(rho,vel,eps):
    """
    Convert to conserved variables
    """
    q = np.zeros((3,len(rho)))
    q[0,:] = rho
    q[1,:] = rho*vel
    q[2,:] = rho*eps+1./2.*rho*vel**2


    return q

def con2prim(q,gamma):
    """
    Convert back to primitive variables
    """
    rho = q[0,:]
    vel = q[1,:]/rho
    eps = (q[2,:]-1./2.*rho*vel**2)/rho
    press = eos_press(rho,eps,gamma)

    return (rho,eps,press,vel)

def apply_bcs(hyd):
    """
    Boundary conditions routine.
    """
    hyd.rho[0:hyd.g-1] = hyd.rho[hyd.g]
    hyd.vel[0:hyd.g-1] = hyd.vel[hyd.g]
    hyd.eps[0:hyd.g-1] = hyd.eps[hyd.g]
    hyd.press[0:hyd.g-1] = hyd.press[hyd.g]

    hyd.rho[hyd.n-hyd.g:hyd.n-1] = hyd.rho[hyd.n-hyd.g-1]
    hyd.vel[hyd.n-hyd.g:hyd.n-1] = hyd.vel[hyd.n-hyd.g-1]
    hyd.eps[hyd.n-hyd.g:hyd.n-1] = hyd.eps[hyd.n-hyd.g-1]
    hyd.press[hyd.n-hyd.g:hyd.n-1] = hyd.press[hyd.n-hyd.g-1]

    return hyd

def minmod(a,b):
    if(a*b < 0):
        mm = 0.0
    elif(abs(a)<abs(b)):
        mm=a
    else:
        mm=b
    return mm

def tvd_minmod_reconstruct(n,g,f,x,xi):
    fp = np.zeros(n)
    fm = np.zeros(n)
    for i in range(g-1,n-g+1):
        dx_up = x[i] - x[i-1]
        dx_down = x[i+1] - x[i]
        dx_m = x[i] - xi[i]
        dx_p = xi[i+1] - x[i]
        df_up = (f[i]-f[i-1]) / dx_up
        df_down = (f[i+1]-f[i]) / dx_down
        delta = minmod(df_up,df_down)
        fp[i] = f[i] + delta*dx_p
        fm[i] = f[i] - delta*dx_m

    return (fp,fm)

def signum(x,y):
    if(y >= 0):
        return abs(x)
    else:
        return -abs(x)

    
def tvd_mc_reconstruct(n,g,f,x,xi):
    fp = np.zeros(n)
    fm = np.zeros(n)
    for i in range(g-1,n-g+1):
        dx_up = x[i] - x[i-1]
        dx_down = x[i+1] - x[i]
        dx_m = x[i] - xi[i]
        dx_p = xi[i+1] - x[i]
        df_up = (f[i]-f[i-1]) / dx_up
        df_down = (f[i+1]-f[i]) / dx_down
        if(df_up*df_down < 0):
            delta = 0.0
        else:
            delta = signum(min(2.0*abs(df_up),2.0*abs(df_down),\
                               0.5*(abs(df_up)+abs(df_down))),\
                               df_up + df_down)

        fp[i] = f[i] + delta*dx_p
        fm[i] = f[i] - delta*dx_m

    return (fp,fm)


def reconstruct(hyd,gamma,type):

    if(type=='pc'):
        # piecewise constant reconstruction 
        for i in range(hyd.g-1,hyd.n-hyd.g+1):
            hyd.rhop[i] = hyd.rho[i]
            hyd.rhom[i] = hyd.rho[i]
            hyd.epsp[i] = hyd.eps[i] 
            hyd.epsm[i] = hyd.eps[i]
            hyd.velp[i] = hyd.vel[i]
            hyd.velm[i] = hyd.vel[i]


    elif(type=='minmod'):
        (hyd.rhop,hyd.rhom) = tvd_minmod_reconstruct(hyd.n,hyd.g,hyd.rho,hyd.x,hyd.xi)
        (hyd.epsp,hyd.epsm) = tvd_minmod_reconstruct(hyd.n,hyd.g,hyd.eps,hyd.x,hyd.xi)
        (hyd.velp,hyd.velm) = tvd_minmod_reconstruct(hyd.n,hyd.g,hyd.vel,hyd.x,hyd.xi)

    elif(type=='mc'):
        (hyd.rhop,hyd.rhom) = tvd_mc_reconstruct(hyd.n,hyd.g,hyd.rho,hyd.x,hyd.xi)
        (hyd.epsp,hyd.epsm) = tvd_mc_reconstruct(hyd.n,hyd.g,hyd.eps,hyd.x,hyd.xi)
        (hyd.velp,hyd.velm) = tvd_mc_reconstruct(hyd.n,hyd.g,hyd.vel,hyd.x,hyd.xi)
                
    else:
        print "reconstruction type not known; abort!"
        sys.exit()


    hyd.pressp = eos_press(hyd.rhop,hyd.epsp,gamma)
    hyd.pressm = eos_press(hyd.rhom,hyd.epsm,gamma)

    hyd.qp = prim2con(hyd.rhop,hyd.velp,hyd.epsp)
    hyd.qm = prim2con(hyd.rhom,hyd.velm,hyd.epsm)


    return hyd

def eos_press(rho,eps,gamma):
    press = (gamma - 1.0) * rho * eps
    return press

def eos_cs2(rho,eps,gamma):
    prs = (gamma - 1.0) * rho *eps
    dpde = (gamma - 1.0) * rho
    dpdrho = (gamma - 1.0) * eps
    cs2 = dpdrho + dpde * prs/(rho+1.0e-30)**2
    return cs2

def calc_dt(hyd,dtp):
    cs = np.sqrt(eos_cs2(hyd.rho,hyd.eps,gamma))
    dtnew = 1.0
    for i in range(hyd.g,hyd.n-hyd.g):
        dtnew = min(dtnew, (hyd.x[i+1]-hyd.x[i]) / \
                    max(abs(hyd.vel[i]+cs[i]), \
                        abs(hyd.vel[i]-cs[i])))

    dtnew = min(cfl*dtnew,1.05*dtp)
    return dtnew


def hlle(hyd):
    fluxdiff = np.zeros((3,hyd.n))
    # compute eigenvalues
    evl  = np.zeros((3,hyd.n))
    evr  = np.zeros((3,hyd.n))
    smin = np.zeros(hyd.n)
    smax = np.zeros(hyd.n)
    csp  = np.sqrt(eos_cs2(hyd.rhop,hyd.epsp,gamma))
    csm  = np.sqrt(eos_cs2(hyd.rhom,hyd.epsm,gamma))
    for i in range(1,hyd.n-2):
        evl[0,i] = hyd.velp[i]
        evr[0,i] = hyd.velm[i+1]
        evl[1,i] = hyd.velp[i] - csp[i]
        evr[1,i] = hyd.velm[i+1] - csm[i+1]
        evl[2,i] = hyd.velp[i] + csp[i]
        evr[2,i] = hyd.velm[i+1] + csm[i+1] 

        # min and max eigenvalues
        smin[i] = min(evl[0,i],evl[1,i],evl[2,i],\
                   evr[0,i],evr[1,i],evr[2,i],0.0)
        smax[i] = max(evl[0,i],evl[1,i],evl[2,i],\
                   evr[0,i],evr[1,i],evr[2,i],0.0)


    # set up flux left L and right R of the interface
    # at i+1/2
    fluxl = np.zeros((3,hyd.n))
    fluxr = np.zeros((3,hyd.n))

    # calculate numerical flux left and right of the
    # interface at i+1/2    
    # for example:
    # fluxl[0,:] corresponds the flux from the left cell
    # so it should be rhop[0,:] * velp[0,i], or, simply
    # hyd.qp[0,1] * hyd.velp[i]
    #
    # Similar expressions must be filled in for all
    # three Euler equations and the stuff coming from the
    # left and coming from the right.
    # Note that the states at the i+1/2 interface from the
    # right are qm[:,i+1]
    for i in range(1,hyd.n-2):
        fluxl[0,i] = hyd.qp[0,i]*hyd.velp[i]
        fluxl[1,i] = hyd.qp[1,i]*hyd.velp[i]+hyd.pressp[i]
        fluxl[2,i] = (hyd.qp[2,i]+hyd.pressp[i])*hyd.velp[i]

        fluxr[0,i] = hyd.qm[0,i+1]*hyd.velm[i+1]
        fluxr[1,i] = hyd.qm[1,i+1]*hyd.velm[i+1]+hyd.pressm[i+1]
        fluxr[2,i] = (hyd.qm[2,i+1]+hyd.pressm[i+1])*hyd.velm[i+1]

    # solve the Riemann problem for the i+1/2 interface
    # with the HLLE solver
    ds = smax - smin
    flux = np.zeros((3,hyd.n))
    for i in range(hyd.g-1,hyd.n-hyd.g+1):
        flux[:,i] = ( smax[i]*fluxl[:,i] \
                      - smin[i]*fluxr[:,i] \
                      + smax[i]*smin[i] \
                      * (hyd.qm[:,i+1] - hyd.qp[:,i])) \
                      / ds[i]

    # flux differences
    for i in range(hyd.g,hyd.n-hyd.g):
        rm = hyd.xi[i]
        rp = hyd.xi[i+1]
        dxi = 1.0/(rp - rm)
        fluxdiff[:,i] = dxi * (flux[:,i] - flux[:,i-1])

    return fluxdiff


def calc_rhs(hyd):
    hyd = reconstruct(hyd,gamma,reconstruction_type)
    fluxdiff = hlle(hyd)
    return -fluxdiff

#########################################################
def run():
    # Global parameters
    gamma = 1.4
    cfl = 0.5
    dt = 1.0e-5
    dtp = dt
    reconstruction_type = 'mc' # minmod, mc, pc
    nzones = 500
    tend = 0.4

    # initialize
    hyd = mydata(nzones)
    
    # set up grid
    xmin = -0.5
    xmax = 0.5
    rchange = 0.0
    hyd.setup_grid(xmin,xmax)

    # set up initial data
    hyd.setup_ID(rchange,gamma)
    
    # get initial timestep
    dt = calc_dt(hyd,dt)
    
    # initial prim2con
    hyd.q = prim2con(hyd.rho,hyd.vel,hyd.eps)

    t = 0.0
    i = 0
    
    mpl.ion()
    mpl.figure()
    mpl.plot(hyd.x,hyd.rho,"r-")
    mpl.xlabel('x',fontsize=20)
    mpl.ylabel('Density',fontsize=20)
    timestring = "t=%5.3f" % (t)
    ax = mpl.gca()
    mpl.text(0.7,0.88,timestring,fontsize=25,
    	    horizontalalignment="left",rotation="horizontal",
    	    transform=ax.transAxes)
    mpl.savefig('t0.pdf')
    mpl.show()
    
    while(t < tend):

    	# some convenience output of the density
    	if(i % 10 == 0):
    	    print "%5d %15.6E %15.6E" % (i,t,dt)
    	    mpl.clf()
    	    mpl.plot(hyd.x,hyd.rho,"r-")
    	    mpl.xlabel('x',fontsize=20)
    	    mpl.ylabel('Density',fontsize=20)
    
    	    timestring = "t=%5.3f" % (t)
    	    ax = mpl.gca()
    	    mpl.text(0.7,0.88,timestring,fontsize=25,
    		    horizontalalignment="left",rotation="horizontal",
    		    transform=ax.transAxes)
    	    mpl.draw()
    
    	if i==70:
    	    mpl.savefig('t70.pdf')
    
    	if i==100:
    	    mpl.savefig('t100.pdf')
    
    	if i==200:
    	    mpl.savefig('t200.pdf')
    
    	if i==350:
    	    mpl.savefig('t350.pdf')
    
    	if i==450:
    	    mpl.savefig('t450.pdf')
    
    	if i==650:
    	    mpl.savefig('t650.pdf')
    
    	if i==700:
	    mpl.savefig('t700.pdf')

    	if i==505:
    	    rho_save=hyd.rho
    	    x_save=hyd.x
	    v_save=hyd.vel
	    p_save=hyd.press
	    eps_save=hyd.eps
    	# calculate new timestep
    	dt = calc_dt(hyd,dt)
    
    	hydold = hyd
    	qold = hyd.q

    	# calc rhs
    	k1 = calc_rhs(hyd)
    	# calculate intermediate step
    	hyd.q = qold + 1.0/2.0 * dt * k1
    	# con2prim
    	(hyd.rho,hyd.eps,hyd.press,hyd.vel) = con2prim(hyd.q,gamma)
    	# boundaries
    	hyd = apply_bcs(hyd)

    	#calc rhs
    	k2 = calc_rhs(hyd)
    	#apply update
    	hyd.q = qold + dt * (0.5 * k1 + 0.5 * k2)
    	# con2prim
    	(hyd.rho,hyd.eps,hyd.press,hyd.vel) = con2prim(hyd.q,gamma)
    	# apply bcs
    	hyd = apply_bcs(hyd)

    	# update time
    	t = t + dt
    	i = i + 1
    
    
    mpl.ioff()
    
    mpl.clf()
    mpl.xlabel('x',fontsize=20)
    mpl.ylabel('Density',fontsize=20)
    mpl.plot(hyd.x,hyd.rho,"r-")
    timestring = "t=%5.3f" % (t)
    ax = mpl.gca()
    mpl.text(0.7,0.88,timestring,fontsize=25,
    	    horizontalalignment="left",rotation="horizontal",
    	    transform=ax.transAxes)
    mpl.savefig('tend.pdf')
    mpl.show()
    
    return x_save,rho_save,v_save,p_save,eps_save
    
    
def plot_things(x15,x16,rho15,rho16,v15,v16,p15,p16,e15,e16):
	mpl.clf()
	mpl.xlabel('x',fontsize=20)
	mpl.ylabel('Density',fontsize=20)
	mpl.xlim(-0.5,0.5)
	a,=mpl.plot(x15,rho15,'k--')
	b,=mpl.plot(x15,rho16,'k-')
	mpl.legend([a,b],['Worksheet 15','Worksheet 16'],loc=3)
	mpl.savefig('rhocomp.pdf')
	
	mpl.clf()
	mpl.xlabel('x',fontsize=20)
	mpl.ylabel('Velocity',fontsize=20)
	mpl.xlim(-0.5,0.5)
	a,=mpl.plot(x15,v15,'k--')
	b,=mpl.plot(x15,v16,'k-')
	mpl.legend([a,b],['Worksheet 15','Worksheet 16'],loc=2)
	mpl.savefig('vcomp.pdf')

	mpl.clf()
	mpl.xlabel('x',fontsize=20)
	mpl.ylabel('Pressure',fontsize=20)
	mpl.xlim(-0.5,0.5)
	a,=mpl.plot(x15,p15,'k--')
	b,=mpl.plot(x15,p16,'k-')
	mpl.legend([a,b],['Worksheet 15','Worksheet 16'],loc=3)
	mpl.savefig('pcomp.pdf')

	mpl.clf()
	mpl.xlabel('x',fontsize=20)
	mpl.ylabel('Energy',fontsize=20)
	mpl.xlim(-0.5,0.5)
	a,=mpl.plot(x15,e15,'k--')
	b,=mpl.plot(x15,e16,'k-')
	mpl.legend([a,b],['Worksheet 15','Worksheet 16'],loc=3)
	mpl.savefig('ecomp.pdf')
    
    
    
