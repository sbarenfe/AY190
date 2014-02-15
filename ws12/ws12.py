import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate

def one():
	r,col1,col2=np.loadtxt(fname='presupernova.dat',usecols=(2,3,4),unpack=True)
	plt.clf()
	a,=plt.plot(r,col1)
	b,=plt.plot(r,col2)
	#plt.xscale('log')
	#plt.yscale('log')
	plt.xlim(0,10**8)
	plt.xlabel('Radius (cm)',fontsize=20)
	plt.ylabel('Column',fontsize=20)	
	plt.legend([a,b],['Column 3','Column 4'])
	plt.savefig('onea.pdf')

	plt.clf()
	plt.plot(r,col2,'k-')
	plt.xscale('log')
	plt.yscale('log')
	plt.xlim(0,10**9)
	plt.ylim(10**4,10**11)
	plt.xlabel('Radius (cm)',fontsize=20)
	plt.ylabel('Density (g/cm$^3$)',fontsize=20)
	plt.savefig('onea2.pdf')

def mk_grid(n):
	rad,den=np.loadtxt(fname='presupernova.dat',usecols=(2,4),unpack=True)
	#Python's built in spline interpolator
	tck=interpolate.splrep(rad,den)
	r=np.linspace(0,10**9,n)
	rho=interpolate.splev(r,tck,der=0)
	return r,rho,rad,den

def two():
	r,rho,rad,den=mk_grid(100)
	plt.clf()
	a,=plt.plot(rad,den,'k-')
	b,=plt.plot(r,rho,'k--')
	plt.xscale('log')
	plt.yscale('log')
	plt.xlim(0,10**9)
	plt.ylim(10**4,10**11)
	plt.xlabel('Radius (cm)',fontsize=20)
	plt.ylabel('Density (g/cm$^3$)',fontsize=20)
	plt.legend([a,b],['Data Table Values','Interpolation'])	
	plt.savefig('two.pdf')

def three():
#Set up grid
	npoints=1000
	#r,rho=mk_grid(npoints)
	G=6.67*10**(-8)
	r=np.linspace(0,10**9,npoints)
	rho=np.ones(len(r))
	dr=r[1]-r[0]
	r+=0.5*dr

#Boundary conditions
	A=0 #phi(0)
	B=0 #phi'(0)

#Integrate
	yy=integrate_FE(r,A,B)

#Shift Phi to match outer BC
	dm=4*np.pi*r**2*dr*rho
	M=np.sum(dm)
	Phi_out=-G*M/r[npoints-1]
	offset=Phi_out-yy[npoints-1,0]
	Phi=yy[:,0]+offset

#Analytical Solution
	phi_ana=2./3.*np.pi*G*rho*(r**2-3*(r[npoints-1]**2))

#plot the result
	plt.clf()
	a,=plt.plot(r,Phi,'k-')
	b,=plt.plot(r,phi_ana,'k--')
	plt.xlim(0,10**9)
	plt.xlabel('Radius (cm)',fontsize=20)
	plt.ylabel('$\Phi$ (cgs)',fontsize=20)
	plt.legend([a,b],['Numerical','Analytical'])	
	plt.savefig('hsphere.pdf')

#Plot error
#	plt.clf()
#	err=np.abs(Phi-phi_ana)/(-Phi)
#	plt.plot(r,err,'k-')
#	plt.xlabel('Radius (cm)',fontsize=20)
#	plt.ylabel('Error',fontsize=20)
#	plt.xlim(0,10**9)
#	plt.savefig('hsphere_err100.pdf')


def integrate_FE(x,A,B):
    # forward-Euler integrator 
        
    # make an array for all points
    # entry 0 contains y
    # entry 1 contains y'
    npoints=len(x)
    #trash,rho=mk_grid(npoints)
    rho=np.ones(len(x))
    dx=x[1]-x[0]
    yy = np.zeros((npoints,2))    

    yy[0,0] = A # boundary value A for y at x=0
    yy[0,1] = B # guessed boundary value for y' at x=0

    for i in range(npoints-1):
        yy[i+1,:] = yy[i,:] + dx*calc_rhs(yy[i,1],x[i],rho[i])
	
    return yy

def calc_rhs(z,x,rho):
    # rhs routine
    # rhs[0] is rhs for y
    # rhs[1] is rhs for u
    G=6.67*10**(-8)
    rhs = np.zeros(2)
    rhs[0] = z 
    rhs[1] = 4*np.pi*G*rho-2.0/x*z

    return rhs

