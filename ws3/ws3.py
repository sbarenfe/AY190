import numpy as np
import matplotlib.pyplot as plt
from scipy import special as sp

def one():
	x=np.linspace(0,np.pi,101)#grid points
	h=x[1]-x[0]#step size
	Q=np.zeros(len(x)-1)
	#f=lambda x:np.sin(x) # part a)
	f=lambda x:x*np.sin(x) # part b)

	print "For %i steps:"%len(Q)
#midpoint rule:
	for i in range(len(Q)):
		Q[i]=h*f((x[i]+x[i+1])/2.)
	Qtot=np.sum(Q)
	print "Midpoint rule gives: %.10f"%Qtot,np.abs(Qtot-np.pi)

#trapezoidal rule:
	for i in range(len(Q)):
		Q[i]=h*(f(x[i])+f(x[i+1]))/2.
	Qtot=np.sum(Q)
	print "Trapezoidel rule gives: %.10f"%Qtot,np.abs(Qtot-np.pi)

#Simpson's Rule:
	for i in range(len(Q)):
		Q[i]=h/6.*(f(x[i])+4*f((x[i]+x[i+1])/2.)+f(x[i+1]))
	Qtot=np.sum(Q)
	print "Simpson's rule gives: %.11f"%Qtot,np.abs(Qtot-np.pi)

def twoa(n):
	kT=3.2*10**(-5)#cgs
	hc=3.16*10**(-17)#cgs
#use Python's built-in Gaussian-Laguerre roots and weights
	[la_r,la_w]=sp.l_roots(n)
	f=lambda x:x**2*np.exp(x)/(1.+np.exp(x))#f(x)
	Q=np.zeros(len(la_r))	
#calculate each Qi
	for i in range(len(Q)):
		Q[i]=la_w[i]*f(la_r[i])
#and the sum
	Qtot=np.sum(Q)
#then put the constants back in
	ne=(Qtot*8.*np.pi*(kT)**3)/(2.*np.pi*hc)**3
	print "The integral is %.11f"%Qtot
	print "ne is %g cm^-3"%ne
	
def twob():
	n=20 #matches the n=20 used for highest accuracy in part a)
	deltaE=5 #MeV
	kT=3.2*10**(-5) #cgs
	hc=3.16*10**(-17) #cgs
#use Python's built-in Gauss-Legendre roots and weights
	[le_r,le_w]=sp.p_roots(n)
#define the substitution I made
	x=lambda u,i:(5.*u+10.*i+5.)/(2.)
#and f(x)
	f=lambda u,i:(x(u,i))**2/(np.exp(x(u,i))+1)
	ne=np.zeros(150./5.)#array of n_e for each energy bin
	Q=np.zeros(len(le_r))

#loop over each energy bin
	for i in range(len(ne)):
#and calculate n_e in each bin
		for j in range(len(Q)):
			Q[j]=le_w[j]*f(le_r[j],i)
		Qtot=np.sum(Q)
		ne[i]=5./2.*(Qtot*8.*np.pi*(kT)**3)/(2.*np.pi*hc)**3
	dndE=ne/deltaE
#print the result
	print dndE
#and check if the total n_e matches part a)
	print np.sum(dndE)*deltaE


