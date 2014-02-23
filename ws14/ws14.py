import sys,math
import numpy as np
import matplotlib.pyplot as mpl
import scipy as sp

def apply_bcs(x,y):
    # apply boundary conditions
    # you need to fill in code
    y[0]=y[1]
    y[len(y)-1]=y[len(y)-2]
    return y

#inintial condition
def analytic(x,t,y):
    return 1./8.*np.sin(2*np.pi*(x-y*t)/np.float64(L))

#use upwind and downwind schemes for positive and negative u regions, respectively
def upwind_update(y):
    ynew=np.zeros(len(y))
    for j in range(1,len(ynew)-1):
	if y[j]>=0:
    	    ynew[j]=y[j]-y[j]*dt/dx*(y[j]-y[j-1])
	if y[j]<0:
	    ynew[j]=y[j]-y[j]*dt/dx*(y[j+1]-y[j])
    return ynew

def FTCS_update(y):
    ynew=np.zeros(len(y))
    for j in range(1,len(ynew)-1):
	ynew[j]=y[j]-y[j]*dt/(2.0*dx)*(y[j+1]-y[j-1])
    return ynew

# set up the grid here. Use a decent number of zones;
# perhaps to get a dx of 0.1
L=100
x = np.arange(0,L,.1)
# parameters
dx = x[1]-x[0]

n = len(x)
y = np.zeros(n)

t = 0.0
#set up initial conditions
y = analytic(x,t,y)

dt = 5*dx

# evolve (and show evolution)
mpl.ion()
mpl.figure()
#mpl.ylim(0,1)
mpl.plot(x,y,'x-') # numerical data
mpl.plot(x,analytic(x,t,y),'r-') # analytic data
mpl.show()

yold2 = y
yold = y
y0=y
ntmax = 500
#err=np.zeros(ntmax)
#yana=analytic(x,0,0)
for it in range(ntmax):
    t = it*dt  
    if t%1==0:    
	print t
    # save previous and previous previous data
    yold2 = yold
    yold = y

    # get new data; ideally just call a function
    #y = ????
    y=upwind_update(y)
    #y=FTCS_update(y)

    # after update, apply boundary conditions
    # apply_bcs(x,y) 
    y=apply_bcs(x,y)

    mpl.clf()
    mpl.ylim(-.2,.2)
    # plot numerical result
    a,=mpl.plot(x,y,'k-')
    # plot analytic results
    b,=mpl.plot(x,y0,'k--')
    mpl.xlabel('x',fontsize=20)
    mpl.ylabel('u',fontsize=20)
    mpl.legend([a,b],['u(x,t)','u$_0$(x,t=0)'],loc=1)
    mpl.draw()
    if t==0:
	mpl.savefig('t0.pdf')
    if t==50:
	mpl.savefig('t50.pdf')
    if t==100:
	mpl.savefig('t100.pdf')
    if t==140:
	mpl.savefig('t140.pdf')
    if t==249:
	mpl.savefig('t249.pdf')


mpl.show()


