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


def analytic(x,t,x0,sigma):
    return np.exp(-((x-v*t)-x0)**2/(2*sigma**2))

def upwind_update(y):
    ynew=np.zeros(len(y))
    for j in range(1,len(ynew)):
	ynew[j]=y[j]-v*dt/dx*(y[j]-y[j-1])
    return ynew

def FTCS_update(y):
    ynew=np.zeros(len(y))
    for j in range(1,len(ynew)-1):
	ynew[j]=y[j]-v*dt/(2.0*dx)*(y[j+1]-y[j-1])
    return ynew

# set up the grid here. Use a decent number of zones;
# perhaps to get a dx of 0.1
x = np.arange(0,100,0.1)
# parameters
dx = x[1]-x[0]
v = 0.1

n = len(x)
y = np.zeros(n)
cfl = 0.5
dt = cfl*dx/v
t = 0.0

# for initial data
sigma = np.sqrt(15.0)
x0 = 30.0

#set up initial conditions
y = analytic(x,t,x0,sigma)

# evolve (and show evolution)
mpl.ion()
mpl.figure()
mpl.ylim(0,1)
mpl.plot(x,y,'x-') # numerical data
mpl.plot(x,analytic(x,t,x0,sigma),'r-') # analytic data
mpl.show()

yold2 = y
yold = y
ntmax = 500
err=np.zeros(ntmax)
for it in range(ntmax):
    t = it*dt
    # save previous and previous previous data
    yold2 = yold
    yold = y

    # get new data; ideally just call a function
    #y = ????
    #y=upwind_update(y)
    y=FTCS_update(y)

    # after update, apply boundary conditions
    # apply_bcs(x,y) 
    y=apply_bcs(x,y)

    # get analytic result for time t
    yana = analytic(x,t,x0,sigma)
    # compute error estimage
    # err = ???
    err[it]=np.abs(np.max(y)-np.max(yana))
    print "it = ",it,err[it]
    mpl.clf()
    mpl.ylim(0,1)
    # plot numerical result
    mpl.plot(x,y,'k-')
    # plot analytic results
    mpl.plot(x,yana,'r-')
    mpl.draw()
    if it==100:
	mpl.savefig('partc100.pdf')
    if it==300:
	mpl.savefig('partc300.pdf')
    if it==350:
	mpl.savefig('partc350.pdf')
    if it==499:
	mpl.savefig('partc499.pdf')


mpl.show()

time=np.zeros(ntmax)
for i in range(ntmax):
    time[i]=i*dt
mpl.clf()
mpl.plot(time,err,'k-')
mpl.xlabel('Time',fontsize=20)
mpl.ylabel('Error',fontsize=20)
mpl.yscale('log')
mpl.show()
mpl.savefig('onecFTCSerr.pdf')


