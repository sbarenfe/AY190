import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate

#problem 1
def one():
	#define x0 and x1
	x0=np.float32(1.)
	x1=np.float32(1./3.)
	#define the recursion relation
	xnp1=lambda xn,xnm1:np.float32(13./3.)*np.float32(xn)-np.float32(4./3.)*np.float32(xnm1)
	#define the array to store the x values in
	x=np.zeros(16,dtype=np.float32)
	x[0]=x0
	x[1]=x1
	for i in range(2,len(x)):
		x[i]=xnp1(x[i-1],x[i-2])
	print np.float32(x)

	#Calculate the exact relation
	y=np.zeros(len(x))
	for i in range(len(y)):
		y[i]=(1./3.)**i
	print y

	#Calculate the error
	abserr=np.abs(x-y)
	relerr=np.abs(x-y)/y
	
	#and print it for n=15
	print abserr[15]
	print relerr[15]

#problem 2
#Calculate forward and central derivatives
def derivative(h):
	x=np.linspace(-2,6,8./h+1)
	#function we are differentiating:
	f=x**3-5*x**2+x
	fprimefwd=np.zeros(len(f))
	#forward differencing scheme:
	for i in range(len(fprimefwd)-1):
		fprimefwd[i]=(f[i+1]-f[i])/h

	#analytical derivative
	fprimetrue=3*x**2-10*x+1

	fprimecen=np.zeros(len(f))
	#central differencing scheme:
	for i in range(1,len(fprimecen)-1):
		fprimecen[i]=(f[i+1]-f[i-1])/(2.*h)
	
	#calculate the numerical errors
	fwderr=np.abs(fprimetrue-fprimefwd)
	cenerr=np.abs(fprimetrue-fprimecen)

	return fwderr,cenerr
		
	
def two():
	h=0.1
	x1=np.linspace(-2,6,8./h+1)
	x2=np.linspace(-2,6,8./(h/2.)+1)

	#use my derivative function above to find the forward and 
	#central errors for the two different h values
	f1,c1=derivative(h)
	f2,c2=derivative(h/2.)

	#plot the results
	plt.clf()
	a,=plt.plot(x1,f1,'k-')
	b,=plt.plot(x2,f2,'k--')
	c,=plt.plot(x1,f1/2.,'k+')
	plt.xlim(-2,6)
	plt.ylim(0,1)
	plt.xlabel('X',fontsize=20)
	plt.ylabel('Absolute Error',fontsize=20)
	plt.legend([a,b,c],["h=0.1","h=0.05","Half of h=0.1"],loc=9)
	plt.savefig('twoi.pdf')	
	plt.show()

	plt.clf()
	d,=plt.plot(x1,c1,'k-')
	e,=plt.plot(x2,c2,'k--')
	f,=plt.plot(x1,c1/4.,'k+')
	plt.xlim(-2,6)
	plt.ylim(0,.02)
	plt.xlabel('X',fontsize=20)
	plt.ylabel('Absolute Error',fontsize=20)
	plt.legend([d,e,f],["h=0.1","h=0.05","1/4 of h=0.1"],loc=9)	
	plt.savefig('twoii.pdf')	
	plt.show()

#Problem 4, Lagrange interpolation
def four():
	#data
	t=np.array([0.0,0.2,0.3,0.4,0.5,0.6,0.7,0.8,1.0])
	m=np.array([0.302,0.185,0.106,0.093,0.24,0.579,0.561,0.468,0.302])

	x=np.linspace(0,1,100)
	#Calculate L_nj given a j and time value
	def L(j,time):
		prod=1.
		for i in range(len(t)):
			if i==j:
				continue
			prod=prod*(time-t[i])/(t[j]-t[i])
		return prod

	p=np.zeros(len(x))
	#Calculate p(x) using the formula in the notes
	for i in range(len(p)):
		for j in range(len(t)):
			p[i]+=m[j]*L(j,x[i])
	#plot the results
	plt.clf()
	plt.scatter(t,m,c='k')
	plt.plot(x,p,'k-')
	plt.xlabel('Time (days)',fontsize=20)
	plt.ylabel('Apparent Magnitude',fontsize=20)
	plt.ylim(0,3)
	plt.savefig('foura.pdf')
	plt.show()

	#Compare to spline interpolation
	tck=interpolate.splrep(t,m)
	yspl=interpolate.splev(x,tck,der=0)

	plt.clf()
	plt.scatter(t,m,c='k')
	a,=plt.plot(x,p,'k-')
	b,=plt.plot(x,yspl,'k--')
	plt.xlabel('Time (days)',fontsize=20)
	plt.ylabel('Apparent Magnitude',fontsize=20)
	plt.ylim(0,3)
	plt.legend([a,b],['Lagrange','Spline'],loc=((0.4,0.8)))
	plt.savefig('five2.pdf')
	plt.show()
	
#Problem 5, spline interpolation
def five():
	#data
	t=np.array([0.0,0.2,0.3,0.4,0.5,0.6,0.7,0.8,1.0])
	m=np.array([0.302,0.185,0.106,0.093,0.24,0.579,0.561,0.468,0.302])
	
	#Python's built in spline interpolator
	tck=interpolate.splrep(t,m)
	x=np.linspace(0,1,100)
	y=interpolate.splev(x,tck,der=0)

	#plot the results
	plt.clf()
	plt.scatter(t,m,c='k')
	plt.plot(x,y,'k-')
	plt.xlabel('Time (days)',fontsize=20)
	plt.ylabel('Apparent Magnitude',fontsize=20)
	plt.ylim(0,1)
	plt.savefig('fiveb.pdf')
	plt.show()
	
