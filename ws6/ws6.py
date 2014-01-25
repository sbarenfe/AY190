import numpy as np
from timeit import timeit
import matplotlib.pyplot as plt

def dft(x,time=False):
	N=len(x)
	W=np.zeros((N,N),dtype=complex)
#Calculate the W matrix
	for k in range(N):
		for j in range(N):
			W[k][j]=np.exp((-2.*np.pi*1J*j*k)/N)
#and dot it with x
	y=np.dot(W,x)
#print out the results, but only if I'm not trying to time the function.
	if time==False:
		print "My DFT function gives:"	
		print y
		y2=np.fft.fft(x)
		print "Python's FFT funciton gives:"
		print y2
		print "The difference is:"
		print y-y2

def time():
	N=np.arange(10,101,1)
	t=np.zeros(len(N))
#Time the speed of my DFT function for N from 10 to 100
	for i in range(len(N)):
		t[i]=timeit("dft(x,True)",number=10,setup="from ws6 import dft; import pylab; x=pylab.randn(%d)" %N[i])
	tfft=np.zeros(len(N))
#and Python's FFT function over the same range of N
	for i in range(len(N)):
		tfft[i]=timeit("fft(x)",number=10,setup="from numpy.fft import fft; import pylab; x=pylab.randn(%d)" %N[i])
	x=np.linspace(10,100,100)
	y=x**2/16000.
#Plot the Results
	plt.clf()
	a,=plt.plot(N,t)
	b,=plt.plot(x,y)
	c,=plt.plot(N,tfft)
	plt.xscale('log')
	plt.yscale('log')
	plt.xlabel('N',fontsize=20)
	plt.ylabel('Time to Run (s)',fontsize=20)
	plt.legend([a,b,c],["My DFT",r"$t=\frac{x^2}{16000}$","Python's FFT"],loc=2)
	plt.savefig('time.pdf')

#So, who's trying for the second fastest DFT?
def fast_dft(x):
	print "The Fourier Transform of x is y."
