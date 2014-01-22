import numpy as np

def find_root(guess,omega,t,e):
	Enew=guess
	ratio=1.
	step=0
#expression I am finding the root of:
	f=lambda x:x-omega*t-e*np.sin(x)
#and it's derivative
	fprime=lambda x:1-e*np.cos(x)
#implement Newton's method until convergence
	while np.abs(ratio) >= 10**(-10):
		Eold=Enew
		Enew=Eold-f(Eold)/fprime(Eold)
		step+=1
		ratio=np.abs(Enew-Eold)/Eold
#return the eccentric anomaly, number of steps required, and ratio between the last two steps
	return Enew,step,ratio
	
def one():
#define constants:
	T=365.25635 #days
	#e=0.0167 #part a)
	e=0.99999 #part b)
	a=1.0 #AU
	b=a*(1-e**2)**(0.5)
	omega=2*np.pi/T

	t=np.array([91,182,273])
	E=np.zeros(len(t))
	steps=np.zeros(len(t))
	ratio=np.zeros(len(t))
#use the find_root fucntion above to find E for each time
	for i in range(len(E)):
		E[i],steps[i],ratio[i]=find_root(2,omega,t[i],e)	
#calculate x and y
	x=a*np.cos(E)
	y=b*np.sin(E)
#print the results
	print E
	print steps
	print ratio	
	print x
	print y

