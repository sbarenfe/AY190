import numpy as np

def find_root(guess,omega,t,e):
	Enew=guess
	ratio=1.
	step=0
	f=lambda x:x-omega*t-e*np.sin(x)
	fprime=lambda x:1-e*np.cos(x)
	while ratio >= 10**(-10):
		Eold=Enew
		Enew=Eold-f(Eold)/fprime(Eold)
		step+=1
		ratio=np.abs(Enew-Eold)/Eold
	return Enew,step,ratio
	
def one():
	T=365.25635
	e=0.0167
	a=1.0
	b=a*(1-e**2)**(0.5)
	omega=2*np.pi/T
	t=np.array([91,182,273])
	E=np.zeros(len(t))
	steps=np.zeros(len(t))
	ratio=np.zeros(len(t))
	for i in range(len(E)):
		E[i],steps[i],ratio[i]=find_root(100,omega,t[i],e)	
	print E
	print steps
	print ratio	
