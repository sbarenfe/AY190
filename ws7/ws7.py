import numpy as np
import matplotlib.pyplot as plt

def one():
#Generate randome (x,y) points, keeping track of how many are inside the circle
	np.random.seed(1)
	N=np.linspace(10,10000,1000)#do this for 10 to 10000 points
	blueberry=np.zeros(len(N))
	for j in range(len(N)):
		n=0.0
		x=np.random.rand(N[j])-0.5
		y=np.random.rand(N[j])-0.5
		r=((x)**2+(y)**2)**0.5
		for i in range(len(r)):
			if r[i]<=0.5:
				n+=1
		blueberry[j]=(4.*n)/N[j]#the aptly name fraction inside the circle

	lemon_meringue=np.pi*np.ones(len(N))#measured value of pi
	apple=np.abs(blueberry-lemon_meringue)#and its error from the real value
#plot the error vs. N to show convergence
	plt.clf()
	a,=plt.plot(N,apple,'k-')
	plt.xlabel('Number of Points Drawn',fontsize=20)
	plt.ylabel('Error',fontsize=20)
	plt.xscale('log')
	plt.yscale('log')
	print '%.4f' %blueberry[-2]
#and check that the convergence goes like 1/N^0.5	
	x=np.linspace(10,10000,100)
	y=1.2/x**(0.5)
	b,=plt.plot(x,y,'k--')
	plt.legend([a,b],['Error in Calculation','$y\propto x^{-1/2}$'],loc=3)
	plt.savefig('ws7a.pdf')

def two():
	N=10000
	num=np.arange(1,51,1)
	frac=np.zeros(len(num))
	for l in range(len(num)):#for groups of 1 to 50 people
		print num[l]
		n=0.0
		for i in range(N):#generate birthdays for the group 10000 times
			people=np.random.randint(1,366,num[l])
			found=False
			for j in range(len(people)):#check if any birthdays are the same
				for k in range(j+1,len(people)):
					if people[j]==people[k]:#and increment the counter if they are
						n+=1
						found=True#but only do this once (don't count two pairs as two events)
						break
				if found:
					break
		frac[l]=n/N#find the total probability for each sized group
	print frac
#and plot the results
	plt.clf()
	plt.scatter(num,frac,c='k')
	plt.xlim(0,50)
	plt.ylim(0,1)
	x=np.linspace(0,50,10)
	y=0.5*np.ones(len(x))
	plt.plot(x,y,'k')
	plt.xlabel('Number of People',fontsize=20)
	plt.ylabel('Probability',fontsize=20)
	plt.savefig('ws7b.pdf')

