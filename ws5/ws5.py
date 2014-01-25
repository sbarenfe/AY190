import numpy as np
import matplotlib.pyplot as plt

def one():
#part a)
#Read in the data
	sigma,sigerr,logM,logMerrH,logMerrL=np.loadtxt(fname='m_sigma_table.dat',usecols=(2,3,9,10,11),unpack=True)
	logsigma=np.log10(sigma)
#and plot it
	plt.subplots_adjust(bottom=0.15)
	plt.scatter(logsigma,logM,c='k')
	plt.xlabel(r"$\log\left(\frac{\sigma_*}{km/s}\right)$",fontsize=20)
	plt.ylabel(r'$\log\left(\frac{M_{BH}}{M_{\odot}}\right)$',fontsize=20)	
	plt.savefig('onea.pdf')
#part  b)
#Set up the variables used to find a1 and a2
	S=len(sigma)
	sigmax=np.sum(logsigma)
	sigmay=np.sum(logM)
	sigmax2=np.sum(logsigma*logsigma)
	sigmaxy=np.sum(logsigma*logM)
#Solve for a1 and a2
	a1=(sigmay*sigmax2-sigmax*sigmaxy)/(S*sigmax2-(sigmax)**2)
	a2=(S*sigmaxy-sigmay*sigmax)/(S*sigmax2-(sigmax)**2)
	print a1,a2
#Plot my fit and the Greene 2006 fit	
	x=np.linspace(1.4,2.6,100)
	y=a1+a2*x
	y2=7.96-4.02*np.log10(200)+4.02*x

	plt.clf()
	plt.scatter(logsigma,logM,c='k')
	a,=plt.plot(x,y,'k-')
	b,=plt.plot(x,y2,'k--')
	plt.xlim(1.4,2.6)
	plt.ylim(4,9)
	plt.xlabel(r"$\log\left(\frac{\sigma_*}{km/s}\right)$",fontsize=20)
	plt.ylabel(r'$\log\left(\frac{M_{BH}}{M_{\odot}}\right)$',fontsize=20)
	plt.legend([a,b],['My Fit','Greene & Ho Fit'],loc=2)
	plt.savefig('oneb.pdf')
#part c)
#Propogate the sigma uncertainties
	logsigerr=sigerr/(sigma*np.log(10))
#Set the logM uncertainty to the max of the given values (if two are given)
	logMerr=np.zeros(len(logMerrL))
	for i in range(len(logMerrL)):
		if logMerrH[i]==999:
			logMerr[i]=logMerrL[i]
		else:
			logMerr[i]=np.max([logMerrH[i],logMerrL[i]])
#Calculate the total vertical error 
	extra_err=a2*logsigerr
	err=(extra_err**2+logMerr**2)**(1./2.)
#Do the fit
	S=np.sum(1./extra_err**2)	
	sigmax=np.sum(logsigma/extra_err**2)
	sigmay=np.sum(logM/extra_err**2)
	sigmax2=np.sum(logsigma*logsigma/extra_err**2)
	sigmaxy=np.sum(logsigma*logM/extra_err**2)

	a1=(sigmay*sigmax2-sigmax*sigmaxy)/(S*sigmax2-(sigmax)**2)
	a2=(S*sigmaxy-sigmay*sigmax)/(S*sigmax2-(sigmax)**2)
	print a1,a2	
	y=a1+a2*x
#and plot it
	plt.clf()
	plt.errorbar(logsigma,logM,logMerr,logsigerr,c='k',fmt='o')
	a,=plt.plot(x,y,'k-')
	b,=plt.plot(x,y2,'k--')
	plt.xlim(1.4,2.6)
	plt.ylim(4,9)
	plt.xlabel(r"$\log\left(\frac{\sigma_*}{km/s}\right)$",fontsize=20)
	plt.ylabel(r'$\log\left(\frac{M_{BH}}{M_{\odot}}\right)$',fontsize=20)
	plt.legend([a,b],['My Fit','Greene & Ho Fit'],loc=2)
	plt.savefig('onec.pdf')

