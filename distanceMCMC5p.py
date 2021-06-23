import numpy as np
import scipy
from progressbar import * 
from scipy import special

nominator=np.sqrt(2) #*sigma
logsqrt2pi=0.5*np.log(2*np.pi)
agmin= 0 #np.min(dataAG)
agmax= 3.609 #np.max(dataAG)
#AGErrorSquare=AGError*AGError
sqrtPi2i=2./np.sqrt(2*np.pi) 


logsqrtPi2i=np.log(sqrtPi2i)


def calProb(disCloud,  mu1,mu1sigma,mu2,mu2sigma,dataDis, dataAG,disError, AGError,AGErrorSquare,Nstars):
	##calculate the probablitydisSigma
	#trueError=np.sqrt(disError**2+disSigma**2)
	#trueError= #np.sqrt(disError**2+disSigma**2)

	weight=scipy.stats.norm.cdf(disCloud,dataDis,disError)
	
	#weight=1-scipy.stats.norm.cdf(dataDis,disCloud, trueError)
 
	#fore ground
	errorVar=   mu1sigma**2+AGErrorSquare
	errorVar_i=1./errorVar
	convolveSTDi=  np.sqrt(errorVar_i)
	a=  special.erf((agmax-mu1)/nominator*convolveSTDi)+ special.erf((mu1-agmin)/nominator*convolveSTDi)  
	
	#common factor,  np.exp(-0.5*(dataAG-mu1)**2/errorVar)*convolveSTDi*sqrtPi2i/a
	
	#pForeground= np.exp(-0.5*(dataAG-mu1)**2/errorVar)*convolveSTDi*sqrtPi2i/a
	pForeground= np.exp(-0.5*np.square(dataAG-mu1)*errorVar_i)*convolveSTDi/a

	
	#pForeground=	  np.sum(-0.5*(dataAG-mu1)**2/errorVar - np.log(a)-0.5*np.log( errorVar)  )
 
	#fore ground
	errorVar2=  mu2sigma**2 +AGErrorSquare
	errorVar2_i=1./errorVar2
	convolveSTDi2=  np.sqrt(errorVar2_i)
	a2=  special.erf((agmax-mu2)/nominator*convolveSTDi2)+ special.erf((mu2-agmin)/nominator*convolveSTDi2)  
	
	
	pBackground= np.exp(-0.5*np.square(dataAG-mu2)*errorVar2_i)*convolveSTDi2/a2
	
	#totalP= pBackground+(pForeground-pBackground)*weight #true norm
	
	totalP= pBackground*(1-weight)+ pForeground*weight #true norm
	
	#print np.sum(np.log(totalP) )  +Nstars*logsqrtPi2i
	return np.sum(np.log(totalP) )  +Nstars*logsqrtPi2i
 



def getNextValues( runSamples,indexPara):
	
 
	 
	RMSs=[100,0.5,0.5,0.5,0.5]
	
	
	
	currentValues=[ ]

	for j in range(len(runSamples)):
		currentValues.append(runSamples[j][-1])


	currentValues[indexPara]=currentValues[indexPara] +np.random.normal(0,RMSs[indexPara])
	
	return currentValues
 




def getDisAndErrorMCMCTrucatedGau5p(  dataDis, dataAG,disError, AGError, processID, returnSampleDic, sampleN=1000,burn_in=50,thinning=10,maxDisCut=1000. ):
	
	"""
	do MCMC in this single function
	"""

	#print "Using Truncated Gaussion"
	
		
	#first dataDis, dataAG,disError, AGError
	#print "Using truncated Gaussion"
	AGErrorSquare=AGError*AGError
	dataAGSquare=dataAG*dataAG
	
	Nstars=  len(dataAG )*1.0  #,dataAGSquare):
	np.random.seed()
	#print "Calculating distances with MCMC...total smaple number is {}....".format(sampleN)
 
	minDis=int(np.min(dataDis))
 
	maxDis=int(np.max(dataDis))  #to avoid touch the edge
 
	#last100pcAG=dataAG[dataDis> maxDisCut-50]
	
	last100pcAG=dataAG[-50:]
	
	
	mu20=np.mean(last100pcAG)
	
	sigma0=np.std(last100pcAG,ddof=1)

 

	aaa=0.5*np.log(2*np.pi)
	
	sqrt2pi=np.sqrt(2*np.pi)
	sqrt2=np.sqrt(2)
 
	agmin= 0 #np.min(dataAG)
	agmax= 3.609 #np.max(dataAG)
 
 
	modeGibbs=0

	disList=[]
	mu1List=[]
	mu2List=[]
	mu1sigmaList=[]
	mu2sigmaList=[]
	disSigmaList=[]


	disCloud0=0
	sumAGError=np.sum(-aaa-np.log(AGError))


 
	disCloud=np.random.uniform(minDis,maxDis-50) # maxDistanceSearch) #np.random.normal(distance0,distanceError) # #choice(pns)
	mu2=  np.random.exponential(mu20)   # np.random.uniform(0,3.609)  #abs( np.random.normal(mu20,0.45) )

	#while mu2> agmax:
		#mu2=  np.random.exponential(mu20)
	
 
 	mu1sigma=  np.random.exponential(0.5) 
	mu2sigma=  np.random.exponential(sigma0)
	mu1=  np.random.exponential(0.5)

	# mu2 should be large than m1
	
	while mu2<mu1 or mu2>agmax:
		
		#print "What The fuck"
		#print mu2,mu1,mu20
		mu1=  np.random.exponential(0.5)

		mu2=  np.random.exponential(mu20)   # np.random.uniform(0,3.609)  #abs( np.random.normal(mu20,0.45) )

	
	

	
 	p0=calProb(disCloud, mu1,mu1sigma,mu2,mu2sigma,dataDis, dataAG,disError, AGError,AGErrorSquare,Nstars)
	

	
		
	#if np.isnan(p0):
		#pass
	 
	runSamples= [[disCloud],[mu1],[mu1sigma],[mu2], [mu2sigma] ]

	widgets = ['MCMCSmapleDistance: ', Percentage(), ' ', Bar(marker='>',left='|',right='|'),
	           ' ', ETA(), ' ', FileTransferSpeed()] #see docs for other options
	
	pbar = ProgressBar(widgets=widgets, maxval=sampleN+burn_in+1)
	pbar.start() 

	searchMax=maxDis-50

	acceptN=0

	 

	for i in range(10000000):
		if i>5000 and np.mean(disList)==disList[-1] :
			
			print p0,p1,"Wrong, restart"
			print "Parameter values:",disCloud, mu1,mu1sigma,mu2,mu2sigma
				
			getDisAndErrorMCMCTrucatedGau5p(  dataDis, dataAG,disError, AGError, processID, returnSampleDic, sampleN=sampleN,burn_in=burn_in,maxDisCut=  maxDisCut )
			return

		paraJ=0
		while paraJ<5:
			
			valueCand= getNextValues(runSamples,paraJ)
			disCloud,  mu1,mu1sigma,mu2,mu2sigma=valueCand
			
			if mu1sigma <0 or mu2sigma<0 or mu2<0 or disCloud<minDis or disCloud>searchMax or mu1<0  or mu1>mu2 or mu2> agmax    :
				#print mu1,mu2
				#paraJ=paraJ-1
				continue
			p1=	 calProb(disCloud, mu1,mu1sigma,mu2,mu2sigma,dataDis, dataAG,disError, AGError,AGErrorSquare,Nstars)
 
			randomR=np.random.uniform(0,1)
			
			
			
			if p1>=p0 or p1-p0>np.log(randomR):
				#print disCloud,  imu1sigma,mu2,imu2sigma,"--->",processID

				p0=p1;
				runSamples[paraJ].append( valueCand[paraJ] )
				acceptN=acceptN+1
			else:
				runSamples[paraJ].append( runSamples[paraJ][-1] )
			paraJ=paraJ+1
		
		runSamples=[ runSamples[0][-1:],runSamples[1][-1:], runSamples[2][-1:], runSamples[3][-1:],runSamples[4][-1:] ]


		if i%thinning==0:
 

			disList.append(runSamples[0][-1] )
			mu1List.append(runSamples[1][-1] )
			mu1sigmaList.append(runSamples[2][-1]  )
			mu2List.append(runSamples[3][-1] )
			mu2sigmaList.append(runSamples[4][-1]  )
 
			

			pbar.update(len(disList)) #this adds a little symbol at each iteration
			if len(disList)>burn_in+sampleN:
			
				break
					
 

	pbar.finish()
	
	
	print "The average accept rate is", acceptN/1./i/5
 

	if 1: #test normal samplig
		
		disArray=np.array(disList[burn_in:-1]) 
		
		#for mu1 t	
		mu1Array=np.array(mu1List[burn_in:-1])
		#print "The modeled mu1 is ",np.mean(mu1ArrayT)
	 

		#mu1Array=np.array(disSigmaList[burn_in:])


			
		mu2Array=np.array(mu2List[burn_in:-1])
			
		mu1SigmaArray=np.array( mu1sigmaList[burn_in:-1])
		mu2SigmaArray=np.array( mu2sigmaList[burn_in:-1])
			
			
		#print "Testing correlation time"
	
		#calAuto(disArray)
		#calAuto(mu1Array)
		#calAuto(mu1SigmaArray)
		#calAuto(mu2Array)	
		#calAuto(mu2SigmaArray)	
			
			
		returnSampleDic[processID]= [disArray, mu1Array, mu1SigmaArray,mu2Array,mu2SigmaArray] 
		return 
	

