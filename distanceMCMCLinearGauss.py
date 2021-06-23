import numpy as np
import scipy
from progressbar import * 
from scipy import special

#first linear increase, then gaussion,
# to eliminae the effect of greadually incrase AG

# in total 6 parameters are need 

"""
Abundent truncated gaussion, 
"""

nominator=np.sqrt(2) #*sigma 
logsqrt2pi=0.5*np.log(2*np.pi) 
agmin= 0 #np.min(dataAG) 
agmax= 3.609 #np.max(dataAG)
#AGErrorSquare=AGError*AGError
sqrtPi2i=2./np.sqrt(2*np.pi) 
 

logsqrtPi2i=np.log(sqrtPi2i)


def calProb6p(disCloud, a, b, imu1sigma,mu2,imu2sigma,dataDis, dataAG,disError, AGError,AGErrorSquare,Nstars):
	##calculate the probablitydisSigma
	#trueError=np.sqrt(disError**2+disSigma**2)
	#trueError= #np.sqrt(disError**2+disSigma**2)

	weight=scipy.stats.norm.cdf(disCloud,dataDis,disError)
	
	#weight=1-scipy.stats.norm.cdf(dataDis,disCloud, trueError)
 
	#fore ground
	mu1=a*dataDis+b
	errorVar= 1./imu1sigma**2+ AGErrorSquare+a**2*disError**2  #AGErrorSquare
	errorVar_i=1./errorVar
	convolveSTDi=  np.sqrt(errorVar_i) 
	#a=  special.erf((agmax-mu1)/nominator*convolveSTDi)+ special.erf((mu1-agmin)/nominator*convolveSTDi)  
	
	#common factor,  np.exp(-0.5*(dataAG-mu1)**2/errorVar)*convolveSTDi*sqrtPi2i/a
	
	#pForeground= np.exp(-0.5*(dataAG-mu1)**2/errorVar)*convolveSTDi*sqrtPi2i/a
	pForeground=  0.5*np.exp(-0.5*np.square(dataAG-mu1)*errorVar_i)*convolveSTDi  #*convolveSTDi #/a

	
	#pForeground=	  np.sum(-0.5*(dataAG-mu1)**2/errorVar - np.log(a)-0.5*np.log( errorVar)  )
 
	#fore ground
	errorVar2= 1./imu2sigma**2 +AGErrorSquare
	errorVar2_i=1./errorVar2
	convolveSTDi2=  np.sqrt(errorVar2_i)
	a2=  special.erf((agmax-mu2)/nominator*convolveSTDi2)+ special.erf((mu2-agmin)/nominator*convolveSTDi2)  
	
	
	pBackground= np.exp(-0.5*np.square(dataAG-mu2)*errorVar2_i)*convolveSTDi2/a2 #+logsqrtPi2i
	
	totalP= pBackground+(pForeground-pBackground)*weight #true norm
	
	return np.sum(np.log(totalP) )  +Nstars*logsqrtPi2i
 



def getNextValues6p( runSamples,indexPara):
	
	
	 
	RMSs=[50, 0.00002,0.1 , 0.3,0.3,0.3]
	
	
	
	currentValues=[ ]

	for j in range(len(runSamples)):
		currentValues.append(runSamples[j][-1])


	currentValues[indexPara]=currentValues[indexPara] +np.random.normal(0,RMSs[indexPara])
	
	return currentValues
 


def getDisAndErrorMCMCTrucatedGau6p(  dataDis, dataAG,disError, AGError, processID, returnSampleDic, sampleN=1000,burn_in=50,thinning=10,maxDisCut=1000. ):
	
	"""
	do MCMC in this single function
	"""

	
 
	
	#print min( dataDis )
	#print max( dataDis)
 
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
	
	last100pcAG=dataAG[-20:]
	
	
	mu20=np.mean(last100pcAG)
	
	isigma0=np.std(last100pcAG,ddof=1)


	#print dataAG[dataDis>maxDis-50]
	
	isigma0= 1./isigma0 

	aaa=0.5*np.log(2*np.pi)
	
	sqrt2pi=np.sqrt(2*np.pi)
	sqrt2=np.sqrt(2)
 
	agmin= 0 #np.min(dataAG)
	agmax= 3.609 #np.max(dataAG)
 
	
	

	modeGibbs=0

	disList=[]
	mu1List=[]
	mu2List=[]
	imu1sigmaList=[]
	imu2sigmaList=[]
	disSigmaList=[]


	disCloud0=0
	sumAGError=np.sum(-aaa-np.log(AGError))


 
	disCloud=np.random.uniform(minDis,maxDis-50) # maxDistanceSearch) #np.random.normal(distance0,distanceError) # #choice(pns)
	mu2=  np.random.exponential(mu20)   # np.random.uniform(0,3.609)  #abs( np.random.normal(mu20,0.45) )
	while mu2> agmax:
		mu2=  np.random.exponential(mu20)
	a=np.random.exponential(0.0001)
	b=np.random.normal(0,0.5)
 	imu1sigma=  np.random.exponential(2) 
	imu2sigma=  np.random.exponential(isigma0)
	mu1=  np.random.exponential(0.5)

 	p0=calProb6p(disCloud, a,b,imu1sigma,mu2,imu2sigma,dataDis, dataAG,disError, AGError,AGErrorSquare,Nstars)
	
	if np.isnan(p0):
		print "Wrong intial value!!!!!!!!!!!!!!!!!!!!!!!!!!"
		for i in range(1000):
			disCloud=np.random.uniform(minDis,maxDis) # maxDistanceSearch) #np.random.normal(distance0,distanceError) # #choice(pns)
			mu1=  np.random.exponential(0.5)

			mu2=  np.random.exponential(mu20)   # np.random.uniform(0,3.609)  #abs( np.random.normal(mu20,0.45) )
		 	imu1sigma=  np.random.exponential(2) 
			imu2sigma=  np.random.exponential(isigma0)
			
			a=np.random.exponential(0.0001)
			b=np.random.normal(-1,0.5)
			
 			p0=calProb6p(disCloud,  a,b,imu1sigma,mu2,imu2sigma,dataDis, dataAG,disError, AGError,AGErrorSquare,Nstars)
			if not np.isnan(p0):
				break
	 
	runSamples= [[disCloud],[a],[b],[imu1sigma],[mu2], [imu2sigma] ]

	widgets = ['MCMCSmapleDistance: ', Percentage(), ' ', Bar(marker='>',left='|',right='|'),
	           ' ', ETA(), ' ', FileTransferSpeed()] #see docs for other options
	
	pbar = ProgressBar(widgets=widgets, maxval=sampleN+burn_in+1)
	pbar.start() 

	searchMax=maxDis-50

	acceptN=0

	 

	for i in range(10000000):
		if i>5000 and np.mean(disList)==disList[-1] :
			
			print p0,p1,"Wrong, restart"
			print "Parameter values:",disCloud, mu1,imu1sigma,mu2,imu2sigma
				
			getDisAndErrorMCMCTrucatedGau6p(  dataDis, dataAG,disError, AGError, processID, returnSampleDic, sampleN=sampleN,burn_in=burn_in,maxDisCut=  maxDisCut )
			return

		paraJ=0
		while paraJ<6:
			
			valueCand= getNextValues6p(runSamples,paraJ)
			disCloud,  a,b,imu1sigma,mu2,imu2sigma=valueCand
			
			if imu1sigma <0 or imu2sigma<0   or disCloud<minDis or disCloud>searchMax or mu2<0  or mu2>agmax :
				#paraJ=paraJ-1
				continue
			p1=	 calProb6p(disCloud, a,b,imu1sigma,mu2,imu2sigma,dataDis, dataAG,disError, AGError,AGErrorSquare,Nstars)
 
			randomR=np.random.uniform(0,1)
			
			
			
			if p1>=p0 or p1-p0>np.log(randomR):
				#print disCloud,  imu1sigma,mu2,imu2sigma,"--->",processID

				p0=p1;
				runSamples[paraJ].append( valueCand[paraJ] )
				acceptN=acceptN+1
			else:
				runSamples[paraJ].append( runSamples[paraJ][-1] )
			paraJ=paraJ+1
		
		runSamples=[ runSamples[0][-1:],runSamples[1][-1:], runSamples[2][-1:], runSamples[3][-1:],runSamples[4][-1:],runSamples[5][-1:]   ]


		if i%thinning==0:
 

			disList.append(runSamples[0][-1] )
			mu1List.append(runSamples[1][-1] )
			imu1sigmaList.append(runSamples[2][-1]  )
			mu2List.append(runSamples[3][-1] )
			imu2sigmaList.append(runSamples[4][-1]  )
 
			#print processID,runSamples[0][-1] ,runSamples[1][-1] ,runSamples[2][-1] ,runSamples[3][-1] ,runSamples[4][-1] ,runSamples[5][-1] 

			pbar.update(len(disList)) #this adds a little symbol at each iteration
			if len(disList)>burn_in+sampleN:
			
				break
					
 

	pbar.finish()
	
	
	print "The accept rate is", acceptN/1./i/5
 

	if 1: #test normal samplig
		
		disArray=np.array(disList[burn_in:-1]) 
		
		#for mu1 t	
		mu1Array=np.array(mu1List[burn_in:-1])
		#print "The modeled mu1 is ",np.mean(mu1ArrayT)
	 

		#mu1Array=np.array(disSigmaList[burn_in:])


			
		mu2Array=np.array(mu2List[burn_in:-1])
			
		mu1SigmaArray=1./np.array( imu1sigmaList[burn_in:-1])
		mu2SigmaArray=1./np.array( imu2sigmaList[burn_in:-1])
			
			
		#print "Testing correlation time"
	
		#calAuto(disArray)
		#calAuto(mu1Array)
		#calAuto(mu1SigmaArray)
		#calAuto(mu2Array)	
		#calAuto(mu2SigmaArray)	
			
			
		returnSampleDic[processID]= [disArray, mu1Array, mu1SigmaArray,mu2Array,mu2SigmaArray] 
		return 
	

