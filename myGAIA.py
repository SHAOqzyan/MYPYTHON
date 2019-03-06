from astropy.table import   Table
import numpy as np
import os
#os.environ['QT_QPA_PLATFORM']='offscreen'
#gaiaDIS.py

import matplotlib as mpl
mpl.use('Agg')
from gaiaTB import  GAIATB
from progressbar import * 
import pymc3 as pm
from astropy.table import Column
from scipy.stats import norm,expon
from scipy import special
import multiprocessing
import scipy
import corner
from astropy.wcs import WCS
import pywcsgrid2
from matplotlib import rc
import os

from astropy.io import fits

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm


nominator=np.sqrt(2) #*sigma
logsqrt2pi=0.5*np.log(2*np.pi)
agmin= 0 #np.min(dataAG)
agmax= 3.609 #np.max(dataAG)
#AGErrorSquare=AGError*AGError
sqrtPi2i=2./np.sqrt(2*np.pi) 


logsqrtPi2i=np.log(sqrtPi2i)


def calProb(disCloud,  mu1,imu1sigma,mu2,imu2sigma,dataDis, dataAG,disError, AGError,AGErrorSquare,Nstars):
	##calculate the probablitydisSigma
	#trueError=np.sqrt(disError**2+disSigma**2)
	#trueError= #np.sqrt(disError**2+disSigma**2)

	weight=scipy.stats.norm.cdf(disCloud,dataDis,disError)
	
	#weight=1-scipy.stats.norm.cdf(dataDis,disCloud, trueError)
 
	#fore ground
	errorVar= 1./imu1sigma**2+AGErrorSquare
	errorVar_i=1./errorVar
	convolveSTDi=  np.sqrt(errorVar_i)
	a=  special.erf((agmax-mu1)/nominator*convolveSTDi)+ special.erf((mu1-agmin)/nominator*convolveSTDi)  
	
	#common factor,  np.exp(-0.5*(dataAG-mu1)**2/errorVar)*convolveSTDi*sqrtPi2i/a
	
	#pForeground= np.exp(-0.5*(dataAG-mu1)**2/errorVar)*convolveSTDi*sqrtPi2i/a
	pForeground= np.exp(-0.5*np.square(dataAG-mu1)*errorVar_i)*convolveSTDi/a

	
	#pForeground=	  np.sum(-0.5*(dataAG-mu1)**2/errorVar - np.log(a)-0.5*np.log( errorVar)  )
 
	#fore ground
	errorVar2= 1./imu2sigma**2 +AGErrorSquare
	errorVar2_i=1./errorVar2
	convolveSTDi2=  np.sqrt(errorVar2_i)
	a2=  special.erf((agmax-mu2)/nominator*convolveSTDi2)+ special.erf((mu2-agmin)/nominator*convolveSTDi2)  
	
	
	pBackground= np.exp(-0.5*np.square(dataAG-mu2)*errorVar2_i)*convolveSTDi2/a2
	
	totalP= pBackground+(pForeground-pBackground)*weight #true norm
	
	return np.sum(np.log(totalP) )+Nstars*logsqrtPi2i
 



def getNextValues( runSamples,indexPara):
	
 
	 
	RMSs=[50,0.3,0.3,0.3,0.3]
	
	
	
	currentValues=[ ]

	for j in range(len(runSamples)):
		currentValues.append(runSamples[j][-1])


	currentValues[indexPara]=currentValues[indexPara] +np.random.normal(0,RMSs[indexPara])
	
	return currentValues

def getDisAndErrorMCMCTrucatedGau5p(  dataDis, dataAG,disError, AGError, processID, returnSampleDic,   sampleN=1000,burn_in=50,thinning=10):
	
	"""
	do MCMC in this single function
	"""
	#first dataDis, dataAG,disError, AGError
	#print "Using truncated Gaussion"
	AGErrorSquare=AGError*AGError
	dataAGSquare=dataAG*dataAG
	
	Nstars=  len(dataAG )*1.0  #,dataAGSquare):
	np.random.seed()
	#print "Calculating distances with MCMC...total smaple number is {}....".format(sampleN)
 
	minDis=int(np.min(dataDis))
 
	maxDis=int(np.max(dataDis))  #to avoid touch the edge
 
	last100pcAG=dataAG[dataDis>maxDis-50]
	mu20=np.mean(last100pcAG)
	isigma0=np.std(last100pcAG,ddof=1)

	
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
 	imu1sigma=  np.random.exponential(2) 
	imu2sigma=  np.random.exponential(isigma0)
	mu1=  np.random.exponential(0.5)

 	p0=calProb(disCloud, mu1,imu1sigma,mu2,imu2sigma,dataDis, dataAG,disError, AGError,AGErrorSquare,Nstars)
	
	if np.isnan(p0):
		print "Wrong intial value!!!!!!!!!!!!!!!!!!!!!!!!!!"
		for i in range(1000):
			disCloud=np.random.uniform(minDis,maxDis) # maxDistanceSearch) #np.random.normal(distance0,distanceError) # #choice(pns)
			mu1=  np.random.exponential(0.5)

			mu2=  np.random.exponential(mu20)   # np.random.uniform(0,3.609)  #abs( np.random.normal(mu20,0.45) )
		 	imu1sigma=  np.random.exponential(2) 
			imu2sigma=  np.random.exponential(isigma0)
 			p0=calProb(disCloud,  mu1,imu1sigma,mu2,imu2sigma,dataDis, dataAG,disError, AGError,AGErrorSquare,Nstars)
			if not np.isnan(p0):
				break
	 
	runSamples= [[disCloud],[mu1],[imu1sigma],[mu2], [imu2sigma] ]

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
				
			getDisAndErrorMCMCTrucatedGau5p(  dataDis, dataAG,disError, AGError, processID, returnSampleDic, sampleN=sampleN,burn_in=burn_in)
			return

		paraJ=0
		while paraJ<5:
			
			valueCand= getNextValues(runSamples,paraJ)
			disCloud,  mu1,imu1sigma,mu2,imu2sigma=valueCand
			
			if imu1sigma <0 or imu2sigma<0 or mu2<0 or disCloud<minDis or disCloud>searchMax or mu1<0:
				#paraJ=paraJ-1
				continue
			p1=	 calProb(disCloud, mu1,imu1sigma,mu2,imu2sigma,dataDis, dataAG,disError, AGError,AGErrorSquare,Nstars)
 
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
			imu1sigmaList.append(runSamples[2][-1]  )
			mu2List.append(runSamples[3][-1] )
			imu2sigmaList.append(runSamples[4][-1]  )
 


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
	


 
 


#this  module is used to calculate the distance with Gaia ag extinction

class GAIADIS:
	#name="GSF_Gaia2_all"
	
	
	coint="coint"
	
	gaiaDo=GAIATB()
	
	nChains=8
	
	chainSample= 500 #1000
	burn_in=50
	
	thinning=10 #get on sample every XX stars.
	

	
	GAIA_source_id ="source_id"


	GAIA_parallax ="parallax"

	GAIA_parallax_err ="parallax_err"
	
	GAIA_l ="l"
	GAIA_b ="b"

	GAIA_distance ="distance"
	GAIA_distanceError ="distance_err"


	GAIA_a_g_val ="a_g_val"

	GAIA_a_g_percentile_lower ="a_g_percentile_lower"
	GAIA_a_g_percentile_upper ="a_g_percentile_upper"

	agError ="agError"

	relative_error ="relative_error"

	def __init__(self,sourceName ):

		self.sourceName=sourceName
	
	def calDisWith5Ps(self,oncloudStars):
		"""
		calculate the distance with oncloud stars, no figure will be provided, and only sample will be returned
		
		5 parameers are used 
		mu1 is assgned to be zero
		"""
		#pass
 
		
		dataDis, dataAG,disError, AGError=self.getDisAndAGFromTB(oncloudStars)

		rawData=[dataDis, dataAG,disError, AGError]

		nChain=self.nChains # =================you could use 8 chain, if your computer allows it
		print "Calculating {} chains and each chain has {} samples, thinned by {}.".format(self.nChains,self.chainSample,self.thinning)  

		procs = []
	
		manager = multiprocessing.Manager()
		returnSampleDic = manager.dict()

		# instantiating process with arguments
		for name in range(nChain):
			# print(name)
			print "Starting process ",name
			
			proc = multiprocessing.Process(target=getDisAndErrorMCMCTrucatedGau5p, args=(dataDis,dataAG,disError,AGError,name, returnSampleDic,self.chainSample,self.burn_in,self.thinning))
			procs.append(proc)
			proc.start()
 
		for proc in procs:
			proc.join()
 
		sampleArray=[]
		
	 	
	 	for j in range(len(returnSampleDic[0])):
			
	 		sampleArrayj=[]
			for i in  range(len( returnSampleDic) ):
				
				sampleArrayj.append(returnSampleDic[i][j])
			
	 		sampleArray.append( np.concatenate(sampleArrayj) )

		return sampleArray,rawData
 

	def getSmoothArrayWithParallaxByMean(self, list1,list2,list1Error=None,list2Error=None,smoothNumber=5):
		"""
		This function is used to sort list1 and list1Error, is parallax
		
		"""
		
		#radius=smoothNumber/2. #not a number ,but pc
		
		list1Copy=list1.copy()
 
		NewList2=[]
		
		minPara=list1.min()
		maxPara=list1.max()
		
		maxDis=1./minPara*1000
		minDis=1./maxPara*1000
		
		minDis=int(minDis)
		maxDis=int(maxDis)
		 
		newDisRange=np.arange(minDis,maxDis+smoothNumber,smoothNumber)
		newParaRange=1000./newDisRange 
		
 
		newDisList=[]
		newExtictionList=[]
		
		error1=[]
		error2=[]
		
 
		
		for countN,eachD in enumerate(newParaRange[0:-1]):
 
			endPara=eachD
			
 
			beginPara=newParaRange[countN+1]
			 
			cut1=list1[ (list1Copy>=beginPara) &  (list1Copy<=endPara)]
 
			cut2=list2[ (list1Copy>=beginPara) &  (list1Copy<=endPara)]

			if len(cut1)==0:
				continue

			cut1Err=list1Error[ (list1Copy>=beginPara) &  (list1Copy<=endPara)]  #this is the error of distanes
			cut2Err=list2Error[ (list1Copy>=beginPara) &  (list1Copy<=endPara)]  #this is the error of distanes

			weight1=1./cut1Err**2/(np.sum(1./cut1Err**2))
			weight2=1./cut2Err**2/(np.sum(1./cut2Err**2))

	
			newDisList.append(1./np.average(cut1,weights=weight1)*1000 )
			
			newExtictionList.append( np.average(cut2,weights=weight2) )
			
			
			
  
		return np.array(newDisList),np.array(newExtictionList),None,None





	def getDisArray(self,lRange,bRange,fitsName, signalLevel,noiseLevel, saveTBName, lowerDisCut=1,upperDisCut=2000,paraErrorCut=0.2, newQuery=False ):
		"""
		
		
		calculate the distance with MCMC
		
		Parameters:
		
		lRange: the galactic longtitude range of the molecular cloud
		lRange: the galactic latitude range of the molecular cloud
		
		fitsName: the background fits file, which is used to classify Gaia stars
		signalLevel,noiseLevel: the signal and noise level cut
		
		saveTBName: the name of TB file, which is saved in the current path.
		
		lwoerDisCut,upperDisCut: the distance range we care about #unit pc
		
		paraErrorCut: the error limit of parallax, default value is 20%
		
		newQuery: upate the saveTBName, otherwise we would use the saved file to avoiding querying the database every time.
		

		"""
		
		
		#read fits file
		

		lowerPara=1./upperDisCut*1000
		upperPara=1./lowerDisCut*1000
		
		
		if newQuery: #save the table
			gaiaOnCloudStars=self.gaiaDo.getByLBRange(lRange,bRange,lowerPara=lowerPara,upperPara=upperPara, paraError=paraErrorCut)
			os.system("rm "+saveTBName)
			gaiaOnCloudStars.write(saveTBName)
			
			
		gaiaAllStars=Table.read(saveTBName)
		
		gaiaAllStars=self.assignOnsourceGaia(gaiaAllStars,fitsName)
		
		
		gaiaAllStars.add_index(self.coint)
		
		
		
		
		#split the gaia stars into on and off stars
		#on-cloud stars, has high emission and off-cloud stars has low emission
		# you have custom this divided levels if you use your own CO fits.
		gaiaOnCloudStars=gaiaAllStars.loc[self.coint,signalLevel:]
		gaiaOffCloudStars=gaiaAllStars.loc[self.coint, :noiseLevel]


		gaiaOnCloudStars.remove_indices(self.coint)


		gaiaOnCloudStars.sort( self.GAIA_parallax )

		gaiaOnCloudStars.reverse() #
		
		
		gaiaOffCloudStars.remove_indices(self.coint)
		gaiaOffCloudStars.sort( self.GAIA_parallax )
		gaiaOffCloudStars.reverse()



		p_coParallax= gaiaOnCloudStars[self.GAIA_parallax] 
		p_coExtention=gaiaOnCloudStars[self.GAIA_a_g_val]
		p_coParallaxErr=gaiaOnCloudStars[self.GAIA_parallax_err] 
		p_coExtentionErr= gaiaOnCloudStars[self.agError]*p_coExtention
		
		p_nocoParallax= gaiaOffCloudStars[self.GAIA_parallax]  #pc
		p_nocoExtention=gaiaOffCloudStars[self.GAIA_a_g_val]
		p_nocoParallaxErr=gaiaOffCloudStars[self.GAIA_parallax_err]  #pc
		p_nocoExtentionErr= gaiaOffCloudStars[self.agError]*p_nocoExtention


		
		smoothNumber=5 #pc
		#get somoothed on and off cloud star
		smsortCOPara,smsortCOExtention,smweigtCOPara,smweightCOExtention=self.getSmoothArrayWithParallaxByMean(p_coParallax.data,p_coExtention.data,list1Error=p_coParallaxErr.data,list2Error=p_coExtentionErr.data,smoothNumber=smoothNumber)
		#smsortCOPara,smsortCOExtention,smweigtCOPara,smweightCOExtention=self.getSmoothArrayWithParallax(p_coParallax.data,p_coExtention.data,list1Error=p_coParallaxErr.data,list2Error=p_coExtentionErr.data,smoothNumber=smoothNumber)

		
		drawNOCOPara,drawNOCOAG,smweightNOCOPara,smweightNOCOExtention=self.getSmoothArrayWithParallaxByMean(p_nocoParallax.data,\
		p_nocoExtention.data,list1Error=p_nocoParallaxErr.data,list2Error=p_nocoExtentionErr.data,smoothNumber=smoothNumber)
		
		
		
		
 
		
		sampleArray,rawData=self.calDisWith5Ps(gaiaOnCloudStars )


		dataDis, dataAG,disError, AGError=rawData 
		

		distance,lowerDis,upperDis=self.getDisAndHPD(sampleArray[0])
		agmu1,loweragmu1,upperagmu1=self.getmuAndHPD(sampleArray[1])
		agmu1Sigma,loweragmu1Sigma,upperagmu1Sigma=self.getmuAndHPD(sampleArray[2])
		agmu2,loweragmu2,upperagmu2=self.getmuAndHPD(sampleArray[3])
		agmu2Sigma,loweragmu2Sigma,upperagmu2Sigma=self.getmuAndHPD(sampleArray[4])


		print "The distance is:",distance 

		#draw figuers

		fig, axs = plt.subplots(ncols=1, nrows=2,figsize=(16.5,8) )  
		rc('text', usetex=True )
		rc('text.latex',  preamble=r'\usepackage{upgreek}')

		for ax in  axs  :
			ax.remove()

		##########draw corner maps
		nn=1
		
		for i in [0,1,2,3,4 ]:
			plt.subplot2grid((10*nn,10*nn), (i*nn*2, 0) ,rowspan=2*nn,colspan=nn)
			
			plt.subplot2grid((10*nn,10*nn), (i*nn*2,  nn),rowspan=2*nn,colspan=nn)
			plt.subplot2grid((10*nn,10*nn), (i*nn*2,  2*nn),rowspan=2*nn,colspan=nn)
			
			plt.subplot2grid((10*nn,10*nn), (i*nn*2,  3*nn),rowspan=2*nn,colspan=nn)
			plt.subplot2grid((10*nn,10*nn), (i*nn*2,  4*nn),rowspan=2*nn,colspan=nn)



		sampleArray=np.array(sampleArray)

		nx,ny=sampleArray.shape
			
		sampleArray=sampleArray.transpose()
		sampleArray=sampleArray.reshape(ny,nx)

		figureCorner = corner.corner(sampleArray,fig=fig,  show_titles=True )

		ndim=5
 		axes = np.array(figureCorner.axes).reshape((ndim, ndim))
		
		
 
		meanValues = [distance,agmu1, agmu1Sigma,agmu2,agmu2Sigma] #np.mean(sampleArray, axis=0)
		lowerPHD90 = [lowerDis,loweragmu1, loweragmu1Sigma,loweragmu2,loweragmu2Sigma] #np.mean(sampleArray, axis=0)
		upperPHD90 = [upperDis, upperagmu1,upperagmu1Sigma,upperagmu2,upperagmu2Sigma] #np.mean(sampleArray, axis=0)

		#agmu1Sigma,loweragmu1Sigma,upperagmu1Sigma

		titleDis=r'$ \rm D  = {}_{{-{}}}^{{+{}}}$ pc'.format(distance,lowerDis,upperDis)
		
		titleMu1=r'$\upmu_2 =  {:.3f}_{{-{:.3f}}}^{{+{:.3f}}}$ mag'.format(agmu1,loweragmu1,upperagmu1)

		titleMu1Sigma=r'$\upsigma_1 =  {:.3f}_{{-{:.3f}}}^{{+{:.3f}}}$ mag'.format(agmu1Sigma,loweragmu1Sigma,upperagmu1Sigma)

		titleMu2=r'$\upmu_2 =  {:.3f}_{{-{:.3f}}}^{{+{:.3f}}}$ mag'.format(agmu2,loweragmu2,upperagmu2)
		
		titleMu2Sigma=r'$ \upsigma_2 =  {:.3f}_{{-{:.3f}}}^{{+{:.3f}}}$ mag'.format(agmu2Sigma,loweragmu2Sigma,upperagmu2Sigma)

		labels=[titleDis,titleMu1, titleMu1Sigma,titleMu2,titleMu2Sigma]



		# Loop over the diagonal
		for i in range(ndim):
			axCorner = axes[i, i]
			#pass
			#ax.axvline(meanValues[i], color="black" )
			#axCorner.title("aaa")

			axCorner.set_title(labels[i],fontsize=12)
			if i==0  :
				axCorner.set_title(labels[i],fontsize=14)

			axCorner.axvline(meanValues[i], color="black" ,linestyle='-',linewidth=1.5)
			axCorner.axvline(meanValues[i]-lowerPHD90[i], color="black" ,linestyle='--',linewidth=1.5)
			axCorner.axvline(meanValues[i]+upperPHD90[i], color="black" ,linestyle='--',linewidth=1.5)
			
 
 

		#draw gia Rawstars
		ax=plt.subplot2grid((10*nn,10*nn), (0,4),colspan=6,rowspan=4)
 
		maxY= max(smsortCOExtention)+ (np.max(dataAG)-max(smsortCOExtention))*0.8

		#ax.scatter(sortNOCOPara,sortNOCOExtention,lw=0.3,facecolors='b',s=5, edgecolors='b',label="Off-cloud stars" )
		ax.fill_betweenx(y=[-0.5, maxY*0.9], x1=distance-lowerDis, x2=upperDis+distance, color='moccasin',lw=0.1 );

		ax.scatter(dataDis, dataAG, edgecolors='darkgray', facecolors='darkgray',s=2,label="Raw on-cloud stars") # raw data
		
		
		ax.scatter(drawNOCOPara,drawNOCOAG,lw=0.3,facecolors='b',s=8, edgecolors='b',label="Binned off-cloud stars" )

		ax.scatter(smsortCOPara, smsortCOExtention, edgecolors="g" ,facecolors='g',s=8,label="Binned on-cloud stars")


		ax.plot([distance,distance],[-0.5,maxY*0.9],lw=1,color="black")
		ax.text(distance-dataDis.max()/100.*5,maxY*0.91, r'${}_{{-{}}}^{{+{}}}$ pc'.format(distance,lowerDis,upperDis),fontsize=15)
		
		ax.plot([min(dataDis),distance],[agmu1,agmu1], 'r--',lw=1.5  ,dashes=(4, 2))
		ax.plot([distance,max(dataDis)],[agmu2,agmu2], 'r--',lw=1.5  ,dashes=(4, 2))




		#draw fits
		import pywcsgrid2
		from mpl_toolkits.axes_grid1.inset_locator import inset_axes
		import mpl_toolkits.axes_grid1.inset_locator as inset_locator
		


		backFISTHDR=fits.open(fitsName)

		backData=backFISTHDR[0].data
		backHead=backFISTHDR[0].header

		axBack=pywcsgrid2.subplot(224,  header= WCS(backHead))
		
		vmin= noiseLevel/2. 
		vmax= signalLevel*1.5
		axBack.imshow(backData, origin="lower",cmap="Greys",norm=LogNorm(vmin=vmin, vmax=vmax),interpolation=None)


		#cmap = plt.cm.winter
		#cmap.set_bad('white',1.)
		#axBack[WCS(backHead)].contour(cropData,  [IRASNoise,IRASSignal],cmap=cmap,vmin=IRASNoise*0.7,vmax=IRASSignal*1.1, linewidths=0.7)
		
		
		axBack["gal"].plot([lRange[0],lRange[0]],bRange,'b--',lw=0.8)
		axBack["gal"].plot([lRange[1],lRange[1]],bRange,'b--',lw=0.8)
		
		axBack["gal"].plot(lRange,[bRange[0],bRange[0]],'b--',lw=0.8)
		axBack["gal"].plot(lRange,[bRange[1],bRange[1]],'b--',lw=0.8)


 

		saveFigureName="{}_{}_{}_{}{}".format(self.sourceName,lowerDisCut,upperDisCut,paraErrorCut,'_extinctionGaiaAg.png')

		#plt.savefig(  self.sourceName+'_extinctionGaiaAg.pdf', bbox_inches="tight")
		plt.savefig( saveFigureName, bbox_inches="tight")

		return distance,lowerDis,upperDis







	def getDisAndAGFromTB(self,oncloudStars):
		print "Calculating distances and their errors...."
		dataDis= oncloudStars["parallax"]*0
		disError=dataDis*0
		dataAG=oncloudStars["a_g_val"]
		AGError= oncloudStars["agError"]* dataAG

		#if 1 :#:#calculate parallax with mcmc
		for i in range(len(dataDis)):
			para=oncloudStars[i]["parallax"]
			paraError=oncloudStars[i]["parallax_err"] #+0.1 #   #
			
			dA=1./np.random.normal(para,paraError,10000)*1000
			
			dataDis[i]=np.mean(dA)
			disError[i]=np.std(dA,ddof=1)
		#no bining
		
		return dataDis, dataAG,disError, AGError
 


	def getDisAndHPD(self,disSample):
		
		

		disSample=np.array(disSample)
		
		meanDis=np.mean(disSample)
		
		#disHPD=pm.stats.hpd(disSample,alpha=0.1)
		disHPD=pm.stats.hpd(disSample,alpha=0.05)

		
		print disHPD,"===95%HPD===="
		#print disSample[-10:]
		lowerDis=meanDis-min(disHPD)
		upperDis= max(disHPD)-meanDis
		
		meanDis,lowerDis,upperDis=map(round,[meanDis,lowerDis,upperDis])
		meanDis,lowerDis,upperDis=map(int,[meanDis,lowerDis,upperDis])
		
		return meanDis, lowerDis,upperDis 


	def getmuAndHPD(self,muSample):
		muSample=np.array(muSample)
		
		meanMu=np.mean(muSample)
		
		muHPD=pm.stats.hpd(muSample,alpha=0.1)
		lowerMu=meanMu-min(muHPD)
		upperMu= max(muHPD)-meanMu
 

		return np.round(meanMu,3), np.round(lowerMu,3),np.round(upperMu,3) 


	def assignOnsourceGaia(self,GaiaTB,bgFITS):
 
		GaiaTB= GaiaTB.copy()

		
		bgData,bgHead= fits.open(bgFITS)[0].data,fits.open(bgFITS)[0].header    #myFITS.readFITS(bgFITS)
		
 
		
		if len(bgData.shape)==3:
			
			bgData=bgData[0]
			
			bgHead["NAXIS"]=2
			del bgHead["CRPIX3"] 
			del bgHead["CDELT3"] 
			
			del bgHead["CRVAL3"] 
			del bgHead["CTYPE3"] 
 
		addCol=Column(name=self.coint,data=np.zeros(len(GaiaTB)))
		GaiaTB.add_column(addCol)
 
		bgWCS=WCS(bgHead)
		
		Ls=GaiaTB["l"]
		Bs=GaiaTB["b"]
 
		Xs,Ys=bgWCS.all_world2pix(Ls,Bs,0)

 		Xs=map(round,Xs)
 		Ys=map(round,Ys)
 		Xs=map(int,Xs)
 		Ys=map(int,Ys)
 
		#[y,x]

		maxY,maxX=bgData.shape


		for i in range(len(Xs)):
			
			if Ys[i] >=maxY-1 or  Xs[i]>=maxX-1:
				GaiaTB[i][self.coint]= -100 #self.getCointByLB(bgData, WCS(bgHead), eachStar["l"],eachStar["b"])

			else:
			
				GaiaTB[i][self.coint]= bgData[Ys[i]][Xs[i]]  #self.getCointByLB(bgData, WCS(bgHead), eachStar["l"],eachStar["b"])



		return  GaiaTB 


	def calWithFITS(self,fitsName,Lrange,Brange,noise,signalLevel,noiseLevel):
		
		"""
		the noiseLevel is only used to compare, not for calculation
		"""
		#selectColName

		#first select Gaia by  LB range
		
		gaiaAllSourceStars=self.gaiaDo.getByLBRange(Lrange,Brange,upperPara= 1/2.2,paraError=0.2) # error<20%
 		gaiaAllSourceStars=self.assignOnsourceGaia(gaiaAllSourceStars,fitsName)

		print len(gaiaAllSourceStars)
		gaiaAllSourceStars.add_index(self.coint)
		gaiaOnSourceStars=gaiaAllSourceStars.loc[self.coint,noise*signalLevel:]
			
		print "The total number of on-cloud gaia stars: ", len(gaiaOnSourceStars)
		gaiaOffSourceStars=gaiaAllSourceStars.loc[self.coint,:noise*noiseLevel]
	
		dataDis, dataAG,disError, AGError=self.getDisAndAGFromTB(gaiaOnSourceStars)
		dataDisOff, dataAGOff,disErrorOff, AGErrorOff=self.getDisAndAGFromTB(gaiaOffSourceStars)

 
		sampleArray=getDisAndErrorMCMCTrucatedGau6p(dataDis, dataAG,disError, AGError)
		#	return  [disArray,disSigmaArray,mu1Array,mu1SigmaArray,mu2Array,mu2SigmaArray] #disEqual,errorEqual,round()
		

		
		distance,lowerDis,upperDis=self.getDisAndHPD(sampleArray[0])
		
		disError,lowerDisError,upperDisError=self.getDisAndHPD(sampleArray[1])
		
		agmu1,loweragmu1,upperagmu1=self.getmuAndHPD(sampleArray[2])

		
		
		agmu1Sigma,loweragmu1Sigma,upperagmu1Sigma=self.getmuAndHPD(sampleArray[3])
		
		
		agmu2,loweragmu2,upperagmu2=self.getmuAndHPD(sampleArray[4])
		agmu2Sigma,loweragmu2Sigma,upperagmu2Sigma=self.getmuAndHPD(sampleArray[5])
		

		if 1: #draw
			
			#fig, ax = plt.subplots(figsize=(8, 6))
			#fig, axs = plt.subplots(nrows=1, ncols=2,  figsize=(8, 6))
			bgData,bgHead=myFITS.readFITS(fitsName)

			axFITS=pywcsgrid2.subplot(122,header=WCS(bgHead) )
			
			axStar=plt.subplot(121 )

 
			
			axStar.scatter(dataDis, dataAG,s=3)
			#ax.scatter(dataDisOff, dataAGOff)
			axFITS.imshow(bgData[0],origin='lower',interpolation=None)
			
			
			axFITS["gal"].plot([Lrange[0],Lrange[0]],Brange,'b--',lw=0.8)
			axFITS["gal"].plot([Lrange[1],Lrange[1]],Brange,'b--',lw=0.8)
			
			axFITS["gal"].plot(Lrange,[Brange[0],Brange[0]],'b--',lw=0.8)
			axFITS["gal"].plot(Lrange,[Brange[1],Brange[1]],'b--',lw=0.8)



			#crop fits
			
			outCropFITS= " FTSIDiscrop.fits"
			myFITS.cropFITS2D(fitsName,outFITS=outCropFITS,Lrange=Lrange,Brange=Brange,overWrite=True)
			
			cropData,cropHead=myFITS.readFITS(outCropFITS)
			cmap = plt.cm.winter
			cmap.set_bad('white',1.)
			axFITS[WCS(cropHead)].contour(cropData,  [noise*noiseLevel,noise*signalLevel],cmap=cmap,vmin=noise*noiseLevel*0.7,vmax=noise*signalLevel*1.1, linewidths=0.1)
 
			
			plt.savefig(  'distance.png', bbox_inches="tight",dpi=300)


	def ZZZ(self):
		pass

	#def testDis



