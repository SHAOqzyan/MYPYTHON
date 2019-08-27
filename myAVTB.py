
#some table of gaia data base

import MySQLdb
from astropy.table import   Table
import numpy as np

from matplotlib.offsetbox import AnchoredText
from distanceMCMC5p import getDisAndErrorMCMCTrucatedGau5p
from sklearn.isotonic import IsotonicRegression
from distanceMCMCLinearGauss import getDisAndErrorMCMCTrucatedGau6p

from distanceMCMCTwoGauss  import getDisAndErrorMCMCG2



class gaiaAV:
	
 	gaiaTBFile=None 
 	gaiaTB=None
	GAIA_source_id ="source_id"


	GAIA_parallax ="parallax"

	GAIA_parallax_err ="parallax_err" 
	
	GAIA_l ="l"
	GAIA_b ="b"



	GAIA_a_g_val ="a_g_val"

	GAIA_a_g_percentile_lower ="a_g_percentile_lower"
	GAIA_a_g_percentile_upper ="a_g_percentile_upper"

	GAIA_av16 ="av16"
	GAIA_av50 ="av50"
	GAIA_av84 ="av84" 

	GAIA_dist16 ="dist16"
	GAIA_dist50  ="dist50"
	GAIA_dist84  ="dist84" 

	coint="coint"
	
	
	def __init__(self,gaiaFile):
		
		"""
		"""
 
		self.gaiaTBFile=  gaiaFile

 
		self.gaiaTB=Table.read( self.gaiaTBFile )
		self.gaiaTB.sort(self.GAIA_dist50)

		
		self.gaiaTB.add_index(self.GAIA_l)
		self.gaiaTB.add_index(self.GAIA_b)
		
		
	def assignOnsourceGaia(self,GaiaTB,bgFITS,colName=None):
 
		GaiaTB= GaiaTB.copy()

		
		bgData,bgHead= fits.open(bgFITS)[0].data,fits.open(bgFITS)[0].header    #myFITS.readFITS(bgFITS)
		
 

		if colName ==None:
			colName=self.coint
			#addCol=Column(name=self.coint,data=np.zeros(len(GaiaTB)))
 
			
		addCol=Column(name=colName,data=np.zeros(len(GaiaTB)))
			
		GaiaTB.add_column(addCol)
 
		bgWCS=WCS(bgHead,naxis=2)
		
		
		
		
		Ls=GaiaTB["l"]
		Bs=GaiaTB["b"]
 
		Xs,Ys=bgWCS.all_world2pix(Ls,Bs,0)

 		Xs=map(round,Xs)
 		Ys=map(round,Ys)
 		Xs=map(int,Xs)
 		Ys=map(int,Ys)
 
		#[y,x]

		sizeData=bgData.shape
		
		if len(sizeData)==3:
			bgData=bgData[0]
		if len(sizeData)==4:
			bgData=bgData[0]
			bgData=bgData[0]

		maxY,maxX=bgData.shape


		for i in range(len(Xs)):
			
			if Ys[i] >maxY-1 or  Xs[i]>maxX-1:
				GaiaTB[i][colName]= -1000 #self.getCointByLB(bgData, WCS(bgHead), eachStar["l"],eachStar["b"])

			else:
			
				GaiaTB[i][colName]= bgData[Ys[i]][Xs[i]]  #self.getCointByLB(bgData, WCS(bgHead), eachStar["l"],eachStar["b"])

				 

		return  GaiaTB 


	def getByLBRange(self,lRange,bRange,lowerDis=None,upperDis=None):
		
		newTB=self.gaiaTB.copy()
		
		newTB=newTB.loc[self.GAIA_l, min(lRange):max(lRange)]
		newTB=newTB.loc[self.GAIA_b, min(bRange):max(bRange)]

		return self.cutByDis(newTB,disLow=lowerDis,disUp=upperDis )
		
	def cutByDis(self,TB,disLow=None,disUp=None):
		
		newTB=TB.copy()
		
		if disLow==None and disUp==None:
			return newTB
			
			
		if disLow!=None:
		
			newTB=newTB[ newTB[self.GAIA_dist50] >disLow/1000. ] 
			
		
		if disUp!=None:

			newTB=newTB[ newTB[self.GAIA_dist50] <disUp/1000. ] 
			
 
		return newTB

	


	def getDisAndAVFromTB(self,TB):
		#print "Calculating distances and their errors...."
		
		distances=TB[self.GAIA_dist50]*1000.
		
		av=TB[self.GAIA_av50]

		disError= (TB[self.GAIA_dist84] -  TB[self.GAIA_dist16] )/2.*1000.

		avError= (TB[self.GAIA_av84] -  TB[self.GAIA_av16] )/2. 

		return distances, av,  disError,avError



	def getSmoothArrayWithTB(self,TB,smScale=10):
		"""
		the unnie of smoothScale is pc
		"""
		TB=TB.copy()
		#TB.sort(self.GAIA_dist50)


		newDisList=[]
		newExtictionList=[]
		
		error1=[]
		error2=[]

		distances, av,  disError,avError =self. getDisAndAVFromTB( TB)


 
 		
		minDis=int( min( distances) )
		maxDis=int( max(distances) ) +1
		
		newDisRange=np.arange(minDis,maxDis+smScale,smScale)

 
  
		for countN,beginD in enumerate(newDisRange[0:-1]):
 
 
			endD=newDisRange[countN+1]
 
			selectCriteria=   (distances>=beginD) &  (distances<=endD)


			cut1=distances[  selectCriteria]
 
			cut2=av[  selectCriteria ]

			if len(cut1)==0:
				continue


 

			cut1Err=disError[  selectCriteria]  #this is the error of distanes
			cut2Err=avError[  selectCriteria ]  #this is the error of distanes

			weight1=1./cut1Err**2/(np.sum(1./cut1Err**2))
			weight2=1./cut2Err**2/(np.sum(1./cut2Err**2))
			
	
			newDisList.append( np.average(cut1,weights=weight1) ) #distance in pc
			
			newExtictionList.append( np.average(cut2,weights=weight2) ) # ag average
			
		#no errors are needed
			
  
		return np.array(newDisList),np.array(newExtictionList),None,None



	def getBaseLine(self,offcloudStars):
		
		
		distances, av,  disError,avError =self. getDisAndAVFromTB( offcloudStars)
		
		
		
		baseLine = IsotonicRegression()

 
		#rawDataOff=[dataDisOff.copy(), dataAGOff.copy(),disErrorOff.copy(), AGErrorOff.copy()] 
		
		
		weight=1./avError**2
		
		#agWeight=  1./offcloudStars[self.agError] 

		baseLine.fit(distances,av,weight)
		
		return baseLine



	def getAVBase(self,dataDis,baseLine):
		
		
		#it is possible that the nearst star is closer than off sttars
		
		
		#print dataDis
		dataDis.mask=None
 		
 
		
		baseAG= baseLine.predict(  dataDis.data  )
		
		
		#print baseLine.predict([  baseLine.X_min_  ])[0], baseLine.predict([  baseLine.X_max_  ])[0]
		
		baseAG[ dataDis<  baseLine.X_min_  ]  = baseLine.predict([  baseLine.X_min_  ])[0]
		
		baseAG[ dataDis>  baseLine.X_max_  ]  =baseLine.predict([  baseLine.X_max_  ])[0]

		return baseAG











	def ZZZ(self):
		pass