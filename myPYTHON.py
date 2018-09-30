#!/usr/bin/python
#this File contains my own defined functions, which is very customed

# put this as the only version

#
from scipy.odr import *
import numpy as np

from astropy.io import fits

from matplotlib import pyplot as plt

from astropy.wcs import WCS
import os

import os.path
import numpy as np
import pywcsgrid2
from matplotlib import rc
from  mpl_toolkits.axes_grid1.axes_rgb import imshow_rgb
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes, mark_inset,inset_axes
from numpy import mean, sqrt, square



from subprocess import call 
import img_scale 

paperPath="/Users/qzyan/Desktop/UrsaMajorApJS/"



def model(p, x):
   a, b = p
   return a + b*x
	

def fitLinearODR(x,y,x_err,y_err):
	# Create a model for fitting.
 
	quad_model = Model(model)
	
	# Create a RealData object using our initiated data from above.
	data = RealData(x, y, sx=x_err, sy=y_err)
	
	# Set up ODR with the model and data.
	odr = ODR(data, quad_model, beta0=[1, 1 ])
 
	# Run the regression.
	out = odr.run()
	#out.pprint()
	
	return out









class myFITS:

	"This is fits is my own class of fits"
	fitsDATA=None
	fitsHeader=None
	fitsName=None
	fitsPath=None
	thisFITS=None


	def __init__(self,FITSname=None):
		if FITSname:#with fits input
			print "An object processing {} created!".format(FITSname)
			self.fitsDATA,self.fitsHeader=self.readFITS(FITSname)
			self.fitsPath,self.fitsName=os.path.split(FITSname)
			if self.fitsPath=="":
				self.fitsPath="./"

			self.thisFITS=FITSname


		else:
			print "No input FITS is provided for myFITS"
		#If there is a paramters input, automatically read this FITS
		#and some paramters will be automatically figured out


	@staticmethod
	def lineFit(X1, Y1, errorX, errorY):
		"""
		Fitting a line
		"""
		
		out1=fitLinearODR(X1, Y1, errorX, errorY)
		
		intercept,slope =out1.beta  #fitobj.params
		
		interceptError,std_err=out1.sd_beta 

		return  [slope,intercept,std_err,interceptError]

	@staticmethod
	def getVoxValue(data,dataWCS,LBV):
		"""
		This function is used to get a voxel value from a data cube
		
		v, km/s
		"""
		l,b,v=LBV
 
		x,y,z=data.shape 
		
		indexX,indexY,indexZ= dataWCS.wcs_world2pix(l,b,v*1000,0)
		
		indexX=int(round(indexX))
		indexY=int(round(indexY))
		indexZ=int(round(indexZ))
		
		#print indexX,indexY,indexZ
		
		#print x,y,z
		#print indexX,indexY,indexZ
		
		
		if indexX>=z or indexY>=y or indexZ>=z:
			return np.NaN
		
		return data[indexZ,indexY,indexX]
		
	@staticmethod
	def downTo2D(fitsFile,outPUT=None,overwrite=True):
		"""
		some files creaseted by Miriad is 3D, and there is only one axis in the third 
		
		this functino is used to transformd 3D to 2d 
		"""
		
		if str.strip(fitsFile)=='':
			print "Empty file, quitting...."
			return 

		if outPUT is None:
			writeName="to2D_"+fitsFile
 
		else:
			writeName=outPUT
		

		fileExist=os.path.isfile(writeName)
 
		
		if overwrite and fileExist:
			os.remove(writeName)
			
		if not overwrite and fileExist:
			print "File exists, quitting..."
			return 

		#read file
		hdu=fits.open(fitsFile)[0]		
		
		head=hdu.header
		#wmap=WCS(hdu.header)
		data=hdu.data
		
		data=data[0]
		
		del head["CRPIX3"]
		del head["CDELT3"]
		del head["CRVAL3"]
		del head["CTYPE3"] 
		
		#writefits

		
		
		
		fits.writeto(writeName,data,header=head)

	@staticmethod
	def roundToInt(someArray):

		"""
		round to nearest integer
		"""
		
		convertInt=np.rint(someArray)
		convertInt=map(int,convertInt)

		return convertInt
		


	@staticmethod
	def getLBVbyIndex(wcs,XYZ):
		X,Y,Z=XYZ
		X=int(round(X));
		Y=int(round(Y));
		Z=int(round(Z));
		
		l,b,v=wcs.all_pix2world(X,Y,Z,0)[0:3]
		
		return [l,b,v]

	def getVoxValueByIndex(self,data,XYZ):
		"""
		This function is used to get a voxel value from a data cube
		
		v, km/s
		
		"""
		X,Y,Z=XYZ
		X=int(round(X));
		Y=int(round(Y));
		Z=int(round(Z));
		
 
		x,y,z=data.shape 
		
		
		#print indexX,indexY,indexZ
		
		#print x,y,z
		#print indexX,indexY,indexZ
		
		
		if X>=z or Y>=y or Z>=z or X<0 or Y<0 or Z<0:
			return np.NaN
		
		return data[Z][Y][X]
		
 

	def smoothSpaceFITS(self,data,dataHeader,rawBeam,resultBeam,outPutName): # arcmin
		
		"""
		This function is used to smooth the fits in order to get the wanted resolution 
		
		The smooth kernel is calculatel withe Raw Beam and ResultBeam
		
		The unit of the Beam should be arcmin
		#At high Galactic latitude, The pixels are still equal size, however, the beam size with pixel, is distorted

		"""
 
		#dataNew=COdata.copy()
		
		
		
		# the pixel resolution
		
		#resolCO=0.125 #deg each pixe
		
		# convert RawBeam in arcmin to degree, calculate the size with pixels
		
		rawBeamDegree=rawBeam/60.
		resultBeamDegree=resultBeam/60.

		resolX=abs(dataHeader["CDELT1"]) #degree
		resolY=abs(dataHeader["CDELT2"]) #degree

		rawPixForX =  rawBeamDegree/resolX/np.cos(np.radians(38)) #
		# dut to the shrink of galactic longitude, the actually beam is a little bit larger
		# use a 38 degree
		# 
		rawPixForY =  rawBeamDegree/resolY
		
		
		 
		resultPixForX =  resultBeamDegree/resolX/np.cos(np.radians(38)) #
		resultPixForY = resultBeamDegree/resolY
		
		
		gaussianBeamX=(resultPixForX**2-rawPixForX**2)**0.5
		gaussianBeamY=(resultPixForY**2-rawPixForY**2)**0.5
		
		# the gaussianBeam is FWHM, the gaussian kernel should use sgima
		sigmaX=gaussianBeamX/2.355
		sigmaY=gaussianBeamY/2.355
		

		#creat an eliptical Gauss
		
		#print rawPixForX,rawPixForY
		#print resultPixForX,resultPixForY
		#print sigmaX,sigmaY



		g2 = Gaussian2D(10,x_stddev=sigmaX,y_stddev=sigmaY)
		smoothKernel=Model2DKernel(g2, x_size=17, y_size=17)
		
		smoothKernel.normalize()
		
		#print smoothKernel.array
		
		#smooth the data using this kernel
 
 		if len(data.shape)==3:
 		
	 		z,y,x=data.shape
 
	 
			for i in range(z):
				
				data[i]=convolve(data[i], smoothKernel)
				
		if len(data.shape)==2:
			data=convolve(data, smoothKernel)
		#written the fits
		
		#CO12[0].header["CTYPE3"]="VELO-LSR"
		#CO12.writeto("CO12SM20.fits")
		#outPutName
		#hdu = pyfits.PrimaryHDU(data)
		os.remove(outPutName)
		fits.writeto(outPutName, data, dataHeader)

	def getAverageSpec(self,fitsFile,path="./",cores=None):
		"""
		This function is dedicated to get the average spectrum for the CO lines
		
		#Basicly, we only care about the spectra at the peak posion of the cores
		
		#The index of peak position given by Duchamp is from 0.
		
		"""
		#read the file
		cores=self.getCoreTable()
		
		COdata,COheader=self.readFITS(path+fitsFile)
		
		#print len(cores)
		avgSpec=0
		
		for eachCore in cores:
		    #l,b= eachCore["GLON_deg"],eachCore["GLAT_deg"]
		    #spectrum,vs =self.getSpectraByLB( COdata,COheader,l,b) 
		
		    X,Y= eachCore["X_peak"],eachCore["Y_peak"]
		    spectrum,vs =self.getSpectraByIndex( COdata,COheader,int(X),int(Y))
		
		
		    avgSpec=avgSpec+spectrum
		    #print l,b,spectrum[0]
		avgSpec=avgSpec/1./len(cores)
		
		if 0:
		    l,b= cores[0]["GLON_deg"],cores[0]["GLAT_deg"]
		
		    avgSpec,vs =self.getSpectraByLB( COdata,COheader,l,b)
		
		if 0:
		    fig, ax = plt.subplots()
		    ax.plot(vs,avgSpec)
		    plt.show()
		
		return avgSpec,vs





	def mytrim(self,d, vmin, vmax):
		dd = (d + vmin) /(vmin + vmax)
		return np.clip(dd, 0, 1)






	def momentFITS(self,FITSFile,Vrange,mom,outFITS=None,cutEdge=False):##Vrange kms
		"""
		Parameters: FITSFile,Vrange,mom,outPutPath=None,outPutName=None,cutEdge=False
		
		This method do the moment operature with Miriad
		
		return the data and header of the moment FITS
		
		the outputPATH is nessary, because miriad needs a path to run

		The unit of Vrange is kms, and cutEdge is not concerned
		
		"""


		#If no outPutName or outPUTPATH is provided then no file is going to be saved

		#Split FITSFile

		
		processPath,FITSname=os.path.split(FITSFile);
 
		if processPath=="":
			processPath="./"

		else:
			processPath=processPath+"/"

		print "Doing moment {} in the velocity range of {} kms".format(mom,Vrange)

		ReadOutName="tempoutFITS"
 
		# if no file name is provided, created one
		#if not outPutName:
		#	outPUTName=FITSFile[0:-5]+"_M{}.fits".format(mom)
		

		if not outFITS:
			outPutName=FITSname[0:-5]+"_M{}.fits".format(mom)
		else:

			aa,outPutName=os.path.split(outFITS);


		#delete this frisrt
		deleteFITS1="rm -rf %s"%ReadOutName
		#step 1 read the file 
		ReadFITS="fits in=%s out=%s op=xyin"%(FITSname,ReadOutName)
	 
		#step 2
		##do moment
		#moment in=GCcont out=GCcont.2d mom=-1 region='kms,images(-50,50)'
		momentTempOut="momenttempout"
		deleteFITS2="rm -rf %s"%momentTempOut
		
		
		momentString="moment in=%s out=%s mom=%s region='kms,images(%s,%s)'"%(ReadOutName,momentTempOut,mom,Vrange[0],Vrange[1])
		
		##step 3 output the fits file
		
		deleteFITS3="rm -rf %s"%outPutName
		
		outPUTFITS="fits in=%s out=%s op=xyout"%(momentTempOut, outPutName )
		##step3 run commonds
  
		#goToPath="cd "+outPutPath
  
		goToPath="cd "+processPath

		saveScriptPath="scriptPath=$PWD"
		backToScriptPath="cd $scriptPath"
		copyFITS="mv {}{}  {}".format(processPath,outPutName,outPutName)
		if outFITS:
			pass
		
 
	
		self.runShellCommonds([saveScriptPath,goToPath,ReadFITS,momentString,deleteFITS3,outPUTFITS,deleteFITS1,deleteFITS2,backToScriptPath,copyFITS],"./")
 
		return self.readFITS(outPutName)



	@staticmethod
	def getSpectraByLB(data,dataHeader,l,b):
		"""
		Parameters: data, dataHeader,l,b
		
		This function is used to get a voxel value from a data cube
		
		
		
		v, km/s
		the unit of returned veloicyt is kms
		"""
		wcs = WCS(dataHeader)
		xindex,yindex=wcs.all_world2pix(l,b,0,0)[0:2]
		xindex=int(round(xindex))
		yindex=int(round(yindex))
 
		
		##it is possible that the yindex,and xindex can exceed the boundary of spectrum
 
		if yindex> data.shape[1]-1 or xindex>data.shape[2]-1:
			return None, None
		spectral=data[:,yindex,xindex]
		##just for test
		#print  w.all_world2pix(0, 0, 0,0)
		#print data.shape[0]
		velocityIndex= range(data.shape[0])
		
		velocities=wcs.all_pix2world(0, 0,velocityIndex,0)[2]/1000.
 
		# 
		return spectral,velocities

 

	def getSpectraByIndex(data,dataHeader,indexX,indexY):
		"""
		paramters: data,dataHeader,indexX,indexY
		This function is used to get a voxel value from a data cube
		
		v, km/s
		"""
		wcs = WCS(dataHeader)
  
		##it is possible that the yindex,and xindex can exceed the boundary of spectrum
 
		spectral=data[:, indexY,indexX]
		##just for test
		#print  w.all_world2pix(0, 0, 0,0)
		#print data.shape[0]
		velocityIndex= range(data.shape[0])
		
		velocities=wcs.all_pix2world(indexX, indexY,velocityIndex,0)[2]/1000.
 
		# 
		return spectral,velocities




	
	@staticmethod
	def downLoadSkyview(survey,centerLB,sizeLB=[0.1,0.1],resolution=None,overWrite=True,savePath="/home/qzyan/WORK/backup/tempskyviewDownload/",saveName=None):
		"""
		all sizes are in degree
		Download IRAS MBM,
		
		@resoltion of wise 22 0.00038189999999999996
		
		"""
		
		print "Download {}, at (l, b)=({}, {}) ----------  sizeL: {} deg; sizeB: {} deg".format(survey,centerLB[0],centerLB[1],sizeLB[0],sizeLB[1])
		
		
		
		
		centerL,centerB=centerLB
		
		sizeL,sizeB=centerLB
		
		if resolution is None:
			sizePix=500
		else:
			sizePix=max(sizeLB)/resolution
			sizePix=int(round(sizePix))
			
 
		
		if saveName is None:
			saveName="{}_{}_{}.fits".format(survey.replace(" ",""),centerLB[0],centerLB[1])
			
			
		
		print "Saving to  {}{}  ....".format(savePath,saveName)
		
		
		
		if not overWrite and os.path.isfile(savePath+saveName):
			print "File exist, returning...."
			return 
		
		#command="skvbatch_wget file=./{}/{} position='{}, {}' Survey='{}'  Projection='Car' Coordinates='Galactic' Pixels={}".format(savePath,tempName,centerl,centerb,survey,sizePix)
		command="skvbatch_wget file='{}{}' position='{}, {}' Survey='{}'  Projection='Car' Coordinates='Galactic' Pixels={}".format(savePath,saveName, centerL,centerB,survey,sizePix)



		os.system(command)
		#first test a small 
		
		
		
		
		
		return 
		
		
		sizePix=1500
		
		#survey="IRIS 100"
		 
		#savePath=self.IRASDownloadFolder
		savePath="MBMFITS" #IRAS 100 umself, #self.MBMFITSPath

		#saveNameRaw=MBMID+"_IRAS100.fits"
		
		if survey.endswith('545'):
		
			saveNameRaw=MBMID+"_planck545.fits"
		if survey.endswith('857'):
		
			saveNameRaw=MBMID+"_planck857.fits"
		if survey.endswith('100'):
		
			saveNameRaw=MBMID+"_IRAS100.fits"
		#survey="Planck 857"
		 
		
		saveName=saveNameRaw.replace(' ', '\ ')
		if not overWrite:
			
			outputFile="./{}/{}".format(savePath,saveName)
			outputFile2="./{}/{}".format(savePath,saveNameRaw)
			
	  
			
			if  os.path.isfile(outputFile2):
				print outputFile,"exits...doing nothing..."
				return 
		tempName="tempIRAS.fits" #in case some filenames contains space
		command="skvbatch_wget file=./{}/{} position='{}, {}' Survey='{}'  Projection='Car' Coordinates='Galactic' Pixels={}".format(savePath,tempName,centerl,centerb,survey,sizePix)
		
		os.system(command)
		copyFile="mv ./{}/{} ./{}/{} ".format(savePath,tempName,savePath,saveName)
		print copyFile
		os.system(copyFile )





	@staticmethod
	def readFITS(fitsFile):
		"""
		parameters: fitsFile
		This file will return the data and header of the fits
		"""
		
		fitsRead=fits.open(fitsFile)

		head=fitsRead[0].header
		try:
			del head["COMMENT"]
			del head["HISTROY"]
		except:
			pass


		return fitsRead[0].data,head


	def downLoadSurveyByRange(self,Survey,LRange,BRange,Original=False,Pixels=None,size=None):

		"""
		This survey is used to download a survey with the best resolution 
		covering the input LBRange

		default size is 0.5 degree
		"""
	

		processFolder=Survey.replace(" ","")+"_Mosaic"

		downLoadPath="./{}/".format(processFolder)
		os.system("mkdir "+downLoadPath)


		centerL,centerB=np.mean(LRange),np.mean(BRange)
		
		if not Pixels and not size:

			print "No download size is assigned, quit"
			return

		if Original:
			#if the largest pixel resolution is wanted

		
			resolution=self.detectSurveyResolution(Survey,[centerL,centerB])
			if Pixels and size:
				print "Can't assign size and pixels simultaneously in Original model, quit"
				return

			if size and not Pixels:
				Pixels=size/resolution
				

		if not Original:

			#in this case, pixels and size must be provided


			if not size:
				size=""
			if not Pixels:
				Pixels=""

		#download with pixels and size assigned




			#download with pixels size

		#a bug is here
		

		extraPixels=50

		downLoadSize=size/Pixels*extraPixels+size
		downLoadPixels=Pixels+extraPixels
		#estimate pixels
			
		tilesL=abs(LRange[1]-LRange[0])/size
		tilesB=abs(BRange[1]-BRange[0])/size

		minL=min(LRange)
		maxL=max(LRange)
		
		minB=min(BRange)
		maxB=max(BRange)

		mosaicFITS=[]

		print tilesL,tilesB

		for i in range(int(tilesL)+1):
			for j in range(int(tilesB)+1):

				centerTileL=minL+size*(i+1)
				centerTileB=minB+size*(j+1)
				print "Downloading (i,j)=({},{}) fits centered:{} {}".format(i+1,j+1,minL+size*(i+1),minB+size/1.*(j+1))
				
				#avoid repeat download
				outputName=downLoadPath+"Mosaic{}{}.fits".format(i+1,j+1)
				if os.path.isfile(outputName):
					continue


				self.getSurvey(Survey,[centerTileL,centerTileB],Pixels=downLoadPixels,size=downLoadSize,outputFITS=outputName)
				#print centerTileL,centerTileB
				#if Original:
					#self.getSurvey(Survey,[centerTileL,centerTileB],Pixels=int(sizePixels+50),outputFITS=outputName)
				#else:
					#download by size

					#self.getSurvey(Survey,[centerTileL,centerTileB],Pixels=downLoadPixels,size=downLoadSize,outputFITS=outputName)

					#self.getSurvey(Survey,[centerTileL,centerTileB],size=0.3,outputFITS=outputName)

		#after download, mergethose files with Montage

	def detectSurveyResolution(self,Survey,LB):
		"""

		to know the resolution of the survey
		"""
		#print Survey,LB

		outputFITS="checkRes.fits"
		self.getSurvey(Survey,LB,outputFITS=outputFITS,Pixels=100)

		data,header=self.readFITS(outputFITS)

		os.system("rm checkRes.fits")

		return abs(header["CDELT1"])# return degree
	@staticmethod
	def getSurvey(Survey,LB,outputFITS=None, Pixels=None,size=None):
		"""
		parameters: Survey,LB,outputFITS=None, Pixels=500,size=None
		
		band list is the wise band to be downloaded
		
		the unit of size is 0.25 degree
		0.25 will keep the original resolution
		size diameter

		"""
		#print "Function works!"
		
		centerL,centerB=LB


		if not outputFITS:
			outputFITS=Survey.replace(" ","")+"_"+str(centerL)+"_"+str(centerB)+".fits"

		if not Pixels and not size:

			#no download paratmers are providided, download with  pixels=500
			
			command="/usr/lib/python2.7/skvbatch_wget file={} position='{},{}' Survey='{}'  Coordinates=Galactic  Projection=Car  Pixels={}".format(outputFITS,centerL,centerB,Survey,500)

			os.system(command)
			return

		if Pixels and size:

			command="/usr/lib/python2.7/skvbatch_wget file={} position='{},{}' Survey='{}'  Coordinates=Galactic  Projection=Car size={} Pixels={}".format(outputFITS,centerL,centerB,Survey,size,Pixels)

			os.system(command)
			return
		if not size:

			command="/usr/lib/python2.7/skvbatch_wget file={} position='{},{}' Survey='{}'  Coordinates=Galactic  Projection=Car  Pixels={}".format(outputFITS,centerL,centerB,Survey,Pixels)
			#print command
		else:
			command="/usr/lib/python2.7/skvbatch_wget file={} position='{},{}' Survey='{}'  Coordinates=Galactic  Projection=Car size={} ".format(outputFITS,centerL,centerB,Survey,size)

		#print command
		os.system(command)

		#call("skvbatch file=example1.fits position='+12 34, -10 23'  Survey='Digitized Sky Survey'")
		#call("./skvbatch_wget file=example2.fits position='0,0' Survey='Digitized Sky Survey' Coordinats=Galactic  Projection=Car size=0.5")




	@staticmethod
	def cropFITS(inFITS,outFITS=None,Vrange=None,Lrange=None,Brange=None,overWrite=False):
		"""
		parameters: inFITS,outFITS=None,Vrange=None,Lrange=None,Brange=None,overWrite=False
		Thiss function is used to create my own function of croping FITS
		Based on Mongate

		In the first version, only Galactic coordinate is supported

		The output fits could be a little bit different as requested

		#no project is concerted in this function		

		# the unit of LBV is degree,degree, kms
		"""

		#read FITS file
			
		hdu=fits.open(inFITS)[0]		
		wmap=WCS(hdu.header)
		data=hdu.data
		if not Vrange and not Lrange and not Brange:
			print "No crop range is provided."
			return
		
		#Examine the maximum number for pixels

		zSize,ySize,xSize=data.shape


		Xrange=[0,xSize-1]  #Galactic Longitude  #
		Yrange=[0,ySize-1]  #Galactic Longitude  #
		Zrange=[0,zSize-1]  #Galactic Longitude  #
		
		firstPoint=wmap.wcs_pix2world(0,0,0,0)
		lastPoint=wmap.wcs_pix2world(xSize-1,ySize-1,zSize-1,0)

		if not Vrange:
			#calculate the range for the 
			Zrange=[firstPoint[2],lastPoint[2]]
		else:
			Zrange=np.array(Vrange)*1000

		if not Lrange:
			Xrange=[firstPoint[0],lastPoint[0]]
		else:
			Xrange=Lrange
		if not Brange:
			Yrange=[firstPoint[1],lastPoint[1]]
		else:
			Yrange=Brange
		
		#revert Galactic longtitude
		if lastPoint[0]<firstPoint[0]:
			Xrange=[max(Xrange),min(Xrange)]
		#print lastPoint[0],firstPoint[0]
	
		#print Xrange,Yrange,Zrange
		cutFIRST=wmap.wcs_world2pix(Xrange[0],Yrange[0],Zrange[0],0)
		cutLAST =wmap.wcs_world2pix(Xrange[1],Yrange[1],Zrange[1],0)

		cutFIRST=map(round,cutFIRST)
		cutLAST=map(round,cutLAST)

		cutFIRST=map(int,cutFIRST)
		cutLAST=map(int,cutLAST)

		cutFIRST[0]=max(0,cutFIRST[0])
		cutFIRST[1]=max(0,cutFIRST[1])
		cutFIRST[2]=max(0,cutFIRST[2])

		cutLAST[0]=min(xSize-1,cutLAST[0])
		cutLAST[1]=min(ySize-1,cutLAST[1])
		cutLAST[2]=min(zSize-1,cutLAST[2])
		#calculate the true pixels according to the input range
		wmapcut=wmap[cutFIRST[2]:cutLAST[2],cutFIRST[1]:cutLAST[1],cutFIRST[0]:cutLAST[0]]
		datacut=data[cutFIRST[2]:cutLAST[2],cutFIRST[1]:cutLAST[1],cutFIRST[0]:cutLAST[0]]
		#datacut=data[1:3,1:5,1:9]

		#hdu = fits.PrimaryHDU(datacut,header=wmapcut)


		if not outFITS:
			"""
			If no output file Name is provide
			"""
			outFITS=inFITS[:-5]+"_C.fits"	

		if not os.path.isfile(outFITS):
			 
			fits.writeto(outFITS,datacut,header=wmapcut.to_header())

		else:

			if overWrite:
				#delete that file
				os.remove(outFITS)
	
				#hdu.data=datacut
				fits.writeto(outFITS,datacut,header=wmapcut.to_header())

			else:
				print "Warring----File ({}) exists and no overwriting!".format(outFITS)


	@staticmethod
	def cropFITS2D(inFITS,outFITS=None,Vrange=None,Lrange=None,Brange=None,overWrite=False):
		"""
		parameters: inFITS,outFITS=None,Vrange=None,Lrange=None,Brange=None,overWrite=False
		Thiss function is used to create my own function of croping FITS
		Based on Mongate

		In the first version, only Galactic coordinate is supported

		The output fits could be a little bit different as requested

		#no project is concerted in this function		

		# the unit of LBV is degree,degree, kms
		"""

		#read FITS file
			
		hdu=fits.open(inFITS)[0]		
		wmap=WCS(hdu.header)
		data=hdu.data
		if not Vrange and not Lrange and not Brange:
			print "No crop range is provided."
			return
		
		#Examine the maximum number for pixels

		ySize,xSize =data.shape

		"Dose this work?"

		Xrange=[0,xSize-1]  #Galactic Longitude  #
		Yrange=[0,ySize-1]  #Galactic Longitude  #
		
		firstPoint=wmap.wcs_pix2world(0,0 ,0)
		lastPoint=wmap.wcs_pix2world(xSize-1,ySize-1 ,0)

 

		if not Lrange:
			Xrange=[firstPoint[0],lastPoint[0]]
		else:
			Xrange=Lrange
		if not Brange:
			Yrange=[firstPoint[1],lastPoint[1]]
		else:
			Yrange=Brange
		
		#revert Galactic longtitude
		#if lastPoint[0]<firstPoint[0]:
		Xrange=[max(Xrange),min(Xrange)]
		Yrange=[min(Yrange),max(Yrange)]
 
		#print lastPoint[0],firstPoint[0]
	
		#print Xrange,Yrange,Zrange
		cutFIRST=wmap.wcs_world2pix(Xrange[0],Yrange[0], 0)
		cutLAST =wmap.wcs_world2pix(Xrange[1],Yrange[1], 0)

 
		cutFIRST=map(round,cutFIRST)
		cutLAST=map(round,cutLAST)

		cutFIRST=map(int,cutFIRST)
		cutLAST=map(int,cutLAST)

		cutFIRST[0]=max(0,cutFIRST[0])
		cutFIRST[1]=max(0,cutFIRST[1])

		cutLAST[0]=min(xSize-1,cutLAST[0])
		cutLAST[1]=min(ySize-1,cutLAST[1])
		
		
 

		#calculate the true pixels according to the input range
		wmapcut=wmap[ cutFIRST[1]:cutLAST[1],cutFIRST[0]:cutLAST[0]]
		datacut=data[ cutFIRST[1]:cutLAST[1],cutFIRST[0]:cutLAST[0]]
		#datacut=data[1:3,1:5,1:9]

		#hdu = fits.PrimaryHDU(datacut,header=wmapcut)


		if not outFITS:
			"""
			If no output file Name is provide
			"""
			outFITS=inFITS[:-5]+"_C.fits"	

		if not os.path.isfile(outFITS):
			 
			fits.writeto(outFITS,datacut,header=wmapcut.to_header())

		else:

			if overWrite:
				#delete that file
				os.remove(outFITS)
	
				#hdu.data=datacut
				fits.writeto(outFITS,datacut,header=wmapcut.to_header())

			else:
				print "Warring----File ({}) exists and no overwriting!".format(outFITS)


 
	@staticmethod
	###############Static functions#######################
	def runShellCommonds(commondList,processPath):
		
		"""
		parameters: commondList,processPath
		run shell commonds
		"""
		##write each commondList into a file and run them and then delete this file
		
		tempSHfile=processPath+"runcommondlistTemp.sh"
		f = open(tempSHfile,'w')
		for eachLine in commondList:
			f.write(eachLine+"\n")
		f.close()
	
		os.system("bash  %s"%( tempSHfile)+"   >/dev/null 2>&1" ) #the end string supress the output of os.system
		# delete file
 
		os.system("rm -rf "+tempSHfile)
		#self.l=l_input
		#self.b=b_input
