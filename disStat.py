import numpy as np
#to draw statistics for a good list of molecular clouds

from astropy.table import Table,vstack,hstack
from  distanceTB import disTB
from astropy.coordinates import SkyCoord
import astropy.units as u
import os
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.offsetbox import AnchoredText
from matplotlib import rc
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from astropy.io import ascii
import sys

from myPYTHON import myFITS
from scipy import optimize
import glob
doFITS = myFITS()
doDis=disTB()
import subprocess

def larsonRelation(x,a,b):
	return a*x+b

def larsonRelationPower(x,a,b):
	return  a*x**b

class disDraw:
	R0=8.178 #34 the background would be changed #  2019A&A...625L..10G

	regionName=None

	goodTBName=None
	dbscanDBName=None
	dbscanTB=None
	parsecToMeter= 3.0857e16 #m

	edgeCloudNameList=None
	edgeCloudIDList=None

	positionGC="GC"
	positionSun="SUN"


	###### 3D position colnames

	colXsun="xSun"
	colYsun="ySun"
	colZsun="zSun"

	colXgc="xGC"
	colYgc="yGC"
	colZgc="zGC"


	## get error position

	colXsunError16="xSun16" #close to the sun
	colYsunError16="ySun16" #close to the sun
	colZsunError16="zSun16" #close to the sun

	colXsunError84="xSun84" #close to the sun
	colYsunError84="ySun84" #close to the sun
	colZsunError84="zSun84" #close to the sun


	besselTBFile="/home/qzyan/WORK/dataDisk/cataLog/apjab4a11t1_mrt.txt"

	besselTB = ascii.read(besselTBFile)

	besselTBFileProcessed="/home/qzyan/WORK/dataDisk/cataLog/apjab4a11t1_mrt.fit"
	besselTBProcessed = Table.read( besselTBFileProcessed )


	massCol="mass"
	massColError="massError"
	lwCol="lineWidth"
	dv=0.167 #km/s
	beamSize=49./50
	disCol= "distance"
	lCol= "l"
	bCol = "b"
	disStdCol = "disStd"
	def __init__(self, goodTBName=None,dbscanTBName=None,regionName=None,calMass=False):


		if regionName!=None:
			self.regionName= regionName
		else:
			print "better provided a region name "

		if goodTBName!=None and dbscanTBName!=None:
			self.goodTBName = goodTBName #Table.read(goodTBName)
			self.dbscanTBName=dbscanTBName

			self.goodTB= Table.read( self.goodTBName  )
			self.dbscanTB=Table.read(self.dbscanTBName)

			self.mergeGoodAndDBSCANTB()
			self.add3DInfo( self.goodTB ) #to goodTB
		else:
			print "please provide a cloud table with distance determined and the DBSCAN Table"



		#calculate mass
		if calMass:
			massArray,massError=self.getMassList()
			self.goodTB[self.massCol]=massArray
			self.goodTB[self.massColError] =  massError


		#self.besselConvert() #add distance ,l, b to the besselConvert



	def besselConvert(self):
		"""
		convert the coordinate to maser from ICRS to Galactic
		:return:
		"""

		lList=[]

		bList=[]

		distanceList=[]
		distanceErrorList=[]


		#remove low precision measurement

		paraErrorRelative =  self.besselTB["e_plx"]/self.besselTB["plx"]
		self.besselTB = self.besselTB[paraErrorRelative<=0.2]
		for eachRow in self.besselTB:
			raStr="{}h{}m{}s".format(eachRow["RAh"], eachRow["RAm"], eachRow["RAs"]  )
			decStr="{}{}d{}m{}s".format( eachRow["DE-"] , eachRow["DEd"] , eachRow["DEm"]  , eachRow["DEs"]    )
			c=SkyCoord( raStr  , decStr , frame='icrs')

			lList.append(    c.galactic.l.deg )
			bList.append( c.galactic.b.deg )

			paraSample=np.random.normal(eachRow["plx"], eachRow["e_plx"],20000)
			disSample=1000./paraSample #pc

			distanceList.append(np.mean(disSample) )

			distanceErrorList.append( np.std(disSample,ddof=1)  )


		#add four colomns to the table

		self.besselTB[disTB.l] = np.asarray( lList)
		self.besselTB[disTB.b] = np.asarray( bList)

		self.besselTB[disTB.distance] = np.asarray(  distanceList )

		self.besselTB[disTB.disStd] = np.asarray (distanceErrorList)
		self.besselTB[disTB.cloudVlsr] = self.besselTB["VLSR"] # in km/s
		self.besselTB["v_cen"] = self.besselTB["VLSR"] # in km/s

		print self.besselTB.write(self.besselTBFileProcessed,overwrite=True)

	def mergeGoodAndDBSCANTB(self):
		"""
		merge the table into to, add DBSCAN info the GoodTable
		:return:
		"""

		goodDBSCAN=self.getGoodTBDBSCAN()
		self.goodTB=hstack([self.goodTB,goodDBSCAN])

		self.goodTB.write(self.regionName+"GoodDBSCANMerge.fit",overwrite=True)





	def compareDisWithA5(self,withReid2019=False ):

		"""
		compare distances with A5 model of
		"""

		#gaia distancs ,may come from different tables;

		gaiaDis=[]
		gaiaDisLow=[]
		gaiaDisUp=[]


		A5Dis=[]
		A5DisLow=[]
		A5DisUp=[]



		disCatTB=self.goodTB



		lS=[]
		vS=[]


		diffGaiaList=[]

		diffPerseus=[]

		for eachC in disCatTB:

			l=eachC[disTB.l]  #cloudVlsr
			b=eachC[disTB.b]
			#v=eachC[disTB.cloudVlsr]*1000


			v = eachC[disTB.cloudVlsr]

			if withReid2019:
				disReid =self.Reid2019(l,b,v)

			else:
				#disReid,rGC,alpha =self.ReidA5(l,b,v)
				disReid  =self.Reid2016(l,b,v)



			A5Dis.append(disReid[0]*1000 )

			A5DisLow.append(disReid[1]*1000 )
			A5DisUp.append(disReid[2]*1000 )

			gaiaDis.append (   eachC[disTB.distance] )
			#gaiaDisLow.append (  eachC[disTB.disHPDLow]  )
			#gaiaDisUp.append ( eachC[disTB.disHPDUp]  )

			# use systematic error here



			disStd =   eachC["disStd"]

			errorWithSys= disStd**2+ ( 0.05* eachC[disTB.distance] )**2
			sysError=np.sqrt(errorWithSys)
			#disErrorWithSys.append( sysError )

			gaiaDisLow.append ( sysError )
			gaiaDisUp.append ( sysError )

			lS.append(l)
			vS.append(v)
			#if v>=10:
			diffGaiaList.append( A5Dis[-1] - gaiaDis[-1] )
			if eachC["vlsr"] <-29.9  :
				diffPerseus.append( A5Dis[-1] - gaiaDis[-1])

		print "Systematic shift",np.nanmean( diffGaiaList )
		print "Systematic shift perseus",np.nanmean( diffPerseus ),len(diffPerseus)

		lS,gaiaDis,gaiaDisLow,gaiaDisUp,A5Dis,A5DisLow,A5DisUp=   zip(*sorted(zip(lS,gaiaDis,gaiaDisLow,gaiaDisUp,A5Dis,A5DisLow,A5DisUp), reverse=True))

		lS,gaiaDis,gaiaDisLow,gaiaDisUp,A5Dis,A5DisLow,A5DisUp=map(np.array, (lS,gaiaDis,gaiaDisLow,gaiaDisUp,A5Dis,A5DisLow,A5DisUp) )

		fig, axs = plt.subplots(nrows=1, ncols=2,  figsize=(8,6),sharex=True)
		#ax1=axs.flat[0]
		rc('text', usetex=True )
		#rc('text.latex',  preamble=r'\usepackage{upgreek}')
		rc('font', **{'family': 'sans-serif',  'size'   :  12,  'serif': ['Helvetica'] })
		#diff=np.array(A5Dis)- np.array(gaiaDis)


		mpl.rcParams['text.latex.preamble'] = [
			r'\usepackage{tgheros}',  # helvetica font
			r'\usepackage{sansmath}',  # math-font matching  helvetica
			r'\sansmath'  # actually tell tex to use it!
			r'\usepackage{siunitx}',  # micro symbols
			r'\sisetup{detect-all}',  # force siunitx to use the fonts
		]

		systematic=np.array(  gaiaDis)*0.05
		systematic=systematic**2+(  (gaiaDisUp+gaiaDisLow)/4.   )**2

		systematic=np.sqrt(systematic)

		gaiaDisLow= systematic #np.array( gaiaDisLow) + systematic
		gaiaDisUp= systematic #np.array( gaiaDisUp) + systematic

		ax=axs[0]
		#ax.tick_params( labelleft=True, labelright=True,    left=True, right=True)
		#ax.set_ylabel(r"Maser\-parallax\-based distance (kpc)")
		ax.set_ylabel( "  Maser-parallax-based distance (kpc)")


		ratio=1./1000.  #to kpc


		ax2=axs[1]

		#print dir(ax2)

		ax.errorbar( gaiaDis*ratio,A5Dis*ratio ,yerr= [abs(A5DisLow)*ratio,  A5DisUp*ratio] , fmt='o',  xerr=[gaiaDisLow*ratio,  gaiaDisUp*ratio],  c='b',marker='o',capsize=1.0,elinewidth=0.5,lw=0.5,markersize=1.5 )

		#ax.scatter( gaiaDis*ratio,A5Dis*ratio,s=3,c='b')


		ax2.errorbar( gaiaDis*ratio,A5Dis*ratio-gaiaDis*ratio ,yerr= [abs(A5DisLow)*ratio,  A5DisUp*ratio] , fmt='o',  xerr=[gaiaDisLow*ratio,  gaiaDisUp*ratio],  c='r',marker='o',capsize=1.0,elinewidth=0.5,lw=0.5,markersize=1.5  )
		#ax2.scatter(gaiaDis*ratio,A5Dis*ratio-gaiaDis*ratio,s=3,c='r')

		#ax.errorbar(gaiaDis,diff,yerr=[np.abs(A5DisLow), A5DisUp],xerr=[np.abs(gaiaDisLow), gaiaDisUp] ,c='b',marker='o',capsize=1.5,elinewidth=0.8,lw=0.4,markersize=3.5 )
		maxX=2

		ax.plot( [0,4], [0,4],'--',color='black',lw=1.0,label="1:1 line" )
		ax2.plot( [0,4], [0,0],'--',color='black',lw= 1.0 )

		#ax.set_xticks([0,0.5,1.0,1.5,2,2.5,3 , 3.5] )
		#ax.set_yticks([0,0.5,1.0,1.5,2,2.5,3, 3.5 ])
		ax.set_xticks([0,1,2,3,4] )
		ax.set_yticks([0,1,2,3,4])

		ax.set_xlim(-0.5,4 )
		ax.set_ylim( -1,4.5 )
		#ax2.set_ylim( -1,1.3)

		#ax2.get_shared_y_axes().join(ax2, ax)

		ax2.set_xlabel(r"Gaia distance (kpc)")
		ax.set_xlabel(r"Gaia distance (kpc)")

		ax.set_aspect('equal' )
		ax2.set_aspect('equal' )
		ax2.set_xticks([0,1,2,3,4 ])
		ax2.set_yticks([  -1,-0,1,2])

		ax2.set_ylabel(r"Distance difference (kpc)")
		plt.subplots_adjust(   hspace=-0.08)

		ax2.set_ylim( -1.5, 3 )
		ax.set_ylim( -0.5,4)


		plt.subplots_adjust(wspace=0.25)
		ax.legend( loc=4, handlelength=1 )
		#plt.gca().invert_xaxis()
		plt.tight_layout()


		#mar a label
		if withReid2019:

			at = AnchoredText('Reid et al. (2019)', loc=2, frameon=False)
			ax.add_artist(at)


		else:

			at = AnchoredText('Reid et al. (2016)', loc=4, frameon=False)
			ax.add_artist(at)

		########


		if self.regionName==None:

			if withReid2019:

				plt.savefig( "compareReid2019G2650.png",bbox_inches='tight',dpi=600)

				plt.savefig( "compareReid2019G2650.pdf", bbox_inches='tight')
			else:
				plt.savefig( "compareA5G2650.png",bbox_inches='tight',dpi=600)

				plt.savefig( "compareA5G2650.pdf", bbox_inches='tight')

		else   :

			if withReid2019:
				plt.savefig( "compareReid2019{}.png".format(self.regionName),bbox_inches='tight',dpi=600)

				plt.savefig( "compareReid2019{}.pdf".format(self.regionName), bbox_inches='tight')


			else:

				plt.savefig( "compareA5{}.png".format(self.regionName),bbox_inches='tight',dpi=600)

				plt.savefig( "compareA5{}.pdf".format(self.regionName), bbox_inches='tight')





	def getCloudNameByLB(self, l,b):

		#if b>=0:

			lStr= str(l)

			bStr="{:+f}".format(b)


			if '.' in lStr:

				lStrPart1,lStrPart2=lStr.split('.')

			else:
				lStrPart1 =lStr
				lStrPart2='0'


			if '.' in bStr:

				bStrPart1,bStrPart2=bStr.split('.')
			else:
				bStrPart1 =bStr
				bStrPart2='0'


			lStr=lStrPart1+'.'+lStrPart2[0:1]


			bStr=bStrPart1+'.'+bStrPart2[0:1]


			lStr=lStr.zfill(5)


			#bStr="{:+.1f}".format(b)
			bStrNumberPart=bStr[1:]
			bStr=bStr[0:1]+  bStrNumberPart.zfill(4)

			cName="G{}{}".format(lStr,bStr)

			return cName


	def getVDlist(self,l=25,b= 0,minV=2,maxV=30,useReid2019=False ):
		"""

		:param l:
		:return:
		"""


		vList= np.arange(minV,maxV)
		dList= []

		for eachV in vList:
			#disAndError,a,a=self.ReidA5( l,b, eachV )
			if useReid2019:
				dis,errorLow,errorUp = self.Reid2019( l,b, eachV )
			else:
				#disAndError,a,a=self.Reid2016( l,b, eachV )
				dis,errorLow,error =self.Reid2016( l,b, eachV )

				#disAndError,a,a=self.ReidA5( l,b, eachV )
				#dis= disAndError[0]
			dList.append( dis )


		return vList,dList

	def ReidA5(self, l,b,v):

		"""
		return the distance with A5 model of Reid 2014

		#what about the errors?

		return distance
		"""


		#first convert l,b to  the string use by the ma

		c = SkyCoord(l=l*u.degree, b=b*u.degree, frame='galactic')


		cor=str(c.icrs.to_string( 'hmsdms'))

		raStr,decStr=cor.split()

		raStr=raStr.replace("h",'')
		raStr=raStr.replace("m",'')
		raStr=raStr.replace("s",'')

		decStr=decStr.replace("d",'')
		decStr=decStr.replace("m",'')
		decStr=decStr.replace("s",'')

		sourceStr="calsource {} {} {} 0".format(raStr,decStr,v)

		#our source

		outSourceStr= 'echo "{}" > /home/qzyan/WORK/projects/maddalena/2014revised_kd/source_file.dat'.format(sourceStr)


		os.system(outSourceStr)

		os.system('currentP="`pwd`";cd /home/qzyan/WORK/projects/maddalena/2014revised_kd; ./a.out > distanceInfo.txt;   cd $currentP')

		with open("/home/qzyan/WORK/projects/maddalena/2014revised_kd/distanceInfo.txt") as f:
			lines=f.readlines()

		disStr=""
		for eachL in lines:

			if eachL[0]=="!":
				continue

			disStr=eachL

		distance= float(disStr.split()[5])

		disErrorUp=  float(disStr.split()[6])
		disErrorLow=  float(disStr.split()[7])

		rGC,alpha=self.getRgcAlpha(distance,l,b)


		return [distance,disErrorLow,disErrorUp  ], rGC,alpha
		#print this to  ./2014revised_kd/source_file.data


	def drawVelLongitudeDiagram(self,useReid2019=False,yRange=[-0.1,5] ):
		"""
		draw a diagram for the distance,v,diagram figures
		:return:
		"""

		import matplotlib as mpl

		fig=plt.figure(figsize=(10,8))
		#fig, axs = plt.subplots(nrows=1, ncols=2,  figsize=(12,6),sharex=True)
		rc('text', usetex=True )
		rc('font', **{'family': 'sans-serif',  'size'   : 20,  'serif': ['Helvetica'] })

		mpl.rcParams['text.latex.preamble'] = [
			r'\usepackage{tgheros}',  # helvetica font
			r'\usepackage{sansmath}',  # math-font matching  helvetica
			r'\sansmath'  # actually tell tex to use it!
			r'\usepackage{siunitx}',  # micro symbols
			r'\sisetup{detect-all}',  # force siunitx to use the fonts
		]
		ax=fig.add_subplot(1,1,1)

		#positiontion 1
		vList25, dList25 = self.getVDlist(l=105,b=-5,minV=-45,maxV=10 ,useReid2019=useReid2019 )
		ax.plot(vList25, dList25,   marker='o', markersize= 4 ,lw=0.8, label=r"$l={}^\circ$, $b={}^\circ$".format(105,-5) )

		#positiontion 2
		vList25, dList25 = self.getVDlist(l=105,b= 5,minV=-45,maxV=10 ,useReid2019=useReid2019 )
		ax.plot(vList25, dList25,   marker='o', markersize=4 ,lw=0.8, label=r"$l={}^\circ$, $b={}^\circ$".format(105, 5) )


		#positiontion 3
		vList25, dList25 = self.getVDlist(l=150,b= -5,minV=-45,maxV=10  ,useReid2019=useReid2019 )
		ax.plot(vList25, dList25,  marker='o', markersize=4 ,lw=0.8, label=r"$l={}^\circ$, $b={}^\circ$".format(150, -5) )

		#positiontion 4
		vList25, dList25 = self.getVDlist(l=150,b= 5,minV=-45,maxV=10 ,useReid2019=useReid2019 )
		ax.plot(vList25, dList25,  marker='o', markersize=4 ,lw=0.8, label=r"$l={}^\circ$, $b={}^\circ$".format(150, 5) )


		#positiontion 5
		vList25, dList25 = self.getVDlist(l=125,b= 0,minV=-45,maxV=10 ,useReid2019=useReid2019 )
		ax.plot(vList25, dList25,  marker='o', markersize=4 ,lw=0.8, label=r"$l={}^\circ$, $b={}^\circ$".format(125, 0) )
		ax.set_xlabel(r"$V_{\rm LSR}$ (km s$^{-1}$)")
		ax.set_ylabel(r"Kinematic distance (kpc)")

		ax.legend()
		ax.set_ylim(yRange)

		if useReid2019:

			at = AnchoredText('Reid et al. (2019)', loc=3, frameon=False)
			ax.add_artist(at)
			plt.savefig("velDisRelation2019.png" , bbox_inches='tight', dpi=600)

		else:
			at = AnchoredText('Reid et al. (2016)', loc=3, frameon=False)
			ax.add_artist(at)
			plt.savefig("velDisRelation2016.png" , bbox_inches='tight', dpi=600)

	def Reid2016(self,l,b,v):

		"""

		:param l:
		:param b:
		:param v:
		:return:
		"""

		os.system("cp ./reid2016KD/EmptySources_info.inp  ./reid2016KD/sources_info.inp")

		if b<0:
			sourceInfo="G{:.2f}-{:.2f}      {:.2f}   {:.2f}   {:.1f}  0.0    0.0    3kN-source@4.8kpc".format(l,abs(b),l,b,v)
		else:
			sourceInfo="G{:.2f}+{:.2f}      {:.2f}   {:.2f}   {:.1f}  0.0    0.0    3kN-source@4.8kpc".format(l,abs(b),l,b,v)


		file_object=open( "./reid2016KD/sources_info.inp" ,'a')
		file_object.write(sourceInfo)
		file_object.close()

		#run the script

		#os.system(  ""  )
		os.system("cd ./reid2016KD; ./a.out  >message.txt; cd ..")

		return self.getDistanceAndError2016()
	def Reid2019(self,l,b,v):

		"""

		:param l:
		:param b:
		:param v:
		:return:
		"""

		os.system("cp /home/qzyan/WORK/projects/maddalena/reid2019KD/EmptySources_info.inp  /home/qzyan/WORK/projects/maddalena/reid2019KD/sources_info.inp")

		if b<0:
			sourceName = "G{:.2f}-{:.2f}".format(l,abs(b))
			sourceInfo="G{:.2f}-{:.2f}      {:.2f}   {:.2f}   {:.1f}  0.0    0.0    3kN-source@4.8kpc".format(l,abs(b),l,b,v)
		else:
			sourceInfo="G{:.2f}+{:.2f}      {:.2f}   {:.2f}   {:.1f}  0.0    0.0    3kN-source@4.8kpc".format(l,abs(b),l,b,v)
			sourceName = "G{:.2f}+{:.2f}".format(l,abs(b))

		file_object=open( "/home/qzyan/WORK/projects/maddalena/reid2019KD/sources_info.inp" ,'a')
		file_object.write(sourceInfo)
		file_object.close()

		#run the script

		#os.system(  ""  )
		#aa=os.system("cd ./reid2019KD; ./a.out  ; cd ..")
		#bbb= os.system('cd ./reid2019KD; ./a.out | grep -m 1 -o "X3333X"  ; rm {}* ;cd ..'.format(sourceName))
		#os.system("cd ./reid2019KD; ./a.out >message.txt  ; cd ..")
		#x = subprocess.run("cd ./reid2019KD; ./a.out | grep -m 1 -o "X3333X ; cd ..", capture_output=True)
		#x = subprocess.run("cd ./reid2019KD; ./a.out /dev/null ; cd ..", capture_output=True)

		#bbb= os.system('./reid2019KD/a.out | grep -m 1 -o "X3333X"   '.format(sourceName))
		os.system('cd /home/qzyan/WORK/projects/maddalena/reid2019KD; ./a.out | grep -m 1 -o "X3333X" ; cd /media/qzyan/maclinux/projects/disAntiGC '.format(sourceName))


		return self.getDistanceAndErrorFromSummary()

		#outSourceStr= 'echo  {}  >> ./reid2019KD/sources_info.inp'.format(sourceInfo)
		#os.system(outSourceStr)
		#print outSourceStr

	def getDistanceAndErrorFromSummary(self):

		with open("/home/qzyan/WORK/projects/maddalena/reid2019KD/summary.prt") as f:
			lines=f.readlines()

		lastLine= lines[-1]
		if lastLine[0]=="!":
			#return [None,None,None]
			return [0,0,0]

		disStr= lastLine

		distance= float(disStr.split()[4])

		disErrorUp=  float(disStr.split()[5])
		disErrorLow=  float(disStr.split()[5])

		return [distance,disErrorLow,disErrorUp  ]

	def getDistanceAndError2016(self):
		"""
		There is no
		:return:
		"""
		with open("./reid2016KD/message.txt") as f:
			lines=f.readlines()


		disStr =  lines[-1]


		distance= float(disStr.split()[3])

		disErrorUp=  float(disStr.split()[4])
		disErrorLow=  float(disStr.split()[4])

		return [distance,disErrorLow,disErrorUp  ]

	def getXYZ(self,d,l,b):
		"""
		:param d:
		:param l:
		:param b:
		:return:
		"""
		radB=np.deg2rad(b)
		radl=np.deg2rad(l)

		z=d*np.sin(radB)
		xyProj = d*np.cos(radB)

		x=self.R0 - xyProj*np.cos(radl)

		y=  xyProj*np.sin(radl)

		return x,y,z


	def get3DPosition(self, d,l,b, withRespectTo="SUN"):
		"""
		3D position with respect to the Galactic center or the sun


		d, in kpc

		l, b in degreee

		:return: return 3D poistions with respect to the GC or the Sun
		"""

		if max(d)>100:
			d=d/1000. #convert to kpc

		bRad= np.deg2rad(b)


		dProj=d*np.cos(  bRad)

		angle0=np.deg2rad(l-90)

		#with respect to the sun
		xSUN=dProj*np.cos(angle0)
		ySUN=dProj*np.sin(angle0)
		zSUN=d*np.sin(bRad)
		#with respect to the Galactic center
		xGC= xSUN
		yGC= ySUN+self.R0
		zGC=zSUN


		if withRespectTo==self.positionGC:
			return  xGC,yGC,zGC

		if withRespectTo==self.positionSun:
			return   xSUN,ySUN,zSUN

		return None #

	def getRgcAlpha(self,d,l,b):

		"""

		"""
		#d=d*np.cos(np.radians(b) )

		R0=self.R0 #8.2 #34 the background would be changed

		Theta0=240.
		dThetadR=-0.2

		#Theta=theta0+dThetadR*(R-R0) #km/s

		d=d*np.cos( np.deg2rad(b) )

		angle0=np.deg2rad(l-90)


		#with respect to the sun
		v1x=d*np.cos(angle0)
		v1y=d*np.sin(angle0)
		v2z=d*np.sin(np.deg2rad(b))
		#with respect to the Galactic center
		vx=v1x
		vy=v1y+R0

		r=np.sqrt( vx**2+vy**2)


		if vy>0:

			alpha=np.arccos(  vx/r)
		if vy<0 and vx>0:
			alpha=np.arcsin(  vy/r)

		if vy<0 and vx<0:
			alpha=np.arcsin(  abs(vy/r) )+np.pi




		return r,alpha*180/np.pi #inradians


	def getDLB_ByVelRelation(self, TB):
		"""

		:param TB:
		:return:
		"""

		vList = TB["v_cen"]

		disList = -0.026* vList +0.564


		return  disList  , TB["x_cen"],TB["y_cen"]
	def getDLB_Reid2019(self,TB): #get distance,l,b,from masers, need another function
		"""

		:param TB:
		:return:
		"""

		#calculate the kinamtic distances with distance error
		disList= []

		disErrorLowList = []


		disErrorUpList = []

		for eachRow in TB:
			l = eachRow["x_cen"]
			b = eachRow["y_cen"]
			v = eachRow["v_cen"]
			disReid  = self.Reid2019(l,b,v)


			disList.append( disReid[0] )

			disErrorLowList.append( disReid[1]    )
			disErrorUpList.append( disReid[2]    )


			#print disReid[0], v

		return  np.asarray(   disList )  , np.asarray(disErrorLowList), np.asarray(disErrorUpList),   TB["x_cen"],TB["y_cen"]


	def getDLB_Reid2014(self,TB): #get distance,l,b,from masers, need another function
		"""

		:param TB:
		:return:
		"""

		#calculate the kinamtic distances with distance error
		disList= []

		disErrorLowList = []


		disErrorUpList = []

		for eachRow in TB:
			l = eachRow["x_cen"]
			b = eachRow["y_cen"]
			v = eachRow["v_cen"]
			disReid,rGC,alpha = self.ReidA5(l,b,v)


			disList.append( disReid[0] )

			disErrorLowList.append( disReid[1]    )
			disErrorUpList.append( disReid[2]    )


			#print disReid[0], v

		return  np.asarray(   disList )  , np.asarray(disErrorLowList), np.asarray(disErrorUpList),   TB["x_cen"],TB["y_cen"]


	def getDLB(self,TB ,disCol=None,lCol=None,bCol=None): #get distance,l,b,from masers, need another function

		if disCol is not  None:
			self.disCol= disCol


		if lCol is not None:
			self.lCol = lCol

		if bCol is not  None:
			self.bCol  = bCol


		print self.disCol,self.lCol, self.bCol ,"??????????"

		return TB[  self.disCol ]/1000,TB[self.lCol],TB[self.bCol ]

	def getDisErrorWithSys(self, withSysError=True):

		if withSysError:

			d=self.goodTB["distance"]/1000.
			disStd = self.goodTB["disStd"] / 1000.
			errorWithSys = disStd ** 2 + (0.05 * d) ** 2
			errorWithSys = np.sqrt(errorWithSys)
			return errorWithSys
		else:
			disStd = self.goodTB["disStd"] / 1000.
			return disStd


	def add3DInfoMWISP(self,processTB):
		"""
		Add the distances by the velcotiy distances
		:param processTB:
		:return:
		"""


		d,    l,b=self.getDLB_ByVelRelation(processTB )

		errorWithSys=  198/1000.  #  ( dErrorUp - dErrorLow )/2.   #dist std, Reid 2014, no distance error of systematic will be considered

		d16= d-errorWithSys
		d84=d+errorWithSys


		xSun,ySun,zSun=self.get3DPosition(d,l,b,withRespectTo=self.positionSun)
		xSun16,ySun16,zSun16=self.get3DPosition(d16,l,b,withRespectTo=self.positionSun)
		xSun84,ySun84,zSun84 =self.get3DPosition(d84,l,b,withRespectTo=self.positionSun)




		processTB[self.colXsun]= xSun
		processTB[self.colYsun]=ySun
		processTB[self.colZsun]=zSun

		processTB[self.colXsunError16]= xSun16
		processTB[self.colYsunError16]=ySun16
		processTB[self.colZsunError16]=zSun16

		processTB[self.colXsunError84]= xSun84
		processTB[self.colYsunError84]=ySun84
		processTB[self.colZsunError84]=zSun84






	def add3DInfoReid2014(self,processTB):
		"""
		add 3D information according to the model of Reid 2014
		:param processTB:
		:return:
		"""
		d, dErrorLow , dErrorUp ,  l,b=self.getDLB_Reid2014(processTB )

		errorWithSys=  ( dErrorUp - dErrorLow )/2.   #dist std, Reid 2014, no distance error of systematic will be considered

		d16= d-errorWithSys
		d84=d+errorWithSys


		xSun,ySun,zSun=self.get3DPosition(d,l,b,withRespectTo=self.positionSun)
		xSun16,ySun16,zSun16=self.get3DPosition(d16,l,b,withRespectTo=self.positionSun)
		xSun84,ySun84,zSun84 =self.get3DPosition(d84,l,b,withRespectTo=self.positionSun)




		processTB[self.colXsun]= xSun
		processTB[self.colYsun]=ySun
		processTB[self.colZsun]=zSun

		processTB[self.colXsunError16]= xSun16
		processTB[self.colYsunError16]=ySun16
		processTB[self.colZsunError16]=zSun16

		processTB[self.colXsunError84]= xSun84
		processTB[self.colYsunError84]=ySun84
		processTB[self.colZsunError84]=zSun84


	def add3DInfoReid2019(self,processTB):
		"""
		add 3D information according to the model of Reid 2014
		:param processTB:
		:return:
		"""
		d, dErrorLow , dErrorUp ,  l,b=self.getDLB_Reid2019(processTB )




		errorWithSys=  ( dErrorUp - dErrorLow )/2.   #dist std, Reid 2014, no distance error of systematic will be considered

		d16= d-errorWithSys
		d84=d+errorWithSys


		xSun,ySun,zSun=self.get3DPosition(d,l,b,withRespectTo=self.positionSun)
		xSun16,ySun16,zSun16=self.get3DPosition(d16,l,b,withRespectTo=self.positionSun)
		xSun84,ySun84,zSun84 =self.get3DPosition(d84,l,b,withRespectTo=self.positionSun)


		processTB[self.colXsun]= xSun
		processTB[self.colYsun]=ySun
		processTB[self.colZsun]=zSun



		processTB[self.colXsunError16]= xSun16
		processTB[self.colYsunError16]=ySun16
		processTB[self.colZsunError16]=zSun16

		processTB[self.colXsunError84]= xSun84
		processTB[self.colYsunError84]=ySun84
		processTB[self.colZsunError84]=zSun84



	def get3DbyLBD(self,l,b,d):
		"""

		:param l:
		:param b:
		:param d:
		:return:
		"""


	def add3DInfo(self, processTB , addDisSys=True ):
		"""
		add 3D information to the goodTB, which will be used to draw face on, or 3D map, test this function
		:return:
		"""
		d,l,b=self.getDLB(processTB )

		#xGC, yGC, zGC=self.get3DPosition(d,l,b,withRespectTo=self.positionGC)
		#xSun,ySun,zSun=self.get3DPosition(d,l,b,withRespectTo=self.positionSun)

		#calculate d16
		disStd = processTB[self.disStdCol] / 1000.
		if addDisSys:
			errorWithSys = disStd ** 2 + (0.05 * d) ** 2
			errorWithSys = np.sqrt(errorWithSys)
		else:
			errorWithSys = disStd

		d16= d-errorWithSys
		d84=d+errorWithSys


		xSun,ySun,zSun=self.get3DPosition(d,l,b,withRespectTo=self.positionSun)
		xSun16,ySun16,zSun16=self.get3DPosition(d16,l,b,withRespectTo=self.positionSun)
		xSun84,ySun84,zSun84 =self.get3DPosition(d84,l,b,withRespectTo=self.positionSun)




		processTB[self.colXsun]= xSun
		processTB[self.colYsun]=ySun
		processTB[self.colZsun]=zSun

		processTB[self.colXsunError16]= xSun16
		processTB[self.colYsunError16]=ySun16
		processTB[self.colZsunError16]=zSun16

		processTB[self.colXsunError84]= xSun84
		processTB[self.colYsunError84]=ySun84
		processTB[self.colZsunError84]=zSun84





	def draw3Dmap(self):
		"""
		draw 3D map of molecular clouds
		:return:
		"""
		fig=plt.figure(figsize=(12,6))

		ax3D=fig.add_subplot(1,2,2,projection='3d')
		ax3D.view_init(10, 35)
		ax3D.scatter3D(0,self.R0, 0 , c='red' )
		X3D = self.goodTB[self.colXsun]
		Y3D = self.goodTB[self.colYsun]+self.R0
		Z3D = self.goodTB[self.colZsun]



		ax3D.scatter3D(X3D,Y3D,Z3D , c='blue', s=15 );
		ax3D.set_xlabel(r"Y (kpc)")
		ax3D.set_ylabel(r"X (kpc)")
		ax3D.set_zlabel(r"Z (kpc)")


		plt.show()



	def selectBesselSub(self,lRange,bRange,vRange,maxDis=4000):
		"""
		:param lRange:
		:param bRange:
		:param vRange:
		:param maxDis:
		:return:
		"""

		subL=doFITS.selectTBByColRange( self.besselTBProcessed , disTB.l, maxV=max(lRange) , minV=min(lRange)   )
		subB=doFITS.selectTBByColRange( subL , disTB.b, maxV=max(bRange) , minV=min(bRange)   )
		subV=doFITS.selectTBByColRange( subB ,  disTB.b, maxV=max(vRange) , minV=min(vRange)   )



		subD=doFITS.selectTBByColRange(  subV  ,  disTB.distance,maxV= maxDis )



		return subD




	def drawArmView(self, drawMaser=True, drawMaserLrange=[90, 180], drawMaserBrange=[-10, 10],
						 drawMaserVrange=[-90, 20], colorVrange=None, drawLrange=None, maxR=2.5, calChains=6):
		"""

		if masers were ddrawn, you need to provide the PPV range of masers

		:return:
		"""

		# simply draw a face on, and on background arm are provided
		xSun = 0.
		ySun = self.R0
		sunRadius = 0.02
		pSun = np.array([xSun, ySun])

		# plot Xs, Ys
		fig = plt.figure(figsize=(12, 6))

		import matplotlib as mpl

		# fig, axs = plt.subplots(nrows=1, ncols=2,  figsize=(12,6),sharex=True)
		rc('text', usetex=True)
		rc('font', **{'family': 'sans-serif', 'size': 12, 'serif': ['Helvetica']})

		mpl.rcParams['text.latex.preamble'] = [
			r'\usepackage{tgheros}',  # helvetica font
			r'\usepackage{sansmath}',  # math-font matching  helvetica
			r'\sansmath'  # actually tell tex to use it!
			r'\usepackage{siunitx}',  # micro symbols
			r'\sisetup{detect-all}',  # force siunitx to use the fonts
		]

		# custom coors with velocity
		ax0 = fig.add_subplot(1, 2, 1)
		axDisVel = fig.add_subplot(1, 2, 2)

		ax0.set_aspect("equal")

		Vs = self.goodTB["v_cen"]
		import matplotlib as mpl
		cmap = plt.cm.jet
		# norm = matplotlib.colors.BoundaryNorm(np.arange(0,30,0.1), cmap.N)
		if colorVrange is None:
			normV = mpl.colors.Normalize(vmin=min(Vs), vmax=max(Vs))

		else:
			normV = mpl.colors.Normalize(vmin=min(colorVrange), vmax=max(colorVrange))

		m = plt.cm.ScalarMappable(norm=normV, cmap=cmap)

		if 1:  # colorbar
			from mpl_toolkits.axes_grid1 import make_axes_locatable
			divider = make_axes_locatable(ax0)
			# cax1 = inset_axes(ax0,
			# width="3%",  # width = 10% of parent_bbox width
			# height="100%",  # height : 50%
			# loc=3,
			# bbox_to_anchor=(1.01, 0, 1, 1),
			# bbox_transform=ax0.transAxes,
			# borderpad=0.
			# )
			cax1 = divider.append_axes("right", size="3%", pad=0.05)
			cb = mpl.colorbar.ColorbarBase(cax1, norm=normV, cmap=cmap)

			cax1.set_ylabel(r"$V_{\rm LSR}$ ($\rm km\ s^{-1}$)")

		ax0.set_xlabel(r"X (kpc)")
		ax0.set_ylabel(r"Y (kpc)")

		if drawMaser:
			print "Drawing masers"

			subD = self.selectBesselSub(drawMaserLrange, drawMaserBrange, drawMaserVrange,maxDis=30000)
			self.drawFaceOnByTB(ax0, subD, m, marker="^")

		if 1:
			axDisVel.scatter( subD["vlsr"] ,subD["distance"],c='blue' ,s=6)
			axDisVel.scatter( self.goodTB["v_cen"] ,self.goodTB["distance"] , c='green' ,s=6)
		if 1:
			print "Draw clouds of  ", self.regionName
			self.drawFaceOnByTB(ax0, self.goodTB, m, marker="o")


		if 0:  #

			print "Draw clouds of  ", self.regionName
			self.drawFaceOnByTB(ax0, self.goodTB, m, marker="o")
			# self.drawVelByTB(ax,self.goodTB,m, lRange=[105,130])
			Vs, Verror, Ds, Derror = self.drawVelByTB(ax, self.goodTB, m)
			doLinearMCMC = linearMCMC3()

			para, paraError = doLinearMCMC.getParameters(Vs, Verror, Ds, Derror, calChains=calChains)

			print "y=ax+b, sigma", "[a, b, sigma]", "aError, bError,sigma"
			print para, paraError

			fitLineXs = np.linspace(min(Vs), max(Vs), 10)

			fitLineYs = fitLineXs * para[0] + para[1]

			formulaStr = r"$\frac{D}{\mathrm{[kpc]}}=\left(" + "{:.3f}".format(
				para[0]) + r"\frac{V_{\rm LSR}}{[\rm km \ s^{-1}]}+" + "{:.3f}".format(para[1]) + r"\right)$"

			ax.plot(fitLineXs, fitLineYs, 'b--', color='black', lw=0.8, alpha=0.8, label=formulaStr)
			ax.legend(loc=1, handlelength=1)

		# draw Line


		# self.drawVelByTB(ax,subD,m,marker="^",vErrorCol="e_VLSR",  lRange=[105,125])
		# draw distance of the first Galacti
		if 0:
			tableQ1 = Table.read("/home/qzyan/WORK/projects/maddalena/dendroDisPath/G2650/G2650goodDisTB.fit")
			self.add3DInfo(tableQ1, addDisSys=True)
			self.drawFaceOnByTB(ax0, tableQ1, m, marker="o")
		###############################################################################################################
		print "Draw lines and Markers"

		# draw circles
		testAngle = 60. / 180 * np.pi  # angle where the distances were put

		drawCircleDistances = np.arange(1, maxR + 1 ,  1 )  # [   0.5, 1.5, 2.5, 3.5]
		for radiusCircle in drawCircleDistances:
			# for radiusCircle in [0.25,  0.5,1,1.5,2,2.5,3,3.5]:
			"""
			"""
			# circle=plt.Circle((xSun,ySun),radius=3 )
			if drawLrange is None:
				drawA = np.linspace(-np.pi / 2, np.pi / 2, 200)

			else:
				drawA = np.linspace(np.deg2rad(min(drawLrange) - 90), np.deg2rad(max(drawLrange) - 90), 200)

			x = radiusCircle * np.cos(drawA)
			y = radiusCircle * np.sin(drawA)
			ax0.plot(x + xSun, y + ySun, '--', color='black', lw=0.5, alpha=0.5)

			if  radiusCircle == 1:
				continue

			if radiusCircle!=maxR:
				ax0.text(radiusCircle * np.cos(testAngle) + xSun + 0.05, radiusCircle * np.sin(testAngle) + ySun + 0.05,  "{}".format(radiusCircle), fontsize=10, ha='center', va="center")
			else:
				ax0.text(radiusCircle * np.cos(testAngle) + xSun + 0.05, radiusCircle * np.sin(testAngle) + ySun + 0.05,  "{} kpc".format(radiusCircle), fontsize=10, ha='center', va="center")


		maxDrawR = max(drawCircleDistances) + 0.5

		ax0.scatter(xSun, ySun, s=5, facecolors="black")

		ax0.scatter(xSun, ySun, s=20, facecolors='none', edgecolors='black', lw=0.3)

		ax0.text(xSun + 0.1, ySun - 0.5, "the Sun", fontsize=9, ha='center')

		# for drawL in [ 180, 190,200,210,220,230,240,250,260,270]:
		# for drawL in [25,30,35,40,45,50  ]:

		if drawLrange is not None:
			lArray = np.arange(min(drawLrange), max(drawLrange) + 10, 10)  #
		else:
			lArray = np.arange(0, 190, 10)
		for drawL in lArray:
			drawAngle = np.radians(drawL - 90)

			unitvector = np.array([np.cos(drawAngle), np.sin(drawAngle)])
			drawp1 = pSun - unitvector * maxDrawR

			drawp1_end = pSun - unitvector * sunRadius

			# ax0.plot( [drawp1_end[0],drawp1[0]],   [drawp1_end[1],drawp1[1]] ,'--',color='black',lw=0.5,alpha=0.5 )

			drawp2_start = pSun + unitvector * sunRadius

			drawp2 = pSun + unitvector * maxDrawR
			ax0.plot([drawp2_start[0], drawp2[0]], [drawp2_start[1], drawp2[1]], '--', color='black', lw=0.5,
					 alpha=0.5)

			unitvector = np.array([np.cos(drawAngle + 0.02), np.sin(drawAngle + 0.025)])

			refPosition = pSun + unitvector * (maxDrawR + 0.07)

			thetaShift = np.radians(drawL)
			# shiftVector=  np.array( [ np.cos(thetaShift), np.sin(thetaShift) ]  )
			# shiftVector=shiftVector *0.5
			drawp2Text = refPosition  # +shiftVector
			ax0.text(drawp2Text[0], drawp2Text[1], r"${}^\circ$".format(drawL), fontsize=8, rotation=drawL - 180,
					 va="center", ha="left", rotation_mode='anchor')

		##############################################################################################################

		# sc=ax0.scatter( Xs,Ys  ,  c =Vs, cmap=cmap,norm=normV, s=13 ,facecolors='none',   lw=0.5,marker="o")

		ax0.set_xlim([-0.7, maxR + 1])
		plt.subplots_adjust(wspace=0.3)
		# plt.show()
		plt.tight_layout(pad=0)

		saveTag = ""
		if drawMaser:
			saveTag = "withMaser"

		plt.savefig("ArmFaceOn{}_{}.png".format(self.regionName, saveTag), bbox_inches='tight', dpi=600)
		plt.savefig("ArmFaceOn{}_{}.pdf".format(self.regionName, saveTag), bbox_inches='tight')




	def drawFaceOnSimple(self, drawMaser=False,drawMaserLrange=[90,180],drawMaserBrange=[-10,10],drawMaserVrange= [-90,20],colorVrange=None, showXrange=None,showYrange=None, drawLrange = None, maxR= 2.5 , calChains=6 ,extraMWISPTB=None ,onlyFaceon=True    ):
		"""

		if masers were ddrawn, you need to provide the PPV range of masers

		:return:
		"""

		#simply draw a face on, and on background arm are provided
		xSun=0.
		ySun=self.R0
		sunRadius=0.02
		pSun=np.array( [xSun,ySun  ])


		#plot Xs, Ys

		if onlyFaceon:
			fig=plt.figure(figsize=(10,8))
		else:
			fig=plt.figure(figsize=(12,6))

		import matplotlib as mpl

		#fig, axs = plt.subplots(nrows=1, ncols=2,  figsize=(12,6),sharex=True)
		rc('text', usetex=True )
		rc('font', **{'family': 'sans-serif',  'size'   : 16,  'serif': ['Helvetica'] })

		mpl.rcParams['text.latex.preamble'] = [
			r'\usepackage{tgheros}',  # helvetica font
			r'\usepackage{sansmath}',  # math-font matching  helvetica
			r'\sansmath'  # actually tell tex to use it!
			r'\usepackage{siunitx}',  # micro symbols
			r'\sisetup{detect-all}',  # force siunitx to use the fonts
		]
 
		#custom coors with velocity



		if onlyFaceon:
			ax0 = fig.add_subplot(1, 1, 1)
			#ax = fig.add_subplot(1, 2, 2)
		else:
			ax0 = fig.add_subplot(1, 2, 1)
			ax = fig.add_subplot(1, 2, 2)



		ax0.set_aspect("equal")

		Vs=self.goodTB["v_cen"]
		import matplotlib as mpl
		cmap=plt.cm.jet
		#norm = matplotlib.colors.BoundaryNorm(np.arange(0,30,0.1), cmap.N)
		if colorVrange is None:
			normV = mpl.colors.Normalize(vmin= min(Vs), vmax= max(Vs) )

		else:
			normV = mpl.colors.Normalize(vmin= min(colorVrange), vmax= max( colorVrange ) )


		m = plt.cm.ScalarMappable(norm= normV , cmap=cmap)



		if 1: #colorbar
			from mpl_toolkits.axes_grid1 import make_axes_locatable
			divider = make_axes_locatable(ax0)
			#cax1 = inset_axes(ax0,
							  #width="3%",  # width = 10% of parent_bbox width
							  #height="100%",  # height : 50%
							  #loc=3,
							  #bbox_to_anchor=(1.01, 0, 1, 1),
							  #bbox_transform=ax0.transAxes,
							  #borderpad=0.
							  #)
			cax1 = divider.append_axes("right", size="3%", pad=0.05)
			cb = mpl.colorbar.ColorbarBase(cax1, norm=normV, cmap=cmap)

			cax1.set_ylabel(r"$V_{\rm LSR}$ ($\rm km\ s^{-1}$)" )



		ax0.set_xlabel(r"X (kpc)")
		ax0.set_ylabel(r"Y (kpc)")

		if not onlyFaceon:
			ax.set_xlabel(r"$V_{\rm LSR}$ ($\rm km\ s^{-1}$)")
			ax.set_ylabel(r"Distance (kpc)")

		print "Draw clouds of  ",self.regionName,"Total number of clouds:", len(self.goodTB )
		self.drawFaceOnByTB(ax0, self.goodTB,m,marker="o")
		#self.drawVelByTB(ax,self.goodTB,m, lRange=[105,130])

		if not onlyFaceon:#


			Vs,Verror, Ds,Derror=self.drawVelByTB(ax,self.goodTB,m   )

			doLinearMCMC=linearMCMC3()
			useReid2019= True
			if useReid2019: #draw lists
				pass
				############## use Reid2019

			# positiontion 1
			vList25, dList25 = self.getVDlist(l=105, b=-5, minV=-45, maxV=10, useReid2019=useReid2019)
			ax.plot(vList25, dList25,   lw=0.8, label= self.getLBLabel(105,-5)  ,zorder=1)

			# positiontion 2
			vList25, dList25 = self.getVDlist(l=105, b=5, minV=-45, maxV=10, useReid2019=useReid2019)
			ax.plot(vList25, dList25,   lw=0.8, label=  self.getLBLabel(105, 5)  ,zorder=1)

			# positiontion 3
			vList25, dList25 = self.getVDlist(l=150, b=-5, minV=-45, maxV=10, useReid2019=useReid2019)
			ax.plot(vList25, dList25,  lw=0.8, label= self.getLBLabel(150, -5)   ,zorder=1)

			# positiontion 4
			vList25, dList25 = self.getVDlist(l=150, b=5, minV=-45, maxV=10, useReid2019=useReid2019)
			ax.plot(vList25, dList25,   lw=0.8, 	label= self.getLBLabel(150,  5)   ,zorder=1)

			# positiontion 5
			vList25, dList25 = self.getVDlist(l=125, b=0, minV=-45, maxV=10, useReid2019=useReid2019)
			ax.plot(vList25, dList25,   lw=0.8, label=  self.getLBLabel(125,  0)  ,zorder=1)
			#ax.set_xlabel(r"$V_{\rm LSR}$ (km s$^{-1}$)")
			#ax.set_ylabel(r"Kinematic distance (kpc)")

			ax.legend(loc=3,handlelength=0.5)
			ax.set_ylim([-0.1,2.5])
			ax.set_xlim([-45,   9  ])

			if useReid2019:

				at = AnchoredText('Reid et al. (2019)', loc=1, frameon=False)
				ax.add_artist(at)

			else:
				at = AnchoredText('Reid et al. (2016)', loc=1, frameon=False)
				ax.add_artist(at)

		if 0:
			para,paraError=doLinearMCMC.getParameters(Vs,Verror,Ds,Derror,calChains= calChains )

			print "y=ax+b, sigma","[a, b, sigma]", "aError, bError,sigma"
			print para,paraError

			fitLineXs= np.linspace(min(Vs),max(Vs), 10)

			fitLineYs= fitLineXs*para[0]+para[1]

			formulaStr=r"$\frac{D}{\mathrm{[kpc]}}=\left("+ "{:.3f}".format(para[0]) +r"\frac{V_{\rm LSR}}{[\rm km \ s^{-1}]}+"+"{:.3f}".format(para[1]) +r"\right)$"

			ax.plot(fitLineXs,fitLineYs, 'b--', color='black',lw=0.8 ,alpha=0.8,label= formulaStr)
			ax.legend(loc=1,handlelength=1 )




		if extraMWISPTB is not None:
			print "draw Extra tables  "

			for eachTB in extraMWISPTB:
				self.drawFaceOnByTB(ax0, eachTB, m, marker="o")
		if 0:
			self.drawSFbook(ax0)


		if drawMaser:
			print "Drawing masers"

			subD=self.selectBesselSub( drawMaserLrange,drawMaserBrange, drawMaserVrange  )
			self.drawFaceOnByTB(ax0, subD,m,marker="^")

			#self.drawVelByTB(ax,subD,m,marker="^",vErrorCol="e_VLSR",  lRange=[105,125])
		#draw distance of the first Galacti		
		if 0:
			tableQ1= Table.read("/home/qzyan/WORK/projects/maddalena/dendroDisPath/G2650/G2650goodDisTB.fit")
			self.add3DInfo(tableQ1,addDisSys=True)
			self.drawFaceOnByTB(ax0, tableQ1,m,marker="o")
		###############################################################################################################
		print "Draw lines and Markers"

		# draw circles


		drawCircleDistances =  np.arange(0.5,maxR+0.5,0.5) # [   0.5, 1.5, 2.5, 3.5]
		for radiusCircle in drawCircleDistances:
		#for radiusCircle in [0.25,  0.5,1,1.5,2,2.5,3,3.5]:
			"""
			"""
			#circle=plt.Circle((xSun,ySun),radius=3 )
			if drawLrange is   None:
				drawA=np.linspace(-np.pi/2, np.pi/2,200)

			else:
				drawA = np.linspace( np.deg2rad(min(drawLrange) -90 ) , np.deg2rad( max(drawLrange)-90  ), 200)

			x=radiusCircle*np.cos(drawA)
			y=radiusCircle*np.sin(drawA)
			ax0.plot( x +xSun ,   y +ySun,'--',color='black',lw=0.5,alpha=0.5 )

			#if radiusCircle!=2:
			testAngle= 90./180*np.pi #angle where the distances were put
			ax0.text( radiusCircle*np.cos(testAngle) +xSun+0.05,  radiusCircle*np.sin(testAngle)+ySun+0.05, "{} kpc".format( radiusCircle) ,fontsize=10 ,ha='center',va="center")

		maxDrawR= max(drawCircleDistances)+0.5

		ax0.scatter(xSun,ySun,s=5,facecolors="black")

		ax0.scatter(xSun,ySun,s=20,facecolors='none', edgecolors='black',lw=0.3)

		ax0.text( xSun  ,ySun-0.15,"the Sun",fontsize=11 ,ha='center'  )


		#for drawL in [ 180, 190,200,210,220,230,240,250,260,270]:
		#for drawL in [25,30,35,40,45,50  ]:


		if drawLrange is not None:
			lArray=  np.arange(min(drawLrange),max(drawLrange)+10,10)  #
		else:
			lArray= np.arange(0,190,10)
		for drawL in  lArray  :

			drawAngle=np.radians(drawL-90)

			unitvector=np.array( [ np.cos(drawAngle), np.sin(drawAngle) ]  )
			drawp1= pSun-unitvector*maxDrawR

			drawp1_end= pSun-unitvector*sunRadius


			#ax0.plot( [drawp1_end[0],drawp1[0]],   [drawp1_end[1],drawp1[1]] ,'--',color='black',lw=0.5,alpha=0.5 )

			drawp2_start= pSun+unitvector*sunRadius

			drawp2= pSun+unitvector*maxDrawR
			ax0.plot( [drawp2_start[0],drawp2[0]],   [drawp2_start[1],drawp2[1]] ,'--',color='black',lw=0.5,alpha=0.5 )

			unitvector=np.array( [ np.cos(drawAngle+0.02), np.sin(drawAngle+0.025 ) ]  )



			refPosition= pSun+unitvector * (maxDrawR+0.07 )

			thetaShift=  np.radians( drawL  )
			#shiftVector=  np.array( [ np.cos(thetaShift), np.sin(thetaShift) ]  )
			#shiftVector=shiftVector *0.5
			drawp2Text=refPosition #+shiftVector
			ax0.text(drawp2Text[0] ,drawp2Text[1]  , r"${}^\circ$".format(drawL),fontsize=10, rotation=drawL-180, va="center", ha="left", rotation_mode='anchor'  )






		##############################################################################################################

		#sc=ax0.scatter( Xs,Ys  ,  c =Vs, cmap=cmap,norm=normV, s=13 ,facecolors='none',   lw=0.5,marker="o")




		ax0.set_xlim([-0.3, maxR+1 ])

		#ax0.set_xlim([-2,  2 ])
		if showXrange is not None:
			ax0.set_xlim( showXrange )
		if showYrange is not None:
			ax0.set_ylim( showYrange )




		plt.subplots_adjust(wspace=0.3 )
		#plt.show()
		plt.tight_layout(pad=0)

		saveTag=""
		if drawMaser:
			saveTag="withMaser"

		plt.savefig( "faceOn{}{}.png".format(self.regionName, saveTag),bbox_inches='tight',dpi=600)
		plt.savefig( "faceOn{}{}.pdf".format(self.regionName, saveTag ), bbox_inches='tight')


	def getLBLabel(self,l,b):

		return r"$\left(\mathit{{l}},\mathit{{b}}\right)= \left({:.0f}^\circ, {:.0f}^\circ\right)$".format(l, b)

	def drawFaceOnByTB_useReid2014(self,ax0, drawTB, marker="o"):
		"""
		#draw the 
		:param ax0:
		:param marker:
		:return:
		"""


		#omit error
		Xs=drawTB[self.colXsun]
		Ys= drawTB[self.colYsun] +self.R0

		#sc=ax0.scatter(Xs ,Ys  , c=Vs, edgecolors = m.to_rgba(Vs), cmap=cmap,norm=normV, s= 8 ,facecolors='none',  lw=0.1,  marker="o", )

		sc=ax0.scatter(Xs ,Ys  ,   edgecolors = "blue",   s= 3 ,facecolors='none',  lw= 0.3 ,  marker= marker  ,alpha=0.8 )

	def drawSFbook(self,ax0):
		"""

		:param ax0:
		:return:
		"""

		from astropy.io import ascii
		tbZucker = ascii.read("/home/qzyan/WORK/dataDisk/cataLog/Handbook_Distances_Zucker2019.dat")


		l = tbZucker["l"]
		b = tbZucker["b"]
		d = tbZucker["d50"] #pc

		Xs, Ys, Zs = self.get3DPosition(d, l, b, withRespectTo=self.positionGC )

		sc=ax0.scatter(Xs ,Ys  ,   edgecolors = 'black' ,   s= 6 ,facecolors='none',  lw= 0.5,  marker= 'o',alpha=0.8   )




	def drawFaceOnByTB(self,ax0,maserSubTB,m,marker="o" ):

		"""
		#draw the face on distribution containing maser spots
		:param maserSubTB:
		:return:
		"""
		self.add3DInfo(maserSubTB,addDisSys=False)

		Xs=maserSubTB[self.colXsun]
		Ys= maserSubTB[self.colYsun] +self.R0

		if m is not None:
			if "v_cen" in maserSubTB.colnames:
				Vs= maserSubTB["v_cen"]
			else:
				Vs= maserSubTB[ disTB.cloudVlsr  ]

		#sc=ax0.scatter(Xs ,Ys  , c=Vs, edgecolors = m.to_rgba(Vs), cmap=cmap,norm=normV, s= 8 ,facecolors='none',  lw=0.1,  marker="o", )

		if "sum" in maserSubTB.colnames:
			massList,massError= self.getMassList(maserSubTB)
			sizeArray=  massList/1000.+1
		else:
			sizeArray=8
		if m is None:
			sc=ax0.scatter(Xs ,Ys  ,   edgecolors =  "blue",   s= sizeArray ,facecolors='none',  lw= 0.5,  marker= marker   )
		else:

			sc = ax0.scatter(Xs, Ys, edgecolors=m.to_rgba(Vs), s=sizeArray, facecolors='none', lw=0.5, marker=marker)

		X16s= maserSubTB[self.colXsunError16]
		Y16s= maserSubTB[self.colYsunError16]+self.R0
		X84s= maserSubTB[self.colXsunError84]
		Y84s= maserSubTB[self.colYsunError84]+self.R0

		for i in range(len(X16s)):
			print i

			if m is not None:
				ax0.plot( [X16s[i],X84s[i]], [Y16s[i],Y84s[i]],color=m.to_rgba(Vs[i]),lw=0.8     )
			else:
				ax0.plot([X16s[i], X84s[i]], [Y16s[i], Y84s[i]], color= "blue" , lw=0.8)

	def drawVelByTB(self,ax,drawTB,m,marker="o",withSysError=True,vErrorCol="v_rms", lRange=None):
		"""
		draw the  velocity distance relationship of draw TB
		:param ax:
		:param drawTB:
		:param m:
		:param marker:
		:return:
		"""
		cmap=plt.cm.jet
		#norm = matplotlib.colors.BoundaryNorm(np.arange(0,30,0.1), cmap.N)

		normV = mpl.colors.Normalize(vmin=  105 , vmax= 150 )


		m = plt.cm.ScalarMappable(norm= normV , cmap=cmap)

		if lRange is not None:
			drawTB=doFITS.selectTBByColRange( drawTB, "l", minV=min(lRange),maxV=max(lRange) )

		Verror = drawTB[ vErrorCol]



		Derror = self.getDisErrorWithSys(withSysError=withSysError)
		Ds = drawTB["distance"] / 1000.
		if "v_cen" in drawTB.colnames:
			Vs= drawTB["v_cen"]
		else:
			Vs= drawTB[ disTB.cloudVlsr  ]

		Ls= drawTB["x_cen"]

		for i in range(len(Vs)):
			#ax.errorbar([Vs[i]], [Ds[i]], yerr=[[Derror[i]], [Derror[i]]], xerr=[[Verror[i]], [Verror[i]]], color=m.to_rgba(Vs[i]), fmt= marker, capsize=1.2, elinewidth=0.7, lw=0.5, markersize=2.5)
			#ax.errorbar([Vs[i]], [Ds[i]], yerr=[[Derror[i]], [Derror[i]]], xerr=[[Verror[i]], [Verror[i]]], color=  m.to_rgba(Ls[i]) , fmt= marker, capsize=1.2, elinewidth=0.7, lw=0.5, markersize=2.5,zorder = 2)
			ax.errorbar([Vs[i]], [Ds[i]], yerr=[[Derror[i]], [Derror[i]]], xerr=[[Verror[i]], [Verror[i]]], color=   "black" , fmt= marker, capsize=1.2, elinewidth=0.7, lw=0.5, markersize=2.5,zorder = 2)

		return Vs,Verror, Ds,Derror


	def getDlistByL(self,l=25):
		"""
		return a list of distance in kpc, with a list of velocity range
		:param l:
		:return:
		"""
		velocityList=  np.linspace(2 ,30,50)

		disList=[]

		for eachV in velocityList:

			dis2014,_,_=self.ReidA5(l,0,eachV)
			disList.append( dis2014[0]  )

		return np.asarray(disList ), velocityList


	def drawVelDis(self,cloudTB):
		"""
		Draw distances, vs velocity
		:return:
		"""

		Vs=[]
		Ds=[]

		Ls=[]

		Bs=[]
		Ms=[]

		cloudTB=Table.read(cloudTB)

		Vrms=[]
		for eachC in self.TB:

			l=eachC[disTB.l]  #cloudVlsr
			b=eachC[disTB.b]
			v=eachC[disTB.cloudVlsr]
			gaiaDis =  eachC[disTB.distance]/1000.

			Vs.append(v)
			Ds.append(gaiaDis)
			Ls.append(l)
			Ms.append( doDis.calMassByRow(eachC)  )

			ID=int( eachC[disTB.sourceName].split("oud")[1] )



			cloudRow= cloudTB[ cloudTB["_idx"]==ID  ][0]


			Vrms.append( cloudRow["v_rms"]   )

		#l25




		fig=plt.figure(figsize=(12,6))
		rc('text', usetex=True )
		rc('font', **{'family': 'sans-serif',  'size'   : 13,  'serif': ['Helvetica'] })
		mpl.rcParams['text.latex.preamble'] = [
			r'\usepackage{tgheros}',  # helvetica font
			r'\usepackage{sansmath}',  # math-font matching  helvetica
			r'\sansmath'  # actually tell tex to use it!
			r'\usepackage{siunitx}',  # micro symbols
			r'\sisetup{detect-all}',  # force siunitx to use the fonts
		]
		ax=fig.add_subplot(1,2,1)
		ax.scatter( Ds,Vs,s=15,c="black")

		ax.set_ylabel(r"Velocity (km/s)")
		ax.set_xlabel(r"Distance (kpc)")


		for draL in [25,30, 35,40,45,50]:


			vList25,dList25 = self.getVDlist(  l=draL  )

			ax.plot(dList25, vList25,  label=r"$l={}^\circ$".format(draL) )

		ax.legend()


		ax1=fig.add_subplot(1,2,2)

		ax1.scatter(np.log10( Ms), np.log10(  Vrms ) )

		ax1.set_xlabel(r"  Log Mass")
		ax1.set_ylabel(r" Log Vrms ")

		plt.savefig( "velDis{}.pdf".format(self.regionName), bbox_inches='tight')

	def getCloudArea(self):
		"""
		use goodTB
		:return:
		"""

		areaList=[] #in pc unite









	def getMassList(self,processTB=None ):
		"""

		#consider the error of flux and distance?

		:param tbRow:
		:return:
		"""


		if processTB==None:
			processTB = self.goodTB

		massList=[]
		massErrorList = []


		for eachC in  processTB :

			l=eachC["l"]
			b=eachC["b"]
			v=eachC["vlsr"]

			#disList,_,_=self.ReidA5(l,b,v)

			COsum = eachC["sum"]

			errorWithSys= np.sqrt(   eachC["disStd"]**2  +  ( eachC["distance"] *0.05 )**2  )

			cloudMass, massError = self.calmassByXfactor(COsum * self.dv, distance=eachC["distance"],errorSys= errorWithSys )  # ...
			#cloudMass, massError = self.calmassByXfactor(COsum * self.dv, distance= disList [0]*1000,errorSys= errorWithSys )  # ...

			massList.append(cloudMass)
			massErrorList.append( massError  )





		return np.asarray(massList) , np.asarray( massErrorList )

	def getCloudID(self,tb):
		idList=[]
		for eachC in tb:
			cloudName= eachC["sourceName"]
			id=cloudName.split("Cloud")[-1]
			idList.append( int(  id ))

		return np.asarray(idList)

	def getCloudIDByRow(self,eachC):

		cloudName= eachC["sourceName"]
		id=cloudName.split("Cloud")[-1]
		return  int(id)





	def getGoodTBDBSCAN(self  ):
		"""
		return a table containing tb
		:param tb:
		:return:
		"""
		tb=self.goodTB.copy()
		tbGoodIndex= self.getCloudID(tb)


		tb["_idx"]=tbGoodIndex
		empTB=myFITS.getTBStructure(self.dbscanTB)

		for eachC in   tb: #self.scimesTB:


			tbID=eachC["_idx"]

			select= self.dbscanTB["_idx"]==tbID

			selectTB=  self.dbscanTB[select]


			empTB.add_row(selectTB[0])

		return empTB


	def test(self):

		farCloud=self.TB[ self.TB[disTB.distance]>1000  ]

		print  farCloud[disTB.sourceName]
		print  farCloud[disTB.cloudVlsr]

	def getDisLatexByRow(self, eachC , emptyRow=None ):

		# disCatTB=Table.read(self.distanceCatalog)
		# use new

		nameStr = eachC[disTB.Note]

		lbvStr = "{:>8.3f} & {:>7.3f} & {:>7.1f}".format(eachC[disTB.l], eachC[disTB.b], eachC[disTB.cloudVlsr])

		# disError_=eachC[disTB.MBM_disError]

		distance = eachC[disTB.distance]
		lowerDisError_ = eachC[disTB.disHPDLow]
		upperDisError_ = eachC[disTB.disHPDUp]




		disStd = eachC[disTB.disStd]

		trueError = np.sqrt(disStd ** 2 + (0.05 * distance) ** 2)

		distanceMCMCStr = "${:>4.0f}_{{-{:>4.0f}}}^{{+{:>4.0f}}}$".format(distance, lowerDisError_, upperDisError_)

		distanceWithSystematicErr = "${:>4.0f}\pm{:>3.0f}$".format(distance, trueError)

		# disLatex=" {:>4s}  -{:>3s}  {:>3s}  ".format( distanceMCMCStr ,   eachC[disTB.onStarN]   )

		onStarNstr = "{:>5.0f}".format(eachC[disTB.onStarN])

		if emptyRow is not None:
 
			emptyRow["N"] =  eachC[disTB.onStarN]

			emptyRow["gaiaDistance"] = distance
			emptyRow["gaiaDistanceHPD95Low"] = lowerDisError_
			emptyRow["gaiaDistanceHPD95Up"] = upperDisError_

		return distanceMCMCStr, onStarNstr  # disLatex  #distanceMCMCStr , distanceWithSystematicErr

		return " {:>4s} {:>3s}  {:>3s}  ".format("", "--", ""), "{:>5s}".format( "--")  # "          --          ","&  --   "

	def calmassByXfactor(self, coInt, distance, errorSys=None, xFactor=2.0e20):

		"""
		The unit of coInt must be K km/s,
		distance pc
		not calDis, calMass

		the errorDis, need to include systematic errors

		"""


		NH2 = coInt * xFactor # cm-2 #

		# distance=2200 # pc

		# parsecToMeter= 3.0857e16 #m
		# length1=np.radians(degSize1)*self.distance*self.parsecToMeter*100. #cm
		# length2=np.radians(degSize2)*self.distance*self.parsecToMeter*100.

		# should use single pix for 12CO

		length1 = np.radians(30. / 60. / 60.) * distance * self.parsecToMeter * 100.  # cm
		# length2=np.radians(degSize2)*self.distance*self.parsecToMeter*100.

		mu = 2.72 # 1.36

		Mh2 = 3.35e-27  # kg
		solarMass = 1.9891e30  # kg
		# s=np.pi*length1*length1
		s = length1 * length1

		coreSolar = s * NH2 * Mh2 * mu / solarMass

		if errorSys is None:

			return coreSolar
		else:

			erorMass=   np.radians(30. / 60. / 60.) *  errorSys *  self.parsecToMeter * 100.*  2* length1* NH2 * Mh2 * mu / solarMass

			return coreSolar, erorMass

	def getLatexName(self, cName):
		"""
		"""

		if "+" in cName:
			return cName.split('+')[0] + "$+$" + cName.split('+')[1].zfill(4)
		if "-" in cName:
			return cName.split('-')[0] + "$-$" + cName.split('-')[1].zfill(4)
		return cName

	def combineline(self, latexName, lbvStr, projectAreaStr, disStr, onStarNstr, signalLevel, disCutStr, massStr,
					 distanceReid2019Str, useAV, Note):
		"""

		"""
		latexRowExtra = " {} & {} & {} & {} & {} & {} & {}  & {} & {} & {} & {:<7s} \\\\ ".format(latexName, lbvStr,
																							 projectAreaStr, disStr,
																							 onStarNstr, signalLevel,
																							 disCutStr, massStr,
																							  distanceReid2019Str ,useAV, Note)
		if signalLevel == "":
			latexRowExtra = " {} & {} & {} & {} & {} & {}  & {} & {} & {} & {:<7s} \\\\ ".format(latexName, lbvStr,
																							projectAreaStr, disStr,
																							onStarNstr, disCutStr,
																							massStr,   distanceReid2019Str, useAV,
																							Note)

		return latexRowExtra

	def printCatlog(self  ):
		"""
		print catalog of good list
		:param goodTBName:
		:return:
		"""

		printTB= self.goodTB #Table.read( self.goodTBName )

		latexFile=self.regionName+"LatexDisTB.txt"

		latexFile=  open( latexFile ,"w")

		printTB.sort("l")



		publishCat= Table()
		publishCat["name"]=  np.array([],dtype= 'S20' )

		publishCat["l"]=  np.array([],dtype=np.float)
		publishCat["l"].unit = "deg"
		publishCat["b"]= np.array([],dtype=np.float)
		publishCat["b"].unit = "deg"

		publishCat["vlsr"]= np.array([],dtype=np.float)
		publishCat["vlsr"].unit="km/s"

		publishCat["area"]= np.array([],dtype=np.float)
		publishCat["area"].unit="square deg"

		publishCat["gaiaDistance"]= np.array([],dtype=np.float)
		publishCat["gaiaDistance"].unit="pc"


		publishCat["gaiaDistanceHPD95Low"]= np.array([],dtype=np.float)
		publishCat["gaiaDistanceHPD95Low"].unit="pc"

		publishCat["gaiaDistanceHPD95Up"]= np.array([],dtype=np.float)
		publishCat["gaiaDistanceHPD95Up"].unit="pc"


		publishCat["N"]= np.array([],dtype=np.int)

		publishCat["COcut"]= np.array([],dtype=np.float)
		publishCat["COcut"].unit="K"


		publishCat["Dcut"]= np.array([],dtype=np.float)
		publishCat["Dcut"].unit="pc"

		publishCat["mass"]= np.array([],dtype=np.float)
		publishCat["mass"].unit="solarMass"




		##############
		publishCat["DReid2019"]= np.array([],dtype=np.float)
		publishCat["DReid2019"].unit="kpc"
		publishCat["DStdLowReid2019"]= np.array([],dtype=np.float)
		publishCat["DStdLowReid2019"].unit="kpc"
		publishCat["DStdUpReid2019"]= np.array([],dtype=np.float)
		publishCat["DStdUpReid2019"].unit="kpc"



		publishCat["extinction"]=  np.array([], dtype= 'S10' )

		publishCat["Note"]=  np.array([],  dtype= 'S20' )


		publishCatCopy=publishCat.copy()

		publishCatCopy.add_row()

		emptyRow= publishCatCopy[0]

		for eachCloud in printTB:
			#Note=""
			tmpEmptyTB= publishCatCopy.copy()
			emptyRow = tmpEmptyTB[0]
			# eachCloud is the #math the cloud with the raw catalog of SCIMES

			sourceName= eachCloud["sourceName"]
			cloudID= eachCloud["_idx"] #int(   sourceName.split("oud")[1]   )
			#eachCloud= self.scimesTB[ self.scimesTB["_idx"]==cloudID ]

			#l = np.float(eachC["l"])
			#b = np.float(eachC["b"])
			#v = np.float(eachC["vlsr"])    # to km/s
			#disReid, rGC, alpha = self.ReidA5(l, b, v)



			l = np.float(eachCloud["x_cen"])
			b = np.float(eachCloud["y_cen"])
			v = np.float(eachCloud["v_cen"])    # to km/s
			#disReid, rGC, alpha = self.ReidA5(l, b, v)
			#disReid  = self.Reid2016(l, b, v) ########


			disReid2019 = self.Reid2019(l, b, v)




			emptyRow["l"] = l
			emptyRow["b"] = b
			emptyRow["vlsr"] = v


			#emptyRow["kinematicDis2016"] =  disReid[0]
			#emptyRow["kinematicStdLow2016"] =   abs(disReid[1])
			#emptyRow["kinematicStdUp2016"] =   disReid[2]

			emptyRow["DReid2019"] =  disReid2019[0]
			emptyRow["DStdLowReid2019"] =   abs(disReid2019[1])
			emptyRow["DStdUpReid2019"] =   disReid2019[2]


			#print l,b,v


			#distanceReidStr = "${:>4.2f}_{{-{:>4.2f}}}^{{+{:>5.2f}}}$".format(disReid[0], abs(disReid[1]), disReid[2])
			distanceReid2019Str = "${:>4.2f}_{{-{:>4.2f}}}^{{+{:>5.2f}}}$".format(disReid2019[0], abs(disReid2019[1]), disReid2019[2])
			
			projectArea = eachCloud["area_exact"] / 3600.  # deg^2
			emptyRow["area"] = projectArea



			projectAreaStr = "{:>5.2f}".format(projectArea)
			cloudName = self.getCloudNameByLB(l, b)
			emptyRow["name"] = cloudName

			latexDis,onStarNstr = self.getDisLatexByRow( eachCloud, emptyRow=emptyRow  )

			lbvStr= " {:>8.3f} & {:>7.3f} & {:>7.1f} ".format( l   , b  , v   )

			#Note="Dendrogram"
			Note=""

			if cloudName=="G134.8+01.4":
				Note="W4"


			#########calculate mass

			disCutStr = eachCloud[disTB.cutDistanceUpper]
			disCutStr = "{:>4.0f}".format(disCutStr)
			emptyRow["Dcut"] = eachCloud[disTB.cutDistanceUpper]

			COsum = eachCloud[disTB.COsum]
			SL = eachCloud[disTB.signalLevel]


			cloudMass =  eachCloud[self.massCol] #...

			# print cloudMass," solar mass"
			emptyRow["mass"] = cloudMass
			massStr = "{:>5.1f}".format(cloudMass / 1000. )

			signalLevel = "{:>2.0f}".format(SL)
			signalLevel = ""#"Do not show "
			emptyRow["COcut"] = SL

			#AV AG

			print cloudName



			path="/home/qzyan/WORK/projects/maddalena/dendroDisPath/{}/goodFigures/".format(self.regionName)
			files=glob.glob(path+cloudName+"*")

			figureFile= files[0]
			if "AV" in figureFile:
				useAV="$A_V$"
				emptyRow['extinction'] = "AV"

			else:
				useAV="$A_G$"
				emptyRow['extinction'] = "AG"


			if cloudName=="G029.6+03.7":
				Note="W40"


			if cloudName=="G115.6-02.7":
				Note="L1265"


			#if cloudName=="G122.4-00.6":
				#Note="L1302"


			if cloudName=="G122.4-00.6":
				Note="L1287,L1302"

			emptyRow['Note'] = Note

			publishCat.add_row(emptyRow)

			latexRow = self.combineline(self.getLatexName(cloudName), lbvStr, projectAreaStr, latexDis, onStarNstr, signalLevel,
										disCutStr, massStr,  distanceReid2019Str, useAV, Note)

			#print latexRow

			latexFile.write(latexRow + " \n")


		latexFile.close()
		publishCat.remove_column("COcut")
		publishCat.write(self.regionName+"PubDisCat.fit",overwrite=True)
	def removeEdgeCriteria(self):


		tfList=[]


		for eachC in self.goodTB:


			tfList.append(   eachC["Note"] not in self.edgeCloudNameList   )


		return np.asarray(tfList)


	def drawSideView(self):
		print "Draing sider views"

		fig=plt.figure(figsize=(8,6))
		rc('text', usetex=True )
		rc('font', **{'family': 'sans-serif',  'size'   : 13,  'serif': ['Helvetica'] })
		mpl.rcParams['text.latex.preamble'] = [
			r'\usepackage{tgheros}',  # helvetica font
			r'\usepackage{sansmath}',  # math-font matching  helvetica
			r'\sansmath'  # actually tell tex to use it!
			r'\usepackage{siunitx}',  # micro symbols
			r'\sisetup{detect-all}',  # force siunitx to use the fonts
		]


		mass=self.getMassList(  )
		ax=fig.add_subplot(1,1,1)

		edgeCloud=self.removeEdgeCriteria()

		if 1:#draw mass

			b=self.goodTB["b"]
			l=self.goodTB["l"]
			dis=self.goodTB["distance"]
			vs=self.goodTB["v_cen"]
			vrms=  self.goodTB["v_rms"]*2.35
			ax.scatter(  vrms[edgeCloud]  ,mass[edgeCloud]  )
			ax.set_xscale("log")
			ax.set_yscale("log")

		if 0:
			b=goodTB["b"]
			l=goodTB["l"]

			centerL=125
			dis=goodTB["distance"]
			bRad=np.deg2rad(b)
			y=np.sin(bRad)*dis
			theta=np.deg2rad( l-centerL )
			projTheta=dis*np.cos(bRad)* np.cos(theta)
			ax.scatter( projTheta , y ,s=10)


		plt.savefig( "siderView{}.pdf".format(self.regionName), bbox_inches='tight')

	def getCloudSize(self, TB,inPC=False, dis=None  ):

		beamSize =  49./50


		if not inPC :
			return np.sqrt(4 * TB["area_exact"] / np.pi - beamSize ** 2)

		if inPC and dis is None:

			print "Error .............. you need to provide distance if you want pc results"
			return

		size =  np.sqrt(4 * TB["area_exact"] / np.pi - beamSize ** 2)
		sizeInPC = np.deg2rad(size / 60.) * dis
		return sizeInPC

	def appendLineWidth(self , TBWithLinewidth ):
		"""

		:param TBWithLinewidth:
		:return:
		"""

		tbLW=Table.read( TBWithLinewidth )

		lwList=[]

		for eachR in self.goodTB:

			id=self.getCloudIDByRow(eachR )

			sourceName="Cloud{}".format(id)

			lwValue=  Table( tbLW[ tbLW["sourceName"] == sourceName  ])[0][ self.lwCol ]
			lwList.append( lwValue )


		self.goodTB[ self.lwCol ] =  lwList


	def getRgc(self,TB):
		"""
		get the dstance to the Galactic center
		:param TB:
		:return:
		"""

		X3D = TB[self.colXsun]
		Y3D = TB[self.colYsun]+self.R0
		Z3D = TB[self.colZsun]

		Rgc= np.sqrt( X3D**2+Y3D**2 + Z3D**2  )
		#Rgc= np.sqrt( X3D**2+Y3D**2   )

		return Rgc

	def drawMassSize(self,drawSizeVel=False,BFFcor=True ):
		"""

		draw the relationship beteen mass and velocity

		:return:
		"""

		#self.appendLineWidth( "disQ2goodDisT_LW.fit")


		#save a copy of the goodTable

		#self.goodTB.write(self.regionName+"mergeGoodTB.fit")

		fig=plt.figure(figsize=(8,6))
		rc('text', usetex=True )
		rc('font', **{'family': 'sans-serif',  'size'   : 16,  'serif': ['Helvetica'] })
		mpl.rcParams['text.latex.preamble'] = [
			r'\usepackage{tgheros}',  # helvetica font
			r'\usepackage{sansmath}',  # math-font matching  helvetica
			r'\sansmath'  # actually tell tex to use it!
			r'\usepackage{siunitx}',  # micro symbols
			r'\sisetup{detect-all}',  # force siunitx to use the fonts
		]



		ax=fig.add_subplot(1,1,1)

		#self.goodTB = self.goodTB[ self.goodTB[self.]  ]

		edgeCloud=self.removeEdgeCriteria()

		#self.goodTB["aaaaaaaaaa"] =self.goodTB["disStd"]


		if 1:#draw mass

			dis=self.goodTB["distance"]
			#vrms= # self.goodTB["pixN"] # self.goodTB["lineWidth"]*2.35
			#using Rgc as color
			Rgc = self.getRgc(self.goodTB)

			print self.goodTB.colnames

			Rgc = Rgc  #np.log(  self.goodTB["sum"] /self.goodTB["pixN"] )  #self.getRgc(self.goodTB)  # np.sqrt( self.goodTB["area_exact"] ) #test other qulities #self.getRgc(self.goodTB)

			meanSNR= self.goodTB["sum"] /self.goodTB["pixN"] /0.5
			BFFsens= (meanSNR-2.329)**2/(  meanSNR-2.329+0.423 )



			cmap = plt.cm.jet
			normV = mpl.colors.Normalize(vmin=np.min(Rgc), vmax=np.max(Rgc))

			m = plt.cm.ScalarMappable(norm=normV, cmap=cmap)

			################################
			if self.lwCol in self.goodTB.colnames:

				vrms=  self.goodTB[self.lwCol]
			else:
				vrms = self.goodTB["v_rms"]

			size=self.getCloudSize(self.goodTB )
			sizeInPC= np.deg2rad(size/60.)*  dis

			sizeInbeam= size/self.beamSize
			if BFFcor:
				fillingFactors= sizeInbeam**2/(sizeInbeam+2.573 )**2
			else:
				fillingFactors=  1

			mass = self.goodTB[self.massCol] / fillingFactors#  /fillingFactors
			massError =self.goodTB[self.massColError ] /  fillingFactors

			ax.set_xscale("log")
			ax.set_yscale("log")





			
			if drawSizeVel:


				Xs =  sizeInPC[edgeCloud]
				Ys =   vrms[edgeCloud]



				#YsError= massError[edgeCloud]
				#ax.scatter(  Xs ,  Ys  ,  s= 8 ,  color="blue"  )
				colors = m.to_rgba(Rgc)
				#ax.errorbar(Xs ,  Ys ,  fmt='o',  c='black', marker='o', capsize=1.0, elinewidth=0.5, lw=0.5, markersize=1.5,label="")
				ax.scatter(Xs ,  Ys ,    c= colors , marker='o', s=10,  label="")
				if 1:  # colorbar
					from mpl_toolkits.axes_grid1 import make_axes_locatable
					divider = make_axes_locatable(ax)
					cax1 = divider.append_axes("right", size="3%", pad=0.05)
					cb = mpl.colorbar.ColorbarBase(cax1, norm=normV, cmap=cmap)

					# cb = fig.colorbar(sc, cax=cax1)
					# cax1.set_ylabel(r"$V_{\rm LSR}$ ($\rm km\ s^{-1}$)")

					cax1.set_ylabel(r"$\mathrm R_{\rm gal}$ (kpc)")
				#for i in range( len(Xs)):
					#ax.errorbar(Xs[i] ,  Ys[i] , yerr= YsError[i] , fmt='o',  c=  colors[i]  , marker='o', capsize=1.0, elinewidth=0.5, lw=0.5, markersize=1.5,label="")

				logX = np.log10( Xs )
				logY = np.log10(  Ys  )

				print np.polyfit( logX, logY,1)
				#logYError= np.log10(np.e)/Ys*YsError


				if 0:
					para,paraError= self.fitLarson(logX,logY,   None )

					print para," (a, b) of ax**b "
					drawXs=np.linspace(  min(Xs)-0.3, max(Xs)+1,10   )

					modelValue= larsonRelation( np.log10( drawXs),*para )
					ax.plot( drawXs,10**modelValue,color= "red",lw= 0.5 ,label=r"$\sigma \propto \mathit L^{{ {:.2f}\pm{:.2f}}}$".format(para[0], paraError[0] )  )

				else: #direct fit

					para,paraError= self.fitLarsonPower( Xs,Ys,None)

					print para," (a, b) of ax**b "
					drawXs=np.linspace(  min(Xs)-0.3, max(Xs)+1,10   )

					modelValue= larsonRelationPower(  drawXs ,*para )
					ax.plot( drawXs, modelValue,color= "black",lw= 0.7 ,label=r"$\sigma_\nu \propto \mathit L^{{ {:.2f}\pm{:.2f}}}$".format(para[1], paraError[1] )  )



				print


				#ax.scatter(  Xs ,  modelValue  ,  s= 8 ,  color="red"  )
				ax.set_xlim(2.5,70)

				ax.set_yticks([2,4,8,16])

				ax.set_xticks([5,10,20,40])
				ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
				ax.get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())

				ax.legend(loc=2 ,handlelength= 0.8 )

				ax.set_ylabel(  r"Equivalent line-width (km s$^{-1}$)" )
				ax.set_xlabel("Size (pc)")

				plt.savefig("velSize{}.pdf".format(self.regionName), bbox_inches='tight')
				plt.savefig("velSize{}.png".format(self.regionName), bbox_inches='tight', dpi=300)


			else:

				Xs =  sizeInPC[edgeCloud]
				Ys =   mass[edgeCloud]
				YsError= massError[edgeCloud]
				#ax.scatter(  Xs ,  Ys  ,  s= 8 ,  color="blue"  )

				#ax.errorbar(Xs ,  Ys , yerr= YsError , fmt='o',  c='b', marker='o', capsize=1.0, elinewidth=0.5, lw=0.5, markersize=1.5,label="")
				colors= m.to_rgba(Rgc)
				for i in range( len(Xs)):
					ax.errorbar(Xs[i] ,  Ys[i] , yerr= YsError[i] , fmt='o',  c=  colors[i]  , marker='o', capsize=1.0, elinewidth=0.5, lw=0.5, markersize=1.5,label="")

				#ax.errorbar(Xs ,  Ys , yerr= YsError , fmt='o',  c=  [0.        , 0.47254902, 1.        , 1.        ]  , marker='o', capsize=1.0, elinewidth=0.5, lw=0.5, markersize=1.5,label="")
				if 1:  # colorbar
					from mpl_toolkits.axes_grid1 import make_axes_locatable
					divider = make_axes_locatable(ax)
					cax1 = divider.append_axes("right", size="3%", pad=0.05)
					cb = mpl.colorbar.ColorbarBase(cax1, norm=normV, cmap=cmap)

					# cb = fig.colorbar(sc, cax=cax1)
					# cax1.set_ylabel(r"$V_{\rm LSR}$ ($\rm km\ s^{-1}$)")

					cax1.set_ylabel(r"$\mathrm R_{\rm gal}$ (kpc)")


				if 0:
					logX = np.log10( Xs )
					logY = np.log10(  Ys  )

					print np.polyfit( logX, logY,1)




					logYError= np.log10(np.e)/Ys*YsError

					para,paraError= self.fitLarson(logX,logY,  logYError )
					print para," (a, b) of ax**b "
					drawXs=np.linspace(  min(Xs)-0.3, max(Xs)+1,10   )

					modelValue= larsonRelation( np.log10( drawXs),*para )

					ax.plot( drawXs,10**modelValue,color= "red",lw= 0.5 ,label=r"$\mathit M\propto \mathit L^{{ {:.2f}\pm{:.2f}}}$".format(para[0], paraError[0] )  )
				else:
					para,paraError=  self.fitLarsonPower( Xs,Ys,YsError)
					print para," (a, b) of ax**b "
					drawXs=np.linspace(  min(Xs)-0.3, max(Xs)+5,50   )

					modelValue= larsonRelationPower(  drawXs ,*para )
					ax.plot( drawXs, modelValue,color= "black",lw= 0.7 ,label=r"$\mathit M\propto \mathit L^{{ {:.2f}\pm{:.2f}}}$".format(para[1], paraError[1] )  )



				#ax.scatter(  Xs ,  modelValue  ,  s= 8 ,  color="red"  )
				ax.set_xlim(2.5,70)

				ax.set_xticks([5,10,20,40])
				ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
				ax.legend(loc=2 ,handlelength= 0.8 )

				ax.set_ylabel(  r"Mass ($M_{\odot}$)" )
				ax.set_xlabel("Size (pc)")


				plt.savefig( "massSize{}.pdf".format(self.regionName), bbox_inches='tight')
				plt.savefig( "massSize{}.png".format(self.regionName), bbox_inches='tight',dpi=300)




	def fitLarsonPower(self,X,Y,Yerror):
		if Yerror is None:



			params, paramas_covariance = optimize.curve_fit(larsonRelationPower, X, Y,  p0=[ 1, 0.5])

			return params, np.sqrt(np.diag(paramas_covariance))
		else:


			params, paramas_covariance =  optimize.curve_fit(larsonRelationPower , X, Y, sigma=Yerror,   absolute_sigma=True, p0=[1, 0.5 ])


			return params,  np.sqrt(np.diag(paramas_covariance))
	def fitLarson(self,X,Y,Yerror):
		"""

		:param X:
		:param Y:
		:param Yerror:
		:return:
		"""

		if Yerror is None:



			params, paramas_covariance = optimize.curve_fit(larsonRelation, X, Y,  p0=[ 1, 0.5])

			return params, np.sqrt(np.diag(paramas_covariance))
		else:


			params, paramas_covariance =  optimize.curve_fit(larsonRelation , X, Y, sigma=Yerror,   absolute_sigma=True, p0=[1, 0.5 ])


			return params,  np.sqrt(np.diag(paramas_covariance))


	def COExtinctionLaw(self):
		"""
		examine the relation between
		:return:
		"""
		searchTablePath="/home/qzyan/WORK/projects/maddalena/dendroDisPath/disQ2/tmpFiles"


		deltaCOList = [] # CO different

		deltaExtinction = []



		for eachR in self.goodTB:

			cloudName = self.getCloudNameByLB( eachR["l"], eachR["b"] )

			figurePath = "/home/qzyan/WORK/projects/maddalena/dendroDisPath/{}/goodFigures/".format(self.regionName)
			files = glob.glob(figurePath + cloudName + "*")

			figureFile = files[0]
			if "AV" in figureFile:
				extictionStr = "AV"

			else:
				extictionStr = "AG"



			ID=self.getCloudIDByRow(eachR)

			tableOnFile= os.path.join( searchTablePath, "Cloud{}OnCloud{}.fit".format(ID,  extictionStr )  )
			tableOffFile= os.path.join( searchTablePath, "Cloud{}OffCloud{}.fit".format(ID,  extictionStr )  )

			TableOn=Table.read( tableOnFile )
			TableOff=Table.read( tableOffFile )


			dCO= np.mean(TableOn["coint"]) -   np.mean(TableOff["coint"])
			dAG = eachR["AG2"] -eachR["AG1"]

			deltaCOList.append( dCO )
			deltaExtinction.append( dAG )




		#draw


		fig=plt.figure(figsize=(8,6))
		rc('text', usetex=True )
		rc('font', **{'family': 'sans-serif',  'size'   : 16,  'serif': ['Helvetica'] })
		mpl.rcParams['text.latex.preamble'] = [
			r'\usepackage{tgheros}',  # helvetica font
			r'\usepackage{sansmath}',  # math-font matching  helvetica
			r'\sansmath'  # actually tell tex to use it!
			r'\usepackage{siunitx}',  # micro symbols
			r'\sisetup{detect-all}',  # force siunitx to use the fonts
		]



		ax=fig.add_subplot(1,1,1)

		ax.scatter( deltaCOList,   deltaExtinction,s=10 )

		plt.savefig("coAG{}.pdf".format(self.regionName), bbox_inches='tight')
		plt.savefig("coAG{}.png".format(self.regionName), bbox_inches='tight', dpi=300)



	def testDraw(self):
		"""
		this function is used to some test
		:return:
		"""




		fig=plt.figure(figsize=(8,6))
		rc('text', usetex=True )
		rc('font', **{'family': 'sans-serif',  'size'   : 16,  'serif': ['Helvetica'] })
		mpl.rcParams['text.latex.preamble'] = [
			r'\usepackage{tgheros}',  # helvetica font
			r'\usepackage{sansmath}',  # math-font matching  helvetica
			r'\sansmath'  # actually tell tex to use it!
			r'\usepackage{siunitx}',  # micro symbols
			r'\sisetup{detect-all}',  # force siunitx to use the fonts
		]



		ax=fig.add_subplot(1,1,1)

		Y=self.goodTB["distance"] #*np.cos( np.deg2rad(b) )
		l=self.goodTB["l"]
		b=self.goodTB["b"]

		if 0: #get the cloud with the largest scale height
			scaleArray =  abs( np.sin( np.deg2rad(b) ) * Y )


			maxScale= np.max( abs(scaleArray))

			print self.goodTB[ scaleArray ==maxScale ]

		if 1:
			pass
		#print

		X= self.goodTB["v_cen"] #*np.cos( np.deg2rad(b) )
		vStd= self.goodTB["v_rms"] #*np.cos( np.deg2rad(b) )



		#import matplotlib as mpl
		cmap=plt.cm.jet
		#norm = matplotlib.colors.BoundaryNorm(np.arange(0,30,0.1), cmap.N)

		normV = mpl.colors.Normalize(vmin= min(l), vmax= max(l) )

		m = plt.cm.ScalarMappable(norm= normV , cmap=cmap)




		ax.scatter(  X, Y    ,s=10  )
		#ax.scatter( b, 1000*np.sin( np.deg2rad(b) ) ,s=10 )






		plt.savefig("testDraw{}.pdf".format(self.regionName), bbox_inches='tight')
		plt.savefig("testDraw{}.png".format(self.regionName), bbox_inches='tight', dpi=300)



	def testQ2PerseusArm(self):
		"""
		Examine the stream motion of the Perseus motion in the second quadrant
		:return:
		"""
		mergedTB=self.goodTB



		perSeus = mergedTB[  mergedTB["vlsr"]<-29 ]
		print len( perSeus  )

		dA5List=[]

		d2019List = []

		for eachPC in perSeus:

			gaiaDis= eachPC["distance"]
			l,b,v= eachPC["l"] ,  eachPC["b"] ,  eachPC["vlsr"]

			v=v+ 16

			reid2014 = self.ReidA5(l,b,v)
			dA5List.append( reid2014[0][0]*1000 -  gaiaDis   )
			reid2019= self.Reid2019(l,b,v)
			d2019List.append( reid2019[0] *1000 -  gaiaDis   )


		print np.mean(dA5List), np.mean(d2019List),"distance different in pc of A5 and 2019 "



	def drawFaceOnSingleTB(self,drawTB,	drawLrange, saveTag="" ):
		"""
		draw faceon view of
		:return:
		"""


		xSun = 0.
		ySun = self.R0
		sunRadius = 0.02
		pSun = np.array([xSun, ySun])

		maxR=3
		fig = plt.figure(figsize=(10, 8))
		# fig, axs = plt.subplots(nrows=1, ncols=2,  figsize=(12,6),sharex=True)
		rc('text', usetex=True)
		rc('font', **{'family': 'sans-serif', 'size': 16, 'serif': ['Helvetica']})

		mpl.rcParams['text.latex.preamble'] = [
			r'\usepackage{tgheros}',  # helvetica font
			r'\usepackage{sansmath}',  # math-font matching  helvetica
			r'\sansmath'  # actually tell tex to use it!
			r'\usepackage{siunitx}',  # micro symbols
			r'\sisetup{detect-all}',  # force siunitx to use the fonts
		]

		# custom coors with velocity


		ax0 = fig.add_subplot(1, 1, 1)
		ax0.set_aspect("equal")

		ax0.set_xlabel(r"X (kpc)")
		ax0.set_ylabel(r"Y (kpc)")


		print "Draw clouds of  ", self.regionName, "Total number of clouds:", len(drawTB)
		self.drawFaceOnByTB(ax0, drawTB, None , marker="o")


		print "Draw lines and Markers"



		drawCircleDistances = np.arange(0.5, maxR + 0.5, 0.5)  # [   0.5, 1.5, 2.5, 3.5]
		for radiusCircle in drawCircleDistances:
			# for radiusCircle in [0.25,  0.5,1,1.5,2,2.5,3,3.5]:
			"""
				"""
			# circle=plt.Circle((xSun,ySun),radius=3 )
			if drawLrange is None:
				drawA = np.linspace(-np.pi / 2, np.pi / 2, 200)

			else:
				drawA = np.linspace(np.deg2rad(min(drawLrange) - 90), np.deg2rad(max(drawLrange) - 90), 200)

			x = radiusCircle * np.cos(drawA)
			y = radiusCircle * np.sin(drawA)
			ax0.plot(x + xSun, y + ySun, '--', color='black', lw=0.5, alpha=0.5)

			# if radiusCircle!=2:
			testAngle = 90. / 180 * np.pi  # angle where the distances were put
			ax0.text(radiusCircle * np.cos(testAngle) + xSun + 0.05, radiusCircle * np.sin(testAngle) + ySun + 0.05,
					 "{} kpc".format(radiusCircle), fontsize=10, ha='center', va="center")

		maxDrawR = max(drawCircleDistances) + 0.5

		ax0.scatter(xSun, ySun, s=5, facecolors="black")

		ax0.scatter(xSun, ySun, s=20, facecolors='none', edgecolors='black', lw=0.3)

		ax0.text(xSun, ySun - 0.15, "the Sun", fontsize=11, ha='center')

		# for drawL in [ 180, 190,200,210,220,230,240,250,260,270]:
		# for drawL in [25,30,35,40,45,50  ]:

		if drawLrange is not None:
			lArray = np.arange(min(drawLrange), max(drawLrange) + 10, 10)  #
		else:
			lArray = np.arange(0, 190, 10)
		for drawL in lArray:
			drawAngle = np.radians(drawL - 90)

			unitvector = np.array([np.cos(drawAngle), np.sin(drawAngle)])
			drawp1 = pSun - unitvector * maxDrawR

			drawp1_end = pSun - unitvector * sunRadius

			# ax0.plot( [drawp1_end[0],drawp1[0]],   [drawp1_end[1],drawp1[1]] ,'--',color='black',lw=0.5,alpha=0.5 )

			drawp2_start = pSun + unitvector * sunRadius

			drawp2 = pSun + unitvector * maxDrawR
			ax0.plot([drawp2_start[0], drawp2[0]], [drawp2_start[1], drawp2[1]], '--', color='black', lw=0.5, alpha=0.5)

			unitvector = np.array([np.cos(drawAngle + 0.02), np.sin(drawAngle + 0.025)])

			refPosition = pSun + unitvector * (maxDrawR + 0.07)

			thetaShift = np.radians(drawL)
			# shiftVector=  np.array( [ np.cos(thetaShift), np.sin(thetaShift) ]  )
			# shiftVector=shiftVector *0.5
			drawp2Text = refPosition  # +shiftVector
			ax0.text(drawp2Text[0], drawp2Text[1], r"${}^\circ$".format(drawL), fontsize=10, rotation=drawL - 180,
					 va="center", ha="left", rotation_mode='anchor')

		##############################################################################################################

		# sc=ax0.scatter( Xs,Ys  ,  c =Vs, cmap=cmap,norm=normV, s=13 ,facecolors='none',   lw=0.5,marker="o")

		ax0.set_xlim([-0.3-maxR, maxR + 1])

		# ax0.set_xlim([-2,  2 ])
		#if showXrange is not None:
			#ax0.set_xlim([-3,2])
		#if showYrange is not None:
			#ax0.set_ylim(showYrange)

		plt.subplots_adjust(wspace=0.3)
		# plt.show()
		plt.tight_layout(pad=0)


		plt.savefig("faceOn{}{}.png".format(self.regionName, saveTag), bbox_inches='tight', dpi=600)
		plt.savefig("faceOn{}{}.pdf".format(self.regionName, saveTag), bbox_inches='tight')

	def drawFaceOnOfAllMCQ2(self):
		"""
		draw faceon view of
		:return:
		"""

		drawTB= "/home/qzyan/WORK/myDownloads/disQ2/goodCloudQ2.fit"
		drawTB = Table.read( drawTB )
		#drawTBReid = drawTB[0:50]

		drawTBReid = drawTB.copy() #drawTB[0:50]


		add3DInfoMWISP = drawTB.copy()

		#disArray, errorLow,errorUp,  lArray,bArray = self.getDLB_Reid2014(drawTB )

		#self.add3DInfoReid2014(drawTBReid )
		self.add3DInfoReid2019(drawTBReid )
		drawTBReid = drawTBReid[drawTBReid[self.colXsun] >0 ]
		print drawTBReid

		self.add3DInfoMWISP(add3DInfoMWISP )



		#drawFaceOnView

		# plot Xs, Ys
		fig = plt.figure(figsize=(12, 6))

		import matplotlib as mpl

		# fig, axs = plt.subplots(nrows=1, ncols=2,  figsize=(12,6),sharex=True)
		rc('text', usetex=True)
		rc('font', **{'family': 'sans-serif', 'size': 14, 'serif': ['Helvetica']})

		mpl.rcParams['text.latex.preamble'] = [
			r'\usepackage{tgheros}',  # helvetica font
			r'\usepackage{sansmath}',  # math-font matching  helvetica
			r'\sansmath'  # actually tell tex to use it!
			r'\usepackage{siunitx}',  # micro symbols
			r'\sisetup{detect-all}',  # force siunitx to use the fonts
		]

		# custom coors with velocity
		axReid = fig.add_subplot(1, 2, 1)
		axMWISP = fig.add_subplot(1, 2, 2) #draw faceon view use the linear relationship of distances

		self.drawFaceOnByTB_useReid2014(axReid, drawTBReid  )


		#

		#####



		self.drawFaceOnByTB_useReid2014(axMWISP, add3DInfoMWISP  )

		axReid.set_xlabel(r"X (kpc)")
		axReid.set_ylabel(r"Y (kpc)")

		axMWISP.set_xlabel(r"X (kpc)")
		axMWISP.set_ylabel(r"Y (kpc)")


		###
		xSun=0.
		ySun=self.R0
		sunRadius=0.02
		###################
		axMWISP.scatter(xSun,ySun,s=5,facecolors="black")
		axMWISP.scatter(xSun,ySun,s=20,facecolors='none', edgecolors='black',lw=0.3)
		axMWISP.text( xSun+0.12  ,ySun-0.15,"the Sun",fontsize=9 ,ha='center'  )

		axReid.scatter(xSun,ySun,s=5,facecolors="black")
		axReid.scatter(xSun,ySun,s=20,facecolors='none', edgecolors='black',lw=0.3)
		axReid.text( xSun+0.12  ,ySun-0.15,"the Sun",fontsize=9 ,ha='center'  )

		##############

		at = AnchoredText('Kinematic distances with Reid et al. (2019)', loc=4, frameon=False)
		axReid.add_artist(at)

		at = AnchoredText(r'Distances with the $V_{\rm LSR}$-distance relationship', loc=4, frameon=False)
		axMWISP.add_artist(at)

		axMWISP.set_aspect('equal' )
		axReid.set_aspect('equal' )

		axReid.set_xlim(-0.2, 5 )
		axMWISP.set_xlim(-0.2,5 )

		axReid.set_ylim( 7.9  , 13 )
		axMWISP.set_ylim( 7.9 , 13  )



		saveTag="reidAndMWIPS"
		plt.savefig( "faceOn{}_{}.png".format(self.regionName, saveTag),bbox_inches='tight',dpi=600)
		plt.savefig( "faceOn{}_{}.pdf".format(self.regionName, saveTag ), bbox_inches='tight'  )


	def ZZZ(self):
			pass



def zzzzzzzz(self):
	pass




#cloudTB=Table.read("/home/qzyan/WORK/myDownloads/testScimes/mosaicV1NewTB.fit")

if 0:

	dbscanTBName=  "/media/qzyan/maclinux/projects/disAntiGC/tmpPath/cleanTBAGC.fit"

	disAGCGoodTB="/home/qzyan/WORK/projects/maddalena/dendroDisPath/disAGC/disAGCgoodDisTB.fit"
	drawDis=disDraw(goodTBName= disAGCGoodTB ,dbscanTBName= dbscanTBName, regionName="disAGC")

	drawDis.edgeCloudNameList=["G149.4+03.1"] #not edge cloud concerned
	print drawDis.ReidA5( 16.103 ,0.896 ,  143.16 )
	print drawDis.Reid2019( 16.103 ,0.896 ,  143.16 )
	print drawDis.Reid2019( 16.995 , -0.647 ,  124.93  )

if 0:#disAGC

	dbscanTBName=  "/media/qzyan/maclinux/projects/disAntiGC/tmpPath/cleanTBAGC.fit"

	disAGCGoodTB="/home/qzyan/WORK/projects/maddalena/dendroDisPath/disAGC/disAGCgoodDisTB.fit"
	drawDis=disDraw(goodTBName= disAGCGoodTB ,dbscanTBName= dbscanTBName, regionName="disAGC")

	drawDis.edgeCloudNameList=["G149.4+03.1"] #not edge cloud concerned
	drawDis.drawVelLongitudeDiagram()
	drawDis.drawVelLongitudeDiagram(useReid2019=True)

if 0: #test Reid 2016
	dbscanTBName=  "/media/qzyan/maclinux/projects/disAntiGC/tmpPath/cleanTBAGC.fit"

	disAGCGoodTB="/home/qzyan/WORK/projects/maddalena/dendroDisPath/disAGC/disAGCgoodDisTB.fit"
	drawDis=disDraw(goodTBName= disAGCGoodTB ,dbscanTBName= dbscanTBName, regionName="disAGC")

	drawDis.edgeCloudNameList=["G149.4+03.1"] #not edge cloud concerned
	print drawDis.Reid2016(125,0, -20 )
	print drawDis.Reid2016(150, 5, -30  )
	print drawDis.Reid2016(105, -5, -10 )


if 0:#disG220

	dbscanTBName=  "/media/qzyan/maclinux/projects/disG220/tmpPath/cropCOFITS_disG220dbscanS2P4Con1_Clean.fit"

	disG220GoodTB="/home/qzyan/WORK/projects/maddalena/dendroDisPath/disG220/disG220goodDisTB.fit"
	drawDis=disDraw(goodTBName= disG220GoodTB ,dbscanTBName= dbscanTBName, regionName="disG220")

	drawDis.edgeCloudNameList=["G149.4+03.1"] #not edge cloud concerned

	drawDis.drawMassSize(drawSizeVel= True  )
	drawDis.drawMassSize(drawSizeVel= False , BFFcor= True )




if 0:#disAGC
	dbscanTBName= "/home/qzyan/WORK/diskMWISP/MWISPData/G105150Tmp/mergedCubedbscanS2P4Con1_Clean.fit"

	disQ2GoodTB="/home/qzyan/WORK/projects/maddalena/dendroDisPath/disQ2/disQ2goodDisTB.fit"
	drawDis=disDraw(goodTBName= disQ2GoodTB ,dbscanTBName= dbscanTBName, regionName="disQ2")

	drawDis.edgeCloudNameList=["G149.4+03.1","G106.5+04.0","G145.8+03.4","G150.1-01.4"]




	pass
	#extraTB1 = Table.read("/home/qzyan/WORK/projects/maddalena/dendroDisPath/disG220/disG220goodDisTB.fit")
	#extraTB2 = Table.read("/home/qzyan/WORK/projects/maddalena/dendroDisPath/disQ2/disQ2goodDisTB.fit")
	#extraTB3 = Table.read("/home/qzyan/WORK/projects/maddalena/dendroDisPath/G2650/G2650goodDisTB.fit")
	#extTB=[ extraTB1,extraTB2    ]
	#drawDis.drawFaceOnSimple( drawLrange =[  100,230],calChains=1,drawMaser=False ,drawMaserLrange=[  100,230],drawMaserBrange=[-5.25,5.25], drawMaserVrange=[-20, 70], extraMWISPTB  = None, onlyFaceon=True )
	#drawDis.draw3Dmap()
	#drawDis.drawMassSize(drawSizeVel= True  )
	#drawDis.drawMassSize(drawSizeVel= False , BFFcor= True )

	#drawDis.drawFaceOnSimple( drawLrange =[  80,180],calChains=1,drawMaser=False ,drawMaserLrange=[  100,230],drawMaserBrange=[-5.25,5.25], drawMaserVrange=[-20, 70], extraMWISPTB  = None, onlyFaceon=True )


if 0: #draw distance
	pass
	drawDis.printCatlog( )
if 0:



	#drawDis.compareDisWithA5(withReid2019=False)
	#drawDis.compareDisWithA5(withReid2019=False)

	drawDis.compareDisWithA5(withReid2019=True)
	drawDis.compareDisWithA5(withReid2019=True)
if 0: #draw distance
 
	drawDis.testQ2PerseusArm( )



if 0:
	drawDis.testDraw()

if 0:
	print "local Clouds"
	#print  "105 degrees", drawDis.ReidA5(105 ,0 ,0)



	for disCut in [-79,-6,5,30,70,139]:

		print "25 degress {:.2f} km/s".format(disCut),drawDis.ReidA5(25 ,0 , disCut)
		print "50 degress {:.2f} km/s".format(disCut) , drawDis.ReidA5(50, 0, disCut)

#print drawDis.ReidA5(46.284, -1.660, 7.1)
	########################


if 0:
	pass
	#drawDis.COExtinctionLaw()

	#drawDis.drawSideView()
	#drawDis.draw3Dmap()
	#drawDis.drawFaceOnSimple( drawMaser=False , colorVrange=[-45,5 ] ,drawLrange=[60, 180],maxR= 2.5, calChains =1 )

	#drawDis.drawFaceOnSimple( drawMaser=True , colorVrange=[-50,5] ,drawLrange=[60,180],maxR= 2.5, calChains =1 )
	#drawDis.drawArmView( drawMaser=True , colorVrange=[-50,5] ,drawLrange=[80,180],maxR= 8, calChains =1 )


	drawDis.drawFaceOnOfAllMCQ2()

	#print drawDis.getMassList()
	#drawDis.drawMassSize(drawSizeVel=True)
	#drawDis.drawMassSize(drawSizeVel=False )

#print "The mass of the largest molecular cloud is ",drawDis.calmassByXfactor(80336760.0*drawDis.dv, 1000)

if 0:
	scimesTBName= "/home/qzyan/WORK/myDownloads/MWISPcloud/G2650CloudForDisCat.fit"

	G2650GoodTB="/home/qzyan/WORK/projects/maddalena/dendroDisPath/G2650/G2650goodDisTB.fit"
	drawDis=disDraw( G2650GoodTB ,"G2650")
	drawDis.scimesTB=Table.read( scimesTBName  )


if 0:
	print drawDis.ReidA5(46.284,-1.660, 7.1)

	sys.exit()
if 0:
	print

	sys.exit()
if 0:

	drawDis.compareDisWithA5(withReid2019=False)
	drawDis.compareDisWithA5(withReid2019=False)

	drawDis.compareDisWithA5(withReid2019=True)
	drawDis.compareDisWithA5(withReid2019=True)



#drawDis.test()


if 0:#test Reid 2019

	print drawDis.ReidA5(35.19,-0.74,30)
	print drawDis.Reid2019(35.19,-0.74,30)

	#print drawDis.Reid2019(26 ,0, 16.3 )

	#print drawDis.Reid2019(50 ,0, 16.3 )