
#some table of gaia data base
#this version of Gaia,perform systematic correction

import psycopg2

import MySQLdb
from astropy.table import   Table,vstack
import numpy as np


class GAIATB:

	
	#name='public."G2650AG"' #table name
	name='public."AGall"' #table name

	GAIA_source_id ="source_id"


	GAIA_parallax ="parallax"

	GAIA_parallax_err ="parallax_error"
	
	GAIA_l ="l"
	GAIA_b ="b"

	GAIA_distance ="distance"
	GAIA_distanceError ="distance_err"


	GAIA_a_g_val ="a_g_val"

	GAIA_a_g_percentile_lower ="a_g_percentile_lower"
	GAIA_a_g_percentile_upper ="a_g_percentile_upper"

	agError ="agError" # relative ag error


	relative_error ="relative_error" ##relative parallax error
	phot_g_mean_mag ="phot_g_mean_mag" ##relative parallax error


	#colnames=[GAIA_source_id, GAIA_parallax,GAIA_parallax_err, GAIA_l,  GAIA_b, GAIA_a_g_val,  GAIA_a_g_percentile_lower, GAIA_a_g_percentile_upper,relative_error,agError   ]
 
	colnames=[GAIA_source_id, GAIA_parallax,GAIA_parallax_err, GAIA_l,  GAIA_b, GAIA_a_g_val,  GAIA_a_g_percentile_lower, GAIA_a_g_percentile_upper, phot_g_mean_mag  ]

	dataTypes=[ float,float,float,float, float, float,float,float ,float ]  #all float
	
	GAIA_distance ="distance"
	GAIA_distanceError ="distance_err"




	def getDB(self):
		#read dbInfo.txt
 		
		#fileA = open('dbInfo.txt', 'r')
		#a = fileA.readlines()
		#db = MySQLdb.connect(str.strip(a[0]),str.strip(a[1]),str.strip(a[2]),str.strip(a[3]) )
		#db = MySQLdb.connect("127.0.0.1","root","shao1234","gaia")
		#return db
		#db = MySQLdb.connect(str.strip(a[0]),str.strip(a[1]),str.strip(a[2]),str.strip(a[3]) )
		#return db

		#why not using mysql???
		connection = psycopg2.connect(user = "postgres", password = "100425", host = "127.0.0.1", port = "5432", database = "gaia")


		#db = MySQLdb.connect(str.strip(a[0]),str.strip(a[1]),str.strip(a[2]),str.strip(a[3]) )
		return connection


	def getByLBRange(self,Lrange,Brange,lowerPara=0.2,paraError=0.2,upperPara=None,calDis=False):
		
		"""
		dedicated to find other regions that are overlaping with the current source

		# lowerPara is the samallest parallax, corresponding to the farthest distance
		
		# upperPara, is largest parallax, corresponding to the nearest distance
		
		
		
		By default, paraError less than 20%, and disances less than 5 kpc
		
		"""

		sL= min(Lrange)
		
		eL= max(Lrange)
#
		sB=  min(Brange)
		eB=  max(Brange)


		db = self.getDB() #MySQLdb.connect("localhost","root","shao1234","gaia" )
		# prepare a cursor object using cursor() method
		cursor = db.cursor()
		
		
		
 		
		if upperPara is None:

			sqlCommand="select * from {} where  l > {} and l < {} and b > {} and b < {} and parallax>{} and parallax_error<parallax*{};".format(self.name,sL,eL,sB,eB,lowerPara,paraError)
		
		else:
			sqlCommand="select * from {} where  l > {} and l < {} and b > {} and b < {} and parallax>{} and parallax<{}  and parallax_error<parallax*{};".format(self.name,sL,eL,sB,eB,lowerPara,upperPara,paraError)

 



		# execute SQL query using execute() method.
		cursor.execute(sqlCommand)
		# Fetch a single row using fetchone() method.
 
		data = np.array(cursor.fetchall() )
		db.close()
 
		if len(data)==0  :
			#db.commit()
 
 			print sqlCommand
			print "No Gaia stars found."
			return None
 
		#t=  self.converSqlToTB(data)
		#construct a TB with data
		#db.commit()
		
		queryTB=Table(rows= data , names=self.colnames,dtype=('f8','f8','f8','f8','f8','f8','f8','f8','f8' ))

		queryTB=self.correctParallax(queryTB)
		queryTB=self.correctParallaxError(queryTB)

		#becaue after correction, error tends to be larger
		#need to reject stars with large errors, should we do this first?
		queryTB=  queryTB[ queryTB[self.GAIA_parallax_err ] <0.2*queryTB[self.GAIA_parallax ]   ]
		#since there is an extradis, we do not care parallax


		#add relative_error,agError
		queryTB[self.relative_error] = queryTB[self.GAIA_parallax_err]/ queryTB[ self.GAIA_parallax ]

		queryTB[self.agError] = (  queryTB[self.GAIA_a_g_percentile_upper] - queryTB[self.GAIA_a_g_percentile_lower]  )/2/queryTB[ self.GAIA_a_g_val ]



		if calDis:
			print "Calculating distances...."
			return self.addDisToTB(queryTB )

		return queryTB
		


	def correctParallaxError(self,TB):
		"""
		Could be empty
		:param TB:
		:return:
		"""
		print "Correcting systematic parallax errors..."

		TB.add_index(self.phot_g_mean_mag)

		part1= TB.loc[ self.phot_g_mean_mag,:11  ]

		part1[ self.GAIA_parallax_err ]=   part1[self.GAIA_parallax_err ]*1.2


		part2= TB.loc[ self.phot_g_mean_mag,11: 15  ]
		part2[self.GAIA_parallax_err ]=  part2[self.GAIA_parallax_err]*( 0.22*part2[self.phot_g_mean_mag]-1.22    )

		part3=  TB.loc[ self.phot_g_mean_mag,15:   ]
		part3[ self.GAIA_parallax_err ]=  part3[self.GAIA_parallax_err]*(  1.08+np.exp(15-    part3[self.phot_g_mean_mag] )   )

		newTB=vstack([part1,part2,part3])
		return newTB

	def correctParallax(self,TB):
		"""

		:param TB:
		:return:
		"""

		print "Correcting systematic Parallax..."


		TB.add_index(self.phot_g_mean_mag)

		part1= TB.loc[ self.phot_g_mean_mag,:14  ]

		part1[ self.GAIA_parallax ]=  part1[self.GAIA_parallax ]+0.05

		part2= TB.loc[ self.phot_g_mean_mag,14: 16.5  ]
		part2[self.GAIA_parallax ]=  part2[self.GAIA_parallax]+  (0.1676 -0.0084*part2[self.phot_g_mean_mag] )

		part3=  TB.loc[ self.phot_g_mean_mag,16.5:   ]
		part3[ self.GAIA_parallax ]=  part3[self.GAIA_parallax ]+0.029

		newTB=vstack([part1,part2,part3])
		return newTB

	def addDisToTB(self,gaiaTB ):
		
		"""
		
		calculate distance, and add them to gaiaTB,
		
		usually only do once
		"""
		#copy file
		
		newTB=gaiaTB.copy()
 




		parallaxErrCol="parallax_err"

		if "parallax_error" in newTB.colnames: # for new gaia TB
			parallaxErrCol="parallax_error"


		#newTB[self.GAIA_distance]=
		#newTB.add_columns( [disCol,  disErrCol  ])



		paraCol= newTB["parallax"]
		paraErrorCol=newTB[   parallaxErrCol ]


		totalN=len( paraCol )



		samplePara=np.random.normal( paraCol,paraErrorCol,size=(20000,totalN) )
		sampleDis= 1./samplePara*1000

		disArray= np.mean( sampleDis , axis=0 )
		disArray= np.round( disArray , 2 )

		disStdArray= np.std( sampleDis,  ddof=1, axis=0 )
		disStdArray= np.round( disStdArray,   2 )

		newTB[ self.GAIA_distance ] =  disArray
		newTB[ self.GAIA_distanceError ] = disStdArray

		return newTB


		for eachRow in newTB: 
			para=eachRow["parallax"]
			paraError= eachRow[parallaxErrCol] #+0.1 #   #
			
			dA=1./np.random.normal(para,paraError,20000)*1000
			
			eachRow[self.GAIA_distance]=  round(np.mean(dA), 2) 
			
			eachRow[self.GAIA_distanceError]=   round(np.std(dA,ddof=1), 2) 
		




	def converSqlToTB(self,sqlTB):
		
		"""
		"""
		
		pass
		
	def filterByLB(self,TB):
		
		"""
		"""
		
		pass