
#deal with  dendrotree

from astropy.table import Table,vstack
import numpy as np
import math
def weighted_avg_and_std(values, weights=None):
    """
    Return the weighted average and standard deviation.

    values, weights -- Numpy ndarrays with the same shape.
    """
    average = np.average(values, weights=weights)
    # Fast and numerically precise:
    variance = np.average((values-average)**2, weights=weights)
    return (average, math.sqrt(variance))


class dendroTree:
	
	
	treeDict=None
	
	
	goodList=None # sources that have distances deteremined
	
	#usually star from turunk, usually, we only care about trunk
	
	dendroTB=True

	childDict=None

	def __init__(self,treeFile,dendroCat=None):
		"""
		pass
		"""
		
		#treeFile is the only file needed 
		
		
		
		self.treeFile=treeFile
		self.dendroCat=dendroCat #the dendroCat		
		
		if self.dendroCat != None:
			
			self.dendroTB=Table.read(self.dendroCat)
		
		
		self.getTreeDict()





	def getTreeDict(self):
		treeDict={}

		childDict={}
		
		with open(self.treeFile) as f:
		    lines = f.readlines()
			
		for eachLine in lines:
			indexC,indexP =map(int,eachLine.split())			
			
			treeDict[indexC]=indexP #search Parent

			if indexP in childDict.keys():
				childDict[indexP ] .append(  indexC )
			else:
				childDict[indexP ]=[indexC]


			#a node may have many two children, but only one parent
		self.treeDict=treeDict
		self.childDict=childDict

	def isTrunk(self,ID):
		
		return self.treeDict[ID]==-1


	############################################################## 

	def isLeaf(self,ID):
		
		return  ID not in self.childDict.keys()

 				
	def getAllTrunk(self):
		"""
		print all turnkID
		"""
		turnkIDS=[]
		for k, v in self.treeDict.items():
			if v==-1:
				turnkIDS.append(k)
		turnkIDS.sort()
		return turnkIDS
				#return False



	def getChildren(self,ID):

		if ID not in self.childDict.keys():
			return None

				
		return self.childDict[ID]


	def getAllLeaves(self,ID):


		if ID not in self.childDict.keys():
			return [ID]

		returnList=[]

		self.getAllLeavesRecurse(ID,returnList)


		return returnList






	def getAllLeavesRecurse(self,ID,collectList):
		"""
		Return all leave IDs
		:param ID:
		:return:
		"""

		allChindren=self.getChildren(ID)

		if allChindren!=None:

			#print (self.dendroTB[allChindren[0]]['v_cen']-self.dendroTB[allChindren[1]]['v_cen'])/1000.

			self.getAllLeavesRecurse( allChindren[0] ,collectList )
			self.getAllLeavesRecurse(  allChindren[1] ,collectList )

		else:

			collectList.append(ID)

	def getMaxDV(self,ID):
		"""

		:param ID:
		:return:
		"""
		# get the maximum velocity range for branches

		allLeaf=self.getAllLeaves(ID)



		subTB=self.dendroTB[allLeaf]
		vrms=subTB[0]['v_rms']
		if len(subTB)==1:
			return vrms

		dv= max(subTB['v_cen'].data) -  min(subTB['v_cen'].data)
		#return max(dv/1000.,vrms)

		return  dv/1000.
	def getMaxDVAnLostFlux(self,ID):
		"""

		:param ID:
		:return:
		"""
		# get the maximum velocity range for branches

		allLeaf=self.getAllLeaves(ID)


		subTB=self.dendroTB[allLeaf]
		vrms=subTB[0]['v_rms']
		if len(subTB)==1:
			return vrms

		dv= max(subTB['v_cen'].data) -  min(subTB['v_cen'].data)
		#return max(dv/1000.,vrms)

		#mean,std=weighted_avg_and_std( subTB['v_cen'].data, weights= subTB['flux'].data )
		mean,std=weighted_avg_and_std( subTB['v_cen'].data  )

		#return   6*std/1000.,  np.sum( subTB["flux"]  ),len(subTB)
		return  dv/1000.,  np.sum( subTB["flux"]  ),len(subTB)

 				
	def getTrunkID(self,ID):
		
		"""
		Find the trunk of the ID
		"""
		
		sID=ID		
		
		while sID>=0:
			
			
			
			pID=self.treeDict[sID]
			
			if pID>=0:
				sID=pID
			else:
				break
		return sID
			
			#for k, v in self.treeDict.items():
				
		
	
		 
	def getExtendNodes(self,nodes):
		
		oldNodes=nodes[:]
		newNodes=nodes[:]
			
		for eachNode in oldNodes:
			
			children=self.getChildren(eachNode)
 
 
			if children==[]:
				continue

			else:
				
				for eachChild in children:
					
					if eachChild in newNodes:
						continue
					else:
						newNodes.append( eachChild)
				
 
			
		return newNodes
				
			
 				
 				
	def getAllRelated(self,ID):
		
		
		# get all the nodes that are in the same molecular clouds,
		
		trunkID=self.getTrunkID(ID)
		
 
 
		
		searchNodes=[ trunkID]
		
		while True:
			
			beforeSearch=searchNodes[:]
			
			searchNodes=self.getExtendNodes(searchNodes)
 
			if len(searchNodes)==len(beforeSearch):
				break
		return searchNodes
		#first find the trunck
		
		#for k, v in self.treeDict.items():
			#if v==ID:
				
				#childrenList.append(k)
			
				
	def setCatFile(self,catFile):
		
		self.dendroCat= catFile		
		
		self.dendroTB=Table.read(self.dendroCat)
		
	
			
	def setGoodList(self,goodList):
		
		self.goodList= goodList			

		#the length of CatFile should be the same as tree file


	def getBadTrunk(self):
		
		
		if self.goodList == None:
			print "No good source list provided!"
			return 


		goodTrunk=[]
		
		for eachGID in self.goodList:
			
			goodTID=self.getTrunkID(eachGID )
			
			goodTrunk.append(goodTID)

		allTrunkID=self.getAllTrunk()
		
		badTrunk=[]
		
		for eachTID in allTrunkID:
			
			if eachTID not in goodTrunk:
				badTrunk.append( eachTID)


		return badTrunk
		
	def getDendroRowByID(self,ID):
		
		"""
		"""
		
 
		if  self.dendroTB==None:
			
			print "Please set the dendroCatFile!"
			
			return
			
		IDList=list(  self.dendroTB["_idx"]  )

		TBID=None
		
		try :
			TBID=IDList.index( ID )
			
		except:
			pass
			
		#####
		
		if TBID==None:
			print "No ID found"
			return 
		
		
		
		return self.dendroTB[TBID]
		



	def filterListByVel(self,inputList,startV=None,endV=None):
		
		"""
		"""
		
		newList=[]
		
		if startV!=None and endV!=None:
			
			oldV=[ startV, endV]
			
			startV=min( oldV )
			end=max( oldV )

		
		
		for eachID in inputList:
			
			row=self.getDendroRowByID(eachID)

			if row==None:
				continue 
			print row['_idx']
			
			if startV!=None and row["v_cen"]/1000.< startV  :
				continue
				
			if endV!=None and row["v_cen"]/1000.> endV  :
				continue


			newList.append( eachID )
			
			
		return newList
		


	def ZZZ(self):
		pass




#treeFile="/home/qzyan/WORK/projects/maddalena/G130150Path/G130150Tree.txt"


#doTree=dendroTree(treeFile)

#print doTree.getAllRelated(1)