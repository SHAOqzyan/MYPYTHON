
#deal with  dendrotree

from astropy.table import Table,vstack

class dendroTree:
	
	
	treeDict=None
	
	
	goodList=None # sources that have distances deteremined
	
	#usually star from turunk, usually, we only care about trunk
	
	dendroTB=True
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
		
		
		with open(self.treeFile) as f:
		    lines = f.readlines()
			
		for eachLine in lines:
			indexC,indexP =map(int,eachLine.split())			
			
			treeDict[indexC]=indexP #search Parent
			#a node may have many two children, but only one parent
		self.treeDict=treeDict
		
	def isTrunk(self,ID):
		
		return self.treeDict[ID]==-1


	############################################################## 

	def isLeaf(self,ID):
		
 
		for k, v in self.treeDict.items():
			if v==ID:
				return False
				
		return True
 				
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
		
 
		
		children=[]
		
		if self.isLeaf(ID):
			return []
		
		for k, v in self.treeDict.items():
			if v==ID:
				children.append(k)
				
		return children
 				

 				
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