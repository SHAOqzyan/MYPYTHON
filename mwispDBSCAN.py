
#This script is used reconstruct the script of DBSCAN.
#To make is more usefule to use. and can be imported septeratrely

import numpy as np
from astropy import units as u
from spectral_cube import SpectralCube
from astropy.table import Table,vstack
#need to download myPYTHON.py
from myPYTHON import *
import os
from progressbar import *

import sys
from skimage.morphology import erosion, dilation
from scipy.ndimage import label, generate_binary_structure,binary_erosion,binary_dilation
from sklearn.cluster import DBSCAN
import requests

doFITS=myFITS()

#print doFITS.weighted_avg_and_std

class  MWISPDBSCAN(object):


    processPath="./"#  by default save the product in the current path

    rawCOFITS=None
    rmsFITS=None
    labelFITSName=None
    catFITSName=None

    modelTBFile=None

    rmsData = None #must has the same shape with Ny,and Nz of the data file

    averageRMS = None

    cleanCatName = None
    cleanFITSName = None
    ################# DBSCAN parameters

    cutoff_sigma = 2
    minPts = 4
    connectivity = 1


    ###################### seting of catalog selection

    minVox=16
    minChannel=3
    hasBeam=1
    minPeakSigma=5
    ###



    def __init__(self):
        self.getModelTB()


    def setDBSCANParameters(self,cutoff_sigma=2,minPts=4,connectivity=1):

        self.cutoff_sigma=cutoff_sigma
        self.minPts=minPts
        self.connectivity=connectivity

    def setCatalogSelectionCriteria(self,minVox=16,minChannel=3,hasBeam=1,minPeakSigma=5):
        self.minVox=minVox
        self.minChannel=minChannel
        self.hasBeam=hasBeam
        self.minPeakSigma=minPeakSigma



    def sumEdgeByCon1(self, extendMask):  # 7 in total
        raw = extendMask[1:-1, 1:-1, 1:-1]

        leftShiftZ = extendMask[0:-2, 1:-1, 1:-1]
        rightShiftZ = extendMask[2:, 1:-1, 1:-1]

        leftShiftY = extendMask[1:-1, 0: -2, 1:-1]
        rightShiftY = extendMask[1:-1, 2:, 1:-1]

        leftShiftX = extendMask[1:-1, 1:-1, 0:-2]
        rightShiftX = extendMask[1:-1, 1:-1, 2:]

        sumAll = raw + leftShiftZ + rightShiftZ + leftShiftY + rightShiftY + leftShiftX + rightShiftX

        return sumAll

    def sumEdgeByCon2(self, extendMask):  # 27 in total
        sumAll = extendMask[1:-1, 1:-1, 1:-1] * 0
        Nz, Ny, Nx = sumAll.shape
        for i in [-1, 0, 1]:
            for j in [-1, 0, 1]:
                for k in [-1, 0, 1]:

                    if np.sqrt(abs(i) + abs(j) + abs(k)) > 1.5:
                        continue

                    sumAll = sumAll + extendMask[1 + i:Nz + 1 + i, j + 1:Ny + 1 + j, k + 1: Nx + 1 + k]

        return sumAll

    def sumEdgeByCon3(self, extendMask):  # 27 in total
        raw = extendMask[1:-1, 1:-1, 1:-1]
        Nz, Ny, Nx = raw.shape
        sumAll = raw * 0
        for i in [-1, 0, 1]:
            for j in [-1, 0, 1]:
                for k in [-1, 0, 1]:
                    sumAll = sumAll + extendMask[1 + i:Nz + 1 + i, j + 1:Ny + 1 + j, k + 1: Nx + 1 + k]

        return sumAll


    ############

    def setOutputPath(self,targetPath):

        if os.path.isdir(targetPath):
            self.processPath=targetPath
        else:
            print "The target output path does not exist, please check your path: ", targetPath


    def getModelTB(self):
        """
        provide a TB file that can be used a model of catalog, new columns could be added to this table
        :return:
        """

        #targetFileName="https://raw.githubusercontent.com/SHAOqzyan/MWISPdbscan/master/minV3minP16_dendroCatTrunk.fit"

        targetFileName="minV3minP16_dendroCatTrunk.fit"

        if os.path.isfile(targetFileName):
            self.modelTBFile=targetFileName
        else:
            print "Downloading a modelTB file from the internet..."
            os.system("wget https://raw.githubusercontent.com/SHAOqzyan/MWISPdbscan/master/minV3minP16_dendroCatTrunk.fit")

        if os.path.isfile(targetFileName):
            self.modelTBFile=targetFileName

        else:
            print "Cannot download modelTB onthe internet, you need to manually download a table file named minV3minP16_dendroCatTrunk.fit (https://raw.githubusercontent.com/SHAOqzyan/MWISPdbscan/master/minV3minP16_dendroCatTrunk.fit)"




    def getDBSCANTag(self ):

        """
        return
        :param cutoff: #the sigma to which the overall datacube will be cut
        :param minpts:
        :param conType: #thre are
        :return:
        """

        saveName = "dbscanS{}P{}Con{}".format( self.cutoff_sigma,self.minPts,self.connectivity)

        return saveName



    def getUniformRMSData(self,averageRMS,Ny,Nx):



        return averageRMS+np.zeros(  (Ny,Nx) )

    def computeDBSCAN(  self  ):
        """
        There are two steps of computing DBSCAN, compute the DBSCAN label, 2 extract catlog, select post catalog and produce clean fits, would be doe seperatrely
        :param COdata:
        :param COHead:
        :param min_sigma: cutoff
        :param min_pix:
        :param connectivity:
        :param region:
        :param getMask:
        :param savePath:
        :param mimicDendro:
        :param rmsFITS:
        :param inputRMS:
        :return:
        """
        # pass

        if self.rawCOFITS is None:

            print "rawCOFITS need to be provided, stoping..."
            return

        COdata, COHead= doFITS.readFITS( self.rawCOFITS )

        Nz, Ny, Nx = COdata.shape


        if self.averageRMS is None and self.rmsFITS is None:
            print "You need either use averaged rms, or use rmsFITS, which can be accurate to each spectrum. Please check you input"


        if self.rmsFITS == None: #use average rms, te rmsFITS has a higher priority
            rmsData= self.getUniformRMSData(self.averageRMS,Ny,Nx)

        else:
            print "Using rms fits...", self.rmsFITS
            rmsData, rmsHead = myFITS.readFITS(self.rmsFITS)
        self.rmsData=rmsData #record the rms Data

        rmsCOData = COdata / rmsData
        goodValues = rmsCOData >= self.cutoff_sigma
        extendMask = np.zeros([Nz + 2, Ny + 2, Nx + 2], dtype=int)


        extendMask[1:-1, 1:-1, 1:-1] = goodValues  # [COdata>=minValue]=1

        s = generate_binary_structure(3, self.connectivity )

        if self.connectivity == 1:
            coreArray = self.sumEdgeByCon1(extendMask)

        if self.connectivity == 2:
            coreArray = self.sumEdgeByCon2(extendMask)

        if self.connectivity == 3:
            coreArray = self.sumEdgeByCon3(extendMask)

        coreArray = coreArray >= self.minPts
        coreArray[~goodValues] = False  # nan could be, #remove falsely, there is a possibility that, a bad value may have lots of pixels around and clould be
        #coreArray = coreArray + 0

        labeled_core, num_features = label(coreArray,    structure=s)  # first label core, then expand, otherwise, the expanding would wrongly connected

        selectExpand = np.logical_and(labeled_core == 0, goodValues)
        # expand labeled_core
        # coreLabelCopy=labeled_core.copy()

        expandTry = dilation(labeled_core,   s)  # first try to expand, then only keep those region that are not occupied previously
        # it is possible  that a molecular cloud may have less pixN than 8, because of the closeness of two
        # only expanded, one time
        labeled_core[selectExpand] = expandTry[selectExpand]

        labeled_array = labeled_core


        saveTag = self.getDBSCANTag( )

        print num_features, "features found!"



        baseName  =  os.path.basename( self.rawCOFITS )

        saveLabelFITSName= os.path.join(self.processPath,baseName[0:-5]+saveTag+".fits")

        fits.writeto(  saveLabelFITSName , labeled_array, header=COHead, overwrite=True)
        self.labelFITSName=saveLabelFITSName

        return saveLabelFITSName

    def getIndices(self, Z0, Y0, X0, values1D, choseID):

        cloudIndices = np.where(values1D == choseID)

        cX0 = X0[cloudIndices]
        cY0 = Y0[cloudIndices]
        cZ0 = Z0[cloudIndices]

        return tuple([cZ0, cY0, cX0])

    def getSaveCatName(self,doClean=True):

        if not doClean: # _Clean.fit
            return self.labelFITSName[0:-1]

        else:
            return self.labelFITSName[0:-5]+"_Clean.fit"


    def getCatFromLabelArray(self,doClean=True ):
        """
        Extract catalog from label fits, the minPix and rms is only used for saving
        #need a TB Model, if TB model does not exist, better downlooad from online
        This catalog has already been cleaned
        four criteria minPix"
        if you want to change the catlog selection criteria, please see setCatalogSelectionCariteria
        :param labelArray:
        :param head:
        :return:
        """

        if self.labelFITSName is None:

            print "The DBLABLE FITS is needed  "
            return

        dataCluster, headCluster = myFITS.readFITS(self.labelFITSName)
        Nz, Ny, Nx = dataCluster.shape

        if self.rawCOFITS is None:

            print "A COFITS is needed  "
            return

        dataCO, headCO = myFITS.readFITS(self.rawCOFITS)

        if self.averageRMS is not None and self.rawCOFITS is not None:

            self.getUniformRMSData(self.averageRMS,Ny,Nx)


        if self.rmsData is None:

            print "A rmsData is needed to select peakSgima"
            return



        saveCatName=self.getSaveCatName(doClean=doClean)


        clusterTBOld = Table.read(self.modelTBFile)



        minV = np.nanmin(dataCluster[0]) #noise mask value



        wcsCloud = WCS(headCluster)

        clusterIndex1D = np.where(dataCluster > minV)
        clusterValue1D = dataCluster[clusterIndex1D]

        Z0, Y0, X0 = clusterIndex1D

        newTB = Table(clusterTBOld[0])
        newTB["sum"] = newTB["flux"]


        newTB["l_rms"] = newTB["v_rms"]*0
        newTB["b_rms"] = newTB["v_rms"]*0

        newTB["pixN"] = newTB["v_rms"]*0
        newTB["peak"] = newTB["v_rms"]*0

        newTB["peakL"] = newTB["v_rms"]*0
        newTB["peakB"] = newTB["v_rms"]*0
        newTB["peakV"] = newTB["v_rms"]*0
        newTB["area_accurate"] = newTB["area_exact"]*0  # the column of area_accurate need to cos(b) facor

        # newTB["Nchannel"]=newTB["v_rms"] #number of channels, whis is used to sellect

        newTB["allChannel"] = newTB["v_rms"]*0  # number channel involved
        newTB["has22"] = newTB["v_rms"]*0  # number channel involved


        #clear TB
        newTB["sum"].unit="K"
        del newTB["area_ellipse"]
        del newTB["flux"]
        del newTB["major_sigma"]
        del newTB["minor_sigma"]
        del newTB["radius"]
        del newTB["position_angle"]





        newTB["v_cen"].unit="km/s"
        newTB["x_cen"].unit="deg"
        newTB["y_cen"].unit="deg"
        newTB["l_rms"].unit="deg"
        newTB["b_rms"].unit="deg"
        newTB["peakB"].unit=""
        newTB["peakL"].unit=""
        newTB["peakV"].unit=""
        newTB["peak"].unit="K"

        newTB["allChannel"].unit=""
        newTB["has22"].unit=""






        zeroProjection = np.zeros((Ny, Nx))  # one zero channel, used to get the projection area and
        zeroProjectionExtend = np.zeros((Ny + 1, Nx + 1))

        idCol = "_idx"

        # count all clusters

        # ids,count=np.unique(dataCluster,return_counts=True )
        ids, count = np.unique(clusterValue1D, return_counts=True)
        GoodIDs = ids
        GoodCount = count
        print "Total number of trunks? ", len(ids)
        # print "Total number of Good Trunks? ",len(GoodIDs)

        # dataCO,headCO=doFITS.readFITS( CO12FITS )
        widgets = ['Recalculating cloud parameters: ', Percentage(), ' ', Bar(marker='0', left='[', right=']'), ' ', ETA(),
                   ' ', FileTransferSpeed()]  # see docs for other options

        catTB = newTB.copy()
        catTB.remove_row(0)

        # zeroP
        # remove any cloud with voxels less than minPix
        selectVox = GoodCount >= self.minVox  # #criteria 1

        if  doClean:
            GoodCount = GoodCount[selectVox]
            GoodIDs = GoodIDs[selectVox]
            print len(GoodIDs), "cluster has more than {} voxels.".format(self.minVox)

        pbar = ProgressBar(widgets=widgets, maxval=len(GoodIDs))
        pbar.start()

        ####

        for i in range(len(GoodIDs)):

            # i would be the newID
            newID = GoodIDs[i]
            pixN = GoodCount[i]

            newRow = newTB[0]

            newRow[idCol] = newID

            cloudIndex = self.getIndices(Z0, Y0, X0, clusterValue1D, newID)
            coValues = dataCO[cloudIndex]
            cloudV = cloudIndex[0]
            cloudB = cloudIndex[1]
            cloudL = cloudIndex[2]
            # step1, check peak
            # peak should be done first

            peakIndex = coValues.argmax()
            peakV = cloudV[peakIndex]
            peakB = cloudB[peakIndex]
            peakL = cloudL[peakIndex]
            peak=coValues[peakIndex]
            # #criteria 2

            if peak < self.minPeakSigma * self.rmsData[peakB,peakL] and doClean:  #  accurate to lines
                continue #do not consider the minimum peaks


            # #criteria 3, check consecutive channels

            ############### ###
            diffVs = np.unique(cloudV)

            if len(diffVs) < self.minChannel and doClean:  # reject all cloud that has channels less than 3 channels
                continue



            # get the exact peak position, which would be used to

            projectIndex = tuple([cloudB, cloudL])

            zeroProjection[projectIndex] = 1
            zeroProjectionExtend[0:-1, 0:-1] = zeroProjection
            sum22 = zeroProjectionExtend[0:-1, 0:-1] + zeroProjectionExtend[0:-1, 1:] + zeroProjectionExtend[1:,
                                                                                        0:-1] + zeroProjectionExtend[1:, 1:]

            # if any pixel>4:
            #

            if 4 in sum22:
                newRow["has22"] = 1

            else:  # Serious bug
                newRow["has22"] = 0



            if self.hasBeam and newRow["has22"]==0 and doClean: # #criteria 4
                zeroProjection[projectIndex] = 0
                zeroProjectionExtend[0:-1, 0:-1] = zeroProjection
                continue
            # newRow["has22"] = 0

            # calculate the accurate, this is only used for data at high Galactic latitude, and useless form the MWISP survey
            # to increase the speed, this part is commentted

            # indexB2D,indexL2D=np.where(zeroProjection==1 )
            # _,BS2D, LS2D = wcsCloud.wcs_pix2world(indexL2D,indexB2D,  0, 0)
            # area_accurate=np.sum( np.cos( np.deg2rad(BS2D) )    )*0.25
            # newRow["area_accurate"]= area_accurate

            sumCO = np.sum(coValues)

            Vcen, Vrms = doFITS.weighted_avg_and_std(cloudV, coValues)
            Bcen, Brms = doFITS.weighted_avg_and_std(cloudB, coValues)
            Lcen, Lrms = doFITS.weighted_avg_and_std(cloudL, coValues)

            # calculate the exact area

            # LBcore = zip(cloudB, cloudL)
            # pixelsN= {}.fromkeys(LBcore).keys() #len( set(LBcore) )
            # area_exact=len(pixelsN)*0.25 #arc mins square
            area_exact = np.sum(zeroProjection) * 0.25

            newRow["area_exact"] = area_exact

            # dataClusterNew[cloudIndex] =newID

            # save values
            newRow["pixN"] = pixN
            newRow["peak"] = peak

            newRow["peakV"] = peakV
            newRow["peakB"] = peakB
            newRow["peakL"] = peakL

            # newRow["peak2"]= peak2

            newRow["sum"] = sumCO

            newRow["x_cen"], newRow["y_cen"], newRow["v_cen"] = wcsCloud.wcs_pix2world(Lcen, Bcen, Vcen, 0)
            newRow["v_cen"] = newRow["v_cen"] / 1000.
            dv = headCluster["CDELT3"] / 1000.  # km/s

            dl = abs(headCluster["CDELT1"])  # deg

            newRow["v_rms"] = Vrms * dv
            newRow["l_rms"] = Lrms * dl
            newRow["b_rms"] = Brms * dl

            # _, Nchan=np.unique( cloudV, return_counts=True)

            # newRow["Nchannel"] =    np.max(P3Value)# if there is a three consecutive spectra in the cloud
            newRow["allChannel"] = len(diffVs)

            catTB.add_row(newRow)

            zeroProjection[projectIndex] = 0
            zeroProjectionExtend[0:-1, 0:-1] = zeroProjection

            pbar.update(i)

        pbar.finish()
        # save the clouds

        # fits.writeto(self.regionName+"NewCloud.fits", dataClusterNew,header=headCluster,overwrite=True   )
        catTB.write(saveCatName, overwrite=True)
        if doClean:
            self.cleanCatName=saveCatName
            #self.catFITSName=saveCatName

        else:
            self.catFITSName=saveCatName
        return saveCatName

    def getEmptyTB(self,tb):
        newTB=Table( tb[0] )

        newTB.remove_row(0)

        return newTB


    def cleanTB(self):

        """
        select clouds that satifies the selection criteria
        :return:
        """



        if self.catFITSName is None:
            print "no catalog file provided"

            return

        if self.rmsData is None:
            print "You need to provide the rms FITS data to select peak accordingly"


        saveCleanTBName=self.catFITSName[0:-4]+"_Clean.fit"

        rawTB=Table.read(self.catFITSName )

        #first voxel
        filterTB = rawTB[rawTB["pixN"] >= self.minVox ]


        filterTB = filterTB[filterTB["has22"] ==self.hasBeam ]

        #select by minCkannel
        filterTB = filterTB[filterTB["allChannel"] >= self.minChannel  ]

        #second peak

        peakB = map(int,filterTB["peakB"] )
        peakL = map(int,filterTB["peakL"] )
        peakCoordinate=tuple( [peakB, peakL] )


        rmsArray=self.rmsData[peakCoordinate]
        filterTB = filterTB[filterTB["peak"] >= self.minPeakSigma * rmsArray ]

        filterTB.write(saveCleanTBName, overwrite=True)
        self.cleanCatName=saveCleanTBName

        return  filterTB




    def produceCleanFITS(self):
        """
        remove noise cluster according to the label fits and catFITS name
        :return:
        """
        if self.labelFITSName is None :

            print "You need to provde the label FITS  "

            return

        if self.cleanCatName is not None:

            targetTBFile= self.cleanCatName #table has alread been cleaned
        elif self.catFITSName is not None:
            #catalog has been calculated but not cleaned, need to clean table first, then do clean
            self.cleanTB()
            targetTBFile= self.cleanCatName #table has alread been cleaned

        elif self.catFITSName is None:
            "You need to provde a raw or a cleaned fits first to produce a clean fits"
            return

        saveCleanFITSName=targetTBFile+"s"

        TB=Table.read(targetTBFile)
        dataCluster,headCluster=doFITS.readFITS(self.labelFITSName)
        noiseLabel = np.min(dataCluster[0])
        clusterIndex1D = np.where(dataCluster > 0)
        clusterValue1D = dataCluster[clusterIndex1D]
        Z0, Y0, X0 = clusterIndex1D

        #

        widgets = ['Cleaning label fits:', Percentage(), ' ', Bar(marker='0', left='[', right=']'), ' ', ETA(), ' ',
                   FileTransferSpeed()]  # see docs for other options

        pbar = ProgressBar(widgets=widgets, maxval=len(TB))
        pbar.start()

        indexRun = 0

        returnCluster = np.zeros_like(dataCluster) + noiseLabel

        for eachDBRow in TB:
            indexRun = indexRun + 1
            pbar.update(indexRun)
            cloudID = eachDBRow["_idx"]

            cloudIndex = self.getIndices(Z0, Y0, X0, clusterValue1D, cloudID)
            returnCluster[cloudIndex] = cloudID  # remove this cluster and do not record this cluster

        pbar.finish()
        fits.writeto(saveCleanFITSName,returnCluster,header=headCluster,overwrite=True)

        self.cleanFITSName=saveCleanFITSName
        return saveCleanFITSName

