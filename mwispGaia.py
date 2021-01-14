import configparser
import os
import matplotlib.pyplot as plt
from astropy.table import Column
import matplotlib.patches as mpatches
import pywcsgrid2
import img_scale
from mpl_toolkits.axes_grid1.axes_rgb import imshow_rgb
import numpy as np
#from mpl_toolkits.axes_grid.anchored_artists import AnchoredText
from astropy import units as u
from matplotlib.offsetbox import AnchoredText
from matplotlib import rc
import numpy.ma as ma
from pyds9 import DS9
from astropy.convolution import convolve, Gaussian2DKernel, Tophat2DKernel
import matplotlib.colors as colors
import pyregion
from  distanceTB import disTB
from astropy.modeling.models import Gaussian2D
from astropy.convolution import CustomKernel, Model2DKernel
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes, mark_inset,inset_axes
from numpy import mean, sqrt, square
import pywcsgrid2
from astropy.wcs import WCS
import corner
from astropy.table import Table,vstack
from scipy.spatial import ConvexHull
import matplotlib as mpl
from matplotlib.patches import Rectangle,FancyArrowPatch
import scipy.spatial.distance as calLBDis

from astropy.io import fits
from matplotlib import rcParams

from mpl_toolkits.axes_grid1.axes_grid import ImageGrid
from mpl_toolkits.axes_grid1 import Grid
from matplotlib import gridspec
from myPYTHON import *


from matplotlib.colors import LogNorm
from progressbar import *

from dendroTB import dendroCat


import glob

#from gaiaTBCor import  GAIATB
#from gaiaTBCor import  GAIATB
from gaiaTBCorAGEDR3 import GAIATBAGEDR3 as GAIATB
#import gaiaTBCorAGEDR3.GAIATBAGEDR3 as GAIATB

from  myGAIA import GAIADIS

from spectral_cube import SpectralCube

from gaiaAVFlagTB import GAIAAVFlag
#from gaiaAVFlagTBEDR3 import GAIAAVFlag #Gaia EDR2



#from distanceAnalysis import gaiaMWISP
from myTree import dendroTree
import sys
import pickle

gaiaDis=GAIADIS()
doAV=GAIAAVFlag()
doAG=GAIATB()

doFITS=myFITS()

class dbscanDis: #a aclass dealing with distances to DBSCAN clouds

    rootPath = "/home/qzyan/WORK/projects/maddalena/dendroDisPath/" #a path, that saves distance resolts

    configFile = rootPath + "regions.ini" #a file that colects distance paths

    maskPathName = "maskPath"
    treeFileName = "treeFile"
    CO12FITSName = "CO12FITS"
    dendroCatName = "dendroCat"
    lvPath = "lvPath"

    foregroundFolder = "foregroundFiles"  # the folder used to store



    # signalLevel=4.
    # noiseLevel=1.

    tmpFolder = "tmpFiles"  # folder used to save temperary files

    testFigures = "testFigures"
    distFigures = "distFigures"
    goodFigures = "goodFigures"

    extraDis = 300

    disSources = None  # sources that have distances

    disCuts = None  # sources that have distances
    signalLevs = None

    useAV = False

    dv = 0.15873764455
    beamsize = 49. / 60. / 60.  # degree
    resolutionCOPix = 30. / 60. / 60.  # degre

    offExtendMask = "offExtend"

    wiseBack = None

    dataCO12 = None
    headCO12 = None
    binCol = "bin"

    def __init__(self, sourceName, useAV=False):

        """
        """

        # the sourceName has to be assinged

        self.config = configparser.ConfigParser()
        self.config.read(self.configFile)

        self.sourceName = sourceName

        self.useAV = useAV

        if not sourceName in self.config.sections():

            print "configue of source does not exist!"




        else:  # read configures
            self.maskPath = self.config[sourceName][self.maskPathName]
            #self.treeFile = self.config[sourceName][self.treeFileName]
            self.CO12FITS = self.config[sourceName][self.CO12FITSName]
            self.dendroCat = self.config[sourceName][self.dendroCatName]

            #self.doTree = dendroTree(self.treeFile, dendroCat=self.dendroCat)

            #self.cloudLvPath = self.config[sourceName][self.lvPath]

            self.goodDistanceTB = self.rootPath + sourceName + "/{}goodDisTB.fit".format(sourceName)
        # pvFITS

        # initiallize....
        self.sourcePath = self.rootPath + sourceName + "/"

        self.pvFITS = self.sourcePath + sourceName + "_LV.fits"

        self.dataPath = self.rootPath + sourceName + "/data/"

        self.dendroTB = Table.read(self.dendroCat)

        if not os.path.isdir(self.sourcePath):
            os.mkdir(self.sourcePath)

        self.testFigurePath = self.sourcePath + self.testFigures + "/"

        if not os.path.isdir(self.testFigurePath):
            os.mkdir(self.testFigurePath)

        self.foregroundPath = self.sourcePath + self.foregroundFolder + "/"
        if not os.path.isdir(self.foregroundPath):
            os.mkdir(self.foregroundPath)

        self.tmpPath = self.sourcePath + self.tmpFolder + "/"
        if not os.path.isdir(self.tmpPath):
            os.mkdir(self.tmpPath)

        self.distFigurePath = self.sourcePath + self.distFigures + "/"

        if not os.path.isdir(self.distFigurePath):
            os.mkdir(self.distFigurePath)
        ###########################################################################
        self.distFigureGoodPath = self.sourcePath + self.goodFigures + "/"
        if not os.path.isdir(self.distFigureGoodPath):
            os.mkdir(self.distFigureGoodPath)

        if not self.useAV:

            self.disCat = self.sourcePath + "disCat{}.fit".format(sourceName)
            self.doDis = disTB(self.disCat)

        else:

            self.disCat = self.sourcePath + "disCat{}_AV.fit".format(sourceName)
            self.doDis = disTB(self.disCat)

    def getForegroundFITS(self, cloudID):
        """
        according the velocity range to produce foreground fits files
        :param ID:
        :return:
        """
        cloudID = str(cloudID)

        if "Cloud" not in cloudID:
            cloudID = "Cloud{}".format(cloudID)

        # searchIntFITS=  glob.glob(self.maskPath+"*{}_level*int.fits".format(cloudID))
        # searchMaskFITS=  glob.glob(self.maskPath+"*{}_level*mask.fits".format(cloudID))



        searchStr = os.path.join(  self.maskPath,    "*{}_*fore.fits".format(cloudID ) )
        searchIntFITSFore = glob.glob(  searchStr   )

        if len( searchIntFITSFore ) ==0:
            return None

        return  searchIntFITSFore[0]

    def getAllCandidateID(self,candidatePath  ):
        """
        get the ID list for the candidate
        :return:
        """
        #candidatePath = "./goodCandidate/"

        pngFiles = glob.glob(candidatePath + "Cloud*Test.png")

        idList = []

        for eachFile in pngFiles:
            idStr = eachFile.split("Test")[0]
            idStr = idStr.split("Cloud")[1]

            idList.append(int(idStr))

        print  len(idList)

        return idList


    def fastTestByID( self, ID, NL=1., SL=4., lowerDisCut=1., cutDis=3000., lExpand=0.15, calibrateAngle=1., bExpand=0.5,
                   testAllBins=False, useMask=True, upCO=None, useBin=False, useAllBins=True, offBinSize=3, onBinSize=1,
                   extendRatio=0.5, paraErrorCut=0.2, useAV=True, pureInt=False, extraDis=300., useBaseline=True,
                   pureTest=False, useForegroundFITS=False, foreOn=True, foreOff=True, lRange=None, bRange=None,
                   extendRegion=False, extendSize=1., offQuadrant=None, useHighAv=False,saveTag="",savePath= "" ):

        """
        #can only use to test small molecular clouds, if its distance is detectable
        :param SL:
        :param lowerDisCut:
        :param cutDis:
        :param lExpand:
        :param calibrateAngle:
        :param bExpand:
        :param testAllBins:
        :param useMask:
        :param upCO:
        :param useBin:
        :param useAllBins:
        :param offBinSize:
        :param onBinSize:
        :param extendRatio:
        :param paraErrorCut:
        :param useAV:
        :param pureInt:
        :param extraDis:
        :param useBaseline:
        :param pureTest:
        :param useForegroundFITS:
        :param foreOn:
        :param foreOff:
        :param lRange:
        :param bRange:
        :param extendRegion:
        :param extendSize:
        :param offQuadrant:
        :param useHighAv:
        :param saveTag:
        :return:
        """

        foregroundFITS = None

        foreGroundCut = 1.5

        if useForegroundFITS:  # get the foreground fitsname and get the nose

            foregroundFITS = self.getForegroundFITS(ID)  #

            foregroundRMS = doFITS.getCOFITSNoise(foregroundFITS)

            foreGroundCut = max(foreGroundCut, 3 * foregroundRMS)

            print "Using foreground fits (noise:{})...cut off: {}".format(foregroundRMS, foreGroundCut)

        gaiaDis.useAV = useAV
        noiseLevel = NL  # RMSCO12*NL
        signalLevel = SL  # MSCO12*SL

        cloudID = 'Cloud{}'.format(ID)
        # maskFITS=		glob.glob(self.maskPath+"*"+cloudID  +"_*mask.fits")[0]
        # intFITS=		glob.glob(self.maskPath+"*"+cloudID  +"_*int.fits")[0]

        intFITS, maskFITS = self.getIntMaskByID(ID, pureInt=pureInt)

        intData, intHead = myFITS.readFITS(intFITS)
        maskData, maskHead = myFITS.readFITS(maskFITS)

        newDisRow = self.doDis.getRowByName(cloudID)
        if newDisRow == None:
            newDisRow = self.doDis.getEmptyRow()

        newDisRow[disTB.fitsFile] = intFITS
        newDisRow[disTB.sourceName] = cloudID

        newDisRow[disTB.noiseLevel] = noiseLevel
        newDisRow[disTB.signalLevel] = signalLevel

        newDisRow[disTB.cutDistanceLower] = lowerDisCut
        newDisRow[disTB.cutDistanceUpper] = cutDis

        lRangeMask, bRangeMask = self.getLBrangeWithMask(maskFITS)

        if not useMask:
            lRangeMask=lRange
            bRangeMask=bRange

        if lRange is not None and bRange is not None:
            lRangeMask=lRange
            bRangeMask=bRange


        backWCS = WCS(intFITS, naxis=(1, 2))

        if lExpand is not  None and bExpand is not  None:

            searchlRange, searchbRange = self.expandLBRangeByDeg(lRangeMask, bRangeMask, intData, backWCS, lExpand=lExpand,
                                                                 bExpand=bExpand)
        else:
            searchlRange, searchbRange = self.expandLBRange(lRangeMask, bRangeMask, intData, backWCS, extendRatio)

        lowerPara = 1000. / (cutDis + self.extraDis)
        upperPara = 1000. / lowerDisCut


        if useAV:
            gaiaStars = doAV.getByLBRange(searchlRange, searchbRange, mimicAG=True, lowerPara=lowerPara,
                                          paraError=paraErrorCut, upperPara=upperPara, useHighAv=useHighAv)

        else:
            gaiaStars = doAG.getByLBRange(searchlRange, searchbRange, lowerPara=lowerPara, calDis=True, upperPara= upperPara,
                                          paraError=paraErrorCut)

        if foregroundFITS != None:  # assuming this is foreground CO  fits

            gaiaStars = gaiaDis.assignOnsourceGaia(gaiaStars, foregroundFITS, gaiaDis.foreCOCol)



        onCloudStars, offCloudStars = gaiaDis.getOnAndOffStars(gaiaStars, NL=noiseLevel, SL=signalLevel, maskFITS=maskFITS,       intFITS=intFITS, useMask=useMask)

        # by here, to test if the off-cloud stars would be better if we extended the region of molecular cloud a little big, for exapmle, 1 beam?
        if extendRegion:
            extendMask = self.getExtendMask(ID, extendBeam=extendSize)

            offCloudStars = gaiaDis.assignOnsourceGaia(offCloudStars, extendMask, self.offExtendMask)
            offCloudStars = offCloudStars[offCloudStars[self.offExtendMask] < 0.5]  # only do this on off cloud star

        if foreOn and foregroundFITS is not  None:
            print "Removing foreground stars for on-cloud stars...:before ",len( onCloudStars )
            onCloudStars = onCloudStars[onCloudStars[gaiaDis.foreCOCol] < foreGroundCut]  # only do this on off cloud star
            print "After ",len( onCloudStars )


        if foreOff and foregroundFITS is not  None:

            print "Removing foreground stars for on-cloud stars...:before ",len( offCloudStars )
            offCloudStars = offCloudStars[offCloudStars[gaiaDis.foreCOCol] < foreGroundCut]  # only do this on off cloud star
            print "After ",len( offCloudStars )




        if lRange is not  None:
            onCloudStars = onCloudStars[onCloudStars["l"] < max(lRange)]
            onCloudStars = onCloudStars[onCloudStars["l"] > min(lRange)]

            offCloudStars = offCloudStars[offCloudStars["l"] < max(lRange)]
            offCloudStars = offCloudStars[offCloudStars["l"] > min(lRange)]

        if bRange is not  None:
            onCloudStars = onCloudStars[onCloudStars["b"] < max(bRange)]
            onCloudStars = onCloudStars[onCloudStars["b"] > min(bRange)]

            offCloudStars = offCloudStars[offCloudStars["b"] < max(bRange)]
            offCloudStars = offCloudStars[offCloudStars["b"] > min(bRange)]

        denDroCat = Table.read(self.dendroCat)

        cloudRow = denDroCat[denDroCat["_idx"] == ID]

        l = float(cloudRow["x_cen"])
        b = float(cloudRow["y_cen"])
        v = float(cloudRow["v_cen"]) # / 1000.
        figNameMark = self.getCloudNameByLB(l, b)

        # r= float(  cloudRow["major_sigma"] )  #/3600. )

        newDisRow[disTB.l] = l
        newDisRow[disTB.b] = b
        newDisRow[disTB.cloudVlsr] = v

        # RMSCO12=self.getRMS(intFITS)

        newDisRow[disTB.Note] = figNameMark

        newDisRow[disTB.boxCenterL] = np.mean(searchlRange)
        newDisRow[disTB.boxCenterB] = np.mean(searchbRange)
        newDisRow[disTB.boxSizeL] = abs(searchlRange[0] - searchlRange[1])
        newDisRow[disTB.boxSizeB] = abs(searchbRange[0] - searchbRange[1])

        COsum = np.sum(intData * maskData)
        newDisRow[disTB.COsum] = COsum

        offCloudStars = self.cleanByAGError(offCloudStars)

        # usually foregroundFITS only apply to off cloud stars, but for clouds hat are fore, it applies to both on and off cloud stars for far molecular clouds

        # assign the star with foreground fits

        if upCO != None:
            onCloudStars = onCloudStars[onCloudStars["coint"] <= upCO]

        print " --------- {}: on and off cloud stars are {} and {} ".format(cloudID, len(onCloudStars), len(offCloudStars))
        # save on cloud stars
        if useAV:
            onCloudStars.write(self.tmpPath + str(cloudID) + "OnCloudAV.fit", overwrite=True)
            offCloudStars.write(self.tmpPath + str(cloudID) + "OffCloudAV.fit", overwrite=True)

        else:
            onCloudStars.write(self.tmpPath + str(cloudID) + "OnCloudAG.fit", overwrite=True)
            offCloudStars.write(self.tmpPath + str(cloudID) + "OffCloudAG.fit", overwrite=True)

        #
        onCloudStarsRaw = onCloudStars.copy()
        offCloudStarsRaw = offCloudStars.copy()



        # print neighbourhood.argmin(axis=1) # for each on cloud bin

        # print neighbourhood.argmin(axis=0)

        # print "Tesing small number of on-cloud stars......"

        # onCloudStars=gaiaDis.getRandomRows(onCloudStars,50)

        # print " --------- {}: After selection on and off cloud stars are {} and {} ".format(cloudID, len(onCloudStars), len(offCloudStars)  )

        # get baseline
        baseLine = gaiaDis.getBaseLine(offCloudStars)
        onCloudStars = onCloudStars[onCloudStars[doAG.GAIA_distance] <= cutDis]
        offCloudStars = offCloudStars[offCloudStars[doAG.GAIA_distance] <= cutDis]
        offCloudStars = offCloudStars[   offCloudStars["distance_err"] > 0]  # distance_errr cannot be zeoro, possibly fore extremelly close stars
        offCloudStars.sort("distance")
        #draw

        fig = plt.figure(figsize=(12, 6))
        rc('text', usetex=True)
        rc('font', **{'family': 'sans-serif', 'size': 16, 'serif': ['Helvetica']})
        mpl.rcParams['text.latex.preamble'] = [
            r'\usepackage{tgheros}',  # helvetica font
            r'\usepackage{sansmath}',  # math-font matching  helvetica
            r'\sansmath'  # actually tell tex to use it!
            r'\usepackage{siunitx}',  # micro symbols
            r'\sisetup{detect-all}',  # force siunitx to use the fonts
        ]

        axRaw = fig.add_subplot(1, 2, 1)

        #axRaw.scatter(offCloudStars["distance"] ,  offCloudStars["a_g_val"] , s=7,color='gray')
        axRaw.scatter(onCloudStars["distance"] ,  onCloudStars["a_g_val"] , s=5,color='green',label=None)


        axRaw.plot(offCloudStars["distance"] ,   baseLine.predict(offCloudStars["distance"] )   , lw=2,color='red',alpha = 0.8 ,label=None)

        axBaseline  = fig.add_subplot(1, 2, 2)

        onDis =  onCloudStars["distance"]
        onAV =  onCloudStars["a_g_val"]  -  baseLine.predict(onDis )

        if len(onAV)>0:
            axBaseline.scatter( onDis ,  onAV , s=5,color='gray')


            #draw base line
            axBaseline.set_xlabel(r"distance (kpc)")
            axRaw.set_xlabel(r"distance (kpc)")
            axRaw.set_ylabel(r"AV (mag)")

            mu1List,ratio ,disArray= self.guessGoodOrBad(onDis,onAV )

            #minimum position of ratio

            guessDis=disArray[  ratio.argmin() ]
            axBaseline.plot( [guessDis,guessDis],  [min(onAV), max(onAV)],lw=0.8,color='black'  )
            axRaw.plot( [guessDis,guessDis],  [min(onAV), max(onAV)],lw=0.8,color='black',label="cloud: {}".format(ID )  )

            axRaw.legend(loc= 4,handlelength=0.3)

            #axRaw.plot( disArray, ratio,lw=0.8,color='black'  )

            axBaseline.plot( disArray, ratio,lw=0.8,color='blue'  )
            #axBaseline.plot( disArray,mu2List,lw=0.8,color='red'  )



        plt.savefig(savePath+"Cloud{}Test.png".format(ID), bbox_inches='tight', dpi=300)
        plt.clf()



    def guessGoodOrBad(self,onDis,onAV):
        """
        find a way to test, if this cloud is good
        :param onDis:
        :param onAV:
        :return:
        """
        #
        AVcut=0.5
        examDis =np.linspace(0,3000,300)

        mu1List=[]
        mu2List=[]

        ratioList = []




        for eachD in examDis:
            foreGroundStars= onAV[onDis<=eachD ]
            backgroundStars=  onAV[onDis>eachD ]


            q1N =    len(   backgroundStars[backgroundStars<AVcut] )

            q2N = len( backgroundStars) -q1N 



            q3N =    len(   foreGroundStars[foreGroundStars>=AVcut] )

            q4N = len( foreGroundStars) -q3N
            ratio= (q3N +q1N )/1./len(onAV )*1.5
            ratioList.append( ratio )



        return np.asarray(mu1List),np.asarray(ratioList), examDis
            


    def calDisByID(self, ID, NL=1., SL=4., lowerDisCut=1., cutDis=3500., lExpand=0.5, calibrateAngle=1., bExpand=0.15,
                   testAllBins=False, useMask=True, upCO=None, useBin=False, useAllBins=True, offBinSize=3, onBinSize=1,
                   extendRatio=0.5, paraErrorCut=0.2, useAV=True, pureInt=False, extraDis=300., useBaseline=True,
                   pureTest=False, useForegroundFITS=False, foreOn=True, foreOff=True, lRange=None, bRange=None,
                   extendRegion=False, extendSize=1., offQuadrant=None, useHighAv=False,saveTag="" ,inputName= None,legendCol=1):
        """
        # this function take a cloudName, search mask and in file, extract on and off cloud stars, then calculate distances

        #and save the figure

        is it fitting a polynomial better?

        better try this with the ratio one
        foregroundFITS is the FITS used to trace those sources affected by local molecular clouds

        #lRange and bRange are manually choose regions for the on-clouds stars. Usually these two parameters are only used for checking if there two components in clouds that are
        identified by dendrogram

        offQuadrant: 1,2,3,4
        :param ID:
        :param NL:
        :param SL:
        :param lowerDisCut:
        :param cutDis:
        :param extendRatio:
        :param paraErrorCut:
        :param useAV:
        :param extraDis:
        :param useBaseline:
        :param pureTest:
        :param useForegroundFITS:
        :param foreOn:
        :param foreOff:
        :param lRange:
        :param bRange:
        :param extendRegion: # extend the cloud region of molecular clouds to see if we can get a better baseline? By default this should not be done
        :param extendSize: # for the first try, we extend the cloud by 1 beam
        :return:
        """

        foregroundFITS = None

        foreGroundCut = 1.5

        if useForegroundFITS:  # get the foreground fitsname and get the nose

            foregroundFITS = self.getForegroundFITS(ID)  #

            foregroundRMS = doFITS.getCOFITSNoise(foregroundFITS)

            foreGroundCut = max(foreGroundCut, 3 * foregroundRMS)

            print "Using foreground fits (noise:{})...cut off: {}".format(foregroundRMS, foreGroundCut)

        gaiaDis.useAV = useAV
        noiseLevel = NL  # RMSCO12*NL
        signalLevel = SL  # MSCO12*SL

        cloudID = 'Cloud{}'.format(ID)
        # maskFITS=		glob.glob(self.maskPath+"*"+cloudID  +"_*mask.fits")[0]
        # intFITS=		glob.glob(self.maskPath+"*"+cloudID  +"_*int.fits")[0]

        intFITS, maskFITS = self.getIntMaskByID(ID, pureInt=pureInt)

        intData, intHead = myFITS.readFITS(intFITS)
        maskData, maskHead = myFITS.readFITS(maskFITS)

        newDisRow = self.doDis.getRowByName(cloudID)
        if newDisRow == None:
            newDisRow = self.doDis.getEmptyRow()

        newDisRow[disTB.fitsFile] = intFITS
        newDisRow[disTB.sourceName] = cloudID

        newDisRow[disTB.noiseLevel] = noiseLevel
        newDisRow[disTB.signalLevel] = signalLevel

        newDisRow[disTB.cutDistanceLower] = lowerDisCut
        newDisRow[disTB.cutDistanceUpper] = cutDis

        lRangeMask, bRangeMask = self.getLBrangeWithMask(maskFITS)

        if not useMask:
            lRangeMask=lRange
            bRangeMask=bRange

        if lRange is not None and bRange is not None:
            lRangeMask=lRange
            bRangeMask=bRange


        backWCS = WCS(intFITS, naxis=(1, 2))

        if lExpand is not  None and bExpand is not  None:

            searchlRange, searchbRange = self.expandLBRangeByDeg(lRangeMask, bRangeMask, intData, backWCS, lExpand=lExpand,
                                                                 bExpand=bExpand)
        else:
            searchlRange, searchbRange = self.expandLBRange(lRangeMask, bRangeMask, intData, backWCS, extendRatio)

        lowerPara = 1000. / (cutDis + self.extraDis)
        upperPara = 1000. / lowerDisCut


        if useAV:
            gaiaStars = doAV.getByLBRange(searchlRange, searchbRange, mimicAG=True, lowerPara=lowerPara,
                                          paraError=paraErrorCut, upperPara=upperPara, useHighAv=useHighAv)

        else:
            gaiaStars = doAG.getByLBRange(searchlRange, searchbRange, lowerPara=lowerPara, calDis=True, upperPara= upperPara,
                                          paraError=paraErrorCut)

        if foregroundFITS != None:  # assuming this is foreground CO  fits

            gaiaStars = gaiaDis.assignOnsourceGaia(gaiaStars, foregroundFITS, gaiaDis.foreCOCol)



        onCloudStars, offCloudStars = gaiaDis.getOnAndOffStars(gaiaStars, NL=noiseLevel, SL=signalLevel, maskFITS=maskFITS,       intFITS=intFITS, useMask=useMask)

        # by here, to test if the off-cloud stars would be better if we extended the region of molecular cloud a little big, for exapmle, 1 beam?
        if extendRegion:
            extendMask = self.getExtendMask(ID, extendBeam=extendSize)

            offCloudStars = gaiaDis.assignOnsourceGaia(offCloudStars, extendMask, self.offExtendMask)
            offCloudStars = offCloudStars[offCloudStars[self.offExtendMask] < 0.5]  # only do this on off cloud star

        if foreOn and foregroundFITS is not  None:
            print "Removing foreground stars for on-cloud stars...:before ",len( onCloudStars )
            onCloudStars = onCloudStars[onCloudStars[gaiaDis.foreCOCol] < foreGroundCut]  # only do this on off cloud star
            print "After ",len( onCloudStars )


        if foreOff and foregroundFITS is not  None:

            print "Removing foreground stars for on-cloud stars...:before ",len( offCloudStars )
            offCloudStars = offCloudStars[offCloudStars[gaiaDis.foreCOCol] < foreGroundCut]  # only do this on off cloud star
            print "After ",len( offCloudStars )




        if lRange is not  None:
            onCloudStars = onCloudStars[onCloudStars["l"] < max(lRange)]
            onCloudStars = onCloudStars[onCloudStars["l"] > min(lRange)]

            offCloudStars = offCloudStars[offCloudStars["l"] < max(lRange)]
            offCloudStars = offCloudStars[offCloudStars["l"] > min(lRange)]

        if bRange is not  None:
            onCloudStars = onCloudStars[onCloudStars["b"] < max(bRange)]
            onCloudStars = onCloudStars[onCloudStars["b"] > min(bRange)]

            offCloudStars = offCloudStars[offCloudStars["b"] < max(bRange)]
            offCloudStars = offCloudStars[offCloudStars["b"] > min(bRange)]

        denDroCat = Table.read(self.dendroCat)

        cloudRow = denDroCat[denDroCat["_idx"] == ID]

        l = float(cloudRow["x_cen"])
        b = float(cloudRow["y_cen"])
        v = float(cloudRow["v_cen"]) # / 1000.
        figNameMark = self.getCloudNameByLB(l, b)

        if inputName is not None:
            figNameMark= inputName

        # r= float(  cloudRow["major_sigma"] )  #/3600. )

        newDisRow[disTB.l] = l
        newDisRow[disTB.b] = b
        newDisRow[disTB.cloudVlsr] = v

        # RMSCO12=self.getRMS(intFITS)

        newDisRow[disTB.Note] = figNameMark

        newDisRow[disTB.boxCenterL] = np.mean(searchlRange)
        newDisRow[disTB.boxCenterB] = np.mean(searchbRange)
        newDisRow[disTB.boxSizeL] = abs(searchlRange[0] - searchlRange[1])
        newDisRow[disTB.boxSizeB] = abs(searchbRange[0] - searchbRange[1])

        COsum = np.sum(intData * maskData)
        newDisRow[disTB.COsum] = COsum

        offCloudStars = self.cleanByAGError(offCloudStars)

        #is the folling code necessary?
        #onCloudStars = self.cleanByAGError(onCloudStars)

        # usually foregroundFITS only apply to off cloud stars, but for clouds hat are fore, it applies to both on and off cloud stars for far molecular clouds

        # assign the star with foreground fits

        if upCO != None:
            onCloudStars = onCloudStars[onCloudStars["coint"] <= upCO]

        print " --------- {}: on and off cloud stars are {} and {} ".format(cloudID, len(onCloudStars), len(offCloudStars))
        # save on cloud stars
        if useAV:
            onCloudStars.write(self.tmpPath + str(cloudID) + "OnCloudAV.fit", overwrite=True)
            offCloudStars.write(self.tmpPath + str(cloudID) + "OffCloudAV.fit", overwrite=True)

        else:
            onCloudStars.write(self.tmpPath + str(cloudID) + "OnCloudAG.fit", overwrite=True)
            offCloudStars.write(self.tmpPath + str(cloudID) + "OffCloudAG.fit", overwrite=True)

        #
        onCloudStarsRaw = onCloudStars.copy()
        offCloudStarsRaw = offCloudStars.copy()

        if useBin:  # test use binned data

            rowCloud = self.getRowByID(ID)

            area = rowCloud["area_exact"].to(u.deg ** 2)

            onBinTB, onBinInfo = self.bin2D(ID, onCloudStars, nPerPC10=onBinSize, position="On")
            # offBin=self.bin2D(ID, offCloudStars , nPerPC10=3 ,position="off")

            # onBin=self.bin2D(ID, onCloudStars , area=area ,position="On")
            # offBin=self.bin2D(ID, offCloudStars , area=area*2 ,position="off")
            offBinTB, offBinInfo = self.bin2D(ID, offCloudStars, nPerPC10=offBinSize, position="off")

            bestPair, disEachOther, newoffBin = self.getBestPair(onBinInfo, offBinInfo,
                                                                 maximumOffVolume=1.)  # np.unravel_index(np.argmin(disEachOther, axis=None), disEachOther.shape)

            # bestPair only cares about newoffBin

            if not useAllBins and not testAllBins:  # use the best one
                bestOffIndex = newoffBin[bestPair[1]][binTB.binID]
                onCloudStars = onCloudStars[onBinTB[self.binCol] == bestPair[0]]
                offCloudStars = offCloudStars[offBinTB[self.binCol] == bestOffIndex]

            if useAllBins and not testAllBins:  # use all the stars, all bins has to use different baseline to correct
                collectiveOn = Table(onCloudStars[0])
                collectiveOn.remove_row(0)

                collectiveOnRaw = Table(onCloudStars[0])
                collectiveOnRaw.remove_row(0)

                collectiveOff = Table(offCloudStars[0])
                collectiveOff.remove_row(0)

                minimumPair = np.argmin(disEachOther, axis=1)
                minimumV = np.min(disEachOther, axis=1)

                correctedBins = 0

                for i in range(len(onBinInfo)):

                    if minimumV[
                        i] > calibrateAngle:  # for large pairs, we reject them, because they cannot be correctly calibrated
                        print "No off-cloud stars around, reject this bin"
                        continue

                    correctedBins = correctedBins + 1
                    onCloudPart = onCloudStars[onBinTB[self.binCol] == i]

                    offBinID = newoffBin[minimumPair[i]][binTB.binID]

                    offCloudPart = offCloudStars[offBinTB[self.binCol] == offBinID]

                    dataDisTT, dataAGTT, disErrorTT, AGErrorTT = gaiaDis.getDisAndAGFromTB(onCloudPart)
                    baselineTMP = gaiaDis.getBaseLine(offCloudPart)
                    agBaseTT = gaiaDis.getAgBase(dataDisTT, baselineTMP)

                    collectiveOnRaw = vstack([collectiveOnRaw, onCloudPart.copy()])
                    if useAV:
                        onCloudPart["av50"] = onCloudPart["av50"] - agBaseTT
                    else:
                        onCloudPart["a_g_val"] = onCloudPart["a_g_val"] - agBaseTT

                    #

                    collectiveOn = vstack([collectiveOn, onCloudPart])
                    collectiveOff = vstack([collectiveOff, offCloudPart])

                # aa
                onCloudStars = collectiveOn
                offCloudStars = collectiveOff  # offCloudStars[   offBin["tb"][self.binCol] == bestPair[1]  ]

                if correctedBins == 1:  # if only one  on bins has off bins nearby, use this

                    onCloudStars = collectiveOnRaw

                    offCloudStars = collectiveOff  # #offCloudStars[   offBin["tb"][self.binCol] == bestPair[1]  ]

                    useBin = False
                    testAllBins = False

                    useAllBins = False

                if correctedBins == 0:
                    print "All off bins are far.. take the closest one..{} deg".format(np.min(disEachOther))

                    bestOffIndex = newoffBin[bestPair[1]][binTB.binID]
                    onCloudStars = onCloudStarsRaw[onBinTB[self.binCol] == bestPair[0]]
                    offCloudStars = offCloudStarsRaw[offBinTB[self.binCol] == bestOffIndex]

                    useBin = False
                    testAllBins = False

                    useAllBins = False
        else:

            testAllBins = False

            useAllBins = False

        if useBin and testAllBins:  # examine bins one by one

            minimumPair = np.argmin(disEachOther, axis=1)
            minimumV = np.min(disEachOther, axis=1)

            for i in range(len(onBinInfo)):

                if minimumV[
                    i] > calibrateAngle:  # for large pairs, we reject them, because they cannot be correctly calibrated
                    print "No off-cloud stars around, reject this bin"
                    continue

                offBinID = newoffBin[minimumPair[i]][binTB.binID]
                onCloudPart = onCloudStars[onBinTB[self.binCol] == i]
                offCloudPart = offCloudStars[offBinTB[self.binCol] == offBinID]

                if pureTest:
                    # os.remove( "saveForTestOn.fit" )
                    # onCloudStars.write("saveForTestOn.fit")
                    self.testAVOn(onCloudPart, offCloudPart, ID)

                    return

                baseLine = gaiaDis.getBaseLine(offCloudPart)

                saveFig = self.distFigurePath + figNameMark + "Distance{}_Bin{}.pdf".format(self.getAVAG(useAV), i)

                disResult = gaiaDis.calDisAndrawWithOnAndOffCloudStars(cloudID, intFITS, maskFITS, onCloudPart,
                                                                       offCloudPart, baseLine, searchlRange, searchbRange,
                                                                       saveFig, \
                                                                       noiseLevel=noiseLevel, signalLevel=signalLevel,
                                                                       draw=True, figNameMark=figNameMark,
                                                                       inputRow=newDisRow, correctBaseline=testAllBins,vlsr=v)

                self.doDis.updateRow(newDisRow)

            return  # stop here

        # print neighbourhood.argmin(axis=1) # for each on cloud bin

        # print neighbourhood.argmin(axis=0)

        # print "Tesing small number of on-cloud stars......"

        # onCloudStars=gaiaDis.getRandomRows(onCloudStars,50)

        # print " --------- {}: After selection on and off cloud stars are {} and {} ".format(cloudID, len(onCloudStars), len(offCloudStars)  )

        # get baseline
        baseLine = gaiaDis.getBaseLine(offCloudStars)
        onCloudStars = onCloudStars[onCloudStars[doAG.GAIA_distance] <= cutDis]
        offCloudStars = offCloudStars[offCloudStars[doAG.GAIA_distance] <= cutDis]
        offCloudStars = offCloudStars[   offCloudStars["distance_err"] > 0]  # distance_errr cannot be zeoro, possibly fore extremelly close stars

        # save onCloudStars,

        if pureTest:
            # os.remove( "saveForTestOn.fit" )
            # onCloudStars.write("saveForTestOn.fit")
            self.testAVOn(onCloudStars, offCloudStars, ID)

            return

        if useBaseline:

            saveFig = self.distFigurePath + figNameMark + "Distance{}{}.pdf".format(self.getAVAG(useAV), saveTag)

            disResult = gaiaDis.calDisAndrawWithOnAndOffCloudStars(cloudID, intFITS, maskFITS, onCloudStars, offCloudStars,
                                                                   baseLine, searchlRange, searchbRange, saveFig, \
                                                                   noiseLevel=noiseLevel, signalLevel=signalLevel,
                                                                   draw=True, figNameMark=figNameMark, inputRow=newDisRow,
                                                                   correctBaseline=not useAllBins ,vlsr=v,  legendCol=legendCol  )

            self.doDis.updateRow(newDisRow)
            #### reurn dis and error
            return disResult

        else:

            saveFig = self.distFigurePath + figNameMark + "Distance{}Normal.pdf".format(self.getAVAG(useAV))

            gaiaDis.calDisWithRowNormal(newDisRow, maskFITS=maskFITS, foregroundCOFITS=None, saveFigureName=saveFig)

            self.doDis.updateRow(newDisRow)

            return


    def getIntMaskByID(self, ID, pureInt=False):
        """
        """

        cloudID = str(ID)

        if "Cloud" not in cloudID:
            cloudID = "Cloud{}".format(cloudID)

        # searchIntFITS=  glob.glob(self.maskPath+"*{}_level*int.fits".format(cloudID))
        # searchMaskFITS=  glob.glob(self.maskPath+"*{}_level*mask.fits".format(cloudID))
        if pureInt:
            searchIntFITS = glob.glob( os.path.join( self.maskPath ,  "*{}_pureInt.fits".format(cloudID)) )

        else:

            searchIntFITS = glob.glob( os.path.join( self.maskPath , "*{}_*int.fits".format(cloudID)) )

        searchMaskFITS = glob.glob(  os.path.join( self.maskPath,  "*{}_*mask.fits".format(cloudID)) )


        return [searchIntFITS[0], searchMaskFITS[0]]

    def getLBrangeWithMask(self, maskFITS):
        """
        """

        # get the least region box, that contains the mask fits

        data, head = myFITS.readFITS(maskFITS)

        axisX = np.sum(data, axis=0)
        axisY = np.sum(data, axis=1)

        nonZerosX = np.nonzero(axisX)[0]

        firstX = nonZerosX[0]
        secondX = nonZerosX[-1]

        nonZerosY = np.nonzero(axisY)[0]

        firstY = nonZerosY[0]
        secondY = nonZerosY[-1]

        tempWCS = WCS(head, naxis=(1, 2))

        endL, startB = tempWCS.wcs_pix2world(firstX, firstY, 0)

        startL, endB = tempWCS.wcs_pix2world(secondX, secondY, 0)

        lRange = [min([endL, startL]), max([endL, startL])]

        bRange = [min([startB, endB]), max([startB, endB])]

        return lRange, bRange

    def expandLBRangeByDeg(self, lRange, bRange, backData, backWCS, lExpand=0.5, bExpand=0.25):
        """
        ###use absolute expand, to controle the size of box,
        :param lRange:
        :param bRange:
        :param backData:
        :param backWCS:
        :param ratioEx:
        :return:
        """

        sizeY, sizeX = backData.shape
        maxL, minB = backWCS.wcs_pix2world(0, 0, 0)
        minL, maxB = backWCS.wcs_pix2world(sizeX - 1, sizeY - 1, 0)

        newSmallL = max(min(lRange) - lExpand, minL)

        newSmallB = max(min(bRange) - bExpand, minB)

        newLargeL = min(max(lRange) + lExpand, maxL)

        newLargeB = min(max(bRange) + bExpand, maxB)

        return [newSmallL, newLargeL], [newSmallB, newLargeB]

    def expandLBRange(self, lRange, bRange, backData, backWCS, ratioEx=1):
        """
        expanding the ratio of LB range by 0.5
        """

        centerL = np.mean(lRange)
        centerB = np.mean(bRange)

        sizeL = max(lRange) - min(lRange)

        sizeB = max(bRange) - min(bRange)

        newsizeL = sizeL * (1 + ratioEx)
        newsizeB = sizeB * (1 + ratioEx)
        sizeY, sizeX = backData.shape

        # surpress the range in the fits

        maxL, minB = backWCS.wcs_pix2world(0, 0, 0)
        minL, maxB = backWCS.wcs_pix2world(sizeX - 1, sizeY - 1, 0)

        cutLLeft = min([centerL + newsizeL / 2., maxL])

        cutLRight = max([centerL - newsizeL / 2., minL])

        cutBBottome = max([centerB - newsizeB / 2., minB])

        cutBTop = min([centerB + newsizeB / 2., maxB])

        return [cutLRight, cutLLeft], [cutBBottome, cutBTop]

    def getCloudNameByLB(self, l, b):

        # if b>=0:

        lStr = str(l)

        bStr = "{:+f}".format(b)

        if '.' in lStr:

            lStrPart1, lStrPart2 = lStr.split('.')

        else:
            lStrPart1 = lStr
            lStrPart2 = '0'

        if '.' in bStr:

            bStrPart1, bStrPart2 = bStr.split('.')
        else:
            bStrPart1 = bStr
            bStrPart2 = '0'

        lStr = lStrPart1 + '.' + lStrPart2[0:1]

        bStr = bStrPart1 + '.' + bStrPart2[0:1]

        lStr = lStr.zfill(5)

        # bStr="{:+.1f}".format(b)
        bStrNumberPart = bStr[1:]
        bStr = bStr[0:1] + bStrNumberPart.zfill(4)

        cName = "G{}{}".format(lStr, bStr)

        return cName

    def cleanByAGError(self, TB):
        """

        remove those stars, that have too accurate AG measures, which is outliers

        """
        newTB = TB.copy()
        dataAG = newTB[doAG.GAIA_a_g_val]
        AGError = newTB["agError"] * dataAG

        newTB = newTB[AGError > 0.05]

        return newTB

    def getAVAG(self, useAV):

        if useAV:
            return "AV"

        return "AG"

    def ZZZ(self):
        pass