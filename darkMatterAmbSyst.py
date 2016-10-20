import ROOT as r
from   ROOT import gROOT, TCanvas, TFile, TF1, TH1F, TGraph
import math, sys, optparse, copy, re, array, itertools, os
import darkmatterSys
from multiprocessing import Pool
import include.CutManager as CutManager


def analyse(sys):
    print 'running class', sys
    analysis = darkmatterSys.darkmatterSys('theAnalysis',sys)
    print 'class loaded', sys
    analysis.GetCorrectionFactors()
    print 'we have the correction factors', sys
    analysis.GetDataDriven()
    resultHistos = []; dataDriven = []
    resultHistos.append( copy.deepcopy(analysis.tt    ))
    resultHistos.append( copy.deepcopy(analysis.ttz   ))
    resultHistos.append( copy.deepcopy(analysis.other ))
    resultHistos.append( copy.deepcopy(analysis.sig   ))
    resultHistos.append( copy.deepcopy(analysis.sig300))
    dataDriven  .append( analysis.TTScaleFactor        )
    dataDriven  .append( analysis.TTZScaleFactor       )

    print 'done', sys
    del analysis
    return (sys,resultHistos, dataDriven)



class darkMatterAmbSyst:
    def __init__(self, name):
        self.name = name


    def loadAll(self):
        print 'loading all'
        syst = ['','jecUp','jecDn','jecUp','bHeDn','bLiUp','bLiDn', 'ElUp', 'ElDn', 'MuUp', 'MuDn']
        print 'calling pool'

        # analysis = darkmatterSys.darkmatterSys('theAnalysis','bHeDn')
        # analysis.GetCorrectionFactors()
        # analysis.GetDataDriven()

        pool = Pool(3)
        print 'starting'


        
        analyses    = pool.map(analyse,syst)

        for sys,anal,dd in analyses:
            print sys, anal[0].Integral(), dd[0], dd[1]
            



if __name__ == '__main__':
    a = darkMatterAmbSyst('bla')
    a.loadAll()
