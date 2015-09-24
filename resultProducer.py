#####################################################################
######                                                              #
###### 8=========D         OO                        8==========D   #  
###### OO                  ||                                ,88   #
###### ||                  ||                              ,88'    #  
###### ||O---\     ,adPPYb,|| ,adPPYb,d8  ,adPPYba,      ,88'      #
###### ||O---/    a8'    `Y||a8'    `Y88 a8P_____88    ,88'        #
###### ||         8b       ||8b       88 8PP'''''''  ,88'          #
###### \/         '8a,   ,d||'8a,   ,d88 '8b,   ,aa 88'            #
###### 8=========D `'8bbdP'\/  `'YbbdP'Y8  `'Ybbd8'' 8==========D   #
######                       aa,    ,88                             #
######                         'Y8bbdP'                             #
######                                                              #
#####################################################################

import ROOT as r
import math, sys, optparse
import Canvas, CutManager, Sample

from ROOT import gROOT, TCanvas, TFile, TF1, TPaveStats

import include.helper as helper


if __name__ == '__main__':

    print 'Starting r_SFOF analysis...'
    parser = optparse.OptionParser(usage='usage: %prog [options] FilenameWithSamples FilenameWithIngredients', version='%prog 1.0')
    parser.add_option('-m', '--mode', action='store', dest='mode', default='rsfof', help='Operation mode')
    (options, args) = parser.parse_args()


    if len(args) != 2:
      parser.error('wrong number of arguments')

    sampleFile     = args[0]
    ingredientFile = args[1]

    print 'Going to load DATA and MC trees...'
    mcDatasets = ['TTJets']#, 'DYJetsToLL_M10to50', 'DYJetsToLL_M50']
    daDatasets = ['DoubleMuon_Run2015C', 'DoubleEG_Run2015C', 'MuonEG_Run2015C']

    treeMC = Sample.Tree(Sample.selectSamples(sampleFile, mcDatasets, 'MC'), 'MC'  , 0)
    treeDA = Sample.Tree(Sample.selectSamples(sampleFile, daDatasets, 'DA'), 'DATA', 1)

    print 'Trees successfully loaded...'


    lumi = 0.020
    print 'Running with an integrated luminosity of', lumi,'fb-1'


   
    gROOT.ProcessLine('.L tdrstyle.C')
    gROOT.SetBatch(1)
    r.setTDRStyle() 
    cuts = CutManager.CutManager()


    doClosure = True
    isBlinded = True

    rsMC = helper.ingredients(ingredientFile, False)
    rsDA = helper.ingredients(ingredientFile, True )
        
