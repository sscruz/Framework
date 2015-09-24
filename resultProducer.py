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

def check(test, string):
    return all(i in string for i in test)

class valErrs:
    def __init__(self, val, sys, stat, name):
        self.val  = val
        self.sys  = sys
        self.stat = stat
        self.vals = []
        self.name = name
    def totError(self):
        self.err = math.sqrt(self.sys**2 + self.stat**2)
        self.vals.append(self.err)
    def setVals(self, line):
        self.val  = float(line.split()[-3])
        self.sys  = float(line.split()[-2])
        self.stat = float(line.split()[-1])
        self.vals.extend([self.val, self.sys, self.stat])
        self.totError()
    def printValues(self):
        print '%s %.3f +- %.3f (%.3f stat. %.3f syst.)' %(
                self.name, self.val, self.err, self.stat, self.sys)

class ingredients:
    def __init__(self, infile, isData):
        self.infile = infile
        self.isData = isData
        self.dataMC = 'DATA' if self.isData else 'MC'
        self.rs = []
        ## rmue
        self.rmue_sr_lm     = valErrs(-1., -1., -1., 'r_mue SR lm'    ); self.rs.append(self.rmue_sr_lm    )
        self.rmue_sr_onZ    = valErrs(-1., -1., -1., 'r_mue SR onZ'   ); self.rs.append(self.rmue_sr_onZ   )
        self.rmue_sr_hm     = valErrs(-1., -1., -1., 'r_mue SR hm'    ); self.rs.append(self.rmue_sr_hm    )
        self.rmue_dycr_dym  = valErrs(-1., -1., -1., 'r_mue dycr dym' ); self.rs.append(self.rmue_dycr_dym )
        ## rsfof
        self.rsfof_sr_lm    = valErrs(-1., -1., -1., 'R_sfof SR lm'   ); self.rs.append(self.rsfof_sr_lm   )
        self.rsfof_sr_onZ   = valErrs(-1., -1., -1., 'R_sfof SR onZ'  ); self.rs.append(self.rsfof_sr_onZ  )
        self.rsfof_sr_hm    = valErrs(-1., -1., -1., 'R_sfof SR hm'   ); self.rs.append(self.rsfof_sr_hm   )
        self.rsfof_ttcr_lm  = valErrs(-1., -1., -1., 'R_sfof TTCR lm' ); self.rs.append(self.rsfof_ttcr_lm )
        self.rsfof_ttcr_onZ = valErrs(-1., -1., -1., 'R_sfof TTCR onZ'); self.rs.append(self.rsfof_ttcr_onZ)
        self.rsfof_ttcr_hm  = valErrs(-1., -1., -1., 'R_sfof TTCR hm' ); self.rs.append(self.rsfof_ttcr_hm )
        ## RT
        self.rt_region      = valErrs(-1., -1., -1., 'R_T region'     ); self.rs.append(self.rt_region     )
        ## fill the values from the file
        self.readValues()
        ## check if none is unset. exit if one is.
        self.checkValues()
        print 'loaded all ingredients from %s for %s' %(self.infile, self.dataMC)

    def readValues(self):
        print 'reading values from %s for %s' %(self.infile, self.dataMC)
        f = open(self.infile, 'r')
        lines = f.read().splitlines()
        for line in lines:
            if '#' in line or not len(line.strip()): continue
            dataMC = 'DATA' if self.isData else 'MC'
            #rmue
            if check(['rmue' , dataMC, 'sr_lm'   ], line): self.rmue_sr_lm    .setVals(line)
            if check(['rmue' , dataMC, 'sr_onZ'  ], line): self.rmue_sr_onZ   .setVals(line)
            if check(['rmue' , dataMC, 'sr_hm'   ], line): self.rmue_sr_hm    .setVals(line)
            if check(['rmue' , dataMC, 'dycr_dym'], line): self.rmue_dycr_dym .setVals(line)
            #rsfof
            if check(['rsfof', dataMC, 'sr_lm'   ], line): self.rsfof_sr_lm   .setVals(line)
            if check(['rsfof', dataMC, 'sr_onZ'  ], line): self.rsfof_sr_onZ  .setVals(line)
            if check(['rsfof', dataMC, 'sr_hm'   ], line): self.rsfof_sr_hm   .setVals(line)
            if check(['rsfof', dataMC, 'ttcr_lm' ], line): self.rsfof_ttcr_lm .setVals(line)
            if check(['rsfof', dataMC, 'ttcr_onZ'], line): self.rsfof_ttcr_onZ.setVals(line)
            if check(['rsfof', dataMC, 'ttcr_hm' ], line): self.rsfof_ttcr_hm .setVals(line)
            #RT
            if check(['rt'   , dataMC, 'region'  ], line): self.rt_region     .setVals(line)
        f.close()
    def checkValues(self):
        for thing in self.rs:
            if any(i < 0 for i in thing.vals):
                print 'ERROR: some of the ingredients aren\'t set properly'
                print thing.printValues()
                #sys.exit('exiting...')



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

    rsMC = ingredients(ingredientFile, False)
    rsDA = ingredients(ingredientFile, True )
        
