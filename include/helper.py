import math, sys
import ROOT as r
from   ROOT import TGraphErrors, gROOT, TCanvas, TFile

def selectSamples(inputfile, selList, sType = 'DATA'):
    f = open(inputfile, 'r')
    tmp_file = open('.tmp_sampleFile%s.txt' %sType, 'w')
    checkedList = []
    typeList    = []
    for line in f.readlines():
        if '#' in line: continue
        for _sample in selList:
            if _sample == line.split()[2]:
                tmp_file.write(line)
                checkedList.append(_sample)
                typeList   .append(int(line.split()[-1]))
    for _selSample in selList:
        if _selSample not in checkedList:
            print 'ERROR: some samples weren\'t selected, check all sample names!'
            sys.exit('exiting...')
    if not len(set(typeList)) == 1:
            print 'ERROR: you\'re mixing DATA and MC!'
            sys.exit('exiting...')
            
    return tmp_file.name


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

    def check(self, test, string):
        return all(i in string for i in test)

    def readValues(self):
        print 'reading values from %s for %s' %(self.infile, self.dataMC)
        f = open(self.infile, 'r')
        lines = f.read().splitlines()
        for line in lines:
            if '#' in line or not len(line.strip()): continue
            dataMC = 'DATA' if self.isData else 'MC'
            #rmue
            if self.check(['rmue' , 'sr_lm'   , dataMC], line): self.rmue_sr_lm    .setVals(line)
            if self.check(['rmue' , 'sr_onZ'  , dataMC], line): self.rmue_sr_onZ   .setVals(line)
            if self.check(['rmue' , 'sr_hm'   , dataMC], line): self.rmue_sr_hm    .setVals(line)
            if self.check(['rmue' , 'dycr_dym', dataMC], line): self.rmue_dycr_dym .setVals(line)
            #rsfof
            if self.check(['rsfof', 'sr_lm'   , dataMC], line): self.rsfof_sr_lm   .setVals(line)
            if self.check(['rsfof', 'sr_onZ'  , dataMC], line): self.rsfof_sr_onZ  .setVals(line)
            if self.check(['rsfof', 'sr_hm'   , dataMC], line): self.rsfof_sr_hm   .setVals(line)
            if self.check(['rsfof', 'ttcr_lm' , dataMC], line): self.rsfof_ttcr_lm .setVals(line)
            if self.check(['rsfof', 'ttcr_onZ', dataMC], line): self.rsfof_ttcr_onZ.setVals(line)
            if self.check(['rsfof', 'ttcr_hm' , dataMC], line): self.rsfof_ttcr_hm .setVals(line)
            #RT
            if self.check(['rt'   , 'region'  , dataMC], line): self.rt_region     .setVals(line)
        f.close()
    def checkValues(self):
        for thing in self.rs:
            if any(i < 0 for i in thing.vals):
                print 'ERROR: some of the ingredients aren\'t set properly'
                print thing.printValues()
                #sys.exit('exiting...')

class rsfofRegion:
    def __init__(self, name, cenFwd, cuts, bins, isGraph, var, doData):
        self.name      = name
        self.cuts      = cuts
        self.bins      = bins
        self.isGraph   = isGraph
        self.cenFwd    = cenFwd
        self.isCentral = (cenFwd == 'central')
        self.doData    = doData

    def set_rsfof(self, rsfof, data=0):
        if not data:
            self.rsfof         = rsfof
            self.rsfof_gr      = r.TGraphErrors(rsfof)
            self.rsfof.GetYaxis().SetRangeUser(0, 2)
        else:
            self.rsfof_data    = rsfof
            self.rsfof_data_gr = r.TGraphErrors(rsfof)
            self.rsfof_data.GetYaxis().SetRangeUser(0, 2)

    def printValues(self):
        print 'REGION', self.name, self.cenFwd
        for _bin in range(1,self.rsfof.GetNbinsX()+1):
            print 'r_sfof in [%.0f, %.0f] in %s: %.3f +- %.3f' %(
                  self.rsfof.GetXaxis().GetBinLowEdge(_bin),
                  self.rsfof.GetXaxis().GetBinUpEdge (_bin),
                  self.cenFwd,
                  self.rsfof.GetBinContent(_bin),
                  self.rsfof.GetBinError(_bin))

class rmueRegion:
    def __init__(self, name, cenFwd, cuts, bins, isGraph, mllMet, doData):
        self.name      = name
        self.cuts      = cuts
        self.bins      = bins
        self.isGraph   = isGraph
        self.cenFwd    = cenFwd
        self.isCentral = (cenFwd == 'central')
        self.doMll     = ('mll' in mllMet)
        self.doMet     = ('met' in mllMet)
        self.doData    = doData

    def set_rmue_mll(self, rmue_mll, data=0):
        if not data:
            self.rmue_mll         = rmue_mll
            self.rmue_mll_gr      = TGraphErrors(rmue_mll)
            self.rmue_mll.GetYaxis().SetRangeUser(0, 2)
        else:
            self.rmue_mll_data    = rmue_mll
            self.rmue_mll_data_gr = TGraphErrors(rmue_mll)
            self.rmue_mll_data.GetYaxis().SetRangeUser(0, 2)

    def set_rmue_met(self, rmue_met, data=0):
        if not data:
            self.rmue_met         = rmue_met
            self.rmue_met_gr      = TGraphErrors(rmue_met)
            self.rmue_met.GetYaxis().SetRangeUser(0, 2)
        else:
            self.rmue_met_data    = rmue_met
            self.rmue_met_data_gr = TGraphErrors(rmue_met)
            self.rmue_met_data.GetYaxis().SetRangeUser(0, 2)


    def printValues(self):
        for i in ([1, 2] if self.doData else [1]):
            print 'REGION %s \t %s \t %s' %(self.name, self.cenFwd, 'DATA' if i>1 else 'MC')
            obj = self.rmue_mll_data if i>1 else self.rmue_mll
            for _bin in range(1,self.rmue_mll.GetNbinsX()+1):
                print 'r_mue in [%.0f, %.0f] in %s: %.3f +- %.3f' %(
                      obj.GetXaxis().GetBinLowEdge(_bin),
                      obj.GetXaxis().GetBinUpEdge (_bin),
                      self.cenFwd,
                      obj.GetBinContent(_bin),
                      obj.GetBinError(_bin))

    def saveInFile(self, ingFile, pattern, systErr):
        print 'writing calculated values into file...'
        f = open(ingFile, 'r')
        lines = f.readlines()
        newlines = []
        for line in lines:
            appended = False
            for t in ['MC', 'DATA'] if self.doData else ['MC']:
                if all(s in line for s in pattern+[t]):
                    obj = self.rmue_mll if t =='MC' else self.rmue_mll_data
                    newlines.append('%s \t %s \t %s \t %.4f \t %.4f \t %s \t %.4f \t %.4f \t %s\n' %(
                            line.split()[0], line.split()[1], t,
                            obj.GetBinContent(1) if     self.isCentral else float(line.split()[3]),
                            obj.GetBinError  (1) if     self.isCentral else float(line.split()[4]),
                            str(systErr)         if     self.isCentral else line.split()[5],
                            obj.GetBinContent(1) if not self.isCentral else float(line.split()[6]),
                            obj.GetBinError  (1) if not self.isCentral else float(line.split()[7]),
                            str(systErr)         if not self.isCentral else line.split()[8]))
                    appended = True

            if not appended:
                newlines.append(line)
        f.close()
        g = open(ingFile, 'w')
        g.writelines(newlines)
