import math, sys
import ROOT as r
from   ROOT import TGraphErrors, gROOT, TCanvas, TFile

def selectSamples(inputfile, selList, sType = 'DATA'):
    f = open(inputfile, 'r')
    tmp_file = open('.tmp_sampleFile%s.txt' %sType, 'w')
    checkedList = []
    typeList    = []
    for line in f.readlines():
        if '#' in line or not len(line.rstrip('\r')): continue
        for _sample in selList:
            if _sample == line.split()[2]:
                if not sType == 'SYNCH':
                    tmp_file.write(line)
                else:
                    tmp_splitline = line.split()
                    tmp_splitline[0] = 'synching'
                    tmp_file.write('  '.join(tmp_splitline+['\n']))
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

class color:
       purple = '\033[95m'
       cyan = '\033[96m'
       darkcyan = '\033[36m'
       blue = '\033[94m'
       green = '\033[92m'
       yellow = '\033[93m'
       red = '\033[91m'
       bold = '\033[1m'
       underline = '\033[4m'
       end = '\033[0m'

class valErrs:
    def __init__(self, cen_val, cen_sys, cen_stat, fwd_val, fwd_sys, fwd_stat, name):
        self.cen_val  = cen_val
        self.cen_sys  = cen_sys
        self.cen_stat = cen_stat
        self.fwd_val  = fwd_val
        self.fwd_sys  = fwd_sys
        self.fwd_stat = fwd_stat
        self.vals = []
        self.name = name
    def totError(self):
        self.cen_err = math.sqrt((self.cen_val*self.cen_sys)**2 + self.cen_stat**2)
        self.fwd_err = math.sqrt((self.fwd_val*self.fwd_sys)**2 + self.fwd_stat**2)
        self.vals.append(self.cen_err)
        self.vals.append(self.fwd_err)
    def setVals(self, line):
        self.cen_val  = float(line.split()[-6])
        self.cen_sys  = float(line.split()[-4])
        self.cen_stat = float(line.split()[-5])
        self.vals.extend([self.cen_val, self.cen_sys, self.cen_stat])
        self.fwd_val  = float(line.split()[-3])
        self.fwd_sys  = float(line.split()[-1])
        self.fwd_stat = float(line.split()[-2])
        self.vals.extend([self.fwd_val, self.fwd_sys, self.fwd_stat])
        self.totError()
    def printValues(self):
        print '%s central %.3f +- %.3f (%.3f stat. %.3f syst.)' %(
                self.name, self.cen_val, self.cen_err, self.cen_stat, self.cen_sys)
        print '%s forward %.3f +- %.3f (%.3f stat. %.3f syst.)' %(
                self.name, self.fwd_val, self.fwd_err, self.fwd_stat, self.fwd_sys)

class ingredients:
    def __init__(self, infile, dType):
        self.infile = infile
        self.isData = (dType == 'DATA')
        self.dType  = dType
        self.rs = []
        ## rmue
        self.rmue_sr_lm     = valErrs(-1., -1., -1., -1., -1., -1., 'r_mue SR lm'    ); self.rs.append(self.rmue_sr_lm    )
        self.rmue_sr_onZ    = valErrs(-1., -1., -1., -1., -1., -1., 'r_mue SR onZ'   ); self.rs.append(self.rmue_sr_onZ   )
        self.rmue_sr_hm     = valErrs(-1., -1., -1., -1., -1., -1., 'r_mue SR hm'    ); self.rs.append(self.rmue_sr_hm    )
        self.rmue_dycr_dym  = valErrs(-1., -1., -1., -1., -1., -1., 'r_mue dycr dym' ); self.rs.append(self.rmue_dycr_dym )
        ## rsfof                                     
        self.rsfof_sr_lm    = valErrs(-1., -1., -1., -1., -1., -1., 'R_sfof SR lm'   ); self.rs.append(self.rsfof_sr_lm   )
        self.rsfof_sr_onZ   = valErrs(-1., -1., -1., -1., -1., -1., 'R_sfof SR onZ'  ); self.rs.append(self.rsfof_sr_onZ  )
        self.rsfof_sr_hm    = valErrs(-1., -1., -1., -1., -1., -1., 'R_sfof SR hm'   ); self.rs.append(self.rsfof_sr_hm   )
        self.rsfof_ttcr_lm  = valErrs(-1., -1., -1., -1., -1., -1., 'R_sfof TTCR lm' ); self.rs.append(self.rsfof_ttcr_lm )
        self.rsfof_ttcr_onZ = valErrs(-1., -1., -1., -1., -1., -1., 'R_sfof TTCR onZ'); self.rs.append(self.rsfof_ttcr_onZ)
        self.rsfof_ttcr_hm  = valErrs(-1., -1., -1., -1., -1., -1., 'R_sfof TTCR hm' ); self.rs.append(self.rsfof_ttcr_hm )
        ## rinout                                    
        self.rinout_dy_lm    = valErrs(-1., -1., -1., -1., -1., -1., 'R_inout SR lm'   ); self.rs.append(self.rinout_dy_lm)
        self.rinout_dy_bz    = valErrs(-1., -1., -1., -1., -1., -1., 'R_inout SR oz'   ); self.rs.append(self.rinout_dy_bz)
        self.rinout_dy_oz    = valErrs(-1., -1., -1., -1., -1., -1., 'R_inout SR bz'   ); self.rs.append(self.rinout_dy_oz)
        self.rinout_dy_az    = valErrs(-1., -1., -1., -1., -1., -1., 'R_inout SR az'   ); self.rs.append(self.rinout_dy_az)
        self.rinout_dy_hm    = valErrs(-1., -1., -1., -1., -1., -1., 'R_inout SR hm'   ); self.rs.append(self.rinout_dy_hm)
        ## RT
        self.rt_region      = valErrs(-1., -1., -1., -1., -1., -1., 'R_T region'     ); self.rs.append(self.rt_region     )
        ## fill the values from the file
        self.readValues()
        ## check if none is unset. exit if one is.
        self.checkValues()
        print 'loaded all ingredients from %s for %s' %(self.infile, self.dType)

    def check(self, test, string):
        return all(i in string for i in test)

    def readValues(self):
        print 'reading values from %s for %s' %(self.infile, self.dType)
        f = open(self.infile, 'r')
        lines = f.read().splitlines()
        for line in lines:
            if '#' in line or not len(line.strip()): continue
            #rmue
            if self.check(['rmue' , 'sr_lm'   , self.dType], line): self.rmue_sr_lm    .setVals(line)
            if self.check(['rmue' , 'sr_onZ'  , self.dType], line): self.rmue_sr_onZ   .setVals(line)
            if self.check(['rmue' , 'sr_hm'   , self.dType], line): self.rmue_sr_hm    .setVals(line)
            if self.check(['rmue' , 'dycr_dym', self.dType], line): self.rmue_dycr_dym .setVals(line)
            #rsfof
            if self.check(['rsfof', 'sr_lm'   , self.dType], line): self.rsfof_sr_lm   .setVals(line)
            if self.check(['rsfof', 'sr_onZ'  , self.dType], line): self.rsfof_sr_onZ  .setVals(line)
            if self.check(['rsfof', 'sr_hm'   , self.dType], line): self.rsfof_sr_hm   .setVals(line)
            if self.check(['rsfof', 'ttcr_lm' , self.dType], line): self.rsfof_ttcr_lm .setVals(line)
            if self.check(['rsfof', 'ttcr_onZ', self.dType], line): self.rsfof_ttcr_onZ.setVals(line)
            if self.check(['rsfof', 'ttcr_hm' , self.dType], line): self.rsfof_ttcr_hm .setVals(line)
            #rinout
            if self.check(['rinout', 'dy_lm' , self.dType], line): self.rinout_dy_lm   .setVals(line)
            if self.check(['rinout', 'dy_bz' , self.dType], line): self.rinout_dy_bz   .setVals(line)
            if self.check(['rinout', 'dy_oz' , self.dType], line): self.rinout_dy_oz   .setVals(line)
            if self.check(['rinout', 'dy_az' , self.dType], line): self.rinout_dy_az   .setVals(line)
            if self.check(['rinout', 'dy_hm' , self.dType], line): self.rinout_dy_hm   .setVals(line)
            #RT
            if self.check(['rt'   , 'region'  , self.dType], line): self.rt_region     .setVals(line)
        f.close()
    def checkValues(self):
        for thing in self.rs:
            if any(i < 0 for i in thing.vals):
                print 'ERROR: some of the ingredients aren\'t set properly'
                print thing.printValues()
                #sys.exit('exiting...')

