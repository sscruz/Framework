import math, sys, os
import ROOT as r
from   ROOT import TGraphErrors, gROOT, TCanvas, TFile

def ensureDirectory(_path):
   #d = os.path.dirname(_path)
   #print d
   if not os.path.exists(_path):
      os.makedirs(_path)

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
            print 'check this sample:', _selSample
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
    def __init__(self, val, sys, stat, name):
        self.val  = val
        self.sys  = sys
        self.stat = stat
        self.vals = []
        self.name = name
    def totError(self):
        print 'total error', self.name
        self.err = math.sqrt((self.val*self.sys)**2 + self.stat**2)
        self.vals.append(self.err)
    def setVals(self, line):
        self.val  = float(line.split()[-3])
        self.sys  = float(line.split()[-1])
        self.stat = float(line.split()[-2])
        self.vals.extend([self.val, self.sys, self.stat])
        self.totError()
    def printValues(self):
        print '%s central %.3f +- %.3f (%.3f stat. %.3f syst.)' %(
                self.name, self.val, self.err, self.stat, self.sys)


class ingredients:
    def __init__(self, infile, dType):
        self.infile = infile
        self.isData = (dType == 'DATA')
        self.dType  = dType
        self.rs = []
        ## rmue
        self.rmue_alone   = valErrs(-1., -1., -1., 'r_mue alone'  ); self.rs.append(self.rmue_alone  )
        self.rmue_factor  = valErrs(-1., -1., -1., 'r_mue factor' ); self.rs.append(self.rmue_factor )
        ## rsfof                                   
        self.rsfof_direct = valErrs(-1., -1., -1., 'R_sfof direct'); self.rs.append(self.rsfof_direct)
        self.rsfof_final  = valErrs(-1., -1., -1., 'R_sfof final '); self.rs.append(self.rsfof_final )
        ## RT
        self.rt_region    = valErrs(-1., -1., -1., 'R_T'          ); self.rs.append(self.rt_region   )
        ## rinout                                  
        self.rinout_dy_lm = valErrs(-1., -1., -1., 'R_inout SR lm'); self.rs.append(self.rinout_dy_lm)
        self.rinout_dy_bz = valErrs(-1., -1., -1., 'R_inout SR oz'); self.rs.append(self.rinout_dy_bz)
        self.rinout_dy_oz = valErrs(-1., -1., -1., 'R_inout SR bz'); self.rs.append(self.rinout_dy_oz)
        self.rinout_dy_az = valErrs(-1., -1., -1., 'R_inout SR az'); self.rs.append(self.rinout_dy_az)
        self.rinout_dy_hm = valErrs(-1., -1., -1., 'R_inout SR hm'); self.rs.append(self.rinout_dy_hm)
        ## fill the values from the file
        self.readValues()
        ## check if none is unset. exit if one is.
        self.checkValues()
        ## caluclate final rsfof as well
        #self.calculateFullRSFOF()
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
            if self.check(['rmue'  , 'alone' , self.dType], line): self.rmue_alone  .setVals(line)
            print 'setting rmue factor'
            if self.check(['rmue'  , 'factor', self.dType], line): self.rmue_factor .setVals(line)
            #rsfof
            if self.check(['rsfof' , 'direct', self.dType], line): self.rsfof_direct.setVals(line)
            if self.check(['rsfof' , 'final' , self.dType], line): self.rsfof_final .setVals(line)
            #RT
            if self.check(['rt'    , 'region', self.dType], line): self.rt_region   .setVals(line)
            #rinout
            if self.check(['rinout', 'dy_lm' , self.dType], line): self.rinout_dy_lm.setVals(line)
            if self.check(['rinout', 'dy_bz' , self.dType], line): self.rinout_dy_bz.setVals(line)
            if self.check(['rinout', 'dy_oz' , self.dType], line): self.rinout_dy_oz.setVals(line)
            if self.check(['rinout', 'dy_az' , self.dType], line): self.rinout_dy_az.setVals(line)
            if self.check(['rinout', 'dy_hm' , self.dType], line): self.rinout_dy_hm.setVals(line)
        f.close()

    def calculateFullRSFOF(self):
        self.rsfof_fac   = self.rmue_factor.val * self.rt_region.val
        self.rsfof_fac_e = math.sqrt(self.rmue_factor.err**2 + self.rt_region.err**2)
        self.rsfof_fac_w = 1./(self.rsfof_fac_e**2)

        self.rsfof_dir    = self.rsfof_direct.val
        self.rsfof_dir_e  = self.rsfof_direct.err
        self.rsfof_dir_w  = 1./(self.rsfof_dir_e**2)

        ## calculate the numerator and denominator for the weighted average
        self.rsfof_final_num = self.rsfof_fac * self.rsfof_fac_w + self.rsfof_dir * self.rsfof_dir_w
        self.rsfof_final_den = self.rsfof_fac_w + self.rsfof_dir_w
        ## here come the final value and its uncertainty
        self.rsfof_final_val = self.rsfof_final_num / self.rsfof_final_den
        self.rsfof_final_err = 1./math.sqrt(self.rsfof_final_den)


        f=open(self.infile, 'r')
        newlines = []
        for line in list(f):
            if len(line.split()) > 2 and line.split()[0] == 'rsfof'  and line.split()[1] == 'final' and line.split()[2] == self.dType:
                newlines.append('rsfof       final           %-4s        %.4f      %.4f      %.4f      %.4f      %.4f      %.4f\n'%(
                    str(self.dType), self.rsfof_final_cen_val, self.rsfof_final_cen_err, 0.0))
            else:
                newlines.append(line)
        f.close()
        g=open(self.infile, 'w')
        g.writelines(newlines)
        g.close()


    def checkValues(self):
        for thing in self.rs:
            if any(i < 0 for i in thing.vals):
                print 'ERROR: some of the ingredients aren\'t set properly'
                print thing.printValues()
                #sys.exit('exiting...')



##################################################################################
### These functions were created by Aachen in order to define the proper style ###
##################################################################################
def createMyColors():
    iIndex = 2000

    containerMyColors = []
    for color in defineMyColors.keys():
       	tempColor = r.TColor(iIndex,
       	float(defineMyColors[color][0]) / 255, float(defineMyColors[color][1]) / 255, float(defineMyColors[color][2]) / 255)
       	containerMyColors.append(tempColor)

       	myColors.update({ color: iIndex })
       	iIndex += 1

    return containerMyColors


defineMyColors = {
        'Black' : (0, 0, 0),
        'White' : (255, 255, 255),
        'Red' : (255, 0, 0),
        'DarkRed' : (128, 0, 0),
        'Green' : (0, 255, 0),
        'Blue' : (0, 0, 255),
        'Yellow' : (255, 255, 0),
        'Orange' : (255, 128, 0),
        'DarkOrange' : (255, 64, 0),
        'Magenta' : (255, 0, 255),
        'KDEBlue' : (64, 137, 210),
        'Grey' : (128, 128, 128),
        'DarkGreen' : (0, 128, 0),
        'DarkSlateBlue' : (72, 61, 139),
        'Brown' : (70, 35, 10),

        'MyBlue' : (36, 72, 206),
        'MyDarkBlue' : (18, 36, 103),
        'MyGreen' : (70, 164, 60),
        'AnnBlueTitle' : (29, 47, 126),
        'AnnBlue' : (55, 100, 255),
#        'W11AnnBlue' : (0, 68, 204),
#        'W11AnnBlue' : (63, 122, 240),
    }


myColors = {
            'W11ttbar':  855,
            'W11singlet':  854,
            'W11ZLightJets':  401,
            'W11ZbJets':  400,
            'W11WJets':  842,
            'W11Diboson':  920,
            'W11AnnBlue': 856,
            'W11Rare':  630,
            }


