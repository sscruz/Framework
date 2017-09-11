import ROOT as r
from   ROOT import gROOT, TCanvas, TFile, TF1, TPaveStats
import math, sys, optparse, copy, re, array, os, pickle
import random, string

import include.helper     as helper
import include.Region     as Region
import include.Canvas     as Canvas
import include.CutManager as CutManager
import include.Sample     as Sample
import include.Tables     as Tables
import include.Scans      as Scans
import include.LeptonSF 
import include.FastSimSF
import include.nll

import subprocess

#import include.nll

from multiprocessing import Pool

global _r
_r = r.TRandom3(42)


def makeMCDatacards():
    print 'producing mc datacards'
    othersDatasets = [ 'GluGluToContinToZZTo2e2nu', 'GluGluToContinToZZTo2mu2nu', 'WZTo3LNu', 'WWG', 'WWW', 'WWZ',  'WZZ', 'ZZZ', 'ZGTo2LG',  'WWDouble', 'WpWpJJ',  'ZZTo2L2Nu', 'TTHnobb_pow', 'VHToNonbb', 'TWZ', 'WZTo2L2Q']
    ttDatasets = ['TTTo2L2Nu','WWTo2L2Nu', 'GGHWWTo2L2Nu', 'GGWWTo2L2Nu']
    dyDatasets = ['DYJetsToLL_M10to50_LO', 'DYJetsToLL_M50_LO']
    print 'getting trees'
    treeTT     = Sample.Tree(helper.selectSamples('samplesSlepton.dat', ttDatasets, 'TT'), 'TT', 0, isScan = 0)
    treeDY     = Sample.Tree(helper.selectSamples('samplesSlepton.dat', dyDatasets, 'DY'), 'DY', 0, isScan = 0)
    treeOTHERS = Sample.Tree(helper.selectSamples('samplesSlepton.dat', othersDatasets, 'OTHERS'), 'OTHERS', 0, isScan = 0)
    print 'getting yields'
    print cuts.AddList([scan.cuts_norm,cuts.goodLepton])
    print cuts
    tt = treeTT.getTH1F(lumi, 'FSbkg', scan.srID, scan.srIDMax+1, -0.5, scan.srIDMax+0.5, cuts.AddList([scan.cuts_norm,cuts.goodLepton]), '', '', "1", "noKFactor"); tt.SetName('FSbkg')
    dy = treeDY.getTH1F(lumi, 'dybkg', scan.srID, scan.srIDMax+1, -0.5, scan.srIDMax+0.5, cuts.AddList([scan.cuts_norm,cuts.goodLepton]), '', '', "1", "noKFactor"); dy.SetName('dybkg')
    others = treeOTHERS.getTH1F(lumi, 'otherbkg', scan.srID, scan.srIDMax+1, -0.5, scan.srIDMax+0.5, cuts.AddList([scan.cuts_norm,cuts.goodLepton]), '', '', "1", "noKFactor"); others.SetName('otherbkg')
    data = tt.Clone('data_obs'); data.Add(dy);data.Add(others)
    helper.ensureDirectory('datacards/datacards_%s/%s'%(scan.name,scan.name))
    for SR, label in scan.shortLabels.items():
        if scan.hasOther:
            datacard = '''imax 1 number of bins
jmax 3 number of processes minus 1
kmax *  number of nuisance parameters
----------------------------------------------------------------------------------------------------------------------------------
bin          {label}
observation  {obs}   
----------------------------------------------------------------------------------------------------------------------------------
bin          {label}      {label}        {label}      {label}
process      XXSIGNALXX     FSbkg          DY         Other
process      0             1              2           3
rate         XXSIGRATEXX    {fs}           {DY}       {other}
----------------------------------------------------------------------------------------------------------------------------------
fs_stat_{label} gmN  {fs_int}  -            1.0             -         - 
fs_unc lnN             -            1.05             -          -
jec   lnN           XXjecXX            -             -          -
El    lnN           XXElXX             -             -          -
Mu    lnN           XXMuXX             -             -          -
FastSimEl   lnN     XXFastSimElXX      -             -          -
FastSimMu   lnN     XXFastSimMuXX      -             -          -
PU    lnN             1.0              -             -          - 
SigTrig  lnN          1.05             -             -          -
bHe   lnN           XXbHeXX            -             -          -
bLi   lnN           XXbLiXX            -             -          -
genMet lnU              XXgenMetXX     -             -          - 
signalMCstats_{label} lnN    XXmcStatXX  -           -          -
dySys  lnN           -                  -            1.4        - 
otherSys lnN         -                  -            -          1.5
lumi   lnN             1.026              -            -          - 
'''.format(obs = int(tt.GetBinContent(tt.FindBin(SR))) + int(dy.GetBinContent(dy.FindBin(SR))) + int(others.GetBinContent(others.FindBin(SR))), 
           fs  = tt.GetBinContent(tt.FindBin(SR)), DY = dy.GetBinContent(dy.FindBin(SR)), other = others.GetBinContent(others.FindBin(SR)), 
           label = label, fs_int = int(tt.GetBinContent(tt.FindBin(SR))))

        else:
            datacard = '''imax 1 number of bins
jmax 2 number of processes minus 1
kmax *  number of nuisance parameters
----------------------------------------------------------------------------------------------------------------------------------
bin          {label}
observation  {obs}   
----------------------------------------------------------------------------------------------------------------------------------
bin          {label}      {label}        {label}
process      XXSIGNALXX     FSbkg          DY
process      0             1              2   
rate         XXSIGRATEXX    {fs}           {DY}
----------------------------------------------------------------------------------------------------------------------------------
fs_stat_{label} gmN  {fs_int}  -            1.0             - 
fs_unc lnN             -            1.05            - 
jec   lnN           XXjecXX            -               -
El    lnN           XXElXX            -               -
Mu    lnN           XXMuXX            -               -
FastSimEl   lnN     XXFastSimElXX     -              - 
FastSimMu   lnN     XXFastSimMuXX     -              -
PU    lnN           1.0            -             -          - 
SigTrig  lnN          1.05             -               -
bHe   lnN           XXbHeXX            -               -
bLi   lnN           XXbLiXX            -               -
genMet lnU              XXgenMetXX         -               - 
signalMCstats_{label} lnN    XXmcStatXX       -                -
dySys  lnN           -                  -            1.4
lumi   lnN             1.026            -               -
'''.format(obs = tt.GetBinContent(tt.FindBin(SR)) + dy.GetBinContent(dy.FindBin(SR)), 
           fs  = tt.GetBinContent(tt.FindBin(SR)), DY = dy.GetBinContent(dy.FindBin(SR)), 
           label = label, fs_int = int(tt.GetBinContent(tt.FindBin(SR))))
        outputFile = open('datacards/datacards_{scan}/{scan}/template_{sr}.txt'.format(scan=scan.name,sr=label),'w')
        outputFile.write(datacard)
        outputFile.close()

    # no shapes anymore :(
    # shapes = r.TFile('datacards/datacards_%s/shapes.root'%scan.name,'recreate'); 
    # tt.Write(); dy.Write(); data.Write();
    # shapes.Close()
    

def makeDataCardsFromRootFileSlepton():                                                                                                                                            
    print "Making datacard from rootfile"
    rootfile = TFile.Open('datacards/forDatacards_slepton_2017_.root')
    da_SF    = rootfile.Get('da_SF'   )
    da_OF    = rootfile.Get('da_OF'   )
    tf_CR_SR = rootfile.Get('tf_CR_SR')
    rareMC   = rootfile.Get('rare' )
    wz       = rootfile.Get('wz_SF' )
    zz       = rootfile.Get('zz_SF' )
    zz_jecUp = rootfile.Get('zz_jecUp')
    zz_jecDn = rootfile.Get('zz_jecDn')
    wz_jecUp = rootfile.Get('wz_jecUp')
    wz_jecDn = rootfile.Get('wz_jecDn')
    zz_qcdUp = rootfile.Get('zz_qcdUp')
    zz_qcdDn = rootfile.Get('zz_qcdDn')
    wz_qcdUp = rootfile.Get('wz_qcdUp')
    wz_qcdDn = rootfile.Get('wz_qcdDn')
    zz_pdf   = rootfile.Get('pdfZZ')
    wz_pdf   = rootfile.Get('pdfWZ')
    print "opened ", rootfile    
    #puFile = TFile.Open('datacards/PuSysts_{name}.root'.format(name ='CharNeu_Moriond2017' if scan.name=='ChiZZ_Moriond2017' else scan.name), 'read')
    #puSystUp = puFile.Get('puUp')
    #puSystDn = puFile.Get('puDn')
#fs_stat_{label} gmN  {fs_int}  -  {tf}                -                  -
#fs_unc  lnN         -                 {tf_e}          -                  -                  -
#fs   lnN         -                 {fs_int}          -                  -                  -

    for SR, label in scan.shortLabels.items():
        #print "SR", SR
        #print "zz_pdf.GetBinContent(SR+1) ", zz_pdf.GetBinContent(SR+1)
        #print "wz_pdf.GetBinContent(SR+1) ", wz_pdf.GetBinContent(SR+1)
        if scan.hasOther:
            datacard = '''imax 1 number of bins
jmax 4 number of processes minus 1
kmax *  number of nuisance parameters
----------------------------------------------------------------------------------------------------------------------------------
bin          {label}
observation  {obs}   
----------------------------------------------------------------------------------------------------------------------------------
bin                                {label}            {label}       {label}            {label}            {label}
process                            XXSIGNALXX         FSbkg         ZZ                 WZ                 rare
process                            0                  1             2                  2                  3
rate                               XXSIGRATEXX        {fs}  {ZZ}  {WZ}  {rare}
--------------------------------------------------------------------------------------------------------------------------------------
fs_stat_{label} gmN  {fs_int}  -  {tf}                -                  -                  -              - 
fs_unc  lnN         -                 {tf_e}          -                  -                  -
jec  lnN         XXjecXX            -             {zz_jUp}/{zz_jDn}                  {wz_jUp}/{wz_jDn}                  -
El  lnN          XXElXX             -             1.03               1.03               -
Mu  lnN          XXMuXX             -             1.06               1.06               -
SigTrig  lnN          1.05               -             -                  -                  -
Trig  lnN          -                  -             1.03               1.03               -
bHe  lnN          XXbHeXX            -             -                  -                  -
bLi  lnN          XXbLiXX            -             -                  -                  -
signalMCstats_{label} lnN          XXmcStatXX         -             -                  -                  -
zzSys  lnN          -                  -             1.1                -                  - 
wzSys  lnN          -                  -             -                  1.07               - 
rareSys  lnN          -                  -             -                  -                  1.5
scale  lnN          1.03               -             -                  -                  - 
pdf  lnN          1.03               -             {zz_p}      {wz_p}      - 
lumi  lnN          1.025              -             1.025              1.025              -               
'''.format(label = label, 
           obs = da_OF.GetBinContent(SR+1), # this is for blinding!!!!!!!!!!!!!!!
           #obs = da_SF.GetBinContent(SR+2), # the plots for ewk start in [50,100]
           fs  = da_OF.GetBinContent(SR+1)*tf_CR_SR.GetBinContent(SR+1),
           ZZ = zz.GetBinContent(SR+1),
           WZ = wz.GetBinContent(SR+1),
           rare = rareMC.GetBinContent(SR+1),
           fs_int = int(da_OF.GetBinContent(SR+1)),
           tf = tf_CR_SR.GetBinContent(SR+1),
           tf_e = 1 + tf_CR_SR.GetBinError(SR+1)/tf_CR_SR.GetBinContent(SR+1),
           zz_jUp = 1 + zz_jecUp.GetBinContent(SR+1),
           zz_jDn = 1 + zz_jecUp.GetBinContent(SR+1),
           wz_jUp = 1 + wz_jecUp.GetBinContent(SR+1),
           wz_jDn = 1 + wz_jecDn.GetBinContent(SR+1),
           #rare_e = 1 + rareMC.GetBinError(SR+2) / rareMC.GetBinContent(SR+2),
           zz_p = 1 + zz_pdf.GetBinContent(SR+1),
           wz_p = 1 + wz_pdf.GetBinContent(SR+1))
           #dy_e = 1.5 if dy_shape.GetBinContent(SR+2) == 0 else 1 + dy_shape.GetBinError(SR+2)/dy_shape.GetBinContent(SR+2),
           #pu_Up = 1 + puSystUp.GetBinContent(SR+1),
           #pu_Dn = 1 + puSystDn.GetBinContent(SR+1),
           #zz_e = 1 + ZZ.GetBinError(SR+2) / ZZ.GetBinContent(SR+2),
           #wz_e = 1 + WZ.GetBinError(SR+2) / WZ.GetBinContent(SR+2),
        print da_OF.GetBinContent(SR+1)
        #print da_SF.GetBinContent(SR+2)

        outputFile = open('datacards/datacards_{scan}/{scan}/datacard_{sr}.txt'.format(scan=scan.name,sr=label),'w')
        outputFile.write(datacard)
        outputFile.close()                                                                                                                                              


def makeDataCardsFromRootFile():
    print "Making datacard from rootfile"
    rootfile = TFile.Open('datacards/forDatacards_%s.root'%('CharNeu_Moriond2017' if scan.name=='ChiZZ_Moriond2017' else scan.name))
    da_SF    = rootfile.Get('da_SF'   )
    da_OF    = rootfile.Get('da_OF'   )
    tf_CR_SR = rootfile.Get('tf_CR_SR')
    dy_shape = rootfile.Get('dy_shape')
    mc_full  = rootfile.Get('mc_full' )
    print "opened ", rootfile    
    puFile = TFile.Open('datacards/PuSysts_{name}.root'.format(name ='CharNeu_Moriond2017' if scan.name=='ChiZZ_Moriond2017' else scan.name), 'read')

    puSystUp = puFile.Get('puUp')
    puSystDn = puFile.Get('puDn')

    for SR, label in scan.shortLabels.items():
        if scan.hasOther:
            datacard = '''imax 1 number of bins
jmax 3 number of processes minus 1
kmax *  number of nuisance parameters
----------------------------------------------------------------------------------------------------------------------------------
bin          {label}
observation  {obs}   
----------------------------------------------------------------------------------------------------------------------------------
bin          {label}      {label}        {label}      {label}
process      XXSIGNALXX     FSbkg          DY         Other
process      0             1              2           3
rate         XXSIGRATEXX    {fs}           {DY}       {other}
----------------------------------------------------------------------------------------------------------------------------------
fs_stat_{label} gmN  {fs_int}  -            {tf}             -         - 
fs_unc lnN             -            {tf_e}             -          -
jec   lnN           XXjecXX            -             -          -
El    lnN           XXElXX             -             -          -
Mu    lnN           XXMuXX             -             -          -
FastSimEl   lnN     XXFastSimElXX      -             -          -
FastSimMu   lnN     XXFastSimMuXX      -             -          -
SigTrig  lnN          1.05             -             -          -
bHe   lnN           XXbHeXX            -             -          -
bLi   lnN           XXbLiXX            -             -          -
PU    lnN           {pu_Up}/{pu_Dn}            -             -          - 
genMet lnU              XXgenMetXX     -             -          -
ISR   lnN              XXISRXX     -             -          - 
signalMCstats_{label} lnN    XXmcStatXX  -           -          -
dySys  lnN           -                  -            {dy_e}        - 
otherSys lnN         -                  -            -          {other_e}
scale    lnN         1.03               -            -            - 
lumi   lnN             1.026              -            -          - 
'''.format(label = label, 
           obs = da_SF.GetBinContent(SR+2), # the plots for ewk start in [50,100]
           fs  = da_OF.GetBinContent(SR+2)*tf_CR_SR.GetBinContent(SR+2),
           DY = dy_shape.GetBinContent(SR+2),
           other = mc_full.GetBinContent(SR+2),
           fs_int = int(da_OF.GetBinContent(SR+2)),
           tf = tf_CR_SR.GetBinContent(SR+2),
           tf_e = 1 + tf_CR_SR.GetBinError(SR+2)/tf_CR_SR.GetBinContent(SR+2),
           dy_e = 1.5 if dy_shape.GetBinContent(SR+2) == 0 else 1 + dy_shape.GetBinError(SR+2)/dy_shape.GetBinContent(SR+2),
           pu_Up = 1 + puSystUp.GetBinContent(SR+1),
           pu_Dn = 1 + puSystDn.GetBinContent(SR+1),
           other_e = 1 + mc_full.GetBinError(SR+2) / mc_full.GetBinContent(SR+2))
        print da_SF.GetBinContent(SR+2)

        outputFile = open('datacards/datacards_{scan}/{scan}/datacard_{sr}.txt'.format(scan=scan.name,sr=label),'w')
        outputFile.write(datacard)
        outputFile.close()                                                                                                                                              

def makeDataCardsFromRootFileForEdge():
    print "Making datacard from rootfile"
    rootfileBelow = TFile.Open('datacards/forDatacards_%s_%s.root'%(scan.name,'nllBelow21'))
    rootfileAbove = TFile.Open('datacards/forDatacards_%s_%s.root'%(scan.name,'nllAbove21'))

    da_SFBelow    = copy.deepcopy( rootfileBelow.Get('da_SF'   )); da_SFBelow   .SetName('da_SFBelow'   )
    da_OFBelow    = copy.deepcopy( rootfileBelow.Get('da_OF'   )); da_OFBelow   .SetName('da_OFBelow'   )
    tf_CR_SRBelow = copy.deepcopy( rootfileBelow.Get('tf_CR_SR')); tf_CR_SRBelow.SetName('tf_CR_SRBelow')
    dy_shapeBelow = copy.deepcopy( rootfileBelow.Get('dy_shape')); dy_shapeBelow.SetName('dy_shapeBelow')
    mc_fullBelow  = copy.deepcopy( rootfileBelow.Get('mc_full' )); mc_fullBelow .SetName('mc_fullBelow' )

    da_SFAbove    = copy.deepcopy( rootfileAbove.Get('da_SF'   )); da_SFAbove   .SetName('da_SFAbove'   )
    da_OFAbove    = copy.deepcopy( rootfileAbove.Get('da_OF'   )); da_OFAbove   .SetName('da_OFAbove'   )
    tf_CR_SRAbove = copy.deepcopy( rootfileAbove.Get('tf_CR_SR')); tf_CR_SRAbove.SetName('tf_CR_SRAbove')
    dy_shapeAbove = copy.deepcopy( rootfileAbove.Get('dy_shape')); dy_shapeAbove.SetName('dy_shapeAbove')
    mc_fullAbove  = copy.deepcopy( rootfileAbove.Get('mc_full' )); mc_fullAbove .SetName('mc_fullAbove' )

    puFile = TFile.Open('datacards/PuSysts_{name}.root'.format(name=scan.name), 'read')

    puSystUp = puFile.Get('puUp')
    puSystDn = puFile.Get('puDn')

    for SR, label in scan.shortLabels.items():
        bin = SR + 1 if SR < 2 else SR+2 if SR < 7 else SR-6 if SR < 9 else SR-5
        print SR, bin
        if scan.hasOther:
            datacard = '''imax 1 number of bins
jmax 3 number of processes minus 1
kmax *  number of nuisance parameters
----------------------------------------------------------------------------------------------------------------------------------
bin          {label}
observation  {obs}   
----------------------------------------------------------------------------------------------------------------------------------
bin          {label}      {label}        {label}      {label}
process      XXSIGNALXX     FSbkg          DY         Other
process      0             1              2           3
rate         XXSIGRATEXX    {fs}           {DY}       {other}
----------------------------------------------------------------------------------------------------------------------------------
fs_stat_{label} gmN  {fs_int}  -            {tf}             -         - 
fs_unc lnN             -            {tf_e}             -          -
jec   lnN           XXjecXX            -             -          -
El    lnN           XXElXX             -             -          -
Mu    lnN           XXMuXX             -             -          -
FastSimEl   lnN     XXFastSimElXX      -             -          -
FastSimMu   lnN     XXFastSimMuXX      -             -          -
SigTrig  lnN          1.05             -             -          -
bHe   lnN           XXbHeXX            -             -          -
bLi   lnN           XXbLiXX            -             -          -
PU    lnN           {pu_Up}/{pu_Dn}            -             -          - 
genMet lnU              XXgenMetXX     -             -          -
ISR   lnN              XXISRXX     -             -          - 
signalMCstats_{label} lnN    XXmcStatXX  -           -          -
dySys  lnN           -                  -            {dy_e}        - 
otherSys lnN         -                  -            -          {other_e}
lumi   lnN             1.026              -            -          - 
'''.format(label = label, 
           obs = da_SFBelow.GetBinContent(bin) if SR < 7 else da_SFAbove.GetBinContent(bin),
           fs  = da_OFBelow.GetBinContent(bin)*tf_CR_SRBelow.GetBinContent(bin) if SR < 7 else da_OFAbove.GetBinContent(bin)*tf_CR_SRAbove.GetBinContent(bin),
           DY = dy_shapeBelow.GetBinContent(bin) if SR < 7 else dy_shapeAbove.GetBinContent(bin),
           other = mc_fullBelow.GetBinContent(bin) if SR < 7 else mc_fullAbove.GetBinContent(bin),
           fs_int = int(da_OFBelow.GetBinContent(bin)) if SR < 7 else int(da_OFAbove.GetBinContent(bin)),
           tf = tf_CR_SRBelow.GetBinContent(bin) if SR < 7 else tf_CR_SRAbove.GetBinContent(bin),
           tf_e = (1 + tf_CR_SRBelow.GetBinError(bin)/tf_CR_SRBelow.GetBinContent(bin)) if SR < 7  else (1 + tf_CR_SRAbove.GetBinError(bin)/tf_CR_SRAbove.GetBinContent(bin)),
           dy_e = ( 1.5 if dy_shapeBelow.GetBinContent(bin) == 0 else 1 + dy_shapeBelow.GetBinError(bin)/dy_shapeBelow.GetBinContent(bin) ) if SR < 7 else ( 1.5 if dy_shapeAbove.GetBinContent(bin) == 0 else 1 + dy_shapeAbove.GetBinError(bin)/dy_shapeAbove.GetBinContent(bin) ),
           pu_Up = 1 + puSystUp.GetBinContent(bin),
           pu_Dn = 1 + puSystDn.GetBinContent(bin),
           other_e = (1 + mc_fullBelow.GetBinError(bin) / mc_fullBelow.GetBinContent(bin)) if SR < 7 else (1 + mc_fullAbove.GetBinError(bin) / mc_fullBelow.GetBinContent(bin)))

        outputFile = open('datacards/datacards_{scan}/{scan}/datacard_{sr}.txt'.format(scan=scan.name,sr=label),'w')
        outputFile.write(datacard)
        outputFile.close()


# def runCmd(cmd):
#     print cmd[2]
#     os.chdir (cmd[2])
#     print cmd[0]
#     os.system(cmd[0])
#     print cmd[1]
#     os.system(cmd[1])
#     return True
def runCmd(cmd):
    if os.path.exists('/pool/'):
        command = '''#!/bin/sh                                                                             
export SCRAM_ARCH=slc6_amd64_gcc530                                                                                 
export VO_CMS_SW_DIR=/cvmfs/cms.cern.ch                                                                             
source ${VO_CMS_SW_DIR}/cmsset_default.sh                                                                           
export CMS_PATH=${VO_CMS_SW_DIR}                                                                                    
cmsenv                                                                                                              
cd /afs/cern.ch/user/m/mvesterb/edgeSW/combine/CMSSW_7_4_7/src/
cmsenv\n'''                                                                                            
    else:
        command = '''#!/bin/sh                                                                             
export SCRAM_ARCH=slc6_amd64_gcc530                                                                                 
export VO_CMS_SW_DIR=/cvmfs/cms.cern.ch                                                                             
source ${VO_CMS_SW_DIR}/cmsset_default.sh                                                                           
export CMS_PATH=${VO_CMS_SW_DIR}                                                                                    
cmsenv                                                                                                              
cd /afs/cern.ch/work/m/mvesterb/private/Edge/produceLimits/CMSSW_7_4_7/src/                                         
cmsenv\n'''                                                                                            
    command = command+'cd ' + cmd[2] + '\n'
    command = command+cmd[0] + '\n'
    command = command+cmd[1] + '\n'
    randomHash = ''.join(random.choice(string.ascii_uppercase+string.digits) for _ in range(15))
    fil = open('command'+randomHash + '.sh','w')
    fil.write(command)
    fil.close()
    # running command
    cmd2 = '%s/%s'%(os.getcwd(),'command'+randomHash + '.sh')
    os.chmod(cmd2, 0o755)
    subprocess.call(cmd2,shell=True)
#    os.system('rm command'+randomHash + '.sh')


    return True


def produceLimits( njobs ):
    basedir = os.path.abspath('.')+'/datacards/datacards_{name}/{name}/'.format(name=scan.name)
    subdirs = [i for i in os.listdir(basedir) if os.path.isdir(basedir+i)]
    pool = Pool(njobs)
    tasks = []
    for ind,d in enumerate(subdirs):
        mass = ''.join(i for i in d.split('_') if i.isdigit() )
        fd = basedir+d
        # srCards = ''
        # for SR, label in scan.shortLabels.items():
        #     srCards += ' {fd}/datacard_{mass}_{label}.txt'.format(fd=fd,mass=d,label=label)

        # some SRs are empty so, better this way
        srCards = ' {fd}/datacard_{mass}_*.txt'.format(fd=fd,mass=d)
        dc_name    = '{fd}/datacard_{suffix}.txt'.format(fd=fd,suffix=d)
        #runcmd = 'combineCards.py -S 0'
        runcmd = 'combineCards.py -S {bstr} > {final_dc}'.format(bstr=srCards,final_dc=dc_name)
        combinecmd = 'combine -m {mass} -M Asymptotic {dc_name}'.format(mass=mass,dc_name=dc_name)
        combinecmd = combinecmd + '\n' + 'combine -m {mass} -M  ProfileLikelihood  --uncapped 1 --significance --rMin -5 {dc_name}'.format(mass=mass,dc_name=dc_name) # for significance maps
#        print runcmd
        tasks.append([runcmd,combinecmd, fd])
    pool.map(runCmd, tasks)
    haddcmd = 'hadd -f {bd}/{name}_allLimits.root {bd}/*/higgs*Asymptotic*.root'.format(name=scan.name,bd=basedir)
    os.system(haddcmd)
    haddcmd = 'hadd -f {bd}/{name}_allSigs.root {bd}/*/higgs*ProfileLikelihood*.root'.format(name=scan.name,bd=basedir)
    os.system(haddcmd)


def adaptBinning(target, current):
    final = copy.deepcopy(target)
    final.Reset()
    final = final.Project3D('yx')
    current_yx = current.Project3D('yx')
    ret_histo = copy.deepcopy(target)
    ret_histo.Reset()
    for i in range(1,final.GetXaxis().GetNbins()+1):
        for j in range(1, final.GetYaxis().GetNbins()+1):
            xval = final.GetXaxis().GetBinCenter(i)
            yval = final.GetYaxis().GetBinCenter(j)
            cont = current_yx.GetBinContent(current_yx.FindBin(xval, yval))
            err  = current_yx.GetBinError  (current_yx.FindBin(xval, yval))
            final.SetBinContent(final.FindBin(xval, yval), cont)
            final.SetBinError  (final.FindBin(xval, yval), err )
            for k in range(1,ret_histo.GetZaxis().GetNbins()+1):
                ret_histo.SetBinContent(ret_histo.FindBin(xval, yval, ret_histo.GetZaxis().GetBinCenter(k)), cont)
                ret_histo.SetBinError  (ret_histo.FindBin(xval, yval, ret_histo.GetZaxis().GetBinCenter(k)), err )
    return final, ret_histo


def getEffMapsSys(sys):
    print ('getting scan for sys', sys) if len(sys) > 0 else 'getting nominal scan'
    theCuts = scan.cuts_norm
    srId = scan.srID
    for rpl in replaceCutsForSys[sys]:
        print 'replacing',rpl[0],rpl[1]
        theCuts = theCuts.replace(rpl[0],rpl[1])
        srId    = srId.replace(rpl[0],rpl[1])
    print '##############'
    print "sys ", sys
    print "the Cuts ", theCuts
    print "srId ", srId
    print '##############'
    effMap = scan.tree.getTH3F(1., 'nPass_norm'+sys, srId+':'+scan.yvar+':'+scan.xvar,
                               scan.xbins.n+1, scan.xbins._min-scan.xbins.w/2.,
                               scan.xbins._max+scan.xbins.w/2., scan.ybins.n+1,
                               scan.ybins._min-scan.ybins.w/2., 
                               scan.ybins._max+scan.ybins.w/2., 
                               scan.srIDMax+1, -0.5, scan.srIDMax+0.5, theCuts, '',
                               scan.xtitle, scan.ytitle, scan.ztitle, 
                               extraWeightsForSys[sys])
#    if sys == '':
#        kk = TFile('map.root','recreate')
#        effMap.Write()
#        scan.ngen_3d.Write()
#        kk.Close()
        
    #     print mierda
        
    effMap.GetXaxis().GetBinCenter(1)
    print effMap.Integral()
    effMap.Divide(scan.ngen_3d)
    print 'kk is ', effMap.Integral()
    return effMap

def getSREffMaps():
    if hasattr(scan,'maps'): 
        if '' in scan.maps: effmap =  scan.maps[''] 
    else: effmap = getEffMapsSys('')

    tmp = effmap.Project3D('yx'); tmp.Reset()

    scan.effMapSR = {}
    effmaps = r.TFile('effmaps.root','recreate')
    for SR, label in scan.SRLabels.items():
        scan.effMapSR[SR] = tmp.Clone()
        for i in range(1, scan.effMapSR[SR].GetXaxis().GetNbins()+1):
            for j in range(1, scan.effMapSR[SR].GetYaxis().GetNbins()+1):
                bin = effmap.FindBin(effmap.GetXaxis().GetBinCenter(i),
                                     effmap.GetYaxis().GetBinCenter(j),
                                     SR)
                scan.effMapSR[SR].SetBinContent(i, j, effmap.GetBinContent(bin))
        scan.effMapSR[SR].SetName(label)
        scan.effMapSR[SR].SetTitle(label)
        scan.effMapSR[SR].Write()
    effmaps.Close()

def PutHistosIntoRootFiles():
#    print 'everything into datacards'
#    kk = TFile.Open('map2.root','recreate')
#    scan.maps[''].Write()
#    kk.Close()
    for i in range(1, scan.norm.GetXaxis().GetNbins()+1):
        for j in range(1, scan.norm.GetYaxis().GetNbins()+1):
            xval = scan.norm.GetXaxis().GetBinCenter(i)
            yval = scan.norm.GetYaxis().GetBinCenter(j)
            if (yval > xval): continue
            massString = 'mSlepton_%.0f_mchi1_%.0f'%(xval, yval)
            #massString = 'mSbottom_%.0f_mchi2_%.0f'%(xval, yval)
            print xval, yval
           
            sysHistos = {}
            for sys, histo in scan.maps.items():
                out = histo.ProjectionZ('%4.0f_%4.0f'%(xval,yval), i, i, j, j)
                out.Scale(scan.br*lumi*scan.xsecs[xval][0])
                out.SetName(massString + sys)
                sysHistos[sys] = out
#                if sys == '':
#                    tfile = TFile.Open(massString+'.root','recreate')
#                    out.Write()
#                    tfile.Close()
            print "sysHistos[''].Integral() ", sysHistos[''].Integral()
            if sysHistos[''].Integral() == 0: 
                print 'no sensitivty for', massString
                continue # if no sensitivity
            
            helper.ensureDirectory('datacards/datacards_{scan}/{scan}/{mass}/'.format(scan=scan.name,
                                                                                      mass=massString))
            for SR, label in scan.shortLabels.items():
                itsOk = True
                if scan.makeMCDatacards:
                    template = open('datacards/datacards_{scan}/{scan}/template_{sr}.txt'.format(scan=scan.name,sr=label),'r').read()
                else:
                    template = open('datacards/datacards_{scan}/{scan}/datacard_{sr}.txt'.format(scan=scan.name,sr=label),'r').read()   
                template = template.replace('XXSIGNALXX',massString)
                # central value is the average of nominal and genMET -.-
                nom = sysHistos['']        .GetBinContent(sysHistos['']        .FindBin(SR))
                if nom == 0: 
                    print 'no events in', label, 'for', i, j
                    continue
                var = sysHistos['genMet']  .GetBinContent(sysHistos['genMet']  .FindBin(SR))
                nom = (nom+var) / 2 
                template = template.replace('XXSIGRATEXX', '%4.4f'%nom)
                

                for sys in scan.SysStringUpDown.split():
                    up = sysHistos[sys+'Up']  .GetBinContent(sysHistos[sys+'Up']  .FindBin(SR))
                    dn = sysHistos[sys+'Down'].GetBinContent(sysHistos[sys+'Down'].FindBin(SR))
                    nom= sysHistos['']        .GetBinContent(sysHistos['']        .FindBin(SR))
                    scan.SysForTableMax[sys] = max( abs(up/nom-1), abs(dn/nom-1), scan.SysForTableMax[sys])
                    scan.SysForTableMin[sys] = min( abs(up/nom-1), abs(dn/nom-1), scan.SysForTableMin[sys])

                    template = template.replace('XX'+sys+'XX', '%4.4f/%4.4f'%(dn/nom,up/nom))
                    if dn/nom < 1e-1000 or up/nom < 1e-1000: 
                    #if dn/nom < 1e-4 or up/nom < 1e-4: 
                        print 'theres an issue in', xval, yval, label, sys, ' probably related to not having enough mc in that region. Skiping that region'
                        itsOk = False

                if not itsOk: continue
                for sys in scan.SysString.split():
                    var = sysHistos[sys]  .GetBinContent(sysHistos[sys]  .FindBin(SR))
                    nom= sysHistos['']    .GetBinContent(sysHistos['']   .FindBin(SR))
                    scan.SysForTableMax[sys] = max( abs(var/nom-1), scan.SysForTableMax[sys])
                    scan.SysForTableMin[sys] = min( abs(var/nom-1), scan.SysForTableMin[sys])

                    # gen met variation is different...
                    if sys == 'genMet':
                        var = (nom+var) / 2 
                    template = template.replace('XX'+sys+'XX', '%4.4f'%(var/nom))
                mcStat = sysHistos[''].GetBinError(sysHistos[''].FindBin(SR)) / sysHistos[''].GetBinContent(sysHistos[''].FindBin(SR))
                print sysHistos[''].GetBinError(sysHistos[''].FindBin(SR)), sysHistos[''].GetBinContent(sysHistos[''].FindBin(SR))
                template = template.replace('XXmcStatXX', '%f'%(1 + mcStat))
                card = open('datacards/datacards_{scan}/{scan}/{mass}/datacard_{mass}_{label}.txt'.format(scan=scan.name,
                                                                                                          mass=massString,
                                                                                                          label=label),'w')
                card.write(template)
                card.close()

def makeTable():
    print 'Source of uncertainty                 &    Uncertainty (\%)'
    print 'Luminosity                            &    6.2 \\'
    print 'Pileup                                &    %4.0f-%4.0f\\'%(scan.SysForTableMin['PU'],scan.SysForTableMax['PU'])
    print 'b tag modelling                       &    %4.0f-%4.0f\\'%(math.sqrt( scan.SysForTableMin['bHe']**2 + scan.SysForTableMin['bLi']**2), math.sqrt( scan.SysForTableMax['bHe']**2 + scan.SysForTableMax['bLi']**2))
    print 'Lepton reconstruction and isolation   &    %4.0f-%4.0f\\'%(math.sqrt( scan.SysForTableMin['Mu']**2 + scan.SysForTableMin['El']**2), math.sqrt( scan.SysForTableMax['Mu']**2 + scan.SysForTableMax['El']**2))
    print 'Fast simulation scale factors         &    %4.0f-%4.0f\\'%(math.sqrt( scan.SysForTableMin['FastSimMu']**2 + scan.SysForTableMin['FastSimEl']**2), math.sqrt( scan.SysForTableMax['FastSimMu']**2 + scan.SysForTableMax['FastSimEl']**2))
    print 'Trigger modelling                     &    5 \\'
    print 'Jet energy scale                      &    %4.0f-%4.0f\\'%(scan.SysForTableMin['jec'],scan.SysForTableMax['jec'])
    print 'ISR modelling TODO \\'
    print 'Statistical uncertainty               &    %4.0f-%4.0f\\'%()
    
     
if __name__ == "__main__":

    r.gStyle.SetPaintTextFormat(".2f")
    r.gStyle.SetOptStat(0)

    parser = optparse.OptionParser(usage="usage: %prog [opts] FilenameWithSamples", version="%prog 1.0")
    parser.add_option('-s', '--samples', action='store', type=str, dest='sampleFile', default='samples.dat', help='the samples file. default \'samples.dat\'')
    parser.add_option('-r', '--reloadScan' , action='store_true', dest='reloadScan', help='reload scan. default %default')
    parser.add_option('-l', '--reloadLimits' , action='store_true', dest='reloadLimits', help='reload limits. default %default')
    parser.add_option('-n', '--scanName'     , action='store', type=str, dest='scanName', default='T6bbslepton', help='scan name. e.g. T6bbslepton (this is the default)')
    (opts, args) = parser.parse_args()

    print 'running with these options'
    for key, value in opts.__dict__.items():
        print '%-20s : %-20s' %(key, value)
    print ' \n\n'
    print 'Going to load the tree(s)...'

    global lumi, scan

    cuts = CutManager.CutManager()
    lumi = 35.9
    global replaceCutsForSys, extraWeightsForSys
    replaceCutsForSys = {'':      [],
                         'jecUp'       : [['nBJetMedium35_Edge','nBJetMedium35_jecUp_Edge'],
                                          ['nBJetMedium25_Edge','nBJetMedium25_jecUp_Edge'],
                                          ['mt2_Edge','mt2_jecUp_Edge'],
                                          ['mbb_Edge','mbb_jecUp_Edge'],
                                          ['nll_Edge', 'nll_jecUp_Edge'],
                                          ['met_Edge', 'met_jecUp_Edge'],
                                          ['dphiMjj_Edge', 'dphiMjj_jecUp_Edge'],
                                          ['mt2bb_Edge', 'mt2bb_jecUp_Edge'],
                                          ['mbb_Edge', 'mbb_jecUp_Edge'],
                                          ['htJet35j_Edge', 'htJet35j_jecUp_Edge'],
                                          ['nJetSel_Edge','nJetSel_jecUp_Edge'],
                                          ['FS_central_jets_Edge','FS_central_jets_jecUp_Edge']],
                         'jecDown'     : [['nBJetMedium35_Edge','nBJetMedium35_jecDn_Edge'],
                                          ['nBJetMedium25_Edge','nBJetMedium25_jecDn_Edge'],
                                          ['mt2_Edge','mt2_jecDn_Edge'],
                                          ['mbb_Edge','mbb_jecDn_Edge'],
                                          ['nll_Edge', 'nll_jecDn_Edge'],
                                          ['met_Edge', 'met_jecDn_Edge'],
                                          ['dphiMjj_Edge', 'dphiMjj_jecDn_Edge'],
                                          ['mt2bb_Edge', 'mt2bb_jecDn_Edge'],
                                          ['mbb_Edge', 'mbb_jecDn_Edge'],
                                          ['htJet35j_Edge', 'htJet35j_jecDn_Edge'],
                                          ['nJetSel_Edge','nJetSel_jecDn_Edge'],
                                          ['FS_central_jets_Edge','FS_central_jets_jecDn_Edge']],
                         'bHeUp'        : [],
                         'bHeDown'      : [],
                         'bLiUp'        : [],
                         'bLiDown'      : [],
                         'ElUp'         : [],
                         'ElDown'       : [],
                         'FastSimElUp'  : [],
                         'FastSimElDown': [],
                         'FastSimMuUp'  : [],
                         'FastSimMuDown': [],
                         'MuUp'   : [],
                         'MuDown' : [],
                         'PUUp'   : [],
                         'PUDown' : [],
                         'ISRUp'  : [],
                         'ISRDown': [],
                         'genMet' : [['met_Edge','genMet_Edge'],
                                     ['nll_Edge','nll_genMet_Edge']]}
    
    extraWeightsForSys = {''          : '1',
                          'jecUp'     : '1',
                          'jecDown'   : '1',
                          'bHeUp'     : 'weight_btagsf_heavy_UP_Edge / weight_btagsf_Edge',
                          'bHeDown'   : 'weight_btagsf_heavy_DN_Edge / weight_btagsf_Edge',
                          'bLiUp'     : 'weight_btagsf_light_UP_Edge / weight_btagsf_Edge',
                          'bLiDown'   : 'weight_btagsf_light_DN_Edge / weight_btagsf_Edge',
                          'ElUp'      : 'LepSFElUp(Lep1_pt_Edge,Lep1_eta_Edge,Lep1_pdgId_Edge)*LepSFElUp(Lep2_pt_Edge,Lep2_eta_Edge,Lep2_pdgId_Edge) / ( LepSF(Lep1_pt_Edge,Lep1_eta_Edge,Lep1_pdgId_Edge)*LepSF(Lep2_pt_Edge,Lep2_eta_Edge,Lep2_pdgId_Edge) )',
                          'ElDown'      : 'LepSFElDn(Lep1_pt_Edge,Lep1_eta_Edge,Lep1_pdgId_Edge)*LepSFElDn(Lep2_pt_Edge,Lep2_eta_Edge,Lep2_pdgId_Edge) / ( LepSF(Lep1_pt_Edge,Lep1_eta_Edge,Lep1_pdgId_Edge)*LepSF(Lep2_pt_Edge,Lep2_eta_Edge,Lep2_pdgId_Edge) )',
                          'FastSimElUp'  : 'LepSFFastSimElUp(Lep1_pt_Edge,Lep1_eta_Edge,Lep1_pdgId_Edge)*LepSFFastSimElUp(Lep2_pt_Edge,Lep2_eta_Edge,Lep2_pdgId_Edge) / ( LepSFFastSim(Lep1_pt_Edge,Lep1_eta_Edge,Lep1_pdgId_Edge)*LepSFFastSim(Lep2_pt_Edge,Lep2_eta_Edge,Lep2_pdgId_Edge) )',
                          'FastSimElDown': 'LepSFFastSimElDn(Lep1_pt_Edge,Lep1_eta_Edge,Lep1_pdgId_Edge)*LepSFFastSimElDn(Lep2_pt_Edge,Lep2_eta_Edge,Lep2_pdgId_Edge) / ( LepSFFastSim(Lep1_pt_Edge,Lep1_eta_Edge,Lep1_pdgId_Edge)*LepSFFastSim(Lep2_pt_Edge,Lep2_eta_Edge,Lep2_pdgId_Edge) )',
                          'FastSimMuUp'  : 'LepSFFastSimMuUp(Lep1_pt_Edge,Lep1_eta_Edge,Lep1_pdgId_Edge)*LepSFFastSimMuUp(Lep2_pt_Edge,Lep2_eta_Edge,Lep2_pdgId_Edge) / ( LepSFFastSim(Lep1_pt_Edge,Lep1_eta_Edge,Lep1_pdgId_Edge)*LepSFFastSim(Lep2_pt_Edge,Lep2_eta_Edge,Lep2_pdgId_Edge) )',
                          'FastSimMuDown': 'LepSFFastSimMuDn(Lep1_pt_Edge,Lep1_eta_Edge,Lep1_pdgId_Edge)*LepSFFastSimMuDn(Lep2_pt_Edge,Lep2_eta_Edge,Lep2_pdgId_Edge) / ( LepSFFastSim(Lep1_pt_Edge,Lep1_eta_Edge,Lep1_pdgId_Edge)*LepSFFastSim(Lep2_pt_Edge,Lep2_eta_Edge,Lep2_pdgId_Edge) )',
                          'MuUp'      : 'LepSFMuUp(Lep1_pt_Edge,Lep1_eta_Edge,Lep1_pdgId_Edge)*LepSFMuUp(Lep2_pt_Edge,Lep2_eta_Edge,Lep2_pdgId_Edge) / ( LepSF(Lep1_pt_Edge,Lep1_eta_Edge,Lep1_pdgId_Edge)*LepSF(Lep2_pt_Edge,Lep2_eta_Edge,Lep2_pdgId_Edge) )',
                          'MuDown'      : 'LepSFMuDn(Lep1_pt_Edge,Lep1_eta_Edge,Lep1_pdgId_Edge)*LepSFMuDn(Lep2_pt_Edge,Lep2_eta_Edge,Lep2_pdgId_Edge) / ( LepSF(Lep1_pt_Edge,Lep1_eta_Edge,Lep1_pdgId_Edge)*LepSF(Lep2_pt_Edge,Lep2_eta_Edge,Lep2_pdgId_Edge) )',
                          'PUUp'      : 'PileupW_Up_Edge / PileupW_Edge',
                          'PUDown'    : 'PileupW_Dn_Edge / PileupW_Edge',
                          'ISRUp'     : 'ISRweight_Up_Edge / ISRweight_Edge',
                          'ISRDown'   : 'ISRweight_Dn_Edge / ISRweight_Edge',
                          'genMet'    : '1'}



    ## ==============================================
    ## == try loading stuff from the pickled file ===
    ## ==============================================

    
    pickfile = 'datacards/datacards_{name}/{name}/{name}.pkl'.format(name=opts.scanName)
    if os.path.isfile(pickfile) and not opts.reloadScan:
        print 'getting the scan object from the pickled file'
        scan = pickle.load(open(pickfile,'r'))
    else:
        ## everything that takes long should be done here!
        scan = Scans.Scan(opts.scanName)
        if scan.makeMCDatacards:
            print "doing makeDataCards"
            makeMCDatacards()
        if scan.makeMCDatacardsSlepton:
            makeDataCardsFromRootFileSlepton()
        else:
            if not 'Edge' in scan.name:
                makeDataCardsFromRootFile()
            else:
                makeDataCardsFromRootFileForEdge()

        scan.tree = Sample.Tree(helper.selectSamples("samplesSlepton.dat", scan.datasets, 'SIG'), 'SIG'  , 0, isScan = True)
        #scan.tree = Sample.Tree(helper.selectSamples(opts.sampleFile, scan.datasets, 'SIG'), 'SIG'  , 0, isScan = True)
        ## Load the number of generated events to produce efficiency maps per systematic
        print "getting treeeees     ____________________________________________________________________________________________________"
        scan.dummy = scan.tree.getTH3F(1., 'dummy', '1:1:1',
                                       scan.xbins.n+1, scan.xbins._min-scan.xbins.w/2.,
                                       scan.xbins._max+scan.xbins.w/2., scan.ybins.n+1,
                                       scan.ybins._min-scan.ybins.w/2., 
                                       scan.ybins._max+scan.ybins.w/2., 
                                       scan.srIDMax+1, -0.5, scan.srIDMax+0.5, '1', '',
                                       scan.xtitle, scan.ytitle, scan.ztitle, 
                                       '1')
        if True:
            scan.ngen = scan.tree.blocks[0].samples[0].smsCount ## take the first slice's ngen histo
            print "scan.ngen ",  scan.ngen
            for ind,i in enumerate(scan.tree.blocks[0].samples):
                print "ind ", ind
                print "i ", i
                if ind: ## do not add the first one twice
                    scan.ngen.Add(i.smsCount, 1.) ## add all others
        else:
            scan.make3DGen()

        for i in range(1, scan.ngen.GetXaxis().GetNbins()+1):
            for j in range(1, scan.ngen.GetYaxis().GetNbins()+1):
                for k in range(1, scan.ngen.GetZaxis().GetNbins()+1):
                    scan.ngen.SetBinError(scan.ngen.GetBin(i,j,k),0.)

        print 'this is the type of scan.ngen before', type(scan.ngen)
        newbinning   = adaptBinning(scan.dummy, scan.ngen)
        scan.ngen_2d = newbinning[0]
        scan.ngen_3d = newbinning[1] ## this one has ngen in every single bin. for every SR. and it's 3D, so that's cool
        print 'this is the type of scan.ngen after', type(scan.ngen)
        getSREffMaps()

        scan.SysStringUpDown = 'El Mu jec bHe bLi ISR'
        #scan.SysStringUpDown = 'El Mu jec bHe bLi FastSimEl FastSimMu ISR'
        scan.SysString = 'genMet'
#        scan.SysStringUpDown = ''
        #scan.SysString = ''
        scan.SysForTableMax = {}
        scan.SysForTableMin = {}
        for sys in scan.SysStringUpDown.split()+scan.SysString.split():
            scan.SysForTableMax[sys] = 0.
            scan.SysForTableMin[sys] = -1.

        scan.maps = {}
        scan.maps['']   = getEffMapsSys('')
        print 'nominal is ', scan.maps[''].Integral()
        print 'normali systematics'
        for sys in scan.SysString.split():
            scan.maps[sys] = getEffMapsSys(sys)
            print 'integral is', scan.maps[sys].Integral()
        print 'up/down systematics'
        for sys in scan.SysStringUpDown.split():
            scan.maps[sys+'Up']   = getEffMapsSys(sys+'Up')
            scan.maps[sys+'Down'] = getEffMapsSys(sys+'Down')

        scan.norm = scan.maps['']
        print 'getting serveral maps'
        scan.xy = scan.norm.Project3D('xy') # this means x versus y. so x is on the y-axis
        scan.xz = scan.norm.Project3D('xz')
        scan.yx = scan.norm.Project3D('yx') # that's the inclusive efficiency map
        scan.yz = scan.norm.Project3D('yz')
        scan.zx = scan.norm.Project3D('zx')
        scan.zy = scan.norm.Project3D('zy')
 
        PutHistosIntoRootFiles()
        print 'done. now calculating limits'
        

    if opts.reloadLimits:
        print 'reloading limits and datacards'
 #        fillAndSaveDatacards(dobs)
        produceLimits(40 if os.path.exists('/pool/') else 8)

    os.system('mkdir -p mkdir -p makeExclusionPlot/config/%s/'%scan.paper)
    scan.makeSignificanceMap()
    scan.makeExclusion()
    scan.makePrettyPlots()

    ## save the scan object in a pickle file to save time the second time around.
    #pickle.dump(scan, open('datacards/datacards_{name}/{name}/{name}.pkl'.format(name=scan.name),'w') )
    print 'marc is stupid'
