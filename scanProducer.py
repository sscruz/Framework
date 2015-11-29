import ROOT as r
from   ROOT import gROOT, TCanvas, TFile, TF1, TPaveStats
import math, sys, optparse, copy, re, array


import include.helper     as helper
import include.Region     as Region
import include.Canvas     as Canvas
import include.CutManager as CutManager
import include.Sample     as Sample
import include.Tables     as Tables

def getSRYield(eta, nb, mll):
    if   eta == 'central':
        etaid = 1
    elif eta == 'forward':
        etaid = 2

    if 'inc' in nb:
        nbid = [0,1,2,3,4,5,6,7]
    elif '0' in nb:
        nbid = [0]
    elif '1' in nb:
        nbid = [1,2,3,4,5,6,7]
    elif '2' in nb:
        nbid = [2,3,4,5,6,7]

    if   'low'   in mll:
        mllid = 1
    elif 'below' in mll:
        mllid = 2
    elif 'on'    in mll:
        mllid = 3
    elif 'above' in mll:
        mllid = 4
    elif high    in mll:
        mllid = 5

    allSRs = []
    for i in nbid:
        allSRs.append(100*etaid + 10*i + mllid)
    
    return allSRs

 
if __name__ == "__main__":


    parser = optparse.OptionParser(usage="usage: %prog [opts] FilenameWithSamples", version="%prog 1.0")
    parser.add_option('-s', '--samples', action='store', type=str, dest='sampleFile', default='samples.dat', help='the samples file. default \'samples.dat\'')
    (opts, args) = parser.parse_args()

    print 'running with these options'
    for key, value in opts.__dict__.items():
        print '%-20s : %-20s' %(key, value)
    print ' \n\n'
    print 'Going to load DATA and MC trees...'

    sigDatasets = ['TTLep_pow']
    treeSIG = Sample.Tree(helper.selectSamples(opts.sampleFile, sigDatasets, 'SIG'), 'SIG'  , 0, isScan = True)
    cuts = CutManager.CutManager()

    signalRegion     = Region.region('signalRegion', 
                                     [cuts.METJetsSignalRegion],
                                     ['mll'],
                                     [ [20., 70., 81., 101., 120., 13000.] ],
                                     False)

    finalCuts = cuts.AddList([cuts.METJetsSignalRegion, cuts.GoodLeptonSF(), cuts.trigger]) ## the trigger will be a problem here
    ## have to think a way of reweighting the events with trigger and lepton SFs.

    lumi = 2.1
    ## this should be then sr-ID:m_slepton:m_sbottom for the final scan
    xvar = 't.srID_Edge'
    yvar = '(abs(t.Lep1_eta_Edge)<1.4 && abs(t.Lep2_eta_Edge)<1.4)'
    zvar = 't.nBJetMedium35_Edge'

    xvar_title = 'SR-ID'
    yvar_title = 'is central'
    zvar_title = 'n_{b-jets}'

    global scanHisto
    scanHisto = treeSIG.getTH3F(lumi, 'signalScan', zvar+':'+yvar+':'+xvar, 160, 100, 260, 2, -0.5, 1.5, 6, 0, 6, finalCuts, '', 'sr-ID', 'eta', 'nb')
    ## we also need the number of generated events here!!
    ##ngenHisto = treeSIG.getTH3F(lumi, 'signalScan', zvar+':'+yvar+':'+xvar, 160, 100, 260, 2, -0.5, 1.5, 6, 0, 6, finalCuts, '', 'sr-ID', 'eta', 'nb')

    effHisto = scanHisto.Clone('efficiencyMap')
    #effHisto.Divide(ngenHisto)

    ## so now: var1 is z-axis, var2 is y-axis, and var3 is x-axis

    xy = effHisto.Project3D('xy')
    xz = effHisto.Project3D('xz')
    yx = effHisto.Project3D('yx')
    yz = effHisto.Project3D('yz')
    zx = effHisto.Project3D('zx')
    zy = effHisto.Project3D('zy')

    ## so for scans, the lumiweight is set to unity. i.e. the efficiency histo is really that

    ## we need a translation between SR and central/fwd/etc.

    central_lowMass_incb_bins = getSRYield('central', 'inc', 'lowMass')

    ## now we need something that takes the default datacard and produces one for each point


