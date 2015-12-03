import ROOT as r
from   ROOT import gROOT, TCanvas, TFile, TF1, TPaveStats
import math, sys, optparse, copy, re, array


import include.helper     as helper
import include.Region     as Region
import include.Canvas     as Canvas
import include.CutManager as CutManager
import include.Sample     as Sample
import include.Tables     as Tables

##def runLimits():

def fillAndSaveDatacards(nb):
    for eta in ['central', 'forward']:
        for mass in ['lowMass', 'belowZ', 'onZ', 'aboveZ', 'highMass']:
            tmp_file = open('datacards/%s_%s_%s.txt'%(eta, mass, nb) ,'r')
            tmp_dc  = tmp_file.read()
            tmp_histo = globals()['eff_%s_%s_%s'%(eta, mass, nb)]
            for i in range(1, tmp_histo.GetXaxis().GetNbins()+1):
                for j in range(1, tmp_histo.GetYaxis().GetNbins()+1):
                    xmass = tmp_histo.GetXaxis().GetBinLowEdge(i)
                    ymass = tmp_histo.GetYaxis().GetBinLowEdge(j)
                    if tmp_histo.GetBinContent(tmp_histo.GetBin(i,j)) == 0. or ymass > xmass:
                        continue
                    mass_string = 'mSbottom_%.0f_mchi2_%.0f'%(xmass, ymass)
                    helper.ensureDirectory('datacards/T6bbslepton/%s'%(mass_string))
                    tmp_rate = lumi*xsecs[xmass][0]*tmp_histo.GetBinContent(tmp_histo.GetBin(i,j))
                    tmp_out = tmp_dc.replace('XXRATEXX', '%.3f'%(tmp_rate))
                    tmp_new = open('datacards/T6bbslepton/%s/%s'%(mass_string, tmp_file.name.split('/')[-1]), 'w')
                    tmp_new.write(tmp_out)
                    tmp_new.close()
            tmp_file.close()


def adaptBinning(target, current):
    final = copy.deepcopy(target)
    final.Reset()
    final = final.Project3D('yx')
    current_yx = current.Project3D('yx')
    ret_histo = copy.deepcopy(target)
    ret_histo.Reset()
    for i in range(1,final.GetXaxis().GetNbins()+1):
        for j in range(1, final.GetYaxis().GetNbins()+1):
            xval = final.GetXaxis().GetBinLowEdge(i)
            yval = final.GetYaxis().GetBinLowEdge(j)
            cont = current_yx.GetBinContent(current_yx.FindBin(xval, yval))
            err  = current_yx.GetBinError  (current_yx.FindBin(xval, yval))
            final.SetBinContent(final.FindBin(xval, yval), cont)
            final.SetBinError  (final.FindBin(xval, yval), err )
            for k in range(1,ret_histo.GetZaxis().GetNbins()+1):
                ret_histo.SetBinContent(ret_histo.FindBin(xval, yval, ret_histo.GetZaxis().GetBinLowEdge(k)), cont)
                ret_histo.SetBinError  (ret_histo.FindBin(xval, yval, ret_histo.GetZaxis().GetBinLowEdge(k)), err )
    return final, ret_histo

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
        mllid = [1]
    elif 'below' in mll:
        mllid = [2]
    elif 'on'    in mll:
        mllid = [3]
    elif 'above' in mll:
        mllid = [4]
    elif 'high'  in mll:
        mllid = [5]
    elif 'all'   in mll:
        mllid = [1,2,3,4,5]

    allSRs = []
    for i in nbid:
        for j in mllid:
            allSRs.append(100*etaid + 10*i + j)

    #print 'all srs for %s in %s and %s are %s'%(nb, eta, mll, str(allSRs))
    
    ret_histo = scan_eff_norm.Clone('eff_%s_%s_%s'%(eta, nb, mll))
    ret_histo = ret_histo.Project3D('yx')
    ret_histo.Reset()
    for i in range(1, scan_eff_norm.GetNbinsX()+1):
        for j in range(1, scan_eff_norm.GetNbinsY()+1):
            tmp_eff  = 0.
            tmp_err2 = 0.
            for k in allSRs:
                tmp_eff  +=  scan_eff_norm.GetBinContent(scan_eff_norm.GetBin(i, j, scan_eff_norm.GetZaxis().FindBin(k)))
                tmp_err2 += (scan_eff_norm.GetBinError  (scan_eff_norm.GetBin(i, j, scan_eff_norm.GetZaxis().FindBin(k))))**2
            ret_histo.SetBinContent(ret_histo.GetBin(i,j), tmp_eff)
            ret_histo.SetBinError  (ret_histo.GetBin(i,j), math.sqrt(tmp_err2))
    return ret_histo

 
if __name__ == "__main__":


    parser = optparse.OptionParser(usage="usage: %prog [opts] FilenameWithSamples", version="%prog 1.0")
    parser.add_option('-s', '--samples', action='store', type=str, dest='sampleFile', default='samples.dat', help='the samples file. default \'samples.dat\'')
    (opts, args) = parser.parse_args()

    print 'running with these options'
    for key, value in opts.__dict__.items():
        print '%-20s : %-20s' %(key, value)
    print ' \n\n'
    print 'Going to load the tree(s)...'

    sigDatasets = ['SMS_T6bbllslepton_mSbottom-600To900_mLSP-200To800']
    treeSIG = Sample.Tree(helper.selectSamples(opts.sampleFile, sigDatasets, 'SIG'), 'SIG'  , 0, isScan = True)
    cuts = CutManager.CutManager()

    signalRegion     = Region.region('signalRegion', 
                                     [cuts.METJetsSignalRegion],
                                     ['mll'],
                                     [ [20., 70., 81., 101., 120., 13000.] ],
                                     False)

    ## have to think a way of reweighting the events with trigger and lepton SFs.
    ## weighting now done with isScan=True flag in samples.py

    global lumi, scan_norm, scan_eff_norm, xsecs
    lumi = 2.1
    ## this should be then sr-ID:m_slepton:m_sbottom for the final scan
    xvar = 'GenSusyMScan1'
    yvar = 'GenSusyMScan2'
    zvar = 't.srID_Edge'

    xvar_title = 'm_{sbottom}'
    yvar_title = 'm_{neu2}'
    zvar_title = 'SR-ID'

    cuts_norm = cuts.AddList([cuts.METJetsSignalRegion, cuts.GoodLeptonSFNoTrigger()]) ## trigger not available in fastsim
    cuts_norm = cuts_norm.replace(cuts.twoLeptons, 't.nPairLep_Edge > 0') ## remove the filters, ugly, but it's a bit intricate in the samples

    scan_norm = treeSIG.getTH3F(lumi, 'nPass_norm', zvar+':'+yvar+':'+xvar,  32, 200, 1000, 32, 200, 1000, 200, 100, 300, cuts_norm, '', xvar_title, yvar_title, zvar_title)
    #do all the systematics

    ## we also need the number of generated events here!!
    scan_ngen = treeSIG.blocks[0].samples[0].smsCount ## careful here!!

    newbinning = adaptBinning(scan_norm, scan_ngen)
    scan_ngen_copy = newbinning[0]
    scan_ngen_3d   = newbinning[1] ## this one has ngen in every single bin. for every SR. andit's 3D, so that's cool

    scan_eff_norm = scan_norm.Clone('efficiency_norm')
    scan_eff_norm.Divide(scan_ngen_3d)

    xy = scan_eff_norm.Project3D('xy') # this means x versus y. so x is on the y-axis
    xz = scan_eff_norm.Project3D('xz')
    yx = scan_eff_norm.Project3D('yx') # that's the inclusive efficiency map
    yz = scan_eff_norm.Project3D('yz')
    zx = scan_eff_norm.Project3D('zx')
    zy = scan_eff_norm.Project3D('zy')


    ## here we get the 2D efficiencies for all signal regions

    for eta in ['central', 'forward']:
        for mass in ['allMass', 'lowMass', 'belowZ', 'onZ', 'aboveZ', 'highMass']:
            for nb in ['incb', '0b', '1b', '2b']:
                globals()['eff_%s_%s_%s'%(eta, mass, nb)] = getSRYield(eta, nb, mass)


    ## histogram and dictionary with the cross-sections and errors
    xsec_histo = r.TH1F('x-sections in fb for sbottom production','xsec_histo', 380, 100, 2000) 
    xsec_histo.Sumw2()
    xsecf = open('datacards/sbottomXsec.txt', 'r')
    xsecs = eval(xsecf.read())
    xsecf.close()
    for key, value in xsecs.items():
        xsecs[key][0] = xsecs[key][0]*1000.
        xsecs[key][1] = xsecs[key][0]*0.01*xsecs[key][1]

    for key,value in xsecs.items():
        xsec_histo.SetBinContent(xsec_histo.FindBin(key), value[0])
        xsec_histo.SetBinError  (xsec_histo.FindBin(key), value[1]) ## it's a percent value


    # now we take the default datacards and save them for each point into a subdirectory of the scan

    fillAndSaveDatacards('incb')

