import ROOT as r
from   ROOT import gROOT, TCanvas, TFile, TF1, TPaveStats
import math, sys, optparse, copy, re, array, os, pickle


import include.helper     as helper
import include.Region     as Region
import include.Canvas     as Canvas
import include.CutManager as CutManager
import include.Sample     as Sample
import include.Tables     as Tables
import include.Scans      as Scans

from multiprocessing import Pool

def runCmd(cmd):
    os.system(cmd[0])
    os.chdir (cmd[1])
    os.system(cmd[2])
    return True


def produceLimits(whichb, njobs):
    basedir = os.path.abspath('.')+'/datacards/{name}/'.format(name=scan.name)
    subdirs = [i for i in os.listdir(basedir) if os.path.isdir(basedir+i)]
    pool = Pool(njobs)
    tasks = []
    for ind,d in enumerate(subdirs):
        mass = ''.join(i for i in d.split('_') if i.isdigit() )
        #if ind > 10: continue
        fd = basedir+d
        bstring = ''; b= ''
        for i in whichb:
            bstring += ' {fd}/{d}*_{b}.txt'.format(fd=fd,b=i,d=d)
            b+=i
        dc_name    = 'datacard_{massp}_full_{b}.txt'.format(massp=d,b=b)
        final_dc   = '{fd}/{dc_name}'.format(fd=fd,dc_name=dc_name)
        runcmd     = 'combineCards.py -S {bstr} > {final_dc}'.format(bstr=bstring,final_dc=final_dc)
        combinecmd = 'combine -m {mass} -M Asymptotic {dc_name}'.format(mass=mass,dc_name=dc_name)
        tasks.append([runcmd, fd, combinecmd])
    pool.map(runCmd, tasks)
    haddcmd = 'hadd -f {bd}/{name}_allLimits.root {bd}/*/higgs*.root'.format(name=scan.name,bd=basedir)
    #print haddcmd
    os.system(haddcmd)

def fillAndSaveDatacards(nbs):
    for region in scan.regions:
        eta = region[0]; mass = region[1]; nb = region[2]
        tmp_file = open('datacards/%s_%s_%s.txt'%(eta, mass, nb) ,'r')
        tmp_dc  = tmp_file.read()
        tmp_histo = getattr(scan, 'eff_%s_%s_%s'%(eta, mass, nb))
        for i in range(1, tmp_histo.GetXaxis().GetNbins()+1):
            for j in range(1, tmp_histo.GetYaxis().GetNbins()+1):
                xmass = tmp_histo.GetXaxis().GetBinCenter(i)
                ymass = tmp_histo.GetYaxis().GetBinCenter(j)
                if tmp_histo.GetBinContent(tmp_histo.GetBin(i,j)) == 0. or ymass > xmass:
                    continue
                mass_string = 'mSbottom_%.0f_mchi2_%.0f'%(xmass, ymass)
                helper.ensureDirectory('datacards/{name}'.format(name=scan.name))
                helper.ensureDirectory('datacards/{name}/{mass}'.format(name=scan.name,mass=mass_string))
                tmp_rate = lumi*scan.xsecs[xmass][0]*tmp_histo.GetBinContent(tmp_histo.GetBin(i,j)) ## LUMI WEIGHTING HAPPENS HERE!!!!
                tmp_out = tmp_dc.replace('XXRATEXX', '%.3f'%(tmp_rate))
                tmp_new = open('datacards/{name}/{mass}/{mass}_{fn}'.format(name=scan.name,mass=mass_string,fn=tmp_file.name.split('/')[-1]), 'w')
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

def getSRYield(eta, mll, nb):
    if   eta == 'central'  : etaids = [1]
    elif eta == 'forward'  : etaids = [2]
    elif eta == 'inclusive': etaids = [1,2]

    if 'inc' in nb: nbid = [0,1,2,3,4,5,6,7]
    elif '0' in nb: nbid = [0              ]
    elif '1' in nb: nbid = [  1,2,3,4,5,6,7]
    elif '2' in nb: nbid = [    2,3,4,5,6,7]

    if   'low'   in mll: mllids = [1]
    elif 'below' in mll: mllids = [2]
    elif 'on'    in mll: mllids = [3]
    elif 'above' in mll: mllids = [4]
    elif 'high'  in mll: mllids = [5]
    elif 'all'   in mll: mllids = [1,2,3,4,5]

    allSRs = []
    for etaid in etaids:
        for bid in nbid:
            for mllid in mllids:
                allSRs.append(100*etaid + 10*bid + mllid)

    ## print 'at signal region',eta, mll, nb
    ## print 'combining signal regions:',allSRs

    scan_sr_yield = scan.ngen_2d.Clone('yield_%s_%s_%s'%(eta, mll, nb))
    scan_sr_yield.SetTitle('yield_%s_%s_%s'%(eta, mll, nb))
    scan_sr_yield.SetName ('yield_%s_%s_%s'%(eta, mll, nb))
    scan_sr_yield.Reset()
    for x in range(1, scan.norm.GetNbinsX()+1):
        for y in range(1, scan.norm.GetNbinsY()+1):
            if y > x: continue
            tmp_yield  = 0.
            tmp_yield2 = 0.
            for z in allSRs:
                zbin = scan.norm.GetZaxis().FindBin(z)
                tmp_yield  +=  scan.norm.GetBinContent(scan.norm.GetBin(x, y, zbin))
                tmp_yield2 += (scan.norm.GetBinError  (scan.norm.GetBin(x, y, zbin)))**2
            scan_sr_yield.SetBinContent(scan_sr_yield.GetBin(x,y), tmp_yield)
            scan_sr_yield.SetBinError  (scan_sr_yield.GetBin(x,y), math.sqrt(tmp_yield2))
    scan_sr_eff = scan_sr_yield.Clone('eff_%s_%s_%s'%(eta, mll, nb))
    scan_sr_eff.SetTitle('eff_%s_%s_%s'%(eta, mll, nb))
    scan_sr_eff.GetZaxis().SetRangeUser(0., scan.zmaxEff)
    scan_sr_eff.Divide(scan.ngen_2d)
    setattr(scan, 'yield_%s_%s_%s'%(eta, mll, nb), scan_sr_yield)
    setattr(scan, 'eff_%s_%s_%s'  %(eta, mll, nb), scan_sr_eff  )

 
if __name__ == "__main__":

    r.gStyle.SetPaintTextFormat(".2f")

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

    ttDatasets = ['TTLep_pow']
    treeTT = Sample.Tree(helper.selectSamples(opts.sampleFile, ttDatasets, 'TT'), 'TT'  , 0, isScan = False)
    cuts = CutManager.CutManager()
    lumi = 2.3

    ## have to think a way of reweighting the events with trigger and lepton SFs.
    ## weighting now done with isScan=True flag in samples.py

    ## this should be then sr-ID:m_slepton:m_sbottom for the final scan
    xvar = 'GenSusyMScan1'
    yvar = 'GenSusyMScan2'
    zvar = 't.srID_Edge'

    #cuts_norm = cuts.AddList([cuts.METJetsSignalRegion, cuts.GoodLeptonSFNoTrigger()]) ## trigger not available in fastsim
    #cuts_norm = cuts_norm.replace(cuts.twoLeptons, 't.nPairLep_Edge > 0') ## remove the filters, ugly, but it's a bit intricate in the samples

    ## ==============================================
    ## == try loading stuff from the pickled file ===
    ## ==============================================
    
    pickfile = 'datacards/{name}/{name}.pkl'.format(name=opts.scanName)
    if os.path.isfile(pickfile) and not opts.reloadScan:
        print 'getting the scan object from the pickled file'
        scan = pickle.load(open(pickfile,'r'))
    else:
        ## everything that takes long should be done here!
        scan = Scans.Scan(opts.scanName)
        scan.tree = Sample.Tree(helper.selectSamples(opts.sampleFile, scan.datasets, 'SIG'), 'SIG'  , 0, isScan = True)
        scan.norm = scan.tree.getTH3F(1., 'nPass_norm', zvar+':'+yvar+':'+xvar,  scan.xbins.n+1, scan.xbins._min-scan.xbins.w/2., scan.xbins._max+scan.xbins.w/2.,  ## lumi set later for scans!!
                                                                                 scan.ybins.n+1, scan.ybins._min-scan.ybins.w/2., scan.ybins._max+scan.ybins.w/2., 
                                                                                 200, 100, 300, scan.cuts_norm, '', scan.xtitle, scan.ytitle, scan.ztitle)

        print '=================================================='
        print '=================================================='
        print '===== this is the cut ============================\n', scan.cuts_norm
        print '=================================================='
        print '=================================================='

        ## we also need the number of generated events here!!
        scan.ngen = scan.tree.blocks[0].samples[0].smsCount ## take the first slice's ngen histo
        for ind,i in enumerate(scan.tree.blocks[0].samples):
            if ind: ## do not add the first one twice
                scan.ngen.Add(i.smsCount, 1.) ## add all others

        ## =====================================================
        ## == if anything takes a long time, load it with pickle
        ## =====================================================


        ## ## look at the min_mlb distribution for a point or so
        ## ## eventually write a more flexible function to do some control plots
        ## sum_mlb_sms = scan.tree.getTH1F(lumi, 'min_mlb1_sms', 't.sum_mlb_Edge',  50, 0, 800, cuts.AddList([cuts_norm, 'GenSusyMScan1 == 750 && GenSusyMScan2 == 300']), '', 'sum_mlb')
        ## sum_mlb_tt  = treeTT   .getTH1F(lumi, 'min_mlb2_tt' , 't.sum_mlb_Edge',  50, 0, 800, cuts_norm, '', 'sum_mlb')


        print 'this is the type of scan.ngen before', type(scan.ngen)
        newbinning   = adaptBinning(scan.norm, scan.ngen)
        scan.ngen_2d = newbinning[0]
        scan.ngen_3d = newbinning[1] ## this one has ngen in every single bin. for every SR. and it's 3D, so that's cool
        print 'this is the type of scan.ngen after', type(scan.ngen)

        tmp_ngen = copy.deepcopy(scan.ngen)

        #print adfasdf

        scan.eff_norm = scan.norm.Clone('efficiency_norm')
        scan.eff_norm.Divide(scan.ngen_3d)

        scan.xy = scan.eff_norm.Project3D('xy') # this means x versus y. so x is on the y-axis
        scan.xz = scan.eff_norm.Project3D('xz')
        scan.yx = scan.eff_norm.Project3D('yx') # that's the inclusive efficiency map
        scan.yz = scan.eff_norm.Project3D('yz')
        scan.zx = scan.eff_norm.Project3D('zx')
        scan.zy = scan.eff_norm.Project3D('zy')


        ## here we get the 2D efficiencies for all signal regions
        for region in scan.regions:
            eta = region[0]; mass = region[1]; nb = region[2]
            getSRYield(eta, mass, nb)

        ## set the scan tree to 0 before pickling. it's huge (that's what she said)
        scan.tree = 0
        scan.ngen = tmp_ngen

    # now we take the default datacards and save them for each point into a subdirectory of the scan

    if opts.reloadLimits:
        print 'reloading limits and datacards'
        dobs = ['0b','1b']
        fillAndSaveDatacards(dobs)
        produceLimits(dobs,5)

    scan.makeExclusion()
    scan.makePrettyPlots()

    ## save the scan object in a pickle file to save time the second time around.
    pickle.dump(scan, open('datacards/{name}/{name}.pkl'.format(name=scan.name),'w') )
