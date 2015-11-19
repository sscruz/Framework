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
from   ROOT import gROOT, TCanvas, TFile, TF1, TPaveStats
import math, sys, optparse, copy

import include.helper     as helper
import include.Region     as Region
import include.Canvas     as Canvas
import include.CutManager as CutManager
import include.Sample     as Sample

def make_rsfof_table(reg):
    lines = []
    for eta in ['central', 'forward']:
        da_nof, da_nsf, da_rsfof, da_rsfof_err = getattr(reg, 'rsfof_yield_DATA_'+eta+'_of'), getattr(reg, 'rsfof_yield_DATA_'+eta+'_sf'), getattr(reg, 'rsfof_DATA_'+eta), getattr(reg, 'rsfof_DATA_'+eta+'_err') 
        mc_nof, mc_nsf, mc_rsfof, mc_rsfof_err = getattr(reg, 'rsfof_yield_MC_'+eta+'_of')  , getattr(reg, 'rsfof_yield_MC_'+eta+'_sf')  , getattr(reg, 'rsfof_MC_'+eta)  , getattr(reg, 'rsfof_MC_'+eta+'_err') 
        mc_nof_err, mc_nsf_err = getattr(reg, 'rsfof_yield_MC_'+eta+'_of_err')  , getattr(reg, 'rsfof_yield_MC_'+eta+'_sf_err')
        lines.append(' \\multicolumn{4}{c}{ \\textbf{%s}} \\\\' %(eta))
        lines.append(' \\multirow{2}{*}{Data}   & N$_{SF}$   & %5d   & \multirow{2}{*}{%.3f \\pm %.3f \\pm %.3f}  \\\\ ' %(da_nsf, da_rsfof, da_rsfof_err, 0.1) )
        lines.append('                          & N$_{OF}$   & %5d   &                                            \\\\ ' %(da_nof) )
        lines.append(' \\multirow{2}{*}{MC}     & N$_{SF}$   & %.0f \\pm %.0f & \multirow{2}{*}{%.3f \\pm %.3f \\pm %.3f}  \\\\ ' %(mc_nsf, mc_nsf_err, mc_rsfof, mc_rsfof_err, 0.1 ) )
        lines.append('                          & N$_{OF}$   & %.0f \\pm %.0f &                                            \\\\ ' %(mc_nof, mc_nof_err))
    for line in lines:
        print line
    return lines


def make_rsfof(histo_sf, histo_of, dataMC):

    ratio = histo_sf.Clone('rsfof_' + histo_sf.GetName())
    ratio.Divide(histo_of)
    ratio.GetYaxis().SetTitle('r_{SFOF}')

    doFit = 0
    if doFit:
        fit = TF1('myfit','pol0', ratio.GetXaxis().GetXmin(), ratio.GetXaxis().GetXmax())
        fit.SetLineColor(r.kBlack if dataMC  == 'DATA' else r.kRed+1)
        fit.SetLineStyle(2)
        
        ratio.Fit('myfit','0')

    ratio.GetYaxis().SetRangeUser(0.5,1.5)

    return ratio



if __name__ == '__main__':

    print 'Starting r_SFOF analysis...'
    parser = optparse.OptionParser(usage='usage: %prog [opts] FilenameWithSamples', version='%prog 1.0')
    parser.add_option('-s', '--samples', action='store', type=str, dest='sampleFile', default='samples.dat', help='the samples file. default \'samples.dat\'')
    parser.add_option('-m', '--mcstudies', action='store_true', dest='mcStudies', default=False, help='do MC studies. only loads one region and ttbar only')
    (opts, args) = parser.parse_args()

    print 'Going to load DATA and MC trees...'
    mcDatasets = ['TTLep_pow']# + ([] if opts.mcStudies else ['DYJetsToLL_M10to50', 'DYJetsToLL_M50'])
    daDatasets = ['DoubleMuon_Run2015C_25ns-05Oct_v1_runs_246908_260627' , 'DoubleEG_Run2015C_25ns-05Oct_v1_runs_246908_260627' , 'MuonEG_Run2015C_25ns-05Oct_v1_runs_246908_260627' ,
                  'DoubleMuon_Run2015D-05Oct_v1_runs_246908_260627'      , 'DoubleEG_Run2015D-05Oct_v1_runs_246908_260627'      , 'MuonEG_Run2015D-05Oct_v2_runs_246908_260627'      ,
                  'DoubleMuon_Run2015D_v4_runs_246908_260627'            , 'DoubleEG_Run2015D_v4_runs_246908_260627'            , 'MuonEG_Run2015D_v4_runs_246908_260627'            ]
    treeMC = Sample.Tree(helper.selectSamples(opts.sampleFile, mcDatasets, 'MC'), 'MC'  , 0)
    treeDA = Sample.Tree(helper.selectSamples(opts.sampleFile, daDatasets, 'DA'), 'DATA', 1)
    print 'Trees successfully loaded...'


    lumi = 2.1 ; maxrun = 999999

    lumi_str = 'lumi'+str(lumi).replace('.', 'p')+'_MCstudies'
    print 'Running with an integrated luminosity of', lumi,'fb-1'

    saveValues = True
   
    gROOT.ProcessLine('.L include/tdrstyle.C')
    gROOT.SetBatch(1)
    r.setTDRStyle() 
    cuts = CutManager.CutManager()

    regions = []
    ttjets_meas    = Region.region('ttjets_meas', 
                                   [cuts.METJetsControlRegion],
                                   ['mll'],
                                   [[20, 70, 81, 101, 120, 300]],
                                   True)
    regions.append(ttjets_meas)
    ttjets_meas_noM= Region.region('ttjets_meas_noM', 
                                   [cuts.METJetsControlRegion],
                                   ['mll'],
                                   [[20, 13000]],
                                   True)
    regions.append(ttjets_meas_noM)
    ttjets_meas_noZ= Region.region('ttjets_meas_noZ', 
                                   [cuts.METJetsControlRegion, cuts.DYmassVeto],
                                   ['mll'],
                                   [[20, 13000]],
                                   True)
    regions.append(ttjets_meas_noZ)
    ttjets_sig_lm  = Region.region('ttjets_sig_lm',
                                   [cuts.METJetsSignalRegion, cuts.lowmass],
                                   ['mll'],
                                   [[20, 70]],
                                   False)
    regions.append(ttjets_sig_lm)
    ttjets_sig_onZ = Region.region('ttjets_sig_onZ', 
                                   [cuts.METJetsSignalRegion, cuts.Zmass],
                                   ['mll'],
                                   [[81, 101]],
                                   False)
    regions.append(ttjets_sig_onZ)
    ttjets_sig_hm  = Region.region('ttjets_sig_hm', 
                                   [cuts.METJetsSignalRegion, cuts.highmass],
                                   ['mll'],
                                   [[120, 300]],
                                   False)
    regions.append(ttjets_sig_hm)
    ttjets_inc     = Region.region('ttjets_inc', 
                                   [cuts.nj2],
                                   ['mll', 'nj', 'met'],
                                   [[20, 70, 81, 101, 120, 300], range(2,9), range(0, 210, 10)],
                                   False)
    regions.append(ttjets_inc)

    if opts.mcStudies:
        regions = []
        regions.append(ttjets_inc)

    for reg in regions:
        print 'i am at region', reg.name
        for eta in ['central', 'forward']:
            print '... in %s' %(eta)

            for var in reg.rvars:

                if   var == 'mll':
                    treevar = 't.lepsMll_Edge'
                    varname = "m_{ll} (GeV)"
                elif var == 'nj':
                    treevar = 't.nJetSel_Edge'
                    varname = 'n_{jets}'
                elif var == 'met':
                    treevar = 'met_pt'
                    varname = 'ME_{T} (GeV)'


                #cuts_sf = cuts.AddList([cuts.MaxRun(maxrun), cuts.GoodLeptonSF(),cuts.trigger]+[cuts.Central() if eta == 'central' else cuts.Forward()]+reg.cuts)
                cuts_of = cuts.AddList([cuts.MaxRun(maxrun), cuts.GoodLeptonOF(),cuts.trigger]+[cuts.Central() if eta == 'central' else cuts.Forward()]+reg.cuts)
                cuts_ee = cuts.AddList([cuts.MaxRun(maxrun), cuts.GoodLeptonee(),cuts.trigger]+[cuts.Central() if eta == 'central' else cuts.Forward()]+reg.cuts)
                cuts_mm = cuts.AddList([cuts.MaxRun(maxrun), cuts.GoodLeptonmm(),cuts.trigger]+[cuts.Central() if eta == 'central' else cuts.Forward()]+reg.cuts)

                for tree in ([treeMC, treeDA] if reg.doData else [treeMC]):

                    dataMC = 'DATA' if tree == treeDA else 'MC'
                    setattr(reg, '%s_ee'%var, tree.getTH1F(lumi, '%s_ee_%s%s%s'%(var, eta, reg.name, dataMC), treevar, reg.bins[reg.rvars.index(var)], 1, 1, cuts_ee, '', varname))
                    setattr(reg, '%s_mm'%var, tree.getTH1F(lumi, '%s_mm_%s%s%s'%(var, eta, reg.name, dataMC), treevar, reg.bins[reg.rvars.index(var)], 1, 1, cuts_mm, '', varname))
                    tmp_ee =   getattr(reg, '%s_ee'%(var))
                    tmp_mm =   getattr(reg, '%s_mm'%(var))
                    tmp_sf =   getattr(reg, '%s_ee'%(var)).Clone(getattr(reg, '%s_ee'%(var)).GetName().replace('ee','sf'))
                    tmp_sf.Add(getattr(reg, '%s_mm'%(var)), 1.)
                    setattr(reg, '%s_sf'%(var), tmp_sf)

                    setattr(reg, '%s_of'%var, tree.getTH1F(lumi, '%s_of_%s%s%s'%(var, eta, reg.name, dataMC), treevar, reg.bins[reg.rvars.index(var)], 1, 1, cuts_of, '', varname))
                    tmp_of = getattr(reg, '%s_of'%(var))

                    tmp_rsfof_histo = make_rsfof(tmp_sf, tmp_of, dataMC)
                    tmp_reeof_histo = make_rsfof(tmp_ee, tmp_of, dataMC)
                    tmp_rmmof_histo = make_rsfof(tmp_mm, tmp_of, dataMC)
                    if var == 'mll': ## this stuff is not needed for any other than mll
                        setattr(reg, "%s_%s_%s_%s"    %("rsfof_yield", dataMC, eta, "of"), reg.mll_of.GetBinContent( reg.mll_of.FindBin(45) ) )
                        setattr(reg, "%s_%s_%s_%s"    %("rsfof_yield", dataMC, eta, "sf"), reg.mll_sf.GetBinContent( reg.mll_sf.FindBin(45) ) )
                        setattr(reg, "%s_%s_%s_%s_err"%("rsfof_yield", dataMC, eta, "of"), reg.mll_of.GetBinError  ( reg.mll_of.FindBin(45) ) )
                        setattr(reg, "%s_%s_%s_%s_err"%("rsfof_yield", dataMC, eta, "sf"), reg.mll_sf.GetBinError  ( reg.mll_sf.FindBin(45) ) )
                        setattr(reg, "%s_%s_%s"    %("rsfof", dataMC, eta), tmp_rsfof_histo.GetBinContent( tmp_rsfof_histo.FindBin(45) ) )
                        setattr(reg, "%s_%s_%s_err"%("rsfof", dataMC, eta), tmp_rsfof_histo.GetBinError  ( tmp_rsfof_histo.FindBin(45) ) )

                    ## set the flavor divided histograms into the region. generic, really
                    setattr(reg, '%s_reeof_%s_%s'%(var, dataMC[:2].lower(), eta), tmp_reeof_histo)
                    setattr(reg, '%s_rmmof_%s_%s'%(var, dataMC[:2].lower(), eta), tmp_rmmof_histo)

                    getattr(reg ,var).setHisto(tmp_rsfof_histo, dataMC, eta)

                    #del reg.mll_sf, reg.mll_of
    

    if opts.mcStudies:
        for var in ['mll', 'nj', 'met']:
            for eta in ['central', 'forward']:
                ## =================
                ## MAKE THE Mll PLOT
                ## =================
                # hard coded values in here.
                val = 1.044 if eta == 'central' else 1.082
                err = math.sqrt( (0.015 if eta == 'central' else 0.023)**2 + (val*0.1)**2 )
                plot_var = Canvas.Canvas('rsfof/%s/plot_rsfof_%s_%s_%s'%(lumi_str, var, eta, 'ttbar'), 'png,pdf', 0.5, 0.2, 0.75, 0.4)
                plot_var.addHisto(getattr(ttjets_inc, var).getHisto('MC'  , eta), 'PE'     , 'ttjets inclusive'  , 'PL', r.kRed+1 , 1, 0)
                plot_var.addBand (getattr(ttjets_inc, var).getHisto('MC', eta).GetXaxis().GetXmin(), val-err, getattr(ttjets_inc, var).getHisto('MC', eta).GetXaxis().GetXmax(), val+err, r.kBlue+2, 0.1)
                plot_var.addLine (getattr(ttjets_inc, var).getHisto('MC', eta).GetXaxis().GetXmin(), val    , getattr(ttjets_inc, var).getHisto('MC', eta).GetXaxis().GetXmax(), val    , r.kBlue+2)
                plot_var.addLatex(0.2, 0.2, eta)
                plot_var.save(1, 1, 0, lumi, 0.2, 1.8)


    else:
        for eta in ['central', 'forward']:
            ## =================
            ## MAKE THE Mll PLOT
            ## =================
            ttjets_meas   .mll.getHisto('MC'  , eta).GetYaxis().SetRangeUser(0.5, 1.5)
            val = ttjets_meas_noZ.mll.getHisto('DATA', eta).GetBinContent(1)
            err = ttjets_meas_noZ.mll.getHisto('DATA', eta).GetBinError(1)
            plot_rsfof = Canvas.Canvas('rsfof/%s/plot_rsfof_mll_%s'%(lumi_str,eta), 'png,pdf', 0.5, 0.2, 0.75, 0.4)
            plot_rsfof.addHisto(ttjets_meas   .mll.getHisto('MC'  , eta), 'PE'     , 'ttjets region - MC'  , 'PL', r.kRed+1 , 1, 0)
            plot_rsfof.addHisto(ttjets_meas   .mll.getHisto('DATA', eta), 'PE,SAME', 'ttjets region - DATA', 'PL', r.kBlack , 1, 1)
            plot_rsfof.addGraph(ttjets_sig_lm .mll.getGraph('MC'  , eta), 'PZ,SAME', 'signal lowmass - MC' , 'PL', r.kBlue-8, 1, 2)
            plot_rsfof.addGraph(ttjets_sig_onZ.mll.getGraph('MC'  , eta), 'PZ,SAME', 'signal onZ - MC'     , 'PL', r.kBlue-7, 1, 3)
            plot_rsfof.addGraph(ttjets_sig_hm .mll.getGraph('MC'  , eta), 'PZ,SAME', 'signal highmass - MC', 'PL', r.kBlue-9, 1, 4)
            #plot_rsfof.addHisto(ttjets_inc .rsfof     , 'E,SAME', 'incl.  region - MC'  , 'PL', r.kBlack , 1, 1)
            plot_rsfof.addBand (ttjets_meas.mll.getHisto('MC', eta).GetXaxis().GetXmin(), val-err, ttjets_meas.mll.getHisto('MC', eta).GetXaxis().GetXmax(), val+err, r.kGreen, 0.2)
            plot_rsfof.addLine (ttjets_meas.mll.getHisto('MC', eta).GetXaxis().GetXmin(), val    , ttjets_meas.mll.getHisto('MC', eta).GetXaxis().GetXmax(), val    , r.kGreen)
            plot_rsfof.addLatex(0.2, 0.2, eta)
            plot_rsfof.save(1, 1, 0, lumi, 0.2, 1.8)
            
        ## =================
        ## PRINT AND SAVE ==
        ## =================
        if saveValues:
            val_mc_cen = ttjets_meas_noZ.mll.getHisto('MC', 'central').GetBinContent(1)
            err_mc_cen = ttjets_meas_noZ.mll.getHisto('MC', 'central').GetBinError(1)
            val_mc_fwd = ttjets_meas_noZ.mll.getHisto('MC', 'forward').GetBinContent(1)
            err_mc_fwd = ttjets_meas_noZ.mll.getHisto('MC', 'forward').GetBinError(1)
            ttjets_meas_noZ.mll.saveInFile(['rsfof', 'direct'], [err_mc_cen/val_mc_cen, err_mc_fwd/val_mc_fwd])


        rsfof_table_fullMass = make_rsfof_table(ttjets_meas_noM)
        rsfof_table_lowMass  = make_rsfof_table(ttjets_meas)
        rsfof_table_noZ      = make_rsfof_table(ttjets_meas_noZ)
