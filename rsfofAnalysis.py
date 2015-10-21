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

    # f = open('txts/'+ratio.GetName()+'_values.txt', 'w')
    # for i in range(1, ratio.GetNbinsX()+1):
    #     min, max = ratio.GetBinLowEdge(i), ratio.GetBinLowEdge(i)+ratio.GetBinWidth(i)
    #     print    '%10s : R_SFOF in [%.2f, %.2f] GeV:\t%.3f +- %.3f'    %(dataMC, min, max, ratio.GetBinContent(i), ratio.GetBinError(i) )
    #     f.write( '%10s : R_SFOF in [%.2f, %.2f] GeV:\t%.3f +- %.3f \n' %(dataMC, min, max, ratio.GetBinContent(i), ratio.GetBinError(i) ) )
    # f.close()


    return ratio



if __name__ == '__main__':

    print 'Starting r_SFOF analysis...'
    parser = optparse.OptionParser(usage='usage: %prog [options] FilenameWithSamples', version='%prog 1.0')
    parser.add_option('-m', '--mode', action='store', dest='mode', default='rsfof', help='Operation mode')
    parser.add_option('-t', '--trigger', action='store', type='int', dest='triggerFlag', default='1', help='Trigger cut. Set to 0 if you want to run without trigger')
    (options, args) = parser.parse_args()


    if len(args) != 1:
      parser.error('wrong number of arguments')

    inputFileName = args[0]

    print 'Going to load DATA and MC trees...'
    mcDatasets = ['TTLep_pow', 'DYJetsToLL_M10to50', 'DYJetsToLL_M50']
    daDatasets = ['DoubleMuon_Run2015D_05Oct_v1_runs_246908_258751', 'DoubleEG_Run2015D_05Oct_v1_runs_246908_258751', 'MuonEG_Run2015D_05Oct_v2_runs_246908_258751',
                  'DoubleMuon_Run2015D_v4_runs_246908_258751', 'DoubleEG_Run2015D_v4_runs_246908_258751', 'MuonEG_Run2015D_v4_runs_246908_258751']
    treeMC = Sample.Tree(helper.selectSamples(inputFileName, mcDatasets, 'MC'), 'MC'  , 0)
    treeDA = Sample.Tree(helper.selectSamples(inputFileName, daDatasets, 'DA'), 'DATA', 1)
    print 'Trees successfully loaded...'


    ##lumi = 0.592
    ## lumi = 1.28
    lumi = 0.58
    maxrun = 258159
    lumi_str = 'lumi'+str(lumi).replace('.', 'p')
    print 'Running with an integrated luminosity of', lumi,'fb-1'
   
    gROOT.ProcessLine('.L tdrstyle.C')
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
                                   [cuts.METJetsControlRegion, cuts.ZmassVeto],
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
                                   ['mll'],
                                   [[20, 70, 81, 101, 120, 300]],
                                   False)
    regions.append(ttjets_inc)

    for reg in regions:
        print 'i am at region', reg.name
        for eta in ['central', 'forward']:
            print '... in %s' %(eta)

            cuts_sf = cuts.AddList([cuts.MaxRun(maxrun), cuts.GoodLeptonSF(),cuts.trigger]+[cuts.Central() if eta == 'central' else cuts.Forward()]+reg.cuts)
            cuts_of = cuts.AddList([cuts.MaxRun(maxrun), cuts.GoodLeptonOF(),cuts.trigger]+[cuts.Central() if eta == 'central' else cuts.Forward()]+reg.cuts)

            for tree in ([treeMC, treeDA] if reg.doData else [treeMC]):

                dataMC = 'DATA' if tree == treeDA else 'MC'
                reg.mll_sf = tree.getTH1F(lumi, 'mll_sf_'+eta+reg.name+dataMC, 't.lepsMll_Edge', reg.bins[reg.rvars.index('mll')], 1, 1, cuts_sf, '', "m_{ll} (GeV)")
                reg.mll_of = tree.getTH1F(lumi, 'mll_of_'+eta+reg.name+dataMC, 't.lepsMll_Edge', reg.bins[reg.rvars.index('mll')], 1, 1, cuts_of, '', "m_{ll} (GeV)")

                tmp_rsfof_histo = make_rsfof(reg.mll_sf, reg.mll_of, dataMC)
                setattr(reg, "%s_%s_%s_%s"    %("rsfof_yield", dataMC, eta, "of"), reg.mll_of.GetBinContent( reg.mll_of.FindBin(45) ) )
                setattr(reg, "%s_%s_%s_%s"    %("rsfof_yield", dataMC, eta, "sf"), reg.mll_sf.GetBinContent( reg.mll_sf.FindBin(45) ) )
                setattr(reg, "%s_%s_%s_%s_err"%("rsfof_yield", dataMC, eta, "of"), reg.mll_of.GetBinError  ( reg.mll_of.FindBin(45) ) )
                setattr(reg, "%s_%s_%s_%s_err"%("rsfof_yield", dataMC, eta, "sf"), reg.mll_sf.GetBinError  ( reg.mll_sf.FindBin(45) ) )
                setattr(reg, "%s_%s_%s"    %("rsfof", dataMC, eta), tmp_rsfof_histo.GetBinContent( tmp_rsfof_histo.FindBin(45) ) )
                setattr(reg, "%s_%s_%s_err"%("rsfof", dataMC, eta), tmp_rsfof_histo.GetBinError  ( tmp_rsfof_histo.FindBin(45) ) )

                reg.mll.setHisto(tmp_rsfof_histo, dataMC, eta)

                #del reg.mll_sf, reg.mll_of
    

    for eta in ['central', 'forward']:
        ## =================
        ## MAKE THE Mll PLOT
        ## =================
        ttjets_meas   .mll.getHisto('MC'  , eta).GetYaxis().SetRangeUser(0.5, 1.5)
        plot_rsfof = Canvas.Canvas('rsfof/%s/plot_rsfof_mll_%s'%(lumi_str,eta), 'png,pdf', 0.5, 0.2, 0.75, 0.4)
        plot_rsfof.addHisto(ttjets_meas   .mll.getHisto('MC'  , eta), 'PE'     , 'ttjets region - MC'  , 'PL', r.kRed+1 , 1, 0)
        plot_rsfof.addHisto(ttjets_meas   .mll.getHisto('DATA', eta), 'PE,SAME', 'ttjets region - DATA', 'PL', r.kBlack , 1, 1)
        plot_rsfof.addGraph(ttjets_sig_lm .mll.getGraph('MC'  , eta), 'PZ,SAME', 'signal lowmass - MC' , 'PL', r.kBlue-8, 1, 2)
        plot_rsfof.addGraph(ttjets_sig_onZ.mll.getGraph('MC'  , eta), 'PZ,SAME', 'signal onZ - MC'     , 'PL', r.kBlue-7, 1, 3)
        plot_rsfof.addGraph(ttjets_sig_hm .mll.getGraph('MC'  , eta), 'PZ,SAME', 'signal highmass - MC', 'PL', r.kBlue-9, 1, 4)
        #plot_rsfof.addHisto(ttjets_inc .rsfof     , 'E,SAME', 'incl.  region - MC'  , 'PL', r.kBlack , 1, 1)
        plot_rsfof.addBand (ttjets_meas.mll.getHisto('MC', eta).GetXaxis().GetXmin(), 0.9, ttjets_meas.mll.getHisto('MC', eta).GetXaxis().GetXmax(), 1.1, r.kGreen, 0.2)
        plot_rsfof.addLine (ttjets_meas.mll.getHisto('MC', eta).GetXaxis().GetXmin(), 1.0, ttjets_meas.mll.getHisto('MC', eta).GetXaxis().GetXmax(), 1.0, r.kGreen)
        plot_rsfof.addLatex(0.2, 0.2, eta)
        plot_rsfof.save(1, 1, 0, lumi)
        
    ## =================
    ## PRINT AND SAVE ==
    ## =================
    ttjets_meas   .mll.saveInFile(['rsfof', 'ttcr_lm' ], 0.1,  50)
    ttjets_meas   .mll.saveInFile(['rsfof', 'ttcr_onZ'], 0.1,  91)
    ttjets_meas   .mll.saveInFile(['rsfof', 'ttcr_hm' ], 0.1, 150)
    ttjets_sig_lm .mll.saveInFile(['rsfof', 'sr_lm'   ], 0.1)
    ttjets_sig_onZ.mll.saveInFile(['rsfof', 'sr_onZ'  ], 0.1)
    ttjets_sig_hm .mll.saveInFile(['rsfof', 'sr_hm'   ], 0.1)


    rsfof_table_fullMass = make_rsfof_table(ttjets_meas_noM)
    rsfof_table_lowMass  = make_rsfof_table(ttjets_meas)
    rsfof_table_noZ      = make_rsfof_table(ttjets_meas_noZ)
