#####################################################################
######                                                              #
###### 88888888888         88                        888888888888   #  
###### 88                  88                                 ,88   #
###### 88                  88                               ,88"    #  
###### 88aaaaa     ,adPPYb,88  ,adPPYb,d8  ,adPPYba,      ,88"      #
###### 88"""""    a8"    `Y88 a8"    `Y88 a8P_____88    ,88"        #
###### 88         8b       88 8b       88 8PP"""""""  ,88"          #
###### 88         "8a,   ,d88 "8a,   ,d88 "8b,   ,aa 88"            #
###### 88888888888 `"8bbdP"Y8  `"YbbdP"Y8  `"Ybbd8"' 888888888888   #
######                       aa,    ,88                             #
######                         "Y8bbdP"                             #
######                                                              #
#####################################################################

import ROOT as r
from   ROOT import gROOT, TCanvas, TFile, TGraphErrors
import math, sys, optparse, array, copy

import include.helper     as helper
import include.Region     as Region
import include.Canvas     as Canvas
import include.CutManager as CutManager
import include.Sample     as Sample
import include.Rounder    as Rounder
import include.Tables     as Tables

def calc_frac(Nin, Nout, Ein, Eout):

    val = [0, 0]
    if(Nin != 0):
        val[0] = Nout/Nin
        val[1] = math.sqrt((Eout * Eout) / (Nin * Nin) + (Ein * Ein) * (Nout*Nout) / (Nin * Nin * Nin * Nin))

    return val


def calc_rinout(NinSF, NinOF, NoutSF, NoutOF, EinSF, EinOF, EoutSF, EoutOF):

    Nin  = NinSF  - NinOF
    Nout = NoutSF - NoutOF
    Ein  = sqrt(EinSF*EinSF   + EinOF*EinOF)
    Eout = sqrt(EoutSF*EoutSF + EoutOF*EoutOF)

    return calc_frac(Nin, Nout, Ein, Eout)


def make_rinout(histo_in_SF, histo_in_OF, histo_out_SF, histo_out_OF):

    histo_in_SF .Add(histo_in_OF , -1.0)
    histo_out_SF.Add(histo_out_OF, -1.0)

    ratio = histo_out_SF.Clone("rinout_" + histo_out_SF.GetName())
    ratio.GetYaxis().SetTitle("r_{out/in}")
    ratio.Divide(histo_in_SF)
  
    return ratio


def make_rinout_histo(histo_SF, histo_OF, scale):

    finalHisto = copy.deepcopy(histo_SF)
    print ' i have %.3f events in the SF sample' %(histo_SF.Integral())
    print ' i have %.3f events in the OF sample' %(histo_OF.Integral())
    finalHisto.Add(histo_OF, -scale)
    if dataMC == 'MC':
        print ' after subtraction i have %.3f events in the SF sample' %(finalHisto.Integral())

    ratio = finalHisto.Clone("rinout_" + finalHisto.GetName())

    # set every bin to the value at the Z
    for _bin in range(1,finalHisto.GetNbinsX()+1):
        finalHisto.SetBinContent(_bin, finalHisto.GetBinContent(finalHisto.FindBin(91.)))
        finalHisto.SetBinError  (_bin, finalHisto.GetBinError  (finalHisto.FindBin(91.)))

    ratio.GetYaxis().SetTitle("r_{out/in}")
    ratio.Divide(finalHisto)
  
    return ratio


if __name__ == "__main__":


    parser = optparse.OptionParser(usage="usage: %prog [opts] FilenameWithSamples", version="%prog 1.0")
    parser.add_option('-s', '--samples', action='store', type=str, dest='sampleFile', default='samples.dat', help='the samples file. default \'samples.dat\'')
    parser.add_option('-i', '--ingredients', action='store', type=str, dest='ingredientFile', default='ingredients.dat', help='the ingredients file. default \'ingredients.dat\'')
    (opts, args) = parser.parse_args()


    plotOnly = False
    if 'plot' in args:
        plotOnly = True

    if not plotOnly:
        print 'Going to load DATA and MC trees...'
        mcDatasets = ['TTLep_pow', 'DYJetsToLL_M10to50', 'DYJetsToLL_M50']
        ## daDatasets = ['DoubleMuon_Run2015C_25ns_05Oct_v1_runs_246908_258714' , 'DoubleEG_Run2015C_25ns_05Oct_v1_runs_246908_258714' , 'MuonEG_Run2015C_25ns_05Oct_v1_runs_246908_258714' ,
        ##               'DoubleMuon_Run2015D_05Oct_v1_runs_246908_258751'      , 'DoubleEG_Run2015D_05Oct_v1_runs_246908_258751'      , 'MuonEG_Run2015D_05Oct_v2_runs_246908_258751'      ,
        ##               'DoubleMuon_Run2015D_v4_runs_246908_258751'            , 'DoubleEG_Run2015D_v4_runs_246908_258751'            , 'MuonEG_Run2015D_v4_runs_246908_258751'            ]
        daDatasets = ['DoubleMuon_Run2015C_25ns-05Oct_v1_runs_246908_260628' , 'DoubleEG_Run2015C_25ns-05Oct_v1_runs_246908_260628' , 'MuonEG_Run2015C_25ns-05Oct_v1_runs_246908_260628' ,
                      'DoubleMuon_Run2015D-05Oct_v1_runs_246908_260628'      , 'DoubleEG_Run2015D-05Oct_v1_runs_246908_260628'      , 'MuonEG_Run2015D-05Oct_v2_runs_246908_260628'      ,
                      'DoubleMuon_Run2015D_v4_runs_246908_260628'            , 'DoubleEG_Run2015D_v4_runs_246908_260628'            , 'MuonEG_Run2015D_v4_runs_246908_260628'            ]
        treeMC = Sample.Tree(helper.selectSamples(opts.sampleFile, mcDatasets, 'MC'), 'MC'  , 0)
        treeDA = Sample.Tree(helper.selectSamples(opts.sampleFile, daDatasets, 'DA'), 'DATA', 1)
        #tree = treeMC
        print 'Trees successfully loaded...'

        ingMC = helper.ingredients(opts.ingredientFile, 'MC'  )
        ingDA = helper.ingredients(opts.ingredientFile, 'DATA')


        gROOT.ProcessLine('.L include/tdrstyle.C')
        gROOT.SetBatch(1)
        r.setTDRStyle()

        cuts = CutManager.CutManager()

        lumi = 2.2
        lumi_str = 'lumi'+str(lumi).replace('.', 'p')+'_forApproval'
 
        regions = []
        dy_nomass = Region.region('DY_nomass',
                           [cuts.DYControlRegion],
                           ['mll', 'met', 'nj'],
                           [ range(20,125,5),
                             range(0,70,10),
                             range(0,8,1) ],
                           True)
        regions.append(dy_nomass)

        ## ===================================================
        ## calculate the rinout values for central and forward
        ## ===================================================

        rinout_meas = Region.region('rinout_meas',
                           [cuts.DYControlRegion],
                           ['mll'],
                           [ [ 20., 70., 81., 101., 120., 210. ] ],
                           True)
         ## ==================================================

        doetas = ['central', 'forward']

        for reg in regions:
            print 'i am at region', reg.name
            for eta in doetas:
                print '... in %s' %(eta)
                etacut = cuts.Central() if eta == 'central' else cuts.Forward()

                thecutsSF_in  = cuts.AddList([cuts.GoodLeptonSF()] + [etacut] + reg.cuts + [cuts.Zmass])
                thecutsOF_in  = cuts.AddList([cuts.GoodLeptonOF()] + [etacut] + reg.cuts + [cuts.Zmass])
                thecutsSF_out = cuts.AddList([cuts.GoodLeptonSF()] + [etacut] + reg.cuts + [cuts.lowmass])
                thecutsOF_out = cuts.AddList([cuts.GoodLeptonOF()] + [etacut] + reg.cuts + [cuts.lowmass])

                for tree in ([treeMC, treeDA] if reg.doData else [treeMC]):
                    dataMC = 'DATA' if tree == treeDA else 'MC'

                    if 'met' in reg.rvars:
                        reg.met_SF_in  = tree.getTH1F(lumi, "met_SF_in" +eta+reg.name+dataMC , "met_pt", reg.bins[reg.rvars.index('met')], 1, 1, thecutsSF_in , "", "E_{T}^{miss.} (GeV)")
                        reg.met_OF_in  = tree.getTH1F(lumi, "met_OF_in" +eta+reg.name+dataMC , "met_pt", reg.bins[reg.rvars.index('met')], 1, 1, thecutsOF_in , "", "E_{T}^{miss.} (GeV)")

                        reg.met_SF_out = tree.getTH1F(lumi, "met_SF_out"+eta+reg.name+dataMC , "met_pt", reg.bins[reg.rvars.index('met')], 1, 1, thecutsSF_out, "", "E_{T}^{miss.} (GeV)")
                        reg.met_OF_out = tree.getTH1F(lumi, "met_OF_out"+eta+reg.name+dataMC , "met_pt", reg.bins[reg.rvars.index('met')], 1, 1, thecutsOF_out, "", "E_{T}^{miss.} (GeV)")
                        reg.met.setHisto(make_rinout(reg.met_SF_in, reg.met_OF_in, reg.met_SF_out, reg.met_OF_out), dataMC, eta)

                    if 'nj' in reg.rvars:
                        reg.nj_SF_in  = tree.getTH1F(lumi, "nj_SF_in" +eta+reg.name+dataMC , "t.nJetSel_Edge", reg.bins[reg.rvars.index('nj')], 1, 1, thecutsSF_in , "", "N_{jets}")
                        reg.nj_OF_in  = tree.getTH1F(lumi, "nj_OF_in" +eta+reg.name+dataMC , "t.nJetSel_Edge", reg.bins[reg.rvars.index('nj')], 1, 1, thecutsOF_in , "", "N_{jets}")

                        reg.nj_SF_out = tree.getTH1F(lumi, "nj_SF_out"+eta+reg.name+dataMC , "t.nJetSel_Edge", reg.bins[reg.rvars.index('nj')], 1, 1, thecutsSF_out, "", "N_{jets}")
                        reg.nj_OF_out = tree.getTH1F(lumi, "nj_OF_out"+eta+reg.name+dataMC , "t.nJetSel_Edge", reg.bins[reg.rvars.index('nj')], 1, 1, thecutsOF_out, "", "N_{jets}")
                        reg.nj.setHisto(make_rinout(reg.nj_SF_in, reg.nj_OF_in, reg.nj_SF_out, reg.nj_OF_out), dataMC, eta)


                    ## do the actual calculus only once
                    if not regions.index(reg):
                        meas_cuts_SF = cuts.AddList([cuts.GoodLeptonSF()] + [etacut] + rinout_meas.cuts)
                        meas_cuts_OF = cuts.AddList([cuts.GoodLeptonOF()] + [etacut] + rinout_meas.cuts)
                        rinout_meas.mll_SF   = tree.getTH1F(lumi, "mll_SF" +eta+rinout_meas.name+dataMC , "t.lepsMll_Edge", rinout_meas.bins[rinout_meas.rvars.index('mll')], 1, 1, meas_cuts_SF , "", "m_{ll} (GeV)")
                        rinout_meas.mll_OF   = tree.getTH1F(lumi, "mll_OF" +eta+rinout_meas.name+dataMC , "t.lepsMll_Edge", rinout_meas.bins[rinout_meas.rvars.index('mll')], 1, 1, meas_cuts_OF , "", "m_{ll} (GeV)")
                        scale = ingDA.rsfof_final_cen_val  if doetas == 'central' else ingDA.rsfof_final_fwd_val
                        s_err = ingDA.rsfof_final_cen_err  if doetas == 'central' else ingDA.rsfof_final_fwd_err
                        if dataMC == 'MC':
                            scale = ingMC.rsfof_final_cen_val  if doetas == 'central' else ingMC.rsfof_final_fwd_val
                            s_err = ingMC.rsfof_final_cen_err  if doetas == 'central' else ingMC.rsfof_final_fwd_err
                        rinout_meas.mll.setHisto(make_rinout_histo(rinout_meas.mll_SF, rinout_meas.mll_OF, scale), dataMC, eta)


        
    for eta in doetas:

        ## ==============================================================
        ## get the mll distribution from the measurement region of rinout
        ## ==============================================================
        etacut = cuts.Central() if eta == 'central' else cuts.Forward()
        thecutsSF_nomass = cuts.AddList([cuts.GoodLeptonSF()] + [etacut] + dy_nomass.cuts)
        thecutsOF_nomass = cuts.AddList([cuts.GoodLeptonOF()] + [etacut] + dy_nomass.cuts)

        dy_nomass.mll_SF_data  = treeDA.getTH1F(lumi, "mll_SF_in" +eta+dy_nomass.name+'data' , "t.lepsMll_Edge", dy_nomass.bins[dy_nomass.rvars.index('mll')], 1, 1, thecutsSF_nomass, "", "m_{ll} (GeV)")
        dy_nomass.mll_OF_data  = treeDA.getTH1F(lumi, "mll_OF_in" +eta+dy_nomass.name+'data' , "t.lepsMll_Edge", dy_nomass.bins[dy_nomass.rvars.index('mll')], 1, 1, thecutsOF_nomass, "", "m_{ll} (GeV)")
        dy_nomass.mll_SF_mc    = treeMC.getTH1F(lumi, "mll_SF_in" +eta+dy_nomass.name+'mc'   , "t.lepsMll_Edge", dy_nomass.bins[dy_nomass.rvars.index('mll')], 1, 1, thecutsSF_nomass, "", "m_{ll} (GeV)")
        dy_nomass.mll_OF_mc    = treeMC.getTH1F(lumi, "mll_OF_in" +eta+dy_nomass.name+'mc'   , "t.lepsMll_Edge", dy_nomass.bins[dy_nomass.rvars.index('mll')], 1, 1, thecutsOF_nomass, "", "m_{ll} (GeV)")

        dy_nomass.mll_SF_mc.GetYaxis().SetRangeUser(1., dy_nomass.mll_SF_mc.GetMaximum()*1.2)
        plot_rinout_mll = Canvas.Canvas("rinout/%s/plot_rinout_mll_%s"%(lumi_str, eta), "png,pdf", 0.2, 0.6, 0.4, 0.8)
        plot_rinout_mll.addHisto(dy_nomass.mll_SF_mc  , "HIST,SAME", "MC - SF", "L" , r.kBlue-4 , 1, 0)
        plot_rinout_mll.addHisto(dy_nomass.mll_OF_mc  , "HIST,SAME", "MC - OF", "L" , r.kRed-4  , 1, 1)
        plot_rinout_mll.addHisto(dy_nomass.mll_SF_data, "PE,SAME", "DATA - SF", "PL", r.kBlue+2 , 1, 2)
        plot_rinout_mll.addHisto(dy_nomass.mll_OF_data, "PE,SAME", "DATA - OF", "PL", r.kRed+2  , 1, 3)
        plot_rinout_mll.addLine ( 81., dy_nomass.mll_SF_mc.GetYaxis().GetXmin(),  81., dy_nomass.mll_SF_mc.GetMaximum(), r.kGreen-3 , 2)
        plot_rinout_mll.addLine (101., dy_nomass.mll_SF_mc.GetYaxis().GetXmin(), 101., dy_nomass.mll_SF_mc.GetMaximum(), r.kGreen-3 , 2)
        plot_rinout_mll.addBand ( 81., dy_nomass.mll_SF_mc.GetYaxis().GetXmin(), 101., dy_nomass.mll_SF_mc.GetMaximum(), r.kGreen-3 , 0.1)
        plot_rinout_mll.addLine ( 70., dy_nomass.mll_SF_mc.GetYaxis().GetXmin(),  70., dy_nomass.mll_SF_mc.GetMaximum(), r.kPink-3, 2)
        plot_rinout_mll.addArrow( 60., 10.,  70., 10., r.kPink-3, "<", 2)
        plot_rinout_mll.addLatex(0.2, 0.2, eta)
        plot_rinout_mll.save(1, 0, 1, lumi)

        ## ============================
        ## get the measured values here
        ## ============================

        rinout_meas.mll.getHisto('MC', eta).GetYaxis().SetRangeUser(0., 0.15)
        meas_histo = rinout_meas.mll.getHisto('DATA', eta)
        meas_value = meas_histo.GetBinContent(meas_histo.FindBin(45.))
        meas_err   = meas_histo.GetBinError  (meas_histo.FindBin(45.))
        full_err = math.sqrt(meas_err**2 + (0.25*meas_value)**2 )
        minus = meas_value - full_err
        plus  = meas_value + full_err
        plot_rinout_meas = Canvas.Canvas("rinout/%s/plot_rinout_measured_%s"%(lumi_str, eta), "png,pdf", 0.6, 0.70, 0.8, 0.80)
        plot_rinout_meas.addHisto(rinout_meas.mll.getHisto('MC'  , eta), "PE,SAME", "MC"       , "PL", r.kRed+1 , 1, 0)
        plot_rinout_meas.addHisto(rinout_meas.mll.getHisto('DATA', eta), "PE,SAME", "data"     , "PL", r.kBlack , 1, 1)
        ##plot_rinout_meas.addBand (rinout_meas.mll.getHisto('MC', eta).GetXaxis().GetXmin(), minus     , rinout_meas.mll.getHisto('MC', eta).GetXaxis().GetXmax(), plus      , r.kGreen, 0.2)
        ##plot_rinout_meas.addLine (rinout_meas.mll.getHisto('MC', eta).GetXaxis().GetXmin(), meas_value, rinout_meas.mll.getHisto('MC', eta).GetXaxis().GetXmax(), meas_value, r.kGreen)
        plot_rinout_meas.addLatex(0.2, 0.2, eta)
        plot_rinout_meas.save(1, 0, 0, lumi)

        table = Tables.makeRinoutTable(rinout_meas)
        ## ===========================
        ## make the depency plots here
        ## ===========================
        
        dy_nomass.met.getHisto('MC'  , eta).GetYaxis().SetRangeUser(0., 0.15)
        plot_rinout_met = Canvas.Canvas("rinout/%s/plot_rinout_met_%s"%(lumi_str,eta), "png,pdf", 0.6, 0.70, 0.8, 0.80)
        plot_rinout_met.addHisto(dy_nomass.met.getHisto('MC'  , eta), "PE,SAME", "MC"       , "PL", r.kRed+1 , 1, 0)
        plot_rinout_met.addHisto(dy_nomass.met.getHisto('DATA', eta), "PE,SAME", "data"     , "PL", r.kBlack , 1, 0)
        plot_rinout_met.addBand (dy_nomass.met.getHisto('MC', eta).GetXaxis().GetXmin(), minus     , dy_nomass.met.getHisto('MC', eta).GetXaxis().GetXmax(), plus      , r.kGreen, 0.2)
        plot_rinout_met.addLine (dy_nomass.met.getHisto('MC', eta).GetXaxis().GetXmin(), meas_value, dy_nomass.met.getHisto('MC', eta).GetXaxis().GetXmax(), meas_value, r.kGreen)
        plot_rinout_met.addLatex(0.2, 0.2, eta)
        plot_rinout_met.save(1, 0, 0, lumi)

        dy_nomass.nj.getHisto('MC'  , eta).GetYaxis().SetRangeUser(0., 0.15)
        plot_rinout_nj = Canvas.Canvas("rinout/%s/plot_rinout_nj_%s"%(lumi_str,eta), "png,pdf", 0.6, 0.70, 0.8, 0.80)
        plot_rinout_nj.addHisto(dy_nomass.nj.getHisto('MC'  , eta), "PE,SAME", "MC"       , "PL", r.kRed+1 , 1, 0)
        plot_rinout_nj.addHisto(dy_nomass.nj.getHisto('DATA', eta), "PE,SAME", "data"     , "PL", r.kBlack , 1, 0)
        plot_rinout_nj.addBand (dy_nomass.nj.getHisto('MC', eta).GetXaxis().GetXmin(), minus     , dy_nomass.nj.getHisto('MC', eta).GetXaxis().GetXmax(), plus      , r.kGreen, 0.2)
        plot_rinout_nj.addLine (dy_nomass.nj.getHisto('MC', eta).GetXaxis().GetXmin(), meas_value, dy_nomass.nj.getHisto('MC', eta).GetXaxis().GetXmax(), meas_value, r.kGreen)
        plot_rinout_nj.addLatex(0.2, 0.2, eta)
        plot_rinout_nj.save(1, 0, 0, lumi)

    ## =================
    ## PRINT AND SAVE ==
    ## =================
    ## rinout_meas.mll.saveInFile(['rinout', 'dy_cr' ], 0.25,  45.)
    ## rinout_meas.mll.saveInFile(['rinout', 'dy_lm' ], 0.25,  45.)
    ## rinout_meas.mll.saveInFile(['rinout', 'dy_bz' ], 0.25,  75.)
    ## rinout_meas.mll.saveInFile(['rinout', 'dy_oz' ], 0.25,  91.)
    ## rinout_meas.mll.saveInFile(['rinout', 'dy_az' ], 0.25, 105.)
    ## rinout_meas.mll.saveInFile(['rinout', 'dy_hm' ], 0.25, 125.)
