############################################################
############################################################
##          _          _            _              _      ##
##         /\ \       /\ \         /\ \           /\ \    ##
##        /  \ \     /  \ \____   /  \ \         /  \ \   ##
##       / /\ \ \   / /\ \_____\ / /\ \_\       / /\ \ \  ##
##      / / /\ \_\ / / /\/___  // / /\/_/      / / /\ \_\ ##
##     / /_/_ \/_// / /   / / // / / ______   / /_/_ \/_/ ##
##    / /____/\  / / /   / / // / / /\_____\ / /____/\    ##
##   / /\____\/ / / /   / / // / /  \/____ // /\____\/    ##
##  / / /______ \ \ \__/ / // / /_____/ / // / /______    ##
## / / /_______\ \ \___\/ // / /______\/ // / /_______\   ##
## \/__________/  \/_____/ \/___________/ \/__________/   ##
############################################################
############################################################
                                                       

import ROOT as r
from ROOT import gROOT, TCanvas, TFile, TF1, TGraphAsymmErrors
import math,sys,optparse
import math



import include.helper     as helper
import include.Region     as Region
import include.Canvas     as Canvas
import include.CutManager as CutManager
import include.Sample     as Sample
import include.Rounder    as rounder


############################################################
def getRT(eff_ee, unc_ee, eff_mm, unc_mm, eff_em, unc_em):
    RT = math.sqrt(eff_ee*eff_mm)/eff_em
    uncRTee = (math.sqrt(eff_mm)/eff_em)*0.5*(unc_ee/math.sqrt(eff_ee)) 
    uncRTmm = (math.sqrt(eff_ee)/eff_em)*0.5*(unc_mm/math.sqrt(eff_mm)) 
    uncRTem = math.sqrt(eff_ee*eff_mm)*((unc_em)/(eff_em*eff_em)) 
    uncRT = math.sqrt(uncRTee * uncRTee + uncRTmm * uncRTmm + uncRTem * uncRTem)
    uncsysRTee = (math.sqrt(eff_mm)/eff_em)*0.5*(0.05*eff_ee/math.sqrt(eff_ee)) 
    uncsysRTmm = (math.sqrt(eff_ee)/eff_em)*0.5*(0.05*eff_mm/math.sqrt(eff_mm)) 
    uncsysRTem = math.sqrt(eff_ee*eff_mm)*((0.05*eff_em)/(eff_em*eff_em)) 
    uncsysRT = math.sqrt(uncsysRTee * uncsysRTee + uncsysRTmm * uncsysRTmm + uncsysRTem * uncsysRTem)

    return [RT, uncRT, uncsysRT]


def RT(eff_ee, eff_mm, eff_em):

    RTratio = eff_ee.Clone("RT_" + eff_ee.GetName())

    RTratio.GetYaxis().SetTitle("RT")
    RTratio.GetXaxis().SetTitle(eff_mm.GetXaxis().GetTitle())

    for i in range(0, eff_mm.GetNbinsX()+1):

        if(eff_em.GetBinContent(i) != 0 and eff_ee.GetBinContent(i) != 0 and eff_mm.GetBinContent(i) != 0):
            [rt, uncrt, uncsys] = getRT(eff_ee.GetBinContent(i), eff_ee.GetBinError(i), eff_mm.GetBinContent(i), eff_mm.GetBinError(i), eff_em.GetBinContent(i), eff_em.GetBinError(i)) 
            RTratio.SetBinContent(i, rt)
            RTratio.SetBinError(i, math.sqrt(uncrt*uncrt+uncsys*uncsys))

    return RTratio


def getTriggerEffs(num, den):

    trig = den.Clone("trig_" + den.GetName())
    trig.Reset()



    errs = TGraphAsymmErrors(num, den, 'v')
    x = errs.GetX()
    y = errs.GetY()
    eyh = errs.GetEYhigh()
    eyl = errs.GetEYlow()
    for i in range(0, num.GetNbinsX()+1):
        if(x[i] != 0 and y[i] != 0):
            trig.SetBinContent(i, y[i])
            trig.SetBinError(i, max(eyh[i], eyl[i]))

    return trig


if __name__ == '__main__':

    parser = optparse.OptionParser(usage='usage: %prog [opts] FilenameWithSamples', version='%prog 1.0')
    parser.add_option('-m', '--mode', action='store', dest='mode', default='rsfof', help='Operation mode')
    parser.add_option('-s', '--samples', action='store', type=str, dest='sampleFile', default='samples.dat', help='the samples file. default \'samples.dat\'')
    (opts, args) = parser.parse_args()

    print 'Going to load DATA and MC trees...'
    mcDatasets = ['TTLep_pow', 'DYJetsToLL_M10to50', 'DYJetsToLL_M50']
    daDatasets = ['HTMHT_Run2015C_25ns-05Oct_v1_runs_246908_260627', 'HTMHT_Run2015D-05Oct_v1_runs_246908_260627', 'HTMHT_Run2015D_v4_runs_246908_260627',
                  'JetHT_Run2015C_25ns-05Oct_v1_runs_246908_260627', 'JetHT_Run2015D-05Oct_v1_runs_246908_260627', 'JetHT_Run2015D_v4_runs_246908_260627']
    treeMC = Sample.Tree(helper.selectSamples(opts.sampleFile, mcDatasets, 'MC'), 'MC'  , 0)
    treeDA = Sample.Tree(helper.selectSamples(opts.sampleFile, daDatasets, 'DA'), 'DATA', 1)
    #tree = treeMC
    print 'Trees successfully loaded...'

    lumi = 2.1
    lumi_str = 'lumi'+str(lumi).replace('.', 'p')
  
    gROOT.ProcessLine('.L include/tdrstyle.C')
    gROOT.SetBatch(1)
    r.setTDRStyle() 
    cuts = CutManager.CutManager()
 
    vetoSignalCR = "(!(" + cuts.METJetsSignalRegion + ")&& !(" + cuts.METJetsControlRegion + "))"
    denominator_EE = cuts.AddList([cuts.GoodLeptonNoTriggeree(), cuts.triggerHT, cuts.HT, vetoSignalCR])
    denominator_MM = cuts.AddList([cuts.GoodLeptonNoTriggermm(), cuts.triggerHT, cuts.HT, vetoSignalCR])
    denominator_EM = cuts.AddList([cuts.GoodLeptonNoTriggerOF(), cuts.triggerHT, cuts.HT, vetoSignalCR])
    numerator_EE = cuts.AddList([denominator_EE, cuts.trigEEc])
    numerator_MM = cuts.AddList([denominator_MM, cuts.trigMMc])
    numerator_EM = cuts.AddList([denominator_EM, cuts.trigEMc])


    regions = []
    mll_inc = Region.region('mll_inc',
                       [''],
                       ['mll'],
                       [[20, 2000]],
                       True)
    regions.append(mll_inc)
    regular = Region.region('regular',
                       ['', '', ''],
                       ['mll', 'l1pt','l2pt'],
                       [[20, 40, 60, 80, 100, 120, 140, 180, 220, 260, 300], [20,40,60,100,140,200,240,280,320], [20,40,60,100,140,200,240,280,320]],
                       True)
    regions.append(regular)


    for reg in regions:
        print 'i am at region', reg.name
        for eta in ['central', 'forward']:
            print '... in %s' %(eta)
            for tree in ([treeMC, treeDA] if reg.doData else [treeMC]):

                dataMC = 'DATA' if tree == treeDA else 'MC'
                cuts_num_ee = cuts.AddList([numerator_EE] + [cuts.Central() if eta == 'central' else cuts.Forward()])
                cuts_den_ee = cuts.AddList([denominator_EE] + [cuts.Central() if eta == 'central' else cuts.Forward()])
                cuts_num_mm = cuts.AddList([numerator_MM] + [cuts.Central() if eta == 'central' else cuts.Forward()])
                cuts_den_mm = cuts.AddList([denominator_MM] + [cuts.Central() if eta == 'central' else cuts.Forward()])
                cuts_num_em = cuts.AddList([numerator_EM] + [cuts.Central() if eta == 'central' else cuts.Forward()])
                cuts_den_em = cuts.AddList([denominator_EM] + [cuts.Central() if eta == 'central' else cuts.Forward()])
                if dataMC == 'MC':
                    cuts_num_ee = cuts.AddList([cuts_num_ee, "(genWeight>0)"])
                    cuts_den_ee = cuts.AddList([cuts_den_ee, "(genWeight>0)"])
                    cuts_num_mm = cuts.AddList([cuts_num_mm, "(genWeight>0)"])
                    cuts_den_mm = cuts.AddList([cuts_den_mm, "(genWeight>0)"])
                    cuts_num_em = cuts.AddList([cuts_num_em, "(genWeight>0)"])
                    cuts_den_em = cuts.AddList([cuts_den_em, "(genWeight>0)"])
                             
 
                for theRvar in reg.rvars:

                    if(theRvar == 'mll'):
                      var = "t.lepsMll_Edge"
                      label = "m_{ll} (GeV)"          
                    if(theRvar == 'l1pt'):
                      var = "t.Lep1_pt_Edge"
                      label = "p_{t,leading} (GeV)"          
                    if(theRvar == 'l2pt'):
                      var = "t.Lep2_pt_Edge"
                      label = "p_{t,trailing} (GeV)"          

 
                    reg.eff_num_ee = tree.getTH1F(lumi, "eff_num_ee"+eta+reg.name+dataMC+theRvar, var, reg.bins[reg.rvars.index(theRvar)], 1, 1, cuts_num_ee, "", label)
                    reg.eff_den_ee = tree.getTH1F(lumi, "eff_den_ee"+eta+reg.name+dataMC+theRvar, var, reg.bins[reg.rvars.index(theRvar)], 1, 1, cuts_den_ee, "", label)
                    reg.eff_num_mm = tree.getTH1F(lumi, "eff_num_mm"+eta+reg.name+dataMC+theRvar, var, reg.bins[reg.rvars.index(theRvar)], 1, 1, cuts_num_mm, "", label)
                    reg.eff_den_mm = tree.getTH1F(lumi, "eff_den_mm"+eta+reg.name+dataMC+theRvar, var, reg.bins[reg.rvars.index(theRvar)], 1, 1, cuts_den_mm, "", label)
                    reg.eff_num_em = tree.getTH1F(lumi, "eff_num_em"+eta+reg.name+dataMC+theRvar, var, reg.bins[reg.rvars.index(theRvar)], 1, 1, cuts_num_em, "", label)
                    reg.eff_den_em = tree.getTH1F(lumi, "eff_den_em"+eta+reg.name+dataMC+theRvar, var, reg.bins[reg.rvars.index(theRvar)], 1, 1, cuts_den_em, "", label)
                    h_eff_ee = getTriggerEffs(reg.eff_num_ee, reg.eff_den_ee)
                    h_eff_mm = getTriggerEffs(reg.eff_num_mm, reg.eff_den_mm)
                    h_eff_em = getTriggerEffs(reg.eff_num_em, reg.eff_den_em)
                    h_RT     = RT(h_eff_ee, h_eff_mm, h_eff_em)
                    setattr(reg, "%s_%s_%s_%s_%s"        %("eff", dataMC, eta, theRvar, "ee"), h_eff_ee.GetBinContent(1))
                    setattr(reg, "%s_%s_%s_%s_%s_err"    %("eff", dataMC, eta, theRvar, "ee"), h_eff_ee.GetBinError(1))
                    setattr(reg, "%s_%s_%s_%s_%s"        %("eff", dataMC, eta, theRvar, "mm"), h_eff_mm.GetBinContent(1))
                    setattr(reg, "%s_%s_%s_%s_%s_err"    %("eff", dataMC, eta, theRvar, "mm"), h_eff_mm.GetBinError(1))
                    setattr(reg, "%s_%s_%s_%s_%s"        %("eff", dataMC, eta, theRvar, "em"), h_eff_em.GetBinContent(1))
                    setattr(reg, "%s_%s_%s_%s_%s_err"    %("eff", dataMC, eta, theRvar, "em"), h_eff_em.GetBinError(1))
                  
                    if(theRvar == 'mll'):
                      reg.mll.setHisto(h_eff_ee, dataMC, eta)
                      reg.mll.setHisto(h_eff_mm, dataMC, eta)
                      reg.mll.setHisto(h_eff_em, dataMC, eta)
                      reg.mll.setHisto(h_RT, dataMC, eta)
                    if(theRvar == 'l1pt'):
                      reg.l1pt.setHisto(h_eff_ee, dataMC, eta)
                      reg.l1pt.setHisto(h_eff_mm, dataMC, eta)
                      reg.l1pt.setHisto(h_eff_em, dataMC, eta)
                      reg.l1pt.setHisto(h_RT, dataMC, eta)
                    if(theRvar == 'l2pt'):
                      reg.l2pt.setHisto(h_eff_ee, dataMC, eta)
                      reg.l2pt.setHisto(h_eff_mm, dataMC, eta)
                      reg.l2pt.setHisto(h_eff_em, dataMC, eta)
                      reg.l2pt.setHisto(h_RT, dataMC, eta)
 
                    if(len(reg.bins[reg.rvars.index(theRvar)]) > 2):
                    	plot_rt = Canvas.Canvas("rt/%s/plot_rt_%s_%s_%s"%(lumi_str, theRvar, eta, dataMC), "png,pdf", 0.6, 0.15, 0.8, 0.35)
                        h_RT.GetYaxis().SetRangeUser(0.5, 1.5)
                        plot_rt.addHisto(h_RT, "E", ""       , "", r.kRed+1 , 1, -1)
                        plot_rt.save(1, 0, 0, lumi)
                    	plot_eff = Canvas.Canvas("rt/%s/plot_eff_%s_%s_%s"%(lumi_str, theRvar, eta, dataMC), "png,pdf", 0.6, 0.15, 0.8, 0.35)
                        h_eff_ee.GetYaxis().SetRangeUser(0.5, 1.1)
                        h_eff_mm.GetYaxis().SetRangeUser(0.5, 1.1)
                        h_eff_em.GetYaxis().SetRangeUser(0.5, 1.1)
                        plot_eff.addHisto(h_eff_ee, "E", "DoubleElectron"       , "PL", r.kRed+1 , 1, 0)
                        plot_eff.addHisto(h_eff_mm, "E,SAME", "DoubleMuon"       , "PL", r.kGreen+1 , 1, 0)
                        plot_eff.addHisto(h_eff_em, "E,SAME", "MuonElectron"       , "PL", r.kBlue+1 , 1, 0)
                        plot_eff.save(1, 0, 0, lumi)



    mll_inc.mll.saveInFile(['rt', 'region'], 0, 1)



