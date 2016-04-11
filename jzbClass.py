import ROOT as r
from   ROOT import gROOT, TCanvas, TFile, TF1, TH1F, TGraph
import math, sys, optparse, copy, re, array, itertools


import include.helper     as helper
import include.Region     as Region
import include.Canvas     as Canvas
import include.CutManager as CutManager
import include.Sample     as Sample
import include.Tables     as Tables

class jzbAnalysis:
    def __init__(self, name):
        self.name = name
        self.loadIngs()

    def loadIngs(self):
        self.ingMC = helper.ingredients(opts.ingredientFile, 'MC'  )
        self.ingDA = helper.ingredients(opts.ingredientFile, 'DATA')

    def loadTrees(self):
        print 'loading trees'
        ttDatasets = ['TTLep_pow']
#        dyDatasetsHT = ['DYJetsToLL_M50_HT100to200', 'DYJetsToLL_M50_HT200to400',
#                        'DYJetsToLL_M50_HT400to600', 'DYJetsToLL_M50_HT600toInf']
#        mcDatasetsHT = ['DYJetsToLL_M50_HT100to200', 'DYJetsToLL_M50_HT200to400',
#                        'DYJetsToLL_M50_HT400to600', 'DYJetsToLL_M50_HT600toInf','TTLep_pow',
#                        'WWTo2L2Nu', 'WZTo2L2Q', 'ZZTo2L2Q','TToLeptons_tch_amcatnlo',
#                        'TTZToLLNuNu', 'WZTo3L1Nu','TBar_tWch','T_tWch','VHToNobb_M125', 'TTHToNobb_M125']
#        mcDatasetsHT = ['DYJetsToLL_M50_HT100to200', 'DYJetsToLL_M50_HT200to400',
#                        'DYJetsToLL_M50_HT400to600', 'DYJetsToLL_M50_HT600toInf','TTLep_pow']
#        dyDatasets   = ['DYJetsToLL_M50']
#        mcDatasets   = [ 'DYJetsToLL_M50','TTLep_pow', 'WWTo2L2Nu', 'WZTo2L2Q', 'ZZTo2L2Q','TToLeptons_tch_amcatnlo',
#                        'TTZToLLNuNu', 'WZTo3L1Nu','TBar_tWch','T_tWch','VHToNobb_M125', 'TTHToNobb_M125']
        mcDatasets   = [ 'DYJetsToLL_M50','TTLep_pow']
        siDatasets   = ['SMS_T5ZZ_mGluino1000To1250_mLSP100To1000', 'SMS_T5ZZ_mGluino1200To1350_mLSP100To1200',
                        'SMS_T5ZZ_mGluino1400To1550_mLSP100To1400', 'SMS_T5ZZ_mGluino600To700_mLSP100To500',
                        'SMS_T5ZZ_mGluino800To950_mLSP100To800']

#        daDatasets = ['MuonEG_Run2015C_25ns-05Oct_v1_runs_246908_260628',
#                      #'HTMHT_Run2015C_25ns-05Oct_v1_runs_246908_260628',
#                      #'DoubleEG_Run2015C_25ns-05Oct_v1_runs_246908_260628',
#                      #'JetHT_Run2015D-05Oct_v1_runs_246908_260628',
#                      #'DoubleEG_Run2015D-05Oct_v1_runs_246908_260628',
#                      'MuonEG_Run2015D-05Oct_v2_runs_246908_260628',
#                      #'HTMHT_Run2015D_v4_runs_246908_260628',
#                      #'DoubleMuon_Run2015D-05Oct_v1_runs_246908_260628',
#                      #'DoubleMuon_Run2015D_v4_runs_246908_260628',
#                      #'HTMHT_Run2015D-05Oct_v1_runs_246908_260628',
#                      #'JetHT_Run2015D_v4_runs_246908_260628',
#                      #'DoubleEG_Run2015D_v4_runs_246908_260628',
#                      #'JetHT_Run2015C_25ns-05Oct_v1_runs_246908_260628',
#                      #'DoubleMuon_Run2015C_25ns-05Oct_v1_runs_246908_260628',
#                      'MuonEG_Run2015D_v4_runs_246908_260628']
        self.treeTT   = Sample.Tree(helper.selectSamples(opts.sampleFile, ttDatasets,   'TT'), 'TT', 0, isScan = 0)
#        self.treeDA   = Sample.Tree(helper.selectSamples(opts.sampleFile, daDatasets,   'DA'), 'DA', 1, isScan = 0)
#        self.treeDY   = Sample.Tree(helper.selectSamples(opts.sampleFile, dyDatasets,   'DY'), 'DY', 0, isScan = 0)
        self.treeMC   = Sample.Tree(helper.selectSamples(opts.sampleFile, mcDatasets,   'MC'), 'MC', 0, isScan = 0)
#        self.treeDYHT = Sample.Tree(helper.selectSamples(opts.sampleFile, dyDatasetsHT, 'DY'), 'DY', 0, isScan = 0)
#        self.treeMCHT = Sample.Tree(helper.selectSamples(opts.sampleFile, mcDatasetsHT, 'MC'), 'MC', 0, isScan = 0)

        self.treeSI = Sample.Tree(helper.selectSamples(opts.sampleFile, siDatasets, 'SI'), 'SI', 0, isScan = 1)

    def loadStack(self):
        print 'loading stacks'
        self.loadTrees()
        cuts = CutManager.CutManager()
        for key in ['2jets','3jets']:
            #for ht in ['','HT']:
            for ht in ['']:
                for fl in ['OF','SF']:
                    for eta in ['','Central','Forward']:
                        setattr(self,'stackMC%s%s%s%s'%(key,ht,fl,eta),copy.deepcopy(getattr(self,'treeMC%s'%ht).getStack(lumi, "jzb_%s_%s_%s%s"%(key,ht,fl,eta), "t.lepsJZB_Edge", 50, -300., 300, cuts.AddList([getattr(cuts,'GoodLepton%s'%fl)(), 't.nJetSel_Edge == 2' if key == '2jets' else 't.nJetSel_Edge > 2', cuts.Zmass, '' if eta == '' else getattr(cuts,eta)()]), "", "JZB [GeV]")))


    def loadHistos(self):
        print 'loading histos'
        cuts = CutManager.CutManager()
        self.loadTrees()
        for key in ['2jets','3jets']:
            #for ht in ['','HT']:
            for ht in ['']:
                for fl in ['SF','OF']:
                    for eta in ['Central','Forward']:
                        setattr(self,'histoDY%s%s%s%s'%(key,ht,fl,eta),copy.deepcopy(getattr(self,'treeDY%s'%ht).getTH1F(lumi,'dy%s%s%s%s'%(key,ht,fl,eta),"t.lepsJZB_Edge", 50, -300., 300, cuts.AddList([getattr(cuts,'GoodLepton%s'%fl)(), 't.nJetSel_Edge == 2' if key == '2jets' else 't.nJetSel_Edge > 2', cuts.Zmass,'' if eta == '' else getattr(cuts,eta)()]), "", "JZB [GeV]")))
                        setattr(self,'center%s%s%s'%(key,ht,eta),self.getCentering(getattr(self,'histoDY%s%s%s%s'%(key,ht,'SF',eta))))
                        
                        setattr(self,'histoDY%s%s%s%s_center'%(key,ht,fl,eta),copy.deepcopy(getattr(self,'treeDY%s'%ht).getTH1F(lumi,'dy%s%s%s%s_center'%(key,ht,fl,eta),"t.lepsJZB_Edge - %f"%getattr(self,'center%s%s%s'%(key,ht,eta)), 50, -300., 300, cuts.AddList([getattr(cuts,'GoodLepton%s'%fl)(), 't.nJetSel_Edge == 2' if key == '2jets' else 't.nJetSel_Edge > 2', cuts.Zmass,'' if eta == '' else getattr(cuts,eta)()]), "", "JZB [GeV]")))

    def loadClosure(self):
        print 'loading histos'
        self.loadTrees()
        cuts = CutManager.CutManager()
        for key in ['2jets','3jets']:
            #for ht in ['','HT']:
            for ht in ['']:
                for fl in ['OF','SF']: 
                    for eta in ['Central','Forward']:
                        setattr(self,'histoDY%s%s%s%s_pos'%(key,ht,fl,eta),copy.deepcopy(getattr(self,'treeDY%s'%ht).getTH1F(lumi,'dy%s%s%s%s_pos'%(key,ht,fl,eta),"t.lepsJZB_Edge - %f"%(getattr(self,'center%s%s%s'%(key,ht,eta))), 25, 0., 200, cuts.AddList([getattr(cuts,'GoodLepton%s'%fl)(), 't.nJetSel_Edge == 2' if key == '2jets' else 't.nJetSel_Edge > 2', cuts.Zmass, "t.lepsJZB_Edge - %f > 0"%(getattr(self,'center%s%s%s'%(key,ht,eta))),'' if eta == '' else getattr(cuts,eta)() ]), "", "JZB [GeV]"))) 
                        setattr(self,'histoDY%s%s%s%s_neg'%(key,ht,fl,eta),copy.deepcopy(getattr(self,'treeDY%s'%ht).getTH1F(lumi,'dy%s%s%s%s_neg'%(key,ht,fl,eta),"TMath::Abs(t.lepsJZB_Edge - %f)"%getattr(self,'center%s%s%s'%(key,ht,eta)), 25, 0., 200, cuts.AddList([getattr(cuts,'GoodLepton%s'%fl)(), 't.nJetSel_Edge == 2' if key == '2jets' else 't.nJetSel_Edge > 2', cuts.Zmass, "t.lepsJZB_Edge - %f < 0"%getattr(self,'center%s%s%s'%(key,ht,eta)),'' if eta == '' else getattr(cuts,eta)() ]), "", "JZB [GeV]"))) 
                        

    def plotStack(self):
        self.loadStack()
        for key in ['2jets','3jets']:
            #for ht in ['','HT']:
            for ht in ['']:
                for fl in ['OF','SF']:
                    for eta in ['','Central','Forward']:
                        plot = Canvas.Canvas('finalPlots/stack_%s_%s_%s_%s'%(key,ht,fl,eta),'png,pdf,cxx', 0.3, 0.60, 0.5, 0.8)  
                        plot.addStack(getattr(self,'stackMC%s%s%s%s'%(key,ht,fl,eta)),'HIST',True,0,['t#bar{t}','Drell-Yan','Rare'])
                        plot.save(True,False,True,lumi,0.,10*getattr(self,'stackMC%s%s%s%s'%(key,ht,fl,eta)).GetMaximum())
    def plotHistos(self):
        self.loadHistos()
        for key in ['2jets','3jets']:
            #for ht in ['','HT']:
            for ht in ['']:
                for eta in ['Central','Forward']:
                    plot = Canvas.Canvas('finalPlots/jzb_%s_%s_%s'%(key,ht,eta),'png,pdf,cxx', 0.4, 0.20, 0.6, 0.32)
                    SF   = getattr(self,'histoDY%s%s%s%s'%(key,ht,'SF',eta))
                    print key, ht,eta, SF.Integral()
                    OF   = getattr(self,'histoDY%s%s%s%s'%(key,ht,'OF',eta))
                    SF.Add(OF,-1.*self.ingMC.rsfof_final_cen_val if eta == 'Central' else -1.*self.ingMC.rsfof_final_fwd_val)
                    SF_center   = getattr(self,'histoDY%s%s%s%s_center'%(key,ht,'SF',eta))
                    OF_center   = getattr(self,'histoDY%s%s%s%s_center'%(key,ht,'OF',eta))
                    SF_center.Add(OF_center,-1.)
                    plot.addHisto(SF, "E,SAME","Drell-Yan","P", r.kRed, True,0)
                    plot.addHisto(SF_center, "E,SAME","Drell-Yan (corrected)","P", r.kBlue, True,1)
                    plot.save(True,False,True,lumi)
                plot = Canvas.Canvas('finalPlots/jzb_%s_%s_%s'%(key,ht,eta),'png,pdf,cxx', 0.4, 0.20, 0.6, 0.32)
                SF        = getattr(self,'histoDY%s%s%s%s'%(key,ht,'SF','Central')).Clone()
                SF.     Add(getattr(self,'histoDY%s%s%s%s'%(key,ht,'SF','Forward')))
                SF_center = getattr(self,'histoDY%s%s%s%s_center'%(key,ht,'SF','Central')).Clone()
                SF_center.Add(getattr(self,'histoDY%s%s%s%s_center'%(key,ht,'SF','Forward')))
                plot.addHisto(SF, "E,SAME","Drell-Yan","P", r.kRed, True,0)
                plot.addHisto(SF_center, "E,SAME","Drell-Yan (corrected)","P", r.kBlue, True,1)
                plot.save(True,False,True,lumi)

    def plotClosure(self):
        self.loadClosure()
        for key in ['2jets','3jets']:
            #for ht in ['','HT']:
            for ht in ['']:
                for eta in ['Central','Forward']:
                    plot = Canvas.Canvas('finalPlots/closure_%s_%s_%s'%(key,ht,eta),'png,pdf,cxx', 0.65, 0.50, 0.90, 0.72)
                    SF_pos = getattr(self,'histoDY%s%s%s%s_pos'%(key,ht,'SF',eta))
                    OF_pos = getattr(self,'histoDY%s%s%s%s_pos'%(key,ht,'OF',eta))
                    SF_pos.Add(OF_pos,-1.*self.ingMC.rsfof_final_cen_val if eta == 'Central' else -1.*self.ingMC.rsfof_final_fwd_val)
                    SF_neg = getattr(self,'histoDY%s%s%s%s_neg'%(key,ht,'SF',eta))
                    OF_neg = getattr(self,'histoDY%s%s%s%s_neg'%(key,ht,'OF',eta))
                    SF_neg.Add(OF_neg,-1.*self.ingMC.rsfof_final_cen_val if eta == 'Central' else -1.*self.ingMC.rsfof_final_fwd_val)
                    plot.addHisto(SF_pos,"E,SAME", "JZB > 0","P",r.kRed,True,0)
                    plot.addHisto(SF_neg,"HIST,SAME", "JZB < 0","L",r.kBlue,True,1)
                    plot.addBandRatio( 0., 0.7, 250., 1.3, r.kRed, 0.4)
                    plot.saveRatio(True,False,True,lumi,SF_pos,SF_neg)

                plot = Canvas.Canvas('finalPlots/closure_%s_%s'%(key,ht),'png,pdf,cxx', 0.65, 0.50, 0.90, 0.72) 
                SF_pos = getattr(self,'histoDY%s%s%s%s_pos'%(key,ht,'SF','Central')).Clone()
                SF_pos.Add(getattr(self,'histoDY%s%s%s%s_pos'%(key,ht,'SF','Forward')))
                SF_neg = getattr(self,'histoDY%s%s%s%s_neg'%(key,ht,'SF','Central')).Clone()
                SF_neg.Add(getattr(self,'histoDY%s%s%s%s_neg'%(key,ht,'SF','Forward')))
                plot.addHisto(SF_pos,"E,SAME", "JZB > 0","P",r.kRed,True,0)
                plot.addHisto(SF_neg,"HIST,SAME", "JZB < 0","L",r.kBlue,True,1)
                plot.addBandRatio( 0., 0.7, 250., 1.3, r.kRed, 0.4)
                plot.saveRatio(True,False,True,lumi,SF_pos,SF_neg)

    def getCentering(self,histo):
        histo.Fit('gaus','0Q','',-30.,30.)
        return histo.GetFunction('gaus').GetParameter(1)
    
    def loadEstimation(self,jzbCut):
        self.loadHistos()
        cuts = CutManager.CutManager()
  
        for key in ['2jets','3jets']:
            #for ht in ['','HT']:
            for ht in ['']:
                for eta in ['Central','Forward']:
                    for fl in ['OF','SF']:
                        for sig in ['pos','neg']:
                            #setattr(self, 'mll%s%s%s%s%s'%(key,ht,eta,fl,sig), getattr(self,'treeMC%s'%ht).getTH1F(lumi,'mll%s%s%s%s%s'%(key,ht,eta,fl,sig),"t.lepsMll_Edge", 50, 0., 200, cuts.AddList([getattr(cuts,'GoodLepton%s'%fl)(), 't.nJetSel_Edge == 2' if key == '2jets' else 't.nJetSel_Edge > 2', cuts.Zmass, "t.lepsJZB_Edge - %f %s%f"%(getattr(self,'center%s%s%s'%(key,ht,eta)),'> ' if sig == 'pos' else '< -',jzbCut),'' if eta == '' else getattr(cuts,eta)() ]), "", "M_{ll} [GeV]"))
                            setattr(self, 'jzb%s%s%s%s%s'%(key,ht,eta,fl,sig), getattr(self,'treeMC%s'%ht).getTH1F(lumi,'jzb%s%s%s%s%s'%(key,ht,eta,fl,sig),"TMath::Abs(t.lepsJZB_Edge-%f)"%(getattr(self,'center%s%s%s'%(key,ht,eta))), 50, 0., 200, cuts.AddList([getattr(cuts,'GoodLepton%s'%fl)(), 't.nJetSel_Edge == 2' if key == '2jets' else 't.nJetSel_Edge > 2', cuts.Zmass, "t.lepsJZB_Edge - %f %s 0"%(getattr(self,'center%s%s%s'%(key,ht,eta)),'> ' if sig == 'pos' else '< -'),'' if eta == '' else getattr(cuts,eta)() ]), "", "JZB [GeV]"))

    def loadDataCard(self,jzbCut):
        self.loadHistos()
        cuts = CutManager.CutManager()
  
        for key in ['2jets','3jets']:
            for eta in ['Central','Forward']:
                for bjet in ['inc','0','1','2']:
                    for fl in ['OF','SF']:
                        for sig in ['pos','neg']:
                            if 'inc' in bjet: bcut = ''
                            else:             bcut = 'nBJetMedium35_Edge %s %s'%('==' if bjet == '0' else '>=', bjet)
                            print 'a'
#setattr(self, 'datacard%s%s%s%s%s'%(key,eta,bjet,fl,sig), getattr(self,'treeMC').getTH1F(lumi,'datacard%s%s%s%s%s'%(key,eta,bjet,fl,sig),"t.lepsMll_Edge", [20., 70., 81., 101., 120., 13000.],0,0 cuts.AddList([getattr(cuts,'GoodLepton%s'%fl)(), 't.nJetSel_Edge == 2' if key == '2jets' else 't.nJetSel_Edge > 2', cuts.Zmass, "t.lepsJZB_Edge - %f %s%f"%(getattr(self,'center%s%s%s'%(key,ht,eta)),'> ' if sig == 'pos' else '< -',jzbCut),'' if eta == '' else getattr(cuts,eta)() , bcut]), "", "M_{ll} [GeV]"))
                            
    def makeEstimation(self,jzbCut):
        fil = open('JZByields.txt','w')
        self.loadEstimation(jzbCut)
        print 'Getting yields'
        fil.write('            %s      %s  %s   %s\n'%('SR','FS Bkg','Drell-Yan','Total'))
        for key in ['2jets','3jets']:
            #for ht in ['','HT']:
            for ht in ['']:
                 for eta in ['Central','Forward']:
                    for obs in ['mll','jzb']:
#                     for obs in ['jzb']:
                        SF_pos    = getattr(self,'%s%s%s%s%s%s'%(obs,key,ht,eta,'SF','pos'))
                        SF_neg    = getattr(self,'%s%s%s%s%s%s'%(obs,key,ht,eta,'SF','neg'))
                        OF_pos    = getattr(self,'%s%s%s%s%s%s'%(obs,key,ht,eta,'OF','pos'))
                        OF_neg    = getattr(self,'%s%s%s%s%s%s'%(obs,key,ht,eta,'OF','neg'))
                        bkg_stack = r.THStack()
                        OF_pos.Scale(self.ingMC.rsfof_final_cen_val if eta == 'Central' else self.ingMC.rsfof_final_fwd_val)                    
                        for _bin in range(1,OF_pos.GetNbinsX()+1):
                            errStat = OF_pos.GetBinError(_bin)
                            errSyst = (self.ingMC.rsfof_final_cen_err if eta == 'Central' else self.ingMC.rsfof_final_fwd_err)*OF_pos.GetBinContent(_bin)
                            OF_pos.SetBinError(_bin, math.sqrt(errStat*errStat + errSyst*errSyst))
                        OF_neg.Scale(self.ingMC.rsfof_final_cen_val if eta == 'Central' else self.ingMC.rsfof_final_fwd_val)
                        SF_neg.Add(OF_neg,-1.)
                        for _bin in range(1,SF_neg.GetNbinsX()+1):
                            errStat = SF_neg.GetBinError(_bin)
                            errSyst = 0.2*SF_neg.GetBinContent(_bin)
                            SF_neg.SetBinError(_bin,math.sqrt(errStat*errStat + errSyst*errSyst))
                        OF_pos.SetFillColor(r.kWhite)
                        SF_neg.SetFillColor(r.kBlue)
                        bkg_stack.Add(OF_pos)
                        bkg_stack.Add(SF_neg)
                        bkg_stack.Draw()
                        totalBkg = OF_pos.Clone()
                        totalBkg.Add(SF_neg)
                        plot = Canvas.Canvas('finalPlots/signalregions/estimation_%s_%s_%s_%s'%(key,ht,eta,obs),'png,pdf,cxx', 0.65, 0.50, 0.90, 0.72)
                        plot.addStack(bkg_stack,'HIST',True,0,['FS background','Drell-Yan'])
                        plot.addHisto(SF_pos,"E,SAME", "Expected","P",r.kBlack,True,0)
                        plot.saveRatio(True,False,True,lumi,SF_pos,totalBkg)
                        if obs == 'jzb':
                            fil.write('%s %s %s\n' %(key,ht,eta))
                            for jzb in [100,150,200]:
                                errSF_pos = r.Double(0.); errOF_pos = r.Double(0.); errSF_neg = r.Double(0.); errTlb = r.Double(0.)
                                SFpos = SF_pos.IntegralAndError(SF_pos.FindBin(jzb),SF_pos.GetNbinsX(),errSF_pos)
                                OFpos = OF_pos.IntegralAndError(OF_pos.FindBin(jzb),OF_pos.GetNbinsX(),errOF_pos)
                                SFneg = SF_neg.IntegralAndError(SF_neg.FindBin(jzb),SF_neg.GetNbinsX(),errSF_neg)
                                totlb = totalBkg.IntegralAndError(totalBkg.FindBin(jzb),totalBkg.GetNbinsX(),errTlb)
                                fil.write('%4f            %4.2f   %4.2f +/- %4.2f   %4.2f +/- %4.2f   %4.2f +/- %4.2f\n'%(jzb, SFpos, OFpos,errOF_pos,SFneg,errSF_neg,totlb,errTlb))
        fil.close()

    def dataCard(self):
        self.loadDataCard(jzbCut)
        for eta in ['Central','Forward']:
            for key in ['2jets','3jets']:
                for bjet in ['inc','0','1','2']:
                    SF_pos    = getattr(self,'datacard%s%s%s%s%s'%(key,eta,bjet,'SF','pos'))
                    SF_neg    = getattr(self,'datacard%s%s%s%s%s'%(key,eta,bjet,'SF','neg'))
                    OF_pos    = getattr(self,'datacard%s%s%s%s%s'%(key,eta,bjet,'OF','pos'))
                    OF_neg    = getattr(self,'datacard%s%s%s%s%s'%(key,eta,bjet,'OF','neg'))
                    OF_pos.Scale(self.ingMC.rsfof_final_cen_val if eta == 'Central' else self.ingMC.rsfof_final_fwd_val)                    
                    OF_neg.Scale(self.ingMC.rsfof_final_cen_val if eta == 'Central' else self.ingMC.rsfof_final_fwd_val)
                    SF_neg.Add(OF_neg,-1.)
                    for _bin in range(1,OF_pos.GetNbinsX()+1):
                        errStat = OF_pos.GetBinError(_bin)
                        errSyst = (self.ingMC.rsfof_final_cen_err if eta == 'Central' else self.ingMC.rsfof_final_fwd_err)*OF_pos.GetBinContent(_bin)
                        OF_pos.SetBinError(_bin, math.sqrt(errStat*errStat + errSyst*errSyst))
                        
                        errStat = SF_neg.GetBinError(_bin)
                        errSyst = 0.2*SF_neg.GetBinContent(_bin)
                        SF_neg.SetBinError(_bin,math.sqrt(errStat*errStat + errSyst*errSyst))

                        if _bin == SF_pos.GetNbinsX(): continue
                        mr   = 'lm'      if i == 1 else 'bz'     if i == 2 else 'oz'  if i == 3 else 'az'     if i == 4 else 'hm'
                        mass = 'lowMass' if i == 1 else 'belowZ' if i == 2 else 'onZ' if i == 3 else 'aboveZ' if i == 4 else 'highMass' if i == 5 else 'tooHigh'
                        obs  = SF_pos.GetBinContent(_bin)
                        fs   = OF_pos.GetBinContent(_bin)
                        fs_e = OF_pos.GetBinError(_bin)
                        tmp_rinout   = getattr(getattr(ingDA, 'rinout_dy_%s'%mr), '%s_val'%_eta)
                        tmp_rinout_e = getattr(getattr(ingDA, 'rinout_dy_%s'%mr), '%s_err'%_eta) if mr != 'oz' else 0.
                        dy   = SF_neg.GetBinContent(_bin)*tmp_rinout
                        dy_e = math.sqrt( (SF_neg.GetBinError(_bin)*tmp_rinout)**2 + (SF_neg.GetBinContent(_bin)*tmp_rinout_e)**2)
                        bin_name = '%s %s %s bjet'%(eta,key,bjet)
                        dc = '''# this is the datacard for bin {bin_name}
imax 1  number of channels
jmax 2  number of backgrounds
kmax *  number of nuisance parameters (sources of systematical uncertainties)
------------
# only one channel here with observation {obs}
bin            {bin_name}
observation    {obs}
------------
bin        {bin_name}     {bin_name}     {bin_name}
process    sig            FS             DY    
process    0              1              2
rate       XXRATEXX         {fs_bkg:.2f}       {dy_bkg:.2f}
------------
deltaS       lnN              1.15      -           -       
lumi         lnN              1.12      -           -       
FS_unc       lnN              -         {fs_unc:.2f}    -       
{bin_name}_fs_stat      gmN {of_yield}   -         {rsfof:.3f}        -       
DY_unc       lnN              -         -           {dy_unc:.2f}'''.format(bin_name=bin_name, obs=obs, fs_bkg=fs, fs_unc=1+rsfof_e, dy_bkg=dy, dy_unc=1+dy_e/dy, of_yield=int(of_yield), rsfof=rsfof)


        
    def RSFOF(self):
        cuts = CutManager.CutManager()
        fil = open('jzb_rsfof.txt','w')
        self.loadHistos()
        for key in ['2jets','3jets']:
            fil.write('%s \n'%key)
            for eta in ['Central','Forward']:
                fil.write('    %s\n'%eta)
                SF = self.treeTT.getTH1F(lumi,'tt%s%s%s'%('SF',eta,key),"t.lepsJZB_Edge - %f"%getattr(self,'center%s%s%s'%(key,'HT',eta)), 50, 0., 300, cuts.AddList([getattr(cuts,'GoodLepton%s'%'SF')(), 't.nJetSel_Edge == 2' if key == '2jets' else 't.nJetSel_Edge > 2', cuts.Zmass, '' if eta == '' else getattr(cuts,eta)() ]), "", "JZB [GeV]")
                OF = self.treeTT.getTH1F(lumi,'tt%s%s%s'%('OF',eta,key),"t.lepsJZB_Edge - %f"%getattr(self,'center%s%s%s'%(key,'HT',eta)), 50, 0., 300, cuts.AddList([getattr(cuts,'GoodLepton%s'%'OF')(), 't.nJetSel_Edge == 2' if key == '2jets' else 't.nJetSel_Edge > 2', cuts.Zmass, '' if eta == '' else getattr(cuts,eta)() ]), "", "JZB [GeV]")
                sf = r.Double(0.); sf_err = r.Double(0.)
                of = r.Double(0.); of_err = r.Double(0.) 
                for jzb in [100,150,200]:
                    sf = SF.IntegralAndError(SF.FindBin(jzb),SF.GetNbinsX(),sf_err)
                    of = OF.IntegralAndError(OF.FindBin(jzb),OF.GetNbinsX(),of_err)
                    fil.write('        R(SF/OF)  JZB cut: %4.0f  => %4.2f +/- %4.2f'%(jzb,sf/of,math.sqrt((sf_err/of)*(sf_err/of) + (sf*of_err/(of*of))*(sf*of_err/(of*of)))))
                    fil.write('\n')

if __name__ == "__main__":
    parser = optparse.OptionParser(usage="usage: %prog [opts] FilenameWithSamples", version="%prog 1.0")
    parser.add_option('-s', '--samples', action='store', type=str, dest='sampleFile', default='samples_skimmed.dat', help='the samples file. default \'samples_skimmed.dat\'')
    parser.add_option('-b', '--batch'  , action='store', type=int, dest='batchMode' , default=1            , help='set batch mode. default true')
    parser.add_option('-d', '--dody'   , action='store', type=int, dest='doDY'      , default=0            , help='do dy sample as well. default is false')
    parser.add_option('-g', '--dosig'   , action='store', type=int, dest='doSig'     , default=1            , help='do signal sample as well. default is true')
    parser.add_option('-q', '--quick'  , action='store', type=int, dest='quick'     , default=0            , help='quick test, not full thing. default false')
    parser.add_option('-i', '--ingredients', action='store', type=str, dest='ingredientFile', default='ingredients.dat', help='the ingredients file. default \'ingredients.dat\'')
    (opts, args) = parser.parse_args()


    if opts.batchMode:
        gROOT.SetBatch(True)
    gROOT.ProcessLine('.L include/tdrstyle.C')
    r.setTDRStyle() 

    print 'running with these options'
    for key, value in opts.__dict__.items():
        print '%-20s : %-20s' %(key, value)

    global lumi
    lumi = 2.3


    

    a = jzbAnalysis('a')
#    a.plotStack()
#    a.plotHistos()
#    a.plotClosure()
    a.loadTrees()
    cuts = CutManager.CutManager()
    histo = a.treeSI.getTH1F(lumi,'tt%s%s'%('OF',key),"t.lepsJZB_Edge", 50, 0., 300, cuts.AddList(['1']), "", "JZB [GeV]")
    c = TCanvas()
    histo.Draw()
    c.SaveAs('mierda.pdf')
#    a.makeEstimation(100.)
#    a.RSFOF()
