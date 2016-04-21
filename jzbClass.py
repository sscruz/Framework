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
    def __init__(self, name, data=False):
        self.name = name
        self.doData = data
        self.loadIngs()
        self.binning = range(0,260,10); self.binning.extend([300])
    def loadIngs(self):
        self.ingMC = helper.ingredients(opts.ingredientFile, 'MC'  )
        self.ingDA = helper.ingredients(opts.ingredientFile, 'DATA')

    def loadTrees(self):
        print 'loading trees'
        ttDatasets = ['TTLep_pow']
        dyDatasetsHT = ['DYJetsToLL_M50_HT100to200', 'DYJetsToLL_M50_HT200to400',
                        'DYJetsToLL_M50_HT400to600', 'DYJetsToLL_M50_HT600toInf']
        mcDatasetsHT = ['DYJetsToLL_M50_HT100to200', 'DYJetsToLL_M50_HT200to400',
                        'DYJetsToLL_M50_HT400to600', 'DYJetsToLL_M50_HT600toInf','TTLep_pow',
                        'WWTo2L2Nu', 'WZTo2L2Q', 'ZZTo2L2Q','TToLeptons_tch_amcatnlo']
        dyDatasets   = ['DYJetsToLL_M50']
        mcDatasets   = [ 'DYJetsToLL_M50','TTLep_pow', 'WWTo2L2Nu', 'WZTo2L2Q', 'ZZTo2L2Q','TToLeptons_tch_amcatnlo']
        
        siDatasets   = ['SMS_T5ZZ_mGluino1000To1250_mLSP100To1000', 'SMS_T5ZZ_mGluino1200To1350_mLSP100To1200',
                        'SMS_T5ZZ_mGluino1400To1550_mLSP100To1400', 'SMS_T5ZZ_mGluino600To700_mLSP100To500',
                        'SMS_T5ZZ_mGluino800To950_mLSP100To800']
        daDatasets = ['MuonEG_Run2015C_25ns-05Oct_v1_runs_246908_260628',
                      'DoubleEG_Run2015C_25ns-05Oct_v1_runs_246908_260628',
                      'DoubleEG_Run2015D-05Oct_v1_runs_246908_260628',
                      'MuonEG_Run2015D-05Oct_v2_runs_246908_260628',
                      'DoubleMuon_Run2015D-05Oct_v1_runs_246908_260628',
                      'DoubleMuon_Run2015D_v4_runs_246908_260628',
                      'DoubleEG_Run2015D_v4_runs_246908_260628',
                      'DoubleMuon_Run2015C_25ns-05Oct_v1_runs_246908_260628',
                      'MuonEG_Run2015D_v4_runs_246908_260628']
        self.treeTT   = Sample.Tree(helper.selectSamples(opts.sampleFile, ttDatasets,   'TT'), 'TT', 0, isScan = 0)
        self.treeDA   = Sample.Tree(helper.selectSamples(opts.sampleFile, daDatasets,   'DA'), 'DA', 1, isScan = 0)
        self.treeDY   = Sample.Tree(helper.selectSamples(opts.sampleFile, dyDatasets,   'DY'), 'DY', 0, isScan = 0)
        self.treeMC   = Sample.Tree(helper.selectSamples(opts.sampleFile, mcDatasets,   'MC'), 'MC', 0, isScan = 0)
        self.treeDYHT = Sample.Tree(helper.selectSamples(opts.sampleFile, dyDatasetsHT, 'DY'), 'DY', 0, isScan = 0)
        self.treeMCHT = Sample.Tree(helper.selectSamples(opts.sampleFile, mcDatasetsHT, 'MC'), 'MC', 0, isScan = 0)

        self.treeSI = Sample.Tree(helper.selectSamples(opts.sampleFile, siDatasets, 'SI'), 'SI', 0, isScan = 1)

    def loadStack(self):
        print 'loading stacks'
        self.loadTrees()
        cuts = CutManager.CutManager()
        for key in ['2jets','3jets']:
            for ht in ['','HT']:
            #for ht in ['']:
                for fl in ['OF','SF']:
                    for eta in ['','Central','Forward']:
                        setattr(self,'stackMC%s%s%s%s'%(key,ht,fl,eta),copy.deepcopy(getattr(self,'treeMC%s'%ht).getStack(lumi, "jzb_%s_%s_%s%s"%(key,ht,fl,eta), "t.lepsJZB_Edge", 50, -300., 300., cuts.AddList([getattr(cuts,'GoodLepton%s'%fl)(), 't.nJetSel_Edge == 2' if key == '2jets' else 't.nJetSel_Edge > 2', cuts.Zmass, '' if eta == '' else getattr(cuts,eta)()]), "", "JZB [GeV]")))


    def loadHistos(self):
        print 'loading histos'
        cuts = CutManager.CutManager()
        self.loadTrees()
        for key in ['2jets','3jets']:
#            for ht in ['','HT']:
            for ht in ['HT']:
                for fl in ['SF','OF']:
                    for eta in ['Central','Forward']:
                        setattr(self,'histoDY%s%s%s%s'%(key,ht,fl,eta),copy.deepcopy(getattr(self,'treeDY%s'%ht).getTH1F(lumi,'dy%s%s%s%s'%(key,ht,fl,eta),"t.lepsJZB_Edge", 50, -300., 300., cuts.AddList([getattr(cuts,'GoodLepton%s'%fl)(), 't.nJetSel_Edge == 2' if key == '2jets' else 't.nJetSel_Edge > 2', cuts.Zmass,'' if eta == '' else getattr(cuts,eta)()]), "", "JZB [GeV]")))
                        setattr(self,'center%s%s%s'%(key,ht,eta),self.getCentering(getattr(self,'histoDY%s%s%s%s'%(key,ht,'SF',eta))))
                        
                        setattr(self,'histoDY%s%s%s%s_center'%(key,ht,fl,eta),copy.deepcopy(getattr(self,'treeDY%s'%ht).getTH1F(lumi,'dy%s%s%s%s_center'%(key,ht,fl,eta),"t.lepsJZB_Edge - %f"%getattr(self,'center%s%s%s'%(key,ht,eta)), 50, -300., 300., cuts.AddList([getattr(cuts,'GoodLepton%s'%fl)(), 't.nJetSel_Edge == 2' if key == '2jets' else 't.nJetSel_Edge > 2', cuts.Zmass,'' if eta == '' else getattr(cuts,eta)()]), "", "JZB [GeV]")))

    def loadClosure(self):
        print 'loading histos'
        self.loadTrees()
        cuts = CutManager.CutManager()
        for key in ['2jets','3jets']:
            for ht in ['HT']:
                for fl in ['OF','SF']: 
                    for eta in ['Central','Forward']:
                        setattr(self,'histoDY%s%s%s%s_pos'%(key,ht,fl,eta),copy.deepcopy(getattr(self,'treeDY%s'%ht).getTH1F(lumi,'dy%s%s%s%s_pos'%(key,ht,fl,eta),"t.lepsJZB_Edge - %f"%(getattr(self,'center%s%s%s'%(key,ht,eta))), self.binning, 0, 0, cuts.AddList([getattr(cuts,'GoodLepton%s'%fl)(), 't.nJetSel_Edge == 2' if key == '2jets' else 't.nJetSel_Edge > 2', cuts.Zmass, "t.lepsJZB_Edge - %f > 0"%(getattr(self,'center%s%s%s'%(key,ht,eta))),'' if eta == '' else getattr(cuts,eta)() ]), "", "JZB [GeV]",ofBin = True)))
                        setattr(self,'histoDY%s%s%s%s_neg'%(key,ht,fl,eta),copy.deepcopy(getattr(self,'treeDY%s'%ht).getTH1F(lumi,'dy%s%s%s%s_neg'%(key,ht,fl,eta),"TMath::Abs(t.lepsJZB_Edge - %f)"%getattr(self,'center%s%s%s'%(key,ht,eta)), self.binning, 0, 0, cuts.AddList([getattr(cuts,'GoodLepton%s'%fl)(), 't.nJetSel_Edge == 2' if key == '2jets' else 't.nJetSel_Edge > 2', cuts.Zmass, "t.lepsJZB_Edge - %f < 0"%getattr(self,'center%s%s%s'%(key,ht,eta)),'' if eta == '' else getattr(cuts,eta)() ]), "", "JZB [GeV]",ofBin = True))) 
                        

    def plotStack(self):
        self.loadStack()
        for key in ['2jets','3jets']:
            for ht in ['','HT']:
                for fl in ['OF','SF']:
                    for eta in ['','Central','Forward']:
                        plot = Canvas.Canvas('finalPlots/stack_%s_%s_%s_%s'%(key,ht,fl,eta),'png,pdf,cxx', 0.2, 0.60, 0.4, 0.8)  
                        plot.addStack(getattr(self,'stackMC%s%s%s%s'%(key,ht,fl,eta)),'HIST',True,0,['t#bar{t}','Drell-Yan','Rare'])
                        plot.save(True,False,True,lumi,0.,10*getattr(self,'stackMC%s%s%s%s'%(key,ht,fl,eta)).GetMaximum())
    def plotHistos(self):
        self.loadHistos()
        for key in ['2jets','3jets']:
            #for ht in ['','HT']:
            for ht in ['HT']:
                for eta in ['Central','Forward']:
                    plot = Canvas.Canvas('finalPlots/jzb_%s_%s_%s'%(key,ht,eta),'png,pdf,cxx', 0.4, 0.20, 0.6, 0.32)
                    SF   = getattr(self,'histoDY%s%s%s%s'%(key,ht,'SF',eta))
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
        self.loadHistos()
        self.loadClosure()
        for key in ['3jets','2jets']:
            for ht in ['HT']:
                for eta in ['Central','Forward']:
                    #plot = Canvas.Canvas('finalPlots/closure_%s_%s_%s'%(key,ht,eta),'png,pdf,cxx', 0.65, 0.50, 0.90, 0.72)
                    SF_pos = getattr(self,'histoDY%s%s%s%s_pos'%(key,ht,'SF',eta))
                    OF_pos = getattr(self,'histoDY%s%s%s%s_pos'%(key,ht,'OF',eta))
                    SF_pos.Add(OF_pos,-1.*self.ingMC.rsfof_final_cen_val if eta == 'Central' else -1.*self.ingMC.rsfof_final_fwd_val)
                    SF_neg = getattr(self,'histoDY%s%s%s%s_neg'%(key,ht,'SF',eta))
                    OF_neg = getattr(self,'histoDY%s%s%s%s_neg'%(key,ht,'OF',eta))
                    SF_neg.Add(OF_neg,-1.*self.ingMC.rsfof_final_cen_val if eta == 'Central' else -1.*self.ingMC.rsfof_final_fwd_val)
                    SF_pos.GetXaxis().SetRangeUser(0.,1000000.)
                    #plot.addHisto(SF_pos,"E,SAME", "JZB > 0","P",r.kRed,True,0)
                    #plot.addHisto(SF_neg,"HIST,SAME", "JZB < 0","L",r.kBlue,True,1)
                    #plot.addBandRatio( 0., 0.7, 400., 1.3, r.kRed, 0.4)
                    #plot.saveRatio(True,False,True,lumi,SF_pos,SF_neg)
                plot = Canvas.Canvas('finalPlots/closure_%s_%s'%(key,ht),'png,pdf,cxx', 0.65, 0.50, 0.90, 0.72) 
                SF_pos = getattr(self,'histoDY%s%s%s%s_pos'%(key,ht,'SF','Central')).Clone()
                SF_pos.Add(getattr(self,'histoDY%s%s%s%s_pos'%(key,ht,'SF','Forward')))
                SF_neg = getattr(self,'histoDY%s%s%s%s_neg'%(key,ht,'SF','Central')).Clone()
                SF_neg.Add(getattr(self,'histoDY%s%s%s%s_neg'%(key,ht,'SF','Forward')))
                plot.addHisto(SF_pos,"E,SAME", "JZB > 0","P",r.kRed,True,0)
                plot.addHisto(SF_neg,"HIST,SAME", "JZB < 0","L",r.kBlue,True,1)
                plot.addBandRatio( 0., 0.7, 350., 1.3, r.kRed, 0.4)
                fil = TFile("output.root","recreate")
                c = TCanvas()
                SF_pos.Write()
                SF_pos.Draw()
                print SF_pos.GetXaxis().GetXmax()
                c.SaveAs("hola2.pdf")
                fil.Close()
                print 'lets save'
                plot.saveRatio(True,False,True,lumi,SF_pos,SF_neg)

    def getCentering(self,histo):
        histo.Fit('gaus','0Q','',-30.,30.)
        return histo.GetFunction('gaus').GetParameter(1)
    
    def loadEstimation(self,jzbCut):
        self.loadHistos()
        cuts = CutManager.CutManager()
        for key in ['3jets']:
            for ht in ['HT']:
                for eta in ['Central','Forward']:
                    for fl in ['OF','SF']:
                        for sig in ['pos','neg']:
                            theCuts = cuts.AddList([getattr(cuts,'GoodLeptonNoTrigger%s'%fl)(), 't.nJetSel_Edge == 2' if key == '2jets' else 't.nJetSel_Edge > 2', cuts.Zmass,'' if eta == '' else getattr(cuts,eta)(),
                                                    '(lepsJZB_Edge -%f %s 0)'%(getattr(self,'center%s%s%s'%(key,ht,eta)),'>' if sig == 'pos' else '<')])
                            if self.doData: theCuts = theCuts.replace(cuts.twoLeptons, 'nPairLep_Edge > 0')
                            tree = getattr(self,'treeDA' if self.doData else 'treeMC%s'%ht )
                            setattr(self, 'jzb%s%s%s%s%s'%(key,ht,eta,fl,sig), tree.getTH1F(lumi,'jzb%s%s%s%s%s'%(key,ht,eta,fl,sig),"TMath::Abs(t.lepsJZB_Edge-%f)"%(getattr(self,'center%s%s%s'%(key,ht,eta))), 25, 0., 300, theCuts, "", "JZB [GeV]"))
                            
    def loadDataCard(self):
        self.loadHistos()
        cuts = CutManager.CutManager()
  
        for eta in ['central','forward']:
            for bjet in ['incb','0b','1b']:
                for fl in ['OF','SF']:
                    for sig in ['pos','neg']:
                        if 'incb' in bjet: bcut = ''
                        else: bcut = 'nBJetMedium35_Edge %s '%('== 0' if '0' in bjet  else '>= 1' )
                        etacut = 'Central' if eta == 'central' else 'Forward'
                        center = getattr(self,'center%s%s%s'%('3jets','HT',etacut))
                        setattr(self, 'datacard%s%s%s%s'%(eta,bjet,fl,sig), getattr(self,'treeMCHT').getTH1F(lumi,'datacard%s%s%s%s'%(eta,bjet,fl,sig),"t.lepsMll_Edge", [20., 70., 81., 101., 120., 13000.],0,0,cuts.AddList([getattr(cuts,'GoodLepton%s'%fl)(), cuts.Zmass, cuts.JZBJetsSignalRegion(center) if sig == 'pos' else cuts.JZBJetsControlRegion(center),'' if eta == '' else getattr(cuts,etacut)() , bcut]), "", "M_{ll} [GeV]"))
                                                
    def makeEstimation(self,jzbCut):
        self.loadEstimation(jzbCut)
        print 'Getting yields'
        for key in ['3jets']:
            for ht in ['HT']:
                 for eta in ['Central','Forward']:
                     for obs in ['jzb']:
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
                         totalBkg = OF_pos.Clone()
                         totalBkg.Add(SF_neg)
                         plot = Canvas.Canvas('finalPlots/signalregions/estimation_%s_%s_%s_%s'%(key,ht,eta,obs),'png,pdf,cxx', 0.65, 0.50, 0.90, 0.72)
                         plot.addStack(bkg_stack,'HIST',True,0,['FS background','Drell-Yan'])
                         plot.addHisto(SF_pos,"E,SAME", "Expected","P",r.kBlack,True,0)
                         plot.saveRatio(True,False,True,lumi,SF_pos,totalBkg)
                         #if obs == 'jzb':
                         #if False:
                          #   print 'writing', '%s %s %s \n' %(key,ht,eta)
                             #fil.write('%s %s %s \n' %(key,ht,eta))
                          #   for jzb in [100,150,200]:
                          #       errSF_pos = r.Double(0.); errOF_pos = r.Double(0.); errSF_neg = r.Double(0.); errTlb = r.Double(0.)
                          #       SFpos = SF_pos.IntegralAndError(SF_pos.FindBin(jzb),SF_pos.GetNbinsX(),errSF_pos)
                          #       OFpos = OF_pos.IntegralAndError(OF_pos.FindBin(jzb),OF_pos.GetNbinsX(),errOF_pos)
                          #       SFneg = SF_neg.IntegralAndError(SF_neg.FindBin(jzb),SF_neg.GetNbinsX(),errSF_neg)
                          #       totlb = totalBkg.IntegralAndError(totalBkg.FindBin(jzb),totalBkg.GetNbinsX(),errTlb)
                          #       print 'writing', '%4f            %4.2f   %4.2f +/- %4.2f   %4.2f +/- %4.2f   %4.2f +/- %4.2f\n'%(jzb, SFpos, OFpos,errOF_pos,SFneg,errSF_neg,totlb,errTlb)
                                 #fil.write('%4f            %4.2f   %4.2f +/- %4.2f   %4.2f +/- %4.2f   %4.2f +/- %4.2f\n'%(jzb, SFpos, OFpos,errOF_pos,SFneg,errSF_neg,totlb,errTlb))
    #fil.close()

    def dataCard(self):
        self.loadDataCard()
        for eta in ['central','forward']:
            _eta = 'cen' if eta == 'central' else 'fwd'
            for bjet in ['incb','0b','1b']:
                rsfof   = getattr(self.ingDA, 'rsfof_final_%s_val'%(_eta))
                rsfof_e = getattr(self.ingDA, 'rsfof_final_%s_err'%(_eta))
                SF_pos    = getattr(self,'datacard%s%s%s%s'%(eta,bjet,'SF','pos'))
                SF_neg    = getattr(self,'datacard%s%s%s%s'%(eta,bjet,'SF','neg'))
                OF_pos    = getattr(self,'datacard%s%s%s%s'%(eta,bjet,'OF','pos'))
                OF_neg    = getattr(self,'datacard%s%s%s%s'%(eta,bjet,'OF','neg'))
                of_yields = []
                for i in range(1,OF_pos.GetNbinsX()+1):
                    of_yields.append(OF_pos.GetBinContent(i)) # before scaling
                OF_pos.Scale(self.ingMC.rsfof_final_cen_val if eta == 'central' else self.ingMC.rsfof_final_fwd_val)  
                OF_neg.Scale(self.ingMC.rsfof_final_cen_val if eta == 'central' else self.ingMC.rsfof_final_fwd_val)
                SF_neg.Add(OF_neg,-1.)
                for i in range(1,OF_pos.GetNbinsX()+1):
                    if not i == 3: continue
                    errStat = OF_pos.GetBinError(i)
                    errSyst = (self.ingMC.rsfof_final_cen_err if eta == 'central' else self.ingMC.rsfof_final_fwd_err)*OF_pos.GetBinContent(i)
                    OF_pos.SetBinError(i, math.sqrt(errStat*errStat + errSyst*errSyst))
                    
                    errStat = SF_neg.GetBinError(i)
                    errSyst = 0.2*SF_neg.GetBinContent(i)
                    SF_neg.SetBinError(i,math.sqrt(errStat*errStat + errSyst*errSyst))
                    
                    if i == SF_pos.GetNbinsX(): continue
                    mr   = 'lm'      if i == 1 else 'bz'     if i == 2 else 'oz'  if i == 3 else 'az'     if i == 4 else 'hm'
                    mass = 'lowMass' if i == 1 else 'belowZ' if i == 2 else 'onZ' if i == 3 else 'aboveZ' if i == 4 else 'highMass' if i == 5 else 'tooHigh'
                    fil     = open('datacards/datacards_%s/%s_%s_%s.txt'%('T5ZZ_JZB',eta,mass,bjet),'w')

                    obs  = SF_pos.GetBinContent(i)
                    fs   = OF_pos.GetBinContent(i)
                    fs_e = OF_pos.GetBinError(i)
                    dy   = max(SF_neg.GetBinContent(i),0.)
                    dy_e = math.sqrt( SF_neg.GetBinError(i)**2 + (SF_neg.GetBinContent(i)*0.3)**2)
                    bin_name = '%s%sbjet'%(eta,bjet)
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
DY_unc       lnN              -         -           {dy_unc:.2f}'''.format(bin_name=bin_name, obs=int(obs), fs_bkg=fs, fs_unc=1+rsfof_e, dy_bkg=dy, dy_unc=1+dy_e/dy, of_yield=int(of_yields[i-1]), rsfof=rsfof)
                    fil.write(dc)
                    fil.close()
                    


        
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
    def fastPlots(self):
        self.loadTrees()
        cuts = CutManager.CutManager()
        self.xmasses = [800, 800,900, 900, 1000, 1000]
        self.ymasses = [100, 250,100, 250, 100,  250]
        color        = [r.kBlue,r.kRed,r.kMagenta, r.kBlack, r.kGreen, r.kOrange]
        self.mpoints = zip(self.xmasses,self.ymasses)
        canvas1 = Canvas.Canvas('fastPlots/met', 'png,pdf', 0.65, 0.40, 0.85, 0.72)
        canvas2 = Canvas.Canvas('fastPlots/njets', 'png,pdf', 0.65, 0.40, 0.85, 0.72)
        canvases_pt    = []
        canvases_phi    = []
        names       = ['gluino','lsp','chi2']
        ID          = [1000021,100022,1000023]
        for name in names:
            canv_pt  = Canvas.Canvas('fastPlots/pt_%s'%name,  'png,pdf', 0.65, 0.40, 0.85, 0.72)
            canv_phi = Canvas.Canvas('fastPlots/phi_%s'%name, 'png,pdf', 0.65, 0.40, 0.85, 0.72)
            canvases_pt.append(canv_pt)
            canvases_phi.append(canv_phi)

        for i,mp in enumerate(self.mpoints):
            cut = cuts.GoodLeptonNoTriggerSF()
            cut = cuts.AddList([cut.replace(cuts.twoLeptons, 't.nPairLep_Edge > 0'), 'GenSusyMScan1_Edge == %.0f && GenSusyMScan2_Edge == %.0f'%(mp[0], mp[1])])

            h = self.treeSI.getTH1F(1., 'si_%s_sr_sf_m%dv%d' %('met',mp[0],mp[1]), 'met_Edge', 50, 0., 2000. , cut, '', 'MET [GeV]')
            g = self.treeSI.getTH1F(1., 'si_%s_sr_sf_m%dv%d' %('njet',mp[0],mp[1]), 'nJetSel_Edge', 10, -0.5, 9.5 , cut, '', 'n_{jet}')
            h.SetLineColor(color[i])
            g.SetLineColor(color[i])
            canvas1.addHisto(h,'same,hist','m_{#tilde{g}} = %d, m_{#chi_{1}^{0}} = %d'%(mp[0],mp[1]),'l',color[i],1,i)
            canvas2.addHisto(g,'same,hist','m_{#tilde{g}} = %d, m_{#chi_{1}^{0}} = %d'%(mp[0],mp[1]),'l',color[i],1,i)
            for j,name in enumerate(names):
                pt = self.treeSI.getTH1F(1., 'si_%s_sr_sf_m%dv%d%s' %('met',mp[0],mp[1],name), 'genPart_pt', 50, 0., 2000. , cuts.AddList([cut,'genPart_pdgId == %d'%ID[j]]), '', 'MET [GeV]')
                canvases_pt[j].addHisto(pt,'same,hist','m_{#tilde{g}} = %d, m_{#chi_{1}^{0}} = %d'%(mp[0],mp[1]),'l',color[i],1,i)
        canvas1.save(1,False,0,2.1)
        canvas2.save(1,False,0,2.1)
        for canvas in canvases_pt:
            canvas.save(1,False,0,2.1)

if __name__ == "__main__":
    parser = optparse.OptionParser(usage="usage: %prog [opts] FilenameWithSamples", version="%prog 1.0")
    parser.add_option('-s', '--samples', action='store', type=str, dest='sampleFile', default='samples.dat', help='the samples file. default \'samples.dat\'')
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
#    a.fastPlots()
#    a.dataCard()
#    a.plotStack()
#    a.plotHistos()
    a.plotClosure()
#    a.makeEstimation(100.)
#    a.RSFOF()
