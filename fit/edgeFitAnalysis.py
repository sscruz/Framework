import ROOT
from CutManager import CutManager
from Sample import Sample

## nS  =  3000
## nFS = 3000
## nDY =  500
## nZ  = 3000
## seed = 37
## 
## ROOT.RooRandom.randomGenerator().SetSeed(seed)
## c1 = ROOT.TCanvas("c1","c1")
## 
## w = ROOT.RooWorkspace()
## 
## w.factory("x[0,300],mll,mll")
## w.var("x").setBins(150)
## w.factory("nS[%d]" %nS)
## w.factory("nFS[%d]"%nFS)
## w.factory("nDY[%d]"%nDY)
## w.factory("nZ[%d]" %nZ)
## 
## 
## 
## # -------------------------------------------------------
## # --- OF CONTRIBUTION: POLY+EXPONENTIAL -----------------
## # -------------------------------------------------------
## w.factory("EXPR::of('a*(x^b) * exp(-c*x-d*x^2)',x,a[1.01212e-05],b[3.72646],c[0.0358278],d[-0.00000506689])")
## 
## # -------------------------------------------------------
## # --- DY CONTRIBUTION EXP+(BW*CB) -----------------------
## # -------------------------------------------------------
## 
## w.factory("BreitWigner::dybw(x, 91.1876,1.0021)") # the real width according to the pdg is 0.0021
## w.factory("DoubleCB::dydscb(x,  0, 1.6, 1.16, 2.9, 2.5, 1.04)")
## 
## #convolute the two
## w.factory("FCONV::dybwcb(x, dybw, dydscb)")
## 
## #make the exponential
## w.factory("EXPR::dyexp('exp(-(g*x))',x,g[0.01])")
## 
## #sum the convolution and the exponential with random fractions
## w.factory("SUM::dyall(0.4*dybwcb, 0.6*dyexp)")
## 
## # -------------------------------------------------------
## # --- SIGNAL: EDGE CONVOLUTED WITH GAUSSIAN -------------
## # -------------------------------------------------------
## w.factory("RooEdge::sigtri(x, 50, 1)")
## w.factory("Gaussian::siggau(x,0,sigma[5])")
## 
## #convolute the two
## w.factory("FCONV:signal(x, sigtri, siggau)")
## 
## # -------------------------------------------------------
## # --- MERGE ALL PDFS INTO ONE ---------------------------
## # -------------------------------------------------------
## w.factory("SUM::model_s(nFS*of, nDY*dyall, nS*signal)")
## w.factory("SUM::model_b(nFS*of, nDY*dyall)")
## 
## obs = ROOT.RooArgSet(w.var("x"))
## data_s = w.pdf("model_s").generate(obs,ROOT.RooFit.Extended())
## data_b = w.pdf("model_b").generate(obs,ROOT.RooFit.Extended())
## 
## 
## 
## 
## frame = w.var("x").frame()
## data_s.plotOn(frame)
## w.pdf("model_s").plotOn(frame, ROOT.RooFit.LineColor(ROOT.kRed))
## w.pdf("model_s").plotOn(frame, ROOT.RooFit.Components("background"))
## frame.Draw()
## c1.Print("data_s.png")
## 
## frame = w.var("x").frame()
## data_b.plotOn(frame)
## w.pdf("model_b").plotOn(frame)
## frame.Draw()
## c1.Print("data_b.png")
## 
## frame = w.var("x").frame()
## w.pdf("signal").plotOn(frame)
## frame.Draw()
## c1.Print("data_signalPDF.png")
## 
## frame = w.var("x").frame()
## w.pdf("dyexp").plotOn(frame)
## frame.Draw()
## c1.Print("data_dyexp.png")
## 
## frame = w.var("x").frame()
## w.pdf("dyall").plotOn(frame)
## frame.Draw()
## c1.Print("data_dyall.png")
## 
## 
## bdata_b = ROOT.RooDataHist("data_obs", "", obs, data_b)
## bdata_s = ROOT.RooDataHist("data_sig", "", obs, data_s)


##mll = ROOT.RooRealVar("t.lepsMll_Edge", "mll", 0, 300., "GeV")

l1p = ROOT.RooRealVar("Lep1_pt_Edge"    , "l1p",   0.,13000., "GeV")
l2p = ROOT.RooRealVar("Lep2_pt_Edge"    , "l2p",   0.,13000., "GeV")
l1e = ROOT.RooRealVar("Lep1_eta_Edge"   , "l1e", -2.5,   2.5, ""   )
l2e = ROOT.RooRealVar("Lep2_eta_Edge"   , "l2e", -2.5,   2.5, ""   )
l1i = ROOT.RooRealVar("Lep1_pdgId_Edge" , "l1i",  -15,    15, ""   )
l2i = ROOT.RooRealVar("Lep2_pdgId_Edge" , "l2i",  -15,    15, ""   )
ldr = ROOT.RooRealVar("lepsDR_Edge"     , "ldr",   0., 10.0 , ""   )
njs = ROOT.RooRealVar("nJetSel_Edge"    , "njs",   0 ,   25 , ""   )
nle = ROOT.RooRealVar("nPairLep_Edge"   , "nle",   0 ,   25 , ""   )
met = ROOT.RooRealVar("met_pt"          , "met",   0.,13000., "GeV")

mll = ROOT.RooRealVar("lepsMll_Edge"    , "mll",  20.,   300., "GeV")

# get the weights
genW = ROOT.RooRealVar("genWeight" , "genW",   -10000., 10000., "")

ntupleVarSet = ROOT.RooArgSet(l1p, l2p, l1e, l2e, l1i, l2i, ldr, mll, njs)
ntupleVarSet.add(nle)
ntupleVarSet.add(met)
ntupleVarSet.add(genW)


## tt_file = ROOT.TFile("/afs/cern.ch/work/m/mdunser/public/edgeTrees/74x_edgeSamples/TTJets/treeProducerSusyEdge/tree.root")
## tt_tree = tt_file.Get("tree")
## tt_tree.AddFriend("sf/t", "/afs/cern.ch/work/m/mdunser/public/edgeTrees/74x_edgeSamples/scalarFriends/evVarFriend_TTJets.root")
## 
## 
## dy50_file = ROOT.TFile("/afs/cern.ch/work/m/mdunser/public/edgeTrees/74x_edgeSamples/DYJetsToLL_M50/treeProducerSusyEdge/tree.root")
## dy50_tree = dy50_file.Get("tree")
## dy50_tree.AddFriend("sf/t", "/afs/cern.ch/work/m/mdunser/public/edgeTrees/74x_edgeSamples/scalarFriends/evVarFriend_DYJetsToLL_M50.root")
## 
## dy10_file = ROOT.TFile("/afs/cern.ch/work/m/mdunser/public/edgeTrees/74x_edgeSamples/DYJetsToLL_M10/treeProducerSusyEdge/tree.root")
## dy10_tree = dy10_file.Get("tree")
## dy10_tree.AddFriend("sf/t", "/afs/cern.ch/work/m/mdunser/public/edgeTrees/74x_edgeSamples/scalarFriends/evVarFriend_DYJetsToLL_M10.root")

ttjets = Sample('TTJets', 
                "/afs/cern.ch/work/m/mdunser/public/edgeTrees/74x_edgeSamples/", 
                "scalarFriends/", 
                831760.,
                14000948.,
                0)

dy50   = Sample('DYJetsToLL_M50', 
                "/afs/cern.ch/work/m/mdunser/public/edgeTrees/74x_edgeSamples/", 
                "scalarFriends/", 
                6025200.,
                19310834.824,
                0)

dy10   = Sample('DYJetsToLL_M10', 
                "/afs/cern.ch/work/m/mdunser/public/edgeTrees/74x_edgeSamples/", 
                "scalarFriends/", 
                18610000.,
                22217467.000,
                0)

samples = [ttjets, dy50, dy10]

for sample in [ttjets]:
    if not sample.name == 'TTJets': continue ## for testing
    print 'at sample', sample.name
    sample.ttree.SetBranchStatus("*"               , 0)
    sample.ttree.SetBranchStatus("Lep1_pt_Edge"    , 1)
    sample.ttree.SetBranchStatus("Lep2_pt_Edge"    , 1)
    sample.ttree.SetBranchStatus("Lep1_eta_Edge"   , 1)
    sample.ttree.SetBranchStatus("Lep2_eta_Edge"   , 1)
    sample.ttree.SetBranchStatus("Lep1_pdgId_Edge" , 1)
    sample.ttree.SetBranchStatus("Lep2_pdgId_Edge" , 1)
    sample.ttree.SetBranchStatus("lepsDR_Edge"     , 1)
    sample.ttree.SetBranchStatus("lepsMll_Edge"    , 1)
    sample.ttree.SetBranchStatus("nJetSel_Edge"    , 1)
    sample.ttree.SetBranchStatus("nPairLep_Edge"   , 1)
    sample.ttree.SetBranchStatus("met_pt"          , 1)
    sample.ttree.SetBranchStatus("genWeight"       , 1)



    print 'getting the roodataset from the tree'
    
    sample.dataSet = ROOT.RooDataSet("ntuple", "ntuple", sample.ttree, ntupleVarSet)
    
    ##tt_w = ROOT.RooRealVar("tt_w", "tt_w", 831760./14000948.)
    sample.lumWeight = ROOT.RooRealVar("lum_weight", "lum_weight", sample.lumWeight)

    sample.wgt_f = ROOT.RooFormulaVar("genS","lum_weight*sign(genWeight)", ROOT.RooArgList(sample.lumWeight,genW) )
    sample.wgt = sample.dataSet.addColumn(sample.wgt_f) ## this is now a RooRealVar
    ## wgt.setMin(-10.)
    ## wgt.setMax(10.)
    ## wgt_frame = wgt.frame()
    ## wgt.plotOn(wgt_frame)
    ## wgt_frame.Draw()
    
    
    cuts = CutManager()
    cut_sigNoMOF_cen = (cuts.Add(cuts.SignalNoMassLeptonOF() , cuts.central)).replace("t.", "")
    cut_sigNoMSF_cen = (cuts.Add(cuts.SignalNoMassLeptonSF() , cuts.central)).replace("t.", "")
    
    
    
    # reduce the dataset with the modified cut
    sample.ds_sigNoMOF_cen = sample.dataSet.reduce(cut_sigNoMOF_cen)
    sample.ds_sigNoMSF_cen = sample.dataSet.reduce(cut_sigNoMSF_cen)
    
    
    # get the weighted datasets
    ds_sigNoMOF_cen_weighted = ROOT.RooDataSet("dsName", "dsTitle", sample.ds_sigNoMOF_cen.get(), ROOT.RooFit.Import(sample.ds_sigNoMOF_cen), ROOT.RooFit.WeightVar(sample.wgt) )
    ds_sigNoMSF_cen_weighted = ROOT.RooDataSet("dsName", "dsTitle", sample.ds_sigNoMSF_cen.get(), ROOT.RooFit.Import(sample.ds_sigNoMSF_cen), ROOT.RooFit.WeightVar(sample.wgt) )
    
    
    mllframe = mll.frame(100)
    ds_sigNoMOF_cen_weighted.plotOn(mllframe, ROOT.RooFit.MarkerColor(1))
    ds_sigNoMSF_cen_weighted.plotOn(mllframe, ROOT.RooFit.MarkerColor(2))
    mllframe.Draw()
    
    
    if sample.name == "TTJets":
        aa = ROOT.RooRealVar("aa"   , "aa",    1.01212e-05, -0.01, 0.001)
        bb = ROOT.RooRealVar("bb"   , "bb",        3.72646,    0.,   10.)
        cc = ROOT.RooRealVar("cc"   , "cc",      0.0358278,    0.,    1.)
        dd = ROOT.RooRealVar("dd"   , "dd", -0.00000506689, -0.01,  0.01)
        
        sample.pdf = ROOT.RooGenericPdf("ofPDF", "ofPDF", "aa*(lepsMll_Edge^bb) * exp(-cc*lepsMll_Edge-dd*lepsMll_Edge^2)", ROOT.RooArgList(mll, aa, bb, cc, dd) )
    
        sample.fitres = sample.pdf.fitTo(sample.ds_sigNoMOF_cen, ROOT.RooFit.Save())##, ROOT.RooFit.Extended() )
        print 'fit performed for %s, these are the parameters' %(sample.name)
        print '=========================='
        sample.fitres.Print()
        print '=========================='
    
    sample.pdf.plotOn(mllframe)
    mllframe.Draw()
    
    

##data_obs = ROOT.RooDataSet("x","an mll spectrum from ttjets", mll , ROOT.RooFit.Import('lepsMll_Edge',mll) )

## x = ROOT.RooRealVar(mll,'x')
## #w.factory("EXPR::of1('m*(x^n) * exp(-o*x-p*x^2)',x,m[0.5*1.01212e-05, 2*1.01212e-05],n[0.5*3.72646, 2.*3.72646],o[0.5*0.0358278,2.*0.0358278],p[0.5*-0.00000506689,2*-0.00000506689])")
## 
## #of1pdf = ROOT.RooGenericPdf("of1pdf","1+sin(0.5*lepsMll_Edge)+abs(exp(0.1*lepsMll_Edge)*cos(-1*lepsMll_Edge))", ROOT.RooArgList(mll) )
## of1pdf = ROOT.RooGenericPdf("of1pdf","1+sin(0.5*lepsMll_Edge)+abs(exp(0.1*lepsMll_Edge)*cos(-1*lepsMll_Edge))", reducedDataset )


## ## ------------ Make histograms ---------------------------
## allHistFile = ROOT.TFile("simple-shapes-TH1.root", "RECREATE")
## ## Signal model
## signal_nominal = w.pdf("signal").createHistogram("x")
## signal_nominal.SetName("signal")
## signal_nominal.Scale(nS/signal_nominal.Integral())
## w.var("sigma").setVal(1.6)
## signal_sigmaUp = w.pdf("signal").createHistogram("x");  
## signal_sigmaUp.SetName("signal_sigmaUp")
## signal_sigmaUp.Scale(nS/signal_sigmaUp.Integral())
## w.var("sigma").setVal(0.7)
## signal_sigmaDown = w.pdf("signal").createHistogram("x");  
## signal_sigmaDown.SetName("signal_sigmaDown")
## signal_sigmaDown.Scale(nS/signal_sigmaDown.Integral())
## w.var("sigma").setVal(1.0)
## c1.Clear()
## signal_sigmaDown.Draw("H")
## signal_sigmaDown.SetLineColor(ROOT.kBlue)
## signal_sigmaDown.SetLineWidth(2)
## signal_sigmaUp.Draw("H SAME")
## signal_sigmaUp.SetLineColor(ROOT.kRed)
## signal_sigmaUp.SetLineWidth(2)
## signal_nominal.Draw("H SAME")
## signal_nominal.SetLineColor(ROOT.kBlack)
## signal_nominal.SetLineWidth(3)
## c1.Print("signal_model_binned.png")
## 
## frame = w.var("x").frame()
## w.pdf("signal").plotOn(frame, ROOT.RooFit.LineColor(ROOT.kBlack), ROOT.RooFit.LineWidth(3))
## w.var("sigma").setVal(1.6)
## w.pdf("signal").plotOn(frame, ROOT.RooFit.LineColor(ROOT.kBlue), ROOT.RooFit.LineWidth(2))
## w.var("sigma").setVal(0.7)
## w.pdf("signal").plotOn(frame, ROOT.RooFit.LineColor(ROOT.kRed), ROOT.RooFit.LineWidth(2))
## frame.Draw()
## c1.Print("signal_model_unbinned.png")
## ## background model
## frame = w.var("x").frame()
## background_nominal = w.pdf("background").createHistogram("x")
## background_nominal.SetName("background")
## background_nominal.Scale(nB/background_nominal.Integral())
## w.var("alpha").setVal(-0.2)
## background_alphaUp = w.pdf("background").createHistogram("x");  
## background_alphaUp.SetName("background_alphaUp")
## background_alphaUp.Scale(nB*1.15/background_alphaUp.Integral())
## w.pdf("background").plotOn(frame, ROOT.RooFit.LineColor(ROOT.kRed), ROOT.RooFit.LineWidth(2), ROOT.RooFit.Normalization(1.15))
## w.var("alpha").setVal(-0.4)
## background_alphaDown = w.pdf("background").createHistogram("x")
## background_alphaDown.SetName("background_alphaDown"); 
## background_alphaDown.Scale(nB*0.90/background_alphaDown.Integral())
## w.pdf("background").plotOn(frame, ROOT.RooFit.LineColor(ROOT.kBlue), ROOT.RooFit.LineWidth(2), ROOT.RooFit.Normalization(0.90))
## w.var("alpha").setVal(-0.3)
## w.pdf("background").plotOn(frame, ROOT.RooFit.LineColor(ROOT.kBlack), ROOT.RooFit.LineWidth(3), ROOT.RooFit.Normalization(1.0))
## frame.Draw(); c1.Print("background_model_unbinned.png")
## background_alphaDown.Draw("H"); 
## background_alphaDown.SetLineColor(ROOT.kBlue); 
## background_alphaDown.SetLineWidth(2)
## background_alphaUp.Draw("H SAME"); 
## background_alphaUp.SetLineColor(ROOT.kRed); 
## background_alphaUp.SetLineWidth(2)
## background_nominal.Draw("H SAME"); 
## background_nominal.SetLineColor(ROOT.kBlack); 
## background_nominal.SetLineWidth(3)
## c1.Print("background_model_binned.png")
## ## data
## hdata_b = bdata_b.createHistogram("x")
## hdata_b.SetName("data_obs")
## hdata_s = bdata_s.createHistogram("x")
## hdata_s.SetName("data_sig")
## ## write to file
## allHistFile.WriteTObject(signal_nominal)
## allHistFile.WriteTObject(signal_sigmaUp)
## allHistFile.WriteTObject(signal_sigmaDown)
## allHistFile.WriteTObject(background_nominal)
## allHistFile.WriteTObject(background_alphaUp)
## allHistFile.WriteTObject(background_alphaDown)
## allHistFile.WriteTObject(hdata_b)
## allHistFile.WriteTObject(hdata_s)
## 
## ## ------------ Make RooFit histograms ----------------------------------
## wB = ROOT.RooWorkspace("w","w")
## hobs = ROOT.RooArgList(w.var("x"))
## getattr(wB,'import')(bdata_b) #wB.import(bdata_b)
## getattr(wB,'import')(bdata_s) #wB.import(bdata_s)
## hsignal               = ROOT.RooDataHist("hsignal","",hobs,signal_nominal)
## hsignal_sigmaUp       = ROOT.RooDataHist("hsignal_sigmaUp","",hobs,signal_sigmaUp)
## hsignal_sigmaDown     = ROOT.RooDataHist("hsignal_sigmaDown","",hobs,signal_sigmaDown)
## hbackground           = ROOT.RooDataHist("hbackground","",hobs,background_nominal)
## hbackground_alphaUp   = ROOT.RooDataHist("hbackground_alphaUp","",hobs,background_alphaUp)
## hbackground_alphaDown = ROOT.RooDataHist("hbackground_alphaDown","",hobs,background_alphaDown)
## getattr(wB,'import')((ROOT.RooHistPdf("signal","",obs,hsignal))                            ) #wB.import(ROOT.RooHistPdf("signal","",obs,hsignal))
## getattr(wB,'import')((ROOT.RooHistPdf("signal_sigmaUp","",obs,hsignal_sigmaUp))            ) #wB.import(ROOT.RooHistPdf("signal_sigmaUp","",obs,hsignal_sigmaUp))
## getattr(wB,'import')((ROOT.RooHistPdf("signal_sigmaDown","",obs,hsignal_sigmaDown))        ) #wB.import(ROOT.RooHistPdf("signal_sigmaDown","",obs,hsignal_sigmaDown))
## getattr(wB,'import')((ROOT.RooHistPdf("background","",obs,hbackground))                    ) #wB.import(ROOT.RooHistPdf("background","",obs,hbackground))
## getattr(wB,'import')((ROOT.RooHistPdf("background_alphaUp","",obs,hbackground_alphaUp))    ) #wB.import(ROOT.RooHistPdf("background_alphaUp","",obs,hbackground_alphaUp))
## getattr(wB,'import')((ROOT.RooHistPdf("background_alphaDown","",obs,hbackground_alphaDown))) #wB.import(ROOT.RooHistPdf("background_alphaDown","",obs,hbackground_alphaDown))
## wB.writeToFile("simple-shapes-RooDataHist.root")
## 
## wBP = ROOT.RooWorkspace("w","w")
## getattr(wBP,'import')(bdata_b) #wBP.import(bdata_b)
## getattr(wBP,'import')(bdata_s) #wBP.import(bdata_s)
## getattr(wBP,'import')(w.pdf("signal")) #wBP.import(w.pdf("signal"))
## getattr(wBP,'import')(w.pdf("background")) #wBP.import(w.pdf("background"))
## wBP.writeToFile("simple-shapes-BinnedParam.root")
## 
## wUP = ROOT.RooWorkspace("w","w")
## wUP.var("x[0,10]")
## getattr(wUP,'import')(data_b, ROOT.RooFit.Rename("data_obs")) #wUP.import(data_b, Rename("data_obs"))
## getattr(wUP,'import')(data_s, ROOT.RooFit.Rename("data_sig")) #wUP.import(data_s, Rename("data_sig"))
## getattr(wUP,'import')(w.pdf("signal")) #wUP.import(w.pdf("signal"))
## getattr(wUP,'import')(w.pdf("background")) #wUP.import(w.pdf("background"))
## wUP.writeToFile("simple-shapes-UnbinnedParam.root")

## ## now we make a version in which the alpha is function of a unit gaussian, 
## ## so that we can do normalization and parametric morphing together
## RooWorkspace *wUPN = new RooWorkspace("w","w")
## wUPN.var("x[0,10]")
## wUPN.import(*data_b, Rename("data_obs"))
## wUPN.import(*data_s, Rename("data_sig"))
## wUPN.import(*w.pdf("signal"))
## RooAbsPdf *bgpdf = (RooAbsPdf *) *w.pdf("background").clone("background_norm")
## wUPN.import(*bgpdf)
## wUPN.factory("sum::param_alpha(-0.3,prod(alphaNorm[0],0.1))"); ## param_alpha = -0.3 + 0.1*alphaNorm, so Gauss(-0.3,1)
## wUPN.factory("EDIT::background(background_norm, alpha=param_alpha)")
## wUPN.Print("V")
## wUPN.writeToFile("simple-shapes-UnbinnedParamNorm.root")
    
## 
## if __name__ == "__main()__"
##     
