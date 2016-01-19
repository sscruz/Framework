import ROOT as r

## import include.helper     as helper
## import include.Region     as Region
## import include.Canvas     as Canvas
import include.CutManager as CutManager
## import include.Sample     as Sample
## import include.Tables     as Tables


def deactivateBranches(tree):
    tree.SetBranchStatus('*', 0)
    tree.SetBranchStatus('nPairLep_Edge'     , 1)
    tree.SetBranchStatus('Lep1_pt_Edge'      , 1)
    tree.SetBranchStatus('Lep2_pt_Edge'      , 1)
    tree.SetBranchStatus('Lep1_eta_Edge'     , 1)
    tree.SetBranchStatus('Lep2_eta_Edge'     , 1)
    tree.SetBranchStatus('lepsDR_Edge'       , 1)
    tree.SetBranchStatus('lepsDPhi_Edge'     , 1)
    tree.SetBranchStatus('metl1DPhi_Edge'    , 1)
    tree.SetBranchStatus('metl2DPhi_Edge'    , 1)
    tree.SetBranchStatus('min_mlb1_Edge'     , 1)
    tree.SetBranchStatus('min_mlb2_Edge'     , 1)
    tree.SetBranchStatus('lepsMll_Edge'      , 1)
    tree.SetBranchStatus('Lep1_pdgId_Edge'   , 1)
    tree.SetBranchStatus('Lep2_pdgId_Edge'   , 1)
    tree.SetBranchStatus('met_pt'            , 1)
    tree.SetBranchStatus('nJetSel_Edge'      , 1)
    tree.SetBranchStatus('nBJetLoose35_Edge' , 1)
    tree.SetBranchStatus('htJet35j_Edge'     , 1)
    return tree


def srEv(ev):
    if ev.nPairLep_Edge < 1 : return False
    if ev.Lep1_pt_Edge < 20.: return False
    if ev.Lep2_pt_Edge < 20.: return False
    if ev.lepsDR_Edge < 0.3 : return False
    if 1.4 < abs(ev.Lep1_eta_Edge) < 1.6: return False
    if 1.4 < abs(ev.Lep2_eta_Edge) < 1.6: return False
    if ev.lepsMll_Edge < 20:  return False
    if (ev.Lep1_pdgId_Edge * ev.Lep2_pdgId_Edge) not in [-121, -169]: return False
    if (ev.met_pt < 150 and ev.nJetSel_Edge < 2): return False
    if (ev.met_pt < 100 and ev.nJetSel_Edge < 3): return False
    return True

r.TMVA.Tools.Instance()
 
# note that it seems to be mandatory to have an
# output file, just passing None to TMVA::Factory(..)
# does not work. Make sure you don't overwrite an
# existing file.

xmass = 700
ymass = 600
outfilename  = 'discrimination/m_%s_%s.root' %( (str(xmass), str(ymass)) if xmass > 0 else ('allMasses','') )
geq = '==' if xmass != 0 else '>='

fout = r.TFile(outfilename,'RECREATE')
 
factory = r.TMVA.Factory('TMVAClassification', fout,
                         ':'.join([
                         '!V',
                         '!Silent',
                         'Color',
                         'DrawProgressBar',
                         'Transformations=I;D;P;G,D',
                         'AnalysisType=Classification'])  )

#factory.AddVariable('t.nPairLep_Edge'     , 'F')
#factory.AddVariable('t.Lep1_pt_Edge'      , 'F')
#factory.AddVariable('t.Lep2_pt_Edge'      , 'F')
#factory.AddVariable('Lep1_eta_Edge'       , 'F')
#factory.AddVariable('Lep2_eta_Edge'       , 'F')
factory.AddVariable('t.lepsDR_Edge'       , 'F')
#factory.AddVariable('t.lepsDPhi_Edge'     , 'F')
#factory.AddVariable('t.metl1DPhi_Edge'    , 'F')
#factory.AddVariable('t.metl2DPhi_Edge'    , 'F')
factory.AddVariable('t.sum_mlb_Edge'    , 'F')
#factory.AddVariable('t.Lep1_eta_Edge*Lep2_eta_Edge'    , 'F')
#factory.AddVariable('t.lepsMll_Edge'      , 'F')
#factory.AddVariable('Lep1_pdgId_Edge'     , 'F')
#factory.AddVariable('Lep2_pdgId_Edge'     , 'F')
#factory.AddVariable('met_pt+t.Lep1_pt_Edge+t.Lep2_pt_Edge', 'F')
factory.AddVariable('met_pt', 'F')
factory.AddVariable('t.st_Edge', 'F')
factory.AddVariable('t.lepsZPt_Edge', 'F')
#factory.AddVariable('t.nJetSel_Edge'      , 'F')
#factory.AddVariable('t.nBJetLoose35_Edge' , 'F')
#factory.AddVariable('t.htJet35j_Edge'     , 'F')
 
## get ttbar tree and friends etc
tt_tfile = r.TFile('/afs/cern.ch/work/p/pablom/public/MC/TTLep_pow/treeProducerSusyEdge/tree.root')
tt_ffile = r.TFile('/afs/cern.ch/work/m/mdunser/public/edgeTrees/friendsForPablosTrees_withSRID/evVarFriend_TTLep_pow.root')
tt_tree = tt_tfile.Get('tree')
tt_tree.AddFriend('sf/t', tt_ffile)
## get signal tree and friends etc
t6_tfile = r.TFile('/afs/cern.ch/work/m/mdunser/public/edgeTrees/T6bbslepton/SMS_T6bbllslepton_mSbottom-600To900_mLSP-200To800/treeProducerSusyEdge/tree.root')
t6_ffile = r.TFile('/afs/cern.ch/work/m/mdunser/public/edgeTrees/T6bbslepton/friends/evVarFriend_SMS_T6bbllslepton_mSbottom-600To900_mLSP-200To800.root')
t6_tree = t6_tfile.Get('tree')
t6_tree.AddFriend('sf/t', t6_ffile)

factory.AddSignalTree    (tt_tree)
factory.AddBackgroundTree(t6_tree)


## tt_tree = deactivateBranches(tt_tree)
## t6_tree = deactivateBranches(t6_tree)

# cuts defining the signal and background sample
cuts = CutManager.CutManager()
cuts_sr_sf = cuts.AddList([cuts.GoodLeptonSFNoTrigger(), cuts.METJetsSignalRegion])
cuts_sr_sf = cuts_sr_sf.replace(cuts.twoLeptons, 't.nPairLep_Edge > 0')
sigCut = r.TCut(cuts_sr_sf)
bg_cut_str = cuts_sr_sf + ' && GenSusyMScan1 %s %.0f && GenSusyMScan2 %s %.0f'%( geq, xmass, geq, ymass)
print bg_cut_str
bgCut  = r.TCut(bg_cut_str)

h_minmlb = r.TH1F('h_minmlb', 'h_minmlb', 50, 0., 1000)
h_drleps = r.TH1F('h_drleps', 'h_drleps', 50, 0.,   5.)
h_met    = r.TH1F('h_met'   , 'h_met'   , 50, 0., 750.)
h_zpt    = r.TH1F('h_zpt'   , 'h_zpt'   , 50, 0., 750.)

tt_tree.Draw('(t.sum_mlb_Edge)>>h_minmlb', cuts_sr_sf)
tt_tree.Draw('(t.lepsDR_Edge)>>h_drleps' , cuts_sr_sf)
tt_tree.Draw('(met_pt)>>h_met' , cuts_sr_sf)
tt_tree.Draw('(t.lepsZPt_Edge)>>h_zpt'   , cuts_sr_sf)
 
h_minmlb.Scale(1./h_minmlb.Integral())
h_drleps.Scale(1./h_drleps.Integral())
h_met   .Scale(1./h_met   .Integral())
h_zpt   .Scale(1./h_zpt   .Integral())
#h_ml1phi.Scale(1./h_ml1phi.Integral())
h_minmlb.GetYaxis().SetRangeUser(0.,0.15)
h_drleps.GetYaxis().SetRangeUser(0.,0.07)
h_met   .GetYaxis().SetRangeUser(0.,0.04)
h_zpt   .GetYaxis().SetRangeUser(0.,0.10)
#h_ml1phi.GetYaxis().SetRangeUser(0.,0.04)


## tt_lh = r.TH1F('tt_lh',  'tt_lh', 100, 0, 0.30)
## t6_lh = r.TH1F('t6_lh',  't6_lh', 100, 0, 0.30)
## 
## i = 0
## for ev in t6_tree:
##     i=i+1
##     if i > 100000: break
##     if not i%10000: print i
##     if not srEv(ev): continue
##     tmp_lh_minmlb = h_minmlb.GetBinContent(h_minmlb.FindBin( (ev.min_mlb1_Edge + ev.min_mlb2_Edge) ))
##     tmp_lh_drleps = h_drleps.GetBinContent(h_drleps.FindBin(  ev.lepsDR_Edge                       ))
##     tmp_lh_ml1phi = h_ml1phi.GetBinContent(h_ml1phi.FindBin(  ev.metl1DPhi_Edge                    ))
##     tmp_lh = tmp_lh_minmlb+tmp_lh_drleps+tmp_lh_ml1phi
##     t6_lh.Fill(tmp_lh)
## 
## i=0
## for ev in tt_tree:
##     i=i+1
##     if i > 500000: break
##     if not i%10000: print i
##     if not srEv(ev): continue
##     tmp_lh_minmlb = h_minmlb.GetBinContent(h_minmlb.FindBin( (ev.min_mlb1_Edge + ev.min_mlb2_Edge) ))
##     tmp_lh_drleps = h_drleps.GetBinContent(h_drleps.FindBin(  ev.lepsDR_Edge                       ))
##     tmp_lh_ml1phi = h_ml1phi.GetBinContent(h_ml1phi.FindBin(  ev.metl1DPhi_Edge                    ))
##     tmp_lh = tmp_lh_minmlb+tmp_lh_drleps+tmp_lh_ml1phi
##     tt_lh.Fill(tmp_lh)
## 
## tt_lh.Scale(1./tt_lh.Integral())
## t6_lh.Scale(1./t6_lh.Integral())
## 
## tt_lh.SetLineColor(r.kBlue+2)
## t6_lh.SetLineColor(r.kRed +2)
## tt_lh.SetLineWidth(2)
## t6_lh.SetLineWidth(2)
## 
## t6_lh.Draw()
## tt_lh.Draw('same')


#print sfasfd
 
factory.PrepareTrainingAndTestTree(sigCut,   # signal events
                                   bgCut,    # background events
                                   ':'.join([
                                   'nTrain_Signal=0',
                                   'nTrain_Background=0',
                                   'SplitMode=Random',
                                   'NormMode=NumEvents',
                                   '!V' ]))

### some sort of bdt
#bdt = factory.BookMethod(r.TMVA.Types.kBDT, 'BDT', ':'.join([ '!H', '!V', 'NTrees=850', 'MinNodeSize=0.05', 'MaxDepth=3', 'BoostType=AdaBoost', 'AdaBoostBeta=0.5', 'SeparationType=GiniIndex', 'nCuts=20', 'PruneMethod=NoPruning', ]))

# Fisher discriminant (same as LD)
fisher = factory.BookMethod( r.TMVA.Types.kFisher, "Fisher", ':'.join(['H','!V','Fisher:CreateMVAPdfs','PDFInterpolMVAPdf=Spline2','NbinsMVAPdf=50','NsmoothMVAPdf=10']) )

## likelihood
lh = factory.BookMethod(r.TMVA.Types.kLikelihood, 'LikelihoodD', ':'.join([ '!H', '!V', '!TRansformOutput', 'PDFInterpol=Spline2', 'NSmoothSig[0]=20', 'NSmooth=5', 'NAvEvtPerBin=50', 'VarTransform=Decorrelate' ]))
 
factory.TrainAllMethods()
factory.TestAllMethods()
factory.EvaluateAllMethods()

fout.Close()

r.TMVA.TMVAGui(outfilename)

# # Logon not automatically loaded through PyROOT (logon loads TMVA library) load also GUI
# gROOT.SetMacroPath( "./" )
# gROOT.Macro       ( "./TMVAlogon.C" )
# r.gROOT.LoadMacro   ( '/afs/cern.ch/sw/lcg/app/releases/ROOT/5.34.18/x86_64-slc6-gcc48-opt/root/tmva/test/TMVAGui.C' )

## open the GUI for the result macros    
#r.gROOT.ProcessLine( 'TMVAGui("%s")' % fout.GetName() )
## keep the ROOT thread running
#r.gApplication.Run()

