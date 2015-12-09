import ROOT as r

## import include.helper     as helper
## import include.Region     as Region
## import include.Canvas     as Canvas
import include.CutManager as CutManager
## import include.Sample     as Sample
## import include.Tables     as Tables

r.TMVA.Tools.Instance()
 
# note that it seems to be mandatory to have an
# output file, just passing None to TMVA::Factory(..)
# does not work. Make sure you don't overwrite an
# existing file.

fout = r.TFile('discrimination/test.root','RECREATE')
 
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
factory.AddVariable('t.lepsDPhi_Edge'     , 'F')
factory.AddVariable('t.metl1DPhi_Edge'    , 'F')
factory.AddVariable('t.metl2DPhi_Edge'    , 'F')
factory.AddVariable('t.min_mlb1_Edge+t.min_mlb2_Edge'    , 'F')
#factory.AddVariable('t.lepsMll_Edge'      , 'F')
#factory.AddVariable('Lep1_pdgId_Edge'     , 'F')
#factory.AddVariable('Lep2_pdgId_Edge'     , 'F')
#factory.AddVariable('met_pt'              , 'F')
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

xmass = 900
ymass = 800
 
# cuts defining the signal and background sample
cuts = CutManager.CutManager()
cuts_sr_sf = cuts.AddList([cuts.GoodLeptonSFNoTrigger(), cuts.METJetsSignalRegion])
cuts_sr_sf = cuts_sr_sf.replace(cuts.twoLeptons, 't.nPairLep_Edge > 0')
sigCut = r.TCut(cuts_sr_sf)
bgCut  = r.TCut(cuts_sr_sf + ' && GenSusyMScan1 == %.0f && GenSusyMScan2 == %.0f'%(xmass, ymass))
 
factory.PrepareTrainingAndTestTree(sigCut,   # signal events
                                   bgCut,    # background events
                                   ':'.join([
                                        'nTrain_Signal=0',
                                        'nTrain_Background=0',
                                        'SplitMode=Random',
                                        'NormMode=NumEvents',
                                        '!V'
                                       ]))

bdt = factory.BookMethod(r.TMVA.Types.kBDT, 'BDT',
                         ':'.join([
                         '!H',                       ## don't print help
                         '!V',                       ## not verbose
                         'NTrees=850',               ## 850 trees used
                         'nEventsMin=150',           ## leaf nodes must contain 150 events (euh?)
                         'MaxDepth=3',               ## depth of tree limited to 3
                         'BoostType=AdaBoost',
                         'AdaBoostBeta=0.5',         ## adaptive boosting used
                         'SeparationType=GiniIndex', ## gini index used for best variable
                         'nCuts=20',                 ## 20 steps when scanning cuts on a variable
                         'PruneMethod=NoPruning',    ## no pruning of trees
                         ]))


lh = factory.BookMethod(r.TMVA.Types.kLikelihood, 'LikelihoodD',
                        ':'.join([
                        '!H',
                        '!V',
                        '!TRansformOutput',
                        'PDFInterpol=Spline2',
                        'NSmoothSig[0]=20',
                        'NSmooth=5',
                        'NAvEvtPerBin=50',
                        'VarTransform=Decorrelate'
                        ]))
 
factory.TrainAllMethods()
factory.TestAllMethods()
factory.EvaluateAllMethods()

fout.Close()
