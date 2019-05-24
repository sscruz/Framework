import ROOT, math, optparse, copy


# Electrons:
# GsfElectronToMVATightIDEmuTightIP2DSIP3D4
# MVAVLooseElectronToMini
# MVATightElectronToConvVetoIHit0
# Reco scale factors

# We need the "Emu" version on the first one for the trigger emulation cuts.

# Muons:
# Medium ID 
# miniIso0.2 vs Medium ID 
# dxy0.05, dz0.1 vs Medium ID 
# SIP3D<4 vs Medium ID 
# (reco scale factors when available)

def LoadLeptonSF():

    ROOT.gROOT.LoadMacro('include/LeptonSF.C+')

    ##### Electrons
    # 2016
    elFile = ROOT.TFile('data/2016/ElectronScaleFactors_Run2016.root')
    ROOT.AddElec(copy.deepcopy( elFile.Get('Run2016_LeptonMvaVTIDEmuTightIP2DSIP3D8miniIso04'))   ,2016)
    ROOT.AddElec(copy.deepcopy( elFile.Get('Run2016_ConvIHit0'))             ,2016)
    elFile.Close()
    elFile = ROOT.TFile('data/2016/EGM2D_BtoH_GT20GeV_RecoSF_Legacy2016.root')
    ROOT.AddElec(copy.deepcopy( elFile.Get('EGamma_SF2D')),2016)
    elFile.Close()


    # 2017 
    elFile = ROOT.TFile('data/2017/leptonSF/ElectronScaleFactors_Run2017.root')
    ROOT.AddElec(copy.deepcopy( elFile.Get('Run2017_MVATightTightIP2D3D'))   ,2017)
    ROOT.AddElec(copy.deepcopy( elFile.Get('Run2017_MVAVLooseTightIP2DMini')),2017) 
    ROOT.AddElec(copy.deepcopy( elFile.Get('Run2017_ConvIHit0'))             ,2017)
    elFile.Close()
    elFile = ROOT.TFile('data/2017/leptonSF/egammaEffi.txt_EGM2D_runBCDEF_passingRECO.root')
    ROOT.AddElec(copy.deepcopy( elFile.Get('EGamma_SF2D')),2017)
    elFile.Close()

    # 2018
    elFile = ROOT.TFile('data/2018/egammaEffi.txt_EGM2D_updatedAll.root')
    ROOT.AddElec(copy.deepcopy( elFile.Get('EGamma_SF2D')),2018)
    elFile.Close()

    elFile = ROOT.TFile('data/2018/ElectronScaleFactors_Run2018.root')
    ROOT.AddElec(copy.deepcopy( elFile.Get('Run2018_MVATightTightIP2D3D'))   ,2018)
    ROOT.AddElec(copy.deepcopy( elFile.Get('Run2018_Mini4')),2018) 
    ROOT.AddElec(copy.deepcopy( elFile.Get('Run2018_ConvIHit0'))             ,2018)
    elFile.Close()


    #### Muons

    # 2016
    muFile = ROOT.TFile('data/2016/TnP_NUM_MediumID_DENOM_generalTracks_VAR_map_pt_eta.root')
    ROOT.AddMuon( copy.deepcopy( muFile.Get('SF')),2016)
    muFile = ROOT.TFile('data/2016/TnP_NUM_MiniIsoLoose_DENOM_MediumID_VAR_map_pt_eta.root')
    ROOT.AddMuon( copy.deepcopy( muFile.Get('SF')),2016)
    muFile = ROOT.TFile('data/2016/TnP_NUM_TightIP2D_DENOM_MediumID_VAR_map_pt_eta.root')
    ROOT.AddMuon( copy.deepcopy( muFile.Get('SF')),2016)
    muFile = ROOT.TFile('data/2016/TnP_NUM_TightIP3D_DENOM_MediumID_VAR_map_pt_eta.root')
    ROOT.AddMuon( copy.deepcopy( muFile.Get('SF')),2016)

    # 2017
    muFile = ROOT.TFile('data/2017/leptonSF/RunBCDEF_SF_ID.root')
    ROOT.AddMuon( copy.deepcopy( muFile.Get('NUM_MediumPromptID_DEN_genTracks_pt_abseta') ),2017)
    muFile.Close()
    muFile = ROOT.TFile('data/2017/leptonSF/SF_num_miniiso_denmediumprompt.root')
    ROOT.AddMuon( copy.deepcopy( muFile.Get('TnP_MC_NUM_MiniIso02Cut_DEN_MediumCutidPromptCut_PAR_pt_eta') ),2017)
    muFile.Close()

    # 2018
    muFile = ROOT.TFile('data/2018/RunABCD_SF_ID.root')
    ROOT.AddMuon( copy.deepcopy( muFile.Get('NUM_MediumPromptID_DEN_TrackerMuons_pt_abseta') ),2018)
    muFile.Close()
    # iso missing....

LoadLeptonSF()

if __name__ == "__main__":
    print 'El', ROOT.LepSF(35,-2.3,11 ,'ElUp')
    print 'El', ROOT.LepSF(35, 2.3,11 ,'ElUp')
    print 'El', ROOT.LepSF(500, 2.3,11,'ElUp')
    print 'Mu', ROOT.LepSF(35,-2.3,13 ,'ElUp')
    print 'Mu', ROOT.LepSF(35,-2.3,13 ,'ElDn')
    print 'Mu', ROOT.LepSF(35,-2.3,13 ,'')
    print 'Mu', ROOT.LepSF(35,-2.3,11 ,'MuUp')
    print 'Mu', ROOT.LepSF(35,-2.3,11 ,'MuDn')
    print 'Mu', ROOT.LepSF(35,-2.3,11 ,'')

