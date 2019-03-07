import ROOT, math, optparse, copy

##########33
## problemas e historias
##
#1) el breikindancing
#2) electrones: no hay mva tight vs miniiso 0.4 (solo loose)
#3) muones:     no hay ip3d < 8, solo ip3d < 4


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

    # Electrons
    elFile = ROOT.TFile('data/2017/leptonSF/ElectronScaleFactors_Run2017.root')
    ROOT.AddElec(copy.deepcopy( elFile.Get('Run2017_MVATightTightIP2D3D')))
    ROOT.AddElec(copy.deepcopy( elFile.Get('Run2017_MVAVLooseTightIP2DMini'))) 
    ROOT.AddElec(copy.deepcopy( elFile.Get('Run2017_ConvIHit0')))
    elFile.Close()

    elFile = ROOT.TFile('data/2017/leptonSF/egammaEffi.txt_EGM2D_runBCDEF_passingRECO.root')
    ROOT.AddElec(copy.deepcopy( elFile.Get('EGamma_SF2D')))
    elFile.Close()

    # Muons
    muFile = ROOT.TFile('data/2017/leptonSF/RunBCDEF_SF_ID.root')
    ROOT.AddMuon( copy.deepcopy( muFile.Get('NUM_MediumPromptID_DEN_genTracks_pt_abseta') ))
    muFile.Close()
    muFile = ROOT.TFile('data/2017/leptonSF/SF_num_miniiso_denmediumprompt.root')
    ROOT.AddMuon( copy.deepcopy( muFile.Get('TnP_MC_NUM_MiniIso02Cut_DEN_MediumCutidPromptCut_PAR_pt_eta') ))
    muFile.Close()


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

