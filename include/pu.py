import ROOT, math, optparse, copy


def loadPUStuff():
    ROOT.gROOT.LoadMacro('include/pu.C+')
    puFile = ROOT.TFile("/afs/cern.ch/user/p/pablom/public/pileup_FULL_nominalUpDown.root","READ")
    puHist   = copy.deepcopy( puFile.Get('weightsNominal') )
#    puHistUp = copy.deepcopy( puFile.Get('weightsUp') )
#    puHistDn = copy.deepcopy( puFile.Get('weightsDown') )
    puFile.Close()

    puFileOld = ROOT.TFile("/afs/cern.ch/work/m/mdunser/public/puWeighting/2016/pileup_jul21_nominalUpDown.root","READ")
    puHistOld = copy.deepcopy( puFileOld.Get('weightsNominal') )
    puFileOld.Close()
    
    ROOT.SetPUhist(puHist)
    ROOT.SetOldPUHist(puHistOld)
#    SetPUhistUp(puHistUp)
#    SetPUhistDn(puHistDn)

loadPUStuff()

