import ROOT, math, optparse, copy


class vParams:
    def __init__(self,name):
        self.name = name
        self.setParams()
    def setParams(self):
        if   self.name == 'met':
            self.xmin = 150; self.xmax =1000; self.nbins = 50; self.ymax = 0.35; self.title = 'ME_{T} (GeV)'       ; self.vtree = 'met_Edge'        ; 
        elif self.name == 'zpt':
            self.xmin = 0.0; self.xmax = 1000; self.nbins = 50; self.ymax = 0.24; self.title = 'p_{T}^{ll} (GeV)'   ; self.vtree = 'lepsZPt_Edge'; 
        elif self.name == 'mlb':
            self.xmin = 0.0; self.xmax = 3000; self.nbins = 50; self.ymax = 0.30; self.title = '#Sigma m_{lb} (Gev)'; self.vtree = 'sum_mlb_Edge'; 
        elif self.name == 'ldp':
            self.xmin = 0.0; self.xmax = 3.15; self.nbins = 50; self.ymax = 0.14; self.title = '#Delta 3D'; self.vtree = 'lepsDPhi_Edge'; 


def loadNLLStuff():
    ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.WARNING) ;
    ROOT.gROOT.LoadMacro('include/nll.C+')
    f_pdfs = ROOT.TFile.Open("/afs/cern.ch/work/s/sesanche/public/stuffForMoriond/pdfs_forMoriond_ver3.root","READ");
    w = f_pdfs.Get('w')
    ROOT.SetW(w)
    for var in 'zpt mlb met ldp'.split():
        getattr(ROOT,'Set%spdf'%var)(w.pdf("%s_analyticalPDF_DA"%var))
        print w.var(vParams(var).vtree)
        getattr(ROOT,'Setobs%s'%var)(ROOT.RooArgSet(w.var(vParams(var).vtree)))

loadNLLStuff()

