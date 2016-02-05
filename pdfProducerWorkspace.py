import ROOT, math, optparse, copy

import include.CutManager as CutManager
import include.Sample     as Sample

def cleanCut(cut):
    cut = cut.replace(cuts.twoLeptons, 't.nPairLep_Edge')
    cut = cut.replace('t.', '')
    return cut

class pdfClass:

    def __init__(self, _name, _ttree, _skimVal=1, isMC=False):
        self.name = _name
        self.tree = _ttree
        self.setTreeBranches()
        self.pdfs     = []
        self.ana_pdfs = []
        self.datasets = []
        self.frames   = []
        self.mcWgt = isMC if isMC else 1
        self.makeVarSet()
        self.skim = _skimVal
        ## setting some variables properties for the fit and the frame
        self.setVarProperties()


    def setTreeBranches(self):
        self.tree.SetBranchStatus('*'               , 0); self.tree.SetBranchStatus('evt'             , 1);
        self.tree.SetBranchStatus('Lep1_pt_Edge'    , 1); self.tree.SetBranchStatus('Lep2_pt_Edge'    , 1);
        self.tree.SetBranchStatus('Lep1_eta_Edge'   , 1); self.tree.SetBranchStatus('Lep2_eta_Edge'   , 1);
        self.tree.SetBranchStatus('Lep1_pdgId_Edge' , 1); self.tree.SetBranchStatus('Lep2_pdgId_Edge' , 1);
        self.tree.SetBranchStatus('lepsDR_Edge'     , 1); self.tree.SetBranchStatus('lepsMll_Edge'    , 1);
        self.tree.SetBranchStatus('nJetSel_Edge'    , 1); self.tree.SetBranchStatus('nPairLep_Edge'   , 1);
        self.tree.SetBranchStatus('met_pt'          , 1); self.tree.SetBranchStatus('min_mlb1_Edge'   , 1); 
        self.tree.SetBranchStatus('sum_mlb_Edge'    , 1); self.tree.SetBranchStatus('st_Edge'         , 1);
        self.tree.SetBranchStatus('lepsZPt_Edge'    , 1);
        
    def makeVarSet(self):
        self.l1p = ROOT.RooRealVar('Lep1_pt_Edge'    , 'l1p',   0., 1000., 'GeV'); self.l2p = ROOT.RooRealVar('Lep2_pt_Edge'    , 'l2p',   0., 1000., 'GeV');
        self.l1e = ROOT.RooRealVar('Lep1_eta_Edge'   , 'l1e', -2.5,   2.5, ''   ); self.l2e = ROOT.RooRealVar('Lep2_eta_Edge'   , 'l2e', -2.5,   2.5, ''   );
        self.l1i = ROOT.RooRealVar('Lep1_pdgId_Edge' , 'l1i',  -15,    15, ''   ); self.l2i = ROOT.RooRealVar('Lep2_pdgId_Edge' , 'l2i',  -15,    15, ''   );
        self.ldr = ROOT.RooRealVar('lepsDR_Edge'     , 'ldr',  0.3,  5.8 , ''   ); self.njs = ROOT.RooRealVar('nJetSel_Edge'    , 'njs',   0 ,   25 , ''   );
        self.nle = ROOT.RooRealVar('nPairLep_Edge'   , 'nle',   -2,    3 , ''   ); self.met = ROOT.RooRealVar('met_pt'          , 'met',   0., 1000., 'GeV');
        self.mlp = ROOT.RooRealVar('metl1DPhi_Edge'  , 'mlp',   0.,3.142 , ''   ); self.mll = ROOT.RooRealVar('lepsMll_Edge'    , 'mll',  20., 1000., 'GeV');
        self.mlb = ROOT.RooRealVar('sum_mlb_Edge'    , 'mlb',   0., 1000., 'GeV'); self.st  = ROOT.RooRealVar('st_Edge'         , 'st' , 100., 2000., 'GeV');
        self.zpt = ROOT.RooRealVar('lepsZPt_Edge'    , 'zpt',   0.,  400., 'GeV'); self.evt = ROOT.RooRealVar('evt'             , 'evt',   0.,1e7   , ''   );
        self.wgt = ROOT.RooRealVar('lumwgt'          , 'wgt', self.mcWgt, self.mcWgt, ''); self.wgt.setConstant(1);

        self.varSet = ROOT.RooArgSet(self.l1p, self.l2p, self.l1e, self.l2e, self.l1i, self.l2i, self.ldr, self.njs, self.nle)
        self.varSet.add(self.met); self.varSet.add(self.mlp); self.varSet.add(self.mll); self.varSet.add(self.mlp);
        self.varSet.add(self.mlb); self.varSet.add(self.st); self.varSet.add(self.zpt); self.varSet.add(self.evt);
        self.varSet.add(self.wgt);

    def makeDatasets(self, cutlist):
        for cut, name in cutlist.items():
            print 'loading a dataset into the class for cut %s'%(name)
            cut = cut + ' && evt%{aa} == 0'.format(aa=self.skim)
            tmp_ds = ROOT.RooDataSet('ds_'+name, 'ds_'+name, self.tree, self.varSet, cut, 'lumwgt')
            setattr(self, name, tmp_ds)
            self.datasets.append(tmp_ds)

    def setVarProperties(self, var=0, _min=0, _max=0, rho=0):
        if not var:
            setattr(self, 'mlb_min',   0.); setattr(self, 'mlb_max', 600.); setattr(self, 'mlb_rho', 1.4);
            setattr(self, 'met_min', 150.); setattr(self, 'met_max', 500.); setattr(self, 'met_rho', 1.5);
            setattr(self, 'ldr_min',  0.3); setattr(self, 'ldr_max',  5.8); setattr(self, 'ldr_rho', 2.0);
            setattr(self, 'zpt_min',   0.); setattr(self, 'zpt_max', 400.); setattr(self, 'zpt_rho', 1.4);
            setattr(self, 'st_min' , 140.); setattr(self, 'st_max' ,1000.); setattr(self, 'st_rho' , 1.8);
            self.propSet = True
        else:
            setattr(self, var+'_min', _min); setattr(self, var+'_max', _max); setattr(self, var+'_rho', rho);

    def addNDPDF(self, variable, dataset, option, alpha):
        ds = getattr(self, dataset)
        pdf_name = 'pdf_%s_%s'%(variable, ds.GetName())
        pdf_histo_name = 'pdf_histo_%s_%s'%(variable, ds.GetName())
        print 'getting NDKeysPDF for variable %s in dataset %s of sample %s' %(getattr(self,variable).getTitle(), ds.GetName(), self.name)
        self.setRange(variable)
        ret_pdf = ROOT.RooNDKeysPdf(pdf_name, pdf_name, ROOT.RooArgList(getattr(self,variable)), ds, option, alpha)
        setattr(self, pdf_name, ret_pdf)
        pdf_histo = ret_pdf.createHistogram(pdf_histo_name, getattr(self,variable), ROOT.RooFit.Binning(1000))
        pdf_histo.SetName(self.name+'_'+pdf_histo_name)
        pdf_histo.Scale(1./pdf_histo.Integral())
        setattr(self, pdf_histo_name, pdf_histo)
        self.pdfs.append(ret_pdf)
        return ret_pdf

    def setRange(self, variable):
        getattr(self, variable).setMin(getattr(self, variable+'_min'))
        getattr(self, variable).setMax(getattr(self, variable+'_max'))
        #if   variable == 'ldr':
        #    getattr(self, variable).setMin(0.); getattr(self, variable).setMax(5.8)
        #elif variable == 'mlb' or variable == 'met':
        #    getattr(self, variable).setMin(0.); getattr(self, variable).setMax(750)
        #elif variable == 'zpt':
        #    getattr(self, variable).setMin(0.); getattr(self, variable).setMax(600)
        #elif variable == 'st':
        #    getattr(self, variable).setMin(100); getattr(self, variable).setMax(1000)
        

    def makeFrame(self, variable, _min=0, _max=0, _nbins=100):
        if _min != _max:
            tmp_frame = getattr(self, variable).frame(_min, _max, _nbins)
        else:
            tmp_frame = getattr(self, variable).frame(getattr(self, variable+'_min'), getattr(self, variable+'_max'), _nbins)
        for i,ds in enumerate(self.datasets):
            ds.plotOn(tmp_frame, ROOT.RooFit.MarkerColor(1+i), ROOT.RooFit.MarkerSize(0.8))
            if hasattr(self, 'pdf_%s_%s'%(variable, ds.GetName()) ):
                getattr(self, 'pdf_%s_%s'%(variable, ds.GetName()) ).plotOn(tmp_frame, ROOT.RooFit.LineColor(1+i)) 
        self.frames.append(tmp_frame)
        setattr(self, 'frame_%s'%(variable), tmp_frame)

    def saveInFile(self):
        for obj in dir(self):
            if not obj.startswith('pdf_histo'): continue
            getattr(self, obj).Write()

def getFrame(var, xmin=0, xmax=0):
    if   var == 'met': vartree = 'met_pt'
    elif var == 'mlb': vartree = "sum_mlb_Edge";
    elif var == 'ldr': vartree = 'lepsDR_Edge';
    elif var == 'zpt': vartree = 'lepsZPt_Edge';
    if xmax <= xmin:
        return w.var(vartree).frame()
    else:
        return w.var(vartree).frame(xmin,xmax,100)

def buildModels(w):
    for var in ['met', 'mlb', 'ldr', 'zpt']:
        if var == 'met':
            vartree = 'met_pt'
            w.factory('c_0[-0.028,-1,0.]')
            w.factory('c_1[-0.01,-0.05,-0.005]')
            w.var(vartree).setMin(150.); w.var(vartree).setMax(1000.);
            w.factory('SUM::{var}_analyticalPDF(met_n1[0,1e4]*Exponential::exp1({vartree},c_0),met_n2[0,1e4]*Exponential::exp2({vartree},c_1))'.format(var=var, vartree=vartree))
        elif var == 'mlb':
            vartree = "sum_mlb_Edge";
            w.factory("mlb_peak[170,120,220]")
            w.factory("mlb_sigma[40.,1.,100.]")
            w.factory("mlb_alpha[-1.,-2.5,0.5]")
            w.factory("mlb_n[0.5,0.,2.5]")
            w.factory("CBShape::{var}_analyticalPDF({vartree},mlb_peak,mlb_sigma,mlb_alpha,mlb_n)".format(var=var, vartree=vartree))
        elif var == 'ldr':
            vartree = 'lepsDR_Edge';
            w.factory("ldr_peak[3.0,2.0,3.5]")
            w.factory("ldr_sigma[3.,0.,10.]")
            w.factory("ldr_alpha[1.,0.,4.5]")
            w.factory("ldr_n[2.5,1.,15.0]")
            #w.factory('CBShape::{var}_analyticalPDF({vartree},ldr_peak,ldr_sigma,ldr_alpha,ldr_n)'.format(var=var, vartree=vartree))
            w.factory('SUM::{var}_analyticalPDF(ldr_n1[0,1]*CBShape::{var}_cb({vartree},ldr_peak,ldr_sigma,ldr_alpha,ldr_n),ldr_n2[0,1]*Gaussian::{var}_gaus({vartree},ldr_gm[0,-3,3],ldr_gs[0.5,0.,0.8]))'.format(var=var, vartree=vartree))
        elif var == 'zpt':
            vartree = 'lepsZPt_Edge';
            w.factory("zpt_peak[60.,30.,100.]")
            w.factory("zpt_sigma[25., 5.,80.]")
            w.factory("zpt_alpha[-1.,-2.5,0.]")
            w.factory("zpt_n[0.5,-5.0,100.0]")
            w.factory('CBShape::{var}_analyticalPDF({vartree},zpt_peak,zpt_sigma,zpt_alpha,zpt_n)'.format(var=var, vartree=vartree))

if __name__ == '__main__':

    print 'Starting r_SFOF analysis...'
    parser = optparse.OptionParser(usage='usage: %prog [opts] FilenameWithSamples', version='%prog 1.0')
    parser.add_option('-s', '--skim' , action='store', type=int , dest='skimVal', default=1, help='take only every nth event into the samples.')
    parser.add_option('-t', '--dott' , action='store_true', dest='dott', help='do ttbar MC as well (takes time). default is 0')
    parser.add_option('-b', '--batch', action='store', type=int , dest='batch'  , default=0, help='set batch mode to not open windows. default false.')
    (opts, args) = parser.parse_args()

    if opts.batch:
        ROOT.gROOT.SetBatch()

    print 'running with these options \n'
    for key, value in opts.__dict__.items():
        print '%-20s : %-20s' %(key, value)

    cuts = CutManager.CutManager()
    cuts_sr        = cuts.AddList([cuts.GoodLeptonNoTriggerOF(), cuts.METJetsSignalRegion      ]); cuts_sr        = cleanCut(cuts_sr)
    cuts_sr_met100 = cuts.AddList([cuts.GoodLeptonNoTriggerOF(), cuts.METJetsSignalRegionMET100]); cuts_sr_met100 = cleanCut(cuts_sr_met100)
    cuts_sr_met150 = cuts.AddList([cuts.GoodLeptonNoTriggerOF(), cuts.METJetsSignalRegionMET150]); cuts_sr_met150 = cleanCut(cuts_sr_met150)
    
    tt_tfile = ROOT.TFile('/afs/cern.ch/work/p/pablom/public/MC/TTLep_pow/treeProducerSusyEdge/tree.root')
    tt_ffile = ROOT.TFile('/afs/cern.ch/work/m/mdunser/public/edgeTrees/friendsForPablosTrees_withSRID/evVarFriend_TTLep_pow.root')
    tt_tree = tt_tfile.Get('tree')
    tt_tree.AddFriend('sf/t', tt_ffile)
    em_tfile = ROOT.TFile('/afs/cern.ch/work/m/mdunser/public/edgeTrees/MuEG_all/MuEG/treeProducerSusyEdge/tree.root')
    em_ffile = ROOT.TFile('/afs/cern.ch/work/m/mdunser/public/edgeTrees/MuEG_all/friends/evVarFriend_MuEG.root')
    em_tree = em_tfile.Get('tree')
    em_tree.AddFriend('sf/t', em_ffile)
    
    em_data = pdfClass('em_data', em_tree)
    dss = [em_data]
    if opts.dott: 
        tt_mc   = pdfClass('tt_mc'   , tt_tree, opts.skimVal, 87314.8*2.1/5e6)
        dss.append(tt_mc)
    
    cutsAndDs = {#cuts_sr       : 'cuts_of_sr'       ,
                 #cuts_sr_met100: 'cuts_of_sr_met100'}
                 cuts_sr_met150: 'cuts_of_sr_met150'}
    w = ROOT.RooWorkspace('w')
    for ds in dss:
        ds.makeDatasets(cutsAndDs)
        getattr(w,'import')(getattr(ds, 'cuts_of_sr_met150'),ROOT.RooFit.Rename(ds.name), ROOT.RooFit.RenameVariable('st_Edge','st'))
    
    buildModels(w)
    
    frs = []; pdf_histos = []; frames = []
    writeobjs = []
    for var in ['met']:#, 'mlb', 'ldr']:
        opt = 'a' if var != 'met' else 'am'
        tmp_frame = getFrame(var, 150, 400); model = var+'_analyticalPDF'
        for i, ds in enumerate(dss):
            color = ROOT.kBlack if i==0 else ROOT.kRed
            linestyle = ROOT.kSolid if i == 0 else ROOT.kDashed
            ds.addNDPDF(var, 'cuts_of_sr_met150', opt, getattr(ds, var+'_rho') )

            tmp_fr=w.pdf(model).fitTo(w.data(ds.name),ROOT.RooFit.Verbose(1),ROOT.RooFit.PrintLevel(1),ROOT.RooFit.NumCPU(1,0),ROOT.RooFit.Save(1),ROOT.RooFit.Strategy(2))
            frs.append(tmp_fr)

            tmp_histo = w.pdf(model).createHistogram(ds.name+'_fit_histo_'+var+'_ds_cuts_of_sr_met150', getattr(ds,var), ROOT.RooFit.Binning(1000) )
            #tmp_histo.SetName(tmp_histo.GetName()[:tmp_histo.GetName().rfind('_ds_cuts_of_sr_met150')+len('_ds_cuts_of_sr_met150')] )
            tmp_histo.SetName(tmp_histo.GetName()[:tmp_histo.GetName().rfind('__')] )
            tmp_histo.SetLineColor(color); tmp_histo.SetLineStyle(linestyle)
            ds.ana_pdfs.append(tmp_histo); pdf_histos.append(tmp_histo)

            w.data(ds.name).plotOn(tmp_frame,ROOT.RooFit.MarkerColor(color), ROOT.RooFit.MarkerSize(0.8))
            w.pdf(model).plotOn(tmp_frame,ROOT.RooFit.LineColor(color), ROOT.RooFit.LineStyle(ROOT.kSolid) )
            getattr(ds, 'pdf_%s_ds_cuts_of_sr_met150'%var).plotOn(tmp_frame,ROOT.RooFit.LineColor(color), ROOT.RooFit.LineStyle(ROOT.kDashed) )

            # damn root takes the objects away if i write them now
            writeobjs.append(copy.deepcopy(tmp_histo))
            writeobjs.append(copy.deepcopy(getattr(ds, 'pdf_histo_%s_ds_%s'%(var, 'cuts_of_sr_met150') ) ))

        frames.append(tmp_frame);
        tmp_frame.Draw()

    outfile = ROOT.TFile('pdfs/pdfs_version14.root', 'RECREATE')
    outfile.cd()
    for i in writeobjs:
        i.Write()
    outfile.Close()

