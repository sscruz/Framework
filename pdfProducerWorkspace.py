import ROOT, math, optparse, copy

import include.CutManager as CutManager
import include.Sample     as Sample

import itertools



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
        self.tree.SetBranchStatus('lepsZPt_Edge'    , 1); self.tree.SetBranchStatus('d3D_Edge'        , 1);
        
    def makeVarSet(self):
        self.l1p = ROOT.RooRealVar('Lep1_pt_Edge'    , 'l1p',   0., 1000., 'GeV'); self.l2p = ROOT.RooRealVar('Lep2_pt_Edge'    , 'l2p',   0., 1000., 'GeV');
        self.l1e = ROOT.RooRealVar('Lep1_eta_Edge'   , 'l1e', -2.5,   2.5, ''   ); self.l2e = ROOT.RooRealVar('Lep2_eta_Edge'   , 'l2e', -2.5,   2.5, ''   );
        self.l1i = ROOT.RooRealVar('Lep1_pdgId_Edge' , 'l1i',  -15,    15, ''   ); self.l2i = ROOT.RooRealVar('Lep2_pdgId_Edge' , 'l2i',  -15,    15, ''   );
        self.ldr = ROOT.RooRealVar('lepsDR_Edge'     , 'ldr',  0.3,  5.8 , ''   ); self.njs = ROOT.RooRealVar('nJetSel_Edge'    , 'njs',   0 ,   25 , ''   );
        self.nle = ROOT.RooRealVar('nPairLep_Edge'   , 'nle',   -2,    3 , ''   ); self.met = ROOT.RooRealVar('met_pt'          , 'met',   0., 1000., 'GeV');
        self.mlp = ROOT.RooRealVar('metl1DPhi_Edge'  , 'mlp',   0.,3.142 , ''   ); self.mll = ROOT.RooRealVar('lepsMll_Edge'    , 'mll',  20., 1000., 'GeV');
        self.mlb = ROOT.RooRealVar('sum_mlb_Edge'    , 'mlb',   0., 3000., 'GeV'); self.st  = ROOT.RooRealVar('st_Edge'         , 'st' , 100., 5000., 'GeV');
        self.zpt = ROOT.RooRealVar('lepsZPt_Edge'    , 'zpt',   0., 1000., 'GeV'); self.evt = ROOT.RooRealVar('evt'             , 'evt',   0.,1e7   , ''   );
        self.a3d = ROOT.RooRealVar('d3D_Edge'        , 'a3d', 0.2 , 3.14 , ''   );
        self.wgt = ROOT.RooRealVar('lumwgt'          , 'wgt', self.mcWgt, self.mcWgt, ''); self.wgt.setConstant(1);

        self.varSet = ROOT.RooArgSet(self.l1p, self.l2p, self.l1e, self.l2e, self.l1i, self.l2i, self.ldr, self.njs, self.nle)
        self.varSet.add(self.met); self.varSet.add(self.mlp); self.varSet.add(self.mll); self.varSet.add(self.mlp);
        self.varSet.add(self.mlb); self.varSet.add(self.st); self.varSet.add(self.zpt); self.varSet.add(self.evt);
        self.varSet.add(self.wgt); self.varSet.add(self.a3d);

    def makeDatasets(self, cutlist):
        for cut, name in cutlist.items():
            print 'loading a dataset into the class for cut %s'%(name)
            cut = cut + ' && evt%{aa} == 0'.format(aa=self.skim)
            tmp_ds = ROOT.RooDataSet('ds_'+name, 'ds_'+name, self.tree, self.varSet, cut, 'lumwgt')
            setattr(self, name, tmp_ds)
            self.datasets.append(tmp_ds)

    def setVarProperties(self, var=0, _min=0, _max=0, rho=0):
        if not var:
            setattr(self, 'mlb_min',   0.); setattr(self, 'mlb_max', 3000.); setattr(self, 'mlb_rho', 1.4);
            setattr(self, 'met_min', 150.); setattr(self, 'met_max', 1000.); setattr(self, 'met_rho', 1.5);
            setattr(self, 'ldr_min',  0.3); setattr(self, 'ldr_max',   5.8); setattr(self, 'ldr_rho', 2.0);
            setattr(self, 'zpt_min',   0.); setattr(self, 'zpt_max', 1000.); setattr(self, 'zpt_rho', 1.4);
            setattr(self, 'a3d_min',  0.2); setattr(self, 'a3d_max',  3.14); setattr(self, 'a3d_rho', 1.4);
            setattr(self, 'st_min' , 140.); setattr(self, 'st_max' , 5000.); setattr(self, 'st_rho' , 1.8);
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
    elif var == 'a3d': vartree = 'd3D_Edge';
    if xmax <= xmin:
        return w.var(vartree).frame()
    else:
        return w.var(vartree).frame(xmin,xmax,100)

def buildModels(w):
    for var,t in itertools.product(['a3d', 'met', 'mlb', 'ldr', 'zpt'],['_DA', '_MC']):
        if var == 'met':
            vartree = 'met_pt'
            w.factory('met_c0{t}[-0.028,-1,0.]'.format(t=t))
            w.factory('met_c1{t}[-0.01,-0.05,-0.005]'.format(t=t))
            w.var(vartree).setMin(150.); w.var(vartree).setMax(1000.);
            w.factory('SUM::{var}_analyticalPDF{t}(met_n1{t}[0,1e4]*Exponential::exp1{t}({vartree},met_c0{t}),met_n2{t}[0,1e4]*Exponential::exp2{t}({vartree},met_c1{t}))'.format(var=var, vartree=vartree, t=t))
        elif var == 'mlb':
            vartree = "sum_mlb_Edge";
            w.factory("mlb_peak{t}[170,120,220]".format(t=t))
            w.factory("mlb_sigma{t}[40.,1.,100.]".format(t=t))
            w.factory("mlb_alpha{t}[-1.,-2.5,-0.5]".format(t=t))
            w.factory("mlb_n{t}[1.5,0.,2.5]".format(t=t))
            w.var(vartree).setMin(0.); w.var(vartree).setMax(3000.);
            w.factory("CBShape::{var}_analyticalPDF{t}({vartree},mlb_peak{t},mlb_sigma{t},mlb_alpha{t},mlb_n{t})".format(var=var, vartree=vartree,t=t))
        elif var == 'ldr':
            vartree = 'lepsDR_Edge';
            w.factory("ldr_peak{t}[3.0,2.0,3.5]".format(t=t))
            w.factory("ldr_sigma{t}[3.,0.,10.]".format(t=t))
            w.factory("ldr_alpha{t}[1.,0.,4.5]".format(t=t))
            w.factory("ldr_n{t}[2.5,1.,15.0]".format(t=t))
            #w.factory('CBShape::{var}_analyticalPDF({vartree},ldr_peak,ldr_sigma,ldr_alpha,ldr_n)'.format(var=var, vartree=vartree))
            w.factory('SUM::{var}_analyticalPDF{t}(ldr_n1{t}[0,1]*CBShape::{var}_cb{t}({vartree},ldr_peak{t},ldr_sigma{t},ldr_alpha{t},ldr_n{t}),ldr_n2{t}[0,1]*Gaussian::{var}_gaus{t}({vartree},ldr_gm{t}[0,-3,3],ldr_gs{t}[0.5,0.,0.8]))'.format(var=var, vartree=vartree, t=t))
        elif var == 'zpt':
            vartree = 'lepsZPt_Edge';
            w.factory("zpt_peak{t}[60.,30.,100.]".format(t=t))
            w.factory("zpt_sigma{t}[25., 5.,80.]".format(t=t))
            w.factory("zpt_alpha{t}[-1.,-2.5,0.]".format(t=t))
            w.factory("zpt_n{t}[5.0,5.0,100.0]".format(t=t))
            #w.factory('CBShape::{var}_analyticalPDF({vartree},zpt_peak,zpt_sigma,zpt_alpha,zpt_n)'.format(var=var, vartree=vartree))
            w.factory('SUM::{var}_analyticalPDF{t}(zpt{t}[0,1]*CBShape::{var}_cb{t}({vartree},zpt_peak{t},zpt_sigma{t},zpt_alpha{t},zpt_n{t}))'.format(var=var, vartree=vartree, t=t))
        elif var == 'a3d':
            vartree = 'd3D_Edge'
            ## w.factory("a3d_offset[1.5,1.4,1.6]")
            ## w.factory("a3d_a0[0.,-10.,10.]")
            ## w.factory("a3d_a1[0.,-10.,10.]")
            ## w.factory("a3d_a2[-1.,-5.,0.]")
            w.factory('EXPR::{var}_analyticalPDF{t}("a3d_a2{t}*({vartree}-a3d_offset{t})*({vartree}-a3d_offset{t}) + a3d_a1{t}*({vartree}-a3d_offset{t}) + a3d_a0{t}",{{{vartree},a3d_offset{t}[1.5,1.2,1.8],a3d_a0{t}[3.,2.5,5.],a3d_a1{t}[0.,-10,10],a3d_a2{t}[-1.,-3.,-1.]}})'.format(var=var, vartree=vartree,t=t))
            ## double gaussian
            ## w.factory("a3d_mean[1.5,1.2,1.7]")
            ## w.factory("a3d_sigma[1.,0.,3.]")
            ## #w.factory('Gaussian::{var}_gauss({vartree},a3d_mean,a3d_sigma)'.format(var=var, vartree=vartree))
            ## w.factory("a3d_mean2[1.5,1.5,2.0]")
            ## w.factory("a3d_sigma2[1.,0.,3.]")
            ## w.factory('Gaussian::{var}_gauss({vartree},a3d_mean,a3d_sigma)'.format(var=var, vartree=vartree))
            ## w.factory('Gaussian::{var}_gauss2({vartree},a3d_mean2,a3d_sigma2)'.format(var=var, vartree=vartree))
            ## w.factory('SUM::{var}_analyticalPDF(n11[0,1000]*a3d_gauss2,n22[0,1000]*a3d_gauss)'.format(var=var, vartree=vartree))

            ##double sided crystal ball approach (not very good)
            ## w.factory("a3d_mean[1.5,1.0,5.0]")
            ## w.factory("a3d_sigma[1.2,1.0,2.0]")
            ## w.factory("a3d_alpha1[1.,0.,2.]")
            ## w.factory("a3d_alpha2[1.,0.,2.]")
            ## w.factory("a3d_n1[1.,0.,10.]")
            ## w.factory("a3d_n2[1.,0.,10.]")
            ## w.factory("DoubleCB::{var}_analyticalPDF({vartree},a3d_mean,a3d_sigma,a3d_alpha1,a3d_n1,a3d_alpha2,a3d_n2)".format(var=var, vartree=vartree))

if __name__ == '__main__':

    print 'Starting pdf producer...'
    parser = optparse.OptionParser(usage='usage: %prog [opts] FilenameWithSamples', version='%prog 1.0')
    parser.add_option('-s', '--skim' , action='store', type=int , dest='skimVal', default=1, help='take only every nth event into the samples.')
    parser.add_option('-t', '--dott' , action='store_true', dest='dott' , help='do ttbar MC as well (takes time). default is 0')
    parser.add_option('-b', '--batch', action='store_true', dest='batch', help='set batch mode to not open windows. default false.')
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
    #tt_ffile = ROOT.TFile('/afs/cern.ch/work/m/mdunser/public/edgeTrees/friendsForPablosTrees_withSRID/evVarFriend_TTLep_pow.root')
    tt_ffile = ROOT.TFile('/afs/cern.ch/work/s/sesanche/public/nTuples/mar02/evVarFriend_TTLep_pow.root')
    tt_tree = tt_tfile.Get('tree')
    tt_tree.AddFriend('sf/t', tt_ffile)
    em_tfile = ROOT.TFile('/afs/cern.ch/work/m/mdunser/public/edgeTrees/MuEG_all/MuEG/treeProducerSusyEdge/tree.root')
    #em_ffile = ROOT.TFile('/afs/cern.ch/work/m/mdunser/public/edgeTrees/MuEG_all/friends/evVarFriend_MuEG.root')
    em_ffile = ROOT.TFile('/afs/cern.ch/work/s/sesanche/public/nTuples/mar02/evVarFriend_MuEG.root')
    em_tree = em_tfile.Get('tree')
    em_tree.AddFriend('sf/t', em_ffile)
    
    dss = []
    em_data = pdfClass('em_data', em_tree)
    dss.append(em_data)
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

    doNDKeys = False
    
    frs = []; pdf_histos = []; frames = []
    writeobjs = []
    for var in ['a3d', 'met', 'mlb', 'ldr', 'zpt']:
        opt = 'a' if var != 'met' else 'am'
        tmp_frame = getFrame(var); 
        for i, ds in enumerate(reversed(dss)):
            model = var+'_analyticalPDF'
            if ds.name == 'tt_mc'  : model += '_MC'
            if ds.name == 'em_data': model += '_DA'
            print '====================================================================='
            print '====================================================================='
            print '====== AT DATASET %s ================================================'%(ds.name)
            print '====================================================================='
            print '====================================================================='
            color = ROOT.kBlack if i==0 else ROOT.kRed
            linestyle = ROOT.kSolid if i == 0 else ROOT.kDashed
            if doNDKeys: ds.addNDPDF(var, 'cuts_of_sr_met150', opt, getattr(ds, var+'_rho') )

            tmp_fr=w.pdf(model).fitTo(w.data(ds.name),ROOT.RooFit.Verbose(1),ROOT.RooFit.PrintLevel(1),ROOT.RooFit.NumCPU(1,0),ROOT.RooFit.Save(1),ROOT.RooFit.Strategy(2))
            #tmp_fr=w.pdf(model).fitTo(w.data(ds.name),ROOT.RooFit.Save(1))
            frs.append(copy.deepcopy(tmp_fr))

            ## tmp_histo = w.pdf(model).createHistogram(ds.name+'_fit_histo_'+var+'_ds_cuts_of_sr_met150', getattr(ds,var), ROOT.RooFit.Binning(1000) )
            ## #tmp_histo.SetName(tmp_histo.GetName()[:tmp_histo.GetName().rfind('_ds_cuts_of_sr_met150')+len('_ds_cuts_of_sr_met150')] )
            ## tmp_histo.SetName(tmp_histo.GetName()[:tmp_histo.GetName().rfind('__')] )
            ## tmp_histo.SetLineColor(color); tmp_histo.SetLineStyle(linestyle)
            ## ds.ana_pdfs.append(tmp_histo); pdf_histos.append(tmp_histo)

            ## w.data(ds.name).plotOn(tmp_frame,ROOT.RooFit.MarkerColor(color), ROOT.RooFit.MarkerSize(0.8))
            ## w.pdf(model).plotOn(tmp_frame,ROOT.RooFit.LineColor(color), ROOT.RooFit.LineStyle(ROOT.kSolid) )
            ## if doNDKeys: getattr(ds, 'pdf_%s_ds_cuts_of_sr_met150'%var).plotOn(tmp_frame,ROOT.RooFit.LineColor(color), ROOT.RooFit.LineStyle(ROOT.kDashed) )

            ## # damn root takes the objects away if i write them now
            ## writeobjs.append(copy.deepcopy(tmp_histo))
            writeobjs.append(copy.deepcopy(w.pdf(model)))
            if doNDKeys: writeobjs.append(copy.deepcopy(getattr(ds, 'pdf_histo_%s_ds_%s'%(var, 'cuts_of_sr_met150') ) ))
            #a = raw_input()

        frames.append(tmp_frame);
        tmp_frame.Draw()

    w.writeToFile("pdfs/pdfs_version22_savingTheWorkspace.root")

    outfile = ROOT.TFile('pdfs/pdfs_version21_extendedFit.root', 'RECREATE')
    outfile.cd()
    for i in writeobjs:
        i.Write()
    for i in frs:
        i.Write()
    outfile.Close()


