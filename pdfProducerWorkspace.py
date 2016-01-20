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
        self.zpt = ROOT.RooRealVar('lepsZPt_Edge'    , 'zpt',   0., 1000., 'GeV'); self.evt = ROOT.RooRealVar('evt'             , 'evt',   0.,1e7   , ''   );
        self.wgt = ROOT.RooRealVar('lumwgt'          , 'wgt', self.mcWgt, self.mcWgt, ''); self.wgt.setConstant(1);

        self.varSet = ROOT.RooArgSet(self.l1p, self.l2p, self.l1e, self.l2e, self.l1i, self.l2i, self.ldr, self.njs, self.nle)
        self.varSet.add(self.met); self.varSet.add(self.mlp); self.varSet.add(self.mll); self.varSet.add(self.mlp);
        self.varSet.add(self.mlb); self.varSet.add(self.st); self.varSet.add(self.zpt); self.varSet.add(self.evt);
        self.varSet.add(self.wgt);
        ## mlb_f = ROOT.RooFormulaVar('mlb','min_mlb1_Edge+min_mlb2_Edge', ROOT.RooArgList(ml1,ml2) )
        ## mlb   = dataSet.addColumn(mlb_f) ## this is now a RooRealVar

    def makeDatasets(self, cutlist):
        for cut, name in cutlist.items():
            print 'loading a dataset into the class for cut %s'%(name)
            cut = cut + ' && evt%{aa} == 0'.format(aa=self.skim)
            tmp_ds = ROOT.RooDataSet('ds_'+name, 'ds_'+name, self.tree, self.varSet, cut, 'lumwgt')
            setattr(self, name, tmp_ds)
            self.datasets.append(tmp_ds)

    def setVarProperties(self, var=0, _min=0, _max=0, rho=0):
        if not var:
            setattr(self, 'mlb_min',   0.); setattr(self, 'mlb_max', 750.); setattr(self, 'mlb_rho', 1.4);
            setattr(self, 'met_min', 150.); setattr(self, 'met_max', 750.); setattr(self, 'met_rho', 1.5);
            setattr(self, 'ldr_min',  0.3); setattr(self, 'ldr_max',  5.8); setattr(self, 'ldr_rho', 2.0);
            setattr(self, 'zpt_min',   0.); setattr(self, 'zpt_max', 600.); setattr(self, 'zpt_rho', 1.4);
            setattr(self, 'st_min' , 140.); setattr(self, 'st_max' ,1000.); setattr(self, 'st_rho' , 1.8);
            ##em_data.addNDPDF('mlp_pdf', 'cuts_of_sr_met150', 'am', 2.) ## have to mirror that one.
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


if __name__ == '__main__':

    print 'Starting r_SFOF analysis...'
    parser = optparse.OptionParser(usage='usage: %prog [opts] FilenameWithSamples', version='%prog 1.0')
    parser.add_option('-s', '--skim', action='store', type=int, dest='skimVal', default=1, help='take only every nth event into the samples.')
    (opts, args) = parser.parse_args()

    cuts = CutManager.CutManager()
    #skimCut = ' 1'
    #skimCut = " evt%{aa} == 0".format(aa = opts.skimVal)
    #print skimCut
    cuts_sr        = cuts.AddList([cuts.GoodLeptonNoTriggerOF(), cuts.METJetsSignalRegion      ]); cuts_sr        = cleanCut(cuts_sr)
    cuts_sr_nj2    = cuts.AddList([cuts.GoodLeptonNoTriggerOF(), cuts.METJetsSignalRegion2J    ]); cuts_sr_nj2    = cleanCut(cuts_sr_nj2)
    cuts_sr_nj3    = cuts.AddList([cuts.GoodLeptonNoTriggerOF(), cuts.METJetsSignalRegion3J    ]); cuts_sr_nj3    = cleanCut(cuts_sr_nj3)
    cuts_sr_met100 = cuts.AddList([cuts.GoodLeptonNoTriggerOF(), cuts.METJetsSignalRegionMET100]); cuts_sr_met100 = cleanCut(cuts_sr_met100)
    cuts_sr_met150 = cuts.AddList([cuts.GoodLeptonNoTriggerOF(), cuts.METJetsSignalRegionMET150]); cuts_sr_met150 = cleanCut(cuts_sr_met150)
    
    if not 'loaded' in globals():
        global loaded
        loaded = False

    dott = False
    
    if not loaded:
        tt_tfile = ROOT.TFile('/afs/cern.ch/work/p/pablom/public/MC/TTLep_pow/treeProducerSusyEdge/tree.root')
        tt_ffile = ROOT.TFile('/afs/cern.ch/work/m/mdunser/public/edgeTrees/friendsForPablosTrees_withSRID/evVarFriend_TTLep_pow.root')
        tt_tree = tt_tfile.Get('tree')
        tt_tree.AddFriend('sf/t', tt_ffile)
        em_tfile = ROOT.TFile('/afs/cern.ch/work/m/mdunser/public/edgeTrees/MuEG_all/MuEG/treeProducerSusyEdge/tree.root')
        em_ffile = ROOT.TFile('/afs/cern.ch/work/m/mdunser/public/edgeTrees/MuEG_all/friends/evVarFriend_MuEG.root')
        em_tree = em_tfile.Get('tree')
        em_tree.AddFriend('sf/t', em_ffile)
        
        em_data = pdfClass('emu_data', em_tree)
        if dott: tt_mc   = pdfClass('tt_mc'   , tt_tree, opts.skimVal, 87314.8*2.1/5e6)
        
        cutsAndDs = {#cuts_sr       : 'cuts_of_sr'       ,
                     cuts_sr_met150: 'cuts_of_sr_met150'}
                     #cuts_sr_nj2: 'cuts_of_sr_nj2',
                     #cuts_sr_nj3: 'cuts_of_sr_nj3'}
        em_data.makeDatasets(cutsAndDs)
        if dott: tt_mc  .makeDatasets(cutsAndDs)
        loaded = True
    
    w = ROOT.RooWorkspace('w')
    getattr(w,'import')(getattr(em_data, 'cuts_of_sr_met150'),ROOT.RooFit.Rename('em_data'), ROOT.RooFit.RenameVariable('st_Edge','st'))#,ROOT.RooFit.RenameVariable('lepsDR_Edge','ldr'))#, ROOT.RooFit.RenameVariable('lepsZPt_Edge','zpt'))
    if dott: 
        getattr(w,'import')(getattr(tt_mc  , 'cuts_of_sr_met150'),ROOT.RooFit.Rename('tt_mc'  ),ROOT.RooFit.RenameVariable('lepsDR_Edge','ldr'), ROOT.RooFit.RenameVariable('st_Edge','st'), ROOT.RooFit.RenameVariable('lepsZPt_Edge','zpt') )

    for var in ['mlb', 'met', 'zpt', 'ldr']:#'zpt', 'met', 'mlb', 'st', 'ldr']:
        opt = 'a' if var != 'met' else 'am'
        em_data.addNDPDF(var, 'cuts_of_sr_met150', opt, getattr(em_data, var+'_rho') )
        em_data.makeFrame(var)
        #if dott: tt_mc  .addNDPDF(var, 'cuts_of_sr_met150', opt, getattr(tt_mc  , var+'_rho') )
        if dott: tt_mc  .makeFrame(var)

    ## def buildMETModel(self, name, var, dataset):
    ##     self.w.factory("c_exp[0.0335,0.0335,0.0335]")
    ##     self.w.factory("c_inv[450.,450.,450.]")
    ##     ## make a lin*exp and a 1/x^n function. sum them up
    ##     metexp = ROOT.RooGenericPdf('metexp'+name, 'metexp'+name, '({var})*exp(-1.*(c_exp*{var}))'.format(var=var), ROOT.RooArgList(w.var(var),w.var('c_exp')) )
    ##     metinv = ROOT.RooGenericPdf('metinv'+name, 'metinv'+name, '1./({var}^c_inv)'.format(var=var), ROOT.RooArgList(w.var(var),w.var('c_inv')) )
    ##     self.w.factory('SUM::metpdf('+'metexp'+name+',metinv'+name+')')


    # ## FIT TO MET
    # w.factory("c_exp[0.,-1.,1.0]")
    # w.factory("c_inv[-1.,-3.,3.]")

    # ## make a lin*exp and a 1/x^n function. sum them up
    # name = 'test'; var = 'met_pt';
    # #metexp = ROOT.RooGenericPdf('metexp'+name, 'metexp'+name, '{var}*exp(-1.*(c_exp*{var}))'.format(var=var), ROOT.RooArgList(w.var(var),w.var('c_exp')) )
    # #metinv = ROOT.RooGenericPdf('metinv'+name, 'metinv'+name, '{var}^(c_inv)'.format(var=var), ROOT.RooArgList(w.var(var),w.var('c_inv')) )
    # #getattr(w,'import')(metexp)
    # #getattr(w,'import')(metinv)
    # #w.factory('SUM::metpdf(NSIG[0,5000]*'+'metexp'+name+',metinv'+name+')')
    # w.factory('EXPR::metpdf({var}*RooFit.RooExponential({var},c_exp))'.format(var=var))
    # fr_data=w.pdf('metpdf').fitTo(w.data('em_data'),ROOT.RooFit.Verbose(1),ROOT.RooFit.PrintLevel(1),ROOT.RooFit.NumCPU(1,0),ROOT.RooFit.Save(1) )
    # w.pdf('metpdf').plotOn(em_data.frame_met,ROOT.RooFit.LineColor(2), ROOT.RooFit.LineStyle(1) )

    ## ## FIT TO SUMMLB ## works nicely!!
    ## var = 'sum_mlb_Edge';
    ## w.factory("peak[170,120,220]")
    ## w.factory("sigma[40.,1.,100.]")
    ## w.factory("alpha[-1.,-2.5,0.5]")
    ## w.factory("n[0.5,0.,2.5]")
    ## #mlbexp = ROOT.RooCBShape('mlbcb'+name, 'mlbcb'+name, w.var(var), w.var('peak'), w.var('sigma'), w.var('alpha'), w.var('n') )
    ## w.factory('CBShape::mlbcbtest({var},peak,sigma,alpha,n)'.format(var=var))
    ## fr_data=w.pdf('mlbcbtest').fitTo(w.data('em_data'),ROOT.RooFit.Verbose(1),ROOT.RooFit.PrintLevel(1),ROOT.RooFit.NumCPU(1,0),ROOT.RooFit.Save(1) )
    ## w.pdf('mlbcbtest').plotOn(em_data.frame_mlb,ROOT.RooFit.LineColor(2), ROOT.RooFit.LineStyle(1) )

    ## FIT TO LDR
    #var = 'ldr';
    var = 'lepsDR_Edge';
    ## w.factory("start[0.3,0.,1.]")
    ## w.factory("slope[10.,5.,15.]")
    ## w.factory("end[3.,2.5,4.]")
    ## w.factory("gausmean[4.,0.,8.]")
    ## w.factory("gaussig[1.0,0.1,10.]")
    ## #w.factory('RooEdge::ldrtest({var},start,slope,end)'.format(var=var))
    ## w.factory('FCONV::ldrtest( {var}, RooEdge::ldredge({var},start,slope,end), Gaussian::ldrgaus({var},gausmean,gaussig) )'.format(var=var))
    ## fr_data=w.pdf('ldrtest').fitTo(w.data('em_data'),ROOT.RooFit.Verbose(1),ROOT.RooFit.PrintLevel(1),ROOT.RooFit.NumCPU(1,0),ROOT.RooFit.Save(1) )
    ## w.pdf('ldrtest').plotOn(em_data.frame_ldr,ROOT.RooFit.LineColor(2), ROOT.RooFit.LineStyle(1) )
    w.factory("peak[2.8,2.0,3.5]")
    w.factory("sigma[3.,0.,10.]")
    w.factory("alpha[1.,0.,4.5]")
    w.factory("n[9.5,-2.,30.0]")
    w.factory("gausmean1[1.,0.,3.]")
    w.factory("gaussig1[1.,0.,4.]")
    w.factory("gausmean2[3.,2.,5.]")
    w.factory("gaussig2[3.,0.,6.]")
    #w.factory('CBShape::ldrtest({var},peak,sigma,alpha,n)'.format(var=var))
    #w.factory('CBShape::ldrtest({var},peak,sigma,alpha,n)'.format(var=var))
    #w.factory('SUM::ldrtest(NSIG[0,10000]*Gaussian({var},gausmean1,gaussig1),NBKG[0,10000]*Gaussian({var},gausmean2,gaussig2))'.format(var=var))
    w.factory('SUM::ldrtest(NSIG[0,10000]*Uniform({var}),NBKG[0,10000]*CBShape({var},peak,sigma,alpha,n))'.format(var=var))
    ## doesn't seem to work w.factory("c[-1.,-2.,3.0]")
    ## doesn't seem to work w.factory("offset[3.,0,5.]")
    ## doesn't seem to work w.factory("width[1,0.,10.]")
    ## doesn't seem to work #w.factory('RooErfExpPdf::ldrtest({var},c,offset,width)'.format(var=var))
    ## doesn't seem to work w.factory('SUM::ldrtest(NSIG[0,100000]*RooErfExpPdf({var},c,offset,width),RooUniform({var}))'.format(var=var))
    fr_data=w.pdf('ldrtest').fitTo(w.data('em_data'),ROOT.RooFit.Verbose(1),ROOT.RooFit.PrintLevel(1),ROOT.RooFit.NumCPU(1,0),ROOT.RooFit.Save(1) )
    w.pdf('ldrtest').plotOn(em_data.frame_ldr,ROOT.RooFit.LineColor(2), ROOT.RooFit.LineStyle(1) )

    ## ## FIT TO Z-PT ## works nicely!!
    ## var = 'lepsZPt_Edge';
    ## w.factory("peak[60.,40.,100.]")
    ## w.factory("sigma[25.,10.,80.]")
    ## w.factory("alpha[-1.,-2.5,0.]")
    ## w.factory("n[0.5,0.,1e4]")
    ## w.factory('CBShape::zptcbtest({var},peak,sigma,alpha,n)'.format(var=var))
    ## fr_data=w.pdf('zptcbtest').fitTo(w.data('em_data'),ROOT.RooFit.Verbose(1),ROOT.RooFit.PrintLevel(1),ROOT.RooFit.NumCPU(1,0),ROOT.RooFit.Save(1) )
    ## w.pdf('zptcbtest').plotOn(em_data.frame_zpt,ROOT.RooFit.LineColor(2), ROOT.RooFit.LineStyle(1) )


