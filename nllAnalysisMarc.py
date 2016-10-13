import ROOT as r
from   ROOT import gROOT, TCanvas, TFile, TF1, TH1F, TGraph
import math, sys, optparse, copy, re, array, itertools


import include.helper     as helper
import include.Region     as Region
import include.Canvas     as Canvas
import include.CutManager as CutManager
import include.Sample     as Sample
import include.Tables     as Tables


class vParams:
    def __init__(self,name):
        self.name = name
        self.setParams()
    def setParams(self):
        if   self.name == 'met':
            self.xmin = 150; self.xmax =1000; self.nbins = 50; self.ymax = 0.35; self.title = 'ME_{T} (GeV)'       ; self.vtree = 't.met_Edge'        ; 
        elif self.name == 'ldr':
            self.xmin = 0.3; self.xmax = 5.8; self.nbins = 50; self.ymax = 0.08; self.title = '#Delta R_{ll}'      ; self.vtree = 't.lepsDR_Edge' ; 
        elif self.name == 'zpt':
            self.xmin = 0.0; self.xmax = 1000; self.nbins = 50; self.ymax = 0.24; self.title = 'p_{T}^{ll} (GeV)'   ; self.vtree = 't.lepsZPt_Edge'; 
        elif self.name == 'mlb':
            self.xmin = 0.0; self.xmax = 3000; self.nbins = 50; self.ymax = 0.30; self.title = '#Sigma m_{lb} (Gev)'; self.vtree = 't.sum_mlb_Edge'; 
        elif self.name == 'a3d':
            self.xmin = 0.2; self.xmax = 3.14; self.nbins = 50; self.ymax = 0.14; self.title = '#Delta 3D'; self.vtree = 't.d3D_Edge'; 
        elif self.name == 'ldp':
            self.xmin = 0.; self.xmax = 3.14; self.nbins = 50; self.ymax = 0.14; self.title = '#Delta #phi'; self.vtree = 't.lepsDPhi_Edge'; 
        elif self.name == 'nll':
            self.xmin = 10; self.xmax =  30; self.nbins = 50; self.ymax = 0.20; self.title = 'NLL'              
            self.vtree = '-1.*TMath::Log(lh_ana_met_data_Edge*lh_ana_mlb_data_Edge*lh_ana_a3d_data_Edge*lh_ana_zpt_data_Edge)'; 
        elif self.name == 'nll_nomet':
            self.xmin = 15; self.xmax =  35; self.nbins = 50; self.ymax = 0.15; self.title = 'NLL-no MET' ; 
            self.vtree = '-1.*TMath::Log(lh_ana_mlb_data_Edge*lh_ana_ldp_data_Edge*lh_ana_zpt_data_Edge)'; 
        elif self.name == 'nll_nomlb':
            self.xmin = 15; self.xmax =  35; self.nbins = 50; self.ymax = 0.15; self.title = 'NLL-no #Sigma m_{lb}' ; 
            self.vtree = '-1.*TMath::Log(lh_ana_met_data_Edge*lh_ana_ldp_data_Edge*lh_ana_zpt_data_Edge)'; 
        elif self.name == 'nll_noldr':
            self.xmin = 15; self.xmax =  35; self.nbins = 50; self.ymax = 0.15; self.title = 'NLL-no #Delta #phi_{ll}' ; 
            self.vtree = '-1.*TMath::Log(lh_ana_met_data_Edge*lh_ana_mlb_data_Edge*lh_ana_zpt_data_Edge)'; 
        elif self.name == 'nll_nozpt':
            self.xmin = 15; self.xmax =  35; self.nbins = 50; self.ymax = 0.15; self.title = 'NLL-no p_{T}^{ll}'   ; 
            self.vtree = '-1.*TMath::Log(lh_ana_met_data_Edge*lh_ana_mlb_data_Edge*lh_ana_ldp_data_Edge)';
        elif self.name == 'nll_a3d':
            self.xmin = 15; self.xmax =  35; self.nbins = 50; self.ymax = 0.15; self.title = 'NLL'                ; 
            self.vtree = '-1.*TMath::Log(lh_ana_met_data_Edge*lh_ana_mlb_data_Edge*lh_ana_a3d_data_Edge*lh_ana_zpt_data_Edge)'; 
        elif self.name == 'nll_ldp':
            self.xmin = 15; self.xmax =  35; self.nbins = 50; self.ymax = 0.15; self.title = 'NLL'                ; 
            self.vtree = '-1.*TMath::Log(lh_ana_met_data_Edge*lh_ana_mlb_data_Edge*lh_ana_ldp_data_Edge*lh_ana_zpt_data_Edge)'; 
        elif self.name == 'mt2':
            self.xmin = 0.0; self.xmax = 300; self.nbins = 50; self.ymax = 0.15; self.title = 'MT_{2} (GeV)'       ; self.vtree = 't.mt2_Edge'    ; 
        elif self.name == 'p.d.f._m_et':
            self.xmin = -5.; self.xmax = 50.; self.nbins = 50; self.ymax = 0.15; self.title = 'f_{met}'; self.vtree = '-TMath::Log(lh_ana_met_data_Edge)'
        elif self.name == 'p.d.f._l_dr':
            self.xmin = 0.; self.xmax = 10.; self.nbins = 50; self.ymax = 0.15; self.title = 'f_{ldr}'; self.vtree = '-TMath::Log(lh_ana_ldr_data_Edge)'
        elif self.name == 'p.d.f._z_pt':
            self.xmin = 0.; self.xmax = 10.; self.nbins = 50; self.ymax = 0.15; self.title = 'f_{ldr}'; self.vtree = '-TMath::Log(lh_ana_zpt_data_Edge)'
        elif self.name == 'p.d.f._m_lb':
            self.xmin = 0.; self.xmax = 10.; self.nbins = 50; self.ymax = 0.15; self.title = 'f_{mlb}'; self.vtree = '-TMath::Log(lh_ana_mlb_data_Edge)'


class nllObject:
    def __init__(self, name, pdffname, doDY=False, doSig = True):
        print 'initializing nllObject'
        self.name = name
        self.pdffile   = r.TFile(pdffname)
        self.loadIngs()
#        self.kerfile   = r.TFile('/afs/cern.ch/work/m/mdunser/public/pdfsForLikelihood/pdfs_version20_fitResultsAndFunctions.root')
#        self.kerfile   = r.TFile('/afs/cern.ch/work/m/mdunser/public/pdfsForLikelihood/pdfs_version23_extendedFit_may2016.root')
        self.w         = self.pdffile.Get('w')
        self.doDY = doDY
        self.doSig = doSig
        self.ivars     = ['ldr', 'met', 'mlb', 'zpt']
#        self.typs      = ['fit', 'ker']
        self.typs      = ['fit']
        self.all_pdfs  = []; self.ana_pdfs  = []; self.ker_pdfs  = []
        self.all_cums  = []; self.ana_cums  = []; self.ker_cums  = []

        self.tt_pdfs  = []; self.tt_cums  = [];
        self.da_pdfs  = []; self.da_cums  = [];
        self.dy_pdfs  = []; self.dy_cums  = [];
        self.si_pdfs  = []; self.si_cums  = [];

        self.tt_pdfs_Zleq_90  = []; self.tt_cums_Zleq_90  = [];
        self.da_pdfs_Zleq_90  = []; self.da_cums_Zleq_90  = [];
        self.dy_pdfs_Zleq_90  = []; self.dy_cums_Zleq_90  = [];
        self.si_pdfs_Zleq_90  = []; self.si_cums_Zleq_90  = [];

        self.tt_pdfs_Zgeq_90  = []; self.tt_cums_Zgeq_90  = [];
        self.da_pdfs_Zgeq_90  = []; self.da_cums_Zgeq_90  = [];
        self.dy_pdfs_Zgeq_90  = []; self.dy_cums_Zgeq_90  = [];
        self.si_pdfs_Zgeq_90  = []; self.si_cums_Zgeq_90  = [];
        print 'easy stuff done'
        self.loadPDFs()
        print 'pdfs loadded'
        self.hasNLLPDFs = False; self.treesLoaded = False; self.treeHistsLoaded = False
        self.loadTrees()
        self.mpoints = []
        self.xmasses = [900, 800, 500, 400]
        self.ymasses = [200, 250, 400, 375]
        self.setMassPoints()

    def loadIngs(self):
        self.ingMC = helper.ingredients(opts.ingredientFile, 'MC'  )
        self.ingDA = helper.ingredients(opts.ingredientFile, 'DATA')

    def setMassPoints(self):
        if len(self.xmasses) is not len(self.ymasses):
            print 'you have more xmasses than ymasses or vice versa'
            print 'you are in trouble now'
        for point in zip(self.xmasses,self.ymasses):
            self.mpoints.append(point)
        
    def loadTrees(self):
        ttDatasets = ['TTLep_pow_ext']
        siDatasets = [#'SMS_T6bbllslepton_mSbottom600To775_mLSP150To725',
                      'SMS_T6bbllslepton_mSbottom400To575_mLSP150To550',
                      'SMS_T6bbllslepton_mSbottom800To950_mLSP150To900']
        mcDatasets = ['TTLep_pow_ext']#,'DYJetsToLL_M50']
        dyDatasets = ['DYJetsToLL_M50_LO','DYJetsToLL_M50_HT100to200','DYJetsToLL_M50_HT200to400','DYJetsToLL_M50_HT400to600',
                      'DYJetsToLL_M50_HT600toInf']
#        dyDatasets = ['DYJetsToLL_M50']

                   #   'SMS_T6bbllslepton_mSbottom-400To550_mLSP-200To500' ]
        daDatasets = [#'HTMHT_Run2015C_25ns-05Oct_v1_runs_246908_260628',
                      'MuonEG_Run2015C_25ns-05Oct_v1_runs_246908_260628',
                      'DoubleEG_Run2015C_25ns-05Oct_v1_runs_246908_260628',
                      #'JetHT_Run2015D-05Oct_v1_runs_246908_260628',
                      'DoubleEG_Run2015D-05Oct_v1_runs_246908_260628',
                      'MuonEG_Run2015D-05Oct_v2_runs_246908_260628',
                      #'HTMHT_Run2015D_v4_runs_246908_260628',
                      'DoubleMuon_Run2015D-05Oct_v1_runs_246908_260628',
                      'DoubleMuon_Run2015D_v4_runs_246908_260628',
                      #'HTMHT_Run2015D-05Oct_v1_runs_246908_260628',
                      #'JetHT_Run2015D_v4_runs_246908_260628',
                      'DoubleEG_Run2015D_v4_runs_246908_260628',
                      #'JetHT_Run2015C_25ns-05Oct_v1_runs_246908_260628',
                      'DoubleMuon_Run2015C_25ns-05Oct_v1_runs_246908_260628',
                      'MuonEG_Run2015D_v4_runs_246908_260628']
        self.treeMC = Sample.Tree(helper.selectSamples(opts.sampleFile, mcDatasets, 'MC'), 'MC', 0, isScan = 0)
        print 'reading from file', helper.selectSamples(opts.sampleFile, ttDatasets, 'TT')
        self.treeTT = Sample.Tree(helper.selectSamples(opts.sampleFile, ttDatasets, 'TT'), 'TT', 0, isScan = 0)
        self.treeDY = Sample.Tree(helper.selectSamples(opts.sampleFile, dyDatasets, 'DY'), 'DY', 0, isScan = 0)
        self.treeSI = Sample.Tree(helper.selectSamples(opts.sampleFile, siDatasets, 'SI'), 'SI', 0, isScan = 1)
        self.treeDA = Sample.Tree(helper.selectSamples(opts.sampleFile, daDatasets, 'DA'), 'DA', 1, isScan = 0)

        if self.doDY:
            dyDatasets = ['DYJetsToLL_M50']
            self.treeDY = Sample.Tree(helper.selectSamples(opts.sampleFile, dyDatasets, 'DY'), 'DY', 0, isScan = 0)
        self.treesLoaded = True

    def loadTreeHists(self, doMassCuts = False):
        if not self.treesLoaded: self.loadTrees()
        for massCut in ( [[], ['lepsMll_Edge < 90'],['lepsMll_Edge > 90'] ] if doMassCuts else [[]]) :
            #for var in self.ivars+['nll']+['ldp']:# + ['nll_noldr', 'nll_a3d']:
            #for var in ['nll','ldp']:
            for var in ['mt2']:

                print 'loading tree histograms for variable', var
                c_sr_sf    = cuts.AddList([cuts.METJetsSignalRegionMET150 , cuts.GoodLeptonSF()]+massCut)
                c_sr_of    = cuts.AddList([cuts.METJetsSignalRegionMET150 , cuts.GoodLeptonOF()]+massCut)
                v = vParams(var)
                xmax = v.xmax-(v.xmax-v.xmin)/v.nbins if var not in ['nll', 'ldr'] else v.xmax
                h_tt_sr_sf  = self.treeTT.getTH1F(1., 'tt_%s_sr_sf' %var, v.vtree, v.nbins, v.xmin, v.xmax , c_sr_sf   , '', v.title, ofBin= False); h_tt_sr_sf.Scale(1./h_tt_sr_sf.Integral())
                h_tt_sr_of  = self.treeTT.getTH1F(1., 'tt_%s_sr_of' %var, v.vtree, v.nbins, v.xmin, v.xmax , c_sr_of   , '', v.title, ofBin= False); h_tt_sr_of.Scale(1./h_tt_sr_of.Integral())
                h_tt_sr_sf.SetMarkerColor(r.kAzure-4); h_tt_sr_of.SetMarkerColor(r.kAzure+4); h_tt_sr_sf.SetMarkerSize(0.5); h_tt_sr_of.SetMarkerSize(0.5);
                print 'reading ttbar'
                if massCut == []:
                    self.tt_pdfs.extend([h_tt_sr_sf, h_tt_sr_of]); 
                    self.tt_cums.extend([h_tt_sr_sf.GetCumulative(), h_tt_sr_of.GetCumulative()])
                    setattr(self, h_tt_sr_sf.GetName(), h_tt_sr_sf);
                    setattr(self, h_tt_sr_of.GetName(), h_tt_sr_of)
                elif massCut ==['lepsMll_Edge < 90']: 
                    self.tt_pdfs_Zleq_90.extend([h_tt_sr_sf, h_tt_sr_of]); 
                    self.tt_cums_Zleq_90.extend([h_tt_sr_sf.GetCumulative(), h_tt_sr_of.GetCumulative()])
                    setattr(self, h_tt_sr_sf.GetName() + 'zleq90', h_tt_sr_sf); 
                    setattr(self, h_tt_sr_of.GetName() + 'zleq90', h_tt_sr_of)
                elif massCut == ['lepsMll_Edge > 90']:
                    self.tt_pdfs_Zgeq_90.extend([h_tt_sr_sf, h_tt_sr_of]);
                    self.tt_cums_Zgeq_90.extend([h_tt_sr_sf.GetCumulative(), h_tt_sr_of.GetCumulative()])
                    setattr(self, h_tt_sr_sf.GetName() + 'zgeq90', h_tt_sr_sf); 
                    setattr(self, h_tt_sr_of.GetName() + 'zgeq90', h_tt_sr_of)
                    
            ## the data
                print 'reading data'
                h_da_sr_of  = self.treeDA.getTH1F(1., 'da_%s_sr_of' %var, v.vtree, v.nbins, v.xmin, v.xmax , c_sr_of   , '', v.title, ofBin= False); h_da_sr_of.Scale(1./h_da_sr_of.Integral())
                h_da_sr_of.SetMarkerColor(r.kBlack); h_da_sr_of.SetMarkerSize(0.5); h_da_sr_of.SetMarkerStyle(34)
                if massCut == []:
                    self.da_pdfs.append(h_da_sr_of); self.da_cums.append(h_da_sr_of.GetCumulative())
                    setattr(self, h_da_sr_of.GetName(), h_da_sr_of)
                elif massCut ==['lepsMll_Edge < 90']: 
                    self.da_pdfs_Zleq_90.append(h_da_sr_of); self.da_cums_Zleq_90.append(h_da_sr_of.GetCumulative())
                    setattr(self, h_da_sr_of.GetName() + 'zleq90', h_da_sr_of)
                elif massCut ==['lepsMll_Edge > 90']: 
                    self.da_pdfs_Zgeq_90.append(h_da_sr_of); self.da_cums_Zgeq_90.append(h_da_sr_of.GetCumulative())
                    setattr(self, h_da_sr_of.GetName() + 'zgeq90', h_da_sr_of)
                print 'We have read', len(self.da_pdfs), ' data pdfs'
            ## the DY
                print 'reading dy'
                if self.doDY:
                    h_dy_sr_sf  = self.treeDY.getTH1F(1., 'dy_%s_sr_sf' %var, v.vtree, v.nbins, v.xmin, v.xmax , c_sr_sf   , '', v.title, ofBin= False); h_dy_sr_sf.Scale(1./h_dy_sr_sf.Integral())
                    h_dy_sr_sf.SetMarkerColor(r.kGreen-4); h_dy_sr_sf.SetMarkerSize(0.5)
                    if massCut == []:
                        self.dy_pdfs.append(h_dy_sr_sf); self.dy_cums.append(h_dy_sr_sf.GetCumulative())
                        setattr(self, h_dy_sr_sf.GetName(), h_dy_sr_sf)
                    elif massCut ==['lepsMll_Edge < 90']: 
                        self.dy_pdfs_Zleq_90.append(h_dy_sr_sf); self.dy_cums_Zleq_90.append(h_dy_sr_sf.GetCumulative())
                        setattr(self, h_dy_sr_sf.GetName() + 'zleq90', h_dy_sr_sf)
                    elif massCut ==['lepsMll_Edge > 90']: 
                        self.dy_pdfs_Zgeq_90.append(h_dy_sr_sf); self.dy_cums_Zgeq_90.append(h_dy_sr_sf.GetCumulative())
                        setattr(self, h_dy_sr_sf.GetName() + 'zgeq90', h_dy_sr_sf)

            ## the signal
                if self.doSig:
                    for i,mp in enumerate(self.mpoints):
                        c_sr_sf    = cuts.AddList([cuts.METJetsSignalRegionMET150 , cuts.GoodLeptonNoTriggerSF()] + massCut)
                        c_sr_of    = cuts.AddList([cuts.METJetsSignalRegionMET150 , cuts.GoodLeptonNoTriggerOF()] + massCut)
                        c_sr_sf_si = cuts.AddList([c_sr_sf.replace(cuts.twoLeptons, 't.nPairLep_Edge > 0'), 'GenSusyMScan1_Edge == %.0f && GenSusyMScan2_Edge == %.0f'%(mp[0], mp[1])])
                        h_si_sr_sf = self.treeSI.getTH1F(1., 'si_%s_sr_sf_m%dv%d' %(var,mp[0],mp[1]), v.vtree, v.nbins, v.xmin, v.xmax , c_sr_sf_si, '', v.title, ofBin = False); h_si_sr_sf.Scale(1./h_si_sr_sf.Integral())
                        h_si_sr_sf.SetMarkerColor(3+i); h_si_sr_sf.SetMarkerStyle(20+i); h_si_sr_sf.SetMarkerSize(0.5);
                        if massCut == []:
                            self.si_pdfs.append(h_si_sr_sf); self.si_cums.append(h_si_sr_sf.GetCumulative())
                            setattr(self, h_si_sr_sf.GetName(), copy.deepcopy(h_si_sr_sf))
                        elif massCut ==['lepsMll_Edge < 90']:
                            self.si_pdfs_Zleq_90.append(h_si_sr_sf); self.si_cums_Zleq_90.append(h_si_sr_sf.GetCumulative())
                            setattr(self, h_si_sr_sf.GetName() + 'zleq90', copy.deepcopy(h_si_sr_sf))
                        elif massCut ==['lepsMll_Edge > 90']:
                            self.si_pdfs_Zgeq_90.append(h_si_sr_sf); self.si_cums_Zgeq_90.append(h_si_sr_sf.GetCumulative())
                            setattr(self, h_si_sr_sf.GetName() + 'zgeq90', copy.deepcopy(h_si_sr_sf))

        self.treeHistsLoaded = True

    def getHisto(self, typ, var, data, normalize = True):
        h_name = ('em_data' if data else 'tt_mc' )+'_%s_histo_%s_ds_cuts_of_sr_met150'%('pdf' if typ == 'ker' else typ,var)
        if typ == 'ker':
            histo = copy.deepcopy(self.kerfile.Get(h_name));
            print histo.GetNbinsX()
        else:
            v = self.w.var( vParams(var).vtree.split('.')[1])
            histo = copy.deepcopy(self.w.pdf('%s_analyticalPDF_%s'%(var,'DA' if data else 'MC')).createHistogram(h_name, v, r.RooFit.Binning(1000) ))
            #histo = copy.deepcopy(w.pdf('%s_analyticalPDF_%s'%(var,'DA' if data else 'MC')).createHistogram(h_name, v, r.RooFit.Binning(50) ))
#        histo.Scale(1./histo.Integral())
        for i in range(1,histo.GetNbinsX()+1): 
            histo.SetBinError(i,0.)
        histo.SetLineColor(r.kBlack if data else r.kRed); histo.SetLineWidth(2)
        histo.SetLineStyle(1 if typ == 'fit' else 2)
        return histo

    def getLabel(self, obj):
        name = obj.GetName()
        ret = ''
        if   '_da'   in name: ret+='data'
        elif '_mc'   in name: ret+='MC'
        else: ret+=''
        if   '_fit_' in name: ret+=' - analytic pdf'
        elif '_pdf_' in name: ret+=' - kernel pdf'
        elif 'tt_'   in name: ret+='ttbar'
        elif 'si_'   in name: ret+='sig.'
        elif 'dy_'   in name: ret+='drell-yan'
        if   'ttbar' in ret and name.endswith('_of'): ret+=' OF'
        elif 'ttbar' in ret and name.endswith('_sf'): ret+=' SF'
        if 'sig'     in ret : ret+=' {x}/{y}'.format(x=name.split('_m')[-1].split('v')[0], y=name.split('_m')[-1].split('v')[1].split('_')[0])
        name.split('_m')[-1].split('v')
        return ret

    def loadPDFs(self):
#        for var,t in itertools.product((self.ivars+['p.d.f._m_et', 'p.d.f._l_dr','p.d.f._z_pt','p.d.f._m_lb']),self.typs):
        for var,t in itertools.product((self.ivars+['a3d','ldp']),self.typs):
            print var
            v = vParams(var)
            if var == 'p.d.f._m_et':
                name_ = 'met'
            elif var == 'p.d.f._l_dr':
                name_ = 'ldr'
            elif var == 'p.d.f._z_pt':
                name_ = 'zpt'
            elif var == 'p.d.f._m_lb':
                name_ = 'mlb'
            else:
                name_ = v.name

            name = v.name
            pdf_hist_da = self.getHisto(t, name_, 1); setattr(self, 'pdf_'+v.name+'_%s_da'%t, pdf_hist_da)
            pdf_hist_mc = self.getHisto(t, name_, 0); setattr(self, 'pdf_'+v.name+'_%s_mc'%t, pdf_hist_mc)
            cum_hist_da = pdf_hist_da.GetCumulative(); setattr(self, 'cum_'+v.name+'_%s_da'%t, cum_hist_da)
            cum_hist_mc = pdf_hist_mc.GetCumulative(); setattr(self, 'cum_'+v.name+'_%s_mc'%t, cum_hist_mc)
            self.all_pdfs.append(pdf_hist_da); self.all_pdfs.append(pdf_hist_mc);
            self.all_cums.append(cum_hist_da); self.all_cums.append(cum_hist_mc);
            if t == 'fit':
                self.ana_pdfs.extend([pdf_hist_da, pdf_hist_mc])
                self.ana_cums.extend([cum_hist_da, cum_hist_mc])
            else:
                self.ker_pdfs.extend([pdf_hist_da, pdf_hist_mc])
                self.ker_cums.extend([cum_hist_da, cum_hist_mc])
            print 'done'
    ## not yet def appendToRightList(self, ls):
    ## not yet     for i in ls:
    ## not yet         if all(x in i.GetName() for x in ['fit','pdf'])

    def makeNLLPDF(self, n=250000):
        print 'producing NLL pdfs. this might take a while depending on n. n is',n
        for t in self.typs:
            pdf_nll_da = self.fillLikelihood(t, 1, n, 'pdf_nll_%s_da'%t); setattr(self, 'pdf_nll_%s_da'%t, pdf_nll_da)
            pdf_nll_mc = self.fillLikelihood(t, 0, n, 'pdf_nll_%s_mc'%t); setattr(self, 'pdf_nll_%s_mc'%t, pdf_nll_mc)
            pdf_nll_da.SetLineColor(r.kBlack); pdf_nll_da.SetLineStyle(1 if t == 'fit' else 2);
            pdf_nll_mc.SetLineColor(r.kRed  ); pdf_nll_mc.SetLineStyle(1 if t == 'fit' else 2);
            if pdf_nll_da.GetNbinsX() > 100: pdf_nll_da.Rebin(20); pdf_nll_da.Scale(1/cum_nll_da.Integral())
            if pdf_nll_mc.GetNbinsX() > 100: pdf_nll_mc.Rebin(20); pdf_nll_mc.Scale(1/cum_nll_mc.Integral())
            cum_nll_da = pdf_nll_da.GetCumulative(); setattr(self, 'cum_nll_%s_da'%t, cum_nll_da)
            cum_nll_mc = pdf_nll_mc.GetCumulative(); setattr(self, 'cum_nll_%s_mc'%t, cum_nll_mc)
            self.all_pdfs.append(pdf_nll_da); self.all_pdfs.append(pdf_nll_mc);
            self.all_cums.append(cum_nll_da); self.all_cums.append(cum_nll_mc);
            if t == 'fit':
                self.ana_pdfs.extend([pdf_nll_da, pdf_nll_mc])
                self.ana_cums.extend([cum_nll_da, cum_nll_mc])
            else:
                self.ker_pdfs.extend([pdf_nll_da, pdf_nll_mc])
                self.ker_cums.extend([cum_nll_da, cum_nll_mc])
        self.hasNLLPDFs = True

    def fillLikelihood(self, typ, data, entries, name):
        print 'filling likelihoods'
        v = vParams('nll')
        histo  = self.treeTT.getTH1F(1., name, v.vtree, v.nbins, v.xmin, v.xmax , '0'   , '', v.title, ofBin= False)
        histo.Reset()
        datasets = []
        observs  = []
        pdfs     = []
        varnames = []
        for var in ['a3d', 'met', 'mlb', 'zpt']:
            pdf = self.w.pdf('%s_analyticalPDF_%s'%(var,'DA' if data else 'MC'))
            v   = self.w.var(vParams(var).vtree.split('.')[1])
            pdf.plotOn( v.frame() )
            obs = r.RooArgSet( v  )
            datasets.append( pdf.generate(obs, entries)  )
            observs .append( copy.deepcopy( obs ))
            pdfs    .append( pdf )
            varnames.append( vParams(var).vtree.split('.')[1])

 
        for i in range(entries):
            lh = 1.
            for dat in range(len(datasets)):
                rand = datasets[dat].get(i).getRealValue(varnames[dat] )
                self.w.var(varnames[dat]).setVal( rand )
                lh   = lh*pdfs[dat].getVal(observs[dat])
                #print varnames[dat], self.w.var(varnames[dat]).getVal(), pdfs[dat].getVal(observs[dat])
            histo.Fill(-1.*math.log(lh))
        histo.Scale(1/histo.Integral())
        for i in range(1,histo.GetNbinsX()+1): histo.SetBinError(i,0.)
        return copy.deepcopy(histo)


    def getLikelihood(self, typ, data):
        lh = 1.
        for histo in listofhistos:
            histo.Scale(1./histo.Integral())
            #if data == True  and 'tt_mc'   in histo.GetName(): continue
            #if data == False and 'em_data' in histo.GetName(): continue
            #random = histo.GetRandom(); 
            #lh = lh*histo.GetBinContent(histo.FindBin(random)) 
            lh *= histo.GetRandom()
        nll = -1.*math.log(lh)
        #print 'lh and nll', lh, nll
        return nll

    def plotZMassCuts(self):
        if not self.treeHistsLoaded: self.loadTreeHists(True)
        pdfs = self.tt_pdfs+self.si_pdfs+  self.dy_pdfs + self.da_pdfs + self.tt_pdfs_Zleq_90+self.si_pdfs_Zleq_90+  self.dy_pdfs_Zleq_90 + self.da_pdfs_Zleq_90 + self.tt_pdfs_Zgeq_90+self.si_pdfs_Zgeq_90+  self.dy_pdfs_Zgeq_90 + self.da_pdfs_Zgeq_90 
        print pdfs
        for sample in ['tt','si','da']:
            print 'I am in ', sample
            for j,mp in enumerate(self.mpoints):
                print 'I am in mp', mp[0]
                if not sample == 'si' and j > 0: continue
                    
        for var in self.ivars+['nll']:
            v = vParams(var)
            print 'I am at variable', v.name
            for chan in ['of','sf']:
                plot     = Canvas.Canvas('zMassCuts/pdfs_%s_tt_%s'%(v.name,chan), 'png,pdf', 0.65, 0.40, 0.85, 0.72)
                pdf_ratios = []
                ind = 0
                for i,pdf in enumerate(pdfs):
                    if not '_'+v.name+'_' in pdf.GetName(): continue
                    if not 'tt' in pdf.GetName(): continue
                    if not chan in pdf.GetName(): continue
                    print pdf.GetName()
                    if ind == 0:
                        label = 'No mass cut'
                        color = r.kBlue
                    elif ind == 1:
                        label = 'M_{ll} < 90 GeV'
                        pdf_ratios.append(pdf)
                        color = r.kRed
                    elif ind == 2:
                        label = 'M_{ll} > 90 GeV'
                        pdf_ratios.append(pdf)
                        color = r.kBlack
                    plot   .addHisto(pdf, 'pe same', label, 'PL', color, 1, ind)
                    ind = ind+1
                ratioName = '%s_%s_sr_%s'%('tt',var,chan)
                plot   .saveRatio(1, 0, 0 , 1., getattr(self, ratioName), pdf_ratios, 0.5, 1.5)
                
#            plot     = Canvas.Canvas('zMassCuts/pdfs_%s_da_%s'%(v.name,'of'), 'png,pdf', 0.65, 0.40, 0.85, 0.72)
#            pdf_ratios = []
#            ind = 0
#            for i,pdf in enumerate(pdfs):
#                if not '_'+v.name+'_' in pdf.GetName(): continue
#                if not 'da' in pdf.GetName(): continue
#                print pdf.GetName()
#                if ind == 0:
#                    label = 'No mass cut'
#                    color = r.kBlue
#                elif ind == 1:
#                    label = 'M_{ll} < 90 GeV'
#                    pdf_ratios.append(pdf)
#                    color = r.kRed
#                elif ind == 2:
#                    label = 'M_{ll} > 90 GeV'
#                    pdf_ratios.append(pdf)
#                    color = r.kBlack
#                plot   .addHisto(pdf, 'pe same', label, 'PL', color, 1, ind)
#                ind = ind+1
#                ratioName = '%s_%s_sr_%s'%('da',var,'of')
#                plot   .saveRatio(1, 0, 0 , 1., getattr(self, ratioName), pdf_ratios, 0.5, 1.5)

                for mp in self.mpoints:
                    plot     = Canvas.Canvas('zMassCuts/pdfs_%s_si_%dv%d'%(v.name,mp[0],mp[1]), 'png,pdf', 0.65, 0.40, 0.85, 0.72)
                    pdf_ratios = []
                    ind = 0
                    for i,pdf in enumerate(pdfs):
                        if not '_'+v.name+'_' in pdf.GetName(): continue
                        if not 'si' in pdf.GetName(): continue
                        if not '%dv%d'%(mp[0],mp[1]) in pdf.GetName(): continue
                        print pdf.GetName()
                        if ind == 0:
                            label = 'No mass cut'
                            color = r.kBlue
                        elif ind == 1:
                            label = 'M_{ll} < 90 GeV'
                            pdf_ratios.append(pdf)
                            color = r.kRed
                        elif ind == 2:
                            label = 'M_{ll} > 90 GeV'
                            pdf_ratios.append(pdf)
                            color = r.kBlack
                        plot   .addHisto(pdf, 'pe same', label, 'PL', color, 1, ind)
                        ind = ind+1
                    ratioName = 'si_%s_sr_sf_m%dv%d'%(v.name, mp[0],mp[1])
                    print getattr(self, ratioName), pdf_ratios
                    plot   .saveRatio(1, 0, 0 , 1., getattr(self, ratioName), pdf_ratios, 0.5, 1.5)

            
    def profile_nll_zmass(self):
        c_sr_of    = cuts.AddList([cuts.METJetsSignalRegionMET150 , cuts.GoodLeptonNoTriggerOF()])
        c_sr_sf    = cuts.AddList([cuts.METJetsSignalRegionMET150 , cuts.GoodLeptonNoTriggerSF()])
        histo_of = self.treeTT.getTH2F(1., 'prof_nll_Zmass_of', '(-1.*TMath::Log(lh_ana_met_data_Edge*lh_ana_mlb_data_Edge*lh_ana_ldr_data_Edge*lh_ana_zpt_data_Edge)):lepsMll_Edge', 18, 20., 200.,50,15.,35., c_sr_of   , '', "M_{ll} GeV","nll")
        histo_sf = self.treeTT.getTH2F(1., 'prof_nll_Zmass_sf', '(-1.*TMath::Log(lh_ana_met_data_Edge*lh_ana_mlb_data_Edge*lh_ana_ldr_data_Edge*lh_ana_zpt_data_Edge)):lepsMll_Edge', 18, 20., 200.,50,15.,35., c_sr_sf   , '', "M_{ll} GeV","nll")
        profile_of = histo_of.ProfileX(); mean_of = profile_of.GetSumOfWeights() / profile_of.GetNbinsX()
        profile_sf = histo_sf.ProfileX(); mean_sf = profile_sf.GetSumOfWeights() / profile_sf.GetNbinsX()
        h_of       = profile_of.ProjectionX()
        h_sf       = profile_sf.ProjectionX()

        profile_of.GetYaxis().SetRangeUser(0.,25.)
        plot     = Canvas.Canvas('zMassCuts/profile_all', 'png,pdf', 0.25, 0.20, 0.45, 0.42)
        plot   .addHisto(profile_of, 'same', 'ttbar of', 'L', r.kRed,  1, 0)
        plot   .addHisto(profile_sf, 'same', 'ttbar sf', 'L', r.kBlue, 1, 1)
        plot   .addLine(20.,mean_of,200.,mean_of, r.kRed)
        plot   .addLine(20.,mean_sf,200.,mean_sf, r.kBlue)
        line_OF = TH1F('line_OF','',18,20.,200.)
        for bin in range(line_OF.GetNbinsX()):
            line_OF.SetBinContent(bin+1, mean_of)
            line_OF.SetBinError(bin+1, 0.)
        plot   .saveRatio(1,0,0,1.,h_of, line_OF,0.9,1.1)
        

    def kernelAnalyticComparison(self, plotMC = False, plotData = False):
        if not self.hasNLLPDFs: self.makeNLLPDF() # call this only once!!!!
        if plotMC and not self.treeHistsLoaded: self.loadTreeHists()
        print 'Len of all_pdfs',len(self.all_pdfs)
        #for var in self.ivars+['nll','ldp']:
        #for var in ['nll','ldp']:
        for var in ['mt2']:
            v = vParams(var)
            print 'I am at variable', v.name
            plot     = Canvas.Canvas('kernelAnalyliticWithData%s/pdfs_%s%s'%('_withSig' if self.doSig else '',v.name,'_withMC' if plotMC else '')       , 'png,pdf', 0.65, 0.40, 0.85, 0.72)
            plotCum  = Canvas.Canvas('kernelAnalyliticWithData%s/cumulatives_%s%s'%('_withSig' if self.doSig else '',v.name,'_withMC' if plotMC else ''), 'png,pdf', 0.65, 0.40, 0.85, 0.72)

            ind=0; histos=[]

            
            for i,_pdf in enumerate(self.all_pdfs):
                #cum = self.all_cums[i]
                _pdf.Scale(1./_pdf.Integral())
                if _pdf.GetNbinsX() > 900:
                    pdf = _pdf.Rebin(20); 
                else:
                    pdf = _pdf
                cum = pdf.GetCumulative()
                #print 'number of bins of the cumulative distribution', cum.GetNbinsX()
                if not '_'+v.name+'_' in pdf.GetName(): continue
                pdf.GetXaxis().SetTitle(v.title); cum.GetXaxis().SetTitle(v.title)
                pdf.GetYaxis().SetRangeUser(0., v.ymax)
                cum.GetYaxis().SetRangeUser(0., 1.05  )
                if var == 'nll':
                    print pdf.GetXaxis().GetXmax(), pdf.GetYaxis().GetXmin(), pdf.GetNbinsX()
                pdf.SetLineWidth(2); cum.SetLineWidth(2)
                pdf.SetDrawOption('hist,same,l'); cum.SetDrawOption('hist,same,l')
                pdf.GetYaxis().SetTitle('Prob. density / %4.2f'%pdf.GetBinWidth(1))
                plot   .addHisto(pdf, 'hist same l', self.getLabel(pdf), 'L', pdf.GetLineColor(), 1, ind, doOF = not var == 'ldp')
                print pdf.GetXaxis().GetTitle()
                plotCum.addHisto(cum, 'hist same l', self.getLabel(cum), 'L', cum.GetLineColor(), 1, ind)
                #print '%s cum from pdf histo range: min %.1f max %.1f nbins %.1f'%(v.name, cum.GetXaxis().GetXmin(), cum.GetXaxis().GetXmax(), cum.GetNbinsX())
                #print '%s pdf from pdf histo range: min %.1f max %.1f nbins %.1f'%(v.name, pdf.GetXaxis().GetXmin(), pdf.GetXaxis().GetXmax(), pdf.GetNbinsX())
                ind+=1
#            pdf_ratios = [getattr(self, 'pdf_'+v.name+'_fit_mc'), getattr(self, 'pdf_'+v.name+'_ker_da'), getattr(self, 'pdf_'+v.name+'_ker_mc')]
#            cum_ratios = [getattr(self, 'cum_'+v.name+'_fit_mc'), getattr(self, 'cum_'+v.name+'_ker_da'), getattr(self, 'cum_'+v.name+'_ker_mc')]
            pdf_ratios = []; cum_ratios = []
            
            pdf_ratios = [getattr(self, 'pdf_'+v.name+'_fit_mc'),getattr(self, 'pdf_'+v.name+'_fit_da')]
            cum_ratios = [getattr(self, 'cum_'+v.name+'_fit_mc'),getattr(self, 'cum_'+v.name+'_fit_da')]

            if plotMC:
                for i,pdf in enumerate(self.tt_pdfs+( [] if not self.doSig else self.si_pdfs) +([] if not self.doDY else self.dy_pdfs)):
                    cum = pdf.GetCumulative(); cum.SetLineColor(pdf.GetMarkerColor()); cum.SetMarkerColor(pdf.GetMarkerColor())
                    if not '_'+v.name+'_' in pdf.GetName(): continue
                    pdf.SetDrawOption('pe same'); cum.SetDrawOption('pe same')
                    pdf.GetYaxis().SetTitle('Prob. density / %4.2f'%pdf.GetBinWidth(1))
                    plot   .addHisto(pdf, 'pe same', self.getLabel(pdf), 'PL', pdf.GetMarkerColor(), 1, ind, doOF = not var == 'ldp')
                    print pdf.GetXaxis().GetTitle()
                    plotCum.addHisto(cum, 'pe same', self.getLabel(cum), 'PL', cum.GetMarkerColor(), 1, ind)
                    ind+=1
                    if 'tt_' in pdf.GetName() and not 'of' in pdf.GetName():
                        pdf_ratios.append(pdf); cum_ratios.append(cum)
            if plotData:
                for i,pdf in enumerate(self.da_pdfs):
                    cum = pdf.GetCumulative(); cum.SetLineColor(pdf.GetMarkerColor()); cum.SetMarkerColor(pdf.GetMarkerColor())
                    if not '_'+v.name+'_' in pdf.GetName(): continue
                    #print 'plotting', pdf.GetName()
                    if var == 'nll':
                        print pdf.GetXaxis().GetXmax(), pdf.GetYaxis().GetXmin(), pdf.GetNbinsX()
                    pdf.SetMarkerStyle(r.kFullCircle); cum.SetMarkerStyle(r.kFullCircle)
                    #pdf.SetMarkerSize(0.5); cum.SetMarkerSize(0.5)
                    pdf.SetDrawOption('pe same'); cum.SetDrawOption('pe same')
                    pdf.GetYaxis().SetTitle('Prob. density / %4.2f'%pdf.GetBinWidth(1))
                    plot   .addHisto(pdf, 'pe same', 'Data', 'P', pdf.GetMarkerColor(), 1, ind, doOF = not var == 'ldp')
                    print pdf.GetXaxis().GetTitle()
                    plotCum.addHisto(cum, 'pe same', 'Data', 'P', cum.GetMarkerColor(), 1, ind)#, doOF = True)
                    if 'da_' in pdf.GetName():
                        pdf_ratios.append(pdf); cum_ratios.append(cum)
                    ind+=1
            for h in cum_ratios:
                h.Scale( 1 / h.GetBinContent( h.GetNbinsX() ))
            getattr(self, 'cum_'+v.name+'_fit_da').Scale(1 / getattr(self, 'cum_'+v.name+'_fit_da').GetBinContent( getattr(self, 'cum_'+v.name+'_fit_da').GetNbinsX()))
                    
            plot   .saveRatio(1, 1, True, 2.3, getattr(self, 'pdf_'+v.name+'_fit_da'), pdf_ratios, 0., 2.)
            plotCum.saveRatio(1, 1, True, 2.3, getattr(self, 'cum_'+v.name+'_fit_da'), cum_ratios, 0., 2.)
#            plot   .saveRatio(1, 0, 0 , 2.3, getattr(self, 'tt_%s_sr_of'%var),                 pdf_ratios, 0., 2.)
#            plotCum.saveRatio(1, 0, 0 , 2.3, getattr(self, 'tt_%s_sr_of'%var).GetCumulative(), cum_ratios, 0., 2.)
            
             

    def getROC(self):
        color = [r.kRed,r.kBlue,r.kGreen,r.kMagenta, r.kBlack, r.kOrange]
        c_sr_sf    = cuts.AddList([cuts.METJetsSignalRegionMET150 , cuts.GoodLeptonNoTriggerSF()])
        npoints = 1000
        for i,mp in enumerate(self.mpoints):
            c_sr_sf_si = cuts.AddList([c_sr_sf.replace(cuts.twoLeptons, 't.nPairLep_Edge > 0'), 'GenSusyMScan1_Edge == %.0f && GenSusyMScan2_Edge == %.0f'%(mp[0], mp[1])])
            ind = 0
            plot     = Canvas.Canvas('rocCurves/roc_%d_%d'%(mp[0],mp[1]), 'png,pdf', 0.2, 0.25, 0.4, 0.48)
            yline90 = 0.
            xline90 = 0.
            nll90   = 0.
#            for var in ['nll', 'nll_noldr', 'nll_a3d']:
            for var in ['nll_noldr', 'nll_nozpt','nll_nomlb','nll_nomet','nll_ldp']:
                v = vParams(var)
                print 'I am at variable', 
                h_tt_sr_sf  = self.treeTT.getTH1F(1., 'tt_%s' %var, v.vtree, npoints, v.xmin, v.xmax , c_sr_sf   , '', v.title, ofBin= True); h_tt_sr_sf.Scale(1./h_tt_sr_sf.Integral())
                h_tt_cum    = h_tt_sr_sf.GetCumulative()                     
                h_si_sr_sf  = self.treeSI.getTH1F(1., 'si_%s_sr_sf_m%dv%d' %(var,mp[0],mp[1]), v.vtree, npoints, v.xmin, v.xmax , c_sr_sf_si, '', v.title, ofBin = True); h_si_sr_sf.Scale(1./h_si_sr_sf.Integral())
                h_si_cum    = h_si_sr_sf.GetCumulative()
                gr          = TGraph(npoints)
                
                gr.SetPoint(0,0.,1.)
                for bin in range(npoints):
                    gr.SetPoint(bin+1, h_tt_cum.GetBinContent(bin+1),1-h_si_cum.GetBinContent(bin+1))
                    if var == 'nll_ldp':
                        if math.pow((xline90-0.95),2) > math.pow(h_tt_cum.GetBinContent(bin+1)-0.95,2):
                            yline90 = 1-h_si_cum.GetBinContent(bin+1)
                            xline90 = h_tt_cum.GetBinContent(bin+1)
                            nll90   = h_tt_cum.GetBinCenter(bin+1)
                    print bin, h_tt_cum.GetBinContent(bin+1),1-h_si_cum.GetBinContent(bin+1)
                    
                gr.GetXaxis().SetTitle('t#bar{t} rejection')
                gr.GetYaxis().SetTitle('Signal efficiency')
                gr.SetLineWidth(3)
                gr.GetXaxis().SetRangeUser(0.,1.)
                print 'going to plot', gr.GetXaxis().GetXmin()
                if ind == 0:
                    plot.addGraph(gr,"al",v.title,"l", color[ind],1,ind)
                else:
                    plot.addGraph(gr,"lsame",v.title,"l", color[ind],1,ind)
                ind = ind+1
            print xline90, yline90, nll90
            plot.addLine(0.95,0.,0.95,yline90,r.kRed,2)
            plot.addLine(0.0,yline90,xline90,yline90,r.kRed,2)
            plot.addLatex(0.2,0.2,'m_{#tilde{b}} = %d GeV, m_{#tilde{#chi}_{2}} = %d GeV'%(mp[0],mp[1]))
            plot.save(1,False,0,2.3)

    def getMT2NllComparison(self):
        color = [r.kRed,r.kBlue]
        c_sr_sf    = cuts.AddList([cuts.METJetsSignalRegionMET150 , cuts.GoodLeptonNoTriggerSF()])
        npoints = 1000
        for i,mp in enumerate(self.mpoints):
            c_sr_sf_si = cuts.AddList([c_sr_sf.replace(cuts.twoLeptons, 't.nPairLep_Edge > 0'), 'GenSusyMScan1_Edge == %.0f && GenSusyMScan2_Edge == %.0f'%(mp[0], mp[1])])
            ind = 0
            plot     = Canvas.Canvas('MT2Nll/roc_%d_%d'%(mp[0],mp[1]), 'png,pdf', 0.2, 0.25, 0.4, 0.48)
            yline90 = 0.
            xline90 = 0.
            nll90   = 0.
            for var in ['nll_ldp','mt2']:
                v = vParams(var)
                print 'I am at variable', 
                h_tt_sr_sf  = self.treeTT.getTH1F(1., 'tt_%s' %var, v.vtree, npoints, v.xmin, v.xmax , c_sr_sf   , '', v.title, ofBin= True); h_tt_sr_sf.Scale(1./h_tt_sr_sf.Integral())
                h_tt_cum    = h_tt_sr_sf.GetCumulative()                     
                h_si_sr_sf  = self.treeSI.getTH1F(1., 'si_%s_sr_sf_m%dv%d' %(var,mp[0],mp[1]), v.vtree, npoints, v.xmin, v.xmax , c_sr_sf_si, '', v.title.split()[0], ofBin = True); h_si_sr_sf.Scale(1./h_si_sr_sf.Integral())
                h_si_cum    = h_si_sr_sf.GetCumulative()
                gr          = TGraph(npoints)
                
                gr.SetPoint(0,0.,1.)
                for bin in range(npoints):
                    gr.SetPoint(bin+1, h_tt_cum.GetBinContent(bin+1),1-h_si_cum.GetBinContent(bin+1))
                    if var == 'nll_a3d':
                        if math.pow((xline90-0.9),2) > math.pow(h_tt_cum.GetBinContent(bin+1)-0.9,2):
                            yline90 = 1-h_si_cum.GetBinContent(bin+1)
                            xline90 = h_tt_cum.GetBinContent(bin+1)
                            nll90   = h_tt_cum.GetBinCenter(bin+1)
                    print bin, h_tt_cum.GetBinContent(bin+1),1-h_si_cum.GetBinContent(bin+1)
                    
                gr.GetXaxis().SetTitle('t#bar{t} rejection')
                gr.GetYaxis().SetTitle('Signal efficiency')
                gr.SetLineWidth(3)
                gr.GetXaxis().SetRangeUser(0.,1.)
                print 'going to plot', gr.GetXaxis().GetXmin()
                if ind == 0:
                    plot.addGraph(gr,"al",v.title.split()[0],"l", color[ind],1,ind)
                else:
                    plot.addGraph(gr,"lsame",v.title.split()[0],"l", color[ind],1,ind)
                ind = ind+1
            print xline90, yline90, nll90
            plot.addLine(0.9,0.,0.9,yline90,r.kRed,2)
            plot.addLine(0.0,yline90,xline90,yline90,r.kRed,2)
            plot.addLatex(0.2,0.2,'m_{#tilde{b}} = %d GeV, m_{#tilde{#chi}_{2}} = %d GeV'%(mp[0],mp[1]))
            plot.save(1,False,0,2.3)

    def compareOFSF(self):
        c_sr_sf    = cuts.AddList([cuts.METJetsSignalRegionMET150 , cuts.GoodLeptonNoTriggerSF()])
        c_sr_of    = cuts.AddList([cuts.METJetsSignalRegionMET150 , cuts.GoodLeptonNoTriggerOF()])
        c_sr_ee    = cuts.AddList([cuts.METJetsSignalRegionMET150 , cuts.GoodLeptonNoTriggeree()])
        c_sr_mm    = cuts.AddList([cuts.METJetsSignalRegionMET150 , cuts.GoodLeptonNoTriggermm()])
        for var in ['nll','met','zpt','mlb','ldr']:
            plot     = Canvas.Canvas('compareOFSF/%s'%var, 'png,pdf', 0.65, 0.40, 0.85, 0.72)
            plot_eemm= Canvas.Canvas('compareOFSF/eemm_%s'%var, 'png,pdf', 0.65, 0.40, 0.85, 0.72)
            v = vParams(var);
            h_tt_sr_sf  = self.treeTT.getTH1F(1., 'tt_sf_%s' %var, v.vtree, v.nbins, v.xmin, v.xmax , c_sr_sf   , '', v.title, ofBin= False); h_tt_sr_sf.Scale(1./h_tt_sr_sf.Integral())
            h_tt_sr_of  = self.treeTT.getTH1F(1., 'tt_of_%s' %var, v.vtree, v.nbins, v.xmin, v.xmax , c_sr_of   , '', v.title, ofBin= False); h_tt_sr_of.Scale(1./h_tt_sr_of.Integral())
            h_tt_sr_ee  = self.treeTT.getTH1F(1., 'tt_of_%s' %var, v.vtree, v.nbins, v.xmin, v.xmax , c_sr_ee   , '', v.title, ofBin= False); h_tt_sr_ee.Scale(1./h_tt_sr_ee.Integral())
            h_tt_sr_mm  = self.treeTT.getTH1F(1., 'tt_of_%s' %var, v.vtree, v.nbins, v.xmin, v.xmax , c_sr_mm   , '', v.title, ofBin= False); h_tt_sr_mm.Scale(1./h_tt_sr_mm.Integral())
            plot      .addHisto(h_tt_sr_sf, 'hist same', 't#bar{t} - SF', 'L', r.kRed, 1, 0)
            plot     .addHisto(h_tt_sr_of, 'hist same', 't#bar{t} - OF', 'L', r.kBlue,1, 1)
            plot_eemm.addHisto(h_tt_sr_ee, 'hist same', 't#bar{b} - EE', 'L', r.kRed, 1,0)
            plot_eemm.addHisto(h_tt_sr_mm, 'hist same', 't#bar{b} - MM', 'L', r.kBlack, 1,0)
            plot_eemm.addHisto(h_tt_sr_of, 'hist same', 't#bar{t} - OF', 'L', r.kBlue,1, 1)
            plot     .saveRatio(1, 0, 0 , 1., h_tt_sr_of, h_tt_sr_sf, 0.5, 1.5)
            plot_eemm.saveRatio(1, 0, 0 , 1., h_tt_sr_of, [h_tt_sr_mm, h_tt_sr_ee], 0.5, 1.5)        

    def varsAfterCut(self):
        c_sr_sf    = cuts.AddList([cuts.METJetsSignalRegionMET150 , cuts.GoodLeptonNoTriggerSF(), '-1.*TMath::Log(lh_ana_met_data_Edge*lh_ana_mlb_data_Edge*lh_ana_a3d_data_Edge*lh_ana_zpt_data_Edge) > 19.9' ])
        c_sr_sf_nocuts    = cuts.AddList([cuts.METJetsSignalRegionMET150 , cuts.GoodLeptonNoTriggerSF()])
        for var in ['met','zpt','mlb','a3d']:
            print 'plotting', var
            plot     = Canvas.Canvas('varsAfterCut/%s'%var, 'png,pdf', 0.65, 0.40, 0.85, 0.72)
            v = vParams(var);
            histo  = self.treeTT.getTH1F(1., 'histo_%s' %var, v.vtree, v.nbins, v.xmin, v.xmax , c_sr_sf   , '', v.title, ofBin= False)
            histo_nocuts  = self.treeTT.getTH1F(1., 'histo_%s' %var, v.vtree, v.nbins, v.xmin, v.xmax , c_sr_sf_nocuts   , '', v.title, ofBin= False)
            plot      .addHisto(histo_nocuts, 'hist same', 't#bar{t} ', 'L', r.kBlue, 1, 0)
            plot      .addHisto(histo, 'hist same', 't#bar{t} - nll Cut', 'L', r.kRed, 1, 0)

            plot.save(1,0,0,2.3)

    def getChi2ForFits(self):
#        if not self.hasNLLPDFs: self.makeNLLPDF() # call this only once!!!!
        if not self.treeHistsLoaded: self.loadTreeHists()

        for var in ['zpt','mlb','a3d','met']:
            v   = self.w.var( vParams(var).vtree.split('.')[1])
            pdf = getattr(self, 'pdf_'+vParams(var).name+'_fit_da'); pdf.Rebin(20); pdf.Scale(1/pdf.Integral())
            dat = getattr(self,'da_%s_sr_of'%var);                   dat.Scale(1/dat.Integral())
            c = TCanvas()
            pdf.Draw()
            dat.Draw('same')
            c.SaveAs(var+'.pdf')
            print pdf.GetNbinsX(),  dat.GetNbinsX(), dat.GetXaxis().GetXmin(), pdf.GetXaxis().GetXmin(), dat.GetXaxis().GetXmax(), pdf.GetXaxis().GetXmax()
            print 'chi2 is', pdf.KolmogorovTest(dat)

    def CheckCorrelations(self):
        cuts = CutManager.CutManager()
        theCuts = cuts.AddList([cuts.METJetsSignalRegionMET150 , cuts.GoodLeptonNoTriggerSF()])
        for var1, var2 in itertools.product(['met','a3d','zpt','mlb'],['met','a3d','zpt','mlb']):
            v1 = vParams(var1); v2 = vParams(var2)
            histo = self.treeTT.getTH2F(1., '%s%s'%(var1,var2), '%s:%s'%(v1.vtree,v2.vtree), v2.nbins, v2.xmin, v2.xmax,v1.nbins,v1.xmin,v1.xmax, theCuts , '', "","")
            print var1, var2, histo.GetCorrelationFactor()
            
    def getDYEff(self):
        cuts = CutManager.CutManager()
        theCuts = cuts.AddList([cuts.METJetsSignalRegionMET150 , cuts.GoodLeptonNoTriggerSF()])
        theCuts = theCuts.replace(cuts.twoLeptons, 't.nPairLep_Edge > 0')
        self.loadTrees()
        (denom,err1) = self.treeDY.getYields(2.3,'1',-5,5,theCuts)
        (num,err2)   = self.treeDY.getYields(2.3,'1',-5,5,theCuts + '&&' + cuts.nLL90Cut())
        print denom, '+/-', err1
        print num, '+/-', err2
        print 'drell yan efficiency for the nll cut is', num/denom
        return num/denom

    def makePrediction(self,of_histo, ing, eta):
        central = (eta == 'Central')
        sf_pred = copy.deepcopy(of_histo)
        scale = ing.rsfof_final_cen_val  if central else ing.rsfof_final_fwd_val
        s_err = ing.rsfof_final_cen_err  if central else ing.rsfof_final_fwd_err
#        print '========================================='
#        print 'using rsfof for %s in %s: %.4f +- %.4f' %(ing.dType, eta, scale, s_err)
#        print '========================================='
        for _bin in range(1,of_histo.GetNbinsX()+1):
            tmp_cont  = of_histo.GetBinContent(_bin)
            tmp_err   = of_histo.GetBinError  (_bin)
            tmp_mass  = of_histo.GetXaxis().GetBinCenter(_bin)
            tmp_pred  = scale*tmp_cont
            tmp_prede = math.sqrt(tmp_err*tmp_err + s_err*s_err*tmp_cont*tmp_cont)
            print 'result', tmp_pred, 'error', tmp_prede
            sf_pred.SetBinContent(_bin, tmp_pred )
            sf_pred.SetBinError  (_bin, tmp_prede)
        return sf_pred


    def dataCard(self):
        self.loadTrees()
        doData = True
        theTree = self.treeDA if doData else self.treeMC 
        cuts = CutManager.CutManager()
        theCuts = [cuts.METJetsSignalRegionMET150, cuts.nLL90Cut()]
        theCutsNoNLL = [cuts.METJetsSignalRegionMET150, cuts.antinLL90Cut()]
        theIng  = self.ingDA if doData else self.ingMC
        dy_eff  = self.getDYEff()
        print dy_eff
        setattr(self,'dy_eff_ttbar', 1 - dy_eff)
        setattr(self,'dy_eff_sig'  , dy_eff)

        ##############################################################################
##############################################################################
        print 'esto hay que uitarlo en algum momento'
        dy_3_150= self.treeDY.getTH1F(lumi, 'pass', '0.5', 1, 0,1 , cuts.AddList([cuts.GoodLeptonSF(),'nJetSel_Edge > 2 && met_Edge > 150']), '', 'M_{ll}', ofBin= False).GetBinContent(1) # lo que quiero
        dy_3_100= self.treeDY.getTH1F(lumi, 'pass', '0.5', 1, 0,1 , cuts.AddList([cuts.GoodLeptonSF(),'nJetSel_Edge > 2 && met_Edge > 100']), '', 'M_{ll}', ofBin= False).GetBinContent(1) #lo que tengo

        jets2_dy_Central = 15.7;       jets2_dy_Central_e = 4.5 
        jets3_dy_Central = 9.7;        jets3_dy_Central_e = 1.6

        jets2_dy_Forward = 5.3;        jets2_dy_Forward_e = 1.4  
        jets3_dy_Forward = 4.0;        jets3_dy_Forward_e = 0.7



        self.dy_onz_Central   = jets2_dy_Central + jets3_dy_Central
        self.dy_onz_Central_e = math.sqrt(jets3_dy_Central_e**2 + jets2_dy_Central_e**2)
        self.dy_onz_Forward   = jets2_dy_Forward + jets3_dy_Forward
        self.dy_onz_Forward_e = math.sqrt(jets3_dy_Forward_e**2 + jets2_dy_Forward_e**2)
##############################################################################
        ##############################################################################


        
        print 'drell yan estimation off shell'

        setattr(self,'dy_low',(getattr(theIng.rinout_dy_lm,'val') + getattr(theIng.rinout_dy_bz,'val'))*getattr(self,'dy_onz'))
        setattr(self,'dy_low_e',math.sqrt(((getattr(theIng.rinout_dy_lm,'val') + getattr(theIng.rinout_dy_bz,'val'))*getattr(self,'dy_onz_e'))**2 + (math.sqrt(getattr(theIng.rinout_dy_lm,'err')**2 + getattr(theIng.rinout_dy_bz,'err')**2)*getattr(self,'dy_onz_e'))**2))
        setattr(self,'dy_high',(getattr(theIng.rinout_dy_az,'%s_val'%et) + getattr(theIng.rinout_dy_hm,'%s_val'%et))*getattr(self,'dy_onz_%s'%eta))
        setattr(self,'dy_high_e',math.sqrt(((getattr(theIng.rinout_dy_az,'val') + getattr(theIng.rinout_dy_hm,'val'))*getattr(self,'dy_onz_e'))**2 + (math.sqrt(getattr(theIng.rinout_dy_az,'err')**2 + getattr(theIng.rinout_dy_hm,'err')**2)*getattr(self,'dy_onz_e'))**2))
        print 'low',  getattr(self,'dy_low' ), '+/-',  getattr(self,'dy_low_e' )
        print 'high', getattr(self,'dy_high'), '+/-',  getattr(self,'dy_high_e')



 
        for nll in ['sig','ttbar']: 
            vH_obs = []; vH_of = []; vH_pred = []; vH_DY = []
                
            theCuts = [cuts.METJetsSignalRegionMET150, cuts.nLL90Cut() if 'sig' in nll else cuts.antinLL90Cut()]
            H_obs  = theTree.getTH1F(lumi, 'obs%s'%(nll), 'lepsMll_Edge', [0.,81,101,5000.], 0,0 , cuts.AddList(theCuts+[cuts.GoodLeptonSF()]), '', 'M_{ll}', ofBin= False)
            H_of   = theTree.getTH1F(lumi, 'of%s' %(nll), 'lepsMll_Edge', [0.,81,101,5000.], 0,0 , cuts.AddList(theCuts+[cuts.GoodLeptonOF()]), '', 'M_{ll}', ofBin= False)
                #print 'dy is being predicted with FULLY MC'
                H_pred= self.makePrediction(H_of,theIng, eta)
                H_dy = H_obs.Clone('hdy%s'%(nll)); H_dy.Reset()
                for mll in ['low','onz','high']:
                    bin_numb = 1 if 'low' in mll else 2 if 'onz' in mll else 3
                    H_dy.SetBinContent( bin_numb, getattr(self,'dy_%s'  %(mll)))
                    H_dy.SetBinError  ( bin_numb, getattr(self,'dy_%s_e'%(mll)))
                H_dy.Scale(getattr(self,'dy_eff_%s'%nll))
                                        
        #        vH_obs.append(H_obs); vH_of.append(H_of); vH_pred.append(H_pred); vH_DY.append(H_dy)

        #    vH_obs [0].Add(vH_obs [1]); H_obs  = vH_obs [0]
        #    vH_of  [0].Add(vH_of  [1]); H_of   = vH_of  [0]
        #    vH_pred[0].Add(vH_pred[1]); H_pred = vH_pred[0]
        #    vH_DY  [0].Add(vH_DY  [1]); H_dy   = vH_DY  [0]
            output = '\t low \t onz \t high\n'
            output = output + nll + '\t'
            dy = []; Obs=[]; Fs = []
            dyE = []; FsE = []
            for mll in ['low','onz','high']:
                bin_name = 'inclusive_'+ mll + '_incb_' + nll
                bin_numb = 1 if 'low' in mll else 2 if 'onz' in mll else 3
                
                obs      = H_obs .GetBinContent( bin_numb )
                fs       = H_pred.GetBinContent( bin_numb )
                fs_unc   = H_pred.GetBinError  ( bin_numb ) / fs
                onz      = H_dy  .GetBinContent( bin_numb )
                onz_e    = H_dy  .GetBinError  ( bin_numb )
                of_yield = H_of  .GetBinContent( bin_numb )
                rsfof    = fs / of_yield
                Obs.append(obs); Fs.append(fs); dy.append(max(0.,onz))
                dyE.append(onz_e); FsE.append(H_pred.GetBinError( bin_numb))
                dc = '''# this is the datacard for bin {bin_name}
imax 1  number of channels
jmax 2  number of backgrounds
kmax *  number of nuisance parameters (sources of systematical uncertainties)
------------
# only one channel here with observation {obs}
bin            {bin_name}
observation    {obs}
------------
bin        {bin_name}     {bin_name}     {bin_name}
process    sig            FS             DY    
process    0              1              2
rate       XXRATEXX         {fs_bkg:.2f}       {dy_bkg:.2f}
------------
deltaS       lnN              1.15      -           -       
lumi         lnN              1.12      -           -       
FS_unc       lnN              -         {fs_unc:.2f}    -       
{bin_name}_fs_stat      gmN {of_yield}   -         {rsfof:.3f}        -       
DY_unc       lnN              -         -           {dy_unc:.2f}'''.format(bin_name=bin_name,
                                                                           obs=int(obs),#
                                                                           fs_bkg=fs,#
                                                                           fs_unc=1+fs_unc,#
                                                                           dy_bkg=max(0.,onz),#
                                                                           dy_unc=1+onz_e/abs(onz),
                                                                           of_yield=int(of_yield), #
                                                                           rsfof=rsfof)
                tmp_file = open('datacards/datacards_Smart_NLL/%s.txt'%(bin_name),'w')

                tmp_file.write(dc)
                tmp_file.close()
            print nll
            output = '''Drell-Yan like &    {lowDY:.1f} $\\pm$ {lowDYe:.1f}    & {onzDY:.1f} $\\pm$ {onzDYe:.1f}    & {highDY:.1f} $\\pm$ {highDYe:.1f} \\\\ 
Flavour symmetric        &    {lowFS:.1f} $\\pm$ {lowFSe:.1f}   & {onzFS:.1f} $\\pm$ {onzFSe:.1f}   & {highFS:.1f} $\\pm$ {highFSe:.1f} \\\\ \\midrule
Total background         &    {lowTot:.1f} $\\pm$  {lowTote:.1f}  & {onzTot:.1f} $\\pm$ {onzTote:.1f}  & {highTot:.1f} $\\pm$ {highTote:.1f} \\\\  \\midrule
Observed                 &    {lowObs}     & {onzObs}     & {highObs}  \\\\ \n'''.format(lowDY=dy[0],        onzDY=dy[1],       highDY=dy[2],
                                                                          lowFS=Fs[0],        onzFS=Fs[1],        highFS=Fs[2],
                                                                          lowTot=dy[0]+Fs[0], onzTot=dy[1]+Fs[1], highTot=dy[2]+Fs[2],
                                                                          lowObs=int(Obs[0]), onzObs=int(Obs[1]), highObs=int(Obs[2]),
                                                                          lowDYe=dyE[0] ,     onzDYe=dyE[1],      highDYe =dyE[2],
                                                                          lowFSe=FsE[0] ,     onzFSe =FsE[1],     highFSe =FsE[2],
                                                                          lowTote=math.sqrt(dyE[0]**2 + FsE[0]**2),
                                                                          onzTote=math.sqrt(dyE[1]**2 + FsE[1]**2),
                                                                          highTote=math.sqrt(dyE[2]**2 + FsE[2]**2))
            print output     

def getL2Norm(histos):
    norm  = 0
    histo = histos[0].Clone()
    histo.Add(histos[1],-1.)
    for i in range(histo.GetNbinsX()+1):
        norm += histo.GetBinContent(i)*histo.GetBinContent(i)
    norm *= (histo.GetXaxis().GetXmax() - histo.GetXaxis().GetXmin()) / histo.GetBinWidth(1)
    norm = math.sqrt(norm)
    norm *= 100/histos[0].Integral()
    return norm

 
if __name__ == "__main__":


    parser = optparse.OptionParser(usage="usage: %prog [opts] FilenameWithSamples", version="%prog 1.0")
    parser.add_option('-s', '--samples', action='store', type=str, dest='sampleFile', default='samples.dat', help='the samples file. default \'samples.dat\'')
    parser.add_option('-b', '--batch'  , action='store', type=int, dest='batchMode' , default=1            , help='set batch mode. default true')
    parser.add_option('-d', '--dody'   , action='store', type=int, dest='doDY'      , default=0            , help='do dy sample as well. default is false')
    parser.add_option('-g', '--dosig'   , action='store', type=int, dest='doSig'     , default=1            , help='do signal sample as well. default is true')
    parser.add_option('-q', '--quick'  , action='store', type=int, dest='quick'     , default=0            , help='quick test, not full thing. default false')
    parser.add_option('-i', '--ingredients', action='store', type=str, dest='ingredientFile', default='ingredients.dat', help='the ingredients file. default \'ingredients.dat\'')

    (opts, args) = parser.parse_args()


    if opts.batchMode:
        print 'i am setting this stupid batch mode now'
        gROOT.SetBatch(True)
    gROOT.ProcessLine('.L include/tdrstyle.C')
    r.setTDRStyle() 

    print 'running with these options'
    for key, value in opts.__dict__.items():
        print '%-20s : %-20s' %(key, value)

    global lumi
    lumi = 2.3

    cuts = CutManager.CutManager()


    if opts.quick:
        rfile = r.TFile('/afs/cern.ch/work/m/mdunser/public/pdfsForLikelihood/pdfs_version22_savingTheWorkspace.root', 'READ')
        pdfmet = rfile.Get('em_data_fit_histo_met_ds_cuts_of_sr_met150');pdfmet.Scale(1./pdfmet.Integral())
        pdfmlb = rfile.Get('em_data_fit_histo_mlb_ds_cuts_of_sr_met150');pdfmlb.Scale(1./pdfmlb.Integral())
        pdfzpt = rfile.Get('em_data_fit_histo_zpt_ds_cuts_of_sr_met150');pdfzpt.Scale(1./pdfzpt.Integral())
        pdfldr = rfile.Get('em_data_fit_histo_ldr_ds_cuts_of_sr_met150');pdfldr.Scale(1./pdfldr.Integral())

        pdfnll = r.TH1F('foo', 'bar', 200, 15, 35)
        for i in range(250000):
            lh  = pdfmet.GetBinContent(pdfmet.FindBin(pdfmet.GetRandom()))
            lh *= pdfmlb.GetBinContent(pdfmlb.FindBin(pdfmlb.GetRandom()))
            lh *= pdfzpt.GetBinContent(pdfzpt.FindBin(pdfzpt.GetRandom()))
            lh *= pdfldr.GetBinContent(pdfldr.FindBin(pdfldr.GetRandom()))
            pdfnll.Fill(-1.*math.log(lh))
        pdfnll.Scale(4./pdfnll.Integral())

            
        
        ttDatasets = ['TTLep_pow']
        treeTT = Sample.Tree(helper.selectSamples(opts.sampleFile, ttDatasets, 'TT'), 'TT', 0, isScan = 0)
        ttDatasets = ['MuEG']
        treeDA = Sample.Tree(helper.selectSamples(opts.sampleFile, ttDatasets, 'DA'), 'DA', 1, isScan = 0)
        siDatasets = ['SMS_T6bbllslepton_mSbottom-600To900_mLSP-200To800', 
                      'SMS_T6bbllslepton_mSbottom-400To550_mLSP-200To500' ]
        treeSI = Sample.Tree(helper.selectSamples(opts.sampleFile, siDatasets, 'SI'), 'SI', 0, isScan = 1)
        c_sr_sf    = cuts.AddList([cuts.METJetsSignalRegionMET150 , cuts.GoodLeptonSF()])
        c_sr_of    = cuts.AddList([cuts.METJetsSignalRegionMET150 , cuts.GoodLeptonOF()])
        var = 'ldr'; v = vParams(var)
        h_tt_sr_sf  = treeTT.getTH1F(1., 'tt_%s_sr_sf' %var, v.vtree, v.nbins, v.xmin, v.xmax , c_sr_sf   , '', v.title, ofBin = True); h_tt_sr_sf.Scale(1./h_tt_sr_sf.Integral())
        h_tt_sr_of  = treeTT.getTH1F(1., 'tt_%s_sr_of' %var, v.vtree, v.nbins, v.xmin, v.xmax , c_sr_of   , '', v.title, ofBin = True); h_tt_sr_of.Scale(1./h_tt_sr_of.Integral())
        h_da_sr_of  = treeDA.getTH1F(1., 'da_%s_sr_of' %var, v.vtree, v.nbins, v.xmin, v.xmax , '('+c_sr_of +' met_pt < 150)'  , '', v.title, ofBin = False); h_da_sr_of.Scale(1./h_da_sr_of.Integral())

        pdfnll.SetLineColor(r.kRed)
        h_tt_sr_sf.SetMarkerColor(r.kBlue)
        h_tt_sr_of.SetMarkerColor(r.kGreen)
        h_da_sr_of.SetMarkerColor(r.kBlack)

        pdfnll.GetYaxis().SetRangeUser(0., 0.1)
        pdfnll.Draw('hist norm')
        h_tt_sr_sf.Draw('same norm pe')
        h_tt_sr_of.Draw('same norm pe')
        h_da_sr_of.Draw('same norm pe')

        # quick mt2 staudies
        var = 'mt2'; v = vParams(var)
        h_tt_mt2  = treeTT.getTH1F(1., 'tt_%s_sr_sf' %var, v.vtree, v.nbins, v.xmin, v.xmax , c_sr_sf  , '', v.title, ofBin = False); h_tt_sr_sf.Scale(1./h_tt_sr_sf.Integral())
        
        sigcuts = cuts.AddList([c_sr_sf.replace(cuts.twoLeptons, 't.nPairLep_Edge > 0'), 'GenSusyMScan1_Edge == %.0f && GenSusyMScan2_Edge == %.0f'%(900, 200)])
        h_si_mt2_900_200  = treeSI.getTH1F(1., 'da_%s_sr_of' %var, v.vtree, v.nbins, v.xmin, v.xmax , sigcuts  , '', v.title, ofBin = False); h_si_sr_of.Scale(1./h_si_sr_of.Integral())
        sigcuts = cuts.AddList([c_sr_sf.replace(cuts.twoLeptons, 't.nPairLep_Edge > 0'), 'GenSusyMScan1_Edge == %.0f && GenSusyMScan2_Edge == %.0f'%(400, 375)])
        h_si_mt2_400_375  = treeSI.getTH1F(1., 'da_%s_sr_of' %var, v.vtree, v.nbins, v.xmin, v.xmax , sigcuts  , '', v.title, ofBin = False); h_si_sr_of.Scale(1./h_si_sr_of.Integral())

        h_tt_mt2.SetMarkerColor(r.kBlue)
        h_si_mt2.SetMarkerColor(r.kGreen)
        h_si_mt2.SetMarkerColor(r.kRed-2)
        


    

    else:
#/afs/cern.ch/work/m/mdunser/public/pdfsForLikelihood/pdfs_version22_savingTheWorkspace.root
        a = nllObject('foobar','pdfs/pdfs_version24_savingTheWorkspace_june2016.root', opts.doDY, opts.doSig)
#        a.CheckCorrelations()
        a.kernelAnalyticComparison(1,1)
#        a.getMT2NllComparison()
#        a.plotZMassCuts()
#        a.profile_nll_zmass()
#        a.kernelAnalyticComparison()
#        a.dataCard()

#        a.getROC()
#        a.getChi2ForFits()
        print 'doing vars after cut'
#        a.varsAfterCut()
#        a.compareOFSF()
#        a.getDYEff()
