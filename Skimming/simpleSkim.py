import include.helper as helper
import ROOT as r 
import os

samples = [#'DYJetsToLL_M50_HT100to200', 
           #'DYJetsToLL_M50_HT1200to2500', 
           #'DYJetsToLL_M50_HT200to400', 
           #'DYJetsToLL_M50_HT400to600', 
           #'DYJetsToLL_M50_HT600to800', 
           #'GluGluToContinToZZTo2e2nu', 
           #'GluGluToContinToZZTo2mu2nu', 
           #'GluGluToContinToZZTo2mu2tau', 
           #'TBar_tch_powheg', 
           #'TTHnobb_pow', 
           #'TTJets',
           #'SMS_TChiWZ',
           #'WZTo2L2Q'

    #'TTLLJets_m1to10']
    #'TTWToLNu']
    'TTWW', 
    'TTWZ', 
    'TTW_LO',
    #'TTZH'
    #       'TTZToLLNuNu']
    # 'TTZZ']
    'T_tWch_noFullHad', 
           #  'T_tch_powheg']
    'WWW_4F', 
    'WZG', 
    'WZTo3LNu_amcatnlo', 
    'WZZ', 
    'ZZTo2L2Nu', 
    'ZZZ']
#'tZq_ll']

#samples = ['SMS_TChiWZ_ZToLL']


selection = 'nFatJetSel_Edge > 0 && (nPairLep_Edge > 0&&Lep1_pt_Edge > 25 && Lep2_pt_Edge > 20.&&lepsDR_Edge > 0.1&&lepsMll_Edge > 20)'

out='kk'

for samp in samples:
    fil = helper.selectSamples('samples.dat', [samp], 'MC')
    subs = open(fil,'r').readlines()
    for sub in subs: 
        print 'Processing', sub
        for i in range(20):
            sub = sub.replace('  ',' ')
            
        tfile = r.TFile.Open(sub.split(' ')[4] + '/' + 'evVarFriend_{k}.root'.format(k=sub.split(' ')[2]))
        ft = tfile.Get('sf/t')
        print 'Friend tree is', ft
        isSMS = bool(tfile.Get('CountSMS'))
        outfile = r.TFile.Open(out + '/evVarFriend_{k}.root'.format(k=sub.split(' ')[2]),'recreate')
        outfile.mkdir('sf').cd()
        nu = ft.CopyTree(selection) 
        nu.AutoSave()
        

        outfile.Close()
        os.system('rootcp {inp}:SumGenWeights {out}:SumGenWeights'.format(inp=(sub.split(' ')[4] + '/' + 'evVarFriend_{k}.root'.format(k=sub.split(' ')[2])),
                                                                          out=  out + '/evVarFriend_{k}.root'.format(k=sub.split(' ')[2])))
        
        if isSMS:
            os.system('rootcp {inp}:CountSMS {out}:CountSMS'.format(inp=(sub.split(' ')[4] + '/' + 'evVarFriend_{k}.root'.format(k=sub.split(' ')[2])),
                                                                    out=  out + '/evVarFriend_{k}.root'.format(k=sub.split(' ')[2])))
