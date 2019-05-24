import ROOT as r
reffile = r.TFile.Open('~/Downloads/btageff__ttbar_powheg_pythia8_25ns_Moriond17.root')

refhist = reffile.Get('h2_BTaggingEff_csv_loose_Eff_b').Clone('ref')
refhist.Reset()

stuff={
    '2016': { 'files' : ['TTTo2L2Nu_2016_part1','TTTo2L2Nu_2016_part2'],
              'discs'   : {
                  'btagDeepB'     : [ 0.2217,  0.6321,  0.8953 ],
                  #'btagDeepFlavB' : [ 0.0614,  0.3093,  0.7221 ],
                  # https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation2016Legacy
              }
          },

    '2017': { 'files' : ['TTTo2L2Nu_2017_part1','TTTo2L2Nu_2017_part2'],
              'discs'   : {
                  #'btagCSVV2'     : [0.5803,  0.8838,  0.9693 ],
                  'btagDeepB'     : [0.1522,  0.4941,  0.8001 ],
                  #'btagDeepFlavB' : [0.0521,  0.3033,  0.7489 ],
                  # https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation94X
              }
          },

    '2018': { 'files' : ['TTTo2L2Nu_2018_part1','TTTo2L2Nu_2018_part2'],
              'discs'   : {
                  'btagDeepB'     : [0.1241,  0.4184,  0.7527 ],
                  #'btagDeepFlavB' : [0.0494,  0.2770,  0.7264 ],
                  # https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation102X
              },
          },
}

flavs = {
    'b'   : 5,
    'c'   : 4,
    'udsg': 0,
}

path = "/pool/ciencias/userstorage/sscruz/NanoAOD/EdgeBtagWeightsv2_merged/"

tfile_out = r.TFile.Open('btagEffs_Nanov4.root','recreate')
for year in stuff: 
    Events = r.TChain('Events')
    for fil in stuff[year]['files']: Events.Add('{path}/{fil}/nanoAODskim/Events.root'.format(path=path,fil=fil)) 
    for flav in flavs:
        den = refhist.Clone('den_%s'%year)
        Events.Project('den_%s'%year,'abs(JetSel_Edge_eta):JetSel_Edge_pt','JetSel_Edge_hadronFlavour == %d'%flavs[flav])
        for disc in stuff[year]['discs']:
            for wp,wpname in zip(stuff[year]['discs'][disc],['loose','med','tight']):
                name = 'h2_Eff_%s_%s_%s_%s'%(year,disc,wpname,flav)
                num = refhist.Clone(name)
                Events.Project(name, 'abs(JetSel_Edge_eta):JetSel_Edge_pt','JetSel_Edge_hadronFlavour == %d && JetSel_Edge_%s > %f'%(flavs[flav],disc,wp))
                num.Divide(den)
                print "b-tagging effiency, MC, {disc} {wpname}, {flav}-jets, |#eta| vs p_{{T}}, {year}".format(disc=disc,wpname=wpname,year=year,flav=flav)
                num.SetTitle("b-tagging effiency, MC, {disc} {wpname}, {flav}-jets, |#eta| vs p_{{T}}, {year}".format(disc=disc,wpname=wpname,year=year,flav=flav))
                num.Write()
