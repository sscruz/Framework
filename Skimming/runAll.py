import os
from multiprocessing import Pool



################################################################################################################
### Please make sure that you are skimming what you want in skimmer.C, there are two skimming modes:         ###
### - For general edge analysis -> 2 leptons and 2 or more jets                                              ###   
### - For trigger calculation -> 2 leptons and HT > 200                                                      ###
################################################################################################################

def runDataset(ins):
    os.system('root -l -b -q runSkim.C+\(\\\"{sample}\\\",\\\"{path}\\\"\);'.format(sample=ins[0],path=ins[1]))

#pathList = ['/tmp/mvesterb/']
#pathList = ['/eos/cms/store/user/mvesterb/edgeFriends_Aug11_Unskimmed/skimThis/']
#pathList = ['/afs/cern.ch/work/p/pablom/public/Edge-Production-Friend-Trees-June-2017-Feb23-ReReco/']
pathList = ['/afs/cern.ch/work/m/mvesterb/public/MC_samples/edgeZ/MC_Aug11/make/']
#pathList = ['/afs/cern.ch/work/m/mvesterb/public/MC_samples/edgeZ/data/']

tasks = []
for path in pathList:
    os.system('root -l -b -q runSkim.C+')
    files = os.listdir(path)
    for fil in files:
        if not 'evVarFriend_' in fil: continue
        dataset = fil.replace('evVarFriend_','').replace('.root','')
        tasks.append([dataset,path])
pool = Pool(8)
pool.map(runDataset, tasks)
