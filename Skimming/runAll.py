import os
from multiprocessing import Pool



################################################################################################################
### Please make sure that you are skimming what you want in skimmer.C, there are two skimming modes:         ###
### - For general edge analysis -> 2 leptons and 2 or more jets                                              ###   
### - For trigger calculation -> 2 leptons and HT > 200                                                      ###
################################################################################################################

def runDataset(ins):
    os.system('root -l -b -q runSkim.C+\(\\\"{sample}\\\",\\\"{path}\\\"\);'.format(sample=ins[0],path=ins[1]))

pathList = ['/afs/cern.ch/work/s/sesanche/public/forEdge/nTuplesForMoriond/jan27/']

tasks = []
for path in pathList:
    os.system('root -l -b -q runSkim.C+')
    files = os.listdir(path)
    for fil in files:
        if not 'evVarFriend_' in fil: continue
        dataset = fil.replace('evVarFriend_','').replace('.root','')
        tasks.append([dataset,path])
pool = Pool(4)
pool.map(runDataset, tasks)
