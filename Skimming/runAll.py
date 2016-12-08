import os
from multiprocessing import Pool

def runDataset(ins):
    os.system('root -l -b -q runSkim.C+\(\\\"{sample}\\\",\\\"{path}\\\"\);'.format(sample=ins[0],path=ins[1]))

pathList = ['/mnt/t3nfs01/data01/shome/pablom/trees-Nov-28/']

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
