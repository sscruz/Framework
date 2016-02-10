import os, math



f_list = ['synchOnFullDataset_0.txt',
          'synchOnFullDataset_1.txt',
          'synchOnFullDataset_2.txt',
          'synchOnFullDataset_3.txt',
          'synchOnFullDataset_4.txt',
          'synchOnFullDataset_5.txt',
          'synchOnFullDataset_6.txt',
          'synchOnFullDataset_7.txt',
          'synchOnFullDataset_8.txt']
               
ofile = open('synch_rt.txt', 'wr')
a = set()
b = set()
for i in f_list:
    f = open(i,'r')
    for ev in list(f):
        if (ev.find("lumi") != -1 or ev.find("***") != -1 or ev.find("==") != -1):
            continue  
        reducedString = ':'.join(ev.replace('*',':').split(':')[2:5])+'\n'
        if not reducedString in a:
            a.add(reducedString)
            b.add(ev.replace('*',':')[12:])
    f.close()

print len(a)
print len(b)
for i in b:
    ofile.writelines(i)
ofile.close()
