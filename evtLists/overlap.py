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

ofile = open('synchFullDataset_onZ_edgeSR_SF.txt', 'wr')
a = set()
for i in f_list:
    f = open(i,'r')
    for ev in list(f):
        a.add(ev.replace('*',':').rstrip('\r')[10:])
    f.close()

print len(a)
for i in a:
    ofile.writelines(i)
ofile.close()
