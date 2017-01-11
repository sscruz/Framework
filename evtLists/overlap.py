import os, math, sys



f_all = sys.argv[-1] #'test0.txt'
               
ofile = open(f_all.split('.')[0]+'_reduced.txt', 'wr')
a = set()
b = set()
f = open(f_all,'r')
for ev in list(f):
    if (ev.find("***") != -1 or ev.find("==") != -1):
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
