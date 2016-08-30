import os, math



f_list = ['edgeSR_met150_2j_0.txt',
          #'edgeSR_met150_2j_1.txt',
          'edgeSR_met150_2j_2.txt',
          #'edgeSR_met150_2j_3.txt',
          'edgeSR_met150_2j_4.txt']
          #'edgeSR_met150_2j_5.txt']
               
ofile = open(f_list[0].split('_')[0]+'.txt', 'wr')
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
