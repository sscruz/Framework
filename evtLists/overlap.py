import os, math, sys

def isNumber(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

def isInteger(s):
    try:
        int(s)
        return True
    except ValueError:
        return False

def convertString(string):
    string = string.split('*')[2:]
    newstring = []
    for sub in string:
        sub = sub.replace(' ','')
        if not isNumber(sub):
            newstring.append('{s:10s}'.format(s=sub.replace('_Edge','').replace('_Ed','').replace('_Edg','')) )
        elif isInteger(sub):
            newstring.append('{val:10d}'.format(val=int(sub)) )
        else:
            newstring.append('{val:10.3f}'.format(val=float(sub)) )
    newstring = ' : '.join(newstring)
    return newstring


f_all = sys.argv[-1] #'test0.txt'
               
ofile = open(f_all.split('.')[0]+'_reduced.txt', 'wr')
a = set()
b = set()
f = open(f_all,'r')
for ev in list(f):
    if (ev.find("***") != -1 or ev.find("==") != -1):
        continue  
    reducedString = ':'.join(ev.replace('*',':').split(':')[2:5])+'\n'
    newstring = convertString(ev)
    if not reducedString in a:
        a.add(reducedString)
        #b.add(ev.replace('*',':')[12:])
        b.add(newstring)
f.close()

print len(a)
print len(b)
for i in b:
    ofile.writelines(i)
ofile.close()
