import sys

if len(sys.argv)==2:
    fn_1 = sys.argv[1]
f2=open('rmsd'+fn_1,'w')
for i in open(fn_1,'r'):
    if 'ATOM' in i:
        if int(i[22:26]) != 1:
            f2.write(i)
f2.close()
