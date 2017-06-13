


import sys

if len(sys.argv)==3:
    fn_1 = sys.argv[1]
    A = sys.argv[2]

f2=open(fn_1+'.pdb','w')
for i in open(fn_1,'r'):
    if 'ATOM' in i:
        if i[21]==A:
            f2.write(i)
