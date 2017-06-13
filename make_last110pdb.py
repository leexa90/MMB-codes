import sys,os

if len(sys.argv)==5:
    f1_n=sys.argv[1]
    A=sys.argv[2]
    f2_n=sys.argv[3]
    B=sys.argv[4]
print A,B
f3=open('last.120.pdb','w')

f1=open(f1_n,'r')
f2=open(f2_n,'r')

for i in f1:
    if 'ATOM' in i:
        f3.write(i[0:21]+A+i[22:])
f3.write('TER\n')

for i in f2:
    if 'ATOM' in i:
        f3.write(i[0:21]+B+i[22:])

f1.close()
f2.close()
f3.close()
