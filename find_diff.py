f1= open('rna_tleap_outtopZZZ.pos_out.pdb2','r')
f2=open('rna_tleap_outtopMMB.pdb2','r')
a={}
for i in f1:
    if 'ATOM' in i:
        if int(i[22:26].strip()) not in a:
            a[int(i[22:26].strip())] = [i,]
        else:
            a[int(i[22:26].strip())] += [i,]

b={}
for i in f2:
    if 'ATOM' in i:
        if int(i[22:26].strip()) not in b:
            b[int(i[22:26].strip())] = [i,]
        else:
            b[int(i[22:26].strip())] += [i,]

for i in range(1,127):
    for j in b[i]:
        inside=0
        for k in a[i]:
            if j[12:16].strip()==k[12:16].strip():
                inside=1
        if inside==0:
            print j,
        

for i in range(1,127):
    for j in a[i]:
        inside=0
        for k in b[i]:
            if j[12:16].strip()==k[12:16].strip():
                inside=1
        if inside==0:
            print j,
            
