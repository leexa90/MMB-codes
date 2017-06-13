import sys
file_n = 'topMMB.pdb'
if len(sys.argv) > 1:
    file_n = sys.argv[1]

file_n_w='new'+file_n
dictt = {'G':'GUA','A':'ADE','U':'URA','C':'CYT'}


#dictt = {'G':'RG','A':'RA','U':'RU','C':'RC'}
#dictt = {'GUA':'GUA','ADE':'ADE','URA':'URA','CYT':'CYT'}
f2=open(file_n_w,'w')
for i in open(file_n,'r'):
    if 'ATOM' in i:
        f2.write(i[0:17]+dictt[i[17:20].strip()]+i[20:])
f2.write('TER\n')
f2.close()