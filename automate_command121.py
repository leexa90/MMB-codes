import sys
sys.path.append ('/cluster/apps/x86_64/packages/lib/python2.4/site-packages')
import numpy as np


'''
automatically makes MMB command121 file

'''

KK=999
sysargv=20
if len(sys.argv) > 1:
    KK =  int(sys.argv[2])
    sysargv =int(sys.argv[1])
#pi_backbone = ['O3',
dict_atom = {" O4'": 77, " O2'": 77,  " C2'": 77, 
 " C4'": 77, " O5'": 77, " O3'": 77, 
 " C1'": 77, " C3'": 77, ' P  ': 77, " C5'": 77}
atom_keys = [" C1'", " C2'", " C3'", " C4'", " C5'", 
             " O2'", " O3'", " O4'", " O5'", ' P  ']

atom_keys = [' P  ']

atom_keys = [" C1'", " C2'", " C3'", " C4'", " C5'", 
             " O2'", " O3'", " O4'", " O5'", ' P  ']

dict_atom = { 'CYT':  [" N1 ", " C2 ", " C5 "],
              'URA':  [" N1 ", " C2 ", " C5 "],
              'ADE':  [" C4 ", " N9 ", " N7 "],
              'GUA':  [" C4 ", " N9 ", " N7 "] }

dict_atom = { 'CYT':  ['   P'," C4'","  N1", "  C2", "  C5"],
              'URA':  ['   P'," C4'","  N1", "  C2", "  C5"],
              'ADE':  ['   P'," C4'","  C4", "  N9", "  N7"],
              'GUA':  ['   P'," C4'","  C4", "  N9", "  N7"] }


#dict_atom = { 'CYT':  ["  C4", "  C2", "  C6"],
#              'URA':  ["  C4", "  C2", "  C6"],
#              'ADE':  ["  C4", "  C2", "  C6"],
#              'GUA':  ["  C4", "  C2", "  C6"] }




data = []
COUNTER =0
##for i in open('ZZZ.pos_out', 'r'):
##    if 'ATOM' in i and i[12:16] in dict_atom[i[17:20]] and int(i[22:26])%1==0:
##        xyz = [i[30:38],i[38:46],i[46:54]]
##        xyz = [float(j) for  j in xyz]
##        temp = [i[6:11],i[12:16],np.array(xyz),int(i[22:26])]
##        data += [temp,]
##        COUNTER += 1

data_model= []
COUNTER =0
##for i in open('rna.pos_out', 'r'):
##    if 'ATOM' in i and i[12:16] in dict_atom[i[17:20]] and int(i[22:26])%1==0:
##        xyz = [i[30:38],i[38:46],i[46:54]]
##        xyz = [float(j) for  j in xyz]
##        temp = [i[6:11],i[12:16],np.array(xyz),int(i[22:26])]
##        data_model += [temp,]
##        COUNTER += 1
        
[' 1637', " C5'", np.array([  8.737,  -9.379, -74.458])]
pair_wise_constrains = []
resi_num = [x[3] for x in data]

D_resnum = {}
for i in range(0,len(data)):
    resNum=data[i][3]
    if resNum not in D_resnum.keys():
        D_resnum[resNum] = [data[i],]
    else:
        D_resnum[resNum] += [data[i],]

result = []
f3=open('alignment.input','r')
for i in f3:
    if '>' not in i and i != '\n':
        result += [i,]
a=result[0]
b=result[1]

A,B='',''
counterA,counterB=0,0
Resi_map,template_resi ={},[]
for i in range(len(b)):
    if a[i].lower() in ['a','u','g','c'] and b[i].lower() in ['a','u','g','c']:
        Resi_map[counterB+1]=counterA+1
        template_resi += [counterB+1]
    if a[i].lower() in ['a','u','g','c'] :
        counterA+=1
        A = A + a[i].lower() 
    if b[i].lower() in ['a','u','g','c'] :
        counterB+=1
        B = B + b[i].lower()  

print 'firstStage 121\nlastStage 121'
print 'reportingInterval .5\nnumReportingIntervals 100'
print 'RNA A 1 ', B.upper()
print 'RNA C 1 ', A.upper()
print 'mobilizer Rigid A 1 %s ' %str(len(B))

counter = 0
for i in open("PA.secstruct",'r'):
    if counter ==0:
        seq = i[0:len(i)-1]
    elif counter == 1:
        ss = i[0:len(i)-1]
    counter  += 1

num=[]
for i in range(len(ss)):
    if ss[i] == '.':
        if num != []:
            if num==range(num[0],num[-1]+1):
                None#print 'mobilizer Rigid C' ,num[0],num[-1] 
            num=[]
    else:
        num += [i+1,]
    


def recursive_get_bp(ss,result=[]):
    if '(' not in ss and ')' not in ss:
        return result
    for i in range(len(ss)):
        if ss[i] == '(':
            for j in range(1+i,len(ss)):
                if ss[j] == '(':
                    break
                if ss[j] == ')':
                    #print i,j
                    ss = ss[0:i]+'.'+ss[i+1:j]+'.'+ss[j+1:]
                    return recursive_get_bp(ss,result+[[i+1,j+1],])
base_paired_target=recursive_get_bp(ss)
print 'alignmentForces forceConstant 300.0'
inside_template=[]
for i in Resi_map:
    inside_template += [Resi_map[i],]
    print ('alignmentForces A %s %s C %s %s' %(str(i) ,str(i) ,str(Resi_map[i]),str(Resi_map[i])))
print 'baseInteractionScaleFactor 1'
for j in base_paired_target:
    #if  j[0] not in inside_template and j[1] not in inside_template:
        print('baseInteraction C %s WatsonCrick C %s WatsonCrick Cis' %(str(j[0]),str(j[1])))
##for i in range(inside_template[0],inside_template[-1]+1):
##    if i not in inside_template:
##        print 'includeResidues C %s %s'  %(i,i)
##
##for i in range(len(ss)):
##    resnum=i+1
    



