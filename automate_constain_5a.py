import sys
sys.path.append ('/cluster/apps/x86_64/packages/lib/python2.4/site-packages')
import numpy as np

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
                print 'mobilizer Rigid C' ,num[0],num[-1] 
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
            
inside_template=[]
for i in Resi_map:
    inside_template += [Resi_map[i],]
    print ('alignmentForces A %s %s C %s %s' %(str(i) ,str(i) ,str(Resi_map[i]),str(Resi_map[i])))
print 'baseInteractionScaleFactor 1'
for j in base_paired_target:
    #if  j[0] not in inside_template and j[1] not in inside_template:
        print('baseInteraction C %s WatsonCrick C %s WatsonCrick Cis' %(str(j[0]),str(j[1])))
for i in range(inside_template[0],inside_template[-1]+1):
    if i not in inside_template:
        print 'includeResidues C %s %s'  %(i,i)

for i in range(len(ss)):
    resnum=i+1
    
chk =[]
result = []
def simillar_atoms(i,j):
    if i[1] == j[1]:
        return True
    if ((j[1] == " N1 " or j[1] == " N9 ") and
        (i[1] == " N1 " or i[1] == " N9 ")) : #takes care of 
        return True
    if ((j[1] == " C2 " or j[1] == " C4 ") and
        (i[1] == " C2 " or i[1] == " C4 ")) :
        return True
    if ((j[1] == " C5 " or j[1] == " N7 ") and
        (i[1] == " C5 " or i[1] == " N7 ")) :
        return True
    if ((j[1] == "  N1" or j[1] == "  N9") and
        (i[1] == "  N1" or i[1] == "  N9")) :
        return True
    if ((j[1] == "  C2" or j[1] == "  C4") and
        (i[1] == "  C2" or i[1] == "  C4")) :
        return True
    if ((j[1] == "  C5" or j[1] == "  N7") and
        (i[1] == "  C5" or i[1] == "  N7")) :
        return True
    

def got_close_by(data_i,data_j,sysargv=sysargv):
    for i in data_i:        
        for j in data_j :
            if simillar_atoms(i,j) is True :
                vec = i[2]-j[2]
                vec = round(sum(vec**2)**.5,2)
                #Res_number of template_data == Resi_map[model_resi_number] (resi_map is a mapping dictionary)
                model_i =  [x for x in data_model if x[3]==Resi_map[i[3]]]  
                model_j =  [x for x in data_model if x[3]==Resi_map[j[3]]]
                # only accept simillar atoms
                model_i =  [x for x in model_i if simillar_atoms(i,x)==True]
                model_j =  [x for x in model_j if simillar_atoms(j,x)==True]
                if len(model_i) ==1 and len(model_j) ==1:
                    i,j= model_i[0],model_j[0]
                    if vec > 3 and vec < sysargv:
                        return True
           
list_of_intr_residues = []
list_of_intr_residues_dd = []
list_of_intr_residues_sssd = []
## Iterate all residue numbers of template,from template find residues close together ###
## Close is defined as having 1 atom < sysargv### 
for resnum in sorted(D_resnum.keys())[0:-1]:
    for resnum2 in range(resnum+1, 1+len(D_resnum.keys())):
        if resnum in Resi_map.keys() and resnum2 in Resi_map.keys():
            if got_close_by(D_resnum[resnum],D_resnum[resnum2],sysargv) is True:
                list_of_intr_residues += [[resnum,resnum2],]
                list_of_intr_residues_dd += [[resnum,resnum2],]

### Find equailvalent from model ### 
for i in list_of_intr_residues:
    resnum,resnum2 = i[0],i[1]
    if resnum in template_resi and resnum2 in template_resi:
        for i in D_resnum[resnum]:        
            for j in D_resnum[resnum2] :
                if simillar_atoms(i,j) is True :
                    vec = i[2]-j[2]
                    vec = round(sum(vec**2)**.5,2)
                    ### see if Model contains the residue ### 
                    model_i =  [x for x in data_model if x[3]==Resi_map[i[3]]]
                    model_j =  [x for x in data_model if x[3]==Resi_map[j[3]]]
                    model_i =  [x for x in model_i if simillar_atoms(i,x)==True]
                    model_j =  [x for x in model_j if simillar_atoms(j,x)==True]
                    if len(model_i) ==1 and len(model_j) ==1:
                        i,j= model_i[0],model_j[0]
                        sort = sorted([int(i[0]),int(j[0])])
                        pair_wise_constrains += [[sort[0],sort[1],vec],]                    

#pair_wise_constrains = pair_wise_constrains[0:5] + pair_wise_constrains[0:5] 
for i in range(0,len(pair_wise_constrains)):
    if pair_wise_constrains[i][0:2] not in chk:
        chk += [pair_wise_constrains[i][0:2],]
        result += [pair_wise_constrains[i],]
    elif pair_wise_constrains[i][0:2] in chk:
        None
        print 'double'
    
#print len(result)
for i in sorted(result):
    None
    print ('~constraint['+"\\"+'atom1{'+str(i[0])+'}'
           +"\\"+'atom2{'+str(i[1])+'}'+"\\"+'fk{'+str(KK)
           +'}'+"\\"+'r0{'+str(i[2])+'}]')


