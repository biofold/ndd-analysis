import sys


def get_file(filename):
    d={}
    f=open(filename)
    for line in f:
        if line.find('Gene_set')>-1: continue
        v=line.split('\t')
        sg=v[1].split('(GO:')
        go='GO:'+sg[1][:-1]
        name=sg[0].strip().replace(' ','_')
        n=v[2].split('/')
        d[go]=[name,n[0],n[1],v[4]]
    return d



# Aggregate 3 file output cancer_setx, noncancer_setx and setx
# Calculate fisher using the columns 3,4 5,6
# adjust the pvalue in column 17 (0-based)

if __name__ == '__main__':
    f1=sys.argv[1]    # cancer_setx
    f2=sys.argv[2]    # noncancer_setx
    f3=sys.argv[3]    # setx
    n1=sys.argv[4]    # number cancer_setx
    n2=sys.argv[5]    # number noncancer_setx
    n3=sys.argv[6]    # number setx
    d1=get_file(f1)
    d2=get_file(f2)
    d3=get_file(f3)
    ng={}
    for go in d1.keys():
        ng[go]=d1[go][0]
    for go in d2.keys():
        ng[go]=d2[go][0]
    for go in d3.keys():
        ng[go]=d3[go][0]  
    for k in ng.keys():
        v1=d1.get(k,[ng[k],'0','0','1.0'])[1:]
        r1=str(int(n1)-int(v1[0]))
        v2=d2.get(k,[ng[k],'0','0','1.0'])[1:]
        r2=str(int(n2)-int(v2[0]))
        v3=d3.get(k,[ng[k],'0','0','1.0'])[1:]
        r3=str(int(n3)-int(v3[0]))
        v=[ng[k],'('+k+')',v1[0],r1,v2[0],r2,v3[0],r3,v1[-1],v2[-1],v3[-1]]        
        print ('\t'.join(v))
        

