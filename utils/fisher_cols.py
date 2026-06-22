import sys, numpy
from scipy import stats

def akashi(mat):
        n=mat[0][0]+mat[0][1]+mat[1][0]+mat[1][1]
        p=(mat[0][0]+mat[0][1])/float(n)
        q=(mat[0][0]+mat[1][0])/float(n)
        E=n*p*q
        V=p*(1-p)*q*(1-q)*(n**2)/float(n-1)
        Z=(mat[0][0]-E)/numpy.sqrt(V)
        return Z


def all_fisher(filename,vpos1,vpos2,c=0):
        f=open(filename,'r')
        for line in f:
                v=line.rstrip().split()
                ct=numpy.array([[float(v[vpos1[0]])+c,float(v[vpos1[1]])+c],[float(v[vpos2[0]])+c,float(v[vpos2[1]])+c]])
                print (line.rstrip()+'\t'+cal_fisher(ct))
        return


def cal_fisher(ct):
        Z=akashi(ct)
        oddsr,pvalueg=stats.fisher_exact(ct,alternative='greater')
        oddsr,pvaluel=stats.fisher_exact(ct,alternative='less')
        oddsr,pvalue2=stats.fisher_exact(ct)
        if pvaluel<pvalueg:
                pvalue=pvaluel
                side='L'
        else:
                pvalue=pvalueg
                side='G'
        pvals=" p-valuel: %10.2e p-valueg: %10.2e" %(pvaluel,pvalueg)
        #return "OR: %8.3f" %oddsr+' '+pvals+' '+side+" p-value: %10.2e" %pvalue+" p-value2: %10.2e" %pvalue2+' Z: %7.2f' %Z
        
        return '\t'.join(str(i) for i in ['OR:',oddsr,'p-valuel:',pvaluel,'p-valueg:',pvalueg,'p-value:',pvalue,'p-value2',pvalue2,'Z:',Z])


if __name__=="__main__":
        filename=sys.argv[1]
        pos1=sys.argv[2].split(',')
        pos2=sys.argv[3].split(',')
        vpos1=[]
        vpos2=[]
        for i in pos1:
                vpos1.append(int(i)-1)
        for i in pos2:
                vpos2.append(int(i)-1)
        all_fisher(filename,vpos1,vpos2)        
