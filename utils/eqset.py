#!/usr/bin/env python
import sys

f1=open(sys.argv[1],'r').readlines()
f2=open(sys.argv[2],'r').readlines()
pos1=sys.argv[3].split(',')
vpos1=[]
pos2=sys.argv[4].split(',')
vpos2=[]
try:
        wall=sys.argv[5]
        if wall.upper()=='TRUE': wall=True
except:
        wall=True

tab=False
try:           
        #if sys.argv[6].upper()=='Y': tab=True
        if sys.argv[6]: tab=sys.argv[6]
except:
        tab=False


for i in pos1:
        vpos1.append(int(i)-1)

for i in pos2:
        vpos2.append(int(i)-1)

def createdict(lista,pos,tab=False):
        df={}
        for i in lista:
                if tab:
                        v=i.rstrip().split(tab)        
                else:
                        v=i.rstrip().split()
                name=[]
                try:
                        for j in pos:
                                name.append(v[j])
                except:
                        print ("WARNING: incorrect line",i.rstrip(), file=sys.stderr)        
                        continue
                name=tuple(name)
                if df.get(name,[])==[]: df[name]=''
                df[name]=df[name]+i
                #df[name]=i
                #print name,df[name]
        return df

def difference(s1,s2):
        v1=sets.Set(s1)
        diff=v1.difference(s2)
        return diff

df1=createdict(f1,vpos1,tab)
df2=createdict(f2,vpos2,tab)
d1=set(df1.keys())
d2=set(df2.keys())
dd=d1.intersection(d2)
#d1=sets.Set(df1)
#d2=sets.Set(df2)
#print len(df1.keys()),len(df2.keys())
#dd=d1 & d2
#dd=d1-d2
#dd=difference(df1.keys(),df2.keys())
for i in dd:
        if wall==True:
                if not tab: tab=' '
                print (df1[i][:-1]+tab+df2[i][:-1])
        else:
                print (df1[i][:-1])
