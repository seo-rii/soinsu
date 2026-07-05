import sys,time,math
from math import gcd
sys.path.append('/mnt/data')
from solve7_compact import build,BASE
from fpylll import IntegerMatrix,LLL
m,t=map(int,sys.argv[1:3])
M,mons,sc=build(BASE,m,t,'y')
G=[]
for j in range(M.ncols):
    g=0
    for i in range(M.nrows): g=gcd(g,abs(int(M[i,j])))
    G.append(g or 1)
R=IntegerMatrix(M.nrows,M.ncols)
for i in range(M.nrows):
    for j in range(M.ncols): R[i,j]=int(M[i,j])//G[j]
print('dim',R.nrows,R.ncols,'maxbit',max(abs(int(R[i,j])).bit_length() for i in range(R.nrows) for j in range(R.ncols) if int(R[i,j])))
t0=time.time(); LLL.reduction(R,delta=float(sys.argv[3]) if len(sys.argv)>3 else .99,method='proved',float_type='mpfr',precision=int(sys.argv[4]) if len(sys.argv)>4 else 256); print('done',time.time()-t0)
print([sum(abs(int(R[i,j])) for j in range(R.ncols)).bit_length() for i in range(min(5,R.nrows))])
