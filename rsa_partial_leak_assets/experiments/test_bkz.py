import sys,time
sys.path.append('/mnt/data')
from solve7_compact import build,BASE
from fpylll import IntegerMatrix, LLL, BKZ, GSO, FPLLL
m,t=map(int,sys.argv[1:3])
M,mons,sc=build(BASE,m,t,'y')
print('dim',M.nrows,'maxbit',max(abs(int(M[i,j])).bit_length() for i in range(M.nrows) for j in range(M.ncols) if int(M[i,j])))
# try wrapper LLL first low delta then BKZ?
t0=time.time()
try:
    LLL.reduction(M,delta=float(sys.argv[3]) if len(sys.argv)>3 else .75,method='proved',float_type='mpfr',precision=128)
    print('lll done',time.time()-t0)
except Exception as e: print('lll err',e)
try:
    par=BKZ.Param(block_size=int(sys.argv[4]) if len(sys.argv)>4 else 10,max_loops=1)
    t1=time.time(); BKZ.reduction(M,par); print('bkz done',time.time()-t1,'total',time.time()-t0)
except Exception as e: print('bkz err',type(e),e)
print([sum(abs(int(M[i,j])) for j in range(M.ncols)).bit_length() for i in range(min(10,M.nrows))])
