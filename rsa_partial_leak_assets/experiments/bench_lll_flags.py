import sys,time
sys.path.append('/mnt/data')
from hm_bivar import iiter, poly_pow, poly_shift
from solve_bits import N, LEAK
from fpylll import IntegerMatrix, LLL

def build(m,t):
 def W(s,l): return ((1<<l)-1)<<s
 p0=LEAK & ~W(265,155) & ~W(600,230); X=1<<155; Y=1<<230; inv=pow(1<<265,-1,N); c0=(p0*inv)%N; cy=((1<<600)*inv)%N
 f={(1,0):1,(0,1):cy,(0,0):c0}; mons=[]; fac=[]; gs=[]
 for ii in iiter(m,2):
  k,j=ii; gs.append(poly_shift(poly_pow(f,k),0,j,pow(N,max(t-k,0)))); mons.append(ii); fac.append(X**k*Y**j)
 n=len(gs); M=[[0]*n for _ in range(n)]
 for i,gi in enumerate(gs):
  for j in range(i+1): M[i][j]=gi.get(mons[j],0)*fac[j]
 return M
B=IntegerMatrix.from_matrix(build(9,4)); st=time.time(); flags=int(sys.argv[1]); delta=float(sys.argv[2]); print('start flags',flags,delta,flush=True)
try:
 LLL.reduction(B,delta=delta,eta=0.51,method='proved',float_type='mpfr',precision=128,flags=flags)
 print('done',time.time()-st,flush=True)
except Exception as e: print('ERR',repr(e),time.time()-st,flush=True)
