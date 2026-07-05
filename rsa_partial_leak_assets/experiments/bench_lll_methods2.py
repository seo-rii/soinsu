import sys,time
sys.path.append('/mnt/data')
from hm_bivar import iiter, poly_pow, poly_shift
from solve_bits import N, LEAK
from fpylll import IntegerMatrix, LLL

def build(m,t):
 def W(s,l): return ((1<<l)-1)<<s
 p0=(LEAK & ~W(265,155) & ~W(600,230)); X=1<<155; Y=1<<230; inv=pow(1<<265,-1,N); c0=(p0*inv)%N; cy=((1<<600)*inv)%N
 f={(1,0):1,(0,1):cy,(0,0):c0}; g=[]; mons=[]; fac=[]
 for ii in iiter(m,2):
  k,j=ii; g.append(poly_shift(poly_pow(f,k),0,j,pow(N,max(t-k,0)))); mons.append(ii); fac.append(X**k*Y**j)
 n=len(g); M=[[0]*n for _ in range(n)]
 for i,gi in enumerate(g):
  for j in range(i+1): M[i][j]=gi.get(mons[j],0)*fac[j]
 return M
m,t=9,4; method=sys.argv[1]; ft=sys.argv[2]; delta=float(sys.argv[3])
B=IntegerMatrix.from_matrix(build(m,t)); print('start',method,ft,delta,flush=True); st=time.time()
try:
 LLL.reduction(B,delta=delta,eta=0.51,method=method,float_type=ft)
 print('done',time.time()-st,flush=True)
except Exception as e: print('ERR',repr(e),time.time()-st,flush=True)
