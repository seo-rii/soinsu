import sys,time
sys.path.append('/mnt/data')
from hm_bivar import iiter, poly_pow, poly_shift
from solve_bits import N, LEAK
from fpylll import IntegerMatrix, LLL

def build(m,t):
 def W(s,l): return ((1<<l)-1)<<s
 lo=0; hi=0; base=LEAK|(lo<<150)|(hi<<920); p0=base & ~W(265,155) & ~W(600,230)
 X=1<<155; Y=1<<230; ainv=pow(1<<265,-1,N); c0=(p0*ainv)%N; cy=((1<<600)*ainv)%N
 f={(1,0):1,(0,1):cy%N,(0,0):c0%N}; g=[]; monoms=[]; Xmuls=[]
 for ii in iiter(m,2):
  k,jy=ii; gi=poly_shift(poly_pow(f,k),0,jy,pow(N,max(t-k,0))); g.append(gi); monoms.append((k,jy)); Xmuls.append(X**k*Y**jy)
 n=len(g); M=[[0]*n for _ in range(n)]
 for i,gi in enumerate(g):
  for j in range(i+1): M[i][j]=gi.get(monoms[j],0)*Xmuls[j]
 return M
m,t=map(int,sys.argv[1:3]); method=sys.argv[3]; delta=float(sys.argv[4]); prec=int(sys.argv[5])
M=build(m,t); B=IntegerMatrix.from_matrix(M); st=time.time(); print('start',m,t,method,delta,prec,'dim',B.nrows,flush=True)
kwargs={}
if method!='none': kwargs['method']=method
if prec:
 kwargs['float_type']='mpfr'; kwargs['precision']=prec
else:
 kwargs['float_type']='double'
LLL.reduction(B, delta=delta, eta=0.51, **kwargs)
print('done',time.time()-st, 'b00bits',abs(int(B[0,0])).bit_length(), flush=True)
