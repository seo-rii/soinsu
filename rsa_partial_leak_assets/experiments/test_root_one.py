import sys,time
sys.path.append('/mnt/data')
from hm_bivar import hm_polys_bivar, recover_roots_crt
from solve_bits import N,LEAK,MASK,ct,e

def W(s,l): return ((1<<l)-1)<<s
lo=0; hi=0
base=LEAK | (lo<<150) | (hi<<920)
p0=base & ~W(265,155) & ~W(600,230)
ainv=pow(1<<265,-1,N); c0=(p0*ainv)%N; cy=((1<<600)*ainv)%N
X=1<<155; Y=1<<230
hs,_=hm_polys_bivar(N,c0,cy,X,Y,m=8,t=4,lll_rows=12,verbose=True)
print('hs ready')
def verify(x,y):
 p=p0+(1<<265)*x+(1<<600)*y
 print('verify candidate',x.bit_length(),y.bit_length(),N%p==0, flush=True)
 if N%p==0 and (p&MASK)==LEAK:
  return p
 return None
st=time.time(); res=recover_roots_crt(hs,X,Y,verify,primes=[10007,10009,10037,10039,10061],max_pairs=6,verbose=True); print('res',res,'time',time.time()-st)
