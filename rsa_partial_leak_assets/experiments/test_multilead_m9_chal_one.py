import sys,time
sys.path.append('/mnt/data')
from hm_bivar import linear_bivar_multilead,recover_roots_crt
from solve_bits import N,LEAK,MASK

def W(s,l): return ((1<<l)-1)<<s
X=1<<155;Y=1<<230;p0=LEAK & ~W(265,155)&~W(600,230)
st=time.time(); hs,_=linear_bivar_multilead(N,p0,1<<265,1<<600,X,Y,m=9,t=4,lll_rows=110,verbose=True); print('time',time.time()-st); hs=[h for h in hs if h]; print('nonzero',len(hs))
def verify(x,y):
 p=p0+(1<<265)*x+(1<<600)*y
 print('verify',N%p==0,flush=True)
 if N%p==0 and (p&MASK)==LEAK: return p
 return None
print(recover_roots_crt(hs[:25],X,Y,verify,primes=[10007,10009,10037,10039,10061],max_pairs=30,verbose=True))
