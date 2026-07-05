import sys,time
sys.path.append('/mnt/data')
from hm_bivar import linear_bivar_multilead,recover_roots_crt
from solve_bits import N,LEAK,MASK

def W(s,l): return ((1<<l)-1)<<s
lo=0;hi=0;X=1<<155;Y=1<<230;p0=(LEAK|(lo<<150)|(hi<<920)) & ~W(265,155) & ~W(600,230)
st=time.time(); hs,_=linear_bivar_multilead(N,p0,1<<265,1<<600,X,Y,m=8,t=4,lll_rows=20,verbose=True); print('LLL time',time.time()-st)
def verify(x,y):
 p=p0+(1<<265)*x+(1<<600)*y
 print('verify cand',x.bit_length(),y.bit_length(),N%p==0, flush=True)
 if N%p==0 and (p&MASK)==LEAK: return p
 return None
st=time.time(); res=recover_roots_crt(hs,X,Y,verify,primes=[10007,10009,10037,10039,10061,10067,10069,10079],max_pairs=30,verbose=True); print('res',res,'time',time.time()-st)
