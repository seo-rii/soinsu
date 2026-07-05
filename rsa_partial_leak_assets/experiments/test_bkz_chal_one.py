import sys,time
sys.path.append('/mnt/data')
from hm_bivar import hm_polys_bivar_bkz,recover_roots_crt
from solve_bits import N,LEAK,MASK

def W(s,l): return ((1<<l)-1)<<s
lo=0;hi=0;X=1<<155;Y=1<<230;p0=(LEAK|(lo<<150)|(hi<<920)) & ~W(265,155) & ~W(600,230)
inv=pow(1<<265,-1,N); c0=(p0*inv)%N; cy=((1<<600)*inv)%N
st=time.time(); hs,_=hm_polys_bivar_bkz(N,c0,cy,X,Y,m=8,t=4,block=10,lll_rows=12,verbose=True); print('time',time.time()-st)
def verify(x,y):
 p=p0+(1<<265)*x+(1<<600)*y
 print('verify',N%p==0,flush=True)
 if N%p==0 and (p&MASK)==LEAK: return p
 return None
print(recover_roots_crt(hs,X,Y,verify,primes=[10007,10009,10037],max_pairs=10,verbose=True))
