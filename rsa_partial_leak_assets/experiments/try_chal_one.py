import sys, time
sys.path.append('/mnt/data')
from hm_bivar import hm_polys_bivar, recover_roots_crt
from solve_bits import N,e,ct,MASK,LEAK

def W(s,l): return ((1<<l)-1)<<s
lo=0; hi=0
base=LEAK | (lo<<150) | (hi<<920)
WX=W(265,155); WY=W(600,230)
p0=base & ~WX & ~WY
ainv=pow(1<<265,-1,N)
c0=(p0*ainv)%N
cy=((1<<600)*ainv)%N
X=1<<155; Y=1<<230
for m,t in [(6,3),(7,3),(8,4),(9,4)]:
    st=time.time(); hs,mons=hm_polys_bivar(N,c0,cy,X,Y,m=m,t=t,lll_rows=8,verbose=True); print('time',time.time()-st)
