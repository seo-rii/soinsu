import sys,time
sys.path.append('/mnt/data')
from try_lattice import get_polys_for_candidate, strip_f_factors_fast, recover_with_pair_fast, verify_xy
from rsa_partial_lib import *
primes=[1009,5003,10007,20011,30011,40009,50021,60013,70001,80021,90001,100003,110017,120011,130003,140009,150001,160001]
start=time.time();
for high in range(16):
  for low in range(16):
    C,polys=get_polys_for_candidate(low,high,m=4,t=3,mode='jm',limit=8)
    qs=[]
    for idx,(bb,c) in enumerate(polys):
      q,k=strip_f_factors_fast(c,C)
      if k==3 and len(q)>=3:
        qs.append((idx,q))
    if len(qs)<2:
      print('noqs',low,high,flush=True); continue
    cand=recover_with_pair_fast(qs[0][1],qs[1][1],primes,verbose=False)
    print('cand',low,high,'num',len(cand),'Mbits',cand[0][2].bit_length() if cand else 0,flush=True)
    for xv,yv,M in cand:
      ok,p=verify_xy(C,xv,yv)
      if ok:
        print('FOUND',low,high,hex(p),decrypt_from_p(p)); sys.exit(0)
print('not found',time.time()-start)
