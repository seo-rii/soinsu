import sys,time
sys.path.append('/mnt/data')
from try_lattice import get_polys_for_candidate, strip_f_factors_fast, recover_with_pair_fast_limited, verify_xy
from rsa_partial_lib import *
# 24 small primes ~ 10 bits => product > 2^230
primes=[1009,1013,1019,1021,1031,1033,1039,1049,1051,1061,1063,1069,1087,1091,1093,1097,1103,1109,1117,1123,1129,1151,1153,1163,1171]
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
    cand=recover_with_pair_fast_limited(qs[0][1],qs[1][1],primes,verbose=False,max_roots_per_prime=20,max_cand=50)
    print('cand',low,high,'num',len(cand),'Mbits',cand[0][2].bit_length() if cand else 0,flush=True)
    for xv,yv,M in cand:
      ok,p=verify_xy(C,xv,yv)
      if ok:
        print('FOUND',low,high,hex(p),decrypt_from_p(p)); sys.exit(0)
print('not found',time.time()-start)
