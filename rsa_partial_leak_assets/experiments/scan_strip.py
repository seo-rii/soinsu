import sys,time,itertools
sys.path.append('/mnt/data')
from try_lattice import get_polys_for_candidate, strip_f_factors, recover_with_pair_fast, verify_xy
from rsa_partial_lib import *
# primes around increasing sizes
primes=[1009,1013,1019,1021,1031,1033,1039,1049,1051,1061,1063,1069,1087,1091,1093,1097,1103,1109,1117,1123,1129,1151,1153,1163,1171]
# use small primes; product bits ~250
start=time.time(); hits=[]
mode='jm'; m=4; t=3
for high in range(16):
  for low in range(16):
    C,polys=get_polys_for_candidate(low,high,m=m,t=t,mode=mode,limit=14)
    qs=[]
    for idx,(b,c) in enumerate(polys):
      q,k=strip_f_factors(c,C)
      # useful: stripped by m-1 or less, nontrivial multivariate, not constant/monomial
      if k<=m-1 and len(q)>=3:
        deg=max(i+j for i,j in q)
        qs.append((idx,b,k,deg,q))
    print('cand',low,high,'qs',[(i,b,k,d,len(q)) for i,b,k,d,q in qs[:8]],flush=True)
    tried=0
    for a,b in itertools.combinations(range(min(len(qs),8)),2):
      ia,ba,ka,da,qa=qs[a]; ib,bb,kb,db,qb=qs[b]
      tried+=1
      cand=recover_with_pair_fast(qa,qb,primes,verbose=False)
      if cand:
        print(' possible',low,high,'pair',ia,ib,'num',len(cand),'Mbits',cand[0][2].bit_length(),flush=True)
      for xv,yv,M in cand:
        ok,p=verify_xy(C,xv,yv)
        if ok:
          print('FOUND low high',low,high,'pair',ia,ib)
          print(hex(p))
          print(decrypt_from_p(p))
          sys.exit(0)
    print(' done cand tried',tried,flush=True)
print('not found time',time.time()-start)
