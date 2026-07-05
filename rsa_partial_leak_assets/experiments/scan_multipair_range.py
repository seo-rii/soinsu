import sys,time,itertools
sys.path.append('/mnt/data')
from try_lattice import get_polys_for_candidate, strip_f_factors_fast, recover_with_pair_fast_limited, verify_xy
from rsa_partial_lib import *
primes=[1009,1013,1019,1021,1031,1033,1039,1049,1051,1061,1063,1069,1087,1091,1093,1097,1103,1109,1117,1123,1129,1151,1153,1163,1171]
h0=int(sys.argv[1]); h1=int(sys.argv[2]); mode=sys.argv[3] if len(sys.argv)>3 else 'jm'; m=int(sys.argv[4]) if len(sys.argv)>4 else 4; t=int(sys.argv[5]) if len(sys.argv)>5 else 3
start=time.time(); totalpairs=0
for high in range(h0,h1):
  for low in range(16):
    C,polys=get_polys_for_candidate(low,high,m=m,t=t,mode=mode,limit=16)
    qs=[]
    for idx,(bb,c) in enumerate(polys):
      q,k=strip_f_factors_fast(c,C)
      if k<=m-1 and len(q)>=3:
        deg=max(i+j for i,j in q)
        # avoid too high degree? keep first 8 useful
        qs.append((idx,bb,k,deg,q))
    qs=qs[:8]
    hits=[]; tried=0
    for a,b in itertools.combinations(range(len(qs)),2):
      ia,ba,ka,da,qa=qs[a]; ib,bb,kb,db,qb=qs[b]
      tried+=1; totalpairs+=1
      # quick: first 5 primes only
      cand=recover_with_pair_fast_limited(qa,qb,primes[:5],max_roots_per_prime=10,max_cand=20)
      if not cand: continue
      # full
      cand=recover_with_pair_fast_limited(qa,qb,primes,max_roots_per_prime=10,max_cand=20)
      if cand:
        hits.append((ia,ib,len(cand),cand[0][2].bit_length()))
        for xv,yv,M in cand:
          ok,p=verify_xy(C,xv,yv)
          if ok:
            print('FOUND',low,high,'pair',ia,ib,hex(p),decrypt_from_p(p)); sys.exit(0)
    print('cand',low,high,'qs',[(i,b,k,d) for i,b,k,d,q in qs],'tried',tried,'hits',hits,flush=True)
print('range done',h0,h1,'time',time.time()-start,'pairs',totalpairs)
