import sys,time,itertools
sys.path.append('/mnt/data')
from try_lattice import get_polys_for_candidate, roots_pair_mod, recover_with_pair, verify_xy
from rsa_partial_lib import *

mode='jm'; m=4; t=3
pr=1000003
hits=[]
start=time.time()
for high in range(16):
  for low in range(16):
    C,polys=get_polys_for_candidate(low,high,m,t,mode,limit=14)
    found_roots=[]
    # try pairs excluding row0? all first 10
    for ia in range(min(10,len(polys))):
      for ib in range(ia+1,min(10,len(polys))):
        roots=roots_pair_mod(polys[ia][1],polys[ib][1],pr)
        if roots:
          found_roots.append((ia,ib,polys[ia][0],polys[ib][0],len(roots),roots[:3]))
          break
      if found_roots: break
    if found_roots:
      print('HIT low high',low,high,found_roots[:3],flush=True)
      hits.append((low,high,found_roots))
    else:
      print('no',low,high,flush=True)
print('done hits',hits,'time',time.time()-start)
