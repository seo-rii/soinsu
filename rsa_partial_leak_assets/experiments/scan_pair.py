import sys,time
sys.path.append('/mnt/data')
from try_lattice import get_polys_for_candidate, roots_pair_mod
from rsa_partial_lib import *
mode='jm'; m=4; t=3; pr=1000003
pairs_to_try=[(1,2),(1,3),(2,3),(3,4),(4,5),(7,8)]
start=time.time(); hits=[]
for high in range(16):
  for low in range(16):
    C,polys=get_polys_for_candidate(low,high,m,t,mode,limit=12)
    for ia,ib in pairs_to_try:
      roots=roots_pair_mod(polys[ia][1],polys[ib][1],pr)
      if roots:
        print('HIT',low,high,ia,ib,len(roots),roots[:5],flush=True); hits.append((low,high,ia,ib,roots[:5])); break
    else:
      print('no',low,high,flush=True)
print('done',hits,'time',time.time()-start)
