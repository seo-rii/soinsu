import sys,time
sys.path.append('/mnt/data')
from rsa_partial_lib import *
configs=[]
for m in [5,6,7,8,9,10]:
  for t in [2,3,4,5]:
    configs.append((f'hm{t}_m{m}',m,hm_pairs(m,t),True))
# also boxes
for m in [6,7,8,9]:
  for tx,ty in [(2,2),(3,2),(3,3),(4,3)]:
    configs.append((f'box{tx}{ty}_m{m}',m,box_pairs(tx,ty),True))
for name,m,pairs,inc in configs:
    rows,mons,sh=build_rows(base0,m,pairs,inc)
    if len(rows)>250 or len(mons)>150: continue
    print('CONFIG',name,'rows',len(rows),'cols',len(mons),flush=True)
    try:
        R=reduce_rows(rows,precision=256)
    except Exception as e:
        print('ERR',e,flush=True); continue
    non=[]
    for i,row in enumerate(R):
        b=vector_bits(row)
        if b: non.append((i,b,sum(1 for v in row if v)))
    print(' first',non[:20],flush=True)
