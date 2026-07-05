import sys, time
sys.path.append('/mnt/data')
from rsa_partial_lib import *

configs=[
('std3',3,std_pairs(3),True),
('std4',4,std_pairs(4),True),
('box11_m4',4,box_pairs(1,1),True),
('box21_m4',4,box_pairs(2,1),True),
('box22_m4',4,box_pairs(2,2),True),
('hm2_m4',4,hm_pairs(4,2),True),
('hm3_m4',4,hm_pairs(4,3),True),
('box11_m5',5,box_pairs(1,1),True),
('box21_m5',5,box_pairs(2,1),True),
]
for name,m,pairs,inc in configs:
    print('CONFIG',name,flush=True)
    rows,mons,sh=build_rows(base0,m,pairs,inc)
    print(' rows cols',len(rows),len(mons),'maxbit',max(abs(v).bit_length() for row in rows for v in row if v),flush=True)
    try:
        R=reduce_rows(rows,precision=256)
    except Exception as e:
        print('ERR',e,flush=True); continue
    non=[]
    for i,row in enumerate(R):
        b=vector_bits(row)
        if b: non.append((i,b,sum(1 for v in row if v)))
    print(' first',non[:15],flush=True)
