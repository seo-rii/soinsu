import sys,time,math
sys.path.append('/mnt/data')
from hm_linear_bivar import build_hm_linear_rows,C_for,row_l1_bits,row_max_bits,unscale_row,N
from fpylll import IntegerMatrix, LLL
from math import gcd
for m,t in [(10,3),(13,4)]:
 rows,mons,scales,meta=build_hm_linear_rows(C_for(0,0),m,t,include_dups=False)
 rows2=[]; gs=[]
 for r in rows:
  g=0
  for v in r:
   if v: g=gcd(g,abs(v))
  if g==0: g=1
  gs.append(g)
  rows2.append([v//g for v in r])
 rows2=sorted(rows2,key=lambda r: row_max_bits(r))
 print('mt',m,t,'rows',len(rows2),'cols',len(rows2[0]),'normrange',row_max_bits(rows2[0]),row_max_bits(rows2[-1]), flush=True)
 A=IntegerMatrix(len(rows2),len(rows2[0]))
 for i,r in enumerate(rows2):
  for j,v in enumerate(r): A[i,j]=int(v)
 st=time.time()
 try:
  LLL.reduction(A,delta=0.99,method='proved',float_type='mpfr',precision=192)
  print('done',time.time()-st,[row_l1_bits([int(A[i,j]) for j in range(A.ncols)]) for i in range(20)],'hbits',int(math.log2(N)*.499*t),flush=True)
 except Exception as e:
  print('err',repr(e),time.time()-st,flush=True)
