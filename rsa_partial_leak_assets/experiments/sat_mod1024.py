import sys,time,argparse
sys.path.append('/mnt/data')
from sat_factor import SatMul, q_low_for_low, q_high_common_for_high
from rsa_partial_lib import N,mask,leak,decrypt_from_p

def build_solve(low,high,tlim=20,threads=1,enum=1):
 from pycryptosat import Solver
 sm=SatMul(); sm.s=Solver(threads=threads, verbose=0)
 nbits=1024
 qhp,pref,_,_=q_high_common_for_high(high)
 pbits=[]
 for i in range(nbits):
  if (mask>>i)&1: pbits.append(bool((leak>>i)&1))
  elif 150<=i<154: pbits.append(bool((low>>(i-150))&1))
  elif 920<=i<924: pbits.append(bool((high>>(i-920))&1))
  else: pbits.append(sm.new())
 qbits=[]
 for j in range(nbits):
  if j in qhp: qbits.append(bool(qhp[j]))
  elif j==0 or j==1023: qbits.append(True)
  else: qbits.append(sm.new())
 # q lower from low guess
 ql=q_low_for_low(low)
 for j,val in enumerate(ql):
  if isinstance(qbits[j],bool):
   if qbits[j]!=bool(val): print('contr qlow'); return None
  else: sm.eq_const(qbits[j],bool(val))
 cols=[[] for _ in range(nbits+1)]
 terms=0;t0=time.time()
 for i,pb in enumerate(pbits):
  if pb is False: continue
  maxj=nbits-1-i
  for j in range(maxj+1):
   qb=qbits[j]
   term=qb if pb is True else sm.and_bit(pb,qb)
   if term is not False: cols[i+j].append(term); terms+=1
 for k in range(nbits):
  col=cols[k]
  while len(col)>=3:
   a=col.pop(); b=col.pop(); c=col.pop(); s,carry=sm.full_adder(a,b,c)
   if s is not False: col.append(s)
   if carry is not False and k+1 < len(cols): cols[k+1].append(carry)
  if len(col)==2:
   a=col.pop(); b=col.pop(); s,carry=sm.half_adder(a,b)
   if s is not False: col.append(s)
   if carry is not False and k+1 < len(cols): cols[k+1].append(carry)
  bit=bool((N>>k)&1)
  if len(col)==0:
   if bit: print('contr at',k); return None
  else: sm.eq_const(col[0],bit)
 # ignore carry beyond 1023 (mod condition)
 print('built',low,high,'pref',pref,'vars',sm.var,'clauses',sm.clauses,'xors',sm.xors,'terms',terms,'time',time.time()-t0,flush=True)
 sols=0
 while sols<enum:
  t=time.time(); sat,sol=sm.s.solve(time_limit=tlim)
  print('solve',low,high,'sat',sat,'dt',time.time()-t,flush=True)
  if sat is not True: return None
  p=0; q=0; block=[]
  for i,b in enumerate(pbits):
   val=b if isinstance(b,bool) else (sol[abs(b)] ^ (b<0))
   if val: p|=1<<i
   if not isinstance(b,bool): block.append(-b if val else b)
  for j,b in enumerate(qbits):
   val=b if isinstance(b,bool) else (sol[abs(b)] ^ (b<0))
   if val: q|=1<<j
  print('model',sols,'N%p',N%p,'qmatch',q==N//p if p else False,'pbitsok',(p&mask)==leak, flush=True)
  if N%p==0 and (p&mask)==leak:
   print('FOUND',hex(p),decrypt_from_p(p)); return p
  sm.add_clause(block) # block exact p assignment? clause of opposite literals
  sols+=1
 return None
if __name__=='__main__':
 ap=argparse.ArgumentParser(); ap.add_argument('low',type=int); ap.add_argument('high',type=int); ap.add_argument('--time',type=float,default=20); ap.add_argument('--threads',type=int,default=1); ap.add_argument('--enum',type=int,default=3)
 a=ap.parse_args(); build_solve(a.low,a.high,a.time,a.threads,a.enum)
