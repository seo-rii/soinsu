from z3 import *
from solve_bits import N,LEAK,MASK,ct,e
import time
p=BitVec('p',1024); q=BitVec('q',1024)
s=Solver(); s.set('timeout',120000)
s.add((p & MASK)==LEAK)
s.add(Extract(1023,1023,p)==1, Extract(1023,1023,q)==1)
s.add(ZeroExt(1024,p)*ZeroExt(1024,q)==BitVecVal(N,2048))
print('check start',flush=True); st=time.time(); r=s.check(); print('res',r,'time',time.time()-st)
if r==sat:
 m=s.model(); pv=m[p].as_long(); qv=m[q].as_long(); print(hex(pv)); print(N%pv)
