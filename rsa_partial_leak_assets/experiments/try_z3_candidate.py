from z3 import *
from solve_bits import N,LEAK,MASK
import time, sys
lo=int(sys.argv[1]) if len(sys.argv)>1 else 0; hi=int(sys.argv[2]) if len(sys.argv)>2 else 0
unknown_blocks=[(150,4),(265,84),(362,58),(600,69),(682,87),(784,46),(920,4)]
def W(s,l): return ((1<<l)-1)<<s
base=LEAK | (lo<<150) | (hi<<920)
# p low known to 265
p_low=base & ((1<<265)-1)
q_low=(N*pow(p_low,-1,1<<265))%(1<<265)
# p range with hi fixed and lo fixed
mask_unknown = W(265,84)|W(362,58)|W(600,69)|W(682,87)|W(784,46)
p_min=base & ~mask_unknown
p_max=p_min | mask_unknown
q_min=N//p_max
q_max=N//p_min
# common prefix length
xor=q_min^q_max; pref=1024-xor.bit_length() if xor else 1024
q_pref=q_min>>(1024-pref) if pref>0 else 0
print('lo hi',lo,hi,'q_low',hex(q_low),'qpref',pref,hex(q_pref),flush=True)
p=BitVec('p',1024); q=BitVec('q',1024)
s=Solver(); s.set('timeout',120000)
s.add((p & MASK)==LEAK)
s.add(Extract(153,150,p)==lo, Extract(923,920,p)==hi)
s.add((q & ((1<<265)-1))==q_low)
if pref>0: s.add(Extract(1023,1024-pref,q)==q_pref)
s.add(UGE(p,BitVecVal(p_min,1024)), ULE(p,BitVecVal(p_max,1024)))
s.add(UGE(q,BitVecVal(q_min,1024)), ULE(q,BitVecVal(q_max,1024)))
s.add(ZeroExt(1024,p)*ZeroExt(1024,q)==BitVecVal(N,2048))
st=time.time(); print('check',flush=True); r=s.check(); print('res',r,'time',time.time()-st,flush=True)
if r==sat:
 pv=s.model()[p].as_long(); print(hex(pv), N%pv)
