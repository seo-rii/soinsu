import z3,sys,time
sys.path.append('/mnt/data')
from rsa_partial_lib import N,mask,leak,ct,e,decrypt_from_p
XS,XL,YS,YL=265,155,600,230
BASE=leak & ~(((1<<XL)-1)<<XS) & ~(((1<<YL)-1)<<YS)
lo=int(sys.argv[1]); hi=int(sys.argv[2]); timeout=int(sys.argv[3]) if len(sys.argv)>3 else 60000
x=z3.BitVec('x',XL); y=z3.BitVec('y',YL)
p = z3.BitVecVal(BASE | (lo<<150) | (hi<<920),1024) + (z3.ZeroExt(1024-XL,x) << XS) + (z3.ZeroExt(1024-YL,y) << YS)
s=z3.SolverFor('QF_BV'); s.set('timeout',timeout)
s.add((p & z3.BitVecVal(mask,1024)) == z3.BitVecVal(leak,1024))
N2048=z3.BitVecVal(N,2048); p2048=z3.ZeroExt(1024,p)
s.add(z3.URem(N2048,p2048)==z3.BitVecVal(0,2048))
q=z3.UDiv(N2048,p2048)
s.add(z3.UGE(q,z3.BitVecVal(1<<1023,2048)), z3.ULT(q,z3.BitVecVal(1<<1024,2048)))
print('start',lo,hi,'assertions',len(s.assertions()),flush=True); t=time.time(); r=s.check(); print('res',r,'time',time.time()-t,flush=True)
if r==z3.sat:
 m=s.model(); xv=m[x].as_long(); yv=m[y].as_long(); pv=(BASE|(lo<<150)|(hi<<920))+(xv<<XS)+(yv<<YS); print(hex(pv),N%pv,(pv&mask)==leak); print(decrypt_from_p(pv))
