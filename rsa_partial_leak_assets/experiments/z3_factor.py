import sys,time,argparse
sys.path.append('/mnt/data')
from rsa_partial_lib import N,mask,leak,decrypt_from_p
from sat_factor import q_low_for_low, q_high_common_for_high
import z3

def solve(low,high,timeout_ms=60000):
    p=z3.BitVec('p',1024); q=z3.BitVec('q',1024)
    s=z3.SolverFor('QF_BV'); s.set('timeout', timeout_ms)
    s.add((p & z3.BitVecVal(mask,1024)) == z3.BitVecVal(leak,1024))
    # low/high guesses
    for b in range(4):
        s.add(z3.Extract(150+b,150+b,p) == z3.BitVecVal((low>>b)&1,1))
        s.add(z3.Extract(920+b,920+b,p) == z3.BitVecVal((high>>b)&1,1))
    ql=q_low_for_low(low)
    # combine q low bits as mask eq
    qlval=sum(int(v)<<j for j,v in enumerate(ql))
    s.add(z3.Extract(264,0,q) == z3.BitVecVal(qlval,265))
    qhp,pref,_,_=q_high_common_for_high(high)
    # top prefix eq as extract from 1023 down to 1024-pref
    if pref>0:
        lo=1024-pref; hi=1023
        val=0
        for j in range(lo,hi+1): val |= qhp[j] << (j-lo)
        s.add(z3.Extract(hi,lo,q) == z3.BitVecVal(val,pref))
    s.add(z3.Extract(1023,1023,p)==z3.BitVecVal(1,1))
    s.add(z3.Extract(1023,1023,q)==z3.BitVecVal(1,1))
    prod=z3.ZeroExt(1024,p)*z3.ZeroExt(1024,q)
    s.add(prod == z3.BitVecVal(N,2048))
    print('built constraints',len(s.assertions()),'pref',pref,flush=True)
    t=time.time(); r=s.check(); print('check',r,'time',time.time()-t,flush=True)
    if r==z3.sat:
        m=s.model(); pv=m[p].as_long(); qv=m[q].as_long(); print(hex(pv)); print('check',pv*qv==N,N%pv==0,(pv&mask)==leak); print(decrypt_from_p(pv))
if __name__=='__main__':
    ap=argparse.ArgumentParser(); ap.add_argument('low',type=int); ap.add_argument('high',type=int); ap.add_argument('--timeout',type=int,default=60000)
    a=ap.parse_args(); solve(a.low,a.high,a.timeout)
