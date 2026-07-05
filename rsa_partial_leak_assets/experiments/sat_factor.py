import sys,time,math
sys.path.append('/mnt/data')
from rsa_partial_lib import N,mask,leak,ct,e,decrypt_from_p
from pycryptosat import Solver

class SatMul:
    def __init__(self):
        self.s=Solver(verbose=0)
        self.var=0
        self.clauses=0; self.xors=0
    def new(self):
        self.var+=1; return self.var
    def add_clause(self,cl):
        # simplify constants? clauses have ints only
        self.s.add_clause(cl); self.clauses+=1
    def add_xor_clause(self,vars,rhs):
        self.s.add_xor_clause(vars, rhs); self.xors+=1
    def is_const(self,b): return isinstance(b,bool)
    def notb(self,b):
        return (not b) if isinstance(b,bool) else -b
    def eq_const(self,b,val):
        if isinstance(b,bool):
            if b!=val: raise ValueError('contradiction')
        else:
            self.add_clause([b if val else -b])
    def xor_reduce(self,bits):
        rhs=False; vars=[]
        for b in bits:
            if isinstance(b,bool): rhs ^= b
            else:
                if b<0:
                    rhs ^= True; vars.append(-b)
                else: vars.append(b)
        # remove duplicate variables mod 2
        # Use dict parity
        if len(vars)>1:
            par={}
            for v in vars: par[v]=par.get(v,0)^1
            vars=[v for v,p in par.items() if p]
        return vars,rhs
    def xor_bit(self,bits):
        vars,rhs=self.xor_reduce(bits)
        if not vars: return rhs
        if len(vars)==1: return -vars[0] if rhs else vars[0]
        s=self.new()
        # xor(vars) xor s = rhs
        self.add_xor_clause(vars+[s], rhs)
        return s
    def and_bit(self,a,b):
        if isinstance(a,bool): return b if a else False
        if isinstance(b,bool): return a if b else False
        # if same literal or complements
        if a==b: return a
        if a==-b: return False
        c=self.new()
        # c <-> a & b, where a,b are literals for true
        self.add_clause([-a, -b, c])
        self.add_clause([a, -c])
        self.add_clause([b, -c])
        return c
    def or_bit(self,a,b):
        if isinstance(a,bool): return True if a else b
        if isinstance(b,bool): return True if b else a
        if a==b: return a
        if a==-b: return True
        c=self.new()
        self.add_clause([a,b,-c])
        self.add_clause([-a,c])
        self.add_clause([-b,c])
        return c
    def maj3(self,a,b,c):
        consts=[x for x in (a,b,c) if isinstance(x,bool)]
        lits=[x for x in (a,b,c) if not isinstance(x,bool)]
        ones=sum(1 for x in consts if x)
        if ones>=2: return True
        if ones==1:
            if len(lits)==0: return False
            if len(lits)==1: return lits[0]
            return self.or_bit(lits[0],lits[1])
        # no true constants
        if len(lits)==0: return False
        if len(lits)==1: return False
        if len(lits)==2: return self.and_bit(lits[0],lits[1])
        a,b,c=lits
        # simplify equal/complement? generic clauses handle tautologies? Need avoid 0 lits? We'll just add; tautological okay?
        # handle duplicates by truth table fallback for simplicity
        if len(set([abs(a),abs(b),abs(c)]))<3:
            # compute majority using and/or: (a&b)|(a&c)|(b&c)
            return self.or_bit(self.or_bit(self.and_bit(a,b), self.and_bit(a,c)), self.and_bit(b,c))
        out=self.new()
        self.add_clause([-a,-b,out])
        self.add_clause([-a,-c,out])
        self.add_clause([-b,-c,out])
        self.add_clause([a,b,-out])
        self.add_clause([a,c,-out])
        self.add_clause([b,c,-out])
        return out
    def half_adder(self,a,b):
        return self.xor_bit([a,b]), self.and_bit(a,b)
    def full_adder(self,a,b,c):
        return self.xor_bit([a,b,c]), self.maj3(a,b,c)

def factor_with_sat(low_guess=None, high_guess=None, q_low_bits=None, q_high_prefix=None, time_limit=None):
    sm=SatMul()
    nbits=1024
    pbits=[]; unknown_p=[]
    for i in range(nbits):
        if (mask>>i)&1:
            pbits.append(bool((leak>>i)&1))
        else:
            v=sm.new(); pbits.append(v); unknown_p.append((i,v))
    # optional fix low/high nibble
    if low_guess is not None:
        for b in range(4): sm.eq_const(pbits[150+b], bool((low_guess>>b)&1))
    if high_guess is not None:
        for b in range(4): sm.eq_const(pbits[920+b], bool((high_guess>>b)&1))
    qbits=[]
    for j in range(nbits):
        v=sm.new(); qbits.append(v)
    # q odd and top bit 1
    sm.eq_const(qbits[0], True); sm.eq_const(qbits[1023], True)
    if q_low_bits is not None:
        for j,val in enumerate(q_low_bits):
            sm.eq_const(qbits[j], bool(val))
    if q_high_prefix is not None:
        # q_high_prefix is dict bit->val
        for j,val in q_high_prefix.items(): sm.eq_const(qbits[j], bool(val))
    maxcols=2100
    cols=[[] for _ in range(maxcols)]
    t0=time.time(); terms=0
    for i,pb in enumerate(pbits):
        if pb is False: continue
        for j,qb in enumerate(qbits):
            term = qb if pb is True else sm.and_bit(pb,qb)
            if term is False: continue
            cols[i+j].append(term); terms+=1
        if i%128==0: print('partial i',i,'vars',sm.var,'clauses',sm.clauses,flush=True)
    print('initial terms',terms,'vars',sm.var,'clauses',sm.clauses,'xors',sm.xors,'time',time.time()-t0,flush=True)
    for k in range(maxcols-1):
        col=cols[k]
        while len(col)>=3:
            a=col.pop(); b=col.pop(); c=col.pop()
            s,carry=sm.full_adder(a,b,c)
            if s is not False: col.append(s)
            # if s false, append nothing? wait sum bit 0 at same weight, reducing by3, okay.
            if carry is not False: cols[k+1].append(carry)
        if len(col)==2:
            a=col.pop(); b=col.pop()
            s,carry=sm.half_adder(a,b)
            if s is not False: col.append(s)
            if carry is not False: cols[k+1].append(carry)
        # now len 0 or 1; constrain bit
        Nbit = bool((N>>k)&1) if k<2048 else False
        if len(col)==0:
            if Nbit: raise ValueError('Nbit 1 but empty col')
        elif len(col)==1:
            sm.eq_const(col[0], Nbit)
        else:
            raise AssertionError('col len >1')
        if k%128==0: print('reduce k',k,'vars',sm.var,'clauses',sm.clauses,'xors',sm.xors,'nextlen',len(cols[k+1]),flush=True)
    # any leftover in last column maxcols-1 must be zero
    for b in cols[-1]: sm.eq_const(b,False)
    print('built vars',sm.var,'clauses',sm.clauses,'xors',sm.xors,flush=True)
    t=time.time(); sat,sol=sm.s.solve()
    print('solve result',sat,'time',time.time()-t,flush=True)
    if sat:
        p=0
        for i,b in enumerate(pbits):
            val=b if isinstance(b,bool) else sol[abs(b)] ^ (b<0)
            if val: p|=1<<i
        q=0
        for j,b in enumerate(qbits):
            val=b if isinstance(b,bool) else sol[abs(b)] ^ (b<0)
            if val: q|=1<<j
        print('p',hex(p)); print('q',hex(q)); print('check',p*q==N,N%p==0,(p&mask)==leak)
        if p*q==N: print('plaintext',decrypt_from_p(p))
        return p,q
    return None,None

# helpers for q low/high guess
def q_low_for_low(low):
    pmod=(leak | (low<<150)) & ((1<<265)-1)
    qmod=(N * pow(pmod,-1,1<<265)) % (1<<265)
    return [(qmod>>j)&1 for j in range(265)]

def q_high_common_for_high(high):
    # compute interval of p for fixed high nibble and mask constraints; min/max
    # lower unknown bits fill 0/1. pmin = leak with high set, pmax = pmin | unknown mask below except low? include all unknown bits except high fixed.
    C=leak | (high<<920)
    unk=((1<<1024)-1) & ~mask
    # clear high nibble unknown from unk
    unk &= ~(((1<<4)-1)<<920)
    pmin=C
    pmax=C | unk
    qmin=N//pmax
    qmax=N//pmin
    # common high bits of all q in [qmin,qmax]
    d=qmin^qmax
    if d==0:
        pref=1024
    else:
        pref=1024-d.bit_length()
    res={}
    for j in range(1024-pref,1024):
        res[j]=(qmin>>j)&1
    return res,pref,qmin,qmax

if __name__=='__main__':
    import argparse
    ap=argparse.ArgumentParser()
    ap.add_argument('--low',type=int)
    ap.add_argument('--high',type=int)
    ap.add_argument('--qlow',action='store_true')
    ap.add_argument('--qhigh',action='store_true')
    args=ap.parse_args()
    qlow=q_low_for_low(args.low) if args.qlow and args.low is not None else None
    qhp=None
    if args.qhigh and args.high is not None:
        qhp,pref,_,_=q_high_common_for_high(args.high); print('q high pref',pref,flush=True)
    factor_with_sat(args.low,args.high,qlow,qhp)
