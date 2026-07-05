import pycryptosat, time, signal, sys

class Timeout(Exception): pass

def alarm_handler(signum, frame): raise Timeout()

class SatBuilder:
    def __init__(self):
        self.s = pycryptosat.Solver(time_limit=30, confl_limit=100000)
        self.nv = 0
        self.clauses = 0
        self.xors = 0
    def new(self):
        self.nv += 1
        return self.nv
    def clause(self, lits):
        # lits may contain bools or ints
        out=[]; seen=set()
        for l in lits:
            if l is True: return
            if l is False: continue
            if -l in seen: return
            if l not in seen: seen.add(l); out.append(l)
        if not out:
            raise ValueError('empty clause')
        self.s.add_clause(out); self.clauses += 1
    def unit(self, lit, val=True):
        if lit is True:
            if not val: raise ValueError('false const')
            return
        if lit is False:
            if val: raise ValueError('false const')
            return
        self.clause([lit if val else -lit])
    def xor_bit(self, args, const=False):
        vars=[]
        c=const
        for a in args:
            if a is False: continue
            if a is True: c = not c
            else: vars.append(a)
        if not vars: return c
        z=self.new()
        # z XOR vars = c
        self.s.add_xor_clause([z]+vars, c); self.xors += 1
        return z
    def and2(self,a,b):
        if a is False or b is False: return False
        if a is True: return b
        if b is True: return a
        z=self.new()
        self.clause([-z,a]); self.clause([-z,b]); self.clause([z,-a,-b])
        return z
    def or2(self,a,b):
        if a is True or b is True: return True
        if a is False: return b
        if b is False: return a
        z=self.new()
        self.clause([z,-a]); self.clause([z,-b]); self.clause([-z,a,b])
        return z
    def majority3(self,a,b,c):
        arr=[a,b,c]
        t=sum(1 for x in arr if x is True)
        rem=[x for x in arr if x is not True and x is not False]
        # Need at least 2 true among t + rem
        need=2-t
        if need <= 0: return True
        if need > len(rem): return False
        if need == 1:
            # OR of rem
            z=rem[0]
            for x in rem[1:]: z=self.or2(z,x)
            return z
        if need == len(rem):
            z=rem[0]
            for x in rem[1:]: z=self.and2(z,x)
            return z
        # only case len=3 need=2
        a,b,c=rem
        z=self.new()
        self.clause([-a,-b,z]); self.clause([-a,-c,z]); self.clause([-b,-c,z])
        self.clause([a,b,-z]); self.clause([a,c,-z]); self.clause([b,c,-z])
        return z
    def full_adder(self,a,b,c):
        s=self.xor_bit([a,b,c], False)
        carry=self.majority3(a,b,c)
        return s,carry
    def constrain_eq_const(self,a,bit):
        self.unit(a, bool(bit))

def build_factor_candidate(N, MASK, LEAK, lo, hi, qprefix=True):
    B=SatBuilder()
    nbits=1024
    def W(s,l): return ((1<<l)-1)<<s
    base=LEAK | (lo<<150) | (hi<<920)
    # p range and q low/prefix
    p_low=base & ((1<<265)-1)
    q_low=(N*pow(p_low,-1,1<<265))%(1<<265)
    mask_unknown = W(265,84)|W(362,58)|W(600,69)|W(682,87)|W(784,46)
    p_min=base & ~mask_unknown
    p_max=p_min | mask_unknown
    q_min=N//p_max; q_max=N//p_min
    xor=q_min^q_max; pref=1024-xor.bit_length() if xor else 1024
    q_pref=q_min>>(1024-pref) if pref>0 else 0
    # variables/constants for p,q bits
    pbits=[]; qbits=[]; pvars=[]; qvars=[]
    for i in range(nbits):
        if ((MASK>>i)&1) or (150<=i<154) or (920<=i<924):
            pbits.append(bool((base>>i)&1))
        else:
            v=B.new(); pbits.append(v); pvars.append((i,v))
    for j in range(nbits):
        if j<265:
            qbits.append(bool((q_low>>j)&1))
        elif qprefix and j>=1024-pref:
            qbits.append(bool((q_pref>>(j-(1024-pref)))&1))
        else:
            v=B.new(); qbits.append(v); qvars.append((j,v))
    # bit heap
    maxcols=2060
    cols=[[] for _ in range(maxcols)]
    const=[0]*maxcols
    for i,pi in enumerate(pbits):
        if pi is False: continue
        for j,qj in enumerate(qbits):
            col=i+j
            if col>=maxcols: raise ValueError('maxcol')
            if pi is True:
                term=qj
            elif qj is True:
                term=pi
            elif qj is False:
                continue
            else:
                term=B.and2(pi,qj)
            if term is True: const[col]+=1
            elif term is False: pass
            else: cols[col].append(term)
    # compress columns up to maxcols-2
    for k in range(maxcols-1):
        # fold true constants pairs
        if const[k]>=2:
            const[k+1]+=const[k]//2
            const[k]%=2
        arr=cols[k]
        if const[k]==1:
            arr.append(True); const[k]=0
        while len(arr)>2:
            a=arr.pop(); b=arr.pop(); c=arr.pop()
            s,car=B.full_adder(a,b,c)
            if s is True: arr.append(True)  # rare; will be handled? could make len high
            elif s is not False: arr.append(s)
            if car is True: const[k+1]+=1
            elif car is not False: cols[k+1].append(car)
            if len(arr)>2 and arr[-1] is True:
                # normalize constants in arr occasionally
                tr=sum(1 for x in arr if x is True)
                if tr>=2:
                    arr[:]=[x for x in arr if x is not True]
                    const[k+1]+=tr//2
                    if tr%2: arr.append(True)
        cols[k]=arr
    # final ripple add two remaining rows
    carry=False
    for k in range(maxcols-1):
        arr=cols[k]
        if const[k]>=2:
            const[k+1]+=const[k]//2; const[k]%=2
        if const[k]==1: arr=arr+[True]
        # ensure <=2; if more, compress (shouldn't)
        while len(arr)>2:
            a=arr.pop(); b=arr.pop(); c=arr.pop(); s,car=B.full_adder(a,b,c)
            if s is not False: arr.append(s)
            if car is True: const[k+1]+=1
            elif car is not False: cols[k+1].append(car)
        a=arr[0] if len(arr)>0 else False
        b=arr[1] if len(arr)>1 else False
        sm,carry=B.full_adder(a,b,carry)
        bit=(N>>k)&1 if k<2048 else 0
        B.constrain_eq_const(sm, bit)
    B.constrain_eq_const(carry, False)
    meta={'pvars':pvars,'qvars':qvars,'pref':pref,'q_low':q_low,'p_min':p_min,'p_max':p_max,'base':base}
    return B,meta

def solve_candidate(N,MASK,LEAK,lo,hi,timeout=120):
    print('build',lo,hi,flush=True); st=time.time(); B,meta=build_factor_candidate(N,MASK,LEAK,lo,hi); print('built vars',B.nv,'clauses',B.clauses,'xors',B.xors,'pref',meta['pref'],'time',time.time()-st,flush=True)
    signal.signal(signal.SIGALRM, alarm_handler); signal.alarm(timeout)
    try:
        st=time.time(); sat, sol = B.s.solve(time_limit=timeout, confl_limit=100000); signal.alarm(0); print('solve res',sat,'time',time.time()-st,flush=True)
    except Timeout:
        print('timeout solve',flush=True); return None
    if not sat: return False
    p=meta['base']
    for i,v in meta['pvars']:
        if sol[v]: p |= (1<<i)
        else: p &= ~(1<<i)
    return p
