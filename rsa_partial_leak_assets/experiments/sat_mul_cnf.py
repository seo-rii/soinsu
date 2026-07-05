from pysat.solvers import Solver
import itertools,time,sys

class SatBuilderCNF:
    def __init__(self, solver_name='cadical153'):
        self.s=Solver(name=solver_name)
        self.nv=0; self.clauses=0
    def new(self): self.nv+=1; return self.nv
    def clause(self,lits):
        out=[]; seen=set()
        for l in lits:
            if l is True: return
            if l is False: continue
            if -l in seen: return
            if l not in seen: seen.add(l); out.append(l)
        if not out: raise ValueError('empty')
        self.s.add_clause(out); self.clauses+=1
    def unit(self,lit,val=True):
        if lit is True:
            if not val: raise ValueError('bad const')
            return
        if lit is False:
            if val: raise ValueError('bad const')
            return
        self.clause([lit if val else -lit])
    def xor_bit(self,args,const=False):
        vars=[]; c=const
        for a in args:
            if a is False: continue
            if a is True: c=not c
            else: vars.append(a)
        if not vars: return c
        z=self.new()
        allv=[z]+vars
        # z xor vars = c. For each assignment with parity != c, forbid.
        n=len(allv)
        for bits in itertools.product([False,True], repeat=n):
            if (sum(bits)%2==1) != c:
                # forbid this assignment
                self.clause([(-v if bit else v) for v,bit in zip(allv,bits)])
        return z
    def and2(self,a,b):
        if a is False or b is False: return False
        if a is True: return b
        if b is True: return a
        z=self.new(); self.clause([-z,a]); self.clause([-z,b]); self.clause([z,-a,-b]); return z
    def or2(self,a,b):
        if a is True or b is True: return True
        if a is False: return b
        if b is False: return a
        z=self.new(); self.clause([z,-a]); self.clause([z,-b]); self.clause([-z,a,b]); return z
    def majority3(self,a,b,c):
        arr=[a,b,c]; t=sum(1 for x in arr if x is True); rem=[x for x in arr if x is not True and x is not False]
        need=2-t
        if need<=0: return True
        if need>len(rem): return False
        if need==1:
            z=rem[0]
            for x in rem[1:]: z=self.or2(z,x)
            return z
        if need==len(rem):
            z=rem[0]
            for x in rem[1:]: z=self.and2(z,x)
            return z
        a,b,c=rem; z=self.new()
        self.clause([-a,-b,z]); self.clause([-a,-c,z]); self.clause([-b,-c,z]); self.clause([a,b,-z]); self.clause([a,c,-z]); self.clause([b,c,-z]); return z
    def full_adder(self,a,b,c): return self.xor_bit([a,b,c]), self.majority3(a,b,c)
    def eqconst(self,a,bit): self.unit(a,bool(bit))

def build(N,MASK,LEAK,lo,hi,solver='cadical153'):
    B=SatBuilderCNF(solver)
    def W(s,l): return ((1<<l)-1)<<s
    base=LEAK|(lo<<150)|(hi<<920); p_low=base&((1<<265)-1); q_low=(N*pow(p_low,-1,1<<265))%(1<<265)
    mask_unknown=W(265,84)|W(362,58)|W(600,69)|W(682,87)|W(784,46); p_min=base&~mask_unknown; p_max=p_min|mask_unknown
    q_min=N//p_max; q_max=N//p_min; xor=q_min^q_max; pref=1024-xor.bit_length() if xor else 1024; q_pref=q_min>>(1024-pref) if pref else 0
    pbits=[]; qbits=[]; pvars=[]
    for i in range(1024):
        if ((MASK>>i)&1) or 150<=i<154 or 920<=i<924: pbits.append(bool((base>>i)&1))
        else:
            v=B.new(); pbits.append(v); pvars.append((i,v))
    for j in range(1024):
        if j<265: qbits.append(bool((q_low>>j)&1))
        elif j>=1024-pref: qbits.append(bool((q_pref>>(j-(1024-pref)))&1))
        else: qbits.append(B.new())
    maxcols=2060; cols=[[] for _ in range(maxcols)]; const=[0]*maxcols
    for i,pi in enumerate(pbits):
        if pi is False: continue
        for j,qj in enumerate(qbits):
            col=i+j
            if pi is True: term=qj
            elif qj is True: term=pi
            elif qj is False: continue
            else: term=B.and2(pi,qj)
            if term is True: const[col]+=1
            elif term is not False: cols[col].append(term)
    for k in range(maxcols-1):
        if const[k]>=2: const[k+1]+=const[k]//2; const[k]%=2
        arr=cols[k]
        if const[k]==1: arr.append(True); const[k]=0
        while len(arr)>2:
            a=arr.pop(); b=arr.pop(); c=arr.pop(); sm,car=B.full_adder(a,b,c)
            if sm is True: arr.append(True)
            elif sm is not False: arr.append(sm)
            if car is True: const[k+1]+=1
            elif car is not False: cols[k+1].append(car)
            tr=sum(1 for x in arr if x is True)
            if tr>=2:
                arr[:]=[x for x in arr if x is not True]; const[k+1]+=tr//2
                if tr&1: arr.append(True)
        cols[k]=arr
    carry=False
    for k in range(maxcols-1):
        if const[k]>=2: const[k+1]+=const[k]//2; const[k]%=2
        arr=cols[k]
        if const[k]==1: arr=arr+[True]
        while len(arr)>2:
            a=arr.pop(); b=arr.pop(); c=arr.pop(); sm,car=B.full_adder(a,b,c)
            if sm is not False: arr.append(sm)
            if car is True: const[k+1]+=1
            elif car is not False: cols[k+1].append(car)
        a=arr[0] if len(arr)>0 else False; b=arr[1] if len(arr)>1 else False
        sm,carry=B.full_adder(a,b,carry)
        B.eqconst(sm, (N>>k)&1 if k<2048 else 0)
    B.eqconst(carry,False)
    return B, {'base':base,'pvars':pvars,'pref':pref}

def solve(N,MASK,LEAK,lo,hi,solver='cadical153',conf_budget=100000):
    st=time.time(); B,meta=build(N,MASK,LEAK,lo,hi,solver); print('built',B.nv,B.clauses,'pref',meta['pref'],'time',time.time()-st,flush=True)
    B.s.conf_budget(conf_budget)
    st=time.time(); r=B.s.solve_limited(expect_interrupt=True); print('solve',r,'time',time.time()-st,flush=True)
    if r is True:
        model=B.s.get_model(); vals=[False]*(B.nv+1)
        for lit in model:
            if lit>0: vals[lit]=True
        p=meta['base']
        for i,v in meta['pvars']:
            if vals[v]: p|=1<<i
            else: p&=~(1<<i)
        return p
    return r
