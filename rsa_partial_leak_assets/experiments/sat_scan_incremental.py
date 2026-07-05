import sys,time
sys.path.append('/mnt/data')
from sat_factor import SatMul, q_low_for_low, q_high_common_for_high
from rsa_partial_lib import N,mask,leak,ct,e,decrypt_from_p

def build_base():
    sm=SatMul(); nbits=1024
    pbits=[]
    for i in range(nbits):
        if (mask>>i)&1: pbits.append(bool((leak>>i)&1))
        else:
            v=sm.new(); pbits.append(v)
    qbits=[sm.new() for _ in range(nbits)]
    sm.eq_const(qbits[0], True); sm.eq_const(qbits[1023], True)
    maxcols=2100
    cols=[[] for _ in range(maxcols)]
    t0=time.time(); terms=0
    for i,pb in enumerate(pbits):
        if pb is False: continue
        for j,qb in enumerate(qbits):
            term = qb if pb is True else sm.and_bit(pb,qb)
            if term is not False:
                cols[i+j].append(term); terms+=1
        if i%256==0: print('partial',i,'vars',sm.var,'clauses',sm.clauses,flush=True)
    print('initial terms',terms,'vars',sm.var,'clauses',sm.clauses,'time',time.time()-t0,flush=True)
    for k in range(maxcols-1):
        col=cols[k]
        while len(col)>=3:
            a=col.pop(); b=col.pop(); c=col.pop()
            s,carry=sm.full_adder(a,b,c)
            if s is not False: col.append(s)
            if carry is not False: cols[k+1].append(carry)
        if len(col)==2:
            a=col.pop(); b=col.pop(); s,carry=sm.half_adder(a,b)
            if s is not False: col.append(s)
            if carry is not False: cols[k+1].append(carry)
        Nbit=bool((N>>k)&1) if k<2048 else False
        if len(col)==0:
            if Nbit: raise RuntimeError('empty col for 1')
        elif len(col)==1: sm.eq_const(col[0],Nbit)
        if k%256==0: print('reduce',k,'vars',sm.var,'clauses',sm.clauses,'xors',sm.xors,'next',len(cols[k+1]),flush=True)
    for b in cols[-1]: sm.eq_const(b,False)
    print('built vars',sm.var,'clauses',sm.clauses,'xors',sm.xors,flush=True)
    return sm,pbits,qbits

def lit_for_bit(bit,val):
    # bit is bool or literal var; return assumption literal or None if already satisfied; raise contradiction
    if isinstance(bit,bool):
        if bit!=bool(val): raise ValueError('assumption contradicts const')
        return None
    return bit if val else -bit

def assumptions_for(low,high,pbits,qbits):
    ass=[]
    for b in range(4):
        lit=lit_for_bit(pbits[150+b], (low>>b)&1)
        if lit is not None: ass.append(lit)
    for b in range(4):
        lit=lit_for_bit(pbits[920+b], (high>>b)&1)
        if lit is not None: ass.append(lit)
    ql=q_low_for_low(low)
    for j,val in enumerate(ql):
        lit=lit_for_bit(qbits[j],val)
        if lit is not None: ass.append(lit)
    qhp,pref,_,_=q_high_common_for_high(high)
    for j,val in qhp.items():
        lit=lit_for_bit(qbits[j],val)
        if lit is not None: ass.append(lit)
    return ass,pref

def extract_solution(sol,pbits,qbits):
    p=0; q=0
    for i,b in enumerate(pbits):
        val=b if isinstance(b,bool) else (sol[abs(b)] ^ (b<0))
        if val: p|=1<<i
    for j,b in enumerate(qbits):
        val=b if isinstance(b,bool) else (sol[abs(b)] ^ (b<0))
        if val: q|=1<<j
    return p,q

if __name__=='__main__':
    sm,pbits,qbits=build_base()
    order=[]
    # try all; optionally skip high 7 etc no
    for high in range(16):
      for low in range(16): order.append((low,high))
    for idx,(low,high) in enumerate(order):
        ass,pref=assumptions_for(low,high,pbits,qbits)
        t=time.time()
        sat,sol=sm.s.solve(assumptions=ass,time_limit=5.0,confl_limit=200000)
        dt=time.time()-t
        print('TRY',idx,'low high',low,high,'pref',pref,'sat',sat,'time',dt,flush=True)
        if sat is True:
            p,q=extract_solution(sol,pbits,qbits)
            print('p',hex(p)); print('q',hex(q)); print('check',p*q==N,N%p==0,(p&mask)==leak)
            if p*q==N:
                print('PLAINTEXT',decrypt_from_p(p)); break
        elif sat is None:
            # retry with more time maybe
            print('timeout candidate',low,high,flush=True)
