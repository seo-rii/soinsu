import sys,time,argparse
sys.path.append('/mnt/data')
from sat_factor import SatMul, q_low_for_low, q_high_common_for_high
from rsa_partial_lib import N,mask,leak,decrypt_from_p

def build_and_solve(low,high,solve_time=60,threads=1,verbose=0):
    # monkey patch Solver? SatMul creates Solver(threads?) no; recreate here by modifying after? Use pycryptosat directly
    from pycryptosat import Solver
    sm=SatMul(); sm.s=Solver(verbose=verbose, threads=threads)
    nbits=1024
    ql=q_low_for_low(low)
    qhp,pref,_,_=q_high_common_for_high(high)
    pbits=[]
    for i in range(nbits):
        if (mask>>i)&1:
            pbits.append(bool((leak>>i)&1))
        elif 150<=i<154:
            pbits.append(bool((low>>(i-150))&1))
        elif 920<=i<924:
            pbits.append(bool((high>>(i-920))&1))
        else:
            pbits.append(sm.new())
    qbits=[]
    for j in range(nbits):
        if j < len(ql):
            qbits.append(bool(ql[j]))
        elif j in qhp:
            qbits.append(bool(qhp[j]))
        elif j==0 or j==1023:
            qbits.append(True)
        else:
            qbits.append(sm.new())
    maxcols=2100
    cols=[[] for _ in range(maxcols)]
    t0=time.time(); terms=0
    for i,pb in enumerate(pbits):
        if pb is False: continue
        for j,qb in enumerate(qbits):
            term=qb if pb is True else sm.and_bit(pb,qb)
            if term is not False:
                cols[i+j].append(term); terms+=1
    print('initial',low,high,'pref',pref,'terms',terms,'vars',sm.var,'clauses',sm.clauses,'time',time.time()-t0,flush=True)
    for k in range(maxcols-1):
        col=cols[k]
        while len(col)>=3:
            a=col.pop(); b=col.pop(); c=col.pop(); s,carry=sm.full_adder(a,b,c)
            if s is not False: col.append(s)
            if carry is not False: cols[k+1].append(carry)
        if len(col)==2:
            a=col.pop(); b=col.pop(); s,carry=sm.half_adder(a,b)
            if s is not False: col.append(s)
            if carry is not False: cols[k+1].append(carry)
        Nbit=bool((N>>k)&1) if k<2048 else False
        if len(col)==0:
            if Nbit:
                print('contrad empty at',k); return None
        elif len(col)==1:
            sm.eq_const(col[0],Nbit)
    for b in cols[-1]: sm.eq_const(b,False)
    print('built',low,high,'vars',sm.var,'clauses',sm.clauses,'xors',sm.xors,flush=True)
    t=time.time(); sat,sol=sm.s.solve(time_limit=solve_time)
    print('solved',low,high,'sat',sat,'time',time.time()-t,flush=True)
    if sat is True:
        p=0; q=0
        for i,b in enumerate(pbits):
            val=b if isinstance(b,bool) else (sol[abs(b)] ^ (b<0))
            if val: p|=1<<i
        for j,b in enumerate(qbits):
            val=b if isinstance(b,bool) else (sol[abs(b)] ^ (b<0))
            if val: q|=1<<j
        print('p',hex(p)); print('q',hex(q)); print('check',p*q==N,N%p==0,(p&mask)==leak)
        if p*q==N: print('PLAINTEXT',decrypt_from_p(p))
        return p
    return None
if __name__=='__main__':
    ap=argparse.ArgumentParser(); ap.add_argument('low',type=int); ap.add_argument('high',type=int); ap.add_argument('--time',type=float,default=60); ap.add_argument('--threads',type=int,default=1)
    args=ap.parse_args(); build_and_solve(args.low,args.high,args.time,args.threads)
