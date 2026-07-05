import sympy as sp, time, sys, signal
from hm_core import *
x,y=sp.symbols('x y')
MOD=2147483647

def timeout_handler(signum, frame): raise TimeoutError
signal.signal(signal.SIGALRM, timeout_handler)

def expr(poly,mod=MOD): return sum((c%mod)*x**a*y**b for (a,b),c in poly.items())

def gb_first(polys, n=5, secs=3, mod=MOD):
    signal.alarm(secs)
    try:
        G=sp.groebner([expr(p,mod) for p in polys[:n]], y,x, modulus=mod, order='lex')
        signal.alarm(0)
    except TimeoutError:
        signal.alarm(0); return None,'timeout'
    except Exception as e:
        signal.alarm(0); return None,'err '+repr(e)[:80]
    if len(G.polys)==1 and G.polys[0].is_ground and int(G.polys[0].as_expr())%mod:
        return [],'unit'
    return G,'nonunit '+str([g.total_degree() for g in G.polys])

if __name__=='__main__':
    m=int(sys.argv[1]) if len(sys.argv)>1 else 6; t=int(sys.argv[2]) if len(sys.argv)>2 else 2
    start=int(sys.argv[3]) if len(sys.argv)>3 else 0; end=int(sys.argv[4]) if len(sys.argv)>4 else 256
    survivors=[]; t0=time.time()
    for idx in range(start,end):
        lo=idx%16; hi=idx//16
        try:
            polys,p0,lt,stats=hm_polys(lo,hi,m,t,8)
            G,info=gb_first(polys,5,3)
        except Exception as e:
            info='fail '+repr(e); G=None; lt=-1
        print(f'{idx:03d} lo={lo:x} hi={hi:x} LLL={lt:.2f} {info}', flush=True)
        if G not in (None,[]): survivors.append((lo,hi,info))
    print('survivors',survivors,'time',time.time()-t0)
