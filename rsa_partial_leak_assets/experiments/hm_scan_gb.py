import sys,time,itertools, sympy as sp
from hm_test import hm_polys, N, MASK, LEAK, A, BB, X, Y
x,y=sp.symbols('x y')
primes=[2147483647,2147483629,2147483587,2147483579,2147483563,2147483549,2147483543,2147483497,2147483489,2147483477]

def poly_expr(poly,mod):
    return sum((c%mod)*x**a*y**b for (a,b),c in poly.items())

def gb_roots_mod(polys, mod, n_polys=6):
    exprs=[poly_expr(p,mod) for p in polys[:n_polys]]
    try:
        G=sp.groebner(exprs, y,x, modulus=mod, order='lex')
    except Exception as e:
        return None, f'gb_err {e}'
    if len(G.polys)==1 and G.polys[0].is_ground:
        c=int(G.polys[0].as_expr())%mod
        if c !=0:
            return [], 'unit'
    # Extract simple lex form. Print repr info
    roots=[]
    # find univariate in x and y relation
    uni_x=None; rel_y=None; uni_y=None; rel_x=None
    for g in G.polys:
        vars=g.as_expr().free_symbols
        if vars <= {x} and g.degree(x)>0:
            if uni_x is None or g.degree(x)<uni_x.degree(x): uni_x=sp.Poly(g.as_expr(), x, modulus=mod)
        if vars <= {y} and g.degree(y)>0:
            if uni_y is None or g.degree(y)<uni_y.degree(y): uni_y=sp.Poly(g.as_expr(), y, modulus=mod)
        # linear in y maybe y + h(x)
        P=sp.Poly(g.as_expr(), y,x, modulus=mod)
        if P.degree(y)==1:
            rel_y=g
        P2=sp.Poly(g.as_expr(), x,y, modulus=mod)
        if P2.degree(x)==1:
            rel_x=g
    if uni_x is not None:
        try:
            fl=sp.factor_list(uni_x, modulus=mod)[1]
            xs=[]
            for fac,exp in fl:
                fac=sp.Poly(fac,x,modulus=mod)
                if fac.degree()==1:
                    a=int(fac.nth(1))%mod; b=int(fac.nth(0))%mod
                    xs.append((-b*pow(a,-1,mod))%mod)
            for xv in xs:
                # solve for y using gcd of first equations substituted
                gg=None
                for expr in exprs[:min(4,len(exprs))]:
                    Py=sp.Poly(expr.subs(x,xv), y, modulus=mod)
                    if gg is None: gg=Py
                    else: gg=sp.gcd(gg,Py)
                if gg is not None and gg.degree()==1:
                    a=int(gg.nth(1))%mod; b=int(gg.nth(0))%mod
                    roots.append((xv,(-b*pow(a,-1,mod))%mod))
        except Exception as e:
            return None, 'root_err '+repr(e)
    elif uni_y is not None:
        # similarly
        pass
    return list(set(roots)), 'gb '+str([g.total_degree() for g in G.polys])

def crt_pair(a1,m1,a2,m2):
    t=((a2-a1)%m2)*pow(m1%m2,-1,m2)%m2
    return a1+m1*t,m1*m2

def try_extract(polys,p0,verbose=False):
    residues=[(0,0,1)]
    for mod in primes:
        rs,info=gb_roots_mod(polys,mod,6)
        if verbose: print('mod',mod,'rs',None if rs is None else len(rs),info)
        if rs is None or not rs or len(rs)>20:
            return None
        new=[]
        for xr,yr,M in residues:
            for a,b in rs:
                nx,nM=crt_pair(xr,M,a,mod); ny,_=crt_pair(yr,M,b,mod)
                # don't prune by X/Y until M > bounds? Actually if M small, residue can be >X but root< X with residue? no root = residue + kM; residue can be >X and k>=0 => >X. so prune.
                if nx < X and ny < Y:
                    new.append((nx,ny,nM))
        residues=new
        if not residues: return None
        if residues[0][2] > Y:
            break
    for xr,yr,M in residues:
        p=p0+A*xr+BB*yr
        if N%p==0 and (p&MASK)==LEAK:
            return p,xr,yr
    return None

if __name__=='__main__':
    m=int(sys.argv[1]) if len(sys.argv)>1 else 6
    t=int(sys.argv[2]) if len(sys.argv)>2 else 2
    start=int(sys.argv[3]) if len(sys.argv)>3 else 0
    end=int(sys.argv[4]) if len(sys.argv)>4 else 256
    t0=time.time(); survivors=[]
    for idx in range(start,end):
        lo=idx%16; hi=idx//16
        print(f'[{idx}/256] lo={lo:x} hi={hi:x}', flush=True)
        polys,p0=hm_polys(lo,hi,m,t,8)
        rs,info=gb_roots_mod(polys,2147483647,6)
        print(' firstmod', None if rs is None else len(rs), info, flush=True)
        if rs:
            r=try_extract(polys,p0,True)
            print(' extract',r,flush=True)
            if r:
                print('FOUND',r[0],hex(r[0])); sys.exit(0)
            survivors.append((lo,hi,len(rs),info))
    print('done time',time.time()-t0,'survivors',survivors)
