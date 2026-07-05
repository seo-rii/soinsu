from hm_test import hm_polys, N, LEAK, MASK, A, BB, X, Y, WX, WY
import sympy as sp, itertools, time, math, sys
x,y=sp.symbols('x y')

def to_poly(poly,p, vars=(x,y)):
    expr=0
    for (a,b),c in poly.items():
        expr += (c % p) * (x**a) * (y**b)
    return sp.Poly(expr, x, y, modulus=p)

def roots_pair_mod(poly1, poly2, p):
    P=to_poly(poly1,p); Q=to_poly(poly2,p)
    try:
        R=sp.resultant(P.as_expr(), Q.as_expr(), y)
        Rx=sp.Poly(R, x, modulus=p)
    except Exception as e:
        #print('res err',e)
        return []
    if Rx.is_zero:
        return []
    try:
        fl=sp.factor_list(Rx, modulus=p)[1]
    except Exception as e:
        #print('factor err',e, 'deg', Rx.degree())
        return []
    xs=[]
    for fac, exp in fl:
        fac=sp.Poly(fac,x,modulus=p)
        if fac.degree()==1:
            a=int(fac.nth(1))%p; b=int(fac.nth(0))%p
            if a%p:
                xs.append((-b*pow(a,-1,p))%p)
    out=[]
    for xv in set(xs):
        # substitute x=xv, get gcd in y
        P1=sp.Poly(P.as_expr().subs(x,xv), y, modulus=p)
        Q1=sp.Poly(Q.as_expr().subs(x,xv), y, modulus=p)
        try:
            G=sp.gcd(P1,Q1)
        except Exception:
            continue
        if G.is_zero: continue
        if G.degree()==1:
            a=int(G.nth(1))%p; b=int(G.nth(0))%p
            if a%p:
                yv=(-b*pow(a,-1,p))%p
                out.append((xv,yv))
        elif G.degree()>1 and G.degree()<20:
            try:
                fl2=sp.factor_list(G, modulus=p)[1]
                for fac, exp in fl2:
                    fac=sp.Poly(fac,y,modulus=p)
                    if fac.degree()==1:
                        a=int(fac.nth(1))%p; b=int(fac.nth(0))%p
                        out.append((xv, (-b*pow(a,-1,p))%p))
            except Exception:
                pass
    return list(set(out))

def crt_pair(a1,m1,a2,m2):
    # x=a1 mod m1, a2 mod m2
    t=((a2-a1)%m2)*pow(m1 % m2,-1,m2)%m2
    return a1 + m1*t, m1*m2

def extract(polys,p0, verbose=True):
    primes=[2147483647,2147483629,2147483587,2147483579,2147483563,2147483549,2147483543,2147483497,2147483489,2147483477,2147483423]
    # try several pairs
    for i,j in itertools.combinations(range(min(len(polys),8)),2):
        residues=[(0,0,1)] # x,y,mod
        if verbose: print('pair',i,j)
        good=True
        for pr in primes:
            rs=roots_pair_mod(polys[i],polys[j],pr)
            if verbose: print(' prime',pr,'roots',len(rs))
            if not rs or len(rs)>50:
                good=False; break
            new=[]
            for xr,yr,mod in residues:
                for a,b in rs:
                    nx,nmod=crt_pair(xr,mod,a,pr)
                    ny,_=crt_pair(yr,mod,b,pr)
                    # keep representative under bounds if possible? no
                    if nx < X and ny < Y:
                        new.append((nx,ny,nmod))
                    else:
                        # also if nx mod nmod can have small representative? CRT product increasing; root positive <X,Y so residue should equal root once mod>X/Y; before can be > bound but actual root residue may be > current? If current mod < X, residue <mod<X so okay. if residue >X then no positive root <X with that residue if mod<X? maybe still root = residue + kmod, so if residue>X but mod<X could root >X; prune. same for Y.
                        pass
            residues=new
            if verbose: print('  combined',len(residues),'modbits',residues[0][2].bit_length() if residues else 0)
            if not residues:
                good=False; break
            if residues[0][2] > Y:
                for xr,yr,mod in residues:
                    if xr < X and yr < Y:
                        p=p0 + A*xr + BB*yr
                        if N%p==0 and (p&MASK)==LEAK:
                            return p,xr,yr,(i,j)
                # maybe continue no
        # final verify all
        for xr,yr,mod in residues:
            p=p0 + A*xr + BB*yr
            if N%p==0 and (p&MASK)==LEAK:
                return p,xr,yr,(i,j)
    return None

if __name__=='__main__':
    lo=int(sys.argv[1],16) if len(sys.argv)>1 else 0
    hi=int(sys.argv[2],16) if len(sys.argv)>2 else 0
    polys,p0=hm_polys(lo,hi,8,3,10)
    print('start extract')
    r=extract(polys,p0,True)
    print('RESULT',r)
