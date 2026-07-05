import sympy as sp
x_sym,y_sym=sp.symbols('x y')

def poly_to_expr(p):
    return sum(int(c)*(x_sym**i)*(y_sym**j) for (i,j),c in p.items())

def eval_poly_mod_dict(poly,xx,yy,p):
    s=0
    # precompute powers maybe not needed
    for (i,j),c in poly.items():
        s=(s + (int(c)%p)*pow(xx,i,p)*pow(yy,j,p))%p
    return s

def roots_univar_bf(P,p):
    # P sympy Poly in one variable over modulus p
    coeffs=[int(c)%p for c in P.all_coeffs()]
    roots=[]
    for a in range(p):
        v=0
        for c in coeffs:
            v=(v*a+c)%p
        if v==0:
            roots.append(a)
    return roots

def roots_mod_pair(polys,p,max_pairs=500,max_poly=10):
    exprs=[poly_to_expr(g) for g in polys if g]
    # Try pairs involving low-degree first; caller should sort.
    L=min(len(exprs), max_poly)
    for a in range(L):
      for b in range(a+1,L):
        f=sp.Poly(exprs[a], x_sym, y_sym, modulus=p)
        g=sp.Poly(exprs[b], x_sym, y_sym, modulus=p)
        if f.is_zero or g.is_zero: continue
        try:
            R=sp.Poly(sp.resultant(f.as_expr(),g.as_expr(),x_sym), y_sym, modulus=p)
        except Exception:
            continue
        if R.is_zero or R.degree()<=0: continue
        if R.degree()>80: continue
        yroots=roots_univar_bf(R,p)
        if not yroots: continue
        pairs=[]
        for yy in yroots:
            # brute x over p, use first polynomial after substitution
            fx=sp.Poly(f.as_expr().subs(y_sym,yy), x_sym, modulus=p)
            if fx.is_zero: fx=sp.Poly(g.as_expr().subs(y_sym,yy), x_sym, modulus=p)
            if fx.is_zero: continue
            xroots=roots_univar_bf(fx,p) if fx.degree()<=80 else []
            for xx in xroots:
                if all(eval_poly_mod_dict(poly,xx,yy,p)==0 for poly in polys[:max_poly]):
                    pairs.append((xx,yy))
                    if len(pairs)>max_pairs: return sorted(set(pairs))
        if pairs:
            return sorted(set(pairs))
    return []

def crt_pair(a1,m1,a2,m2):
    x1,y1=a1; x2,y2=a2
    inv=pow(m1,-1,m2)
    M=m1*m2
    return ((x1+m1*(((x2-x1)%m2)*inv%m2))%M, (y1+m1*(((y2-y1)%m2)*inv%m2))%M), M

def recover_crt(polys,X,Y,primes,verbose=False,max_poly=10):
    cands=[]
    for p in primes:
        pairs=roots_mod_pair(polys,p,max_poly=max_poly)
        if verbose: print('prime',p,'pairs',len(pairs),pairs[:3],flush=True)
        if not pairs: continue
        if not cands:
            cands=[(pair,p) for pair in pairs]
        else:
            new=[]
            for a,m in cands:
                for pair in pairs:
                    new.append(crt_pair(a,m,pair,p))
            # uniqueness + cap
            seen=set(); c=[]
            for a,m in new:
                key=(a[0],a[1],m)
                if key not in seen:
                    seen.add(key); c.append((a,m))
            cands=c[:5000]
        if verbose and cands: print(' combined',len(cands),'modbits',cands[0][1].bit_length(),flush=True)
        if cands and cands[0][1] > max(X,Y)*2:
            outs=[]
            for (xx,yy),m in cands:
                if xx<X and yy<Y: outs.append((xx,yy))
            if outs: return outs
    return []
