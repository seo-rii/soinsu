import sympy as sp
x,y=sp.symbols('x y')

def expr(poly): return sum(int(c)*(x**i)*(y**j) for (i,j),c in poly.items())

def eval_mod(poly,xx,yy,p): return sum((int(c)%p)*pow(xx,i,p)*pow(yy,j,p) for (i,j),c in poly.items())%p

def roots_mod_groebner(polys,p,maxK=8):
    exprs=[]
    for poly in polys[:maxK]:
        e=expr(poly)
        if e!=0:
            exprs.append(e)
    if len(exprs)<2: return []
    # Try growing prefixes
    for K in range(2,len(exprs)+1):
        try:
            G=sp.groebner(exprs[:K], x,y, order='lex', modulus=p)
        except Exception as e:
            continue
        gb=[sp.Poly(g,x,y,modulus=p) for g in G.polys if not g.is_zero]
        if not gb: continue
        # Look for univariate y polynomial
        ypolys=[]
        for g in gb:
            if all(mon[0]==0 for mon in g.monoms()):
                yy=sp.Poly(g.as_expr(), y, modulus=p)
                if yy.degree()>0: ypolys.append(yy)
        if not ypolys: continue
        yy=min(ypolys, key=lambda P:P.degree())
        try:
            fl=sp.factor_list(yy, modulus=p)[1]
        except Exception:
            continue
        yroots=[]
        for fac,ex in fl:
            fac=sp.Poly(fac,y,modulus=p)
            if fac.degree()==1:
                a,b=map(lambda z:int(z)%p, fac.all_coeffs())
                if a: yroots.append((-b*pow(a,-1,p))%p)
        pairs=[]
        for yr in yroots:
            # find x via a polynomial linear in x after y substitution
            xroots=set()
            for g in gb:
                gx=sp.Poly(g.as_expr().subs(y,yr), x, modulus=p)
                if gx.is_zero or gx.degree()<=0: continue
                try: flx=sp.factor_list(gx, modulus=p)[1]
                except Exception: continue
                for fac,ex in flx:
                    fac=sp.Poly(fac,x,modulus=p)
                    if fac.degree()==1:
                        a,b=map(lambda z:int(z)%p, fac.all_coeffs())
                        if a: xroots.add((-b*pow(a,-1,p))%p)
            for xr in xroots:
                if all(eval_mod(poly,xr,yr,p)==0 for poly in polys[:maxK]):
                    pairs.append((xr,yr))
        if pairs:
            return sorted(set(pairs))
    return []

def crt_pair(a1,m1,a2,m2):
    inv=pow(m1,-1,m2); x1,y1=a1; x2,y2=a2
    return ((x1+m1*(((x2-x1)%m2)*inv%m2))%(m1*m2), (y1+m1*(((y2-y1)%m2)*inv%m2))%(m1*m2)), m1*m2

def recover(polys,X,Y,primes,verbose=False):
    cands=[]; mod=1
    for p in primes:
        pairs=roots_mod_groebner(polys,p)
        if verbose: print('p',p,'pairs',len(pairs),pairs[:3],flush=True)
        if not pairs: continue
        if not cands:
            cands=[(pair,p) for pair in pairs]
        else:
            new=[]
            for a,m in cands:
                for pair in pairs:
                    new.append(crt_pair(a,m,pair,p))
            # maybe many; filter by lower bound? Keep all up to 500
            seen=set(); c=[]
            for a,m in new:
                key=(a[0],a[1],m)
                if key not in seen:
                    seen.add(key); c.append((a,m))
            cands=c[:500]
        if verbose and cands: print(' modbits',cands[0][1].bit_length(),'cands',len(cands),flush=True)
        if cands and cands[0][1] > max(X,Y)*2:
            outs=[]
            for (xx,yy),m in cands:
                if xx<X and yy<Y: outs.append((xx,yy))
            if outs: return outs
    return []
