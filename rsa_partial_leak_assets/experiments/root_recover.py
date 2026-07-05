import sympy as sp
from math import prod

x_sym,y_sym = sp.symbols('x y')

def poly_to_expr(p):
    return sum(int(c)*(x_sym**i)*(y_sym**j) for (i,j),c in p.items())

def roots_mod_pair(polys, p, max_pairs=200):
    # polys list dict; use first two nonzero polynomials by default
    exprs=[poly_to_expr(g) for g in polys if g]
    # try multiple pairs until roots found
    for a in range(min(len(exprs),5)):
      for b in range(a+1,min(len(exprs),8)):
        f=sp.Poly(exprs[a], x_sym, y_sym, modulus=p)
        g=sp.Poly(exprs[b], x_sym, y_sym, modulus=p)
        if f.is_zero or g.is_zero: continue
        try:
            R=sp.resultant(f.as_expr(), g.as_expr(), x_sym)
            R=sp.Poly(R, y_sym, modulus=p)
        except Exception as e:
            continue
        if R.is_zero or R.degree()<=0: continue
        try:
            fl=sp.factor_list(R, modulus=p)[1]
        except Exception as e:
            continue
        yroots=[]
        for fac,exp in fl:
            fac=sp.Poly(fac, y_sym, modulus=p)
            if fac.degree()==1:
                coeffs=fac.all_coeffs() # a*y+b
                aa=int(coeffs[0])%p; bb=int(coeffs[1])%p
                if aa:
                    yroots.append((-bb*pow(aa,-1,p))%p)
        pairs=[]
        for yy in yroots:
            fx=sp.Poly(f.as_expr().subs(y_sym, yy), x_sym, modulus=p)
            if fx.is_zero:
                # try g instead
                fx=sp.Poly(g.as_expr().subs(y_sym, yy), x_sym, modulus=p)
            if fx.is_zero: continue
            try:
                flx=sp.factor_list(fx, modulus=p)[1]
            except Exception:
                continue
            for fac,exp in flx:
                fac=sp.Poly(fac,x_sym,modulus=p)
                if fac.degree()==1:
                    aa=int(fac.all_coeffs()[0])%p; bb=int(fac.all_coeffs()[1])%p
                    xx=(-bb*pow(aa,-1,p))%p
                    ok=True
                    for EE in polys[:min(len(polys),8)]:
                        if sum((int(cc)%p)*pow(xx,ii,p)*pow(yy,jj,p) for (ii,jj),cc in EE.items()) % p != 0:
                            ok=False; break
                    if ok:
                        pairs.append((xx,yy))
                        if len(pairs)>max_pairs: return pairs
        if pairs:
            # unique
            return sorted(set(pairs))
    return []

def crt_pair(a1,m1,a2,m2):
    x1,y1=a1; x2,y2=a2
    inv=pow(m1,-1,m2)
    x=x1 + m1*(((x2-x1)%m2)*inv % m2)
    y=y1 + m1*(((y2-y1)%m2)*inv % m2)
    return (x%(m1*m2), y%(m1*m2)), m1*m2

def recover_crt(polys, X, Y, primes, verbose=False):
    cands=[((0,0),1)]
    for p in primes:
        pairs=roots_mod_pair(polys,p)
        if verbose: print('prime',p,'pairs',len(pairs),pairs[:4],flush=True)
        if not pairs: continue
        new=[]
        for a,m in cands:
            if m==1:
                for pair in pairs: new.append((pair,p))
            else:
                for pair in pairs:
                    try: new.append(crt_pair(a,m,pair,p))
                    except ValueError: pass
        # reduce duplicates
        seen=set(); cands2=[]
        for a,m in new:
            key=(a[0]%m,a[1]%m,m)
            if key not in seen:
                seen.add(key); cands2.append((a,m))
        cands=cands2[:2000]
        if verbose: print(' combined',len(cands),'modbits',cands[0][1].bit_length() if cands else 0,flush=True)
        if cands and cands[0][1] > max(X,Y)*2:
            # return candidates within bounds (min repr)
            outs=[]
            for (xx,yy),m in cands:
                # roots expected 0<= <bounds, use residue direct; also try negative not needed
                if xx<X and yy<Y: outs.append((xx,yy))
            if outs: return outs
    return []

if __name__=='__main__':
    pass
