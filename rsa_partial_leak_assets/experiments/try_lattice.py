import sys, time, math, itertools, random
sys.path.append('/mnt/data')
from rsa_partial_lib import *
import sympy as sp
x_sym,y_sym=sp.symbols('x y')

# known masks within x,y for final filter
x_known=(leak>>x_start)&(X-1)
y_known=(leak>>y_start)&(Y-1)
x_mask=(mask>>x_start)&(X-1)
y_mask=(mask>>y_start)&(Y-1)

def coeffs_to_poly_mod(coeffs, mod):
    d={mon:int(c%mod) for mon,c in coeffs.items() if c%mod}
    if not d:
        return sp.Poly(0,x_sym,y_sym, modulus=mod)
    return sp.Poly.from_dict(d, (x_sym,y_sym), modulus=mod)

def roots_univar_mod(poly, var, mod):
    P=sp.Poly(poly, var, modulus=mod)
    if P.is_zero:
        return None  # all roots
    roots=[]
    try:
        fl=sp.factor_list(P, modulus=mod)[1]
    except Exception as e:
        # fallback ground roots
        try:
            gr=P.ground_roots()
            return [int(r)%mod for r,mult in gr.items()]
        except Exception:
            return []
    for fac,mult in fl:
        fac=sp.Poly(fac, var, modulus=mod)
        if fac.degree()==1:
            a=int(fac.nth(1))%mod; b=int(fac.nth(0))%mod
            if a:
                roots.append((-b*pow(a,-1,mod))%mod)
    return sorted(set(roots))

def roots_pair_mod(c1,c2,mod, max_roots=50):
    P1=coeffs_to_poly_mod(c1,mod); P2=coeffs_to_poly_mod(c2,mod)
    if P1.is_zero or P2.is_zero: return []
    try:
        res=sp.resultant(P1.as_expr(), P2.as_expr(), y_sym)
        R=sp.Poly(res, x_sym, modulus=mod)
    except Exception as e:
        return []
    xroots=roots_univar_mod(R,x_sym,mod)
    if xroots is None:
        return []
    out=[]
    for xr in xroots[:max_roots]:
        # substitute x=xr, get univar y polynomials
        py1=sp.Poly(P1.as_expr().subs(x_sym,xr), y_sym, modulus=mod)
        py2=sp.Poly(P2.as_expr().subs(x_sym,xr), y_sym, modulus=mod)
        if py1.is_zero and py2.is_zero:
            continue
        if py1.is_zero:
            G=py2
        elif py2.is_zero:
            G=py1
        else:
            try:
                G=sp.gcd(py1,py2)
            except Exception:
                continue
        yroots=roots_univar_mod(G,y_sym,mod)
        if yroots is None:
            continue
        for yr in yroots:
            # verify both polys mod
            if P1.eval(xr,yr)%mod==0 and P2.eval(xr,yr)%mod==0:
                out.append((xr,yr))
    return sorted(set(out))

def crt_pair(a,m,b,n):
    # coprime
    t=((b-a)%n)*pow(m % n,-1,n)%n
    return a + m*t, m*n

def recover_with_pair(c1,c2, primes, verbose=False):
    cand=[(0,0,1)]
    for pr in primes:
        roots=roots_pair_mod(c1,c2,pr)
        if verbose: print(' prime',pr,'roots',len(roots), roots[:5], flush=True)
        if not roots:
            return []
        new=[]
        for xr,yr in roots:
            for ax,ay,M in cand:
                nx,NM=crt_pair(ax,M,xr,pr)
                ny,_=crt_pair(ay,M,yr,pr)
                # reduce representatives
                nx%=NM; ny%=NM
                # optional if modulus exceeds bound, require residue < bound
                if NM>X and nx>=X: continue
                if NM>Y and ny>=Y: continue
                new.append((nx,ny,NM))
                if len(new)>200: break
            if len(new)>200: break
        cand=new
        if verbose: print(' cand',len(cand),'Mbits',cand[0][2].bit_length() if cand else 0, flush=True)
        if not cand: return []
        # if enough recover, verify exact candidates
        if cand[0][2] > max(X,Y):
            return cand
    return cand

def verify_xy(C,xv,yv):
    if not (0<=xv<X and 0<=yv<Y): return False,None
    p=C + (xv<<x_start) + (yv<<y_start)
    if (p & mask) != leak: return False,None
    if N % p == 0:
        return True,p
    return False,None

def get_polys_for_candidate(low,high,m=4,t=3,mode='jm',limit=20):
    C=C_for(low,high)
    if mode=='jm': pairs=jm_pairs(m,t)
    elif mode=='hm': pairs=hm_pairs(m,t)
    elif mode=='box22': pairs=box_pairs(2,2)
    elif mode=='box21': pairs=box_pairs(2,1)
    elif mode=='std': pairs=std_pairs(m)
    else: raise ValueError(mode)
    rows,mons,sh=build_rows(C,m,pairs,True)
    R=reduce_rows(rows,precision=256)
    polys=[]
    for row in R:
        if vector_bits(row)==0: continue
        coeffs=unscale_row(row,mons)
        if coeffs is None: continue
        # skip f^m? We don't know, include but may not vanish. Use bit threshold and not too sparse? 
        polys.append((vector_bits(row), coeffs))
    polys.sort(key=lambda z:z[0])
    return C,polys[:limit]

# CLI: low high mode m t
if __name__=='__main__':
    if len(sys.argv)<3:
        print('usage low high [mode m t]'); sys.exit()
    low=int(sys.argv[1]); high=int(sys.argv[2]); mode=sys.argv[3] if len(sys.argv)>3 else 'jm'; m=int(sys.argv[4]) if len(sys.argv)>4 else 4; t=int(sys.argv[5]) if len(sys.argv)>5 else 3
    C,polys=get_polys_for_candidate(low,high,m,t,mode,limit=25)
    print('polys bits',[b for b,c in polys],flush=True)
    primes=[1000003,1000033,1000037,1000039,1000081,1000099,1000117,1000121,1000133,1000151,1000159,1000171,1000183]
    # use pairs among first polys excluding first maybe try many
    for ia in range(min(12,len(polys))):
      for ib in range(ia+1,min(12,len(polys))):
        b1,c1=polys[ia]; b2,c2=polys[ib]
        print('TRY pair',ia,ib,b1,b2,flush=True)
        cand=recover_with_pair(c1,c2,primes,verbose=True)
        print(' result cand',len(cand), flush=True)
        for xv,yv,M in cand:
            ok,p=verify_xy(C,xv,yv)
            print('  cand bits M',M.bit_length(),'x<',xv.bit_length(),'y<',yv.bit_length(),'ok',ok, flush=True)
            if ok:
                print('FOUND p',hex(p))
                print('plaintext',decrypt_from_p(p))
                sys.exit(0)
    print('not found')

# Fast small-prime routines; coefficients are dict (ix,iy)->int.
def poly_trim(a):
    while len(a)>0 and a[-1]==0: a.pop()
    return a

def poly_mod_mul(a,b,p):
    if not a or not b: return []
    c=[0]*(len(a)+len(b)-1)
    for i,ai in enumerate(a):
        if ai:
            for j,bj in enumerate(b):
                if bj: c[i+j]=(c[i+j]+ai*bj)%p
    return poly_trim(c)

def poly_mod_rem(a,b,p):
    a=a[:]; b=poly_trim(b[:])
    if not b: raise ZeroDivisionError
    db=len(b)-1; inv_l=pow(b[-1],-1,p)
    while len(a)>=len(b) and a:
        coef=a[-1]*inv_l%p
        if coef:
            off=len(a)-len(b)
            for j in range(db+1):
                a[off+j]=(a[off+j]-coef*b[j])%p
        poly_trim(a)
    return a

def poly_mod_gcd(a,b,p):
    a=poly_trim([x%p for x in a]); b=poly_trim([x%p for x in b])
    while b:
        a,b=b,poly_mod_rem(a,b,p)
    if a:
        inv=pow(a[-1],-1,p); a=[(x*inv)%p for x in a]
    return a

def poly_eval_y(poly,y,p):
    s=0
    for c in reversed(poly):
        s=(s*y+c)%p
    return s

def poly_y_at_x(coeffs,xv,p):
    maxj=max((j for i,j in coeffs.keys()), default=0)
    arr=[0]*(maxj+1)
    # powers x
    max_i=max((i for i,j in coeffs.keys()), default=0)
    xp=[1]*(max_i+1)
    for i in range(1,max_i+1): xp[i]=xp[i-1]*xv%p
    for (i,j),c in coeffs.items():
        arr[j]=(arr[j]+(c%p)*xp[i])%p
    return poly_trim(arr)

def roots_pair_mod_fast(c1,c2,p, max_roots=100):
    roots=[]
    for xv in range(p):
        py1=poly_y_at_x(c1,xv,p)
        py2=poly_y_at_x(c2,xv,p)
        if not py1 and not py2:
            continue
        if not py1: g=py2
        elif not py2: g=py1
        else: g=poly_mod_gcd(py1,py2,p)
        if len(g)>1:
            # brute force y for small p
            for yv in range(p):
                if poly_eval_y(g,yv,p)==0:
                    roots.append((xv,yv))
                    if len(roots)>=max_roots: return roots
    return roots

def recover_with_pair_fast(c1,c2, primes, verbose=False):
    cand=[(0,0,1)]
    for pr in primes:
        roots=roots_pair_mod_fast(c1,c2,pr)
        if verbose: print(' prime',pr,'roots',len(roots),roots[:5],flush=True)
        if not roots: return []
        new=[]
        for xr,yr in roots:
            for ax,ay,M in cand:
                nx,NM=crt_pair(ax,M,xr,pr)
                ny,_=crt_pair(ay,M,yr,pr)
                nx%=NM; ny%=NM
                if NM>X and nx>=X: continue
                if NM>Y and ny>=Y: continue
                new.append((nx,ny,NM))
                if len(new)>500: break
            if len(new)>500: break
        cand=new
        if verbose: print(' cand',len(cand),'Mbits',cand[0][2].bit_length() if cand else 0,flush=True)
        if not cand: return []
        if cand[0][2] > max(X,Y): return cand
    return cand

def strip_f_factors(coeffs,C):
    # exact division over QQ by f=C+2^265*x+2^600*y, clear denominators/content
    expr=0
    for (i,j),c in coeffs.items():
        expr += sp.Integer(c)*(x_sym**i)*(y_sym**j)
    P=sp.Poly(expr,x_sym,y_sym,domain=sp.QQ)
    F=sp.Poly(sp.Integer(C) + sp.Integer(1<<x_start)*x_sym + sp.Integer(1<<y_start)*y_sym, x_sym,y_sym,domain=sp.QQ)
    cnt=0
    while True:
        Q,R=P.div(F)
        if R.is_zero:
            P=Q; cnt+=1
        else:
            break
    # clear denominators and content
    dct=P.as_dict()
    if not dct: return {},cnt
    den=1
    for coef in dct.values():
        den=sp.ilcm(den, int(coef.q))
    ints={mon:int(coef*den) for mon,coef in dct.items()}
    # remove gcd content
    g=0
    for val in ints.values(): g=math.gcd(g, abs(val))
    if g>1: ints={mon:val//g for mon,val in ints.items()}
    return {mon:val for mon,val in ints.items() if val},cnt

def polyx_add(p,q,scale=1):
    r=p.copy()
    for i,c in q.items():
        r[i]=r.get(i,0)+scale*c
        if r[i]==0: del r[i]
    return r

def polyx_mul_C_ax(q,C):
    # (C + 2^265 x)*q(x)
    a=1<<x_start
    r={}
    for i,c in q.items():
        r[i]=r.get(i,0)+C*c
        r[i+1]=r.get(i+1,0)+a*c
    return {i:c for i,c in r.items() if c}

def polyx_div_int(p,d):
    r={}
    for i,c in p.items():
        if c%d: return None
        v=c//d
        if v: r[i]=v
    return r

def divide_by_f_fast(coeffs,C):
    b=1<<y_start
    maxj=max((j for i,j in coeffs), default=0)
    if maxj==0:
        return None, False
    P=[{} for _ in range(maxj+1)]
    for (i,j),c in coeffs.items():
        if c: P[j][i]=P[j].get(i,0)+c
    Q=[{} for _ in range(maxj)]
    # top down
    for j in range(maxj,0,-1):
        num=P[j].copy()
        if j<maxj:
            num=polyx_add(num, polyx_mul_C_ax(Q[j],C), scale=-1)  # P_j - L*Q_j
        q=polyx_div_int(num,b)
        if q is None:
            return None, False
        Q[j-1]=q
    rem=polyx_add(P[0].copy(), polyx_mul_C_ax(Q[0],C), scale=-1)
    if rem:
        return None, False
    out={}
    for j,qx in enumerate(Q):
        for i,c in qx.items():
            if c: out[(i,j)]=c
    return out, True

def strip_f_factors_fast(coeffs,C):
    cnt=0; cur=coeffs
    while True:
        q,ok=divide_by_f_fast(cur,C)
        if not ok: break
        cur=q; cnt+=1
    # content gcd
    g=0
    for v in cur.values(): g=math.gcd(g,abs(v))
    if g>1: cur={mon:v//g for mon,v in cur.items()}
    return cur,cnt

def recover_with_pair_fast_limited(c1,c2, primes, verbose=False, max_roots_per_prime=20, max_cand=50):
    cand=[(0,0,1)]
    for pr in primes:
        roots=roots_pair_mod_fast(c1,c2,pr,max_roots=max_roots_per_prime+1)
        if verbose: print(' prime',pr,'roots',len(roots),roots[:5],flush=True)
        if not roots or len(roots)>max_roots_per_prime: return []
        new=[]
        for xr,yr in roots:
            for ax,ay,M in cand:
                nx,NM=crt_pair(ax,M,xr,pr)
                ny,_=crt_pair(ay,M,yr,pr)
                nx%=NM; ny%=NM
                if NM>X and nx>=X: continue
                if NM>Y and ny>=Y: continue
                new.append((nx,ny,NM))
                if len(new)>max_cand: return []
        cand=new
        if verbose: print(' cand',len(cand),'Mbits',cand[0][2].bit_length() if cand else 0,flush=True)
        if not cand: return []
        if cand[0][2] > max(X,Y): return cand
    return cand
