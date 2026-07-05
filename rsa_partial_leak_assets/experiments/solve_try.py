from math import prod
from itertools import product
from fpylll import IntegerMatrix, LLL
import sympy as sp

hex_clean=lambda s:int(''.join(s.split()),16)
MASK=hex_clean('''
ffffffffffffffff fffffffff0ffffff ffffffffffffffff c00000000000fffe
0000000000000000 000003ffe0000000 0000000000ffffff ffffffffffffffff
ffffffffffffffff fffffff000000000 000003ffe0000000 00000000000001ff
ffffffffffffffff fffffffffc3fffff ffffffffffffffff ffffffffffffffff
''')
PAND=hex_clean('''
ffa360d46885c534 d186538170633faf c2c0548a2e24a2c1 c0000000000039e2
0000000000000000 000000a520000000 00000000003e2de4 c436d2ca740a6246
99e1a1af94045c63 261323c000000000 000003bba0000000 00000000000000e5
0b0bc2461fcbac07 26360c2c0809450a 9a892cbf1d98ceee 48827591ccc593c9
''')
N=hex_clean('''
e505004fb5d34eb7 12d48ff4bbe8d27f c388133c6c0e7340 01061c0ee0a4edc6
37c04fe8dd376185 de8ba04d0ccdbabb 93ab7c371b88d92e 865eec42b028c61d
d7004ebf2ebb5d69 d0a09142be5c9de4 da16e514eea31817 2ecda6cd192073eb
afb1e02d522ec053 34590ea6d75960c4 937bf64f9700db17 7a4aa3da6aae6807
e5e32c0d0e428a0d b68d299f20c235d8 4ef459b0cf118286 59c31663c9ea8204
4b28152c89a9c36c 3ec4303bd36664fd 77fb02c58340bdae 21120326d83fc017
34bc90048dec9fe3 5f08c8fdc523abf8 4a91ec430f495672 37c3153a2035ff62
5613b6dc3e6cb14d 50e18b8a79b25d67 8465b3ad02f5b7d8 18a1e2d635a0baf1
''')
ct=hex_clean('''
8919342826ef3821 5af31e00c9290c4c 50ef9ff9e1afc591 47fab5b096361035
e85f5fc95b73b069 7813b57b831a807d 41bcbecde5b9e663 9e2845b14e395ed0
e5d995e63709ac0c 5ee2337228ee76bc bad857b14904aa2e 8e9997671908a634
d0d1dda1d062ce7f 2e3293ddec8f5cce 26029292d594a062 dcf317d2a8380f43
d72551889efceb87 6c8945a50382272e 76ed6b6fcdff1603 44e9e948e2b6e740
e78bedf25f30e2c7 eeb5f74686c8eadc 29cea04ff08cfd86 dfd3d2a1632bf04a
d5cfa369892a2da4 0f0dc0098ce6b731 d841aab3d0c8b78e b69c4625c47c4ad7
158d49bb5d879581 e02bc525abe47f39 f699864bc5ce1de7 19430dae7aa5480b
''')
# polynomial dict helper {(i,j): coeff}
def add(p,q):
    r=p.copy()
    for m,c in q.items():
        r[m]=r.get(m,0)+c
        if r[m]==0: del r[m]
    return r
def mul(p,q):
    r={}
    for (i,j),a in p.items():
        for (k,l),b in q.items():
            m=(i+k,j+l); r[m]=r.get(m,0)+a*b
    return {m:c for m,c in r.items() if c}
def powpoly(p,e):
    r={(0,0):1}
    while e:
        if e&1: r=mul(r,p)
        e//=2
        if e: p=mul(p,p)
    return r
def shift_poly(f,k,ix,iy,N,m):
    g=powpoly(f,k)
    if ix or iy: g=mul(g,{(ix,iy):1})
    c=N**(m-k)
    return {mon:coef*c for mon,coef in g.items()}

def build_lattice(f, X, Y, N, m, extra=0, shape='tri'):
    polys=[]
    # base shifts total degree <= m-k (triangular); with extra y shifts? 
    for k in range(m+1):
        D=m-k
        for ix in range(D+1):
            for iy in range(D-ix+1):
                polys.append(shift_poly(f,k,ix,iy,N,m))
    # Optional extra y shifts as HM? multiply f^m by y^j for j=1..extra (and x? maybe)
    for j in range(1,extra+1):
        polys.append(shift_poly(f,m,0,j,N,m))
    mons=sorted(set(mon for p in polys for mon in p), key=lambda z:(z[0]+z[1],z[0],z[1]))
    M=IntegerMatrix(len(polys), len(mons))
    for r,p in enumerate(polys):
        for c,mon in enumerate(mons):
            a,b=mon
            M[r,c]=p.get(mon,0)*(X**a)*(Y**b)
    return M, mons

def row_to_poly(row, mons, X, Y):
    d={}
    for c,mon in enumerate(mons):
        a,b=mon
        val=int(row[c])
        if val:
            den=(X**a)*(Y**b)
            if val%den!=0:
                # should not happen maybe after combination yes still divisible since all rows columns have factor den
                pass
            coeff=val//den
            if coeff: d[mon]=coeff
    return d

def poly_eval(d,x,y): return sum(c*(x**i)*(y**j) for (i,j),c in d.items())
def poly_to_sympy(d):
    x,y=sp.symbols('x y')
    return sum(sp.Integer(c)*x**i*y**j for (i,j),c in d.items())

def try_candidate(low, high, m=6, extra=0, maxrows=10, verbose=False):
    # p0 with guessed low/high and all known bits; x covers 265..419 (incl known gap), y covers 600..829 (incl known gaps)
    p0 = PAND | (low<<150) | (high<<920)
    # remove x,y ranges from p0? PAND already has known bits in those ranges; but f = p0 + 2^265 x + 2^600 y would double known bits inside ranges.
    # So x,y should be full block values incl known bits, therefore p0 must zero those ranges.
    xmask=((1<<155)-1)<<265
    ymask=((1<<230)-1)<<600
    const = p0 & ~xmask & ~ymask
    # x/y include known gap bits; bounds full range
    X=1<<155; Y=1<<230
    f={(0,0):const,(1,0):1<<265,(0,1):1<<600}
    M,mons=build_lattice(f,X,Y,N,m,extra)
    if verbose: print('dim',M.nrows,M.ncols)
    LLL.reduction(M, method='proved', float_type='mpfr', precision=256)
    hs=[]
    for r in range(min(maxrows,M.nrows)):
        d=row_to_poly([M[r,c] for c in range(M.ncols)], mons, X,Y)
        # skip multiples? 
        if d: hs.append(d)
    # resultant pairs
    x,y=sp.symbols('x y')
    syms=[poly_to_sympy(h) for h in hs]
    if verbose:
        for i,h in enumerate(hs[:3]):
            print('h',i,'terms',len(h),'deg',max(a+b for a,b in h),'bits coeff max',max(abs(c).bit_length() for c in h.values()))
    # Try all pairs
    for i in range(len(syms)):
        for j in range(i+1,len(syms)):
            # ensure independent
            try:
                R=sp.resultant(syms[i],syms[j],x) # polynomial in y
                if R==0: continue
                R=sp.Poly(R,y)
                # remove content
                fac=sp.factor_list(R.as_expr())
                candidates_y=[]
                # root via nroots? no. factor list check linear factors
                for factor, exp in fac[1]:
                    P=sp.Poly(factor,y)
                    if P.degree()==1:
                        a,b=P.all_coeffs()
                        if (-b)%a==0:
                            yy=int((-b)//a)
                            if 0<=yy<Y: candidates_y.append(yy)
                    elif P.degree()<=4:
                        roots=sp.roots(P)
                        for rr,mult in roots.items():
                            if rr.is_integer:
                                yy=int(rr)
                                if 0<=yy<Y: candidates_y.append(yy)
                if candidates_y:
                    for yy in set(candidates_y):
                        # solve x from h_i(x,yy)=0
                        Pi=sp.Poly(syms[i].subs(y,yy),x)
                        for rr,mult in sp.roots(Pi).items():
                            if rr.is_integer:
                                xx=int(rr)
                                if 0<=xx<X:
                                    p=const+(xx<<265)+(yy<<600)
                                    if p>1 and N%p==0 and (p&MASK)==PAND:
                                        return p,(xx,yy),'res'
            except Exception as e:
                if verbose: print('res err',i,j,e)
                continue
    return None

if __name__=='__main__':
    print('Nbits',N.bit_length())
    for high in range(16):
      for low in range(16):
        print('try',high,low,flush=True)
        ans=try_candidate(low,high,m=5,extra=0,maxrows=12,verbose=False)
        if ans:
            print('FOUND',high,low,ans[0]); raise SystemExit
    print('none')
