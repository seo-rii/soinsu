from fpylll import IntegerMatrix, LLL
from math import prod
import sympy as sp, itertools, time, sys, math, os, pickle

def H(s): return int(''.join(s.split()),16)
N=H('''e505004fb5d34eb7 12d48ff4bbe8d27f c388133c6c0e7340 01061c0ee0a4edc6
37c04fe8dd376185 de8ba04d0ccdbabb 93ab7c371b88d92e 865eec42b028c61d
d7004ebf2ebb5d69 d0a09142be5c9de4 da16e514eea31817 2ecda6cd192073eb
afb1e02d522ec053 34590ea6d75960c4 937bf64f9700db17 7a4aa3da6aae6807
e5e32c0d0e428a0d b68d299f20c235d8 4ef459b0cf118286 59c31663c9ea8204
4b28152c89a9c36c 3ec4303bd36664fd 77fb02c58340bdae 21120326d83fc017
34bc90048dec9fe3 5f08c8fdc523abf8 4a91ec430f495672 37c3153a2035ff62
5613b6dc3e6cb14d 50e18b8a79b25d67 8465b3ad02f5b7d8 18a1e2d635a0baf1''')
MASK=H('''ffffffffffffffff fffffffff0ffffff ffffffffffffffff c00000000000fffe
0000000000000000 000003ffe0000000 0000000000ffffff ffffffffffffffff
ffffffffffffffff fffffff000000000 000003ffe0000000 00000000000001ff
ffffffffffffffff fffffffffc3fffff ffffffffffffffff ffffffffffffffff''')
PM=H('''ffa360d46885c534 d186538170633faf c2c0548a2e24a2c1 c0000000000039e2
0000000000000000 000000a520000000 00000000003e2de4 c436d2ca740a6246
99e1a1af94045c63 261323c000000000 000003bba0000000 00000000000000e5
0b0bc2461fcbac07 26360c2c0809450a 9a892cbf1d98ceee 48827591ccc593c9''')
ct=H('''8919342826ef3821 5af31e00c9290c4c 50ef9ff9e1afc591 47fab5b096361035
e85f5fc95b73b069 7813b57b831a807d 41bcbecde5b9e663 9e2845b14e395ed0
e5d995e63709ac0c 5ee2337228ee76bc bad857b14904aa2e 8e9997671908a634
d0d1dda1d062ce7f 2e3293ddec8f5cce 26029292d594a062 dcf317d2a8380f43
d72551889efceb87 6c8945a50382272e 76ed6b6fcdff1603 44e9e948e2b6e740
e78bedf25f30e2c7 eeb5f74686c8eadc 29cea04ff08cfd86 dfd3d2a1632bf04a
d5cfa369892a2da4 0f0dc0098ce6b731 d841aab3d0c8b78e b69c4625c47c4ad7
158d49bb5d879581 e02bc525abe47f39 f699864bc5ce1de7 19430dae7aa5480b''')

def poly_mul(a,b):
    r={}
    for (i,j),ai in a.items():
        for (k,l),bj in b.items():
            e=(i+k,j+l); r[e]=r.get(e,0)+ai*bj
    return {e:c for e,c in r.items() if c}
def poly_pow(f,k):
    r={(0,0):1}
    for _ in range(k): r=poly_mul(r,f)
    return r

def hm2(low,high,m=6,t=3,keep=8,th_exp=1):
    p0=PM | (low<<150) | (high<<920)
    inv=pow(1<<265,-1,N)
    A=((1<<600)*inv)%N; C=(p0*inv)%N
    f={(1,0):1,(0,1):A,(0,0):C}
    shifts=[]
    for k in range(m+1):
        base=poly_pow(f,k); mult=pow(N,max(t-k,0))
        for iy in range(m+1-k):
            shifts.append({(i,j+iy):c*mult for (i,j),c in base.items()})
    mons=sorted(set(e for g in shifts for e in g), key=lambda e:(e[0]+e[1], e[1], e[0]))
    X,Y=1<<155,1<<230
    B=IntegerMatrix(len(shifts),len(mons))
    scales=[pow(X,i)*pow(Y,j) for i,j in mons]
    for r,g in enumerate(shifts):
        for c,e in enumerate(mons): B[r,c]=g.get(e,0)*scales[c]
    LLL.reduction(B, method='proved', delta=0.99)
    x,y=sp.symbols('x y')
    polys=[]
    for r in range(min(keep,B.nrows)):
        expr=0; ns=0; w=0
        for c,(i,j) in enumerate(mons):
            v=int(B[r,c]);
            if v:
                ns += v*v; w+=1; expr += (v//scales[c])*(x**i)*(y**j)
        if expr and (not expr.is_number):
            # optionally filter norm <? use loose
            polys.append(sp.Poly(expr,x,y,domain=sp.ZZ))
    # Try groebner on incremental subsets
    for s in range(2,min(len(polys),keep)+1):
        try:
            G=sp.groebner([p.as_expr() for p in polys[:s]], x,y, order='lex')
            # print('GB',len(G.polys), [g.total_degree() for g in G.polys[:3]])
            # find univar in y or x
            for g in G.polys:
                if len(g.gens)==2:
                    expr=g.as_expr()
                    if not expr.has(x) and expr.has(y):
                        roots=sp.roots(sp.Poly(expr,y))
                        for ry,mul in roots.items():
                            if ry.is_Integer and 0 <= int(ry) < Y:
                                # solve x from f mod? Need use another poly linear in x
                                yy=int(ry)
                                for gg in G.polys:
                                    ex=sp.Poly(gg.as_expr().subs(y,yy), x)
                                    if ex.degree()==1:
                                        aa=int(ex.nth(1)); bb=int(ex.nth(0))
                                        if aa and (-bb)%aa==0:
                                            xx=(-bb)//aa
                                            if 0<=xx<X:
                                                p=p0+(xx<<265)+(yy<<600)
                                                if (p&MASK)==PM and N%p==0: return p
                    if not expr.has(y) and expr.has(x):
                        roots=sp.roots(sp.Poly(expr,x))
                        for rx,mul in roots.items():
                            if rx.is_Integer and 0<=int(rx)<X:
                                xx=int(rx)
                                for gg in G.polys:
                                    ey=sp.Poly(gg.as_expr().subs(x,xx), y)
                                    if ey.degree()==1:
                                        aa=int(ey.nth(1)); bb=int(ey.nth(0))
                                        if aa and (-bb)%aa==0:
                                            yy=(-bb)//aa
                                            if 0<=yy<Y:
                                                p=p0+(xx<<265)+(yy<<600)
                                                if (p&MASK)==PM and N%p==0: return p
        except Exception as e:
            pass
    # pairwise resultants over ZZ maybe
    for a,b in itertools.combinations(polys[:keep],2):
        try:
            R=sp.resultant(a.as_expr(),b.as_expr(), y)
            P=sp.Poly(R,x,domain=sp.ZZ)
            if P.degree()>0 and P.degree()<=20:
                for rx,mul in sp.roots(P).items():
                    if rx.is_Integer and 0<=int(rx)<X:
                        xx=int(rx)
                        # find y linear
                        for poly in (a,b):
                            Py=sp.Poly(poly.as_expr().subs(x,xx),y)
                            for ry,mul2 in sp.roots(Py).items():
                                if ry.is_Integer and 0<=int(ry)<Y:
                                    yy=int(ry); p=p0+(xx<<265)+(yy<<600)
                                    if (p&MASK)==PM and N%p==0: return p
        except Exception as e:
            pass
    return None

if __name__=='__main__':
    print('bits',N.bit_length(),bin(MASK).count('1'))
    order=[(h,l) for h in range(16) for l in range(16)]
    # maybe scan fast params
    for h,l in order:
        st=time.time();
        p=hm2(l,h,m=int(sys.argv[1]) if len(sys.argv)>1 else 6,t=int(sys.argv[2]) if len(sys.argv)>2 else 3,keep=10)
        print('try',h,l,'dt',round(time.time()-st,2),'hit',p is not None, flush=True)
        if p:
            print(hex(p)); print(hex(N//p)); break
