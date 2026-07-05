from fpylll import IntegerMatrix, LLL
from math import isqrt
import time, itertools, sys, random
from collections import defaultdict

H=lambda s:int(''.join(s.split()),16)
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
LEAK=H('''ffa360d46885c534 d186538170633faf c2c0548a2e24a2c1 c0000000000039e2
0000000000000000 000000a520000000 00000000003e2de4 c436d2ca740a6246
99e1a1af94045c63 261323c000000000 000003bba0000000 00000000000000e5
0b0bc2461fcbac07 26360c2c0809450a 9a892cbf1d98ceee 48827591ccc593c9''')
A=1<<265; BB=1<<600; X=1<<155; Y=1<<230
WX=((1<<155)-1)<<265; WY=((1<<230)-1)<<600

def centered(a, mod):
    a %= mod
    if a > mod//2: a -= mod
    return a

def add(p,q):
    r=p.copy()
    for m,c in q.items():
        v=r.get(m,0)+c
        if v: r[m]=v
        elif m in r: del r[m]
    return r

def mul(p,q):
    r={}
    for (a,b),c in p.items():
        for (d,e),f in q.items():
            m=(a+d,b+e); r[m]=r.get(m,0)+c*f
    return {m:c for m,c in r.items() if c}

def shift_y(p,i):
    return {(a,b+i):c for (a,b),c in p.items()}

def scalar(p,s):
    if s==1: return p.copy()
    return {m:c*s for m,c in p.items() if c*s}

def reduce_mod(p,mod):
    return {m:centered(c,mod) for m,c in p.items() if centered(c,mod)}

def eval_poly(p,x,y):
    s=0
    # not efficient
    xp=[1]; yp=[1]
    maxa=max([a for a,b in p]+[0]); maxb=max([b for a,b in p]+[0])
    for i in range(maxa): xp.append(xp[-1]*x)
    for i in range(maxb): yp.append(yp[-1]*y)
    for (a,b),c in p.items(): s += c*xp[a]*yp[b]
    return s

def hm_polys(lo,hi,m=8,t=3,maxrows=8):
    base=LEAK | (lo<<150) | (hi<<920)
    p0=base & ~WX & ~WY
    invA=pow(A,-1,N)
    cy=centered((BB*invA)%N,N)
    c0=centered((p0*invA)%N,N)
    f={(1,0):1,(0,1):cy,(0,0):c0}
    # f powers reduced mod N^k
    fp=[{(0,0):1}]
    cur={(0,0):1}
    for k in range(1,m+1):
        cur=mul(cur,f)
        cur=reduce_mod(cur, pow(N,k))
        fp.append(cur)
    shifts=[]
    for k in range(m+1):
        Nk=pow(N,max(t-k,0))
        for i in range(m-k+1):
            shifts.append(scalar(shift_y(fp[k],i),Nk))
    mons=sorted(set().union(*[set(s.keys()) for s in shifts]), key=lambda ab:(ab[0]+ab[1],ab[0],ab[1]))
    # print len
    scales=[pow(X,a)*pow(Y,b) for a,b in mons]
    M=IntegerMatrix(len(shifts),len(mons))
    for r,s in enumerate(shifts):
        for c,mon in enumerate(mons):
            M[r,c]=s.get(mon,0)*scales[c]
    t0=time.time()
    LLL.reduction(M, delta=0.99, method='proved', float_type='mpfr', precision=256)
    print('LLL time',time.time()-t0,'dim',M.nrows,M.ncols)
    thresh=pow(N,t) # squared threshold N^t (norm2*w < N^t)
    polys=[]
    for r in range(M.nrows):
        row=[int(M[r,c]) for c in range(M.ncols)]
        norm2=sum(v*v for v in row)
        w=sum(1 for v in row if v)
        ok = norm2 * max(w,1) < thresh
        if r<10:
            print('row',r,'norm bits',norm2.bit_length()/2,'w',w,'ok',ok)
        if len(polys)<maxrows: # include first rows regardless? or ok
            poly={}
            bad=False
            for c,val in enumerate(row):
                if val:
                    sc=scales[c]
                    if val%sc: bad=True
                    poly[mons[c]]=val//sc
            if poly and not (len(poly)==1 and (0,0) in poly):
                polys.append(poly)
    return polys, p0

if __name__=='__main__':
    lo=int(sys.argv[1],16) if len(sys.argv)>1 else 0
    hi=int(sys.argv[2],16) if len(sys.argv)>2 else 0
    polys,p0=hm_polys(lo,hi,8,3,10)
    print('polys',len(polys))
    for idx,p in enumerate(polys[:3]):
        print('poly',idx,'terms',len(p),'deg',max(a+b for a,b in p),'coefbits',max(abs(c).bit_length() for c in p.values()))

