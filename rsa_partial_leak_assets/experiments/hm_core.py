from fpylll import IntegerMatrix, LLL
from collections import defaultdict
import time
H=lambda s:int(''.join(s.split()),16)
N=H('''e505004fb5d34eb7 12d48ff4bbe8d27f c388133c6c0e7340 01061c0ee0a4edc6
37c04fe8dd376185 de8ba04d0ccdbabb 93ab7c371b88d92e 865eec42b028c61d
d7004ebf2ebb5d69 d0a09142be5c9de4 da16e514eea31817 2ecda6cd192073eb
afb1e02d522ec053 34590ea6d75960c4 937bf64f9700db17 7a4aa3da6aae6807
e5e32c0d0e428a0d b68d299f20c235d8 4ef459b0cf118286 59c31663c9ea8204
4b28152c89a9c36c 3ec4303bd36664fd 77fb02c58340bdae 21120326d83fc017
34bc90048dec9fe3 5f08c8fdc523abf8 4a91ec430f495672 37c3153a2035ff62
5613b6dc3e6cb14d 50e18b8a79b25d67 8465b3ad02f5b7d8 18a1e2d635a0baf1''')
e=65537
ct=H('''8919342826ef3821 5af31e00c9290c4c 50ef9ff9e1afc591 47fab5b096361035
e85f5fc95b73b069 7813b57b831a807d 41bcbecde5b9e663 9e2845b14e395ed0
e5d995e63709ac0c 5ee2337228ee76bc bad857b14904aa2e 8e9997671908a634
d0d1dda1d062ce7f 2e3293ddec8f5cce 26029292d594a062 dcf317d2a8380f43
d72551889efceb87 6c8945a50382272e 76ed6b6fcdff1603 44e9e948e2b6e740
e78bedf25f30e2c7 eeb5f74686c8eadc 29cea04ff08cfd86 dfd3d2a1632bf04a
d5cfa369892a2da4 0f0dc0098ce6b731 d841aab3d0c8b78e b69c4625c47c4ad7
158d49bb5d879581 e02bc525abe47f39 f699864bc5ce1de7 19430dae7aa5480b''')
MASK=H('''ffffffffffffffff fffffffff0ffffff ffffffffffffffff c00000000000fffe
0000000000000000 000003ffe0000000 0000000000ffffff ffffffffffffffff
ffffffffffffffff fffffff000000000 000003ffe0000000 00000000000001ff
ffffffffffffffff fffffffffc3fffff ffffffffffffffff ffffffffffffffff''')
LEAK=H('''ffa360d46885c534 d186538170633faf c2c0548a2e24a2c1 c0000000000039e2
0000000000000000 000000a520000000 00000000003e2de4 c436d2ca740a6246
99e1a1af94045c63 261323c000000000 000003bba0000000 00000000000000e5
0b0bc2461fcbac07 26360c2c0809450a 9a892cbf1d98ceee 48827591ccc593c9''')
A=1<<265; B=1<<600; X=1<<155; Y=1<<230
WX=((1<<155)-1)<<265; WY=((1<<230)-1)<<600
invA=pow(A,-1,N)

def centered(a,mod):
    a%=mod
    return a-mod if a>mod//2 else a

def mul(p,q):
    r={}
    for (a,b),c in p.items():
        for (d,e),f in q.items():
            m=(a+d,b+e); r[m]=r.get(m,0)+c*f
    return {m:c for m,c in r.items() if c}

def shift_y(p,i): return {(a,b+i):c for (a,b),c in p.items()}
def scalar(p,s): return p.copy() if s==1 else {m:c*s for m,c in p.items() if c*s}
def reduce_mod(p,mod):
    r={}
    for m,c in p.items():
        cc=centered(c,mod)
        if cc: r[m]=cc
    return r

def hm_polys(lo,hi,m=8,t=3,maxrows=10, delta=0.99, verbose=False):
    base=LEAK | (lo<<150) | (hi<<920)
    p0=base & ~WX & ~WY
    cy=centered((B*invA)%N,N); c0=centered((p0*invA)%N,N)
    f={(1,0):1,(0,1):cy,(0,0):c0}
    fp=[{(0,0):1}]; cur={(0,0):1}
    Npows=[1]
    for k in range(1,m+1): Npows.append(Npows[-1]*N)
    for k in range(1,m+1):
        cur=mul(cur,f); cur=reduce_mod(cur,Npows[k]); fp.append(cur)
    shifts=[]
    for k in range(m+1):
        Nk=Npows[max(t-k,0)]
        for i in range(m-k+1): shifts.append(scalar(shift_y(fp[k],i),Nk))
    mons=sorted(set().union(*(set(s.keys()) for s in shifts)), key=lambda ab:(ab[0]+ab[1],ab[0],ab[1]))
    scales=[(X**a)*(Y**b) for a,b in mons]
    M=IntegerMatrix(len(shifts),len(mons))
    for r,s in enumerate(shifts):
        for c,mon in enumerate(mons): M[r,c]=s.get(mon,0)*scales[c]
    t0=time.time()
    LLL.reduction(M, delta=delta, method='proved', float_type='mpfr', precision=256)
    lll_time=time.time()-t0
    polys=[]; stats=[]
    thresh=Npows[t]
    for r in range(M.nrows):
        row=[int(M[r,c]) for c in range(M.ncols)]
        norm2=sum(v*v for v in row); w=sum(1 for v in row if v)
        stats.append((norm2.bit_length()/2,w,norm2*max(w,1)<thresh))
        if len(polys)<maxrows:
            poly={}
            for c,val in enumerate(row):
                if val:
                    poly[mons[c]]=val//scales[c]
            if poly and not (len(poly)==1 and (0,0) in poly): polys.append(poly)
    if verbose: print('LLL',lll_time,'dim',M.nrows,'stats',stats[:5])
    return polys,p0,lll_time,stats

def decrypt(p):
    q=N//p; phi=(p-1)*(q-1); d=pow(e,-1,phi); m=pow(ct,d,N)
    return m.to_bytes((m.bit_length()+7)//8,'big')
