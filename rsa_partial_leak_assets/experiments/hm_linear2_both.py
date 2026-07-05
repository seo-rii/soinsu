import sys, math, time, itertools, sympy as sp
from fpylll import IntegerMatrix, LLL

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
LEAK=H('''ffa360d46885c534 d186538170633faf c2c0548a2e24a2c1 c0000000000039e2
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
e=65537
XS,XL,YS,YL=265,155,600,230
X,Y=1<<XL,1<<YL; A,B=1<<XS,1<<YS
base=LEAK & ~(((1<<XL)-1)<<XS) & ~(((1<<YL)-1)<<YS)

def mul(p,q,mod=None):
    r={}
    for (i,j),a in p.items():
        for (k,l),b in q.items():
            c=a*b
            if mod: c%=mod
            if c:
                e=(i+k,j+l); r[e]=r.get(e,0)+c
    if mod: r={e:c%mod for e,c in r.items() if c%mod}
    return r

def build(C,m,t,both=True):
    polys=[]
    # base normalized in x and y, coefficients reduced modulo N
    bases=[]
    invA=pow(A,-1,N); bases.append({(1,0):1,(0,1):(B*invA)%N,(0,0):(C*invA)%N}) # monic x; shift y^j
    if both:
        invB=pow(B,-1,N); bases.append({(0,1):1,(1,0):(A*invB)%N,(0,0):(C*invB)%N}) # monic y; shift x^j
    for bi,f in enumerate(bases):
        pows=[{(0,0):1}]
        for k in range(1,m+1):
            # Important: kionactf/sage does not reduce f^k mod N; it changes mod coeffs to ZZ then powers.
            # Coefficients are representatives modulo N; multiplication over ZZ.
            pows.append(mul(pows[-1], f, None))
        for k in range(m+1):
            mult=pow(N, max(t-k,0))
            for j in range(m-k+1):
                d={}
                for (u,v),c in pows[k].items():
                    if bi==0: mon=(u,v+j)      # monic x: exclude x, shift y
                    else:     mon=(u+j,v)      # monic y: exclude y, shift x
                    d[mon]=d.get(mon,0)+c*mult
                polys.append(d)
    # dedup rows by dictionary contents (k=0 rows duplicate between bases)
    uniq=[]; seen=set()
    for p in polys:
        key=tuple(sorted(p.items()))
        if key not in seen:
            seen.add(key); uniq.append(p)
    polys=uniq
    mons=sorted(set(e for p in polys for e in p), key=lambda z:(z[0]+z[1],z[0],z[1]))
    scales=[(X**i)*(Y**j) for i,j in mons]
    col={e:i for i,e in enumerate(mons)}
    M=IntegerMatrix(len(polys), len(mons))
    for r,p in enumerate(polys):
        for mon,c in p.items():
            M[r,col[mon]]=int(c*scales[col[mon]])
    return M,mons,scales

def unscale_row(M,r,mons,scales):
    d={}
    for j,mon in enumerate(mons):
        v=int(M[r,j])
        if v:
            s=scales[j]
            # after LLL, divisibility may fail due to linearly combined rows? It should hold over ZZ if rows all scaled same cols.
            if v % s: return None
            d[mon]=v//s
    return d

def eval_poly(poly,xv,yv):
    return sum(c*pow(xv,i)*pow(yv,j) for (i,j),c in poly.items())

def eval_mod(poly,xv,yv,P):
    return sum((c%P)*pow(xv,i,P)*pow(yv,j,P) for (i,j),c in poly.items())%P

# fast brute roots mod small prime, useful when primes ~ 257..1009
_mod_cache={}
def roots_mod(polys,P):
    key=(P, max(i for p in polys for i,j in p), max(j for p in polys for i,j in p))
    if key not in _mod_cache:
        _,mi,mj=key
        xp=[[1]*(mi+1) for _ in range(P)]; yp=[[1]*(mj+1) for _ in range(P)]
        for a in range(P):
            for i in range(1,mi+1): xp[a][i]=xp[a][i-1]*a%P
        for b in range(P):
            for j in range(1,mj+1): yp[b][j]=yp[b][j-1]*b%P
        _mod_cache[key]=(xp,yp)
    xp,yp=_mod_cache[key]
    res=[]
    for a in range(P):
        for b in range(P):
            ok=True
            for poly in polys:
                s=0
                for (i,j),c in poly.items(): s=(s+(c%P)*xp[a][i]*yp[b][j])%P
                if s:
                    ok=False; break
            if ok: res.append((a,b))
    return res

def crt_states(states, residues, P, cap=10000):
    out=[]
    for x,mx,y,my in states:
        invx=pow(mx%P,-1,P); invy=pow(my%P,-1,P)
        for a,b in residues:
            nx=x+mx*(((a-x)%P)*invx%P); ny=y+my*(((b-y)%P)*invy%P)
            nmx=mx*P; nmy=my*P
            # feasible in [0,X), [0,Y)
            if nx<X or nmx<X:
                if ny<Y or nmy<Y:
                    out.append((nx,nmx,ny,nmy))
                    if len(out)>cap: return []
    # dedup
    return list(dict.fromkeys(out))

def recover_from_polys(polys, C, maxpolys=4):
    # choose subsets of first polys; use small primes bruteforce + CRT. Filter known gap bits early.
    # known bit filters within x,y
    xmask=(MASK>>XS)&(X-1); xleak=(LEAK>>XS)&(X-1)
    ymask=(MASK>>YS)&(Y-1); yleak=(LEAK>>YS)&(Y-1)
    def ok_known_partial(v,mod,mask,leak,bound):
        # if mod is power of odd primes, hard to check all bits; only check when modulus exceeds bound
        return True
    primes=[257,263,269,271,277,281,283,293,307,311,313,317,331,337,347,349,353,359,367,373,379,383,389,397,401,409,419,421,431,433,439,443,449,457,461,463,467,479,487,491,499,503,509,521,523,541,547]
    choices=[]
    L=min(len(polys),10)
    # try combinations of 2,3,4 short polynomials, not too many
    for r in (2,3,4):
        for idx in itertools.combinations(range(L),r):
            choices.append(idx)
            if len(choices)>60: break
        if len(choices)>60: break
    for idx in choices:
        use=[polys[i] for i in idx]
        states=[(0,1,0,1)]
        for P in primes:
            R=roots_mod(use,P)
            if not R:
                states=[]; break
            if len(R)>200: # too weak subset
                states=[]; break
            states=crt_states(states,R,P,cap=20000)
            if not states: break
            # if mod covers bounds, test candidates
            if states[0][1]>=X and states[0][3]>=Y:
                for xr,mx,yr,my in states:
                    if xr<X and yr<Y and (xr&xmask)==xleak and (yr&ymask)==yleak:
                        p=C+(xr<<XS)+(yr<<YS)
                        if (p&MASK)==LEAK and N%p==0: return p
                break
    return None

def one(cid,m,t,both=True,keep=24):
    lo,hi=cid&15,cid>>4
    C=base | (lo<<150) | (hi<<920)
    M,mons,scales=build(C,m,t,both)
    st=time.time()
    LLL.reduction(M,delta=0.99,method='proved',float_type='mpfr',precision=256)
    dt=time.time()-st
    hg=1 << (1023*t) # rough bound for filtering only
    rows=[]
    for r in range(M.nrows):
        l1=sum(abs(int(M[r,j])) for j in range(M.ncols))
        p=unscale_row(M,r,mons,scales)
        if p and len(p)>1:
            # optional: keep vectors roughly HG-bound and also nearest shortest
            rows.append((l1.bit_length(),p))
    rows=sorted(rows,key=lambda z:z[0])[:keep]
    bits=[b for b,_ in rows[:6]]
    print(f'cid={cid:03d} lo={lo:x} hi={hi:x} dim={M.nrows}x{M.ncols} LLL={dt:.2f}s bits={bits}', flush=True)
    p=recover_from_polys([p for _,p in rows], C)
    if p:
        return p
    return None

def dec(p):
    q=N//p; d=pow(e,-1,(p-1)*(q-1)); m=pow(ct,d,N)
    return m.to_bytes((m.bit_length()+7)//8,'big')

if __name__=='__main__':
    m=int(sys.argv[1]) if len(sys.argv)>1 else 8
    t=int(sys.argv[2]) if len(sys.argv)>2 else 3
    a=int(sys.argv[3]) if len(sys.argv)>3 else 0
    b=int(sys.argv[4]) if len(sys.argv)>4 else 256
    both=bool(int(sys.argv[5])) if len(sys.argv)>5 else True
    for cid in range(a,b):
        p=one(cid,m,t,both)
        if p:
            q=N//p
            print('FOUND')
            print('p=',hex(p)); print('q=',hex(q)); print('m=',dec(p)); print('mhex=',dec(p).hex())
            sys.exit(0)
    print('not found')
