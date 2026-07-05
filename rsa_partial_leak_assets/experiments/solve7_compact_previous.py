#!/usr/bin/env python3
# Preferred: sage -python solve7_compact.py --mode cuso
# Fallback:  python3 solve7_compact.py --mode local --m 8 --t 3 --a 0 --b 256
import argparse, time, sys
H=lambda s:int(''.join(s.split()),16)
N=H('''e505004fb5d34eb7 12d48ff4bbe8d27f c388133c6c0e7340 01061c0ee0a4edc6
37c04fe8dd376185 de8ba04d0ccdbabb 93ab7c371b88d92e 865eec42b028c61d
d7004ebf2ebb5d69 d0a09142be5c9de4 da16e514eea31817 2ecda6cd192073eb
afb1e02d522ec053 34590ea6d75960c4 937bf64f9700db17 7a4aa3da6aae6807
e5e32c0d0e428a0d b68d299f20c235d8 4ef459b0cf118286 59c31663c9ea8204
4b28152c89a9c36c 3ec4303bd36664fd 77fb02c58340bdae 21120326d83fc017
34bc90048dec9fe3 5f08c8fdc523abf8 4a91ec430f495672 37c3153a2035ff62
5613b6dc3e6cb14d 50e18b8a79b25d67 8465b3ad02f5b7d8 18a1e2d635a0baf1''')
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
e=65537; XS,XL,YS,YL=265,155,600,230; X,Y=1<<XL,1<<YL; A,B=1<<XS,1<<YS
BASE=LEAK & ~(((1<<XL)-1)<<XS) & ~(((1<<YL)-1)<<YS)
def bb(n): return n.to_bytes((n.bit_length()+7)//8,'big')
def testp(p):
    if 1<p<N and N%p==0 and (p&MASK)==LEAK:
        q=N//p; m=pow(ct,pow(e,-1,(p-1)*(q-1)),N)
        print('FOUND'); print('p =',hex(p)); print('q =',hex(q)); print('m.hex =',bb(m).hex()); print('m =',bb(m)); return True
    return False

def cuso_run(a,b):
    from sage.all import var
    import cuso
    x,y=var('x y')
    def get(r,k):
        for kk in (k,str(k)):
            try:
                if kk in r: return int(r[kk])
            except Exception: pass
        return None
    for cid in range(a,b):
        lo,hi=cid&15,cid>>4; C=BASE|(lo<<150)|(hi<<920)
        print(f'try cid={cid} lo={lo:x} hi={hi:x}',flush=True)
        f=C+(1<<XS)*x+(1<<YS)*y
        for r in cuso.find_small_roots([f],{x:(0,X),y:(0,Y)},modulus='p',modulus_multiple=N,modulus_lower_bound=1<<1023,modulus_upper_bound=1<<1024):
            p=int(r['p']) if isinstance(r,dict) and 'p' in r else C+(get(r,x)<<XS)+(get(r,y)<<YS)
            if testp(p): return
    print('not found')

def mul(P,Q,mod=0):
    R={}
    for (i,j),a in P.items():
        for (k,l),b in Q.items(): R[(i+k,j+l)]=R.get((i+k,j+l),0)+a*b
    return {m:(c%mod if mod else c) for m,c in R.items() if (c%mod if mod else c)}
def build(C,m,t,lead):
    from fpylll import IntegerMatrix
    if lead=='y': f={(0,1):1,(1,0):(A*pow(B,-1,N))%N,(0,0):(C*pow(B,-1,N))%N}; sh=(1,0)
    else:         f={(1,0):1,(0,1):(B*pow(A,-1,N))%N,(0,0):(C*pow(A,-1,N))%N}; sh=(0,1)
    P=[]; pw=[{(0,0):1}]
    for k in range(1,m+1): pw.append(mul(pw[-1],f,N**k))
    for k in range(m+1):
        for i in range(m-k+1): P.append({(a+sh[0]*i,b+sh[1]*i):c*N**max(t-k,0) for (a,b),c in pw[k].items()})
    mons=sorted(set().union(*(p.keys() for p in P)),key=lambda z:(z[0]+z[1],z[0])); col={u:i for i,u in enumerate(mons)}; sc=[X**i*Y**j for i,j in mons]
    M=IntegerMatrix(len(P),len(mons))
    for r,p in enumerate(P):
        for mon,c in p.items(): M[r,col[mon]]=int(c*sc[col[mon]])
    return M,mons,sc
def unscale(row,mons,sc):
    P={}
    for v,mon,s in zip(row,mons,sc):
        v=int(v)
        if v:
            if v%s: return None
            P[mon]=v//s
    return P
def ev(P,x,y,p): return sum((c%p)*pow(x,i,p)*pow(y,j,p) for (i,j),c in P.items())%p
def roots_mod(P,pr):
    return [(x,y) for x in range(pr) for y in range(pr) if all(ev(F,x,y,pr)==0 for F in P)]
def crt_join(st,mod,R,pr,cap=30000):
    out=[]; inv=pow(mod,-1,pr)
    for x,y in st:
        for a,b in R:
            out.append((x+mod*(((a-x)%pr)*inv%pr),y+mod*(((b-y)%pr)*inv%pr)))
            if len(out)>cap: return [],0
    return out,mod*pr
def recover(P,k=4):
    st=[(0,0)]; mod=1
    for pr in [257,263,269,271,277,281,283,293,307,311,313,317,331,337,347,349,353,359,367,373,379,383,389,397,401,409,419,421,431,433,439,443,449,457,461,463,467]:
        R=roots_mod(P[:k],pr)
        if not R: return []
        st,mod=crt_join(st,mod,R,pr)
        if not st: return []
        if mod>X and mod>Y: return [(x,y) for x,y in st if x<X and y<Y]
    return []
def local_run(a,b,m,t,lead,delta,prec):
    from fpylll import LLL
    for cid in range(a,b):
        lo,hi=cid&15,cid>>4; C=BASE|(lo<<150)|(hi<<920)
        M,mons,sc=build(C,m,t,lead); t0=time.time()
        LLL.reduction(M,delta=delta,method='proved',float_type='mpfr',precision=prec)
        rows=[]
        for r in range(M.nrows):
            row=[int(M[r,j]) for j in range(M.ncols)]; P=unscale(row,mons,sc)
            if P and len(P)>1: rows.append((sum(abs(v) for v in row).bit_length(),P))
        rows=sorted(rows,key=lambda z:z[0])[:10]
        print(cid,lo,hi,'lead',lead,'dim',M.nrows,'LLL',round(time.time()-t0,2),'bits',[b for b,_ in rows[:5]],flush=True)
        for x,y in recover([p for _,p in rows]):
            if testp(C+(x<<XS)+(y<<YS)): return
    print('not found')

if __name__=='__main__':
    p=argparse.ArgumentParser(); p.add_argument('--mode',choices=['auto','cuso','local'],default='auto'); p.add_argument('--m',type=int,default=8); p.add_argument('--t',type=int,default=3); p.add_argument('--a',type=int,default=0); p.add_argument('--b',type=int,default=256); p.add_argument('--lead',choices=['x','y'],default='y'); p.add_argument('--delta',type=float,default=.99); p.add_argument('--prec',type=int,default=256); a=p.parse_args()
    if a.mode in ('auto','cuso'):
        try: cuso_run(a.a,a.b)
        except Exception as ex:
            if a.mode=='cuso': raise
            print('cuso unavailable/fail:',type(ex).__name__,ex,'=> local fallback')
            local_run(a.a,a.b,a.m,a.t,a.lead,a.delta,a.prec)
    else: local_run(a.a,a.b,a.m,a.t,a.lead,a.delta,a.prec)
