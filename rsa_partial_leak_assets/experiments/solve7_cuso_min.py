# run: sage -python solve7_cuso_min.py
# deps: Sage>=9.8, cuso, flatter, msolve
from sage.all import var
import itertools, cuso
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
e=65537; x,y=var('x y')
XS,XL,YS,YL=265,155,600,230
WX=((1<<XL)-1)<<XS; WY=((1<<YL)-1)<<YS
base=LEAK & ~WX & ~WY
B=lambda z:int(z).to_bytes((int(z).bit_length()+7)//8,'big')
for lo,hi in itertools.product(range(16),repeat=2):
    p0=base|(lo<<150)|(hi<<920)
    f=p0+(1<<XS)*x+(1<<YS)*y
    print(f'try lo={lo:x} hi={hi:x}', flush=True)
    for r in cuso.find_small_roots([f], {x:(0,1<<XL), y:(0,1<<YL)}, modulus='p', modulus_multiple=N, modulus_lower_bound=1<<1023, use_graph_optimization=True):
        xv=int(r.get(x,r.get(str(x),0))); yv=int(r.get(y,r.get(str(y),0)))
        p=p0+(xv<<XS)+(yv<<YS)
        if N%p==0 and (p&MASK)==LEAK:
            q=N//p; d=pow(e,-1,(p-1)*(q-1)); m=pow(ct,d,N)
            print('FOUND'); print('p=',hex(p)); print('q=',hex(q)); print('m.hex=',B(m).hex()); print('m=',B(m)); raise SystemExit
print('not found')
