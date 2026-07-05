import math, time, itertools, os, sys
from fpylll import IntegerMatrix, LLL, BKZ

N = int('''e505004fb5d34eb712d48ff4bbe8d27fc388133c6c0e734001061c0ee0a4edc6
37c04fe8dd376185de8ba04d0ccdbabb93ab7c371b88d92e865eec42b028c61d
d7004ebf2ebb5d69d0a09142be5c9de4da16e514eea318172ecda6cd192073eb
afb1e02d522ec05334590ea6d75960c4937bf64f9700db177a4aa3da6aae6807
e5e32c0d0e428a0db68d299f20c235d84ef459b0cf11828659c31663c9ea8204
4b28152c89a9c36c3ec4303bd36664fd77fb02c58340bdae21120326d83fc017
34bc90048dec9fe35f08c8fdc523abf84a91ec430f49567237c3153a2035ff62
5613b6dc3e6cb14d50e18b8a79b25d678465b3ad02f5b7d818a1e2d635a0baf1'''.replace('\n',''),16)
ct = int('''8919342826ef38215af31e00c9290c4c50ef9ff9e1afc59147fab5b096361035
e85f5fc95b73b0697813b57b831a807d41bcbecde5b9e6639e2845b14e395ed0
e5d995e63709ac0c5ee2337228ee76bcbad857b14904aa2e8e9997671908a634
d0d1dda1d062ce7f2e3293ddec8f5cce26029292d594a062dcf317d2a8380f43
d72551889efceb876c8945a50382272e76ed6b6fcdff160344e9e948e2b6e740
e78bedf25f30e2c7eeb5f74686c8eadc29cea04ff08cfd86dfd3d2a1632bf04a
d5cfa369892a2da40f0dc0098ce6b731d841aab3d0c8b78eb69c4625c47c4ad7
158d49bb5d879581e02bc525abe47f39f699864bc5ce1de719430dae7aa5480b'''.replace('\n',''),16)
mask = int(''.join([
'ffffffffffffffff','fffffffff0ffffff','ffffffffffffffff','c00000000000fffe',
'0000000000000000','000003ffe0000000','0000000000ffffff','ffffffffffffffff',
'ffffffffffffffff','fffffff000000000','000003ffe0000000','00000000000001ff',
'ffffffffffffffff','fffffffffc3fffff','ffffffffffffffff','ffffffffffffffff']),16)
leak = int(''.join([
'ffa360d46885c534','d186538170633faf','c2c0548a2e24a2c1','c0000000000039e2',
'0000000000000000','000000a520000000','00000000003e2de4','c436d2ca740a6246',
'99e1a1af94045c63','261323c000000000','000003bba0000000','00000000000000e5',
'0b0bc2461fcbac07','26360c2c0809450a','9a892cbf1d98ceee','48827591ccc593c9']),16)
e = 65537
x_start,x_len=265,155
y_start,y_len=600,230
X=1<<x_len; Y=1<<y_len

def clear_range(val,start,length):
    return val & ~(((1<<length)-1)<<start)
base0=clear_range(clear_range(leak,x_start,x_len),y_start,y_len)

def C_for(low, high):
    return base0 | (low<<150) | (high<<920)

def actual_xy_from_p(p):
    return ((p>>x_start)&(X-1), (p>>y_start)&(Y-1))

def poly_shift_terms(C,a,b,k,sx,sy,Npow):
    terms={}
    fact=math.factorial
    for i in range(k+1):
        for j in range(k-i+1):
            cpower=k-i-j
            coeff=fact(k)//(fact(i)*fact(j)*fact(cpower))
            coeff *= pow(C,cpower)*pow(a,i)*pow(b,j)*Npow
            terms[(i+sx,j+sy)] = terms.get((i+sx,j+sy),0)+coeff
    return terms

def build_rows(C,m,shift_pairs,include_k0=True):
    a=1<<x_start; b=1<<y_start
    rows_terms=[]; shifts=[]
    for k in range(0 if include_k0 else 1,m+1):
        Npow=pow(N,m-k)
        for sx,sy in shift_pairs(k):
            rows_terms.append(poly_shift_terms(C,a,b,k,sx,sy,Npow)); shifts.append((k,sx,sy))
    mons=sorted(set().union(*[set(t.keys()) for t in rows_terms]), key=lambda e:(e[0]+e[1],e[0],e[1]))
    col={mon:i for i,mon in enumerate(mons)}
    # precompute scale maybe
    scale={mon: pow(X,mon[0])*pow(Y,mon[1]) for mon in mons}
    rows=[]
    for terms in rows_terms:
        row=[0]*len(mons)
        for mon,coef in terms.items():
            row[col[mon]]=coef*scale[mon]
        rows.append(row)
    return rows,mons,shifts

def box_pairs(tx,ty):
    return lambda k: [(sx,sy) for sx in range(tx+1) for sy in range(ty+1)]

def tri_pairs(T):
    return lambda k: [(sx,sy) for sx in range(T+1) for sy in range(T+1-sx)]

def std_pairs(m):
    return lambda k: [(sx,sy) for sx in range(m-k+1) for sy in range(m-k+1-sx)]

def hm_pairs(m,t):
    # a common rectangular set: all shifts total degree <= t for every f^k
    return lambda k: [(sx,sy) for sx in range(t+1) for sy in range(t+1-sx)]

def reduce_rows(rows, delta=0.99, method='proved', precision=256, flags=0):
    M=len(rows); n=len(rows[0])
    A=IntegerMatrix(M,n)
    for i,row in enumerate(rows):
        for j,v in enumerate(row):
            A[i,j]=int(v)
    t=time.time()
    LLL.reduction(A, delta=delta, method=method, float_type='mpfr', precision=precision)
    print(f'LLL {M}x{n} {time.time()-t:.3f}s', file=sys.stderr, flush=True)
    out=[]
    for i in range(M):
        out.append([int(A[i,j]) for j in range(n)])
    return out

def vector_bits(row):
    vals=[abs(v).bit_length() for v in row if v]
    return max(vals) if vals else 0

def unscale_row(row, mons):
    coeffs={}
    for val,mon in zip(row, mons):
        val=int(val)
        if val:
            den=pow(X,mon[0])*pow(Y,mon[1])
            if val % den != 0:
                # keep rational? shouldn't happen
                return None
            coeffs[mon]=val//den
    return coeffs

def eval_poly(coeffs,x,y):
    s=0
    # powers cache
    xp={0:1}; yp={0:1}
    for i,j in coeffs:
        if i not in xp: xp[i]=pow(x,i)
        if j not in yp: yp[j]=pow(y,j)
    for (i,j),c in coeffs.items():
        s += c*xp[i]*yp[j]
    return s

def decrypt_from_p(p):
    q=N//p
    phi=(p-1)*(q-1)
    d=pow(e,-1,phi)
    m=pow(ct,d,N)
    return m.to_bytes((m.bit_length()+7)//8,'big')

def jm_pairs(m,t):
    # base triangular for k<m: total shift <= m-k; for k=m: extra total shifts <= t
    def pairs(k):
        T = (m-k) if k < m else t
        return [(sx,sy) for sx in range(T+1) for sy in range(T+1-sx)]
    return pairs

def jm_box_pairs(m,tx,ty):
    def pairs(k):
        if k < m:
            return [(sx,sy) for sx in range(m-k+1) for sy in range(m-k+1-sx)]
        else:
            return [(sx,sy) for sx in range(tx+1) for sy in range(ty+1)]
    return pairs
