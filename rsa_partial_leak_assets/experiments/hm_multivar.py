from fpylll import IntegerMatrix, LLL, BKZ, GSO
import math

# sparse bivariate polynomials dict (i,j)->coeff over ZZ

def clean(p): return {m:c for m,c in p.items() if c}
def add(a,b):
    r=a.copy()
    for m,c in b.items(): r[m]=r.get(m,0)+c
    return clean(r)
def mul(a,b):
    r={}
    for (i,j),c in a.items():
        if not c: continue
        for (k,l),d in b.items():
            if not d: continue
            r[(i+k,j+l)] = r.get((i+k,j+l),0)+c*d
    return clean(r)
def pow_poly(f,k):
    r={(0,0):1}
    for _ in range(k): r=mul(r,f)
    return r
def shift_y(p,j): return {(i,jj+j):c for (i,jj),c in p.items()}
def mod_coeffs(p,N,center=False):
    r={}
    for m,c in p.items():
        c%=N
        if center and c>N//2: c-=N
        r[m]=c
    return clean(r)
def eval_poly(p,x,y):
    return sum(c*(x**i)*(y**j) for (i,j),c in p.items())

def hm_shifts(f,N,m,t):
    # f must be normalized monic in x modulo N, coefficients integers usually [0,N)
    shifts=[]
    fpows=[pow_poly(f,k) for k in range(m+1)]
    for k in range(m+1):
        base=fpows[k]
        if k<t:
            fac=pow(N,t-k)
            base={mon:c*fac for mon,c in base.items()}
        for j in range(m+1-k):
            shifts.append(shift_y(base,j))
    return shifts

def create_lattice(shifts,X,Y):
    mons=sorted(set().union(*[set(s) for s in shifts]), key=lambda z:(z[0]+z[1], z[0], z[1]))
    M=IntegerMatrix(len(shifts),len(mons))
    for r,s in enumerate(shifts):
        for c,mon in enumerate(mons):
            coeff=s.get(mon,0)
            if coeff:
                i,j=mon
                M[r,c]=coeff*pow(X,i)*pow(Y,j)
    return M,mons

def reconstruct_row(M,r,mons,X,Y):
    p={}
    for c,mon in enumerate(mons):
        val=int(M[r,c])
        if val:
            i,j=mon
            den=pow(X,i)*pow(Y,j)
            assert val%den==0
            p[mon]=val//den
    return clean(p)

def reduce_lattice(M,method='lll'):
    if method=='lll':
        LLL.reduction(M, method='proved', float_type='mpfr', precision=256)
    elif method=='bkz':
        LLL.reduction(M, method='proved', float_type='mpfr', precision=256)
        par=BKZ.Param(block_size=20, strategies=BKZ.DEFAULT_STRATEGY, max_loops=4)
        BKZ.reduction(M, par)
    return M

def row_norm_bits(M,r):
    s=0
    for c in range(M.ncols):
        v=int(M[r,c]); s+=v*v
    return s.bit_length()/2 if s else 0

# validation toy
if __name__=='__main__':
    import sympy as sp, random
    for bits in [128,192,256,384]:
        p=int(sp.randprime(1<<(bits//2-1),1<<(bits//2))); q=int(sp.randprime(1<<(bits//2-1),1<<(bits//2))); N=p*q
        # bivariate windows: x 20%, y 15% of pbits product maybe under bound
        a=bits//8; bx=bits//8; b=bits//4; by=bits//8
        mask=((1<<bx)-1)<<a | ((1<<by)-1)<<b
        p0=p & ~mask; x0=(p>>a)&((1<<bx)-1); y0=(p>>b)&((1<<by)-1)
        A=1<<a; B=1<<b
        invA=pow(A,-1,N); C=(B*invA)%N; D=(p0*invA)%N
        f={(1,0):1,(0,1):C,(0,0):D}
        X=1<<bx; Y=1<<by
        print('\ntoy bits',bits,'xbits',bx,'ybits',by,'pbits',p.bit_length())
        for m,t in [(3,1),(4,1),(5,1),(4,2),(5,2),(6,2),(8,3),(10,3),(10,4)]:
            shifts=hm_shifts(f,N,m,t); M,mons=create_lattice(shifts,X,Y)
            reduce_lattice(M,'lll')
            good=0; vals=[]
            for r in range(min(M.nrows,12)):
                pol=reconstruct_row(M,r,mons,X,Y); val=eval_poly(pol,x0,y0)
                if val==0: good+=1
                vals.append((0 if val==0 else abs(val).bit_length(), val==0, len(pol), row_norm_bits(M,r)))
            print('m,t,dim',m,t,M.nrows,M.ncols,'good',good,'first',vals[:4])
            if good>=2: break
