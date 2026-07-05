import sys,time,math,random
sys.path.append('/mnt/data')
from hm_linear_bivar import build_hm_linear_rows,C_for,row_max_bits,row_l1_bits,unscale_row,N
from fpylll import IntegerMatrix, LLL

def select_independent(rows, p=2147483647):
    basis=[]  # list (pivot, vector mod p)
    selected=[]
    for idx,r in enumerate(rows):
        v=[x%p for x in r]
        for piv,b in basis:
            if v[piv]:
                coef=v[piv]*pow(b[piv],-1,p)%p
                if coef:
                    for j in range(piv,len(v)):
                        v[j]=(v[j]-coef*b[j])%p
        # find pivot
        piv=None
        for j,a in enumerate(v):
            if a:
                piv=j; break
        if piv is not None:
            inv=pow(v[piv],-1,p)
            for j in range(piv,len(v)): v[j]=v[j]*inv%p
            basis.append((piv,v)); selected.append(idx)
            basis.sort(key=lambda x:x[0])
        if len(selected)==len(rows[0]): break
    return selected

for m,t in [(10,3),(13,4)]:
    rows,mons,scales,meta=build_hm_linear_rows(C_for(0,0),m,t,include_dups=False)
    sel=select_independent(rows)
    print('mt',m,t,'rows',len(rows),'cols',len(rows[0]),'rank sel',len(sel),'sel first last',sel[:5],sel[-5:],flush=True)
    srows=[rows[i] for i in sel]
    A=IntegerMatrix(len(srows),len(srows[0]))
    for i,r in enumerate(srows):
     for j,v in enumerate(r): A[i,j]=int(v)
    print('LLL start',flush=True)
    st=time.time()
    LLL.reduction(A,delta=0.99,method='proved',float_type='mpfr',precision=192)
    print('LLL done',time.time()-st,'bits',[row_l1_bits([int(A[i,j]) for j in range(A.ncols)]) for i in range(min(12,A.nrows))], 'hbits', int(math.log2(N)*.499*t), flush=True)
