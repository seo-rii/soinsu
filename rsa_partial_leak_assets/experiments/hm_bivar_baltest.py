import sys,sympy as sp,time
sys.path.append('/mnt/data')
from hm_bivar import iiter, poly_pow, poly_shift, eval_poly
from fpylll import IntegerMatrix, LLL

def hm_bal(N,c0,cy,X,Y,m,t,balance=False):
    def bal(a):
        a%=N
        return a-N if balance and a>N//2 else a
    f={(1,0):1,(0,1):bal(cy),(0,0):bal(c0)}
    g=[]; mons=[]; facs=[]
    for ii in iiter(m,2):
        k,j=ii; g.append(poly_shift(poly_pow(f,k),0,j,pow(N,max(t-k,0)))); mons.append(ii); facs.append(X**k*Y**j)
    n=len(g); M=[[0]*n for _ in range(n)]
    for i,gi in enumerate(g):
        for j in range(i+1): M[i][j]=gi.get(mons[j],0)*facs[j]
    B=IntegerMatrix.from_matrix(M); LLL.reduction(B,delta=0.99,eta=0.51,method='proved',float_type='mpfr',precision=128)
    hs=[]
    for i in range(n):
        h={}
        for j,mon in enumerate(mons):
            v=int(B[i,j])
            if v: h[mon]=v//facs[j]
        hs.append(h)
    return hs
p=int(sp.randprime(1<<1023,1<<1024)); q=int(sp.randprime(1<<1023,1<<1024)); N=p*q
X=1<<155; Y=1<<230; s1=265; s2=600
x0=(p>>s1)&(X-1); y0=(p>>s2)&(Y-1); mask=((1<<1024)-1)^(((X-1)<<s1)|((Y-1)<<s2)); p0=p&mask
inv=pow(1<<s1,-1,N); c0=(p0*inv)%N; cy=((1<<s2)*inv)%N
for bal in [False,True]:
 print('balance',bal); st=time.time(); hs=hm_bal(N,c0,cy,X,Y,8,4,balance=bal); print('time',time.time()-st)
 zeros=[]; minb=99999; pbits=(p**4).bit_length()
 for i,h in enumerate(hs):
  v=eval_poly(h,x0,y0)
  if v==0: zeros.append(i)
  else: minb=min(minb,abs(v).bit_length())
 print('zeros',zeros,'count',len(zeros),'minb',minb,'pbits',pbits)
