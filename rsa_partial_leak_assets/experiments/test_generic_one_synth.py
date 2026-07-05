import sys,time,random,sympy as sp,gc
sys.path.append('/mnt/data')
from hm_bivar import generic_shift_polys_bivar, eval_poly
p=int(sp.randprime(1<<1023,1<<1024)); q=int(sp.randprime(1<<1023,1<<1024)); N=p*q
X=1<<155; Y=1<<230; s1=265; s2=600
x0=(p>>s1)&(X-1); y0=(p>>s2)&(Y-1); mask=((1<<1024)-1)^(((X-1)<<s1)|((Y-1)<<s2)); p0=p&mask
inv=pow(1<<s1,-1,N); c0=(p0*inv)%N; cy=((1<<s2)*inv)%N
m,d=3,4
st=time.time(); hs,_=generic_shift_polys_bivar(N,c0,cy,X,Y,m=m,d=d,lll_rows=43,verbose=True); print('time',time.time()-st)
zeros=[]; minb=99999; pmb=(p**m).bit_length()
for i,h in enumerate(hs):
 v=eval_poly(h,x0,y0)
 if v==0: zeros.append(i)
 else: minb=min(minb,abs(v).bit_length())
print('zeros',zeros,'count',len(zeros),'minb',minb,'pmb',pmb)
