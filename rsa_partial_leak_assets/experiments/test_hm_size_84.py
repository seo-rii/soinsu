import sys,time,random,sympy as sp
sys.path.append('/mnt/data')
from hm_bivar import hm_polys_bivar, eval_poly
p=int(sp.randprime(1<<1023,1<<1024)); q=int(sp.randprime(1<<1023,1<<1024)); N=p*q
X=1<<155; Y=1<<230; s1=265; s2=600
x0=(p>>s1)&(X-1); y0=(p>>s2)&(Y-1); mask=((1<<1024)-1)^(((X-1)<<s1)|((Y-1)<<s2)); p0=p&mask
ainv=pow(1<<s1,-1,N); c0=(p0*ainv)%N; cy=((1<<s2)*ainv)%N
print('gen')
st=time.time(); hs,_=hm_polys_bivar(N,c0,cy,X,Y,m=8,t=4,lll_rows=10,verbose=True); print('time',time.time()-st,'p^t bits',(p**4).bit_length())
for i,h in enumerate(hs[:10]):
 v=eval_poly(h,x0,y0)
 print(i,v==0,abs(v).bit_length(),v%(p**4)==0,abs(v)<p**4)
