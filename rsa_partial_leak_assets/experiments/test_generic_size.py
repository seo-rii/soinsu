import sys,time,sympy as sp
sys.path.append('/mnt/data')
from hm_bivar import generic_shift_polys_bivar, eval_poly
p=int(sp.randprime(1<<1023,1<<1024)); q=int(sp.randprime(1<<1023,1<<1024)); N=p*q
s1,l1=265,155; s2,l2=600,230; X=1<<l1; Y=1<<l2
x0=(p>>s1)&(X-1); y0=(p>>s2)&(Y-1); mask=((1<<1024)-1)^(((X-1)<<s1)|((Y-1)<<s2)); p0=p&mask
ainv=pow(1<<s1,-1,N); c0=(p0*ainv)%N; cy=((1<<s2)*ainv)%N
print('root bits',x0.bit_length(),y0.bit_length())
for m,d in [(2,6),(3,4),(3,5),(4,4),(4,5),(5,4)]:
 print('\nparam',m,d)
 st=time.time(); hs,_=generic_shift_polys_bivar(N,c0,cy,X,Y,m=m,d=d,lll_rows=8,verbose=True); print('time',time.time()-st,'p^m bits',(p**m).bit_length())
 for i,h in enumerate(hs[:8]):
  v=eval_poly(h,x0,y0); print(i,v==0,abs(v).bit_length(),v%(p**m)==0,abs(v)<p**m)
